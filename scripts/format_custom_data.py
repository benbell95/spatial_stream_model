#!/usr/bin/env python3
"""
scripts/format_custom_data.py

Universal Standalone Custom Monitoring Data Formatter for the S3N Stream Network Pipeline.

Transforms any external, user-supplied, or citizen-science monitoring dataset into the
canonical `sites.csv` and `samples.csv` schemas required by `01_preprocess_orn.py` and
`run_s3n_determinand.R`.

Key Capabilities:
  1. Universal Auto-Detection:
     Scans column headers and candidate aliases for station identifiers, British National
     Grid (BNG) or WGS84 coordinates, dates/times, determinand codes, concentrations,
     censoring flags, and measurement units.
  2. Explicit Column Overrides:
     Allows overriding any or all column mappings via command-line options
     (--id-col, --x-col, --y-col, --date-col, --time-col, --val-col, --det-col, etc.).
  3. "Gold Standard" Interval-Censored Tobit & Truncated Lognormal Model (Default):
     Test-strip and citizen-science datasets frequently report observations as intervals
     or ranges (e.g. '5-10', '0.2-0.5', '1-2', '<0.2', '>10').
     Rather than crude arithmetic midpoints (which overestimate log-normal concentrations),
     the default "gold standard" fits an interval-censored lognormal maximum likelihood model
     (Tobit formulation) per determinand. It computes mathematically exact conditional
     expectations E[Y | a <= Y <= b] for closed ranges, left-censored detection limits (<0.2),
     and right-censored saturation levels (>10), along with analytical binning uncertainties
     sigma_bin = (log(b) - log(a)) / sqrt(12).
  4. Multiple Range Modes Supported:
     - `tobit`     : (Default) Interval-censored MLE with truncated lognormal expectation.
     - `geometric` : Geometric midpoint y* = sqrt(a * b) = exp((log a + log b) / 2).
     - `midpoint`  : Arithmetic midpoint y* = (a + b) / 2.
     - `lower`     : Lower bound a.
     - `upper`     : Upper bound b.
     - `skip`      : Skips range intervals.
  5. Wide & Long Format Support:
     Supports both single-determinand long-format tables and wide multi-determinand tables
     (e.g., separate columns for Nitrate, Phosphate, pH, Turbidity), unpivoting wide
     measurements into standardized long-format `samples.csv` rows.
  6. Built-in Coordinate Projection:
     Converts WGS84 decimal degrees (lat/lon) into British National Grid (EPSG:27700) meters
     in pure Python standard library (Airy 1830 ellipsoid) with zero external dependencies.
  7. Row Filtering & Subsetting:
     Supports targeted filtering on any column (e.g. `--filter="Country=England"`).
  8. Optional JSON Profile Registration:
     Can automatically generate a permanent `configs/sources/<name>.json` profile.

Zero External Dependencies:
  Runs on pure Python 3 standard library (csv, math, re, os, json, argparse, datetime).
"""

import argparse
import csv
import datetime
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict

# Base workspace directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from s3n_py_utils import (
    load_directory_config as _shared_load_directory_config,
    wgs84_to_bng,
    normalize_name,
    find_matching_column,
    normal_cdf,
    normal_pdf,
    parse_raw_measurement_token,
    IntervalTobitModel,
    compute_alternative_range_result,
    parse_datetime_iso,
    extract_unit_from_string,
)


def load_directory_config(base_dir=BASE_DIR):
    """Loads directory configuration from configs/directories.json, paths.json, or config.json."""
    cfg = _shared_load_directory_config(base_dir=base_dir, auto_create=False)
    return {"data_dir": cfg["data_dir"], "output_dir": cfg["output_dir"], "config_file": cfg["config_file"]}


# ==============================================================================
# 1. Chemical Parameter Dictionary & Unit Matching
# ==============================================================================


def match_known_determinand(col_name):
    """Matches a column name against known water quality parameters."""
    c_norm = normalize_name(col_name)
    for key, info in KNOWN_DETERMINANDS.items():
        if key in c_norm:
            unit = extract_unit_from_string(col_name) or info["unit"]
            return {"code": info["code"], "name": info["name"], "unit": unit}

    clean_name = re.sub(r"\s*[\(\[].*?[\)\]]", "", str(col_name)).strip()
    code = re.sub(r"[^a-zA-Z0-9_]", "_", clean_name.lower())[:16]
    unit = extract_unit_from_string(col_name) or "mg/L"
    return {"code": code, "name": clean_name, "unit": unit}


# ==============================================================================
# 9. Core Dataset Formatter Function
# ==============================================================================

def format_custom_dataset(
    input_path,
    output_dir=None,
    source_name=None,
    source_label=None,
    # Column overrides
    site_id_col=None,
    site_name_col=None,
    easting_col=None,
    northing_col=None,
    date_col=None,
    time_col=None,
    val_cols=None,
    det_code_col=None,
    det_name_col=None,
    unit_col=None,
    sign_col=None,
    # Determinands & units overrides
    default_det_code="0117",
    default_det_name="Nitrate as N",
    default_unit="mg/L",
    # Range conversion options
    ranges_to_midpoint=True,
    range_mode="tobit",  # "tobit" is the default gold-standard approach!
    # Coordinate transformation
    convert_to_bng=True,
    input_crs=None,
    # Grouping & filtering
    group_sites_by="auto",
    site_id_prefix=None,
    filter_expr=None,
    bbox=None,
    limit=None,
    generate_config=False,
    verbose=True,
):
    """
    Reads an external monitoring dataset, formats it, and outputs canonical
    `sites.csv` and `samples.csv` files with interval-censored Tobit expectations
    and binning uncertainties.
    """
    # 1. Resolve Input File
    if os.path.isdir(input_path):
        csvs = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.lower().endswith(".csv")]
        if not csvs:
            raise FileNotFoundError(f"No CSV files found in input directory: {input_path}")
        csvs.sort(key=lambda p: os.path.getsize(p), reverse=True)
        csv_file = csvs[0]
        if verbose:
            print(f"Discovered input CSV in directory: {csv_file}")
    elif os.path.isfile(input_path):
        csv_file = input_path
    else:
        raise FileNotFoundError(f"Input path not found: {input_path}")

    # 2. Resolve Source Name and Output Directory
    if not source_name:
        base_cand = os.path.basename(os.path.dirname(csv_file)) if os.path.dirname(csv_file) else "custom_source"
        if base_cand.lower() in ("data", "csv", "raw", "input", "."):
            base_cand = os.path.splitext(os.path.basename(csv_file))[0]
        source_name = re.sub(r"[^a-zA-Z0-9_]", "_", base_cand.lower()).strip("_")
    source_name = source_name.lower().strip()

    raw_source_label = source_label or source_name.replace("_", " ").title()
    source_label = raw_source_label

    dir_cfg = load_directory_config(BASE_DIR)
    data_root = dir_cfg["data_dir"]

    explicit_out_dir = bool(output_dir)
    if not output_dir:
        output_dir = os.path.join(data_root, "custom", source_name)
    os.makedirs(output_dir, exist_ok=True)

    # 3. Read Header and Determine Column Mapping
    with open(csv_file, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        try:
            fieldnames = next(reader)
        except StopIteration:
            raise ValueError(f"Input CSV {csv_file} is empty!")

    # Aliases
    id_aliases = [
        "site_id", "station_id", "monitoring_id", "location_id", "network_id",
        "sample_id", "site_code", "station_code", "loc_id", "stn_id",
        "site_name", "sitename", "site", "objectid", "globalid", "id"
    ]
    name_aliases = [
        "site_name", "station_name", "location_name", "monitoring_location_name",
        "site_description", "water_body_name", "waterbody_name", "river_name",
        "description", "name", "site"
    ]
    x_aliases = [
        "ngr_easting", "easting", "eastings", "bng_easting", "bng_x",
        "x_coord", "xcoord", "point_x", "longitude", "long", "lon", "x"
    ]
    y_aliases = [
        "ngr_northing", "northing", "northings", "bng_northing", "bng_y",
        "y_coord", "ycoord", "point_y", "latitude", "lat", "y"
    ]
    date_aliases = [
        "sample_datetime", "datetime", "date_time", "survey_datetime",
        "sample_date", "survey_date", "date", "timestamp", "collection_date"
    ]
    time_aliases = [
        "sample_time", "survey_time", "time_of_day", "time", "collection_time"
    ]
    sign_aliases = [
        "meas_sign", "less_than_flag", "censored", "is_censored",
        "qualifier", "flag", "sign", "detection_limit_flag"
    ]
    val_aliases = [
        "meas_result", "result", "value", "val", "reading", "score",
        "concentration", "measurement", "amount"
    ]

    mapped_id = site_id_col or find_matching_column(fieldnames, id_aliases)
    mapped_name = site_name_col or find_matching_column(fieldnames, name_aliases)
    mapped_x = easting_col or find_matching_column(fieldnames, x_aliases)
    mapped_y = northing_col or find_matching_column(fieldnames, y_aliases)
    mapped_date = date_col or find_matching_column(fieldnames, date_aliases)
    mapped_time = time_col or find_matching_column(fieldnames, time_aliases)
    mapped_sign = sign_col or find_matching_column(fieldnames, sign_aliases)
    mapped_det_code = det_code_col or find_matching_column(fieldnames, ["det_code", "determinand_code", "param_code", "code"])
    mapped_det_name = det_name_col or find_matching_column(fieldnames, ["det_name", "determinand_name", "det_desc", "parameter", "compound_name"])
    mapped_unit = unit_col or find_matching_column(fieldnames, ["unit_short_desc", "unit", "units", "uom"])

    if not mapped_x or not mapped_y:
        raise ValueError(
            f"Could not automatically detect Easting/Northing or Longitude/Latitude columns in {csv_file}.\n"
            f"Available columns: {fieldnames}\n"
            "Please specify explicitly with --x-col and --y-col."
        )

    val_col_list = []
    if val_cols:
        val_col_list = [c.strip() for c in val_cols.split(",") if c.strip()]
        for vc in val_col_list:
            if vc not in fieldnames:
                closest = find_matching_column(fieldnames, [vc])
                if closest:
                    val_col_list[val_col_list.index(vc)] = closest
                else:
                    raise ValueError(f"Specified value column '{vc}' not found in CSV columns: {fieldnames}")
    else:
        single_val = find_matching_column(fieldnames, val_aliases)
        if single_val:
            val_col_list = [single_val]
        else:
            for f_col in fieldnames:
                c_norm = normalize_name(f_col)
                c_tokens = set(re.split(r"[^a-z0-9]+", c_norm))
                if any(k in c_norm for k in ("nitrate", "phosphate", "ammon", "turbid", "conduct", "oxygen")) or (
                    "ph" in c_tokens or "temp" in c_tokens or "temperature" in c_tokens or "water_temp" in c_norm
                ):
                    val_col_list.append(f_col)

            if not val_col_list:
                raise ValueError(
                    f"Could not automatically identify measurement value column(s) in {csv_file}.\n"
                    f"Available columns: {fieldnames}\n"
                    "Please specify explicitly with --val-col."
                )

    parsed_filters = []
    if filter_expr:
        for f_part in filter_expr.split(","):
            if "=" in f_part:
                k, v = f_part.split("=", 1)
                k = k.strip()
                v = v.strip().strip("'\"")
                matched_k = find_matching_column(fieldnames, [k]) or k
                parsed_filters.append((matched_k, v.lower()))

    parsed_bbox = None
    if bbox:
        try:
            b_parts = [float(x.strip()) for x in bbox.split(",")]
            if len(b_parts) == 4:
                parsed_bbox = (b_parts[0], b_parts[1], b_parts[2], b_parts[3])
        except ValueError:
            print(f"Warning: Invalid bbox '{bbox}'. Expected 'min_x,min_y,max_x,max_y'. Ignoring.")

    effective_mode = range_mode if ranges_to_midpoint else "skip"

    if verbose:
        print("=" * 76)
        print(f"  Formatting Custom Monitoring Dataset: {source_label} ({source_name})")
        print("=" * 76)
        print(f"  Input File        : {csv_file}")
        print(f"  Output Directory  : {output_dir}")
        print("  Column Mappings   :")
        print(f"    - Site ID       : {mapped_id or '(auto-derived by location)'}")
        print(f"    - Site Name     : {mapped_name or '(none)'}")
        print(f"    - Coordinates   : X='{mapped_x}', Y='{mapped_y}'")
        print(f"    - Sample Date   : '{mapped_date}'" + (f", Time='{mapped_time}'" if mapped_time else ""))
        print(f"    - Measurement(s): {', '.join(val_col_list)}")
        print(f"    - Range Mode    : {effective_mode.upper()} " +
              ("(\"The Gold Standard\": Truncated Lognormal Tobit MLE)" if effective_mode == "tobit" else ""))
        if parsed_filters:
            print(f"    - Filter(s)     : {parsed_filters}")
        print("-" * 76)

    # 4. PASS 1: Read, Filter & Collect Measurement Tokens for Model Fitting
    # We aggregate distinct tokens by determinand to fit the Tobit MLE model instantaneously.
    det_token_counts = defaultdict(Counter)  # det_code -> Counter((type, a, b))
    det_info_map = {}                        # det_code -> {"code", "name", "unit"}

    # Temporary storage of parsed rows before expectation mapping
    staged_records = []
    sites_dict = {}
    coord_site_map = {}
    name_site_map = {}

    total_rows = 0
    filtered_rows = 0
    auto_site_counter = 0
    prefix = site_id_prefix or (source_name.upper()[:4] + "_")

    with open(csv_file, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_rows += 1
            if limit and total_rows > limit:
                break

            # Filter
            if parsed_filters:
                match_all = True
                for fk, fv in parsed_filters:
                    row_val = str(row.get(fk, "")).strip().lower()
                    if row_val != fv:
                        match_all = False
                        break
                if not match_all:
                    filtered_rows += 1
                    continue

            # Coordinates
            raw_x = row.get(mapped_x, "").strip()
            raw_y = row.get(mapped_y, "").strip()
            if not raw_x or not raw_y:
                continue
            try:
                px = float(raw_x)
                py = float(raw_y)
            except ValueError:
                continue

            if parsed_bbox:
                if not (parsed_bbox[0] <= px <= parsed_bbox[2] and parsed_bbox[1] <= py <= parsed_bbox[3]):
                    continue

            is_wgs84 = (abs(px) <= 180.0 and abs(py) <= 90.0) or (input_crs == 4326)
            if is_wgs84 and convert_to_bng:
                east, north = wgs84_to_bng(py, px)
            else:
                east, north = round(px, 2), round(py, 2)

            raw_site_name = row.get(mapped_name, "").strip() if mapped_name else ""
            raw_id_val = row.get(mapped_id, "").strip() if mapped_id else ""

            # Site ID assignment
            if group_sites_by == "coords":
                ckey = (round(px, 5), round(py, 5))
                if ckey not in coord_site_map:
                    auto_site_counter += 1
                    coord_site_map[ckey] = f"{prefix}{auto_site_counter:05d}"
                sid = coord_site_map[ckey]
            elif group_sites_by == "name" and raw_site_name:
                if raw_site_name not in name_site_map:
                    auto_site_counter += 1
                    name_site_map[raw_site_name] = f"{prefix}{auto_site_counter:05d}"
                sid = name_site_map[raw_site_name]
            elif group_sites_by == "name_and_coords":
                nckey = (raw_site_name, round(px, 5), round(py, 5))
                if nckey not in coord_site_map:
                    auto_site_counter += 1
                    coord_site_map[nckey] = f"{prefix}{auto_site_counter:05d}"
                sid = coord_site_map[nckey]
            elif raw_id_val:
                sid = raw_id_val
            elif raw_site_name:
                if raw_site_name not in name_site_map:
                    auto_site_counter += 1
                    name_site_map[raw_site_name] = f"{prefix}{auto_site_counter:05d}"
                sid = name_site_map[raw_site_name]
            else:
                ckey = (round(px, 5), round(py, 5))
                if ckey not in coord_site_map:
                    auto_site_counter += 1
                    coord_site_map[ckey] = f"{prefix}{auto_site_counter:05d}"
                sid = coord_site_map[ckey]

            if sid not in sites_dict:
                s_name = raw_site_name or sid
                sites_dict[sid] = {
                    "site_id": sid,
                    "easting": east,
                    "northing": north,
                    "site_name": s_name,
                    "source": source_label,
                }

            # Datetime
            raw_date = row.get(mapped_date, "").strip() if mapped_date else ""
            raw_time = row.get(mapped_time, "").strip() if mapped_time else ""
            dt_iso = parse_datetime_iso(raw_date, raw_time)

            # Measurement Values
            for v_col in val_col_list:
                raw_v = row.get(v_col, "").strip()
                if not raw_v:
                    continue

                if mapped_det_code and row.get(mapped_det_code):
                    d_code = row.get(mapped_det_code).strip()
                    d_name = row.get(mapped_det_name, "").strip() if mapped_det_name else default_det_name
                    d_unit = row.get(mapped_unit, "").strip() if mapped_unit else default_unit
                elif len(val_col_list) == 1 and default_det_code:
                    d_code = default_det_code
                    d_name = default_det_name
                    d_unit = default_unit
                else:
                    d_info = match_known_determinand(v_col)
                    d_code = d_info["code"]
                    d_name = d_info["name"]
                    d_unit = d_info["unit"]

                det_info_map[d_code] = {"code": d_code, "name": d_name, "unit": d_unit}

                token = parse_raw_measurement_token(raw_v)
                if token is None:
                    continue
                if token["type"] == "interval" and effective_mode == "skip":
                    continue

                # Explicit sign override from dataset column if available
                if mapped_sign and row.get(mapped_sign):
                    ex_sign = row.get(mapped_sign).strip()
                    if ex_sign in ("<", ">"):
                        token["sign"] = ex_sign
                        if ex_sign == "<":
                            token["type"] = "left"
                            token["b"] = token.get("val") or token.get("b")
                        else:
                            token["type"] = "right"
                            token["a"] = token.get("val") or token.get("a")

                token_key = (token["type"], token.get("a"), token.get("b"))
                det_token_counts[d_code][token_key] += 1

                staged_records.append({
                    "site_id": sid,
                    "sample_datetime": dt_iso,
                    "det_code": d_code,
                    "token": token,
                })

    # 5. Fit Tobit Model (or Alternative Models) per Determinand
    tobit_models = {}
    if effective_mode == "tobit":
        if verbose:
            print("Fitting Interval-Censored Tobit Maximum Likelihood Models:")
        for d_code, counts in det_token_counts.items():
            d_name = det_info_map[d_code]["name"]
            model = IntervalTobitModel(name=d_name)
            model.fit_from_token_counts(counts)
            tobit_models[d_code] = model
            if verbose:
                med_val = math.exp(model.mu)
                print(f"  [{d_code}] {d_name:<28} : mu = {model.mu:+.3f}, sigma = {model.sigma:.3f} "
                      f"(Est. Median = {med_val:.3g} {det_info_map[d_code]['unit']}) "
                      f"[Intervals: {model.n_intervals:,}, Left: {model.n_left:,}, Right: {model.n_right:,}, Exact: {model.n_points:,}]")
        if verbose:
            print("-" * 76)

    # 6. PASS 2: Map Tokens to Unbiased Values and Uncertainties
    samples_list = []
    ranges_converted = 0
    censored_count = 0

    for rec in staged_records:
        d_code = rec["det_code"]
        token = rec["token"]
        d_info = det_info_map[d_code]

        if effective_mode == "tobit":
            model = tobit_models.get(d_code)
            meas_res, meas_sign, meas_low, meas_upp, meas_unc, meas_meth = model.predict(
                rec_type=token["type"],
                a=token.get("a"),
                b=token.get("b"),
                val=token.get("val")
            )
        else:
            meas_res, meas_sign, meas_low, meas_upp, meas_unc, meas_meth = compute_alternative_range_result(
                token=token,
                range_mode=effective_mode
            )

        if token["type"] == "interval":
            ranges_converted += 1
        if meas_sign in ("<", ">"):
            censored_count += 1

        samples_list.append({
            "site_id": rec["site_id"],
            "sample_datetime": rec["sample_datetime"],
            "det_code": d_code,
            "det_name": d_info["name"],
            "meas_result": f"{meas_res:.5g}" if isinstance(meas_res, float) else str(meas_res),
            "meas_sign": meas_sign,
            "unit_short_desc": d_info["unit"],
            "meas_lower": f"{meas_low:.5g}" if isinstance(meas_low, float) else str(meas_low),
            "meas_upper": f"{meas_upp:.5g}" if isinstance(meas_upp, float) else str(meas_upp),
            "meas_uncertainty": f"{meas_unc:.4f}" if isinstance(meas_unc, float) else str(meas_unc),
            "meas_method": meas_meth,
        })

    # 7. Write Canonical Output Files (Tag output folder with range mode when ranges are used)
    ranges_used = (ranges_converted > 0 and effective_mode in ("tobit", "geometric", "midpoint", "lower", "upper"))
    tagged_source_name = source_name
    if ranges_used:
        mode_suffix = f"_{effective_mode}"
        if not source_name.endswith(mode_suffix):
            tagged_source_name = f"{source_name}{mode_suffix}"
        
        default_dir = os.path.join(data_root, "custom", source_name)
        if (not explicit_out_dir) or (output_dir == default_dir) or (os.path.basename(output_dir) == source_name):
            output_dir = os.path.join(data_root, "custom", tagged_source_name)
        
        if f"[{effective_mode.title()}]" not in source_label:
            source_label = f"{raw_source_label} [{effective_mode.title()}]"
        for s_data in sites_dict.values():
            s_data["source"] = source_label

    os.makedirs(output_dir, exist_ok=True)
    sites_path = os.path.join(output_dir, "sites.csv")
    samples_path = os.path.join(output_dir, "samples.csv")

    with open(sites_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["site_id", "easting", "northing", "site_name", "source"])
        writer.writeheader()
        for s_data in sites_dict.values():
            writer.writerow(s_data)

    with open(samples_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "site_id", "sample_datetime", "det_code", "det_name",
            "meas_result", "meas_sign", "unit_short_desc",
            "meas_lower", "meas_upper", "meas_uncertainty", "meas_method"
        ])
        writer.writeheader()
        for smp in samples_list:
            writer.writerow(smp)

    # 8. Optional: Generate Declarative JSON Configuration Profile
    cfg_path = None
    if generate_config:
        cfg_dir = os.path.join(BASE_DIR, "configs", "sources")
        os.makedirs(cfg_dir, exist_ok=True)
        cfg_path = os.path.join(cfg_dir, f"{tagged_source_name}.json")
        cfg_obj = {
            "name": tagged_source_name,
            "abbr": f"{source_name.upper()[:6]}_{effective_mode.upper()[:3]}" if ranges_used else source_name.upper()[:6],
            "label": source_label,
            "range_mode": effective_mode if ranges_used else None,
            "description": f"Custom monitoring dataset formatted by format_custom_data.py ({effective_mode.upper()} mode)",
            "data_dir": os.path.relpath(output_dir, BASE_DIR),
            "sites_file": "sites.csv",
            "id_col": "site_id",
            "easting_col": "easting",
            "northing_col": "northing",
            "crs_epsg": 27700 if convert_to_bng else (4326 if is_wgs84 else 27700),
            "datasets": {
                "wq": {
                    "file": "samples.csv",
                    "id_col": "site_id",
                    "datetime_col": "sample_datetime",
                    "det_col": "det_code",
                    "det_desc_col": "det_name",
                    "val_col": "meas_result",
                    "cens_col": "meas_sign",
                    "cens_flag": "<",
                    "unit_col": "unit_short_desc"
                }
            }
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(cfg_obj, f, indent=2)

        # Write alias config if tagged_source_name differs from base source_name
        if tagged_source_name != source_name:
            base_cfg_path = os.path.join(cfg_dir, f"{source_name}.json")
            base_cfg_obj = {
                "name": source_name,
                "abbr": source_name.upper()[:6],
                "label": raw_source_label,
                "range_mode": effective_mode,
                "alias_target": tagged_source_name,
                "description": f"Custom monitoring dataset {source_name} (defaulting to {tagged_source_name})",
                "data_dir": os.path.relpath(output_dir, BASE_DIR),
                "sites_file": "sites.csv",
                "id_col": "site_id",
                "easting_col": "easting",
                "northing_col": "northing",
                "crs_epsg": 27700 if convert_to_bng else (4326 if is_wgs84 else 27700),
                "datasets": {
                    "wq": {
                        "file": "samples.csv",
                        "id_col": "site_id",
                        "datetime_col": "sample_datetime",
                        "det_col": "det_code",
                        "det_desc_col": "det_name",
                        "val_col": "meas_result",
                        "cens_col": "meas_sign",
                        "cens_flag": "<",
                        "unit_col": "unit_short_desc"
                    }
                }
            }
            with open(base_cfg_path, "w", encoding="utf-8") as f:
                json.dump(base_cfg_obj, f, indent=2)

    # 9. Print Execution Summary
    if verbose:
        print("Formatting complete!")
        print(f"  Total Rows Read        : {total_rows:,}")
        if filtered_rows > 0:
            print(f"  Rows Filtered Out      : {filtered_rows:,}")
        print(f"  Unique Sites Created   : {len(sites_dict):,} -> {sites_path}")
        print(f"  Valid Samples Written  : {len(samples_list):,} -> {samples_path}")
        print(f"  Ranges Processed       : {ranges_converted:,} intervals converted ({effective_mode.upper()})")
        if ranges_used:
            print(f"  Output Folder Tagged   : {output_dir} [Mode: {effective_mode.upper()}]")
        print(f"  Censored Values (<, >) : {censored_count:,} flagged with bounds")
        if cfg_path:
            print(f"  JSON Profile Generated : {cfg_path}")
        print("-" * 76)
        print("Preview of sites.csv:")
        for s_sample in list(sites_dict.values())[:3]:
            print(f"  {s_sample}")
        print("Preview of samples.csv:")
        for sm_sample in samples_list[:3]:
            print(f"  {sm_sample}")
        print("=" * 76)
        print("Next Steps in Pipeline:")
        print(f"  1. Preprocess river reaches & snapping:")
        print(f"     python3 scripts/01_preprocess_orn.py --source={tagged_source_name}")
        print(f"  2. Inspect station snapping diagnostics:")
        print(f"     Rscript scripts/06_station_snapping_diagnostics.R --source={tagged_source_name} --list-catchments")
        print(f"  3. Fit S3N model & kriging predictions:")
        det_run = default_det_code if default_det_code else (samples_list[0]["det_code"] if samples_list else "0117")
        print(f"     Rscript scripts/run_s3n_determinand.R {det_run} --source={tagged_source_name}")
        print("=" * 76)

    return {
        "source_name": tagged_source_name,
        "base_source_name": source_name,
        "range_mode": effective_mode if ranges_used else None,
        "ranges_used": ranges_used,
        "sites_path": sites_path,
        "samples_path": samples_path,
        "config_path": cfg_path,
        "n_sites": len(sites_dict),
        "n_samples": len(samples_list),
        "ranges_converted": ranges_converted,
        "censored_count": censored_count,
        "tobit_models": tobit_models,
    }


# ==============================================================================
# 10. CLI Entrypoint
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Format any custom water quality dataset to work with the S3N stream network pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 1. Gold standard interval-censored Tobit model (Default):
  python3 scripts/format_custom_data.py data/Waterblitz \\
    --source-name=waterblitz \\
    --filter="Country=England" \\
    --val-col="Nitrate (mg/L),Phosphate (mg/L)" \\
    --group-sites-by=coords \\
    --site-id-prefix="WB_" \\
    --generate-config

  # 2. Geometric midpoint mode:
  python3 scripts/format_custom_data.py data/Waterblitz \\
    --source-name=waterblitz \\
    --range-mode=geometric

  # 3. Arithmetic midpoint mode (e.g. 5-10 -> 7.5):
  python3 scripts/format_custom_data.py data/Waterblitz \\
    --source-name=waterblitz \\
    --range-mode=midpoint

  # 4. Explicit column mappings for an arbitrary sensor or citizen survey file:
  python3 scripts/format_custom_data.py path/to/my_readings.csv \\
    --source-name=my_survey \\
    --id-col="Station_Number" \\
    --name-col="River_Site" \\
    --x-col="East" \\
    --y-col="North" \\
    --date-col="Survey_Date" \\
    --val-col="Concentration" \\
    --det-code=0117 \\
    --generate-config
        """
    )

    parser.add_argument("input", nargs="?", help="Path to input CSV file or folder containing custom monitoring data")
    parser.add_argument("--input", "-i", dest="input_flag", help="Alternative flag for input CSV file or folder")
    parser.add_argument("--source-name", "-s", default=None, help="Identifier name for this source (e.g. waterblitz, riverfly)")
    parser.add_argument("--source-label", default=None, help="Human-readable label for this source")
    parser.add_argument("--out-dir", "--out", "-o", default=None, help="Output directory (defaults to data/custom/<source-name>/)")

    # Column Mapping Overrides
    parser.add_argument("--id-col", "--site-id-col", default=None, help="Column name for station/site identifier")
    parser.add_argument("--name-col", "--site-name-col", default=None, help="Column name for human-readable site name")
    parser.add_argument("--x-col", "--easting-col", "--lon-col", default=None, help="Column name for Easting (m) or Longitude (deg)")
    parser.add_argument("--y-col", "--northing-col", "--lat-col", default=None, help="Column name for Northing (m) or Latitude (deg)")
    parser.add_argument("--date-col", "--datetime-col", default=None, help="Column name for sample date or datetime")
    parser.add_argument("--time-col", default=None, help="Column name for sample time (if separate from date)")
    parser.add_argument("--val-col", "--value-col", "--result-col", default=None, help="Measurement column(s) (comma-separated for wide format)")
    parser.add_argument("--det-col", "--det-code-col", default=None, help="Column name for determinand code (if long format)")
    parser.add_argument("--det-name-col", default=None, help="Column name for determinand name/description")
    parser.add_argument("--unit-col", default=None, help="Column name for measurement units")
    parser.add_argument("--sign-col", "--cens-col", default=None, help="Column name for censoring flag/sign (<, =, >)")

    # Determinand Constants
    parser.add_argument("--det-code", default="0117", help="Default determinand code (default: 0117 for Nitrate)")
    parser.add_argument("--det-name", default="Nitrate as N", help="Default determinand description (default: 'Nitrate as N')")
    parser.add_argument("--unit", default="mg/L", help="Default measurement unit (default: 'mg/L')")

    # Range Conversion Options
    parser.add_argument(
        "--range-mode",
        choices=["tobit", "geometric", "midpoint", "lower", "upper", "skip"],
        default="tobit",
        help="Range handling strategy: 'tobit' (default, Gold Standard interval-censored MLE), 'geometric' (log-midpoint sqrt(a*b)), 'midpoint' (arithmetic (a+b)/2), 'lower', 'upper', or 'skip'"
    )
    parser.add_argument("--ranges-to-midpoint", action="store_true", default=True, help="Enable range conversion (default: True)")
    parser.add_argument("--no-ranges-to-midpoint", action="store_false", dest="ranges_to_midpoint", help="Do not convert ranges (skip interval records)")

    # Coordinate Projection
    parser.add_argument("--convert-to-bng", action="store_true", default=True, help="Convert WGS84 lat/lon to British National Grid meters (default: True)")
    parser.add_argument("--keep-coords", "--no-convert-bng", action="store_false", dest="convert_to_bng", help="Keep coordinates as raw values without projection")
    parser.add_argument("--crs", type=int, default=None, help="Explicit input CRS (e.g. 4326 for WGS84, 27700 for BNG)")

    # Grouping & Filtering
    parser.add_argument("--group-sites-by", choices=["auto", "id", "name", "coords", "name_and_coords"], default="auto", help="Site grouping strategy (default: auto)")
    parser.add_argument("--site-id-prefix", default=None, help="Prefix for auto-generated site IDs (e.g. 'WB_' or 'SITE_')")
    parser.add_argument("--filter", dest="filter_expr", default=None, help="Row filter expression (e.g. 'Country=England' or 'Country=England,Freshwater body type=River')")
    parser.add_argument("--bbox", default=None, help="Bounding box filter 'min_x,min_y,max_x,max_y'")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N rows (for testing)")
    parser.add_argument("--generate-config", action="store_true", default=False, help="Generate declarative JSON profile in configs/sources/<name>.json")

    args = parser.parse_args()

    input_target = args.input or args.input_flag
    if not input_target:
        parser.print_help()
        sys.exit(1)

    format_custom_dataset(
        input_path=input_target,
        output_dir=args.out_dir,
        source_name=args.source_name,
        source_label=args.source_label,
        site_id_col=args.id_col,
        site_name_col=args.name_col,
        easting_col=args.x_col,
        northing_col=args.y_col,
        date_col=args.date_col,
        time_col=args.time_col,
        val_cols=args.val_col,
        det_code_col=args.det_col,
        det_name_col=args.det_name_col,
        unit_col=args.unit_col,
        sign_col=args.sign_col,
        default_det_code=args.det_code,
        default_det_name=args.det_name,
        default_unit=args.unit,
        ranges_to_midpoint=args.ranges_to_midpoint,
        range_mode=args.range_mode,
        convert_to_bng=args.convert_to_bng,
        input_crs=args.crs,
        group_sites_by=args.group_sites_by,
        site_id_prefix=args.site_id_prefix,
        filter_expr=args.filter_expr,
        bbox=args.bbox,
        limit=args.limit,
        generate_config=args.generate_config,
        verbose=True,
    )


if __name__ == "__main__":
    main()
