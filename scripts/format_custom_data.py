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
    load_directory_config,
    load_determinand_catalog,
    wgs84_to_bng,
    normalize_name,
    find_matching_column,
    parse_raw_measurement_token,
    IntervalTobitModel,
    compute_alternative_range_result,
    parse_datetime_iso,
    extract_unit_from_string,
)


# ==============================================================================
# 1. Chemical Parameter Dictionary & Unit Matching
# ==============================================================================

_DET_STEM_TO_CODE = (
    ("orthophosphate", "0180"),
    ("phosphate", "0180"),
    ("nitrate", "0117"),
    ("nitrite", "0118"),
    ("ammon", "0111"),
    ("oxidised", "0116"),
    ("organiccarbon", "0301"),
    ("conductiv", "0077"),
    ("alkalin", "0162"),
    ("dissolvedoxygen", "9924"),
    ("watertemp", "0076"),
    ("temperature", "0076"),
)
_DET_CATALOG_CACHE = None


def match_known_determinand(col_name):
    """Matches a column name against known water quality parameters in configs/determinands.json."""
    global _DET_CATALOG_CACHE
    if _DET_CATALOG_CACHE is None:
        _DET_CATALOG_CACHE = load_determinand_catalog(BASE_DIR)

    c_norm = normalize_name(col_name)
    c_tokens = set(re.findall(r"[a-z0-9]+", str(col_name).lower()))
    matched_code = None
    for stem, code in _DET_STEM_TO_CODE:
        if stem in c_norm:
            matched_code = code
            break
    if matched_code is None:
        for code, info in _DET_CATALOG_CACHE.items():
            short_norm = normalize_name(info.get("short", ""))
            name_norm = normalize_name(info.get("name", ""))
            if c_norm == code:
                matched_code = code
                break
            if short_norm and len(short_norm) >= 4 and short_norm in c_norm:
                matched_code = code
                break
            if name_norm and len(name_norm) >= 4 and name_norm in c_norm:
                matched_code = code
                break
    if matched_code is None and "ph" in c_tokens:
        matched_code = "0061"

    if matched_code and matched_code in _DET_CATALOG_CACHE:
        info = _DET_CATALOG_CACHE[matched_code]
        unit = extract_unit_from_string(col_name) or info.get("unit", "mg/L")
        return {"code": matched_code, "name": info.get("name", matched_code), "unit": unit, "matched": True}
    if matched_code:
        unit = extract_unit_from_string(col_name) or "mg/L"
        return {"code": matched_code, "name": str(col_name).strip(), "unit": unit, "matched": True}

    clean_name = re.sub(r"\s*[\(\[].*?[\)\]]", "", str(col_name)).strip()
    code = re.sub(r"[^a-zA-Z0-9_]", "_", clean_name.lower())[:16]
    unit = extract_unit_from_string(col_name) or "mg/L"
    return {"code": code, "name": clean_name, "unit": unit, "matched": False}


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
                c_tokens = set(re.findall(r"[a-z0-9]+", f_col.lower()))
                if any(k in c_norm for k in ("nitrate", "phosphate", "ammon", "turbid", "conduct", "oxygen", "watertemp", "organiccarbon")) or (
                    "ph" in c_tokens or "temp" in c_tokens or "temperature" in c_tokens
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

    effective_mode = range_mode

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
    saw_wgs84 = False
    prefix = site_id_prefix or (source_name.upper()[:4] + "_")

    bng_cache = {}
    det_match_cache = {}

    with open(csv_file, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader, None) or []
        col_idx = {h: i for i, h in enumerate(header)}
        n_cols = len(header)

        filter_idx = [(col_idx[fk], fv) for fk, fv in parsed_filters if fk in col_idx] if parsed_filters else []
        missing_filter = bool(parsed_filters) and (len(filter_idx) != len(parsed_filters))
        idx_x = col_idx[mapped_x]
        idx_y = col_idx[mapped_y]
        idx_name = col_idx.get(mapped_name) if mapped_name else None
        idx_id = col_idx.get(mapped_id) if mapped_id else None
        idx_date = col_idx.get(mapped_date) if mapped_date else None
        idx_time = col_idx.get(mapped_time) if mapped_time else None
        idx_det_code = col_idx.get(mapped_det_code) if mapped_det_code else None
        idx_det_name = col_idx.get(mapped_det_name) if mapped_det_name else None
        idx_unit = col_idx.get(mapped_unit) if mapped_unit else None
        idx_sign = col_idx.get(mapped_sign) if mapped_sign else None
        val_col_pairs = [(v_col, col_idx[v_col]) for v_col in val_col_list if v_col in col_idx]

        # Pre-resolve wide-table determinand metadata outside the per-row loop
        wide_det_meta = {}
        for v_col, _ in val_col_pairs:
            d_info = match_known_determinand(v_col)
            if d_info.get("matched") or not (len(val_col_list) == 1 and default_det_code):
                wide_det_meta[v_col] = (d_info["code"], d_info["name"], d_info["unit"])
            else:
                wide_det_meta[v_col] = (default_det_code, default_det_name, default_unit)

        for row in reader:
            if not row:
                continue
            total_rows += 1
            if limit and total_rows > limit:
                break
            if len(row) < n_cols:
                row = row + [""] * (n_cols - len(row))

            # Filter
            if missing_filter:
                filtered_rows += 1
                continue
            if filter_idx:
                match_all = True
                for f_i, fv in filter_idx:
                    if row[f_i].strip().lower() != fv:
                        match_all = False
                        break
                if not match_all:
                    filtered_rows += 1
                    continue

            # Coordinates
            raw_x = row[idx_x].strip()
            raw_y = row[idx_y].strip()
            if not raw_x or not raw_y:
                continue
            try:
                px = float(raw_x)
                py = float(raw_y)
            except ValueError:
                continue

            if input_crs == 27700:
                is_wgs84 = False
            elif input_crs == 4326:
                is_wgs84 = True
            else:
                is_wgs84 = (abs(px) <= 180.0 and abs(py) <= 90.0)
            if is_wgs84:
                saw_wgs84 = True
            if is_wgs84 and convert_to_bng:
                xy_key = (px, py)
                bng_res = bng_cache.get(xy_key)
                if bng_res is None:
                    bng_res = wgs84_to_bng(py, px)
                    bng_cache[xy_key] = bng_res
                east, north = bng_res
            else:
                east, north = round(px, 2), round(py, 2)

            if parsed_bbox:
                is_bbox_bng = any(abs(v) > 360.0 for v in parsed_bbox)
                bx, by = (east, north) if is_bbox_bng else (px, py)
                if not (parsed_bbox[0] <= bx <= parsed_bbox[2] and parsed_bbox[1] <= by <= parsed_bbox[3]):
                    continue

            raw_site_name = row[idx_name].strip() if idx_name is not None else ""
            raw_id_val = row[idx_id].strip() if idx_id is not None else ""

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
            raw_date = row[idx_date].strip() if idx_date is not None else ""
            raw_time = row[idx_time].strip() if idx_time is not None else ""
            dt_iso = parse_datetime_iso(raw_date, raw_time)

            row_det_code = row[idx_det_code].strip() if idx_det_code is not None else ""
            row_det_name = row[idx_det_name].strip() if idx_det_name is not None else ""
            row_unit = row[idx_unit].strip() if idx_unit is not None else ""
            row_sign = row[idx_sign].strip().lower() if idx_sign is not None else ""

            # Measurement Values
            for v_col, v_idx in val_col_pairs:
                raw_v = row[v_idx].strip()
                if not raw_v:
                    continue

                if row_det_code:
                    d_code = row_det_code
                    d_name = row_det_name or default_det_name
                    d_unit = row_unit or default_unit
                elif row_det_name:
                    d_info = det_match_cache.get(row_det_name)
                    if d_info is None:
                        d_info = match_known_determinand(row_det_name)
                        det_match_cache[row_det_name] = d_info
                    d_code = d_info["code"]
                    d_name = d_info["name"] if d_info.get("matched") else row_det_name
                    d_unit = row_unit or d_info["unit"]
                else:
                    d_code, d_name, w_unit = wide_det_meta[v_col]
                    d_unit = row_unit or w_unit

                if d_code not in det_info_map:
                    det_info_map[d_code] = {"code": d_code, "name": d_name, "unit": d_unit}
                else:
                    info_entry = det_info_map[d_code]
                    if d_name and not info_entry.get("name"):
                        info_entry["name"] = d_name
                    if d_unit and not info_entry.get("unit"):
                        info_entry["unit"] = d_unit

                token = parse_raw_measurement_token(raw_v)
                if token is None:
                    continue
                tok_type = token["type"]
                if tok_type == "interval" and effective_mode == "skip":
                    continue

                tok_a = token.get("a")
                tok_b = token.get("b")
                tok_val = token.get("val")

                # Explicit sign override from dataset column if available
                if row_sign:
                    ex_sign = "<" if row_sign in ("<", "yes", "true", "1", "y", "lt") else (">" if row_sign in (">", "gt") else row_sign)
                    if ex_sign in ("<", ">"):
                        token["sign"] = ex_sign
                        bound_val = tok_val if tok_val is not None else (tok_b if ex_sign == "<" else tok_a)
                        if ex_sign == "<":
                            tok_type = "left"
                            tok_a = 0.0
                            tok_b = bound_val
                        else:
                            tok_type = "right"
                            tok_a = bound_val
                            tok_b = None
                        token["type"] = tok_type
                        token["a"] = tok_a
                        token["b"] = tok_b

                token_key = (tok_type, tok_a, tok_b)
                det_token_counts[d_code][token_key] += 1

                staged_records.append((sid, dt_iso, d_code, token))

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
    samples_rows = []
    ranges_converted = 0
    censored_count = 0
    pred_cache = {}

    for sid, dt_iso, d_code, token in staged_records:
        d_info = det_info_map[d_code]
        tok_type = token["type"]
        tok_a = token.get("a")
        tok_b = token.get("b")
        tok_val = token.get("val")
        cache_k = (d_code, tok_type, tok_a, tok_b, tok_val, token.get("sign"))
        cached_fmt = pred_cache.get(cache_k)

        if cached_fmt is None:
            if effective_mode == "tobit":
                model = tobit_models.get(d_code)
                meas_res, meas_sign, meas_low, meas_upp, meas_unc, meas_meth = model.predict(
                    rec_type=tok_type, a=tok_a, b=tok_b, val=tok_val
                )
            else:
                meas_res, meas_sign, meas_low, meas_upp, meas_unc, meas_meth = compute_alternative_range_result(
                    token=token, range_mode=effective_mode
                )
            cached_fmt = (
                f"{meas_res:.5g}" if isinstance(meas_res, float) else ("" if meas_res is None else str(meas_res)),
                meas_sign,
                f"{meas_low:.5g}" if isinstance(meas_low, float) else ("" if meas_low is None else str(meas_low)),
                f"{meas_upp:.5g}" if isinstance(meas_upp, float) else ("" if meas_upp is None else str(meas_upp)),
                f"{meas_unc:.4f}" if isinstance(meas_unc, float) else ("" if meas_unc is None else str(meas_unc)),
                meas_meth,
            )
            pred_cache[cache_k] = cached_fmt

        res_s, meas_sign, low_s, upp_s, unc_s, meas_meth = cached_fmt
        if tok_type == "interval":
            ranges_converted += 1
        if meas_sign in ("<", ">"):
            censored_count += 1

        samples_rows.append([
            sid, dt_iso, d_code, d_info["name"],
            res_s, meas_sign, d_info["unit"],
            low_s, upp_s, unc_s, meas_meth,
        ])

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
        writer.writerows(sites_dict.values())

    with open(samples_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "site_id", "sample_datetime", "det_code", "det_name",
            "meas_result", "meas_sign", "unit_short_desc",
            "meas_lower", "meas_upper", "meas_uncertainty", "meas_method"
        ])
        writer.writerows(samples_rows)

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
            "crs_epsg": 27700 if (saw_wgs84 and convert_to_bng) else (4326 if saw_wgs84 else (input_crs or 27700)),
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
                    "unit_col": "unit_short_desc",
                    "unc_col": "meas_uncertainty",
                    "meth_col": "meas_method"
                }
            }
        }
        with open(cfg_path, "w", encoding="utf-8") as f:
            json.dump(cfg_obj, f, indent=2)

    # 9. Print Execution Summary
    if verbose:
        print("Formatting complete!")
        print(f"  Total Rows Read        : {total_rows:,}")
        if filtered_rows > 0:
            print(f"  Rows Filtered Out      : {filtered_rows:,}")
        print(f"  Unique Sites Created   : {len(sites_dict):,} -> {sites_path}")
        print(f"  Valid Samples Written  : {len(samples_rows):,} -> {samples_path}")
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
        for sm_sample in samples_rows[:3]:
            print(f"  {sm_sample}")
        print("=" * 76)
        print("Next Steps in Pipeline:")
        print(f"  1. Preprocess river reaches & snapping:")
        print(f"     python3 scripts/01_preprocess_orn.py --source={tagged_source_name}")
        print(f"  2. Inspect station snapping diagnostics:")
        print(f"     Rscript scripts/06_station_snapping_diagnostics.R --source={tagged_source_name} --list-catchments")
        print(f"  3. Fit S3N model & kriging predictions:")
        det_run = default_det_code if default_det_code else (samples_rows[0][2] if samples_rows else "0117")
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
        "n_samples": len(samples_rows),
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
    parser.add_argument("--out-dir", "-o", default=None, help="Output directory (defaults to data/custom/<source-name>/)")

    # Column Mapping Overrides
    parser.add_argument("--id-col", default=None, help="Column name for station/site identifier")
    parser.add_argument("--name-col", default=None, help="Column name for human-readable site name")
    parser.add_argument("--x-col", "--lon-col", default=None, help="Column name for Easting (m) or Longitude (deg)")
    parser.add_argument("--y-col", "--lat-col", default=None, help="Column name for Northing (m) or Latitude (deg)")
    parser.add_argument("--date-col", default=None, help="Column name for sample date or datetime")
    parser.add_argument("--time-col", default=None, help="Column name for sample time (if separate from date)")
    parser.add_argument("--val-col", default=None, help="Measurement column(s) (comma-separated for wide format)")
    parser.add_argument("--det-col", default=None, help="Column name for determinand code (if long format)")
    parser.add_argument("--det-name-col", default=None, help="Column name for determinand name/description")
    parser.add_argument("--unit-col", default=None, help="Column name for measurement units")
    parser.add_argument("--sign-col", default=None, help="Column name for censoring flag/sign (<, =, >)")

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

    # Coordinate Projection
    parser.add_argument("--convert-to-bng", action="store_true", default=True, help="Convert WGS84 lat/lon to British National Grid meters (default: True)")
    parser.add_argument("--no-convert-bng", action="store_false", dest="convert_to_bng", help="Keep coordinates as raw values without projection")
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
