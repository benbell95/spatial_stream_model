#!/usr/bin/env python3
"""
01_preprocess_orn.py

Network-agnostic stream network and water quality preprocessor for Scalable
Spatial Stream Network (S3N) modeling.

Supports multiple river network datasets via JSON configuration profiles in
`configs/networks/` (or custom `--network-config` / `--network` GeoPackage paths):
  - `orn`            : Ordnance Survey Open Rivers Network (ORN v2, EPSG:27700)
  - `euhydro_thames` : EU-Hydro v1.3 Thames River Network Package (EPSG:3035)
  - `euhydro_tweed`  : EU-Hydro v1.3 Tweed River Network Package (EPSG:3035)
  - `euhydro`        : EU-Hydro v1.3 Combined Thames + Tweed Network (EPSG:3035)

Key Future-Proofing Capabilities:
  1. Universal GeoPackage WKB geometry parser supporting 2D (XY), 3D (XYZ/XYM),
     and 4D (XYZM) LineString & MultiLineString in ISO WKB and OGC EWKB, with
     automatic 3D vertex Z elevation extraction.
  2. Automatic CRS detection (`gpkg_contents.srs_id`) and batch planar coordinate
     transformation (`EPSG:3035` / any EPSG <-> `EPSG:27700` British National Grid)
     via `cs2cs`.
  3. Three-mode directed river tree topology builder:
       - `nodes`    : `(Fnode, Tnode)` node-pair matching (e.g. ORN, NHDPlus)
       - `pointer`  : `NEXTDOWNID` downstream reach pointer with `(Fnode, Tnode)` fallback (e.g. EU-Hydro, HydroRIVERS)
       - `geometry` : Endpoint coordinate snapping within tolerance
     Includes automatic cycle & braided-divergence breaking so every reach has at
     most one downstream `parent_dn` (strict directed tree required by S3N Tail-Up)
     plus automatic reach vertex orientation checking (`upstream -> downstream`).
  4. Automatic tree-based computation of missing hydrological attributes
     (`CatchID`, `Dist2Mth_km`, `US_Accum_km`, `Shreve`, `Strahler`, `Sinuosity`,
     `Elevation`, `SrcElevation`, `Slope`).
  5. Multi-tier covariate extraction (reach attributes, node table lookups, 3D
     geometry Z/slope/sinuosity, and fast spatial-hash BGS geology join).

Outputs (written to both `data/processed/<network_abbr>/` and `data/processed/`):
  - `england_rsn_reaches.csv`   : Reach-level topology, distances, AFV, covariates, and observations
  - `rsn_sites.csv`             : Station-level snapping summary
  - `rsn_reach_periods.csv`     : Multi-period reach summaries (all, year, season, month, year_season)
  - `rsn_catchments_lookup.csv` : Catchment-level lookup table with main river and RBD names
  - `network_metadata.json`     : Active water network metadata for downstream modeling & plots
"""

import argparse
import csv
import ctypes
import ctypes.util
import glob
import json
import math
import os
import re
import sqlite3
import statistics
import struct
import sys
import time
import zipfile
from collections import Counter, defaultdict

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

try:
    BASE_DIR = "." if os.path.samefile(os.getcwd(), _REPO_ROOT) else os.path.relpath(_REPO_ROOT, os.getcwd())
except OSError:
    BASE_DIR = "."

from s3n_py_utils import (
    DEFAULT_BASE_DIR,
    MINIMAL_DEFAULT_NETWORKS,
    MONTH_TO_SEASON,
    MONTH_NAME_MAP,
    RuntimeConfig,
    compute_afv_reach_weight,
    _rel_join,
    apply_forward_transform,
    apply_ros_to_records,
    build_runtime_config,
    compute_polyline_metrics,
    compute_ros_plotting_positions,
    detect_range_mode as _shared_detect_range_mode,
    ensure_directory_config,
    ensure_network_configs,
    ensure_source_configs,
    ensure_determinands_config,
    estimate_boxcox_lambda,
    estimate_yeojohnson_lambda,
    find_col,
    fit_ros_ols,
    impute_censored_by_site,
    load_directory_config,
    load_source_registry as _shared_load_source_registry,
    parse_datetime_iso,
    parse_gpkg_envelope,
    parse_gpkg_geometry,
    parse_gpkg_multilinestring,
    parse_month_arg,
    point_to_polyline,
    print_available_transforms,
    resolve_source_tag as _shared_resolve_source_tag,
    transform_coords_batch,
    wgs84_to_bng,
)


def ensure_all_configs(base_dir=BASE_DIR, data_dir=None, update_determinands=False):
    """Ensures directory, network, source, and determinand configs all exist."""
    ensure_directory_config(base_dir)
    ensure_network_configs(base_dir)
    ensure_source_configs(base_dir)
    ensure_determinands_config(base_dir, data_dir=data_dir, force_update=update_determinands)


DIR_CONFIG = load_directory_config(BASE_DIR, auto_create=False)
RUNTIME_CFG = build_runtime_config(BASE_DIR)
GIS_DIR = RUNTIME_CFG.gis_dir
DATA_DIR = RUNTIME_CFG.data_dir
RSN_DIR = RUNTIME_CFG.rsn_dir
SSN_DIR = RUNTIME_CFG.ssn_dir
CONFIG_DIR = RUNTIME_CFG.config_dir
SOURCE_CONFIG_DIR = RUNTIME_CFG.source_config_dir
OUT_DIR = RUNTIME_CFG.processed_dir


def load_source_registry(base_dir=BASE_DIR, runtime_cfg=None, data_dir=None):
    eff_data = data_dir or (runtime_cfg.data_dir if runtime_cfg is not None else DATA_DIR)
    return _shared_load_source_registry(base_dir=base_dir, data_dir=eff_data)


def resolve_source_tag(source="rsn", registry=None, runtime_cfg=None, data_dir=None):
    eff_data = data_dir or (runtime_cfg.data_dir if runtime_cfg is not None else DATA_DIR)
    return _shared_resolve_source_tag(source=source, registry=registry, data_dir=eff_data, base_dir=BASE_DIR)


def detect_range_mode(source_key, sources_reg, custom_samples=None, custom_dir=None, cli_mode="auto", runtime_cfg=None, data_dir=None):
    eff_data = data_dir or (runtime_cfg.data_dir if runtime_cfg is not None else DATA_DIR)
    return _shared_detect_range_mode(
        source_key=source_key,
        sources_reg=sources_reg,
        custom_samples=custom_samples,
        custom_dir=custom_dir,
        cli_mode=cli_mode,
        base_dir=BASE_DIR,
        data_dir=eff_data,
    )


def print_available_sources(base_dir=BASE_DIR, runtime_cfg=None):
    sources_reg, composites_reg = load_source_registry(base_dir=base_dir, runtime_cfg=runtime_cfg)
    print("Available Monitoring Data Sources:\n")
    print(f"  {'Source':<18} | {'Abbr':<6} | {'Description':<44} | {'Type':<12}")
    print("-" * 88)
    print(f"  {'rsn':<18} | {'RSN':<6} | {'River Surveillance Network (1,382 stations)':<44} | {'YES (Default)':<12}")
    print(f"  {'ssn':<18} | {'SSN':<6} | {'Small Streams Network (313 headwaters)':<44} | {'Registered':<12}")
    print(f"  {'rss':<18} | {'RSS':<6} | {'Combined River Surveillance & Small Streams':<44} | {'Composite':<12}")

    for s_name, s_info in sources_reg.items():
        if s_name in ("rsn", "ssn"):
            continue
        desc = s_info.get("description") or s_info.get("label", s_name)
        abbr = s_info.get("abbr", s_name.upper())
        print(f"  {s_name[:18]:<18} | {abbr[:6]:<6} | {desc[:44]:<44} | {'Configured':<12}")

    print("\nUsage & Multi-Source Syntax:")
    print("  python3 scripts/01_preprocess_orn.py --source=rsn           # River Surveillance Network (default)")
    print("  python3 scripts/01_preprocess_orn.py --source=ssn           # Small Streams Network only")
    print("  python3 scripts/01_preprocess_orn.py --source=rss           # Combined RSN + SSN")
    print("  python3 scripts/01_preprocess_orn.py --source=rsn+ssn       # On-the-fly combination syntax (+ or ,)")
    print("  python3 scripts/01_preprocess_orn.py --source=custom_example # Registered custom dataset")
    print("  python3 scripts/01_preprocess_orn.py --custom-sites=path/to/sites.csv --custom-samples=path/to/samples.csv\n")


def configure_runtime_dirs(gis_dir=None, data_dir=None, out_dir=None):
    """Builds RuntimeConfig and refreshes active module directory paths."""
    global GIS_DIR, DATA_DIR, RSN_DIR, SSN_DIR, OUT_DIR, DIR_CONFIG, RUNTIME_CFG
    DIR_CONFIG = load_directory_config(BASE_DIR, auto_create=False)
    RUNTIME_CFG = build_runtime_config(
        base_dir=BASE_DIR, gis_dir=gis_dir, data_dir=data_dir, out_dir=out_dir, create_out_dir=True
    )
    GIS_DIR = RUNTIME_CFG.gis_dir
    DATA_DIR = RUNTIME_CFG.data_dir
    RSN_DIR = RUNTIME_CFG.rsn_dir
    SSN_DIR = RUNTIME_CFG.ssn_dir
    OUT_DIR = RUNTIME_CFG.processed_dir
    return RUNTIME_CFG


def _opt_float(v):
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _extract_year_month(raw_dt):
    dt_iso = parse_datetime_iso(raw_dt) if raw_dt else ""
    if len(dt_iso) >= 7 and dt_iso[0:4].isdigit() and dt_iso[5:7].isdigit():
        return int(dt_iso[0:4]), int(dt_iso[5:7])
    return 2023, 6


def load_monitoring_stations_and_records(
    data_dir,
    source="rsn",
    dataset="wq",
    det_code="0117",
    custom_sites=None,
    custom_samples=None,
    custom_dir=None,
    custom_name="custom",
    custom_crs=27700
):
    """
    Loads station coordinates (pts), station-to-source mapping (pts_source),
    core water quality sites for catchment discovery (core_wq_sites),
    and observational sample records (all_records) across RSN, SSN, registered sources,
    or custom user datasets.
    """
    sources_reg, composites_reg = load_source_registry(base_dir=BASE_DIR, data_dir=data_dir)
    source_key, constituents, source_abbr, source_label = resolve_source_tag(
        source, (sources_reg, composites_reg), data_dir=data_dir
    )

    rsn_dir = _rel_join(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
    ssn_dir = _rel_join(data_dir, "NCEA_Small_Streams_Network_ARD_Open")

    pts = {}
    pts_source = {}
    core_wq_sites = set()
    all_records = []

    for src_item in constituents:
        src_label = src_item.upper()
        # 1. Custom CLI override (scoped to custom tokens, never hijacking RSN/SSN/registered sources)
        is_cli_custom = src_item in ("custom", custom_name.lower()) or (
            src_item not in ("rsn", "ssn") and src_item not in sources_reg and bool(custom_sites or custom_dir)
        )
        if is_cli_custom:
            sites_p = custom_sites or (os.path.join(custom_dir, "sites.csv") if custom_dir else None)
            samples_p = custom_samples or (os.path.join(custom_dir, "samples.csv") if custom_dir else None)
            if not sites_p or not os.path.exists(sites_p):
                cand_cs = os.path.join(data_dir, "custom", src_item, "sites.csv")
                if os.path.exists(cand_cs):
                    sites_p = cand_cs
                    samples_p = samples_p or os.path.join(data_dir, "custom", src_item, "samples.csv")
                else:
                    raise FileNotFoundError(f"Custom sites file not found for source '{src_item}'. Specify --custom-sites=path/to/sites.csv")

            print(f"   Loading custom stations from: {sites_p}")
            with open(sites_p, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                fnames = reader.fieldnames or []
                id_col = find_col(fnames, ["site_id", "network_id", "station_id", "id", "location_id", "site", "station", "code"])
                x_col = find_col(fnames, ["ngr_easting", "easting", "eastings", "bng_easting", "x", "point_x", "bng_x", "longitude", "lon"])
                y_col = find_col(fnames, ["ngr_northing", "northing", "northings", "bng_northing", "y", "point_y", "bng_y", "latitude", "lat"])
                if not (id_col and x_col and y_col):
                    raise ValueError(f"Could not find ID and Easting/Northing columns in {sites_p}. Found: {fnames}")

                for r in reader:
                    raw_id = r.get(id_col, "").strip()
                    raw_x = r.get(x_col, "").strip()
                    raw_y = r.get(y_col, "").strip()
                    if not raw_id or not raw_x or not raw_y:
                        continue
                    try:
                        px, py = float(raw_x), float(raw_y)
                    except ValueError:
                        continue

                    if (abs(px) <= 180.0 and abs(py) <= 90.0) or (custom_crs == 4326):
                        px, py = wgs84_to_bng(py, px)
                    elif custom_crs not in (27700, 4326):
                        px, py = transform_coords_batch([(px, py)], custom_crs, 27700)[0]

                    site_meta = {"wims_id": raw_id, "wims_region_code": "CUSTOM"}
                    pts[raw_id] = (px, py, site_meta)
                    pts_source[raw_id] = src_label
                    core_wq_sites.add(raw_id)

            if samples_p and os.path.exists(samples_p):
                print(f"   Loading custom samples from: {samples_p}")
                with open(samples_p, "r", encoding="utf-8-sig") as f:
                    s_reader = csv.DictReader(f)
                    s_fnames = s_reader.fieldnames or []
                    s_id_col = find_col(s_fnames, ["site_id", "network_id", "station_id", "id", "location_id", "site"])
                    dt_col = find_col(s_fnames, ["sample_datetime", "datetime", "date_time", "date", "survey_date", "timestamp"])
                    val_col = find_col(s_fnames, ["meas_result", "result", "value", "val", "reading", "score", "concentration"])
                    unc_col = find_col(s_fnames, ["meas_uncertainty", "uncertainty", "se", "std_err", "standard_error", "sigma_error"])
                    meth_col = find_col(s_fnames, ["meas_method", "method", "model_type", "type"])
                    det_col_cand = find_col(s_fnames, ["det_code", "determinand_code", "param_code", "cas_number", "determinand", "code", "metric_name"])
                    sign_col = find_col(s_fnames, ["meas_sign", "less_than_flag", "censored", "is_censored", "qualifier"])
                    if not (s_id_col and val_col):
                        raise ValueError(f"Could not find site ID and measurement value columns in {samples_p}. Found: {s_fnames}")

                    for r in s_reader:
                        s_nid = r.get(s_id_col, "").strip()
                        if not s_nid:
                            continue
                        if det_col_cand and det_code:
                            row_det = r.get(det_col_cand, "").strip()
                            if row_det and row_det.lower() != det_code.lower():
                                continue
                        raw_val_str = r.get(val_col, "").strip()
                        if not raw_val_str:
                            continue
                        try:
                            val = float(raw_val_str)
                        except ValueError:
                            continue
                        is_cens = False
                        if sign_col:
                            s_sign = r.get(sign_col, "").strip().lower()
                            is_cens = (s_sign in ("<", "yes", "true", "1"))
                        yr, mo = _extract_year_month(r.get(dt_col, "").strip() if dt_col else "")
                        seas = MONTH_TO_SEASON.get(mo, "summer")
                        unc_val = _opt_float(r.get(unc_col, "").strip()) if unc_col else None
                        meth_val = r.get(meth_col, "").strip().lower() if meth_col else ""
                        all_records.append((s_nid, val, is_cens, yr, mo, seas, unc_val, meth_val))
            continue

        # 2. Registered Source Profile (from configs/sources/*.json)
        if src_item in sources_reg and src_item not in ("rsn", "ssn"):
            cfg_s = sources_reg[src_item]
            s_dir = cfg_s.get("data_dir", "")
            if not os.path.isabs(s_dir):
                cand_a = _rel_join(data_dir, s_dir)
                cand_b = _rel_join(BASE_DIR, s_dir)
                s_dir = cand_a if os.path.isdir(cand_a) else (cand_b if os.path.isdir(cand_b) else s_dir)
            if not os.path.isdir(s_dir):
                raise FileNotFoundError(f"Data directory for source '{src_item}' not found at: {s_dir}")

            sites_f = cfg_s.get("sites_file", "sites.csv")
            sites_p = os.path.join(s_dir, sites_f)
            if not os.path.exists(sites_p):
                raise FileNotFoundError(f"Sites file '{sites_f}' not found in {s_dir}")

            id_col_c = cfg_s.get("id_col", "site_id")
            x_col_c = cfg_s.get("easting_col", "easting")
            y_col_c = cfg_s.get("northing_col", "northing")
            crs_s = int(cfg_s.get("crs_epsg", 27700))

            with open(sites_p, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                fnames = reader.fieldnames or []
                id_col = id_col_c if id_col_c in fnames else find_col(fnames, ["site_id", "network_id", "station_id", "id"])
                x_col = x_col_c if x_col_c in fnames else find_col(fnames, ["ngr_easting", "easting", "x", "lon"])
                y_col = y_col_c if y_col_c in fnames else find_col(fnames, ["ngr_northing", "northing", "y", "lat"])
                if not (id_col and x_col and y_col):
                    raise ValueError(f"Could not find ID and coordinate columns in {sites_p}. Found: {fnames}")

                for r in reader:
                    raw_id = r.get(id_col, "").strip()
                    raw_x = r.get(x_col, "").strip()
                    raw_y = r.get(y_col, "").strip()
                    if not raw_id or not raw_x or not raw_y:
                        continue
                    try:
                        px, py = float(raw_x), float(raw_y)
                    except ValueError:
                        continue
                    if (abs(px) <= 180.0 and abs(py) <= 90.0) or crs_s == 4326:
                        px, py = wgs84_to_bng(py, px)
                    elif crs_s not in (27700, 4326):
                        px, py = transform_coords_batch([(px, py)], crs_s, 27700)[0]

                    site_meta = {"wims_id": raw_id, "wims_region_code": cfg_s.get("abbr", "SRC")}
                    pts[raw_id] = (px, py, site_meta)
                    pts_source[raw_id] = cfg_s.get("abbr", src_label)
                    core_wq_sites.add(raw_id)

            # Check samples dataset in registry
            ds_cfg = cfg_s.get("datasets", {}).get(dataset, {})
            samp_f = ds_cfg.get("file", "samples.csv")
            samp_p = os.path.join(s_dir, samp_f)
            if os.path.exists(samp_p):
                with open(samp_p, "r", encoding="utf-8-sig") as f:
                    s_reader = csv.DictReader(f)
                    s_fnames = s_reader.fieldnames or []
                    s_id_col = ds_cfg.get("id_col") or find_col(s_fnames, ["network_id", "site_id", "station_id", "id"])
                    dt_col = ds_cfg.get("datetime_col") or find_col(s_fnames, ["sample_datetime", "datetime", "date", "survey_date"])
                    val_col = ds_cfg.get("val_col") or find_col(s_fnames, ["meas_result", "result", "value", "val"])
                    det_col_c = ds_cfg.get("det_col") or find_col(s_fnames, ["det_code", "determinand_code", "determinand", "code"])
                    cens_col_c = ds_cfg.get("cens_col") or find_col(s_fnames, ["meas_sign", "less_than_flag", "is_censored"])
                    cens_flg = ds_cfg.get("cens_flag", "<").lower()
                    unc_col = ds_cfg.get("unc_col") or find_col(s_fnames, ["meas_uncertainty", "uncertainty", "se", "std_err", "standard_error", "sigma_error"])
                    meth_col = ds_cfg.get("meth_col") or find_col(s_fnames, ["meas_method", "method", "model_type", "type"])
                    if not (s_id_col and val_col):
                        raise ValueError(f"Could not find site ID and measurement value columns in {samp_p}. Found: {s_fnames}")

                    for r in s_reader:
                        s_nid = r.get(s_id_col, "").strip()
                        if not s_nid:
                            continue
                        if det_col_c and det_code:
                            row_det = r.get(det_col_c, "").strip()
                            if row_det and row_det.lower() != det_code.lower():
                                continue
                        raw_v = r.get(val_col, "").strip()
                        if not raw_v:
                            continue
                        try:
                            val = float(raw_v)
                        except ValueError:
                            continue
                        is_cens = False
                        if cens_col_c:
                            s_sign = r.get(cens_col_c, "").strip().lower()
                            is_cens = (s_sign == cens_flg or s_sign in ("<", "yes", "true", "1"))
                        yr, mo = _extract_year_month(r.get(dt_col, "").strip() if dt_col else "")
                        seas = MONTH_TO_SEASON.get(mo, "summer")
                        unc_val = _opt_float(r.get(unc_col, "").strip()) if unc_col else None
                        meth_val = r.get(meth_col, "").strip().lower() if meth_col else str(cfg_s.get("meas_method", "")).lower()
                        all_records.append((s_nid, val, is_cens, yr, mo, seas, unc_val, meth_val))
            continue

        # 3. Standard RSN / SSN official datasets or auto-discovered custom folder
        if src_item not in ("rsn", "ssn"):
            cand_cs = os.path.join(data_dir, "custom", src_item)
            if not os.path.isdir(cand_cs):
                cand_cs = os.path.join(BASE_DIR, "data", "custom", src_item)
            if not os.path.isdir(cand_cs):
                valid_srcs = sorted(list(sources_reg.keys()) + ["rsn", "ssn", "rss"])
                raise ValueError(f"Unknown monitoring data source '{src_item}'. Available sources: {valid_srcs}")

        s_dir = rsn_dir if src_item == "rsn" else (ssn_dir if src_item == "ssn" else cand_cs)
        prefix = "RSN" if src_item == "rsn" else ("SSN" if src_item == "ssn" else src_item.upper())
        if not os.path.isdir(s_dir):
            raise FileNotFoundError(f"Directory for {src_label} data not found at: {s_dir}")

        wq_pts_path = os.path.join(s_dir, f"{prefix} Water Quality Sample Points.csv")
        if os.path.exists(wq_pts_path):
            with open(wq_pts_path, "r", encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    if r.get("ngr_easting") and r.get("ngr_northing"):
                        nid = r["network_id"]
                        pts[nid] = (float(r["ngr_easting"]), float(r["ngr_northing"]), r)
                        pts_source[nid] = prefix

        eco_pts_path = os.path.join(s_dir, f"{prefix} Ecology Sample Points.csv")
        if os.path.exists(eco_pts_path):
            with open(eco_pts_path, "r", encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    nid = r["network_id"]
                    if nid not in pts and r.get("ngr_easting") and r.get("ngr_northing"):
                        eco_meta = {"wims_id": nid, "wims_region_code": ""}
                        pts[nid] = (float(r["ngr_easting"]), float(r["ngr_northing"]), eco_meta)
                        pts_source[nid] = prefix

        wq_samples_path = os.path.join(s_dir, f"{prefix} Water Quality Samples.csv")
        if dataset == "wq":
            if os.path.exists(wq_samples_path):
                with open(wq_samples_path, "r", encoding="utf-8-sig") as f:
                    for r in csv.DictReader(f):
                        r_det = r.get("det_code")
                        if r_det == det_code:
                            core_wq_sites.add(r["network_id"])
                            try:
                                val = float(r["meas_result"])
                            except (ValueError, TypeError):
                                continue
                            is_cens = (r.get("meas_sign") == "<")
                            dt = r.get("sample_datetime", "")
                            if len(dt) >= 7:
                                yr = int(dt[0:4])
                                mo = int(dt[5:7])
                                seas = MONTH_TO_SEASON[mo]
                                all_records.append((r["network_id"], val, is_cens, yr, mo, seas, None, ""))
        else:
            if os.path.exists(wq_samples_path):
                with open(wq_samples_path, "r", encoding="utf-8-sig") as f:
                    for r in csv.DictReader(f):
                        if r.get("det_code") == "0117":
                            core_wq_sites.add(r["network_id"])
            if dataset == "lcms":
                lcms_path = os.path.join(s_dir, f"{prefix} LCMS GCMS Samples.csv")
                if not os.path.exists(lcms_path):
                    if src_item == "ssn":
                        raise ValueError("LCMS/GCMS data is only monitored on RSN and is not available in Small Streams Network (SSN).")
                    continue
                with open(lcms_path, "r", encoding="utf-8-sig") as f:
                    for r in csv.DictReader(f):
                        if r.get("cas_number") == det_code or (r.get("compound_name") and det_code.lower() in r["compound_name"].lower()):
                            if not r.get("meas_result"):
                                continue
                            try:
                                val = float(r["meas_result"])
                            except (ValueError, TypeError):
                                continue
                            is_cens = (r.get("less_than_flag", "").strip().lower() == "yes")
                            dt = r.get("sample_datetime", "")
                            if len(dt) >= 7:
                                yr = int(dt[0:4])
                                mo = int(dt[5:7])
                                seas = MONTH_TO_SEASON[mo]
                                all_records.append((r["network_id"], val, is_cens, yr, mo, seas, None, ""))
                                core_wq_sites.add(r["network_id"])
            else:
                ds_files = {
                    "inverts": (f"{prefix} Macroinvertebrate Samples.csv", "sample_datetime"),
                    "diatoms": (f"{prefix} Diatom Samples.csv", "sample_datetime"),
                    "macrophytes": (f"{prefix} Macrophyte Samples.csv", "sample_datetime"),
                    "rhs": (f"{prefix} River Habitat Surveys.csv", "survey_date"),
                }
                csv_name, dt_col = ds_files[dataset]
                csv_path = os.path.join(s_dir, csv_name)
                if not os.path.exists(csv_path):
                    if src_item == "ssn":
                        raise ValueError(f"Dataset '{dataset}' is only monitored on RSN and is not available in Small Streams Network (SSN).")
                    continue
                col_key = "whpt_aspt" if (dataset == "inverts" and det_code == "0117") else det_code
                with open(csv_path, "r", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f)
                    if col_key not in reader.fieldnames:
                        raise ValueError(f"Column '{col_key}' not found in {csv_name}. Available: {reader.fieldnames}")
                    for r in reader:
                        raw_s = r.get(col_key, "").strip()
                        dt = r.get(dt_col, "").strip()
                        if not raw_s or not dt or len(dt) < 7:
                            continue
                        try:
                            val = float(raw_s)
                        except (ValueError, TypeError):
                            continue
                        yr = int(dt[0:4])
                        mo = int(dt[5:7])
                        seas = MONTH_TO_SEASON[mo]
                        all_records.append((r["network_id"], val, False, yr, mo, seas, None, ""))
                        core_wq_sites.add(r["network_id"])

    return pts, pts_source, core_wq_sites, all_records


def resolve_gis_path(rel_p, gis_dir=None):
    """
    Resolves a GeoPackage or ZIP path against the configured `GIS_DIR` (default `gis`).
    If `rel_p` begins with `gis/`, that prefix is replaced by `GIS_DIR`.
    """
    g_dir = gis_dir or GIS_DIR
    if os.path.isabs(rel_p):
        return rel_p
    norm_p = rel_p.replace("\\", "/")
    if norm_p.startswith("gis/"):
        sub_p = norm_p[4:]
        cand = _rel_join(g_dir, sub_p)
        if os.path.exists(cand):
            return cand
        return cand
    cand_direct = rel_p
    if os.path.exists(cand_direct):
        return cand_direct
    cand_gis = _rel_join(g_dir, rel_p)
    if os.path.exists(cand_gis):
        return cand_gis
    return _rel_join(BASE_DIR, rel_p)

# Shared transforms, ROS censoring, and GPB/WKB geometry helpers are imported from s3n_py_utils.


# ==============================================================================
# GEOPARQUET (.parquet / .geoparquet / .pq) READER & METADATA HELPERS
# ==============================================================================

PARQUET_EXTENSIONS = (".parquet", ".geoparquet", ".pq")
_GDAL_LIB = None


def is_parquet_path(path):
    """Returns True if `path` is a GeoParquet file or a directory containing Parquet files."""
    if not path:
        return False
    low = str(path).lower()
    if low.endswith(PARQUET_EXTENSIONS):
        return True
    if os.path.isdir(path):
        for fn in os.listdir(path):
            if fn.lower().endswith(PARQUET_EXTENSIONS):
                return True
    return False


def _extract_epsg_from_crs_obj(crs_obj):
    """Extracts an integer EPSG code from a GeoParquet PROJJSON `crs` object."""
    if crs_obj is None:
        return None
    if isinstance(crs_obj, int):
        return crs_obj
    if isinstance(crs_obj, str):
        s = crs_obj.strip().upper()
        if s.startswith("EPSG:") and s[5:].isdigit():
            return int(s[5:])
        if "CRS84" in s or "WGS84" in s:
            return 4326
        return None
    if isinstance(crs_obj, dict):
        cid = crs_obj.get("id")
        if isinstance(cid, dict):
            auth = str(cid.get("authority", "")).upper()
            code = cid.get("code")
            if auth == "EPSG" and str(code).isdigit():
                return int(code)
            if auth == "OGC" and str(code).upper() == "CRS84":
                return 4326
        # Search nested identifiers if present
        ids = crs_obj.get("ids")
        if isinstance(ids, list):
            for item in ids:
                if isinstance(item, dict) and str(item.get("authority", "")).upper() == "EPSG":
                    if str(item.get("code", "")).isdigit():
                        return int(item["code"])
    return None


def read_geoparquet_footer_metadata(parquet_path):
    """
    Pure-Python reader that extracts the OGC GeoParquet `'geo'` JSON metadata
    from the Parquet file footer (`PAR1` magic + 4-byte footer length) in < 1 ms.
    Returns a dict with `'primary_column'`, `'crs_epsg'`, `'encoding'`, `'bbox'`, and raw `'geo'`.
    """
    info = {
        "primary_column": None,
        "crs_epsg": None,
        "encoding": "WKB",
        "bbox": None,
        "geo": None,
    }
    if not os.path.isfile(parquet_path) or os.path.getsize(parquet_path) < 12:
        return info
    try:
        with open(parquet_path, "rb") as f:
            f.seek(-8, os.SEEK_END)
            flen, magic = struct.unpack("<I4s", f.read(8))
            if magic != b"PAR1" or flen <= 0 or flen > 64 * 1024 * 1024:
                return info
            f.seek(-(8 + flen), os.SEEK_END)
            footer = f.read(flen)

        # Locate the GeoParquet "geo" JSON payload inside the Parquet FileMetaData key-value list
        for marker in (b'"primary_column"', b'"columns"'):
            idx = footer.find(marker)
            if idx != -1:
                start = footer.rfind(b"{", 0, idx)
                while start != -1:
                    depth = 0
                    in_str = False
                    esc = False
                    end_pos = -1
                    for i in range(start, len(footer)):
                        b = footer[i]
                        if in_str:
                            if esc:
                                esc = False
                            elif b == 0x5C:  # backslash
                                esc = True
                            elif b == 0x22:  # quote
                                in_str = False
                        else:
                            if b == 0x22:
                                in_str = True
                            elif b == 0x7B:  # '{'
                                depth += 1
                            elif b == 0x7D:  # '}'
                                depth -= 1
                                if depth == 0:
                                    end_pos = i + 1
                                    break
                    if end_pos != -1:
                        try:
                            cand = json.loads(footer[start:end_pos].decode("utf-8"))
                            if isinstance(cand, dict) and "columns" in cand and isinstance(cand["columns"], dict):
                                info["geo"] = cand
                                break
                        except Exception:
                            pass
                    start = footer.rfind(b"{", 0, start)
                if info["geo"] is not None:
                    break

        geo = info["geo"]
        if isinstance(geo, dict):
            cols_meta = geo.get("columns", {})
            prim = geo.get("primary_column")
            if not prim and cols_meta:
                prim = next(iter(cols_meta.keys()))
            info["primary_column"] = prim
            if prim and prim in cols_meta and isinstance(cols_meta[prim], dict):
                cmeta = cols_meta[prim]
                info["encoding"] = cmeta.get("encoding", "WKB")
                info["bbox"] = cmeta.get("bbox")
                if "crs" in cmeta:
                    if cmeta["crs"] is None:
                        # Per OGC GeoParquet spec, explicit null/omitted crs means OGC:CRS84 (EPSG:4326)
                        info["crs_epsg"] = 4326
                    else:
                        info["crs_epsg"] = _extract_epsg_from_crs_obj(cmeta["crs"])
    except Exception:
        pass
    return info


def _get_gdal_lib():
    """
    Lazily loads `libgdal` via `ctypes` and configures OGR C API signatures for
    reading and writing GeoParquet (`.parquet` / `.geoparquet`) files.
    """
    global _GDAL_LIB
    if _GDAL_LIB is not None:
        return _GDAL_LIB

    candidates = [
        ctypes.util.find_library("gdal"),
        "libgdal.so.37",
        "libgdal.so.36",
        "libgdal.so.35",
        "libgdal.so.34",
        "libgdal.so",
        "gdal309.dll",
        "gdal308.dll",
        "gdal307.dll",
        "gdal306.dll",
        "gdal.dll",
    ]
    gdal = None
    for cand in candidates:
        if not cand:
            continue
        try:
            gdal = ctypes.CDLL(cand)
            break
        except OSError:
            continue

    if gdal is None:
        raise RuntimeError(
            "Could not load libgdal shared library or pyarrow to read GeoParquet file. "
            "Install `gdal-libs` or `pyarrow`."
        )

    gdal.GDALAllRegister()
    gdal.CPLSetConfigOption.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    gdal.CPLSetConfigOption.restype = None
    gdal.CPLSetConfigOption(b"OGR_PARQUET_ALLOW_ALL_DIMS", b"YES")

    gdal.GDALOpenEx.argtypes = [
        ctypes.c_char_p,
        ctypes.c_uint,
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.POINTER(ctypes.c_char_p),
        ctypes.POINTER(ctypes.c_char_p),
    ]
    gdal.GDALOpenEx.restype = ctypes.c_void_p
    gdal.GDALClose.argtypes = [ctypes.c_void_p]
    gdal.GDALClose.restype = None
    gdal.GDALDatasetGetLayerCount.argtypes = [ctypes.c_void_p]
    gdal.GDALDatasetGetLayerCount.restype = ctypes.c_int
    gdal.GDALDatasetGetLayer.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.GDALDatasetGetLayer.restype = ctypes.c_void_p
    gdal.GDALDatasetGetLayerByName.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    gdal.GDALDatasetGetLayerByName.restype = ctypes.c_void_p

    gdal.OGR_L_GetName.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_GetName.restype = ctypes.c_char_p
    gdal.OGR_L_GetGeometryColumn.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_GetGeometryColumn.restype = ctypes.c_char_p
    gdal.OGR_L_GetSpatialRef.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_GetSpatialRef.restype = ctypes.c_void_p
    gdal.OSRGetAuthorityCode.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    gdal.OSRGetAuthorityCode.restype = ctypes.c_char_p

    gdal.OGR_L_GetLayerDefn.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_GetLayerDefn.restype = ctypes.c_void_p
    gdal.OGR_FD_GetFieldCount.argtypes = [ctypes.c_void_p]
    gdal.OGR_FD_GetFieldCount.restype = ctypes.c_int
    gdal.OGR_FD_GetFieldDefn.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.OGR_FD_GetFieldDefn.restype = ctypes.c_void_p
    gdal.OGR_Fld_GetNameRef.argtypes = [ctypes.c_void_p]
    gdal.OGR_Fld_GetNameRef.restype = ctypes.c_char_p
    gdal.OGR_Fld_GetType.argtypes = [ctypes.c_void_p]
    gdal.OGR_Fld_GetType.restype = ctypes.c_int

    gdal.OGR_L_ResetReading.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_ResetReading.restype = None
    gdal.OGR_L_GetNextFeature.argtypes = [ctypes.c_void_p]
    gdal.OGR_L_GetNextFeature.restype = ctypes.c_void_p
    gdal.OGR_F_Destroy.argtypes = [ctypes.c_void_p]
    gdal.OGR_F_Destroy.restype = None
    gdal.OGR_F_GetFID.argtypes = [ctypes.c_void_p]
    gdal.OGR_F_GetFID.restype = ctypes.c_longlong
    gdal.OGR_F_IsFieldSetAndNotNull.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.OGR_F_IsFieldSetAndNotNull.restype = ctypes.c_int
    gdal.OGR_F_GetFieldAsInteger64.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.OGR_F_GetFieldAsInteger64.restype = ctypes.c_longlong
    gdal.OGR_F_GetFieldAsDouble.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.OGR_F_GetFieldAsDouble.restype = ctypes.c_double
    gdal.OGR_F_GetFieldAsString.argtypes = [ctypes.c_void_p, ctypes.c_int]
    gdal.OGR_F_GetFieldAsString.restype = ctypes.c_char_p

    gdal.OGR_F_GetGeometryRef.argtypes = [ctypes.c_void_p]
    gdal.OGR_F_GetGeometryRef.restype = ctypes.c_void_p
    gdal.OGR_G_WkbSize.argtypes = [ctypes.c_void_p]
    gdal.OGR_G_WkbSize.restype = ctypes.c_int
    gdal.OGR_G_ExportToIsoWkb.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
    gdal.OGR_G_ExportToIsoWkb.restype = ctypes.c_int

    gdal.GDALVectorTranslateOptionsNew.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p]
    gdal.GDALVectorTranslateOptionsNew.restype = ctypes.c_void_p
    gdal.GDALVectorTranslateOptionsFree.argtypes = [ctypes.c_void_p]
    gdal.GDALVectorTranslateOptionsFree.restype = None
    gdal.GDALVectorTranslate.argtypes = [
        ctypes.c_char_p,
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
    ]
    gdal.GDALVectorTranslate.restype = ctypes.c_void_p

    _GDAL_LIB = gdal
    return gdal


def inspect_geoparquet_schema(parquet_path):
    """
    Inspects a GeoParquet (`.parquet` / `.geoparquet` / `.pq`) file to determine:
      - `layer`    : layer/table name
      - `geom_col` : geometry column name
      - `crs_epsg` : integer EPSG code
      - `cols`     : dict mapping lowercase column name -> original column name
    """
    target_file = parquet_path
    if os.path.isdir(parquet_path):
        pq_files = sorted(
            os.path.join(parquet_path, fn)
            for fn in os.listdir(parquet_path)
            if fn.lower().endswith(PARQUET_EXTENSIONS)
        )
        if not pq_files:
            raise FileNotFoundError(f"No .parquet files found inside directory: {parquet_path}")
        target_file = pq_files[0]

    footer_meta = read_geoparquet_footer_metadata(target_file)
    geom_col = footer_meta.get("primary_column")
    crs_epsg = footer_meta.get("crs_epsg")
    layer_name = os.path.splitext(os.path.basename(target_file))[0]
    cols = {}

    # Try pyarrow first if installed
    try:
        import pyarrow.parquet as pq  # type: ignore

        pf = pq.ParquetFile(target_file)
        for name in pf.schema.names:
            cols[name.lower()] = name
        if not geom_col:
            for cand in ("shape", "geometry", "geom", "wkb_geometry"):
                if cand in cols:
                    geom_col = cols[cand]
                    break
        return {
            "layer": layer_name,
            "geom_col": geom_col or "geometry",
            "crs_epsg": crs_epsg or 27700,
            "cols": cols,
        }
    except ImportError:
        pass

    # Fallback to libgdal C OGR API (Parquet driver)
    gdal = _get_gdal_lib()
    ds = gdal.GDALOpenEx(os.fsencode(target_file), 0, None, None, None)
    if not ds:
        raise RuntimeError(f"GDAL failed to open GeoParquet file: {target_file}")
    try:
        lyr = gdal.GDALDatasetGetLayer(ds, 0)
        if lyr:
            lnm = gdal.OGR_L_GetName(lyr)
            if lnm:
                layer_name = lnm.decode("utf-8", errors="ignore")
            gcol_ptr = gdal.OGR_L_GetGeometryColumn(lyr)
            if gcol_ptr:
                gcol_s = gcol_ptr.decode("utf-8", errors="ignore")
                if gcol_s:
                    geom_col = geom_col or gcol_s
                    cols[gcol_s.lower()] = gcol_s
            srs = gdal.OGR_L_GetSpatialRef(lyr)
            if srs and not crs_epsg:
                auth_code = gdal.OSRGetAuthorityCode(srs, None)
                if auth_code:
                    code_s = auth_code.decode("utf-8", errors="ignore")
                    if code_s.isdigit():
                        crs_epsg = int(code_s)
            fdefn = gdal.OGR_L_GetLayerDefn(lyr)
            nflds = gdal.OGR_FD_GetFieldCount(fdefn)
            for i in range(nflds):
                fld = gdal.OGR_FD_GetFieldDefn(fdefn, i)
                fnm = gdal.OGR_Fld_GetNameRef(fld).decode("utf-8", errors="ignore")
                cols[fnm.lower()] = fnm
    finally:
        gdal.GDALClose(ds)

    if not geom_col:
        for cand in ("shape", "geometry", "geom", "wkb_geometry"):
            if cand in cols:
                geom_col = cols[cand]
                break

    return {
        "layer": layer_name,
        "geom_col": geom_col or "geometry",
        "crs_epsg": crs_epsg or 27700,
        "cols": cols,
    }


def iter_geoparquet_rows(parquet_path, geom_col_hint=None):
    """
    Streams rows from a GeoParquet file (or directory of Parquet files), yielding
    `(rdict, geom_wkb_bytes)` for each feature.
    Uses `pyarrow.parquet` if installed, or `libgdal` C OGR API (`Parquet` driver).
    """
    if os.path.isdir(parquet_path):
        pq_files = sorted(
            os.path.join(parquet_path, fn)
            for fn in os.listdir(parquet_path)
            if fn.lower().endswith(PARQUET_EXTENSIONS)
        )
        for pq_f in pq_files:
            yield from iter_geoparquet_rows(pq_f, geom_col_hint=geom_col_hint)
        return

    footer_meta = read_geoparquet_footer_metadata(parquet_path)
    primary_geom_col = footer_meta.get("primary_column") or geom_col_hint

    # Tier 1: pyarrow.parquet (if installed and encoding is WKB)
    if footer_meta.get("encoding", "WKB").upper() == "WKB":
        try:
            import pyarrow.parquet as pq  # type: ignore

            pf = pq.ParquetFile(parquet_path)
            schema_names = pf.schema.names
            gcol = primary_geom_col
            if not gcol or gcol not in schema_names:
                for cand in (geom_col_hint, "Shape", "geometry", "geom", "wkb_geometry"):
                    if cand and cand in schema_names:
                        gcol = cand
                        break
            for batch in pf.iter_batches(batch_size=8192):
                cols_py = {name: batch.column(i).to_pylist() for i, name in enumerate(batch.schema.names)}
                nrows = batch.num_rows
                for r_idx in range(nrows):
                    rdict = {k: v[r_idx] for k, v in cols_py.items()}
                    if "FID" not in rdict:
                        rdict["FID"] = r_idx + 1
                    blob = rdict.get(gcol) if gcol else None
                    yield rdict, blob
            return
        except ImportError:
            pass

    # Tier 2: libgdal C OGR API (supports all GeoParquet 1.0/1.1 encodings & compressions)
    gdal = _get_gdal_lib()
    ds = gdal.GDALOpenEx(os.fsencode(parquet_path), 0, None, None, None)
    if not ds:
        raise RuntimeError(f"Failed to open GeoParquet file: {parquet_path}")
    try:
        n_layers = gdal.GDALDatasetGetLayerCount(ds)
        for l_idx in range(n_layers):
            lyr = gdal.GDALDatasetGetLayer(ds, l_idx)
            if not lyr:
                continue
            gcol_ptr = gdal.OGR_L_GetGeometryColumn(lyr)
            gcol = (
                gcol_ptr.decode("utf-8", errors="ignore")
                if gcol_ptr
                else (primary_geom_col or geom_col_hint or "geometry")
            )
            fdefn = gdal.OGR_L_GetLayerDefn(lyr)
            nflds = gdal.OGR_FD_GetFieldCount(fdefn)
            fields = []
            for i in range(nflds):
                fld = gdal.OGR_FD_GetFieldDefn(fdefn, i)
                fnm = gdal.OGR_Fld_GetNameRef(fld).decode("utf-8", errors="ignore")
                ftyp = gdal.OGR_Fld_GetType(fld)
                fields.append((i, fnm, ftyp))

            gdal.OGR_L_ResetReading(lyr)
            row_num = 0
            while True:
                feat = gdal.OGR_L_GetNextFeature(lyr)
                if not feat:
                    break
                row_num += 1
                try:
                    fid = gdal.OGR_F_GetFID(feat)
                    rdict = {"FID": int(fid) if fid >= 0 else row_num}
                    for i, fnm, ftyp in fields:
                        if not gdal.OGR_F_IsFieldSetAndNotNull(feat, i):
                            rdict[fnm] = None
                        elif ftyp in (0, 12):  # OFTInteger, OFTInteger64
                            rdict[fnm] = int(gdal.OGR_F_GetFieldAsInteger64(feat, i))
                        elif ftyp == 2:  # OFTReal
                            rdict[fnm] = float(gdal.OGR_F_GetFieldAsDouble(feat, i))
                        else:
                            s_ptr = gdal.OGR_F_GetFieldAsString(feat, i)
                            rdict[fnm] = s_ptr.decode("utf-8", errors="ignore") if s_ptr else None

                    geom = gdal.OGR_F_GetGeometryRef(feat)
                    wkb_bytes = None
                    if geom:
                        sz = gdal.OGR_G_WkbSize(geom)
                        if sz > 0:
                            buf = ctypes.create_string_buffer(sz)
                            if gdal.OGR_G_ExportToIsoWkb(geom, 1, buf) == 0:
                                wkb_bytes = buf.raw
                    if gcol:
                        rdict[gcol] = wkb_bytes
                    yield rdict, wkb_bytes
                finally:
                    gdal.OGR_F_Destroy(feat)
    finally:
        gdal.GDALClose(ds)


def convert_gpkg_to_geoparquet(src_gpkg_path, dst_parquet_path, layers=None):
    """
    Converts a GeoPackage (`.gpkg`) layer (or merges multiple layers) into an
    OGC GeoParquet (`.parquet` / `.geoparquet`) file using `libgdal`.
    """
    gdal = _get_gdal_lib()
    src_ds = gdal.GDALOpenEx(os.fsencode(src_gpkg_path), 0, None, None, None)
    if not src_ds:
        raise FileNotFoundError(f"Could not open source GeoPackage: {src_gpkg_path}")
    try:
        os.makedirs(os.path.dirname(os.path.abspath(dst_parquet_path)), exist_ok=True)
        if os.path.exists(dst_parquet_path):
            os.remove(dst_parquet_path)

        opt_tokens = [b"-f", b"Parquet"]
        if layers and len(layers) == 1:
            opt_tokens.append(os.fsencode(layers[0]))
        elif layers and len(layers) > 1:
            # Inspect columns across layers and build a column-aligned UNION ALL query
            conn = sqlite3.connect(src_gpkg_path)
            c = conn.cursor()
            layer_cols = []
            all_cols_ordered = []
            for lyr in layers:
                c_list = [r[1] for r in c.execute(f'PRAGMA table_info("{lyr}")').fetchall()]
                layer_cols.append((lyr, set(c_list)))
                for cn in c_list:
                    if cn not in all_cols_ordered:
                        all_cols_ordered.append(cn)
            conn.close()
            sql_parts = []
            for lyr, c_set in layer_cols:
                sel_exprs = [f'"{cn}"' if cn in c_set else f'NULL AS "{cn}"' for cn in all_cols_ordered]
                sql_parts.append(f'SELECT {", ".join(sel_exprs)} FROM "{lyr}"')
            union_sql = " UNION ALL ".join(sql_parts)
            opt_tokens.extend([b"-sql", os.fsencode(union_sql)])

        c_opts = (ctypes.c_char_p * (len(opt_tokens) + 1))(*opt_tokens, None)
        vt_opts = gdal.GDALVectorTranslateOptionsNew(c_opts, None)
        src_arr = (ctypes.c_void_p * 1)(src_ds)
        err = ctypes.c_int(0)
        dst_ds = gdal.GDALVectorTranslate(
            os.fsencode(dst_parquet_path), None, 1, src_arr, vt_opts, ctypes.byref(err)
        )
        gdal.GDALVectorTranslateOptionsFree(vt_opts)
        if dst_ds:
            gdal.GDALClose(dst_ds)
        if err.value != 0 or not os.path.exists(dst_parquet_path):
            raise RuntimeError(f"GDALVectorTranslate to GeoParquet failed (error={err.value}).")
    finally:
        gdal.GDALClose(src_ds)


def _get_config_network_paths(cfg):
    """Returns configured network file paths from `network_paths`, `parquet_paths`, or `gpkg_paths`."""
    paths = []
    for key in ("network_paths", "parquet_paths", "gpkg_paths"):
        val = cfg.get(key)
        if isinstance(val, list):
            for item in val:
                if item and item not in paths:
                    paths.append(item)
        elif isinstance(val, str) and val.strip():
            if val.strip() not in paths:
                paths.append(val.strip())
    return paths


def list_available_networks():
    ensure_network_configs(BASE_DIR)
    print(f"Available Water Network Profiles (`--network=<name>`, GIS folder: `{GIS_DIR}`):\n")
    print(f"{'Network Name':<18} | {'Short Label':<20} | {'CRS':<10} | {'Network File Path(s) (.gpkg / .parquet)'}")
    print("-" * 105)
    for p in sorted(glob.glob(os.path.join(CONFIG_DIR, "*.json"))):
        with open(p, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        net_files = ", ".join(resolve_gis_path(gp) for gp in _get_config_network_paths(cfg))
        print(f"{cfg['abbr']:<18} | {cfg.get('short_label', cfg['abbr']):<20} | EPSG:{cfg.get('crs_epsg', 27700):<5} | {net_files}")
    print(
        "\nYou can also pass a custom JSON config (`--network-config=path/to/net.json`), "
        "a `.gpkg` or `.parquet`/`.geoparquet` file path (`--network=path/to/network.parquet`), or `--gis-dir=<dir>`."
    )


def _build_auto_network_config(file_path, layer, geom_col, crs_epsg, cols, is_parquet=False):
    """Builds a complete network configuration dict from an inspected GeoPackage or GeoParquet schema."""
    stem = os.path.splitext(os.path.basename(file_path.rstrip("/\\")))[0]
    if stem.lower().endswith(".geo"):
        stem = stem[:-4]
    abbr = stem.lower()
    path_key = "parquet_paths" if is_parquet else "gpkg_paths"
    return {
        "name": abbr,
        "abbr": abbr,
        "label": abbr,
        "short_label": abbr,
        path_key: [file_path],
        "gpkg_paths": [file_path],
        "layers": [layer],
        "geom_col": geom_col,
        "crs_epsg": crs_epsg or 27700,
        "topology": {
            "mode": "pointer" if "nextdownid" in cols else ("nodes" if "fnode" in cols else "geometry"),
            "id_col": cols.get(
                "rivid",
                cols.get("object_id", cols.get("objectid", cols.get("id", cols.get("hyriv_id", "FID")))),
            ),
            "from_node_col": cols.get("fnode", cols.get("startnode", cols.get("fromnode"))),
            "to_node_col": cols.get("tnode", cols.get("endnode", cols.get("tonode"))),
            "next_down_col": cols.get("nextdownid", cols.get("next_down", cols.get("tonode_id"))),
            "catch_id_col": cols.get("catchid", cols.get("catch_id", cols.get("hybas_id"))),
            "rbd_id_col": cols.get("rbd_id"),
        },
        "attributes": {
            "dist_to_mouth_col": cols.get("dist2mth", cols.get("dist_dn_km")),
            "dist_to_mouth_unit": "km" if "dist_dn_km" in cols else "m",
            "afv_weight_col": cols.get("us_accum", cols.get("cum_len", cols.get("upland_skm"))),
            "afv_weight_unit": "km" if "upland_skm" in cols else "m",
            "shreve_col": cols.get("shreve"),
            "strahler_col": cols.get("strahler", cols.get("ord_stra")),
            "slope_col": cols.get("slope", cols.get("pente")),
            "elevation_col": cols.get("elevation", cols.get("elev")),
            "river_name_cols": [c_nm for k, c_nm in cols.items() if "name" in k],
        },
        "covariates": {"mode": "geometry_and_spatial_join"},
        "max_snap_dist_m": 250.0,
    }


def load_network_config(network_arg, config_path=None):
    if config_path:
        p = config_path if os.path.isabs(config_path) else _rel_join(BASE_DIR, config_path)
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)

    ensure_network_configs(BASE_DIR)

    net_clean = network_arg.strip().lower().replace("-", "_")
    cand_json = os.path.join(CONFIG_DIR, f"{net_clean}.json")
    if os.path.exists(cand_json):
        with open(cand_json, "r", encoding="utf-8") as f:
            return json.load(f)

    if net_clean in MINIMAL_DEFAULT_NETWORKS:
        return MINIMAL_DEFAULT_NETWORKS[net_clean]

    if network_arg.strip().endswith(".json") and os.path.exists(network_arg.strip()):
        with open(network_arg.strip(), "r", encoding="utf-8") as f:
            return json.load(f)

    # Direct GeoParquet or GeoPackage path fallback (or stem inside GIS_DIR): auto-inspect schema
    cand_paths = [resolve_gis_path(net_clean)]
    if not any(net_clean.lower().endswith(ext) for ext in (".gpkg",) + PARQUET_EXTENSIONS):
        for ext in PARQUET_EXTENSIONS + (".gpkg",):
            cand_paths.append(resolve_gis_path(net_clean + ext))

    for cand_p in cand_paths:
        if not os.path.exists(cand_p):
            continue
        if is_parquet_path(cand_p):
            pq_info = inspect_geoparquet_schema(cand_p)
            return _build_auto_network_config(
                file_path=cand_p,
                layer=pq_info["layer"],
                geom_col=pq_info["geom_col"],
                crs_epsg=pq_info["crs_epsg"],
                cols=pq_info["cols"],
                is_parquet=True,
            )
        if cand_p.lower().endswith(".gpkg"):
            conn = sqlite3.connect(cand_p)
            c = conn.cursor()
            row = c.execute("SELECT table_name, srs_id FROM gpkg_contents WHERE data_type='features' LIMIT 1").fetchone()
            layer, srs_id = row[0], row[1]
            cols = {r[1].lower(): r[1] for r in c.execute(f'PRAGMA table_info("{layer}")').fetchall()}
            conn.close()
            geom_col = cols.get("shape", cols.get("geom", cols.get("geometry", "Shape")))
            return _build_auto_network_config(
                file_path=cand_p,
                layer=layer,
                geom_col=geom_col,
                crs_epsg=srs_id or 27700,
                cols=cols,
                is_parquet=False,
            )

    raise ValueError(f"Unknown network '{network_arg}'. Run with --list-networks to see available options.")


def ensure_gpkg_paths(cfg):
    """Ensures all configured GeoPackage or GeoParquet files exist, extracting zip_fallbacks if needed."""
    resolved = []
    zip_list = cfg.get("zip_fallbacks", [])
    raw_paths = _get_config_network_paths(cfg)
    for idx, rel_p in enumerate(raw_paths):
        full_p = resolve_gis_path(rel_p)
        if not os.path.exists(full_p) and idx < len(zip_list):
            zip_rel = zip_list[idx]
            zip_p = resolve_gis_path(zip_rel)
            if os.path.exists(zip_p):
                print(f"   Extracting {zip_p}...")
                with zipfile.ZipFile(zip_p, "r") as zf:
                    zf.extractall(os.path.dirname(zip_p))
        if not os.path.exists(full_p):
            raise FileNotFoundError(f"Water network file not found: {full_p}")
        resolved.append(full_p)
    return resolved


def load_orn_spatial_reference(runtime_cfg=None):
    """
    Loads the pre-extracted ORN England reaches (`<processed_dir>/orn_rsn/england_rsn_reaches.csv`)
    into a fast 2 km spatial hash grid in EPSG:27700 so non-ORN networks in England
    can inherit BGS bedrock/superficial geology (`Is_Chalk_Limestone`, `Is_Sandstone`,
    `Is_Aquifer`, `Is_Peat`), `RBD_ID`, and fallback `Elevation`/`Slope` in < 0.2s.
    """
    cfg = runtime_cfg if runtime_cfg is not None else RUNTIME_CFG
    orn_csv = os.path.join(cfg.processed_dir, "orn_rsn", "england_rsn_reaches.csv")
    if not os.path.exists(orn_csv):
        orn_csv = os.path.join(cfg.data_dir, "processed", "orn_rsn", "england_rsn_reaches.csv")
    if not os.path.exists(orn_csv):
        return None, {}
    grid = defaultdict(list)
    with open(orn_csv, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mx = float(r["mid_X"])
            my = float(r["mid_Y"])
            rec = (
                mx,
                my,
                int(r["RBD_ID"] or 0),
                float(r["Elevation"] or 50.0),
                float(r["Slope"] or 0.5),
                float(r["SrcElevation"] or 100.0),
                int(r["Is_Chalk_Limestone"] or 0),
                int(r["Is_Sandstone"] or 0),
                int(r["Is_Aquifer"] or 0),
                int(r["Is_Peat"] or 0),
            )
            grid[(int(mx // 2000), int(my // 2000))].append(rec)

    rbd_names = {}
    orn_lu = os.path.join(os.path.dirname(orn_csv), "rsn_catchments_lookup.csv")
    if os.path.exists(orn_lu):
        with open(orn_lu, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["RBD_ID"] and r["RBD_Name"]:
                    rbd_names[int(r["RBD_ID"])] = r["RBD_Name"]
    return grid, rbd_names


def query_orn_spatial_grid(grid, mx, my):
    """Finds the nearest ORN reach midpoint in EPSG:27700 within 6 km."""
    if grid is None:
        return None
    gx0, gy0 = int(mx // 2000), int(my // 2000)
    best_d2 = 1e18
    best_rec = None
    for radius in (1, 3):
        for gx in range(gx0 - radius, gx0 + radius + 1):
            for gy in range(gy0 - radius, gy0 + radius + 1):
                for rec in grid.get((gx, gy), ()):
                    d2 = (mx - rec[0]) ** 2 + (my - rec[1]) ** 2
                    if d2 < best_d2:
                        best_d2 = d2
                        best_rec = rec
        if best_rec is not None:
            return best_rec
    return best_rec


REACH_OBS_HEADER_COLS = [
    "has_obs", "site_id", "snap_dist_m", "obs_up_dist_km",
    "obs_mean", "obs_median", "trans_obs_mean",
    "meas_uncertainty", "meas_method", "n_samples",
]


def summarize_reach_site_obs(s_list, d2m_km, trans_type="log", eff_lambda=None):
    """Computes reach-level observed summary fields for `england_rsn_reaches.csv`."""
    all_vals = [v for s in s_list for v in s[3]]
    all_uncs = [u for s in s_list for u in s[4]]
    all_meths = {m for s in s_list for m in s[5] if m}
    along_m = statistics.mean(s[2] for s in s_list)
    snap_d = min(s[1] for s in s_list)
    obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
    m_val = statistics.mean(all_vals)
    med_val = statistics.median(all_vals)
    trans_val = apply_forward_transform(m_val, ttype=trans_type, lam=eff_lambda)
    trans_str = round(trans_val, 5) if trans_val is not None else ""
    r_unc = round(math.sqrt(sum(u * u for u in all_uncs) / len(all_uncs)), 6) if all_uncs else ""
    r_meth = sorted(all_meths)[0] if all_meths else ""
    return [
        1, s_list[0][0], round(snap_d, 2), round(obs_up_km, 5),
        round(m_val, 4), round(med_val, 4), trans_str,
        r_unc, r_meth, len(all_vals),
    ]


def write_snapped_sites_and_periods(
    out_dir,
    all_station_ids,
    pts,
    pts_source,
    snap_info,
    imputed_records,
    obs_by_site,
    unc_by_site,
    meth_by_site,
    trans_type="log",
    eff_lambda=None,
):
    """Writes `rsn_sites.csv` and `rsn_reach_periods.csv` and returns `(site_by_reach, sites_csv_path, periods_csv_path)`."""
    site_by_reach = defaultdict(list)
    sites_csv_path = os.path.join(out_dir, "rsn_sites.csv")
    with open(sites_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "site_id", "wims_id", "wims_region_code",
            "easting", "northing", "snap_dist_m", "RivID", "CatchID",
            "obs_up_dist_km", "n_samples", "obs_mean", "obs_median",
            "obs_sd", "trans_obs_mean", "mean_trans_obs",
            "meas_uncertainty", "meas_method", "included_250m",
            "source",
        ])
        for nid in all_station_ids:
            if nid not in snap_info:
                continue
            px, py, meta = pts[nid]
            d, riv_id, catch_id, along_m, d2m_km, inc = snap_info[nid]
            obs_up_km = max(0.0, d2m_km - (along_m / 1000.0))
            src_str = pts_source.get(nid, "RSN")
            vals = obs_by_site.get(nid, [])
            u_vals = unc_by_site.get(nid, [])
            m_set = meth_by_site.get(nid, set())
            site_u = round(math.sqrt(sum(u * u for u in u_vals) / len(u_vals)), 6) if u_vals else ""
            site_m = sorted(m_set)[0] if m_set else ""
            if vals:
                m_val = statistics.mean(vals)
                med_val = statistics.median(vals)
                sd_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                trans_m_val = apply_forward_transform(m_val, ttype=trans_type, lam=eff_lambda)
                trans_m = round(trans_m_val, 5) if trans_m_val is not None else ""
                trans_clean = [
                    tv for tv in (apply_forward_transform(v, ttype=trans_type, lam=eff_lambda) for v in vals)
                    if tv is not None
                ]
                m_trans = round(statistics.mean(trans_clean), 5) if trans_clean else ""
                if inc:
                    site_by_reach[riv_id].append((nid, d, along_m, vals, u_vals, m_set))
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), len(vals), round(m_val, 4), round(med_val, 4),
                    round(sd_val, 4), trans_m, m_trans, site_u, site_m, inc, src_str,
                ])
            else:
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), 0, "", "", "", "", "", "", "", inc, src_str,
                ])

    reach_period_vals = defaultdict(list)
    reach_period_unc = defaultdict(list)
    reach_period_meth = defaultdict(set)
    reach_along_map = defaultdict(list)
    reach_meta_map = {}
    for rec in imputed_records:
        nid, imp_val, yr, mo, seas = rec[0], rec[1], rec[2], rec[3], rec[4]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
        if nid not in snap_info:
            continue
        _d, riv_id, catch_id, along_m, d2m_km, inc = snap_info[nid]
        if inc == 1:
            reach_along_map[riv_id].append(along_m)
            reach_meta_map[riv_id] = (catch_id, d2m_km)
            for p_type, p_val in (
                ("all", "all"),
                ("year", str(yr)),
                ("season", seas),
                ("month", f"{mo:02d}"),
                ("year_season", f"{yr}_{seas}"),
            ):
                reach_period_vals[(riv_id, p_type, p_val)].append(imp_val)
                if unc_val is not None:
                    reach_period_unc[(riv_id, p_type, p_val)].append(unc_val)
                if meth_val:
                    reach_period_meth[(riv_id, p_type, p_val)].add(meth_val)

    periods_csv_path = os.path.join(out_dir, "rsn_reach_periods.csv")
    with open(periods_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "RivID", "CatchID", "period_type", "period_value",
            "obs_up_dist_km", "n_samples", "obs_mean", "obs_median", "trans_obs_mean",
            "meas_uncertainty", "meas_method",
        ])
        for (riv_id, p_type, p_val), vals in sorted(reach_period_vals.items()):
            catch_id, d2m_km = reach_meta_map[riv_id]
            along_m = statistics.mean(reach_along_map[riv_id])
            obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
            m_v = statistics.mean(vals)
            med_v = statistics.median(vals)
            trans_mv_val = apply_forward_transform(m_v, ttype=trans_type, lam=eff_lambda)
            trans_mv = round(trans_mv_val, 5) if trans_mv_val is not None else ""
            u_list = reach_period_unc.get((riv_id, p_type, p_val), [])
            m_set = reach_period_meth.get((riv_id, p_type, p_val), set())
            p_unc = round(math.sqrt(sum(u * u for u in u_list) / len(u_list)), 6) if u_list else ""
            p_meth = sorted(m_set)[0] if m_set else ""
            w.writerow([
                riv_id, catch_id, p_type, p_val,
                round(obs_up_km, 5), len(vals), round(m_v, 4), round(med_v, 4), trans_mv,
                p_unc, p_meth,
            ])

    return site_by_reach, sites_csv_path, periods_csv_path


def write_catchment_lookup(out_dir, rows):
    """Writes `rsn_catchments_lookup.csv` from pre-sorted `[CatchID, Main_River, RBD_ID, RBD_Name, n_reaches, n_obs_reaches]` rows."""
    catch_lu_path = os.path.join(out_dir, "rsn_catchments_lookup.csv")
    with open(catch_lu_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["CatchID", "Main_River", "RBD_ID", "RBD_Name", "n_reaches", "n_obs_reaches"])
        w.writerows(rows)
    return catch_lu_path


def write_reaches_and_finalize(
    out_dir,
    net_abbr,
    out_reaches_path,
    reaches_header,
    reaches_rows,
    meta_obj,
    sites_csv_path,
    periods_csv_path,
    catch_lu_path,
    source_key="rsn",
    runtime_cfg=None,
):
    """Writes `england_rsn_reaches.csv` and `network_metadata.json` directly into `out_dir`."""
    with open(out_reaches_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(reaches_header)
        w.writerows(reaches_rows)

    meta_path = os.path.join(out_dir, "network_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_obj, f, indent=2)
    return meta_path


def preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, obs_by_site, period_label, t0,
                               trans_type="log", eff_lambda=None, pts_source=None, source_key="rsn", range_mode=None,
                               runtime_cfg=None, out_dir=None):
    """
    Network-agnostic preprocessing pipeline for any configured GeoPackage network
    (such as EU-Hydro Thames, EU-Hydro Tweed, combined EU-Hydro, OS Water Network, etc.).
    """
    eff_cfg = runtime_cfg if runtime_cfg is not None else RUNTIME_CFG
    eff_out_dir = out_dir if out_dir is not None else eff_cfg.processed_dir
    if pts_source is None:
        pts_source = {}
    source_key, constituents, source_abbr, source_label = resolve_source_tag(source_key, runtime_cfg=eff_cfg)
    net_abbr = cfg["abbr"]
    net_label = cfg.get("label", net_abbr)
    short_label = cfg.get("short_label", net_label)
    crs_epsg = int(cfg.get("crs_epsg", 27700))
    max_snap_m = float(args.max_snap_dist if args.max_snap_dist is not None else cfg.get("max_snap_dist_m", 250.0))
    gpkg_paths = ensure_gpkg_paths(cfg)
    layers = cfg.get("layers", ["ORN"])
    geom_col = cfg.get("geom_col", "Shape")

    topo_cfg = cfg.get("topology", {})
    topo_mode = args.topology_mode if args.topology_mode != "auto" else topo_cfg.get("mode", "nodes")
    id_col = topo_cfg.get("id_col", "RivID")
    fnode_col = topo_cfg.get("from_node_col", "Fnode")
    tnode_col = topo_cfg.get("to_node_col", "Tnode")
    next_down_col = topo_cfg.get("next_down_col")
    catch_col = topo_cfg.get("catch_id_col", "CatchID")
    rbd_col = topo_cfg.get("rbd_id_col")

    attr_cfg = cfg.get("attributes", {})
    d2m_col = attr_cfg.get("dist_to_mouth_col")
    d2m_scale = 0.001 if attr_cfg.get("dist_to_mouth_unit", "m") == "m" else 1.0
    afv_col = args.afv_col if args.afv_col else attr_cfg.get("afv_weight_col")
    afv_scale = 0.001 if attr_cfg.get("afv_weight_unit", "m") == "m" else 1.0
    shreve_col = attr_cfg.get("shreve_col")
    strahler_col = attr_cfg.get("strahler_col")
    slope_col = attr_cfg.get("slope_col")
    elev_col = attr_cfg.get("elevation_col")
    name_cols = attr_cfg.get("river_name_cols", [])
    extra_num_cols = attr_cfg.get("extra_numeric_cols", [])

    print(
        f"2. Loading network '{net_label}' (abbr='{net_abbr}', CRS=EPSG:{crs_epsg}, "
        f"topology='{topo_mode}', max_snap={max_snap_m:.0f}m)..."
    )

    # Transform RSN stations from EPSG:27700 to network CRS for accurate planar snapping
    all_station_ids = sorted(pts.keys())
    station_bng_xy = [(pts[nid][0], pts[nid][1]) for nid in all_station_ids]
    station_net_xy = transform_coords_batch(station_bng_xy, 27700, crs_epsg)

    site_grid = defaultdict(list)
    for nid, (px, py) in zip(all_station_ids, station_net_xy):
        site_grid[(int(px // 2000), int(py // 2000))].append((nid, px, py))

    best_snap = {nid: (1e18, None, None, 0.0) for nid in all_station_ids}

    # Scan GeoPackage(s) and/or GeoParquet file(s) to snap stations and load reaches
    raw_reaches = []
    node_lookup = {}  # (file_idx, node_id) or ("*", node_id) -> (elev, len_tom)
    nodes_cfg = cfg.get("nodes_table")

    # If nodes_table points to an external GeoParquet or CSV file, load it into node_lookup
    if isinstance(nodes_cfg, dict):
        ext_nodes_rel = nodes_cfg.get("path") or nodes_cfg.get("parquet_path")
        if not ext_nodes_rel:
            cand_tbl = str(nodes_cfg.get("table", ""))
            if cand_tbl.lower().endswith(PARQUET_EXTENSIONS + (".csv", ".gpkg")):
                ext_nodes_rel = cand_tbl
        if ext_nodes_rel:
            ext_nodes_p = resolve_gis_path(ext_nodes_rel)
            if os.path.exists(ext_nodes_p):
                nid_c = nodes_cfg.get("id_col", "OBJECT_ID")
                nel_c = nodes_cfg.get("elev_col", "ELEV")
                nlt_c = nodes_cfg.get("dist_to_mouth_col", "LEN_TOM")
                if is_parquet_path(ext_nodes_p):
                    for nr_dict, _ in iter_geoparquet_rows(ext_nodes_p):
                        if nr_dict.get(nid_c) is not None:
                            node_lookup[("*", str(nr_dict[nid_c]))] = (nr_dict.get(nel_c), nr_dict.get(nlt_c))
                elif ext_nodes_p.lower().endswith(".csv"):
                    with open(ext_nodes_p, "r", encoding="utf-8-sig") as nf:
                        for nr_dict in csv.DictReader(nf):
                            if nr_dict.get(nid_c):
                                el_v = float(nr_dict[nel_c]) if nr_dict.get(nel_c) not in (None, "") else None
                                lt_v = float(nr_dict[nlt_c]) if nr_dict.get(nlt_c) not in (None, "") else None
                                node_lookup[("*", str(nr_dict[nid_c]))] = (el_v, lt_v)

    def _parse_cid_val(raw_val):
        if raw_val is None:
            return 0
        s = str(raw_val).strip()
        if not s:
            return 0
        try:
            return int(float(s))
        except ValueError:
            return (abs(hash(s)) % 90000) + 1

    def _ingest_reach_row(file_idx, rdict, blob):
        raw_id_val = rdict.get(id_col) if id_col else None
        if raw_id_val is None:
            raw_id_val = rdict.get("FID", len(raw_reaches) + 1)
        orig_id = str(raw_id_val)
        raw_cid = rdict.get(catch_col) if catch_col else None
        cid_key = (file_idx, _parse_cid_val(raw_cid))

        env = parse_gpkg_envelope(blob)
        poly = None
        zs = ()
        if env is not None:
            minx, maxx, miny, maxy = env
        else:
            poly, zs = parse_gpkg_geometry(blob)
            if len(poly) < 2:
                return
            xs = [p[0] for p in poly]
            ys = [p[1] for p in poly]
            minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        snap_pad = max(1500.0, max_snap_m + 500.0)
        gx0, gx1 = int((minx - snap_pad) // 2000), int((maxx + snap_pad) // 2000)
        gy0, gy1 = int((miny - snap_pad) // 2000), int((maxy + snap_pad) // 2000)
        cand_sites = []
        for gx in range(gx0, gx1 + 1):
            for gy in range(gy0, gy1 + 1):
                if (gx, gy) in site_grid:
                    cand_sites.extend(site_grid[(gx, gy)])

        if cand_sites:
            if poly is None:
                poly, zs = parse_gpkg_geometry(blob)
                if len(poly) < 2:
                    return
            for nid, px, py in cand_sites:
                d, along = point_to_polyline(px, py, poly)
                if d < best_snap[nid][0]:
                    best_snap[nid] = (d, (file_idx, orig_id), cid_key, along)

        raw_reaches.append((file_idx, orig_id, cid_key, poly, zs, rdict, blob))

    for file_idx, gpkg_p in enumerate(gpkg_paths):
        if is_parquet_path(gpkg_p):
            for rdict, blob in iter_geoparquet_rows(gpkg_p, geom_col_hint=geom_col):
                if blob is None and geom_col in rdict:
                    blob = rdict.get(geom_col)
                _ingest_reach_row(file_idx, rdict, blob)
        else:
            conn = sqlite3.connect(gpkg_p)
            conn.row_factory = sqlite3.Row
            c = conn.cursor()

            if isinstance(nodes_cfg, dict):
                ntbl = nodes_cfg.get("table", "Nodes")
                nid_c = nodes_cfg.get("id_col", "OBJECT_ID")
                nel_c = nodes_cfg.get("elev_col", "ELEV")
                nlt_c = nodes_cfg.get("dist_to_mouth_col", "LEN_TOM")
                try:
                    for nr in c.execute(f'SELECT "{nid_c}", "{nel_c}", "{nlt_c}" FROM "{ntbl}"'):
                        node_lookup[(file_idx, str(nr[0]))] = (nr[1], nr[2])
                except sqlite3.Error:
                    pass

            for layer in layers:
                for row in c.execute(f'SELECT * FROM "{layer}"'):
                    rdict = dict(row)
                    blob = rdict.get(geom_col)
                    _ingest_reach_row(file_idx, rdict, blob)
            conn.close()

    # Determine monitored catchments (or keep all if network has <= 5 catchments)
    monitored_cid_keys = {
        best_snap[nid][2]
        for nid in core_wq_sites
        if nid in best_snap and best_snap[nid][0] <= max_snap_m and best_snap[nid][2] is not None
    }
    if not monitored_cid_keys:
        monitored_cid_keys = {
            best_snap[nid][2]
            for nid in best_snap
            if best_snap[nid][0] <= max_snap_m and best_snap[nid][2] is not None
        }
    if not monitored_cid_keys:
        monitored_cid_keys = {r[2] for r in raw_reaches}

    has_explicit_catch_id = any(r[2][1] > 0 for r in raw_reaches)

    # Filter reaches to monitored catchments
    filtered_reaches = [r for r in raw_reaches if r[2] in monitored_cid_keys]

    # Assign unique integer CatchIDs (preserving original CatchID if single file, else offsetting file_idx > 0)
    cid_map = {}
    used_cids = set()
    for file_idx, orig_cid in sorted(monitored_cid_keys):
        cand_cid = orig_cid if (file_idx == 0 and orig_cid > 0 and orig_cid not in used_cids) else (file_idx * 10000 + orig_cid)
        while cand_cid <= 0 or cand_cid in used_cids:
            cand_cid += 1
        cid_map[(file_idx, orig_cid)] = cand_cid
        used_cids.add(cand_cid)

    # Assign unique positive integer RivIDs (lazily decoding geometry for reaches without nearby sites)
    key_to_rivid = {}
    global_orig_to_rivid = {}
    reaches_by_id = {}
    riv_idx = 0
    for file_idx, orig_id, cid_key, poly, zs, rdict, blob in filtered_reaches:
        if poly is None:
            poly, zs = parse_gpkg_geometry(blob)
        if len(poly) < 2:
            continue
        riv_idx += 1
        if len(gpkg_paths) == 1 and len(layers) == 1 and str(orig_id).isdigit():
            riv_id = int(orig_id)
        else:
            riv_id = riv_idx
        key_to_rivid[(file_idx, orig_id)] = riv_id
        if orig_id not in global_orig_to_rivid:
            global_orig_to_rivid[orig_id] = riv_id
        reaches_by_id[riv_id] = {
            "riv_id": riv_id,
            "file_idx": file_idx,
            "orig_id": orig_id,
            "catch_id": cid_map[cid_key],
            "poly": poly,
            "zs": zs,
            "row": rdict,
        }

    src_fmt_label = (
        "GeoParquet"
        if all(is_parquet_path(p) for p in gpkg_paths)
        else ("GeoPackage + GeoParquet" if any(is_parquet_path(p) for p in gpkg_paths) else "GeoPackage")
    )
    print(
        f"   Loaded {len(reaches_by_id)} reaches across {len(monitored_cid_keys)} monitored catchments "
        f"(from {len(raw_reaches)} total reaches in {src_fmt_label})."
    )

    # 3. Build directed river tree topology (supporting pointer, nodes, or geometry)
    print(f"3. Building directed river trees (mode='{topo_mode}') and verifying reach orientation...")
    fnode_index = {}
    start_pt_grid = defaultdict(list)
    for riv_id, info in reaches_by_id.items():
        f_idx = info["file_idx"]
        rdict = info["row"]
        fn = rdict.get(fnode_col) if fnode_col else None
        if fn:
            # Prefer main channel (MC == 1) when multiple reaches share an FNODE
            mc = rdict.get("MC", 1)
            key = (f_idx, info["catch_id"], str(fn))
            if key not in fnode_index or mc == 1:
                fnode_index[key] = riv_id
            g_key = ("*", info["catch_id"], str(fn))
            if g_key not in fnode_index or mc == 1:
                fnode_index[g_key] = riv_id
        p0 = info["poly"][0]
        start_pt_grid[(f_idx, int(p0[0] // 50), int(p0[1] // 50))].append((riv_id, p0[0], p0[1]))
        start_pt_grid[("*", int(p0[0] // 50), int(p0[1] // 50))].append((riv_id, p0[0], p0[1]))

    raw_parent_dn = {}
    for riv_id, info in reaches_by_id.items():
        f_idx = info["file_idx"]
        rdict = info["row"]
        dn_id = None

        # Mode A: Pointer (NEXTDOWNID)
        if topo_mode == "pointer" and next_down_col:
            nd_val = rdict.get(next_down_col)
            if nd_val:
                # Handle potential semicolon/comma list by taking primary target
                first_tok = str(nd_val).split(";")[0].split(",")[0].strip()
                dn_id = key_to_rivid.get((f_idx, first_tok)) or global_orig_to_rivid.get(first_tok)

        # Mode B / Fallback: Nodes (Tnode -> Fnode)
        if dn_id is None and topo_mode in ("nodes", "pointer") and tnode_col:
            tn_val = rdict.get(tnode_col)
            if tn_val:
                dn_id = fnode_index.get((f_idx, info["catch_id"], str(tn_val))) or fnode_index.get(
                    ("*", info["catch_id"], str(tn_val))
                )

        # Mode C / Fallback: Geometry endpoint snapping within 5m
        if dn_id is None and topo_mode == "geometry":
            pe = info["poly"][-1]
            gx0, gy0 = int(pe[0] // 50), int(pe[1] // 50)
            best_d = 5.0
            for gx in range(gx0 - 1, gx0 + 2):
                for gy in range(gy0 - 1, gy0 + 2):
                    for cand_id, sx, sy in start_pt_grid.get((f_idx, gx, gy), ()):
                        if cand_id != riv_id:
                            d = math.hypot(pe[0] - sx, pe[1] - sy)
                            if d < best_d:
                                best_d = d
                                dn_id = cand_id

        if dn_id is not None and dn_id != riv_id and dn_id in reaches_by_id:
            # Keep within same CatchID (or allow cross-file link when no explicit CatchID column exists)
            if not has_explicit_catch_id or reaches_by_id[dn_id]["catch_id"] == info["catch_id"]:
                raw_parent_dn[riv_id] = dn_id

    # Break any cycles so the directed graph is a strict forest of trees
    parent_dn = {}
    visited_state = {}  # 0 = visiting, 1 = done
    cycles_broken = 0
    for riv_id in sorted(reaches_by_id.keys()):
        curr = riv_id
        chain = []
        while curr is not None and curr not in visited_state:
            visited_state[curr] = 0
            chain.append(curr)
            nxt = raw_parent_dn.get(curr)
            if nxt is not None and visited_state.get(nxt) == 0:
                # Cycle detected: break link at curr
                raw_parent_dn.pop(curr, None)
                cycles_broken += 1
                break
            curr = nxt
        for u in chain:
            visited_state[u] = 1
            if u in raw_parent_dn:
                parent_dn[u] = raw_parent_dn[u]

    parents = defaultdict(list)
    mouths = []
    for riv_id in sorted(reaches_by_id.keys()):
        dn = parent_dn.get(riv_id)
        if dn is None:
            mouths.append(riv_id)
        else:
            parents[dn].append(riv_id)

    # If the input dataset lacked an explicit CatchID column, assign integer CatchIDs
    # automatically from weakly connected river tree components rooted at each mouth
    if not has_explicit_catch_id and len(mouths) > 1:
        tree_members = {}
        for m_id in mouths:
            comp = []
            stack = [m_id]
            while stack:
                u = stack.pop()
                comp.append(u)
                for p in parents[u]:
                    stack.append(p)
            tree_members[m_id] = comp

        snapped_rivids = {
            key_to_rivid[best_snap[nid][1]]
            for nid in best_snap
            if best_snap[nid][0] <= max_snap_m and best_snap[nid][1] in key_to_rivid
        }
        monitored_mouths = [
            m_id for m_id in mouths if any(u in snapped_rivids for u in tree_members[m_id])
        ]
        active_mouths = monitored_mouths if monitored_mouths else mouths
        active_mouths.sort(key=lambda m_id: (-len(tree_members[m_id]), m_id))

        keep_rivids = set()
        for new_cid, m_id in enumerate(active_mouths, 1):
            for u in tree_members[m_id]:
                reaches_by_id[u]["catch_id"] = new_cid
                keep_rivids.add(u)

        if len(keep_rivids) < len(reaches_by_id):
            reaches_by_id = {u: info for u, info in reaches_by_id.items() if u in keep_rivids}
            parent_dn = {u: dn for u, dn in parent_dn.items() if u in keep_rivids and dn in keep_rivids}
            parents = defaultdict(list)
            for u, dn in parent_dn.items():
                parents[dn].append(u)
            key_to_rivid = {k: v for k, v in key_to_rivid.items() if v in keep_rivids}
        mouths = active_mouths

    # Verify polyline vertex orientation (upstream -> downstream) and compute reach lengths
    flipped_cnt = 0
    for riv_id, info in reaches_by_id.items():
        poly = info["poly"]
        dn = parent_dn.get(riv_id)
        if dn is not None:
            dn_poly = reaches_by_id[dn]["poly"]
            d_end_to_dn = min(
                math.hypot(poly[-1][0] - dn_poly[0][0], poly[-1][1] - dn_poly[0][1]),
                math.hypot(poly[-1][0] - dn_poly[-1][0], poly[-1][1] - dn_poly[-1][1]),
            )
            d_start_to_dn = min(
                math.hypot(poly[0][0] - dn_poly[0][0], poly[0][1] - dn_poly[0][1]),
                math.hypot(poly[0][0] - dn_poly[-1][0], poly[0][1] - dn_poly[-1][1]),
            )
            if d_start_to_dn + 10.0 < d_end_to_dn:
                poly = list(reversed(poly))
                info["poly"] = poly
                if info["zs"]:
                    info["zs"] = list(reversed(info["zs"]))
                flipped_cnt += 1

        tot_m, sinu, coords_net = compute_polyline_metrics(poly)
        info["length_km"] = max(1e-4, tot_m / 1000.0)
        info["sinuosity"] = sinu
        info["coords_net"] = coords_net

    # Downstream-to-upstream tree traversal: topological order & Dist2Mth_km
    topo_order = []
    for m_id in sorted(mouths, key=lambda x: (reaches_by_id[x]["catch_id"], x)):
        stack = [m_id]
        while stack:
            u = stack.pop()
            topo_order.append(u)
            u_info = reaches_by_id[u]
            dn = parent_dn.get(u)
            raw_d2m = _opt_float(u_info["row"].get(d2m_col)) if d2m_col else None
            if raw_d2m is not None and raw_d2m > 0:
                u_info["d2m_km"] = raw_d2m * d2m_scale
            else:
                base_d2m = reaches_by_id[dn]["d2m_km"] if dn is not None else 0.0
                u_info["d2m_km"] = base_d2m + u_info["length_km"]
            for p in parents[u]:
                stack.append(p)

    # Upstream-to-downstream pass: compute Shreve, Strahler, and US_Accum_km where needed
    for u in reversed(topo_order):
        u_info = reaches_by_id[u]
        rdict = u_info["row"]
        ps = parents[u]

        # Shreve order
        raw_shr = _opt_float(rdict.get(shreve_col)) if shreve_col else None
        if raw_shr is not None and int(raw_shr) > 0:
            u_info["shreve"] = int(raw_shr)
        else:
            u_info["shreve"] = sum(reaches_by_id[p]["shreve"] for p in ps) if ps else 1

        # Strahler order
        raw_str = _opt_float(rdict.get(strahler_col)) if strahler_col else None
        if raw_str is not None and int(raw_str) > 0:
            u_info["strahler"] = int(raw_str)
        elif not ps:
            u_info["strahler"] = 1
        else:
            child_s = [reaches_by_id[p]["strahler"] for p in ps]
            mx_s = max(child_s)
            u_info["strahler"] = mx_s + 1 if child_s.count(mx_s) >= 2 else mx_s

        # US_Accum_km (AFV weight)
        raw_afv = _opt_float(rdict.get(afv_col)) if afv_col else None
        if raw_afv is not None and raw_afv > 0:
            u_info["us_accum_km"] = max(1e-4, raw_afv * afv_scale)
        else:
            u_info["us_accum_km"] = u_info["length_km"] + sum(reaches_by_id[p]["us_accum_km"] for p in ps)

    # Top-down pass: compute log_AFV and propagate headwater SrcElevation
    afv_scheme = cfg.get("afv_scheme", "length")
    log_afv = {}
    for m_id in sorted(mouths, key=lambda x: (reaches_by_id[x]["catch_id"], x)):
        log_afv[m_id] = 0.0
        stack = [m_id]
        while stack:
            u = stack.pop()
            ps = parents[u]
            if ps:
                p_weights = {
                    p: compute_afv_reach_weight(
                        afv_scheme,
                        reaches_by_id[p]["us_accum_km"],
                        reaches_by_id[p].get("shreve", 1),
                        reaches_by_id[p].get("elev", 0.0),
                        reaches_by_id[p].get("is_aquifer", 0),
                    )
                    for p in ps
                }
                sum_us = sum(p_weights.values())
                for p in ps:
                    log_afv[p] = log_afv[u] + math.log(p_weights[p] / sum_us)
                    stack.append(p)

    # Transform reach (start, mid, end) coordinates to EPSG:27700 so all maps & Easting/Northing are in BNG
    if crs_epsg != 27700:
        print(f"4. Transforming {len(topo_order) * 3} reach coordinates from EPSG:{crs_epsg} to EPSG:27700...")
        flat_pts = []
        for u in topo_order:
            sx, sy, mx, my, ex, ey = reaches_by_id[u]["coords_net"]
            flat_pts.extend([(sx, sy), (mx, my), (ex, ey)])
        bng_pts = transform_coords_batch(flat_pts, crs_epsg, 27700)
        for idx, u in enumerate(topo_order):
            p_s = bng_pts[3 * idx]
            p_m = bng_pts[3 * idx + 1]
            p_e = bng_pts[3 * idx + 2]
            reaches_by_id[u]["coords_bng"] = (p_s[0], p_s[1], p_m[0], p_m[1], p_e[0], p_e[1])
    else:
        for u in topo_order:
            reaches_by_id[u]["coords_bng"] = reaches_by_id[u]["coords_net"]

    # Extract elevation, slope, and BGS geology (via spatial join with ORN reference grid if available)
    orn_grid, orn_rbd_names = load_orn_spatial_reference(runtime_cfg=eff_cfg)
    for u in topo_order:
        u_info = reaches_by_id[u]
        rdict = u_info["row"]
        f_idx = u_info["file_idx"]
        zs = u_info["zs"]
        _, _, mx, my, _, _ = u_info["coords_bng"]
        orn_ref = query_orn_spatial_grid(orn_grid, mx, my)

        # Elevation: 1. explicit attribute, 2. 3D WKB Z mean, 3. Nodes table FNODE/TNODE, 4. ORN spatial join
        elev = _opt_float(rdict.get(elev_col)) if elev_col else None
        if elev is None and zs:
            elev = statistics.mean(zs)
        elif elev is None and node_lookup:
            fn_key = str(rdict.get(fnode_col))
            tn_key = str(rdict.get(tnode_col))
            fn_el = (node_lookup.get((f_idx, fn_key)) or node_lookup.get(("*", fn_key)) or (None, None))[0]
            tn_el = (node_lookup.get((f_idx, tn_key)) or node_lookup.get(("*", tn_key)) or (None, None))[0]
            vals_el = [x for x in (_opt_float(fn_el), _opt_float(tn_el)) if x is not None]
            if vals_el:
                elev = statistics.mean(vals_el)
        if elev is None:
            elev = orn_ref[3] if orn_ref else 50.0
        u_info["elevation"] = elev

        # Slope (%): 1. explicit attribute (e.g. PENTE), 2. 3D WKB Z drop, 3. ORN spatial join
        raw_slope = _opt_float(rdict.get(slope_col)) if slope_col else None
        slope = raw_slope if (raw_slope is not None and raw_slope >= 0) else None
        if slope is None and len(zs) >= 2 and u_info["length_km"] > 0:
            slope = max(0.0, (zs[0] - zs[-1]) / (u_info["length_km"] * 10.0))
        if slope is None:
            slope = orn_ref[4] if orn_ref else 0.5
        u_info["slope"] = slope

        # Geology & RBD_ID from ORN spatial reference (if in England)
        raw_rbd = _opt_float(rdict.get(rbd_col)) if rbd_col else None
        if orn_ref is not None:
            u_info["rbd_id"] = int(raw_rbd) if raw_rbd is not None else orn_ref[2]
            u_info["src_elevation"] = max(elev, orn_ref[5])
            u_info["is_chalk_lime"] = orn_ref[6]
            u_info["is_sand"] = orn_ref[7]
            u_info["is_aq"] = orn_ref[8]
            u_info["is_peat"] = orn_ref[9]
        else:
            u_info["rbd_id"] = int(raw_rbd) if raw_rbd is not None else 1
            u_info["src_elevation"] = elev
            u_info["is_chalk_lime"] = 0
            u_info["is_sand"] = 0
            u_info["is_aq"] = 0
            u_info["is_peat"] = 0

    # Write station-level snapping summary (rsn_sites.csv) and multi-period reach table (rsn_reach_periods.csv)
    unc_by_site = defaultdict(list)
    meth_by_site = defaultdict(set)
    for rec in imputed_records:
        nid = rec[0]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
        if unc_val is not None:
            unc_by_site[nid].append(unc_val)
        if meth_val:
            meth_by_site[nid].add(meth_val)

    snap_info = {}
    for nid in all_station_ids:
        d, r_key, cid_key, along_m = best_snap[nid]
        if r_key not in key_to_rivid:
            continue
        riv_id = key_to_rivid[r_key]
        u_info = reaches_by_id[riv_id]
        snap_info[nid] = (d, riv_id, u_info["catch_id"], along_m, u_info["d2m_km"], 1 if d <= max_snap_m else 0)

    site_by_reach, sites_csv_path, periods_csv_path = write_snapped_sites_and_periods(
        out_dir=eff_out_dir,
        all_station_ids=all_station_ids,
        pts=pts,
        pts_source=pts_source,
        snap_info=snap_info,
        imputed_records=imputed_records,
        obs_by_site=obs_by_site,
        unc_by_site=unc_by_site,
        meth_by_site=meth_by_site,
        trans_type=trans_type,
        eff_lambda=eff_lambda,
    )

    # Build catchment lookup table (rsn_catchments_lookup.csv)
    catch_reaches = defaultdict(list)
    for u in topo_order:
        catch_reaches[reaches_by_id[u]["catch_id"]].append(reaches_by_id[u])

    sorted_cids = sorted(
        catch_reaches.keys(),
        key=lambda cid: (-sum(1 for r in catch_reaches[cid] if r["riv_id"] in site_by_reach), -len(catch_reaches[cid])),
    )
    catch_rows = []
    for cid in sorted_cids:
        rlist = catch_reaches[cid]
        n_obs_c = sum(1 for r in rlist if r["riv_id"] in site_by_reach)
        rbd_counts = Counter(r["rbd_id"] for r in rlist if r["rbd_id"] > 0)
        rbd_id = rbd_counts.most_common(1)[0][0] if rbd_counts else 1
        rbd_name = orn_rbd_names.get(rbd_id, f"RBD {rbd_id}")

        best_name = ""
        best_score = -1.0
        for r in rlist:
            for nc in name_cols:
                val = r["row"].get(nc)
                if val and str(val).strip():
                    sc = r["us_accum_km"]
                    if sc > best_score:
                        best_score = sc
                        best_name = str(val).strip()
                    break
        if not best_name:
            best_name = f"Catchment {cid}"
        catch_rows.append([cid, best_name, rbd_id, rbd_name, len(rlist), n_obs_c])
    catch_lu_path = write_catchment_lookup(eff_out_dir, catch_rows)

    # Write england_rsn_reaches.csv (including any configured extra_numeric_cols)
    out_reaches_path = os.path.join(eff_out_dir, "england_rsn_reaches.csv")
    print(f"5. Writing {out_reaches_path} ({len(topo_order)} reaches, {len(site_by_reach)} observed)...")
    reaches_header = [
        "RivID", "CatchID", "RBD_ID", "parent_dn",
        "Dist2Mth_km", "up_dist_mid_km", "length_km", "US_Accum_km",
        "log_AFV", "Shreve", "Strahler",
        "start_X", "start_Y", "mid_X", "mid_Y", "end_X", "end_Y",
        "Elevation", "Slope", "Sinuosity", "SrcElevation",
        "Is_Chalk_Limestone", "Is_Sandstone", "Is_Aquifer", "Is_Peat",
    ] + REACH_OBS_HEADER_COLS + extra_num_cols

    reaches_rows = []
    for riv_id in topo_order:
        u_info = reaches_by_id[riv_id]
        d2m_km = u_info["d2m_km"]
        l_km = u_info["length_km"]
        up_mid_km = max(0.0, d2m_km - 0.5 * l_km)
        sx, sy, mx, my, ex, ey = u_info["coords_bng"]
        extra_vals = [round(float(u_info["row"].get(ec) or 0.0), 4) for ec in extra_num_cols]
        base_cols = [
            riv_id, u_info["catch_id"], u_info["rbd_id"], parent_dn.get(riv_id, 0),
            round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(u_info["us_accum_km"], 5),
            round(log_afv[riv_id], 6), u_info["shreve"], u_info["strahler"],
            round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
            round(u_info["elevation"], 3), round(u_info["slope"], 4), round(u_info["sinuosity"], 4), round(u_info["src_elevation"], 3),
            u_info["is_chalk_lime"], u_info["is_sand"], u_info["is_aq"], u_info["is_peat"],
        ]
        if riv_id in site_by_reach:
            obs_cols = summarize_reach_site_obs(site_by_reach[riv_id], d2m_km, trans_type=trans_type, eff_lambda=eff_lambda)
        else:
            obs_cols = [0, "", "", round(up_mid_km, 5), "", "", "", "", "", 0]
        reaches_rows.append(base_cols + obs_cols + extra_vals)

    all_meths = sorted({m for s_list in site_by_reach.values() for s in s_list for m in s[5] if m and m != "exact_point"})
    active_meas_method = "gam" if any("gam" in m for m in all_meths) else (all_meths[0] if all_meths else None)
    meta_obj = {
        "name": cfg.get("name", net_abbr),
        "abbr": net_abbr,
        "label": net_label,
        "short_label": short_label,
        "source": source_key,
        "source_abbr": source_abbr,
        "source_label": source_label,
        "range_mode": range_mode,
        "meas_method": active_meas_method,
        "has_meas_uncertainty": bool(active_meas_method),
        "crs_epsg": crs_epsg,
        "topology_mode": topo_mode,
        "max_snap_dist_m": max_snap_m,
        "n_reaches": len(topo_order),
        "n_obs_reaches": len(site_by_reach),
        "n_catchments": len(catch_reaches),
        "core_basins": cfg.get("core_basins", {}),
        "extra_numeric_cols": extra_num_cols,
        "transform_type": trans_type,
        "transform_lambda": eff_lambda,
    }
    write_reaches_and_finalize(
        out_dir=eff_out_dir,
        net_abbr=net_abbr,
        out_reaches_path=out_reaches_path,
        reaches_header=reaches_header,
        reaches_rows=reaches_rows,
        meta_obj=meta_obj,
        sites_csv_path=sites_csv_path,
        periods_csv_path=periods_csv_path,
        catch_lu_path=catch_lu_path,
        source_key=source_key,
        runtime_cfg=eff_cfg,
    )
    print(f"Done in {time.time() - t0:.2f} seconds (saved to {eff_out_dir}/).")


def main():
    global OUT_DIR
    parser = argparse.ArgumentParser(
        description="Preprocess stream network GeoPackages (ORN v2, EU-Hydro, etc.) and EA RSN samples for S3N modeling."
    )
    parser.add_argument(
        "--gis-dir",
        dest="gis_dir",
        default=None,
        help="Custom GIS directory containing water network GeoPackages (default: 'gis').",
    )
    parser.add_argument(
        "--data-dir",
        dest="data_dir",
        default=None,
        help="Custom data directory containing RSN sample folders and processed outputs (default: 'data').",
    )
    parser.add_argument(
        "--out-dir",
        "--output-dir",
        dest="out_dir",
        default=None,
        help="Custom output directory for processed network CSVs (default: '<processed-dir>/<net>_<source>').",
    )
    parser.add_argument(
        "--network",
        default="orn",
        help="Water network preset ('orn', 'euhydro_thames', 'euhydro_tweed', 'euhydro') or path to GeoPackage (.gpkg), GeoParquet (.parquet/.geoparquet), or JSON config.",
    )
    parser.add_argument(
        "--network-config",
        default=None,
        help="Optional path to a custom network JSON configuration file.",
    )
    parser.add_argument(
        "--convert-to-geoparquet",
        dest="convert_to_geoparquet",
        default=None,
        help="Optional output .parquet/.geoparquet path: converts the selected --network GeoPackage layer(s) to an OGC GeoParquet file and exits.",
    )
    parser.add_argument(
        "--list-networks",
        action="store_true",
        help="List all configured water network profiles and exit.",
    )
    parser.add_argument(
        "--topology-mode",
        choices=["auto", "nodes", "pointer", "geometry"],
        default="auto",
        help="Override how river reach connectivity is constructed ('auto' uses network profile).",
    )
    parser.add_argument(
        "--afv-col",
        default=None,
        help="Override attribute column used for Additive Function Value (AFV) confluence weights.",
    )
    parser.add_argument(
        "--max-snap-dist",
        type=float,
        default=None,
        help="Maximum snapping distance in meters from RSN monitoring stations to stream reaches (default: 250m).",
    )
    parser.add_argument(
        "--censor-method",
        choices=["ros", "ros_site", "half_lod"],
        default="ros",
        help=(
            "Method for handling left-censored ('<') observations below detection limit: "
            "'ros' = Pooled Robust Lognormal Regression on Order Statistics (default), "
            "'ros_site' = Station-level Robust Lognormal ROS with global ROS fallback, "
            "'half_lod' = 0.5 * limit of detection substitution."
        ),
    )
    parser.add_argument(
        "--dataset",
        choices=["wq", "inverts", "diatoms", "macrophytes", "rhs", "lcms"],
        default="wq",
        help=(
            "RSN dataset to extract: 'wq' = Water Quality Samples (default), "
            "'inverts' = Macroinvertebrate Samples, 'diatoms' = Diatom Samples, "
            "'macrophytes' = Macrophyte Samples, 'rhs' = River Habitat Surveys, "
            "'lcms' = LCMS/GCMS Organic Pollutant Samples."
        ),
    )
    parser.add_argument(
        "--det-code",
        default="0117",
        help=(
            "4-digit EA RSN determinand code (for --dataset wq, default: '0117') "
            "or column/metric name (for --dataset inverts/diatoms/macrophytes/rhs, e.g. 'whpt_aspt')."
        ),
    )
    parser.add_argument(
        "--year",
        type=int,
        default=None,
        help="Optional calendar year filter (e.g. 2021, 2022, 2023, 2024). Default: all years.",
    )
    parser.add_argument(
        "--season",
        type=str,
        default=None,
        help="Optional meteorological season filter: winter (DJF), spring (MAM), summer (JJA), autumn/fall (SON).",
    )
    parser.add_argument(
        "--month",
        type=str,
        default=None,
        help="Optional calendar month filter (1..12 or Jan..Dec). Default: all months.",
    )
    parser.add_argument(
        "--transform",
        choices=["log", "log10", "log1p", "sqrt", "cbrt", "fourthrt", "reciprocal", "boxcox", "yeojohnson", "asin", "logit", "none", "raw"],
        default="log",
        help="Data transformation applied to determinand values (default: 'log').",
    )
    parser.add_argument(
        "--transform-lambda",
        dest="transform_lambda",
        type=float,
        default=None,
        help="Power parameter lambda for Box-Cox or Yeo-Johnson transformations (default: auto-estimated via Profile MLE).",
    )
    parser.add_argument(
        "--list-transforms",
        action="store_true",
        help="List all supported data transformations and their mathematical properties, then exit.",
    )
    parser.add_argument(
        "--afv-scheme",
        dest="afv_scheme",
        choices=["length", "shreve", "hydro_yield"],
        default="length",
        help=(
            "Physical confluence weighting scheme for Tail-Up Additive Function Values (log_AFV): "
            "'length' = cumulative upstream channel length US_Accum_km (default), "
            "'shreve' = Shreve hydrological link magnitude, "
            "'hydro_yield' = elevation/orographic + Principal Aquifer baseflow hydrological yield proxy."
        ),
    )
    parser.add_argument(
        "--source",
        default="rsn",
        help=(
            "Monitoring network data source: 'rsn' (default), 'ssn', 'rss' (combined RSN+SSN), "
            "any configured source in configs/sources/<name>.json, data/custom/<name>/, "
            "or delimiter combination (e.g. 'rsn+ssn', 'rsn+custom')."
        ),
    )
    parser.add_argument(
        "--custom-sites",
        default=None,
        help="Path to custom monitoring station points CSV (e.g., with site_id, easting, northing).",
    )
    parser.add_argument(
        "--custom-samples",
        default=None,
        help="Path to custom sample observations CSV (e.g., with site_id, sample_datetime, meas_result).",
    )
    parser.add_argument(
        "--custom-dir",
        default=None,
        help="Path to directory containing custom 'sites.csv' and/or 'samples.csv'.",
    )
    parser.add_argument(
        "--custom-name",
        default="custom",
        help="Short name/tag for the custom dataset (default: 'custom').",
    )
    parser.add_argument(
        "--custom-crs",
        type=int,
        default=27700,
        help="EPSG code for custom station coordinates (default: 27700 BNG; use 4326 for WGS84 lat/lon).",
    )
    parser.add_argument(
        "--range-mode",
        choices=["auto", "tobit", "geometric", "midpoint", "lower", "upper", "none"],
        default="auto",
        help=(
            "Range handling strategy tag to apply to output directory: 'auto' (detect from config/samples, default), "
            "'tobit', 'geometric', 'midpoint', 'lower', 'upper', or 'none'."
        ),
    )
    parser.add_argument(
        "--list-sources",
        action="store_true",
        help="List available monitoring data sources (RSN, SSN, RSS, custom) and their configuration, then exit.",
    )
    parser.add_argument(
        "--update-determinands",
        action="store_true",
        help="Rescan raw RSN/SSN CSVs to add any missing determinands to configs/determinands.json.",
    )
    args = parser.parse_args()
    runtime_cfg = configure_runtime_dirs(gis_dir=args.gis_dir, data_dir=args.data_dir, out_dir=args.out_dir)
    ensure_all_configs(BASE_DIR, data_dir=runtime_cfg.data_dir, update_determinands=args.update_determinands)
    if args.list_transforms:
        print_available_transforms()
        return
    if args.list_sources:
        print_available_sources()
        return

    # Auto-route source to custom-name if custom paths supplied
    if (args.custom_sites or args.custom_dir) and args.source == "rsn":
        args.source = args.custom_name

    sources_reg, composites_reg = load_source_registry(base_dir=BASE_DIR, runtime_cfg=runtime_cfg)
    source_key, constituents, source_abbr, source_label = resolve_source_tag(args.source, (sources_reg, composites_reg), runtime_cfg=runtime_cfg)

    active_range_mode = detect_range_mode(
        source_key=source_key,
        sources_reg=sources_reg,
        custom_samples=args.custom_samples,
        custom_dir=args.custom_dir,
        cli_mode=args.range_mode,
        runtime_cfg=runtime_cfg,
    )

    if active_range_mode:
        # Only append range mode if none of the constituent tags already carry it (e.g. waterblitz_tobit+rsn -> rsn+waterblitz_tobit)
        has_rm_suffix = ("+" in source_key and any(t.endswith(f"_{active_range_mode}") for t in source_key.split("+"))) or source_key.endswith(f"_{active_range_mode}")
        if source_key != "rsn" and not has_rm_suffix:
            # Re-check if an explicit <source>_<mode> alias exists in registry (e.g. waterblitz + --range-mode=midpoint -> waterblitz_midpoint)
            if args.range_mode not in ("auto", "none") and "+" not in source_key:
                base_cand = re.sub(r"_(tobit|geometric|midpoint|lower|upper)$", "", source_key)
                mode_cand = f"{base_cand}_{active_range_mode}"
                if mode_cand in sources_reg or os.path.isdir(_rel_join(DATA_DIR, "custom", mode_cand)):
                    source_key, constituents, source_abbr, source_label = resolve_source_tag(mode_cand, (sources_reg, composites_reg))
                else:
                    source_key = mode_cand
            else:
                source_key = f"{source_key}_{active_range_mode}"
            if f"_{active_range_mode.upper()[:3]}" not in source_abbr and f"[{active_range_mode.title()}]" not in source_abbr:
                source_abbr = f"{source_abbr}_{active_range_mode.upper()[:3]}"
            if f"[{active_range_mode.title()}]" not in source_label:
                source_label = f"{source_label} [{active_range_mode.title()}]"
        print(f"   [Range Handling] Detected range mode: {active_range_mode.upper()} (Processed folder tagged with _{active_range_mode})")

    if args.list_networks:
        list_available_networks()
        return

    cfg = load_network_config(args.network, args.network_config)
    cfg["afv_scheme"] = getattr(args, "afv_scheme", "length") or "length"
    net_abbr = cfg["abbr"]

    if not args.out_dir:
        sub_name = f"{net_abbr}_{source_key}"
        out_dir = _rel_join(runtime_cfg.processed_dir, sub_name)
        OUT_DIR = out_dir
    else:
        out_dir = runtime_cfg.processed_dir
    os.makedirs(out_dir, exist_ok=True)

    if args.convert_to_geoparquet:
        src_files = ensure_gpkg_paths(cfg)
        layers = cfg.get("layers")
        print(f"Converting {src_files[0]} (layers={layers}) to GeoParquet: {args.convert_to_geoparquet}...")
        convert_gpkg_to_geoparquet(src_files[0], args.convert_to_geoparquet, layers=layers)
        print(f"Successfully wrote GeoParquet file: {args.convert_to_geoparquet}")
        return

    season_filter = None
    if args.season is not None:
        s_low = args.season.strip().lower()
        if s_low == "fall":
            s_low = "autumn"
        if s_low not in ("winter", "spring", "summer", "autumn"):
            raise ValueError(f"Invalid --season '{args.season}'. Choose from: winter, spring, summer, autumn.")
        season_filter = s_low
    month_filter = parse_month_arg(args.month)

    period_label = "all"
    if args.year or season_filter or month_filter:
        parts = []
        if args.year:
            parts.append(f"year={args.year}")
        if season_filter:
            parts.append(f"season={season_filter}")
        if month_filter:
            parts.append(f"month={month_filter:02d}")
        period_label = ",".join(parts)

    t0 = time.time()
    print(
        f"1. Reading Monitoring Sample Points and Samples [Source: {source_label}] "
        f"(network='{net_abbr}', dataset='{args.dataset}', variable='{args.det_code}', "
        f"censor_method='{args.censor_method}', period='{period_label}')..."
    )
    pts, pts_source, core_wq_sites, all_records = load_monitoring_stations_and_records(
        DATA_DIR,
        source=source_key,
        dataset=args.dataset,
        det_code=args.det_code,
        custom_sites=args.custom_sites,
        custom_samples=args.custom_samples,
        custom_dir=args.custom_dir,
        custom_name=args.custom_name,
        custom_crs=args.custom_crs,
    )

    raw_all_by_site = defaultdict(list)
    for rec in all_records:
        nid, val, is_cens = rec[0], rec[1], rec[2]
        meth_val = rec[7] if len(rec) > 7 else ""
        # Records already imputed via Tobit conditional expectation E[Y | Y <= b] must not be double-imputed
        needs_cens_imp = bool(is_cens) and ("tobit" not in meth_val)
        raw_all_by_site[nid].append((val, needs_cens_imp))
    imp_all_by_site = impute_censored_by_site(raw_all_by_site, method=args.censor_method)

    site_cursors = defaultdict(int)
    imputed_records = []
    for rec in all_records:
        nid, val, is_cens, yr, mo, seas = rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
        idx = site_cursors[nid]
        site_cursors[nid] += 1
        imp_val = imp_all_by_site[nid][idx]
        imputed_records.append((nid, imp_val, yr, mo, seas, is_cens, unc_val, meth_val))

    obs_by_site = defaultdict(list)
    unc_by_site = defaultdict(list)
    meth_by_site = defaultdict(set)
    n_slice_samples = 0
    n_slice_cens = 0
    for rec in imputed_records:
        nid, imp_val, yr, mo, seas, is_cens = rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
        if args.year is not None and yr != args.year:
            continue
        if season_filter is not None and seas != season_filter:
            continue
        if month_filter is not None and mo != month_filter:
            continue
        obs_by_site[nid].append(imp_val)
        if unc_val is not None:
            unc_by_site[nid].append(unc_val)
        if meth_val:
            meth_by_site[nid].add(meth_val)
        n_slice_samples += 1
        if is_cens:
            n_slice_cens += 1

    print(
        f"   Found {len(obs_by_site)} active stations in period '{period_label}' "
        f"({n_slice_samples} samples; {n_slice_cens} censored '<' "
        f"[{100.0 * n_slice_cens / max(1, n_slice_samples):.2f}%] imputed via '{args.censor_method}')."
    )

    # Resolve transformation type and Profile Maximum Likelihood lambda
    trans_type = (args.transform or "log").lower().strip()
    all_sample_vals = [v for s_vals in obs_by_site.values() for v in s_vals if v is not None]
    if trans_type in ("boxcox", "bc"):
        eff_lambda = args.transform_lambda if args.transform_lambda is not None else estimate_boxcox_lambda(all_sample_vals)
        print(f"   Estimated Profile Maximum Likelihood lambda = {eff_lambda:.4f} for Box-Cox transformation.")
    elif trans_type in ("yeojohnson", "yj"):
        eff_lambda = args.transform_lambda if args.transform_lambda is not None else estimate_yeojohnson_lambda(all_sample_vals)
        print(f"   Estimated Profile Maximum Likelihood lambda = {eff_lambda:.4f} for Yeo-Johnson transformation.")
    else:
        eff_lambda = args.transform_lambda

    # Dispatch to generic network pipeline if not using standard ORN fast path
    if cfg.get("covariates", {}).get("mode") != "orn_100m_table" or args.topology_mode != "auto" or args.afv_col:
        preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, obs_by_site, period_label, t0,
                                   trans_type=trans_type, eff_lambda=eff_lambda,
                                   pts_source=pts_source, source_key=source_key,
                                   range_mode=active_range_mode,
                                   runtime_cfg=runtime_cfg, out_dir=out_dir)
        return

    # Standard ORN v2 fast path
    max_snap_m = float(args.max_snap_dist if args.max_snap_dist is not None else cfg.get("max_snap_dist_m", 250.0))
    gpkg_path = ensure_gpkg_paths(cfg)[0]
    print(f"2. Scanning {os.path.basename(gpkg_path)} geometry to snap {source_abbr} stations (<= {max_snap_m:.0f}m)...")
    conn = sqlite3.connect(gpkg_path)
    c = conn.cursor()

    all_station_ids = sorted(pts.keys())
    site_grid = defaultdict(list)
    for nid in all_station_ids:
        px, py, _ = pts[nid]
        site_grid[(int(px // 2000), int(py // 2000))].append((nid, px, py))

    best_snap = {nid: (1e18, None, None, 0.0, 0.0) for nid in all_station_ids}

    for row in c.execute("SELECT RivID, CatchID, Dist2Mth, Shape FROM ORN"):
        riv_id, catch_id, d2m, blob = row
        env = parse_gpkg_envelope(blob)
        if env is not None:
            minx, maxx, miny, maxy = env
        else:
            poly = parse_gpkg_multilinestring(blob)
            if not poly:
                continue
            xs = [p[0] for p in poly]
            ys = [p[1] for p in poly]
            minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        snap_pad_orn = max(1500.0, max_snap_m + 500.0)
        gx0, gx1 = int((minx - snap_pad_orn) // 2000), int((maxx + snap_pad_orn) // 2000)
        gy0, gy1 = int((miny - snap_pad_orn) // 2000), int((maxy + snap_pad_orn) // 2000)
        cand_sites = []
        for gx in range(gx0, gx1 + 1):
            for gy in range(gy0, gy1 + 1):
                if (gx, gy) in site_grid:
                    cand_sites.extend(site_grid[(gx, gy)])
        if not cand_sites:
            continue

        poly = parse_gpkg_multilinestring(blob)
        for nid, px, py in cand_sites:
            d, along = point_to_polyline(px, py, poly)
            if d < best_snap[nid][0]:
                best_snap[nid] = (d, riv_id, catch_id, along, d2m)

    active_monitored_sites = set(obs_by_site.keys()) | core_wq_sites
    monitored_catchments = {
        best_snap[nid][2]
        for nid in active_monitored_sites
        if nid in best_snap and best_snap[nid][0] <= max_snap_m
    }
    snap_info = {}
    for nid in all_station_ids:
        d, riv_id, catch_id, along, d2m = best_snap[nid]
        if catch_id is None or (catch_id not in monitored_catchments and nid not in core_wq_sites):
            continue
        inc = 1 if (d <= max_snap_m and catch_id in monitored_catchments) else 0
        snap_info[nid] = (d, riv_id, catch_id, along, d2m / 1000.0, inc)

    site_by_reach, sites_csv_path, periods_csv_path = write_snapped_sites_and_periods(
        out_dir=out_dir,
        all_station_ids=all_station_ids,
        pts=pts,
        pts_source=pts_source,
        snap_info=snap_info,
        imputed_records=imputed_records,
        obs_by_site=obs_by_site,
        unc_by_site=unc_by_site,
        meth_by_site=meth_by_site,
        trans_type=trans_type,
        eff_lambda=eff_lambda,
    )

    print(
        f"   Snapped {sum(1 for nid in best_snap if best_snap[nid][0] <= max_snap_m)} sites (<={max_snap_m:.0f}m) "
        f"({len(site_by_reach)} active reaches in period '{period_label}') across {len(monitored_catchments)} catchments."
    )

    print("3. Loading ORN reaches and ORN_Sampled_100m covariates for all monitored catchments...")
    bedrock_lu = {r[1]: (r[5] or "") for r in c.execute("SELECT * FROM tblBedrockGeologyLookUp")}
    superficial_lu = {r[1]: (r[5] or "") for r in c.execute("SELECT * FROM tblSuperficialGeologyLookUp")}
    rbd_lu = {r[1]: (r[2] or "") for r in c.execute("SELECT * FROM tblRiverBasinLookUp")}

    catch_list_sql = ",".join(str(cid) for cid in sorted(monitored_catchments))

    reach_100m = {}
    for r in c.execute(f"""
        SELECT PolylineID,
               AVG(OS_T50_Ele),
               AVG(Slope),
               AVG(Sinuosity),
               MAX(LEX_RCS_BEDROCK),
               MAX(LEX_RCS_SUPERFICIAL),
               AVG(OS_T50_SrcEle)
        FROM ORN_Sampled_100m
        WHERE CatchID IN ({catch_list_sql})
        GROUP BY PolylineID
    """):
        reach_100m[r[0]] = r[1:]

    reaches = {}
    reach_len = {}
    reach_coords = {}
    fnode_map = {}
    catch_main_river = {}
    catch_reach_cnt = defaultdict(int)

    for row in c.execute(f"""
        SELECT RivID, Fnode, Tnode, CatchID, RBD_ID, Dist2Mth, US_Accum, Shreve, Strahler, Hack, Shape, name1
        FROM ORN
        WHERE CatchID IN ({catch_list_sql})
    """):
        riv_id = row[0]
        catch_id = row[3]
        rbd_id = row[4]
        shreve = row[7] or 0
        nm1 = (row[11] or "").strip()
        catch_reach_cnt[catch_id] += 1
        if nm1 and (catch_id not in catch_main_river or shreve > catch_main_river[catch_id][0]):
            catch_main_river[catch_id] = (shreve, nm1, rbd_id)
        reaches[riv_id] = row[:10]
        fnode_map[(catch_id, row[1])] = riv_id
        poly = parse_gpkg_multilinestring(row[10])
        tot, _, coords_tup = compute_polyline_metrics(poly)
        reach_len[riv_id] = tot
        reach_coords[riv_id] = coords_tup

    catch_obs_cnt = defaultdict(int)
    for riv_id in site_by_reach:
        if riv_id in reaches:
            catch_obs_cnt[reaches[riv_id][3]] += 1

    catch_rows = []
    for cid in sorted(monitored_catchments, key=lambda k: (-catch_obs_cnt[k], -catch_reach_cnt[k])):
        _, nm1, rbd_id = catch_main_river.get(cid, (0, f"Catchment {cid}", 0))
        catch_rows.append([cid, nm1, rbd_id, rbd_lu.get(rbd_id, ""), catch_reach_cnt[cid], catch_obs_cnt[cid]])
    catch_lu_path = write_catchment_lookup(out_dir, catch_rows)

    print(f"   Loaded {len(reaches)} reaches across {len(monitored_catchments)} catchments.")

    print("4. Building directed river trees, computing log_AFV, and imputing covariates for short reaches...")
    parents = defaultdict(list)
    mouths = []
    for riv_id, r in reaches.items():
        catch_id = r[3]
        dn = fnode_map.get((catch_id, r[2]))
        if dn is None:
            mouths.append(riv_id)
        else:
            parents[dn].append(riv_id)

    parent_dn = {}
    log_afv = {}
    topo_order = []
    visited = set()
    for m_id in sorted(mouths, key=lambda x: reaches[x][3]):
        if m_id in visited:
            continue
        visited.add(m_id)
        log_afv[m_id] = 0.0
        stack = [m_id]
        while stack:
            u = stack.pop()
            topo_order.append(u)
            ps = [p for p in parents[u] if p not in visited]
            if ps:
                for p in ps:
                    visited.add(p)
                    parent_dn[p] = u
                    stack.append(p)

    for _ in range(3):
        for u in topo_order:
            if u not in reach_100m or reach_100m[u][0] is None:
                dn = parent_dn.get(u)
                if dn and dn in reach_100m and reach_100m[dn][0] is not None:
                    reach_100m[u] = reach_100m[dn]
        for u in reversed(topo_order):
            if u not in reach_100m or reach_100m[u][0] is None:
                for p in parents[u]:
                    if p in reach_100m and reach_100m[p][0] is not None:
                        reach_100m[u] = reach_100m[p]
                        break

    def _classify_bedrock(b_desc):
        is_cl = 1 if ("CHALK" in b_desc or "LIMESTONE" in b_desc or "DOLOSTONE" in b_desc) else 0
        is_sd = 1 if ("SANDSTONE" in b_desc and "MUDSTONE" not in b_desc and "SILTSTONE" not in b_desc) else 0
        return is_cl, is_sd, (1 if (is_cl or is_sd) else 0)

    # Compute log_AFV using selected physical confluence weighting scheme (--afv-scheme)
    afv_scheme = cfg.get("afv_scheme", "length")
    for u in topo_order:
        ps = [p for p in parents[u] if parent_dn.get(p) == u]
        if ps:
            p_weights = {}
            for p in ps:
                r100 = reach_100m.get(p)
                elev_p = (r100[0] or 0.0) if r100 else 0.0
                b_desc_p = bedrock_lu.get(r100[3], "") if (r100 and len(r100) > 3) else ""
                _, _, aq_p = _classify_bedrock(b_desc_p)
                p_weights[p] = compute_afv_reach_weight(
                    afv_scheme,
                    reaches[p][6] or 0.0,
                    reaches[p][7] or 1,
                    elev_p,
                    aq_p,
                )
            sum_us = sum(p_weights.values())
            for p in ps:
                log_afv[p] = log_afv[u] + math.log(p_weights[p] / sum_us)

    missing_cnt = sum(1 for u in topo_order if u not in reach_100m or reach_100m[u][0] is None)
    print(f"   Tree traversal complete ({len(topo_order)} reaches, {missing_cnt} missing covariates).")

    out_reaches_path = os.path.join(out_dir, "england_rsn_reaches.csv")
    print(f"5. Writing {out_reaches_path}...")
    reaches_header = [
        "RivID", "CatchID", "RBD_ID", "parent_dn",
        "Dist2Mth_km", "up_dist_mid_km", "length_km", "US_Accum_km",
        "log_AFV", "Shreve", "Strahler",
        "start_X", "start_Y", "mid_X", "mid_Y", "end_X", "end_Y",
        "Elevation", "Slope", "Sinuosity", "SrcElevation",
        "Is_Chalk_Limestone", "Is_Sandstone", "Is_Aquifer", "Is_Peat",
    ] + REACH_OBS_HEADER_COLS

    reaches_rows = []
    for riv_id in topo_order:
        r = reaches[riv_id]
        rc = reach_100m.get(riv_id, (50.0, 0.5, 1.05, 0, 0, 100.0))
        b_desc = bedrock_lu.get(rc[3], "")
        s_desc = superficial_lu.get(rc[4], "")
        is_chalk_lime, is_sand, is_aq = _classify_bedrock(b_desc)
        is_peat = 1 if "PEAT" in s_desc else 0

        d2m_km = r[5] / 1000.0
        l_km = reach_len[riv_id] / 1000.0
        up_mid_km = max(0.0, d2m_km - 0.5 * l_km)
        sx, sy, mx, my, ex, ey = reach_coords[riv_id]

        base_cols = [
            riv_id, r[3], r[4], parent_dn.get(riv_id, 0),
            round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(r[6] / 1000.0, 5),
            round(log_afv[riv_id], 6), r[7], r[8],
            round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
            round(rc[0] or 0.0, 3), round(rc[1] or 0.0, 4), round(rc[2] or 1.0, 4), round(rc[5] or rc[0] or 0.0, 3),
            is_chalk_lime, is_sand, is_aq, is_peat,
        ]
        if riv_id in site_by_reach:
            obs_cols = summarize_reach_site_obs(site_by_reach[riv_id], d2m_km, trans_type=trans_type, eff_lambda=eff_lambda)
        else:
            obs_cols = [0, "", "", round(up_mid_km, 5), "", "", "", "", "", 0]
        reaches_rows.append(base_cols + obs_cols)

    all_meths = sorted({m for s_list in site_by_reach.values() for s in s_list for m in s[5] if m and m != "exact_point"})
    active_meas_method = "gam" if any("gam" in m for m in all_meths) else (all_meths[0] if all_meths else None)
    meta_obj = {
        "name": cfg.get("name", "orn"),
        "abbr": net_abbr,
        "label": cfg.get("label", "OS Open Rivers Network (ORN v2)"),
        "short_label": cfg.get("short_label", "ORN v2"),
        "source": source_key,
        "source_abbr": source_abbr,
        "source_label": source_label,
        "range_mode": active_range_mode,
        "meas_method": active_meas_method,
        "has_meas_uncertainty": bool(active_meas_method),
        "crs_epsg": 27700,
        "topology_mode": "nodes",
        "max_snap_dist_m": max_snap_m,
        "n_reaches": len(topo_order),
        "n_obs_reaches": len(site_by_reach),
        "n_catchments": len(monitored_catchments),
        "core_basins": cfg.get("core_basins", {"humber": 80, "thames": 3055, "severn": 2786}),
        "extra_numeric_cols": [],
        "transform_type": trans_type,
        "transform_lambda": eff_lambda,
    }
    write_reaches_and_finalize(
        out_dir=out_dir,
        net_abbr=net_abbr,
        out_reaches_path=out_reaches_path,
        reaches_header=reaches_header,
        reaches_rows=reaches_rows,
        meta_obj=meta_obj,
        sites_csv_path=sites_csv_path,
        periods_csv_path=periods_csv_path,
        catch_lu_path=catch_lu_path,
        source_key=source_key,
        runtime_cfg=runtime_cfg,
    )
    print(f"Done in {time.time() - t0:.2f} seconds (saved to {out_dir}/).")


if __name__ == "__main__":
    main()
