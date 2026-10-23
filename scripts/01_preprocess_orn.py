#!/usr/bin/env python3
"""
01_preprocess_orn.py

Extracts and preprocesses the Ordnance Survey Open Rivers Network (ORN v2)
and Environment Agency River Surveillance Network (RSN) water quality
dataset for Scalable Spatial Stream Network (S3N) modeling.

Outputs:
  - data/processed/england_rsn_reaches.csv : All 90,876 reaches across the 141
    monitored river catchments in England, with topology (parent_dn, CatchID),
    upstream distances (Dist2Mth_km, up_dist_mid_km, obs_up_dist_km), additive
    function values (log_AFV), environmental/geological covariates, and observed
    water quality summaries.
  - data/processed/rsn_sites.csv : Station-level summary of all RSN
    monitoring sites and their snapping distances to ORN reaches.
  - data/processed/rsn_reach_periods.csv : Multi-period reach summaries across
    all data, by year, by season, and by month.
"""

import argparse
import csv
import math
import os
import sqlite3
import statistics
import struct
import time
from collections import defaultdict

BASE_DIR = "/home/ben/code/r/spatial_stream_model"
WQ_DIR = os.path.join(BASE_DIR, "data/NCEA_River_Surveillance_Network_ARD_Open")
GPKG_PATH = os.path.join(BASE_DIR, "gis/ORN_v2.gpkg")
OUT_DIR = os.path.join(BASE_DIR, "data/processed")
os.makedirs(OUT_DIR, exist_ok=True)

_STD_NORM = statistics.NormalDist()


def compute_ros_plotting_positions(records):
    """
    Computes Weibull / Kaplan-Meier plotting positions P_i in (0, 1) for a list
    of (value, is_censored) tuples with single or multiple censoring thresholds,
    following Hirsch & Stedinger (1987) and Helsel & Cohn (1988).
    """
    n = len(records)
    if n == 0:
        return []
    c_limits = sorted(set(v for v, is_c in records if is_c))
    if not c_limits:
        order = sorted(range(n), key=lambda i: records[i][0])
        pp = [0.0] * n
        for r_idx, i in enumerate(order, 1):
            pp[i] = r_idx / (n + 1.0)
        return pp

    thresholds = [0.0] + c_limits + [float("inf")]
    J = len(c_limits)
    A = [[] for _ in range(J + 1)]
    C = [[] for _ in range(J + 1)]
    limit_to_j = {lim: idx + 1 for idx, lim in enumerate(c_limits)}

    for idx, (v, is_c) in enumerate(records):
        if is_c:
            C[limit_to_j[v]].append(idx)
        else:
            placed = False
            for j in range(J, -1, -1):
                if v >= thresholds[j]:
                    A[j].append(idx)
                    placed = True
                    break
            if not placed:
                A[0].append(idx)

    for j in range(J + 1):
        A[j].sort(key=lambda i: records[i][0])

    # Exceedance probabilities p[j] = P(Y >= d_j) from highest threshold down to lowest
    p = [0.0] * (J + 2)
    p[J + 1] = 0.0
    for j in range(J, 0, -1):
        d_j = thresholds[j]
        A_j = len(A[j])
        B_j = sum(1 for v, is_c in records if (is_c and v <= d_j) or (not is_c and v < d_j))
        if A_j + B_j > 0:
            p[j] = p[j + 1] + (1.0 - p[j + 1]) * (A_j / (A_j + B_j))
        else:
            p[j] = p[j + 1]
    p[0] = 1.0 if len(A[0]) > 0 else p[1]

    pp = [0.0] * n
    for j in range(J + 1):
        A_j = len(A[j])
        for r_idx, i in enumerate(A[j], 1):
            pp[i] = (1.0 - p[j]) + (p[j] - p[j + 1]) * (r_idx / (A_j + 1.0))
    for j in range(1, J + 1):
        C_j = len(C[j])
        for k_idx, i in enumerate(C[j], 1):
            pp[i] = (1.0 - p[j]) * (k_idx / (C_j + 1.0))

    return [min(1.0 - 1e-6, max(1e-6, x)) for x in pp]


def fit_ros_ols(records, pp):
    """
    Fits OLS linear regression of log(y_uncensored) on normal quantiles z = Phi^{-1}(P_i).
    Returns (intercept_mu, slope_sigma) or None if fewer than 3 uncensored values exist.
    """
    zs = [_STD_NORM.inv_cdf(pp[i]) for i, (v, is_c) in enumerate(records) if not is_c and v > 0]
    ly = [math.log(v) for v, is_c in records if not is_c and v > 0]
    if len(zs) < 3:
        return None
    mz = statistics.mean(zs)
    ml = statistics.mean(ly)
    vz = sum((z - mz) ** 2 for z in zs)
    if vz < 1e-9:
        return None
    cov_zl = sum((z - mz) * (l - ml) for z, l in zip(zs, ly))
    slope = max(0.05, cov_zl / vz)
    intercept = ml - slope * mz
    return (intercept, slope)


def apply_ros_to_records(records, fallback_params=None):
    """
    Applies Robust Lognormal ROS (Helsel & Cohn 1988) to a list of (value, is_censored)
    tuples: uncensored detections keep their observed values, while censored ('<')
    observations are imputed via exp(mu_hat + sigma_hat * Phi^{-1}(P_i)), capped at LOD.
    """
    if not any(is_c for _, is_c in records):
        return [v for v, _ in records]
    pp = compute_ros_plotting_positions(records)
    n_unc = sum(1 for v, is_c in records if not is_c and v > 0)
    params = None
    if n_unc >= 3 and (n_unc / len(records)) >= 0.20:
        params = fit_ros_ols(records, pp)
    if params is None:
        params = fallback_params
    if params is None:
        # Fallback to 0.5 * LOD if neither local nor global ROS fit is available
        return [0.5 * v if is_c else v for v, is_c in records]

    mu_hat, sigma_hat = params
    out = []
    for i, (v, is_c) in enumerate(records):
        if is_c:
            z_i = _STD_NORM.inv_cdf(pp[i])
            imp = math.exp(mu_hat + sigma_hat * z_i)
            out.append(min(v, max(1e-6, imp)))
        else:
            out.append(v)
    return out


def impute_censored_by_site(raw_by_site, method="half_lod"):
    """
    Imputes left-censored ('<') water quality records across stations using:
      - 'half_lod' : 0.5 * detection limit substitution
      - 'ros'      : Pooled Robust Lognormal Regression on Order Statistics (ROS)
      - 'ros_site' : Station-level Robust Lognormal ROS with global ROS fallback
    """
    if method == "half_lod":
        return {
            nid: [0.5 * v if is_c else v for v, is_c in recs]
            for nid, recs in raw_by_site.items()
        }

    # Build global ROS model across all stations
    flat_nids = []
    flat_recs = []
    for nid, recs in sorted(raw_by_site.items()):
        for rec in recs:
            flat_nids.append(nid)
            flat_recs.append(rec)

    pp_glob = compute_ros_plotting_positions(flat_recs)
    glob_params = fit_ros_ols(flat_recs, pp_glob)

    if method == "ros":
        flat_imp = apply_ros_to_records(flat_recs, fallback_params=glob_params)
        out = defaultdict(list)
        for nid, val in zip(flat_nids, flat_imp):
            out[nid].append(val)
        return dict(out)
    elif method == "ros_site":
        return {
            nid: apply_ros_to_records(recs, fallback_params=glob_params)
            for nid, recs in raw_by_site.items()
        }
    else:
        raise ValueError(f"Unknown censor_method: {method}")


def parse_gpkg_multilinestring(blob):
    """Fast parser for GeoPackage WKB MultiLineString / LineString in EPSG:27700."""
    env_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
    flags = blob[3]
    env_ind = (flags >> 1) & 0x07
    wkb_offset = 8 + env_sizes[env_ind]
    wkb_type = struct.unpack_from("<I", blob, wkb_offset + 1)[0]
    pts = []
    if wkb_type == 5:  # MultiLineString
        num_geoms = struct.unpack_from("<I", blob, wkb_offset + 5)[0]
        offset = wkb_offset + 9
        for _ in range(num_geoms):
            npts = struct.unpack_from("<I", blob, offset + 5)[0]
            offset += 9
            coords = struct.unpack_from(f"<{2 * npts}d", blob, offset)
            offset += 16 * npts
            for i in range(0, 2 * npts, 2):
                pts.append((coords[i], coords[i + 1]))
    elif wkb_type == 2:  # LineString
        npts = struct.unpack_from("<I", blob, wkb_offset + 5)[0]
        offset = wkb_offset + 9
        coords = struct.unpack_from(f"<{2 * npts}d", blob, offset)
        for i in range(0, 2 * npts, 2):
            pts.append((coords[i], coords[i + 1]))
    return pts


def point_to_polyline(px, py, poly):
    """Returns (min_dist_m, distance_along_from_upstream_m) from point (px, py) to poly."""
    best_d2 = 1e30
    best_along = 0.0
    cum_len = 0.0
    for i in range(len(poly) - 1):
        x1, y1 = poly[i]
        x2, y2 = poly[i + 1]
        dx, dy = x2 - x1, y2 - y1
        seg_len2 = dx * dx + dy * dy
        seg_len = math.sqrt(seg_len2)
        if seg_len2 == 0.0:
            t = 0.0
        else:
            t = max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / seg_len2))
        projx = x1 + t * dx
        projy = y1 + t * dy
        d2 = (px - projx) ** 2 + (py - projy) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_along = cum_len + t * seg_len
        cum_len += seg_len
    return math.sqrt(best_d2), best_along


MONTH_TO_SEASON = {
    12: "winter", 1: "winter", 2: "winter",
    3: "spring",  4: "spring", 5: "spring",
    6: "summer",  7: "summer", 8: "summer",
    9: "autumn", 10: "autumn", 11: "autumn",
}

MONTH_NAME_MAP = {
    "jan": 1, "january": 1, "feb": 2, "february": 2, "mar": 3, "march": 3,
    "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6,
    "jul": 7, "july": 7, "aug": 8, "august": 8, "sep": 9, "september": 9,
    "oct": 10, "october": 10, "nov": 11, "november": 11, "dec": 12, "december": 12,
}


def parse_month_arg(m_str):
    if m_str is None:
        return None
    s = str(m_str).strip().lower()
    if s.isdigit():
        m = int(s)
        if 1 <= m <= 12:
            return m
    if s in MONTH_NAME_MAP:
        return MONTH_NAME_MAP[s]
    raise ValueError(f"Invalid --month '{m_str}'. Use 1..12 or Jan..Dec.")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess ORN_v2.gpkg and EA RSN Water Quality Samples for S3N modeling."
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
    args = parser.parse_args()

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
        f"1. Reading RSN Sample Points and Samples "
        f"(dataset='{args.dataset}', variable='{args.det_code}', censor_method='{args.censor_method}', period='{period_label}')..."
    )
    pts = {}
    with open(os.path.join(WQ_DIR, "RSN Water Quality Sample Points.csv"), "r", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["ngr_easting"] and r["ngr_northing"]:
                pts[r["network_id"]] = (float(r["ngr_easting"]), float(r["ngr_northing"]), r)
    eco_pts_path = os.path.join(WQ_DIR, "RSN Ecology Sample Points.csv")
    if os.path.exists(eco_pts_path):
        with open(eco_pts_path, "r", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                nid = r["network_id"]
                if nid not in pts and r["ngr_easting"] and r["ngr_northing"]:
                    r_compat = {"wims_id": nid, "wims_region_code": ""}
                    pts[nid] = (float(r["ngr_easting"]), float(r["ngr_northing"]), r_compat)

    # Always determine the 141 core monitored catchments from Nitrate (0117) stations so the 90,876-reach network is invariant
    core_wq_sites = set()
    with open(os.path.join(WQ_DIR, "RSN Water Quality Samples.csv"), "r", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["det_code"] == "0117":
                core_wq_sites.add(r["network_id"])

    # Read all samples for the requested dataset & variable along with their timestamps
    all_records = []  # list of (network_id, val, is_cens, year, month, season)
    if args.dataset == "wq":
        with open(os.path.join(WQ_DIR, "RSN Water Quality Samples.csv"), "r", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                if r["det_code"] == args.det_code:
                    val = float(r["meas_result"])
                    is_cens = (r["meas_sign"] == "<")
                    dt = r["sample_datetime"]
                    yr = int(dt[0:4])
                    mo = int(dt[5:7])
                    seas = MONTH_TO_SEASON[mo]
                    all_records.append((r["network_id"], val, is_cens, yr, mo, seas))
    elif args.dataset == "lcms":
        with open(os.path.join(WQ_DIR, "RSN LCMS GCMS Samples.csv"), "r", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                if r["cas_number"] == args.det_code or args.det_code.lower() in r["compound_name"].lower():
                    if not r["meas_result"]:
                        continue
                    val = float(r["meas_result"])
                    is_cens = (r["less_than_flag"].strip().lower() == "yes")
                    dt = r["sample_datetime"]
                    yr = int(dt[0:4])
                    mo = int(dt[5:7])
                    seas = MONTH_TO_SEASON[mo]
                    all_records.append((r["network_id"], val, is_cens, yr, mo, seas))
    else:
        ds_files = {
            "inverts": ("RSN Macroinvertebrate Samples.csv", "sample_datetime"),
            "diatoms": ("RSN Diatom Samples.csv", "sample_datetime"),
            "macrophytes": ("RSN Macrophyte Samples.csv", "sample_datetime"),
            "rhs": ("RSN River Habitat Surveys.csv", "survey_date"),
        }
        csv_name, dt_col = ds_files[args.dataset]
        col_key = "whpt_aspt" if (args.dataset == "inverts" and args.det_code == "0117") else args.det_code
        with open(os.path.join(WQ_DIR, csv_name), "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            if col_key not in reader.fieldnames:
                raise ValueError(f"Column '{col_key}' not found in {csv_name}. Available: {reader.fieldnames}")
            for r in reader:
                raw_s = r[col_key].strip()
                dt = r[dt_col].strip()
                if not raw_s or not dt:
                    continue
                val = float(raw_s)
                yr = int(dt[0:4])
                mo = int(dt[5:7])
                seas = MONTH_TO_SEASON[mo]
                all_records.append((r["network_id"], val, False, yr, mo, seas))

    # Impute censored values across the full time series first so every record has a consistent value
    raw_all_by_site = defaultdict(list)
    for nid, val, is_cens, yr, mo, seas in all_records:
        raw_all_by_site[nid].append((val, is_cens))
    imp_all_by_site = impute_censored_by_site(raw_all_by_site, method=args.censor_method)

    # Re-attach imputed values to (nid, imp_val, yr, mo, seas, is_cens)
    site_cursors = defaultdict(int)
    imputed_records = []
    for nid, val, is_cens, yr, mo, seas in all_records:
        idx = site_cursors[nid]
        site_cursors[nid] += 1
        imp_val = imp_all_by_site[nid][idx]
        imputed_records.append((nid, imp_val, yr, mo, seas, is_cens))

    # Filter to active temporal slice (--year, --season, --month)
    nitrate_by_site = defaultdict(list)
    n_slice_samples = 0
    n_slice_cens = 0
    for nid, imp_val, yr, mo, seas, is_cens in imputed_records:
        if args.year is not None and yr != args.year:
            continue
        if season_filter is not None and seas != season_filter:
            continue
        if month_filter is not None and mo != month_filter:
            continue
        nitrate_by_site[nid].append(imp_val)
        n_slice_samples += 1
        if is_cens:
            n_slice_cens += 1

    print(
        f"   Found {len(nitrate_by_site)} active stations in period '{period_label}' "
        f"({n_slice_samples} samples; {n_slice_cens} censored '<' "
        f"[{100.0 * n_slice_cens / max(1, n_slice_samples):.2f}%] imputed via '{args.censor_method}')."
    )

    print("2. Scanning ORN_v2.gpkg geometry to snap RSN stations to stream reaches...")
    conn = sqlite3.connect(GPKG_PATH)
    c = conn.cursor()

    # Snap all RSN stations (water quality + ecology) so rsn_sites.csv covers every dataset
    all_station_ids = sorted(pts.keys())
    site_grid = defaultdict(list)
    for nid in all_station_ids:
        px, py, _ = pts[nid]
        site_grid[(int(px // 2000), int(py // 2000))].append((nid, px, py))

    best_snap = {nid: (1e18, None, None, 0.0, 0.0) for nid in all_station_ids}

    for row in c.execute("SELECT RivID, CatchID, Dist2Mth, Shape FROM ORN"):
        riv_id, catch_id, d2m, blob = row
        # Quick bounding box check from GeoPackage header
        flags = blob[3]
        env_ind = (flags >> 1) & 0x07
        if env_ind >= 1:
            minx, maxx, miny, maxy = struct.unpack_from("<4d", blob, 8)
        else:
            poly = parse_gpkg_multilinestring(blob)
            xs = [p[0] for p in poly]
            ys = [p[1] for p in poly]
            minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        gx0, gx1 = int((minx - 1000) // 2000), int((maxx + 1000) // 2000)
        gy0, gy1 = int((miny - 1000) // 2000), int((maxy + 1000) // 2000)
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

    # Save station-level snapping summary (keeping the 141 monitored catchments constant)
    monitored_catchments = {
        best_snap[nid][2]
        for nid in core_wq_sites
        if nid in best_snap and best_snap[nid][0] <= 250.0
    }
    site_by_reach = defaultdict(list)
    sites_csv_path = os.path.join(OUT_DIR, "rsn_sites.csv")
    with open(sites_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "site_id", "wims_id", "wims_region_code",
            "easting", "northing", "snap_dist_m", "RivID", "CatchID",
            "obs_up_dist_km", "n_samples", "nitrate_mean", "nitrate_median",
            "nitrate_sd", "log_nitrate_mean", "mean_log_nitrate", "included_250m"
        ])
        for nid in all_station_ids:
            px, py, meta = pts[nid]
            d, riv_id, catch_id, along, d2m = best_snap[nid]
            if catch_id not in monitored_catchments and nid not in core_wq_sites:
                continue
            obs_up_km = max(0.0, (d2m - along) / 1000.0)
            inc = 1 if (d <= 250.0 and catch_id in monitored_catchments) else 0
            vals = nitrate_by_site.get(nid, [])
            if vals:
                m_val = statistics.mean(vals)
                med_val = statistics.median(vals)
                sd_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                log_m = round(math.log(m_val), 5) if m_val > 0 else ""
                m_log = round(statistics.mean(math.log(v) for v in vals if v > 0), 5) if any(v > 0 for v in vals) else ""
                if inc:
                    site_by_reach[riv_id].append((nid, d, along, vals))
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), len(vals), round(m_val, 4), round(med_val, 4),
                    round(sd_val, 4), log_m, m_log, inc
                ])
            else:
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), 0, "", "", "", "", "", inc
                ])

    # Also export multi-period reach aggregation table (all, year, season, month, year_season)
    reach_period_vals = defaultdict(list)
    reach_along_map = defaultdict(list)
    reach_catch_map = {}
    reach_d2m_map = {}
    for nid, imp_val, yr, mo, seas, _ in imputed_records:
        if nid not in best_snap:
            continue
        d, riv_id, catch_id, along, d2m = best_snap[nid]
        if d <= 250.0 and catch_id in monitored_catchments:
            reach_along_map[riv_id].append(along)
            reach_catch_map[riv_id] = catch_id
            reach_d2m_map[riv_id] = d2m
            reach_period_vals[(riv_id, "all", "all")].append(imp_val)
            reach_period_vals[(riv_id, "year", str(yr))].append(imp_val)
            reach_period_vals[(riv_id, "season", seas)].append(imp_val)
            reach_period_vals[(riv_id, "month", f"{mo:02d}")].append(imp_val)
            reach_period_vals[(riv_id, "year_season", f"{yr}_{seas}")].append(imp_val)

    periods_csv_path = os.path.join(OUT_DIR, "rsn_reach_periods.csv")
    with open(periods_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "RivID", "CatchID", "period_type", "period_value",
            "obs_up_dist_km", "n_samples", "nitrate_mean", "nitrate_median", "log_nitrate_mean"
        ])
        for (riv_id, p_type, p_val), vals in sorted(reach_period_vals.items()):
            d2m_km = reach_d2m_map[riv_id] / 1000.0
            along_m = statistics.mean(reach_along_map[riv_id])
            obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
            m_v = statistics.mean(vals)
            med_v = statistics.median(vals)
            log_mv = round(math.log(m_v), 5) if m_v > 0 else ""
            w.writerow([
                riv_id, reach_catch_map[riv_id], p_type, p_val,
                round(obs_up_km, 5), len(vals), round(m_v, 4), round(med_v, 4), log_mv
            ])

    print(
        f"   Snapped {sum(1 for nid in best_snap if best_snap[nid][0] <= 250.0)} sites (<=250m) "
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
    catch_main_river = {}  # catch_id -> (max_shreve, name1, rbd_id)
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
        tot = 0.0
        segs = []
        for i in range(len(poly) - 1):
            d = math.hypot(poly[i + 1][0] - poly[i][0], poly[i + 1][1] - poly[i][1])
            segs.append((d, poly[i], poly[i + 1]))
            tot += d
        reach_len[riv_id] = tot
        half = tot * 0.5
        cur = 0.0
        mx, my = poly[0]
        for d, p1, p2 in segs:
            if cur + d >= half and d > 0:
                t = (half - cur) / d
                mx = p1[0] + t * (p2[0] - p1[0])
                my = p1[1] + t * (p2[1] - p1[1])
                break
            cur += d
        reach_coords[riv_id] = (poly[0][0], poly[0][1], mx, my, poly[-1][0], poly[-1][1])

    catch_obs_cnt = defaultdict(int)
    for riv_id in site_by_reach:
        if riv_id in reaches:
            catch_obs_cnt[reaches[riv_id][3]] += 1

    catch_lu_path = os.path.join(OUT_DIR, "rsn_catchments_lookup.csv")
    with open(catch_lu_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["CatchID", "Main_River", "RBD_ID", "RBD_Name", "n_reaches", "n_obs_reaches"])
        for cid in sorted(monitored_catchments, key=lambda k: (-catch_obs_cnt[k], -catch_reach_cnt[k])):
            _, nm1, rbd_id = catch_main_river.get(cid, (0, f"Catchment {cid}", 0))
            w.writerow([cid, nm1, rbd_id, rbd_lu.get(rbd_id, ""), catch_reach_cnt[cid], catch_obs_cnt[cid]])

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
    topo_order = []  # Downstream-to-upstream order within each catchment
    for m_id in sorted(mouths, key=lambda x: reaches[x][3]):
        log_afv[m_id] = 0.0
        stack = [m_id]
        while stack:
            u = stack.pop()
            topo_order.append(u)
            ps = parents[u]
            if ps:
                sum_us = sum(max(0.0001, reaches[p][6]) for p in ps)
                for p in ps:
                    parent_dn[p] = u
                    log_afv[p] = log_afv[u] + math.log(max(0.0001, reaches[p][6]) / sum_us)
                    stack.append(p)

    # Multi-pass tree neighbor imputation for short reaches (<100m) without a 100m point
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

    missing_cnt = sum(1 for u in topo_order if u not in reach_100m or reach_100m[u][0] is None)
    print(f"   Tree traversal complete ({len(topo_order)} reaches, {missing_cnt} missing covariates).")

    out_reaches_path = os.path.join(OUT_DIR, "england_rsn_reaches.csv")
    print(f"5. Writing {out_reaches_path}...")
    with open(out_reaches_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "RivID", "CatchID", "RBD_ID", "parent_dn",
            "Dist2Mth_km", "up_dist_mid_km", "length_km", "US_Accum_km",
            "log_AFV", "Shreve", "Strahler",
            "start_X", "start_Y", "mid_X", "mid_Y", "end_X", "end_Y",
            "Elevation", "Slope", "Sinuosity", "SrcElevation",
            "Is_Chalk_Limestone", "Is_Sandstone", "Is_Aquifer", "Is_Peat",
            "has_obs", "site_id", "snap_dist_m", "obs_up_dist_km",
            "nitrate_mean", "nitrate_median", "log_nitrate_mean", "n_samples"
        ])
        for riv_id in topo_order:
            r = reaches[riv_id]
            rc = reach_100m.get(riv_id, (50.0, 0.5, 1.05, 0, 0, 100.0))
            b_desc = bedrock_lu.get(rc[3], "")
            s_desc = superficial_lu.get(rc[4], "")
            is_chalk_lime = 1 if ("CHALK" in b_desc or "LIMESTONE" in b_desc or "DOLOSTONE" in b_desc) else 0
            is_sand = 1 if ("SANDSTONE" in b_desc and "MUDSTONE" not in b_desc and "SILTSTONE" not in b_desc) else 0
            is_aq = 1 if (is_chalk_lime or is_sand) else 0
            is_peat = 1 if "PEAT" in s_desc else 0

            d2m_km = r[5] / 1000.0
            l_km = reach_len[riv_id] / 1000.0
            up_mid_km = max(0.0, d2m_km - 0.5 * l_km)
            sx, sy, mx, my, ex, ey = reach_coords[riv_id]

            if riv_id in site_by_reach:
                s_list = site_by_reach[riv_id]
                all_vals = [v for s in s_list for v in s[3]]
                along_m = statistics.mean(s[2] for s in s_list)
                snap_d = min(s[1] for s in s_list)
                obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
                m_no3 = statistics.mean(all_vals)
                med_no3 = statistics.median(all_vals)
                log_no3 = round(math.log(m_no3), 5) if m_no3 > 0 else ""
                w.writerow([
                    riv_id, r[3], r[4], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(r[6] / 1000.0, 5),
                    round(log_afv[riv_id], 6), r[7], r[8],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(rc[0] or 0.0, 3), round(rc[1] or 0.0, 4), round(rc[2] or 1.0, 4), round(rc[5] or rc[0] or 0.0, 3),
                    is_chalk_lime, is_sand, is_aq, is_peat,
                    1, s_list[0][0], round(snap_d, 2), round(obs_up_km, 5),
                    round(m_no3, 4), round(med_no3, 4), log_no3, len(all_vals)
                ])
            else:
                w.writerow([
                    riv_id, r[3], r[4], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(r[6] / 1000.0, 5),
                    round(log_afv[riv_id], 6), r[7], r[8],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(rc[0] or 0.0, 3), round(rc[1] or 0.0, 4), round(rc[2] or 1.0, 4), round(rc[5] or rc[0] or 0.0, 3),
                    is_chalk_lime, is_sand, is_aq, is_peat,
                    0, "", "", round(up_mid_km, 5),
                    "", "", "", 0
                ])

    print(f"Done in {time.time() - t0:.2f} seconds.")


if __name__ == "__main__":
    main()
