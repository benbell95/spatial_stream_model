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
     via `/usr/bin/cs2cs`.
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
import glob
import json
import math
import os
import shutil
import sqlite3
import statistics
import struct
import subprocess
import time
import zipfile
from collections import Counter, defaultdict

BASE_DIR = "/home/ben/code/r/spatial_stream_model"
WQ_DIR = os.path.join(BASE_DIR, "data/NCEA_River_Surveillance_Network_ARD_Open")
CONFIG_DIR = os.path.join(BASE_DIR, "configs/networks")
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
    if method == "half_lod":
        return {
            nid: [0.5 * v if is_c else v for v, is_c in recs]
            for nid, recs in raw_by_site.items()
        }

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


def parse_gpkg_geometry(blob):
    """
    Universal parser for GeoPackage WKB LineString / MultiLineString geometries.
    Supports 2D (XY), 3D (XYZ / XYM), and 4D (XYZM) in both ISO WKB (1000/2000/3000)
    and OGC EWKB (0x80000000 / 0x40000000), returning `(pts_xy, zs)` where:
      - `pts_xy` : list of (x, y) tuples
      - `zs`     : list of valid vertex Z elevations (if 3D/4D Z is present, else [])
    """
    if not blob or len(blob) < 9:
        return [], []
    env_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
    flags = blob[3]
    env_ind = (flags >> 1) & 0x07
    offset = 8 + env_sizes.get(env_ind, 0)
    if offset >= len(blob):
        return [], []

    bo = "<" if blob[offset] == 1 else ">"
    wkb_type = struct.unpack_from(bo + "I", blob, offset + 1)[0]
    low16 = wkb_type & 0xFFFF
    has_z = (wkb_type & 0x80000000) != 0 or (1000 <= low16 < 4000 and (low16 // 1000) in (1, 3))
    has_m = (wkb_type & 0x40000000) != 0 or (1000 <= low16 < 4000 and (low16 // 1000) in (2, 3))
    ndim = 2 + (1 if has_z else 0) + (1 if has_m else 0)
    base_type = low16 % 1000

    pts = []
    zs = []
    if base_type == 5:  # MultiLineString
        num_geoms = struct.unpack_from(bo + "I", blob, offset + 5)[0]
        offset += 9
        for _ in range(num_geoms):
            g_bo = "<" if blob[offset] == 1 else ">"
            g_type = struct.unpack_from(g_bo + "I", blob, offset + 1)[0]
            g_low16 = g_type & 0xFFFF
            g_z = (g_type & 0x80000000) != 0 or (1000 <= g_low16 < 4000 and (g_low16 // 1000) in (1, 3))
            g_m = (g_type & 0x40000000) != 0 or (1000 <= g_low16 < 4000 and (g_low16 // 1000) in (2, 3))
            g_dim = 2 + (1 if g_z else 0) + (1 if g_m else 0)
            npts = struct.unpack_from(g_bo + "I", blob, offset + 5)[0]
            offset += 9
            coords = struct.unpack_from(g_bo + f"{g_dim * npts}d", blob, offset)
            offset += 8 * g_dim * npts
            for i in range(0, g_dim * npts, g_dim):
                pts.append((coords[i], coords[i + 1]))
                if g_z:
                    zv = coords[i + 2]
                    if not math.isnan(zv) and -500.0 < zv < 9000.0:
                        zs.append(zv)
    elif base_type == 2:  # LineString
        npts = struct.unpack_from(bo + "I", blob, offset + 5)[0]
        offset += 9
        coords = struct.unpack_from(bo + f"{ndim * npts}d", blob, offset)
        for i in range(0, ndim * npts, ndim):
            pts.append((coords[i], coords[i + 1]))
            if has_z:
                zv = coords[i + 2]
                if not math.isnan(zv) and -500.0 < zv < 9000.0:
                    zs.append(zv)
    return pts, zs


def parse_gpkg_multilinestring(blob):
    """Backward-compatible helper returning only 2D (x, y) vertices."""
    pts, _ = parse_gpkg_geometry(blob)
    return pts


def transform_coords_batch(xy_list, from_epsg, to_epsg):
    """
    Transforms a list of (x, y) coordinates from `from_epsg` to `to_epsg` using
    `/usr/bin/cs2cs` with traditional GIS easting/northing (`+init=epsg:...`) order.
    """
    if not xy_list or from_epsg == to_epsg:
        return list(xy_list)
    inp = "\n".join(f"{x:.4f} {y:.4f}" for x, y in xy_list)
    cmd = [
        "cs2cs",
        "-f",
        "%.3f",
        f"+init=epsg:{int(from_epsg)}",
        "+to",
        f"+init=epsg:{int(to_epsg)}",
    ]
    proc = subprocess.run(cmd, input=inp, text=True, capture_output=True, check=True)
    out = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out.append((float(parts[0]), float(parts[1])))
    if len(out) != len(xy_list):
        raise RuntimeError(f"cs2cs returned {len(out)} coordinates for {len(xy_list)} inputs.")
    return out


def point_to_polyline(px, py, poly):
    """Returns (min_dist_m, distance_along_from_start_m) from point (px, py) to poly."""
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


def compute_polyline_metrics(poly):
    """Returns (total_length_m, sinuosity, (start_X, start_Y, mid_X, mid_Y, end_X, end_Y))."""
    if not poly:
        return 0.0, 1.0, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    if len(poly) == 1:
        x0, y0 = poly[0]
        return 0.0, 1.0, (x0, y0, x0, y0, x0, y0)
    tot = 0.0
    segs = []
    for i in range(len(poly) - 1):
        d = math.hypot(poly[i + 1][0] - poly[i][0], poly[i + 1][1] - poly[i][1])
        segs.append((d, poly[i], poly[i + 1]))
        tot += d
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
    straight = math.hypot(poly[-1][0] - poly[0][0], poly[-1][1] - poly[0][1])
    sinuosity = max(1.0, tot / straight) if straight > 1e-3 else 1.0
    return tot, sinuosity, (poly[0][0], poly[0][1], mx, my, poly[-1][0], poly[-1][1])


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


def list_available_networks():
    print("Available Water Network Profiles (`--network=<name>`):\n")
    print(f"{'Network Name':<18} | {'Short Label':<20} | {'CRS':<10} | {'GeoPackage Path(s)'}")
    print("-" * 105)
    for p in sorted(glob.glob(os.path.join(CONFIG_DIR, "*.json"))):
        with open(p, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        gpkgs = ", ".join(cfg.get("gpkg_paths", []))
        print(f"{cfg['abbr']:<18} | {cfg.get('short_label', cfg['abbr']):<20} | EPSG:{cfg.get('crs_epsg', 27700):<5} | {gpkgs}")
    print("\nYou can also pass a custom JSON config (`--network-config=path/to/net.json`) or a `.gpkg` path directly.")


def load_network_config(network_arg, config_path=None):
    if config_path:
        p = config_path if os.path.isabs(config_path) else os.path.join(BASE_DIR, config_path)
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)

    net_clean = network_arg.strip()
    cand_json = os.path.join(CONFIG_DIR, f"{net_clean.lower()}.json")
    if os.path.exists(cand_json):
        with open(cand_json, "r", encoding="utf-8") as f:
            return json.load(f)

    # Check common aliases
    alias_map = {
        "orn_v2": "orn",
        "os_orn": "orn",
        "thames": "euhydro_thames",
        "euhydro-thames": "euhydro_thames",
        "eu_hydro_thames": "euhydro_thames",
        "tweed": "euhydro_tweed",
        "euhydro-tweed": "euhydro_tweed",
        "eu_hydro_tweed": "euhydro_tweed",
        "eu_hydro": "euhydro",
        "eu-hydro": "euhydro",
    }
    if net_clean.lower() in alias_map:
        mapped = alias_map[net_clean.lower()]
        with open(os.path.join(CONFIG_DIR, f"{mapped}.json"), "r", encoding="utf-8") as f:
            return json.load(f)

    if net_clean.endswith(".json") and os.path.exists(net_clean):
        with open(net_clean, "r", encoding="utf-8") as f:
            return json.load(f)

    # Direct GeoPackage path fallback: auto-inspect schema
    gpkg_p = net_clean if os.path.isabs(net_clean) else os.path.join(BASE_DIR, net_clean)
    if os.path.exists(gpkg_p) and gpkg_p.endswith(".gpkg"):
        conn = sqlite3.connect(gpkg_p)
        c = conn.cursor()
        row = c.execute("SELECT table_name, srs_id FROM gpkg_contents WHERE data_type='features' LIMIT 1").fetchone()
        layer, srs_id = row[0], row[1]
        cols = {r[1].lower(): r[1] for r in c.execute(f'PRAGMA table_info("{layer}")').fetchall()}
        conn.close()
        abbr = os.path.splitext(os.path.basename(gpkg_p))[0].lower()
        return {
            "name": abbr,
            "abbr": abbr,
            "label": abbr,
            "short_label": abbr,
            "gpkg_paths": [net_clean],
            "layers": [layer],
            "geom_col": cols.get("shape", cols.get("geom", "Shape")),
            "crs_epsg": srs_id or 27700,
            "topology": {
                "mode": "pointer" if "nextdownid" in cols else ("nodes" if "fnode" in cols else "geometry"),
                "id_col": cols.get("rivid", cols.get("object_id", cols.get("objectid", "FID"))),
                "from_node_col": cols.get("fnode"),
                "to_node_col": cols.get("tnode"),
                "next_down_col": cols.get("nextdownid"),
                "catch_id_col": cols.get("catchid"),
                "rbd_id_col": cols.get("rbd_id"),
            },
            "attributes": {
                "dist_to_mouth_col": cols.get("dist2mth"),
                "dist_to_mouth_unit": "m",
                "afv_weight_col": cols.get("us_accum", cols.get("cum_len")),
                "afv_weight_unit": "m",
                "shreve_col": cols.get("shreve"),
                "strahler_col": cols.get("strahler"),
                "slope_col": cols.get("slope", cols.get("pente")),
                "elevation_col": cols.get("elevation"),
                "river_name_cols": [c_nm for k, c_nm in cols.items() if "name" in k],
            },
            "covariates": {"mode": "geometry_and_spatial_join"},
            "max_snap_dist_m": 250.0,
        }

    raise ValueError(f"Unknown network '{network_arg}'. Run with --list-networks to see available options.")


def ensure_gpkg_paths(cfg):
    """Ensures all configured GeoPackage files exist, extracting zip_fallbacks if needed."""
    resolved = []
    zip_list = cfg.get("zip_fallbacks", [])
    for idx, rel_p in enumerate(cfg.get("gpkg_paths", [])):
        full_p = rel_p if os.path.isabs(rel_p) else os.path.join(BASE_DIR, rel_p)
        if not os.path.exists(full_p) and idx < len(zip_list):
            zip_rel = zip_list[idx]
            zip_p = zip_rel if os.path.isabs(zip_rel) else os.path.join(BASE_DIR, zip_rel)
            if os.path.exists(zip_p):
                print(f"   Extracting {zip_p}...")
                with zipfile.ZipFile(zip_p, "r") as zf:
                    zf.extractall(os.path.dirname(zip_p))
        if not os.path.exists(full_p):
            raise FileNotFoundError(f"GeoPackage file not found: {full_p}")
        resolved.append(full_p)
    return resolved


def write_outputs_to_dirs(net_abbr, files_dict):
    """
    Writes generated output files into `data/processed/<net_abbr>/` AND copies them
    into `data/processed/` so downstream scripts work with or without `--network`.
    """
    net_dir = os.path.join(OUT_DIR, net_abbr)
    os.makedirs(net_dir, exist_ok=True)
    for fname, src_path in files_dict.items():
        dst_net = os.path.join(net_dir, fname)
        if os.path.abspath(src_path) != os.path.abspath(dst_net):
            shutil.copy2(src_path, dst_net)


def load_orn_spatial_reference():
    """
    Loads the pre-extracted ORN England reaches (`data/processed/orn/england_rsn_reaches.csv`)
    into a fast 2 km spatial hash grid in EPSG:27700 so non-ORN networks in England
    can inherit BGS bedrock/superficial geology (`Is_Chalk_Limestone`, `Is_Sandstone`,
    `Is_Aquifer`, `Is_Peat`), `RBD_ID`, and fallback `Elevation`/`Slope` in < 0.2s.
    """
    orn_csv = os.path.join(OUT_DIR, "orn", "england_rsn_reaches.csv")
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
    orn_lu = os.path.join(OUT_DIR, "orn", "rsn_catchments_lookup.csv")
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


def preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, nitrate_by_site, period_label, t0):
    """
    Network-agnostic preprocessing pipeline for any configured GeoPackage network
    (such as EU-Hydro Thames, EU-Hydro Tweed, combined EU-Hydro, OS Water Network, etc.).
    """
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

    # Scan GeoPackage(s) to snap stations and load reaches
    raw_reaches = []
    node_lookup = {}  # (file_idx, node_id) -> (elev, len_tom)
    nodes_cfg = cfg.get("nodes_table")

    for file_idx, gpkg_p in enumerate(gpkg_paths):
        conn = sqlite3.connect(gpkg_p)
        conn.row_factory = sqlite3.Row
        c = conn.cursor()

        if nodes_cfg:
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
                poly, zs = parse_gpkg_geometry(blob)
                if len(poly) < 2:
                    continue

                orig_id = str(rdict[id_col])
                raw_cid = rdict.get(catch_col) if catch_col else None
                cid_key = (file_idx, int(raw_cid) if raw_cid is not None else 0)

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

                for nid, px, py in cand_sites:
                    d, along = point_to_polyline(px, py, poly)
                    if d < best_snap[nid][0]:
                        best_snap[nid] = (d, (file_idx, orig_id), cid_key, along)

                raw_reaches.append((file_idx, orig_id, cid_key, poly, zs, rdict))
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

    # Assign unique positive integer RivIDs
    key_to_rivid = {}
    reaches_by_id = {}
    for idx, (file_idx, orig_id, cid_key, poly, zs, rdict) in enumerate(filtered_reaches, 1):
        if len(gpkg_paths) == 1 and len(layers) == 1 and str(orig_id).isdigit():
            riv_id = int(orig_id)
        else:
            riv_id = idx
        key_to_rivid[(file_idx, orig_id)] = riv_id
        reaches_by_id[riv_id] = {
            "riv_id": riv_id,
            "file_idx": file_idx,
            "orig_id": orig_id,
            "catch_id": cid_map[cid_key],
            "poly": poly,
            "zs": zs,
            "row": rdict,
        }

    print(
        f"   Loaded {len(reaches_by_id)} reaches across {len(monitored_cid_keys)} monitored catchments "
        f"(from {len(raw_reaches)} total reaches in GeoPackage)."
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
        p0 = info["poly"][0]
        start_pt_grid[(f_idx, int(p0[0] // 50), int(p0[1] // 50))].append((riv_id, p0[0], p0[1]))

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
                dn_id = key_to_rivid.get((f_idx, first_tok))

        # Mode B / Fallback: Nodes (Tnode -> Fnode)
        if dn_id is None and topo_mode in ("nodes", "pointer") and tnode_col:
            tn_val = rdict.get(tnode_col)
            if tn_val:
                dn_id = fnode_index.get((f_idx, info["catch_id"], str(tn_val)))

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
            # Keep within same CatchID
            if reaches_by_id[dn_id]["catch_id"] == info["catch_id"]:
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
            raw_d2m = u_info["row"].get(d2m_col) if d2m_col else None
            if raw_d2m is not None and float(raw_d2m) > 0:
                u_info["d2m_km"] = float(raw_d2m) * d2m_scale
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
        raw_shr = rdict.get(shreve_col) if shreve_col else None
        if raw_shr is not None and int(float(raw_shr)) > 0:
            u_info["shreve"] = int(float(raw_shr))
        else:
            u_info["shreve"] = sum(reaches_by_id[p]["shreve"] for p in ps) if ps else 1

        # Strahler order
        raw_str = rdict.get(strahler_col) if strahler_col else None
        if raw_str is not None and int(float(raw_str)) > 0:
            u_info["strahler"] = int(float(raw_str))
        elif not ps:
            u_info["strahler"] = 1
        else:
            child_s = [reaches_by_id[p]["strahler"] for p in ps]
            mx_s = max(child_s)
            u_info["strahler"] = mx_s + 1 if child_s.count(mx_s) >= 2 else mx_s

        # US_Accum_km (AFV weight)
        raw_afv = rdict.get(afv_col) if afv_col else None
        if raw_afv is not None and float(raw_afv) > 0:
            u_info["us_accum_km"] = max(1e-4, float(raw_afv) * afv_scale)
        else:
            u_info["us_accum_km"] = u_info["length_km"] + sum(reaches_by_id[p]["us_accum_km"] for p in ps)

    # Top-down pass: compute log_AFV and propagate headwater SrcElevation
    log_afv = {}
    for m_id in sorted(mouths, key=lambda x: (reaches_by_id[x]["catch_id"], x)):
        log_afv[m_id] = 0.0
        stack = [m_id]
        while stack:
            u = stack.pop()
            ps = parents[u]
            if ps:
                sum_us = sum(max(1e-4, reaches_by_id[p]["us_accum_km"]) for p in ps)
                for p in ps:
                    log_afv[p] = log_afv[u] + math.log(max(1e-4, reaches_by_id[p]["us_accum_km"]) / sum_us)
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
    orn_grid, orn_rbd_names = load_orn_spatial_reference()
    for u in topo_order:
        u_info = reaches_by_id[u]
        rdict = u_info["row"]
        f_idx = u_info["file_idx"]
        zs = u_info["zs"]
        _, _, mx, my, _, _ = u_info["coords_bng"]
        orn_ref = query_orn_spatial_grid(orn_grid, mx, my)

        # Elevation: 1. explicit attribute, 2. 3D WKB Z mean, 3. Nodes table FNODE/TNODE, 4. ORN spatial join
        elev = None
        if elev_col and rdict.get(elev_col) is not None:
            elev = float(rdict[elev_col])
        elif zs:
            elev = statistics.mean(zs)
        elif node_lookup:
            fn_el = node_lookup.get((f_idx, str(rdict.get(fnode_col))), (None, None))[0]
            tn_el = node_lookup.get((f_idx, str(rdict.get(tnode_col))), (None, None))[0]
            vals_el = [x for x in (fn_el, tn_el) if x is not None]
            if vals_el:
                elev = statistics.mean(vals_el)
        if elev is None:
            elev = orn_ref[3] if orn_ref else 50.0
        u_info["elevation"] = elev

        # Slope (%): 1. explicit attribute (e.g. PENTE), 2. 3D WKB Z drop, 3. ORN spatial join
        slope = None
        if slope_col and rdict.get(slope_col) is not None and float(rdict[slope_col]) >= 0:
            slope = float(rdict[slope_col])
        elif len(zs) >= 2 and u_info["length_km"] > 0:
            slope = max(0.0, (zs[0] - zs[-1]) / (u_info["length_km"] * 10.0))
        if slope is None:
            slope = orn_ref[4] if orn_ref else 0.5
        u_info["slope"] = slope

        # Geology & RBD_ID from ORN spatial reference (if in England)
        if orn_ref is not None:
            u_info["rbd_id"] = int(rdict[rbd_col]) if (rbd_col and rdict.get(rbd_col)) else orn_ref[2]
            u_info["src_elevation"] = max(elev, orn_ref[5])
            u_info["is_chalk_lime"] = orn_ref[6]
            u_info["is_sand"] = orn_ref[7]
            u_info["is_aq"] = orn_ref[8]
            u_info["is_peat"] = orn_ref[9]
        else:
            u_info["rbd_id"] = int(rdict[rbd_col]) if (rbd_col and rdict.get(rbd_col)) else 1
            u_info["src_elevation"] = elev
            u_info["is_chalk_lime"] = 0
            u_info["is_sand"] = 0
            u_info["is_aq"] = 0
            u_info["is_peat"] = 0

    # Write station-level snapping summary (rsn_sites.csv) and multi-period reach table (rsn_reach_periods.csv)
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
            d, r_key, cid_key, along_m = best_snap[nid]
            if r_key not in key_to_rivid:
                continue
            riv_id = key_to_rivid[r_key]
            u_info = reaches_by_id[riv_id]
            catch_id = u_info["catch_id"]
            d2m_km = u_info["d2m_km"]
            obs_up_km = max(0.0, d2m_km - (along_m / 1000.0))
            inc = 1 if d <= max_snap_m else 0
            vals = nitrate_by_site.get(nid, [])
            if vals:
                m_val = statistics.mean(vals)
                med_val = statistics.median(vals)
                sd_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                log_m = round(math.log(m_val), 5) if m_val > 0 else ""
                m_log = round(statistics.mean(math.log(v) for v in vals if v > 0), 5) if any(v > 0 for v in vals) else ""
                if inc:
                    site_by_reach[riv_id].append((nid, d, along_m, vals))
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

    reach_period_vals = defaultdict(list)
    reach_along_map = defaultdict(list)
    for nid, imp_val, yr, mo, seas, _ in imputed_records:
        if nid not in best_snap:
            continue
        d, r_key, cid_key, along_m = best_snap[nid]
        if d <= max_snap_m and r_key in key_to_rivid:
            riv_id = key_to_rivid[r_key]
            reach_along_map[riv_id].append(along_m)
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
            u_info = reaches_by_id[riv_id]
            d2m_km = u_info["d2m_km"]
            along_m = statistics.mean(reach_along_map[riv_id])
            obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
            m_v = statistics.mean(vals)
            med_v = statistics.median(vals)
            log_mv = round(math.log(m_v), 5) if m_v > 0 else ""
            w.writerow([
                riv_id, u_info["catch_id"], p_type, p_val,
                round(obs_up_km, 5), len(vals), round(m_v, 4), round(med_v, 4), log_mv
            ])

    # Build catchment lookup table (rsn_catchments_lookup.csv)
    catch_reaches = defaultdict(list)
    for u in topo_order:
        catch_reaches[reaches_by_id[u]["catch_id"]].append(reaches_by_id[u])

    catch_lu_path = os.path.join(OUT_DIR, "rsn_catchments_lookup.csv")
    with open(catch_lu_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["CatchID", "Main_River", "RBD_ID", "RBD_Name", "n_reaches", "n_obs_reaches"])
        sorted_cids = sorted(
            catch_reaches.keys(),
            key=lambda cid: (-sum(1 for r in catch_reaches[cid] if r["riv_id"] in site_by_reach), -len(catch_reaches[cid]))
        )
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
            w.writerow([cid, best_name, rbd_id, rbd_name, len(rlist), n_obs_c])

    # Write england_rsn_reaches.csv (including any configured extra_numeric_cols)
    out_reaches_path = os.path.join(OUT_DIR, "england_rsn_reaches.csv")
    print(f"5. Writing {out_reaches_path} ({len(topo_order)} reaches, {len(site_by_reach)} observed)...")
    base_header = [
        "RivID", "CatchID", "RBD_ID", "parent_dn",
        "Dist2Mth_km", "up_dist_mid_km", "length_km", "US_Accum_km",
        "log_AFV", "Shreve", "Strahler",
        "start_X", "start_Y", "mid_X", "mid_Y", "end_X", "end_Y",
        "Elevation", "Slope", "Sinuosity", "SrcElevation",
        "Is_Chalk_Limestone", "Is_Sandstone", "Is_Aquifer", "Is_Peat",
        "has_obs", "site_id", "snap_dist_m", "obs_up_dist_km",
        "nitrate_mean", "nitrate_median", "log_nitrate_mean", "n_samples"
    ]
    with open(out_reaches_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(base_header + extra_num_cols)
        for riv_id in topo_order:
            u_info = reaches_by_id[riv_id]
            d2m_km = u_info["d2m_km"]
            l_km = u_info["length_km"]
            up_mid_km = max(0.0, d2m_km - 0.5 * l_km)
            sx, sy, mx, my, ex, ey = u_info["coords_bng"]
            extra_vals = [
                round(float(u_info["row"].get(ec) or 0.0), 4)
                for ec in extra_num_cols
            ]

            if riv_id in site_by_reach:
                s_list = site_by_reach[riv_id]
                all_vals = [v for s in s_list for v in s[3]]
                along_m = statistics.mean(s[2] for s in s_list)
                snap_d = min(s[1] for s in s_list)
                obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
                m_no3 = statistics.mean(all_vals)
                med_no3 = statistics.median(all_vals)
                log_no3 = round(math.log(m_no3), 5) if m_no3 > 0 else ""
                row_out = [
                    riv_id, u_info["catch_id"], u_info["rbd_id"], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(u_info["us_accum_km"], 5),
                    round(log_afv[riv_id], 6), u_info["shreve"], u_info["strahler"],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(u_info["elevation"], 3), round(u_info["slope"], 4), round(u_info["sinuosity"], 4), round(u_info["src_elevation"], 3),
                    u_info["is_chalk_lime"], u_info["is_sand"], u_info["is_aq"], u_info["is_peat"],
                    1, s_list[0][0], round(snap_d, 2), round(obs_up_km, 5),
                    round(m_no3, 4), round(med_no3, 4), log_no3, len(all_vals)
                ]
            else:
                row_out = [
                    riv_id, u_info["catch_id"], u_info["rbd_id"], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(u_info["us_accum_km"], 5),
                    round(log_afv[riv_id], 6), u_info["shreve"], u_info["strahler"],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(u_info["elevation"], 3), round(u_info["slope"], 4), round(u_info["sinuosity"], 4), round(u_info["src_elevation"], 3),
                    u_info["is_chalk_lime"], u_info["is_sand"], u_info["is_aq"], u_info["is_peat"],
                    0, "", "", round(up_mid_km, 5),
                    "", "", "", 0
                ]
            w.writerow(row_out + extra_vals)

    meta_path = os.path.join(OUT_DIR, "network_metadata.json")
    meta_obj = {
        "name": cfg.get("name", net_abbr),
        "abbr": net_abbr,
        "label": net_label,
        "short_label": short_label,
        "crs_epsg": crs_epsg,
        "topology_mode": topo_mode,
        "max_snap_dist_m": max_snap_m,
        "n_reaches": len(topo_order),
        "n_obs_reaches": len(site_by_reach),
        "n_catchments": len(catch_reaches),
        "core_basins": cfg.get("core_basins", {}),
        "extra_numeric_cols": extra_num_cols,
    }
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_obj, f, indent=2)

    write_outputs_to_dirs(net_abbr, {
        "england_rsn_reaches.csv": out_reaches_path,
        "rsn_sites.csv": sites_csv_path,
        "rsn_reach_periods.csv": periods_csv_path,
        "rsn_catchments_lookup.csv": catch_lu_path,
        "network_metadata.json": meta_path,
    })
    print(f"Done in {time.time() - t0:.2f} seconds (saved to data/processed/{net_abbr}/ and data/processed/).")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess stream network GeoPackages (ORN v2, EU-Hydro, etc.) and EA RSN samples for S3N modeling."
    )
    parser.add_argument(
        "--network",
        default="orn",
        help="Water network preset ('orn', 'euhydro_thames', 'euhydro_tweed', 'euhydro') or path to GeoPackage/JSON.",
    )
    parser.add_argument(
        "--network-config",
        default=None,
        help="Optional path to a custom network JSON configuration file.",
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
    args = parser.parse_args()

    if args.list_networks:
        list_available_networks()
        return

    cfg = load_network_config(args.network, args.network_config)
    net_abbr = cfg["abbr"]

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
        f"(network='{net_abbr}', dataset='{args.dataset}', variable='{args.det_code}', "
        f"censor_method='{args.censor_method}', period='{period_label}')..."
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

    core_wq_sites = set()
    with open(os.path.join(WQ_DIR, "RSN Water Quality Samples.csv"), "r", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r["det_code"] == "0117":
                core_wq_sites.add(r["network_id"])

    all_records = []
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

    raw_all_by_site = defaultdict(list)
    for nid, val, is_cens, yr, mo, seas in all_records:
        raw_all_by_site[nid].append((val, is_cens))
    imp_all_by_site = impute_censored_by_site(raw_all_by_site, method=args.censor_method)

    site_cursors = defaultdict(int)
    imputed_records = []
    for nid, val, is_cens, yr, mo, seas in all_records:
        idx = site_cursors[nid]
        site_cursors[nid] += 1
        imp_val = imp_all_by_site[nid][idx]
        imputed_records.append((nid, imp_val, yr, mo, seas, is_cens))

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

    # Dispatch to generic network pipeline if not using standard ORN fast path
    if cfg.get("covariates", {}).get("mode") != "orn_100m_table" or args.topology_mode != "auto" or args.afv_col:
        preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, nitrate_by_site, period_label, t0)
        return

    # Standard ORN v2 fast path
    max_snap_m = float(args.max_snap_dist if args.max_snap_dist is not None else cfg.get("max_snap_dist_m", 250.0))
    gpkg_path = ensure_gpkg_paths(cfg)[0]
    print(f"2. Scanning {os.path.basename(gpkg_path)} geometry to snap RSN stations (<= {max_snap_m:.0f}m)...")
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

    monitored_catchments = {
        best_snap[nid][2]
        for nid in core_wq_sites
        if nid in best_snap and best_snap[nid][0] <= max_snap_m
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
            inc = 1 if (d <= max_snap_m and catch_id in monitored_catchments) else 0
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

    reach_period_vals = defaultdict(list)
    reach_along_map = defaultdict(list)
    reach_catch_map = {}
    reach_d2m_map = {}
    for nid, imp_val, yr, mo, seas, _ in imputed_records:
        if nid not in best_snap:
            continue
        d, riv_id, catch_id, along, d2m = best_snap[nid]
        if d <= max_snap_m and catch_id in monitored_catchments:
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
    topo_order = []
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

    meta_path = os.path.join(OUT_DIR, "network_metadata.json")
    meta_obj = {
        "name": cfg.get("name", "orn"),
        "abbr": net_abbr,
        "label": cfg.get("label", "OS Open Rivers Network (ORN v2)"),
        "short_label": cfg.get("short_label", "ORN v2"),
        "crs_epsg": 27700,
        "topology_mode": "nodes",
        "max_snap_dist_m": max_snap_m,
        "n_reaches": len(topo_order),
        "n_obs_reaches": len(site_by_reach),
        "n_catchments": len(monitored_catchments),
        "core_basins": cfg.get("core_basins", {"humber": 80, "thames": 3055, "severn": 2786}),
        "extra_numeric_cols": [],
    }
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_obj, f, indent=2)

    write_outputs_to_dirs(net_abbr, {
        "england_rsn_reaches.csv": out_reaches_path,
        "rsn_sites.csv": sites_csv_path,
        "rsn_reach_periods.csv": periods_csv_path,
        "rsn_catchments_lookup.csv": catch_lu_path,
        "network_metadata.json": meta_path,
    })
    print(f"Done in {time.time() - t0:.2f} seconds.")


if __name__ == "__main__":
    main()
