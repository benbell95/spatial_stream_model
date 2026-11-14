#!/usr/bin/env python3
"""
scripts/s3n_py_utils.py

Shared Python utilities for the S3N Stream Network Modeling toolkit:
  - Configuration & RuntimeConfig resolution (configs/directories.json, networks, sources)
  - Coordinate Reference System transforms (pure-Python WGS84->BNG + pyproj/osr/cs2cs batching)
  - Interval-Censored Tobit MLE, Range Parsers, ROS Censoring Imputation, and Power Transforms
  - Universal GeoPackage / WKB / EWKB Geometry & Fast Envelope Readers
  - Source Registry & DateTime Parsers
"""

from dataclasses import dataclass
import bisect
import csv
from functools import lru_cache
import glob
import json
import math
import os
import re
import statistics
import struct
import subprocess
from collections import defaultdict

try:
    from pyproj import Transformer as _PyprojTransformer
except ImportError:
    _PyprojTransformer = None

try:
    from osgeo import osr as _OsgeoOsr
except ImportError:
    _OsgeoOsr = None

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
try:
    DEFAULT_BASE_DIR = "." if os.path.samefile(os.getcwd(), _REPO_ROOT) else os.path.relpath(_REPO_ROOT, os.getcwd())
except OSError:
    DEFAULT_BASE_DIR = "."

_TRANSFORMER_CACHE = {}

MONTH_TO_SEASON = {
    12: "winter", 1: "winter", 2: "winter",
    3: "spring", 4: "spring", 5: "spring",
    6: "summer", 7: "summer", 8: "summer",
    9: "autumn", 10: "autumn", 11: "autumn",
}

MONTH_NAME_MAP = {
    "jan": 1, "january": 1, "feb": 2, "february": 2, "mar": 3, "march": 3,
    "apr": 4, "april": 4, "may": 5, "jun": 6, "june": 6,
    "jul": 7, "july": 7, "aug": 8, "august": 8, "sep": 9, "september": 9,
    "oct": 10, "october": 10, "nov": 11, "november": 11, "dec": 12, "december": 12,
}

TRANSFORM_CATALOG = [
    ("log", "ln", "ln(y)", "y > 0", "Default for water chemistry & LCMS. Multiplicative concentrations spanning orders of magnitude."),
    ("sqrt", "square_root", "y^(1/2)", "y >= 0", "Strongly recommended for moderately skewed concentrations (e.g. UK Nitrate) and taxa counts."),
    ("cbrt", "cube_root", "sign(y)*|y|^(1/3)", "All real y", "Wilson-Hilferty transformation for Gamma-distributed hydrological variables."),
    ("log1p", "log_1p", "ln(1 + y)", "y >= 0", "Accommodates true zeros (y=0) without arbitrary LOD substitution."),
    ("log10", "log_10", "log10(y)", "y > 0", "Order-of-magnitude chemical concentrations (e.g. pathogen counts, conductivity)."),
    ("fourthrt", "qtrt", "y^(1/4)", "y >= 0", "Intermediate between square root and log; classic ecological standard."),
    ("boxcox", "bc, box_cox", "(y^lambda-1)/lam", "y > 0", "Analytically optimizes residual normality via Profile MLE lambda in [-2, 2]."),
    ("yeojohnson", "yj", "Yeo-Johnson(lam)", "All real y", "Extension of Box-Cox supporting zero and negative values."),
    ("reciprocal", "inv, inverse", "-1 / y", "y > 0", "Monotonic rate-dilution processes preserving positive spatial gradient direction."),
    ("asin", "arcsin, arcsine", "asin(sqrt(y/sc))", "0 <= y <= scale", "Angular transformation for bounded proportions and percentages."),
    ("logit", "logit", "ln(p / (1 - p))", "0 < y < scale", "S-shaped logistic transformation for bounded proportions."),
    ("none", "raw, identity", "y", "All real y", "Default for ecological indices (WHPT ASPT, LIFE index, RHS scores)."),
]

# Minimal programmatic defaults used ONLY if configs/ JSON files are absent on disk
MINIMAL_DEFAULT_NETWORKS = {
    "orn": {
        "name": "orn", "abbr": "orn", "label": "OS Open Rivers Network (ORN v2)", "short_label": "ORN v2",
        "gpkg_paths": ["gis/ORN_v2.gpkg"], "zip_fallbacks": ["gis/ORN_v2_GeoPackage.zip"],
        "layers": ["ORN"], "geom_col": "Shape", "crs_epsg": 27700,
        "topology": {"mode": "nodes", "id_col": "RivID", "from_node_col": "Fnode", "to_node_col": "Tnode",
                     "next_down_col": None, "catch_id_col": "CatchID", "rbd_id_col": "RBD_ID"},
        "attributes": {"dist_to_mouth_col": "Dist2Mth", "dist_to_mouth_unit": "m",
                       "afv_weight_col": "US_Accum", "afv_weight_unit": "m",
                       "shreve_col": "Shreve", "strahler_col": "Strahler",
                       "slope_col": None, "elevation_col": None, "river_name_cols": ["name1", "name2"]},
        "covariates": {"mode": "orn_100m_table", "sampled_table": "ORN_Sampled_100m",
                       "bedrock_lookup": "tblBedrockGeologyLookUp",
                       "superficial_lookup": "tblSuperficialGeologyLookUp", "rbd_lookup": "tblRiverBasinLookUp"},
        "core_basins": {"humber": 80, "thames": 3055, "severn": 2786}, "max_snap_dist_m": 250.0,
    }
}

MINIMAL_DEFAULT_SOURCES = {
    "rsn": {
        "name": "rsn", "abbr": "RSN", "label": "River Surveillance Network (RSN)",
        "description": "Environment Agency River Surveillance Network (1,382 baseline stations)",
        "data_dir": "data/NCEA_River_Surveillance_Network_ARD_Open",
        "sites_file": "RSN Water Quality Sample Points.csv",
        "eco_sites_file": "RSN Ecology Sample Points.csv",
        "id_col": "network_id", "easting_col": "ngr_easting", "northing_col": "ngr_northing", "crs_epsg": 27700,
        "datasets": {
            "wq": {"file": "RSN Water Quality Samples.csv", "id_col": "network_id",
                   "datetime_col": "sample_datetime", "det_col": "det_code", "det_desc_col": "det_desc",
                   "val_col": "meas_result", "cens_col": "meas_sign", "cens_flag": "<", "unit_col": "unit_short_desc"},
            "lcms": {"file": "RSN LCMS GCMS Samples.csv", "id_col": "network_id",
                     "datetime_col": "sample_datetime", "det_col": "cas_number", "det_desc_col": "compound_name",
                     "val_col": "meas_result", "cens_col": "less_than_flag", "cens_flag": "yes", "unit_col": "unit_short_desc"},
            "inverts": {"file": "RSN Macroinvertebrate Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
            "diatoms": {"file": "RSN Diatom Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
            "macrophytes": {"file": "RSN Macrophyte Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
            "rhs": {"file": "RSN River Habitat Surveys.csv", "id_col": "network_id", "datetime_col": "survey_date"},
        },
    },
    "ssn": {
        "name": "ssn", "abbr": "SSN", "label": "Small Streams Network (SSN)",
        "description": "Environment Agency Small Streams Network (313 dedicated headwater stations)",
        "data_dir": "data/NCEA_Small_Streams_Network_ARD_Open",
        "sites_file": "SSN Water Quality Sample Points.csv",
        "eco_sites_file": "SSN Ecology Sample Points.csv",
        "id_col": "network_id", "easting_col": "ngr_easting", "northing_col": "ngr_northing", "crs_epsg": 27700,
        "datasets": {
            "wq": {"file": "SSN Water Quality Samples.csv", "id_col": "network_id",
                   "datetime_col": "sample_datetime", "det_col": "det_code", "det_desc_col": "det_desc",
                   "val_col": "meas_result", "cens_col": "meas_sign", "cens_flag": "<", "unit_col": "unit_short_desc"},
            "inverts": {"file": "SSN Macroinvertebrate Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
            "diatoms": {"file": "SSN Diatom Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
            "macrophytes": {"file": "SSN Macrophyte Samples.csv", "id_col": "network_id", "datetime_col": "sample_datetime"},
        },
    },
    "composites": {
        "composites": {
            "rss": {"name": "rss", "abbr": "RSS", "label": "Combined River Surveillance & Small Streams (RSS)",
                    "description": "Integrated RSN and SSN", "sources": ["rsn", "ssn"]},
            "all": {"name": "all", "abbr": "ALL", "label": "All Registered Monitoring Networks",
                    "description": "Union of all active monitoring networks", "sources": ["rsn", "ssn"]},
        }
    },
}


# ==============================================================================
# 1. Directory Configuration & RuntimeConfig
# ==============================================================================

def _rel_join(base, *parts):
    p = os.path.join(base, *parts)
    if p.startswith("./") and len(p) > 2:
        p = p[2:]
    return p


def find_directory_config_file(base_dir=DEFAULT_BASE_DIR):
    """Returns canonical path to configs/directories.json."""
    return _rel_join(base_dir, "configs", "directories.json")


def ensure_directory_config(base_dir=DEFAULT_BASE_DIR, default_data="data", default_gis="gis", default_out="output"):
    """Creates configs/directories.json if it does not exist yet. Returns (cfg_file_path, was_created)."""
    cfg_file = find_directory_config_file(base_dir)
    created = False
    if not os.path.exists(cfg_file):
        os.makedirs(os.path.dirname(os.path.abspath(cfg_file)), exist_ok=True)
        default_cfg = {
            "data_dir": default_data,
            "gis_dir": default_gis,
            "output_dir": default_out,
        }
        with open(cfg_file, "w", encoding="utf-8") as f:
            json.dump(default_cfg, f, indent=2)
            f.write("\n")
        created = True
        print(f"   [Config] Created directory configuration file: {cfg_file}")
    return cfg_file, created


def load_directory_config(base_dir=DEFAULT_BASE_DIR, auto_create=False):
    """Loads directory paths from configs/directories.json."""
    cfg_file = find_directory_config_file(base_dir)
    if not os.path.exists(cfg_file) and auto_create:
        ensure_directory_config(base_dir)

    raw_cfg = {}
    if os.path.exists(cfg_file):
        try:
            with open(cfg_file, "r", encoding="utf-8") as f:
                raw_cfg = json.load(f)
        except Exception as e:
            print(f"   [Warning] Could not parse directory config {cfg_file}: {e}")

    def _resolve(val, default):
        p = str(val).strip() if (val is not None and str(val).strip()) else default
        p = os.path.expanduser(p)
        if not os.path.isabs(p):
            p = _rel_join(base_dir, p)
        return p

    data_dir = _resolve(raw_cfg.get("data_dir"), "data")
    gis_dir = _resolve(raw_cfg.get("gis_dir"), "gis")
    output_dir = _resolve(raw_cfg.get("output_dir"), "output")
    proc_val = raw_cfg.get("processed_dir")
    processed_dir = _resolve(proc_val, _rel_join(data_dir, "processed")) if proc_val else _rel_join(data_dir, "processed")

    return {
        "config_file": cfg_file if os.path.exists(cfg_file) else None,
        "data_dir": data_dir,
        "gis_dir": gis_dir,
        "output_dir": output_dir,
        "processed_dir": processed_dir,
    }


def load_determinand_catalog(base_dir=DEFAULT_BASE_DIR):
    """Loads known determinands metadata from configs/determinands.json."""
    cat_path = _rel_join(base_dir, "configs", "determinands.json")
    if not os.path.isfile(cat_path):
        return {}
    try:
        with open(cat_path, "r", encoding="utf-8") as f:
            raw = json.load(f)
        return {str(k).lower(): v for k, v in raw.items() if isinstance(v, dict)}
    except Exception as e:
        print(f"   [Warning] Could not parse determinands catalog {cat_path}: {e}")
        return {}


@dataclass
class RuntimeConfig:
    """Explicit runtime path configuration passed through preprocessing functions."""
    base_dir: str
    gis_dir: str
    data_dir: str
    processed_dir: str
    output_dir: str
    rsn_dir: str
    ssn_dir: str
    config_dir: str
    source_config_dir: str


def build_runtime_config(base_dir=DEFAULT_BASE_DIR, gis_dir=None, data_dir=None, out_dir=None, create_out_dir=False):
    """Builds a RuntimeConfig from CLI args, environment variables, and configs/directories.json."""
    dir_cfg = load_directory_config(base_dir, auto_create=False)
    eff_gis = gis_dir or os.environ.get("S3N_GIS_DIR") or dir_cfg["gis_dir"]
    eff_data = data_dir or os.environ.get("S3N_DATA_DIR") or dir_cfg["data_dir"]

    rsn_dir = _rel_join(eff_data, "NCEA_River_Surveillance_Network_ARD_Open")
    ssn_dir = _rel_join(eff_data, "NCEA_Small_Streams_Network_ARD_Open")
    eff_proc = out_dir or os.environ.get("S3N_PROCESSED_DIR") or dir_cfg["processed_dir"]
    if create_out_dir:
        os.makedirs(eff_proc, exist_ok=True)

    return RuntimeConfig(
        base_dir=base_dir,
        gis_dir=eff_gis,
        data_dir=eff_data,
        processed_dir=eff_proc,
        output_dir=dir_cfg["output_dir"],
        rsn_dir=rsn_dir,
        ssn_dir=ssn_dir,
        config_dir=_rel_join(base_dir, "configs", "networks"),
        source_config_dir=_rel_join(base_dir, "configs", "sources"),
    )


def ensure_network_configs(base_dir=DEFAULT_BASE_DIR):
    """Ensures configs/networks exists and has at least the default ORN config if missing."""
    net_dir = _rel_join(base_dir, "configs", "networks")
    os.makedirs(net_dir, exist_ok=True)
    created = []
    for net_name, net_cfg in MINIMAL_DEFAULT_NETWORKS.items():
        dst = os.path.join(net_dir, f"{net_name}.json")
        if not os.path.exists(dst):
            try:
                with open(dst, "w", encoding="utf-8") as f:
                    json.dump(net_cfg, f, indent=2)
                    f.write("\n")
                created.append(dst)
            except Exception as e:
                print(f"   [Warning] Failed to write network config {dst}: {e}")
    return created


def ensure_source_configs(base_dir=DEFAULT_BASE_DIR):
    """Ensures configs/sources exists and has core default configs (RSN, SSN, composites) if missing."""
    src_dir = _rel_join(base_dir, "configs", "sources")
    os.makedirs(src_dir, exist_ok=True)
    created = []
    for src_name, src_cfg in MINIMAL_DEFAULT_SOURCES.items():
        dst = os.path.join(src_dir, f"{src_name}.json")
        if not os.path.exists(dst):
            try:
                with open(dst, "w", encoding="utf-8") as f:
                    json.dump(src_cfg, f, indent=2)
                    f.write("\n")
                created.append(dst)
            except Exception as e:
                print(f"   [Warning] Failed to write source config {dst}: {e}")
    return created


_ECOLOGICAL_DETERMINAND_DEFAULTS = {
    "whpt_aspt": ("WHPT ASPT", "WHPT Average Score Per Taxon", "ASPT"),
    "whpt_ntaxa": ("WHPT N-Taxa", "WHPT Taxon Richness (N-Taxa)", "taxa"),
    "life_family_index": ("LIFE Index", "LIFE Flow Velocity Stress Index (Family)", "LIFE"),
    "psi_family_score": ("PSI Score", "PSI Fine Sediment Sensitivity Score", "%"),
    "epsi_family_score": ("ePSI Score", "Empirically Weighted PSI Sediment Score", "%"),
    "cci": ("CCI", "Community Conservation Index", "CCI"),
    "dehli": ("DEHLI", "DEHLI Drought Stress Index", "DEHLI"),
    "armi_score": ("ARMI Score", "Anglers' Riverfly Monitoring Score", "score"),
    "percent_organic_tolerant": ("Organic-Tolerant Diatoms", "Organic-Tolerant Diatom Relative Abundance", "%"),
    "percent_motile": ("Motile Diatoms", "Motile Diatom Guild Relative Abundance", "%"),
    "rmni": ("RMNI", "River Macrophyte Nutrient Index (RMNI)", "RMNI"),
    "rn_a_taxa": ("Macrophyte Taxa", "Macrophyte Aquatic Taxon Richness", "taxa"),
    "hqa_score": ("RHS HQA", "River Habitat Survey Habitat Quality (HQA)", "HQA"),
    "hms_score": ("RHS HMS", "River Habitat Survey Modification Score (HMS)", "HMS"),
}

_FALLBACK_CORE_DETERMINANDS = {
    "0117": {"short": "Nitrate", "name": "Nitrate as N", "unit": "mg/L"},
    "0180": {"short": "Orthophosphate", "name": "Orthophosphate, reactive as P", "unit": "mg/L"},
    "0061": {"short": "pH", "name": "pH : In Situ", "unit": "pH units"},
}


def _normalize_det_unit(unit_str):
    u = str(unit_str or "").strip()
    if not u or u.lower() in ("na", "n/a", "none", "null"):
        return ""
    return re.sub(r"/l\b", "/L", u)


def _derive_short_det_name(full_name):
    s = str(full_name or "").strip()
    if not s:
        return ""
    s = re.sub(r"\s*:\s*In Situ\s*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r",\s*reactive as P\s*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+as\s+[A-Za-z0-9]+$", "", s)
    if "," in s and len(s) > 28:
        s = s.split(",", 1)[0].strip()
    return s


def _scan_wq_determinands_csv(csv_path, catalog):
    if not os.path.isfile(csv_path):
        return
    try:
        with open(csv_path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames:
                return
            fields = {c.strip(): c for c in reader.fieldnames}
            code_col = fields.get("det_code") or fields.get("cas_number")
            name_col = fields.get("det_desc") or fields.get("compound_name")
            unit_col = fields.get("unit_short_desc")
            if not code_col or not name_col:
                return
            for row in reader:
                code = str(row.get(code_col) or "").strip()
                if not code or code in catalog:
                    continue
                name = str(row.get(name_col) or "").strip() or code
                unit = _normalize_det_unit(row.get(unit_col, "")) if unit_col else ""
                catalog[code] = {
                    "short": _derive_short_det_name(name) or code,
                    "name": name,
                    "unit": unit,
                }
    except Exception as e:
        print(f"   [Warning] Could not scan determinands from {csv_path}: {e}")


def ensure_determinands_config(base_dir=DEFAULT_BASE_DIR, data_dir=None, force_update=False):
    """Generates configs/determinands.json from available RSN/SSN data files if missing."""
    cfg_dir = _rel_join(base_dir, "configs")
    os.makedirs(cfg_dir, exist_ok=True)
    cat_path = os.path.join(cfg_dir, "determinands.json")
    if os.path.isfile(cat_path) and not force_update:
        return cat_path

    catalog = {}
    if os.path.isfile(cat_path):
        try:
            with open(cat_path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
            if isinstance(loaded, dict):
                catalog.update(loaded)
        except Exception:
            pass

    eff_data = data_dir or load_directory_config(base_dir, auto_create=False)["data_dir"]
    rsn_dir = _rel_join(eff_data, "NCEA_River_Surveillance_Network_ARD_Open")
    ssn_dir = _rel_join(eff_data, "NCEA_Small_Streams_Network_ARD_Open")

    for d_dir, prefix in ((rsn_dir, "RSN"), (ssn_dir, "SSN")):
        if not os.path.isdir(d_dir):
            continue
        _scan_wq_determinands_csv(os.path.join(d_dir, f"{prefix} Water Quality Samples.csv"), catalog)
        _scan_wq_determinands_csv(os.path.join(d_dir, f"{prefix} LCMS GCMS Samples.csv"), catalog)

    for key, (short_n, full_n, unit_n) in _ECOLOGICAL_DETERMINAND_DEFAULTS.items():
        if key not in catalog:
            catalog[key] = {"short": short_n, "name": full_n, "unit": unit_n}

    if not catalog:
        catalog.update(_FALLBACK_CORE_DETERMINANDS)

    sorted_catalog = {k: catalog[k] for k in sorted(catalog.keys())}
    try:
        with open(cat_path, "w", encoding="utf-8") as f:
            json.dump(sorted_catalog, f, indent=2)
            f.write("\n")
        print(f"   [Config] Wrote {len(sorted_catalog)} determinand entries to {cat_path}")
    except Exception as e:
        print(f"   [Warning] Failed to write determinands catalog {cat_path}: {e}")
    return cat_path


# ==============================================================================
# 2. Coordinate Reference System (CRS) & Projection Helpers
# ==============================================================================

def wgs84_to_bng(lat, lon):
    """
    Converts WGS84 decimal latitude and longitude into British National Grid
    (OSGB36 / EPSG:27700) Eastings and Northings in meters using the standard
    7-parameter Helmert datum shift (WGS84 -> OSGB36) followed by the
    Airy 1830 Transverse Mercator projection (pure Python, zero dependencies).
    """
    a_w, b_w = 6378137.0, 6356752.3141
    e2_w = 1.0 - (b_w * b_w) / (a_w * a_w)
    phi_w = math.radians(lat)
    lam_w = math.radians(lon)
    sin_phi_w = math.sin(phi_w)
    cos_phi_w = math.cos(phi_w)
    nu_w = a_w / math.sqrt(max(1e-12, 1.0 - e2_w * sin_phi_w * sin_phi_w))
    x_w = nu_w * cos_phi_w * math.cos(lam_w)
    y_w = nu_w * cos_phi_w * math.sin(lam_w)
    z_w = (1.0 - e2_w) * nu_w * sin_phi_w

    tx, ty, tz = -446.448, 125.157, -542.060
    s = 1.0 + 20.4894e-6
    sec2rad = math.pi / (180.0 * 3600.0)
    rx = -0.1502 * sec2rad
    ry = -0.2470 * sec2rad
    rz = -0.8421 * sec2rad
    x_o = tx + s * x_w - rz * y_w + ry * z_w
    y_o = ty + rz * x_w + s * y_w - rx * z_w
    z_o = tz - ry * x_w + rx * y_w + s * z_w

    a, b = 6377563.396, 6356256.909
    e2 = 1.0 - (b * b) / (a * a)
    p = math.hypot(x_o, y_o)
    lam = math.atan2(y_o, x_o)
    phi = math.atan2(z_o, p * (1.0 - e2))
    for _ in range(8):
        nu_i = a / math.sqrt(max(1e-12, 1.0 - e2 * (math.sin(phi) ** 2)))
        phi_new = math.atan2(z_o + e2 * nu_i * math.sin(phi), p)
        if abs(phi_new - phi) < 1e-11:
            phi = phi_new
            break
        phi = phi_new

    F0 = 0.9996012717
    lat0 = math.radians(49.0)
    lon0 = math.radians(-2.0)
    N0, E0 = -100000.0, 400000.0
    n = (a - b) / (a + b)

    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    nu = a * F0 / math.sqrt(max(1e-12, 1.0 - e2 * sin_phi * sin_phi))
    rho = a * F0 * (1.0 - e2) / (math.pow(max(1e-12, 1.0 - e2 * sin_phi * sin_phi), 1.5))
    eta2 = max(0.0, nu / rho - 1.0)

    M = b * F0 * (
        (1.0 + n + 1.25 * (n ** 2) + 1.25 * (n ** 3)) * (phi - lat0)
        - (3.0 * n + 3.0 * (n ** 2) + 2.625 * (n ** 3)) * math.sin(phi - lat0) * math.cos(phi + lat0)
        + (1.875 * (n ** 2) + 1.875 * (n ** 3)) * math.sin(2.0 * (phi - lat0)) * math.cos(2.0 * (phi + lat0))
        - (35.0 / 24.0 * (n ** 3)) * math.sin(3.0 * (phi - lat0)) * math.cos(3.0 * (phi + lat0))
    )

    I = M + N0
    II = (nu / 2.0) * sin_phi * cos_phi
    III = (nu / 24.0) * sin_phi * (cos_phi ** 3) * (5.0 - (math.tan(phi) ** 2) + 9.0 * eta2)
    IIIA = (nu / 720.0) * sin_phi * (cos_phi ** 5) * (61.0 - 58.0 * (math.tan(phi) ** 2) + (math.tan(phi) ** 4))
    IV = nu * cos_phi
    V = (nu / 6.0) * (cos_phi ** 3) * (nu / rho - (math.tan(phi) ** 2))
    VI = (nu / 120.0) * (cos_phi ** 5) * (5.0 - 18.0 * (math.tan(phi) ** 2) + (math.tan(phi) ** 4) + 14.0 * eta2 - 58.0 * (math.tan(phi) ** 2) * eta2)

    d_lam = lam - lon0
    N = I + II * (d_lam ** 2) + III * (d_lam ** 4) + IIIA * (d_lam ** 6)
    E = E0 + IV * d_lam + V * (d_lam ** 3) + VI * (d_lam ** 5)
    return round(E, 2), round(N, 2)


def transform_coords_batch(xy_list, from_epsg, to_epsg):
    """
    Transforms a list of (x, y) coordinates from `from_epsg` to `to_epsg`.
    Uses in-process `pyproj` or `osgeo.osr` when installed, falling back to `cs2cs`.
    """
    if not xy_list or int(from_epsg) == int(to_epsg):
        return list(xy_list)
    src_e, dst_e = int(from_epsg), int(to_epsg)

    if _PyprojTransformer is not None:
        key = ("pyproj", src_e, dst_e)
        tr = _TRANSFORMER_CACHE.get(key)
        if tr is None:
            tr = _PyprojTransformer.from_crs(f"EPSG:{src_e}", f"EPSG:{dst_e}", always_xy=True)
            _TRANSFORMER_CACHE[key] = tr
        xs = [p[0] for p in xy_list]
        ys = [p[1] for p in xy_list]
        txs, tys = tr.transform(xs, ys)
        return [(round(float(x), 3), round(float(y), 3)) for x, y in zip(txs, tys)]

    if _OsgeoOsr is not None:
        key = ("osr", src_e, dst_e)
        ct = _TRANSFORMER_CACHE.get(key)
        if ct is None:
            srs_in = _OsgeoOsr.SpatialReference()
            srs_in.ImportFromEPSG(src_e)
            srs_in.SetAxisMappingStrategy(_OsgeoOsr.OAMS_TRADITIONAL_GIS_ORDER)
            srs_out = _OsgeoOsr.SpatialReference()
            srs_out.ImportFromEPSG(dst_e)
            srs_out.SetAxisMappingStrategy(_OsgeoOsr.OAMS_TRADITIONAL_GIS_ORDER)
            ct = _OsgeoOsr.CoordinateTransformation(srs_in, srs_out)
            _TRANSFORMER_CACHE[key] = ct
        pts3 = ct.TransformPoints([(float(x), float(y), 0.0) for x, y in xy_list])
        return [(round(p[0], 3), round(p[1], 3)) for p in pts3]

    inp = "\n".join(f"{x:.8f} {y:.8f}" for x, y in xy_list)
    fmt = "%.8f" if dst_e in (4326, 4258, 4277) else "%.3f"
    cmd = ["cs2cs", "-f", fmt, f"+init=epsg:{src_e}", "+to", f"+init=epsg:{dst_e}"]
    proc = subprocess.run(cmd, input=inp, text=True, capture_output=True, check=True)
    out = []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out.append((float(parts[0]), float(parts[1])))
    if len(out) != len(xy_list):
        raise RuntimeError(f"cs2cs returned {len(out)} coordinates for {len(xy_list)} inputs.")
    return out


# ==============================================================================
# 3. String, Column, Date, and Source Registry Helpers
# ==============================================================================

_RE_NORM_NAME = re.compile(r"[^a-z0-9]")
_RE_TIME_PART = re.compile(r"(\d{1,2}):(\d{2})(?::(\d{2}))?\s*(am|pm)?", re.IGNORECASE)
_RE_TIME_STRIP = re.compile(r"\s+\d{1,2}:\d{2}(?::\d{2})?\s*(?:am|pm)?.*$", re.IGNORECASE)
_RE_DATE_ISO = re.compile(r"^(\d{4})[-/.](\d{1,2})[-/.](\d{1,2})")
_RE_DATE_SLASH = re.compile(r"^(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})")


def normalize_name(s):
    """Strips punctuation, underscores, and extra spaces for resilient matching."""
    if not s:
        return ""
    return _RE_NORM_NAME.sub("", str(s).lower())


def find_col(fieldnames, candidate_names):
    """Case-insensitive matching of column names against candidate variants."""
    if not fieldnames:
        return None
    cols_map = {str(c).lower().strip(): c for c in fieldnames if c}
    for cand in candidate_names:
        cand_low = cand.lower().strip()
        if cand_low in cols_map:
            return cols_map[cand_low]
    return None


def find_matching_column(fieldnames, candidate_aliases):
    """Exact case-insensitive, then normalized alphanumeric, then short substring match."""
    if not fieldnames or not candidate_aliases:
        return None
    lower_map = {f.lower().strip(): f for f in fieldnames}
    for cand in candidate_aliases:
        c_low = cand.lower().strip()
        if c_low in lower_map:
            return lower_map[c_low]
    norm_map = {normalize_name(f): f for f in fieldnames}
    for cand in candidate_aliases:
        c_norm = normalize_name(cand)
        if c_norm in norm_map:
            return norm_map[c_norm]
    for cand in candidate_aliases:
        c_norm = normalize_name(cand)
        if len(c_norm) >= 3:
            for f_norm, original_f in norm_map.items():
                if c_norm == f_norm or (c_norm in f_norm and len(f_norm) <= len(c_norm) + 4):
                    return original_f
    return None


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
    raise ValueError(f"Invalid --month '{m_str}'. Specify 1-12 or month name (e.g. 'jan', 'july').")


def compute_afv_reach_weight(afv_scheme, us_accum_km, shreve=1, elevation=0.0, is_aquifer=0):
    """Compute positive accumulated hydrological mass proxy for Tail-Up AFV splitting."""
    scheme = str(afv_scheme or "length").strip().lower().replace("-", "_")
    us_km = max(1e-4, float(us_accum_km or 0.0))
    if scheme in ("shreve", "order", "magnitude"):
        return max(1.0, float(shreve or 1.0))
    if scheme in ("hydro_yield", "yield", "hydro"):
        elev = max(0.0, float(elevation or 0.0))
        aq = max(0.0, float(is_aquifer or 0))
        return max(1e-4, us_km * (1.0 + 0.002 * elev) * (1.0 + 0.50 * aq))
    return us_km


@lru_cache(maxsize=16384)
def _parse_datetime_iso_cached(ds, ts):
    # Fast path for standard ISO-8601 strings: YYYY-MM-DD[T| ]HH:MM:SS
    if not ts and len(ds) >= 10 and ds[4] == "-" and ds[7] == "-" and ds[0:4].isdigit() and ds[5:7].isdigit() and ds[8:10].isdigit():
        if len(ds) == 10:
            return f"{ds}T12:00:00"
        if len(ds) >= 19 and (len(ds) == 19 or ds[19] in "Zz.+-") and ds[10] in ("T", " ") and ds[13] == ":" and ds[16] == ":" and ds[11:13].isdigit() and ds[14:16].isdigit() and ds[17:19].isdigit():
            return f"{ds[:10]}T{ds[11:19]}"

    time_part = None
    for src_t in (ts, ds):
        if time_part:
            break
        if src_t and src_t.lower() not in ("none", "null", "na", "n/a", ""):
            m_t = _RE_TIME_PART.search(src_t)
            if m_t:
                hr, mn, sc, ampm = m_t.groups()
                hr, mn = int(hr), int(mn)
                sc = int(sc) if sc else 0
                if ampm:
                    ampm_low = ampm.lower()
                    if ampm_low == "pm" and hr < 12:
                        hr += 12
                    elif ampm_low == "am" and hr == 12:
                        hr = 0
                time_part = f"{hr:02d}:{mn:02d}:{sc:02d}"

    date_clean = _RE_TIME_STRIP.sub("", ds).strip()
    m_iso = _RE_DATE_ISO.match(date_clean)
    if m_iso:
        y, m, d = int(m_iso.group(1)), int(m_iso.group(2)), int(m_iso.group(3))
        return f"{y:04d}-{m:02d}-{d:02d}T{time_part or '12:00:00'}"

    m_slash = _RE_DATE_SLASH.match(date_clean)
    if m_slash:
        p1, p2, y = int(m_slash.group(1)), int(m_slash.group(2)), int(m_slash.group(3))
        if p1 > 12:
            d, m = p1, p2
        elif p2 > 12:
            m, d = p1, p2
        else:
            d, m = p1, p2
        return f"{y:04d}-{m:02d}-{d:02d}T{time_part or '12:00:00'}"

    if len(ds) >= 10 and ds[4] == "-" and ds[7] == "-" and ds[0:4].isdigit() and ds[5:7].isdigit() and ds[8:10].isdigit():
        return ds.replace(" ", "T")
    return ds


def parse_datetime_iso(date_str, time_str=None):
    """Normalizes arbitrary date and time strings into ISO-8601 'YYYY-MM-DDTHH:MM:SS'."""
    if not date_str:
        return "2023-06-01T12:00:00"
    ds = str(date_str).strip()
    if not ds:
        return "2023-06-01T12:00:00"
    ts = str(time_str).strip() if time_str else ""
    return _parse_datetime_iso_cached(ds, ts)


def load_source_registry(base_dir=DEFAULT_BASE_DIR, data_dir=None):
    """Loads source configs from configs/sources/*.json and auto-discovers data/custom/*."""
    ensure_source_configs(base_dir)
    cfg_dir = os.path.join(base_dir, "configs", "sources")
    eff_data_dir = data_dir or os.path.join(base_dir, "data")
    sources = {
        "rsn": dict(MINIMAL_DEFAULT_SOURCES["rsn"]),
        "ssn": dict(MINIMAL_DEFAULT_SOURCES["ssn"]),
    }
    composites = dict(MINIMAL_DEFAULT_SOURCES["composites"]["composites"])

    if os.path.isdir(cfg_dir):
        for f in sorted(glob.glob(os.path.join(cfg_dir, "*.json"))):
            fname = os.path.basename(f)
            try:
                with open(f, "r", encoding="utf-8") as jf:
                    data = json.load(jf)
                if fname == "composites.json" or "composites" in data:
                    comps = data.get("composites", data)
                    for c_name, c_info in comps.items():
                        composites[c_name.lower()] = c_info
                else:
                    s_name = data.get("name", os.path.splitext(fname)[0]).lower()
                    sources[s_name] = data
            except Exception as e:
                print(f"Warning: Failed to parse source config {f}: {e}")

    cust_dir = os.path.join(eff_data_dir, "custom")
    if os.path.isdir(cust_dir):
        for sd in sorted(os.listdir(cust_dir)):
            full_sd = os.path.join(cust_dir, sd)
            if os.path.isdir(full_sd):
                s_key = sd.lower()
                if s_key not in sources:
                    sources[s_key] = {
                        "name": s_key, "abbr": s_key.upper(),
                        "label": f"Custom Dataset: {sd}",
                        "description": f"User custom dataset in data/custom/{sd}",
                        "data_dir": full_sd, "sites_file": "sites.csv",
                        "id_col": "site_id", "easting_col": "easting", "northing_col": "northing",
                        "crs_epsg": 27700,
                    }
    return sources, composites


def resolve_source_tag(source="rsn", registry=None, data_dir=None, base_dir=DEFAULT_BASE_DIR):
    """
    Resolves source tags, aliases, composite groups, and delimiter combinations
    (e.g., 'rsn+ssn', 'waterblitz-tobit+rsn').
    Returns (source_key, constituent_list, source_abbr, source_label).
    """
    if not source:
        return "rsn", ["rsn"], "RSN", "River Surveillance Network (RSN)"

    if registry is None:
        sources_reg, composites_reg = load_source_registry(base_dir=base_dir, data_dir=data_dir)
    else:
        sources_reg, composites_reg = registry
    eff_data_dir = data_dir or os.path.join(base_dir, "data")

    raw_s = str(source).lower().strip()
    tokens = [t.strip().lower() for t in re.split(r"[,+ ]+", raw_s) if t.strip()]
    if not tokens:
        tokens = ["rsn"]

    def _canon_token(tk):
        tk_under = tk.replace("-", "_")
        if tk in composites_reg or tk in sources_reg or os.path.isdir(os.path.join(eff_data_dir, "custom", tk)):
            return tk
        if tk_under in composites_reg or tk_under in sources_reg or os.path.isdir(os.path.join(eff_data_dir, "custom", tk_under)):
            return tk_under
        return tk_under

    if len(tokens) > 1:
        expanded = []
        for tk in tokens:
            cand = _canon_token(tk)
            if cand in composites_reg:
                expanded.extend(composites_reg[cand].get("sources", [cand]))
            elif cand in sources_reg and sources_reg[cand].get("alias_target"):
                expanded.append(sources_reg[cand]["alias_target"])
            elif f"{cand}_tobit" in sources_reg:
                expanded.append(f"{cand}_tobit")
            else:
                expanded.append(cand)
        unique_tokens = sorted(set(expanded))
        if set(unique_tokens) == {"rsn", "ssn"}:
            return "rss", ["rsn", "ssn"], "RSS", "Combined River Surveillance & Small Streams (RSS)"
        combo_key = "+".join(unique_tokens)
        combo_abbr = "+".join(sources_reg[t].get("abbr", t.upper()) if t in sources_reg else t.upper() for t in unique_tokens)
        combo_label = f"Combined Monitoring Sources ({combo_abbr})"
        return combo_key, unique_tokens, combo_abbr, combo_label

    s = _canon_token(tokens[0])
    if s in composites_reg:
        c_info = composites_reg[s]
        c_srcs = c_info.get("sources", ["rsn", "ssn"])
        return s, c_srcs, c_info.get("abbr", s.upper()), c_info.get("label", f"Composite: {s}")
    if s in sources_reg:
        s_info = sources_reg[s]
        if s_info.get("alias_target") and s_info["alias_target"] in sources_reg:
            t_key = s_info["alias_target"]
            t_info = sources_reg[t_key]
            return t_key, [t_key], t_info.get("abbr", t_key.upper()), t_info.get("label", f"Monitoring Source: {t_key}")
        return s, [s], s_info.get("abbr", s.upper()), s_info.get("label", f"Monitoring Source: {s}")

    for rm in ("tobit", "geometric", "midpoint", "lower", "upper"):
        cand_key = f"{s}_{rm}"
        if cand_key in sources_reg:
            t_info = sources_reg[cand_key]
            return cand_key, [cand_key], t_info.get("abbr", cand_key.upper()), t_info.get("label", f"Monitoring Source: {cand_key}")

    for rm in ("tobit", "geometric", "midpoint", "lower", "upper"):
        if s.endswith(f"_{rm}"):
            base_s = s[: -len(rm) - 1]
            has_s_dir = os.path.isdir(os.path.join(eff_data_dir, "custom", s))
            if base_s in sources_reg:
                b_info = sources_reg[base_s]
                b_label = b_info.get("label", base_s.title())
                if f"[{rm.title()}]" not in b_label:
                    b_label = f"{b_label} [{rm.title()}]"
                b_abbr = f"{b_info.get('abbr', base_s.upper())}_{rm.upper()[:3]}"
                return s, ([s] if has_s_dir else [base_s]), b_abbr, b_label
            if not has_s_dir and base_s and os.path.isdir(os.path.join(eff_data_dir, "custom", base_s)):
                return s, [base_s], f"{base_s.upper()}_{rm.upper()[:3]}", f"Custom Source: {base_s} [{rm.title()}]"

    clean_tag = re.sub(r"[^a-z0-9_]+", "_", s)
    return clean_tag, [clean_tag], clean_tag.upper(), f"Custom Source: {clean_tag}"


def detect_range_mode(source_key, sources_reg, custom_samples=None, custom_dir=None, cli_mode="auto", base_dir=DEFAULT_BASE_DIR, data_dir=None):
    """Detects active range-handling mode from CLI override, source metadata, or samples CSV."""
    if cli_mode:
        cm = cli_mode.lower().strip()
        if cm == "none":
            return None
        if cm != "auto":
            return cm
    _, constituents, _, _ = resolve_source_tag(source_key, registry=(sources_reg, {}), data_dir=data_dir, base_dir=base_dir)
    found_modes = []
    for tk in constituents:
        tk_mode = None
        if tk in sources_reg and sources_reg[tk].get("range_mode"):
            tk_mode = str(sources_reg[tk]["range_mode"]).lower()
        else:
            tk_low = tk.lower()
            for m in ("tobit", "geometric", "midpoint", "lower", "upper"):
                if tk_low.endswith(f"_{m}") or tk_low.endswith(f"-{m}"):
                    tk_mode = m
                    break
        if tk_mode and tk_mode not in found_modes:
            found_modes.append(tk_mode)
    if len(found_modes) > 1:
        print(f"   [Warning] Mixed range modes across combined sources ({', '.join(found_modes)}); using '{found_modes[0]}'.")
    if found_modes:
        return found_modes[0]

    cand_paths = []
    if custom_samples:
        cand_paths.append(custom_samples)
    elif custom_dir:
        cand_paths.append(os.path.join(custom_dir, "samples.csv"))
    for tk in constituents:
        if tk in sources_reg:
            s_dir = sources_reg[tk].get("data_dir", "")
            if s_dir:
                cand_paths.append(os.path.join(s_dir, "samples.csv"))
                cand_paths.append(os.path.join(base_dir, s_dir, "samples.csv"))

    for cand_samples in cand_paths:
        if cand_samples and os.path.exists(cand_samples):
            try:
                with open(cand_samples, "r", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f)
                    if reader.fieldnames and "meas_method" in reader.fieldnames:
                        for i, r in enumerate(reader):
                            if i > 500:
                                break
                            mm = (r.get("meas_method") or "").strip().lower()
                            for m in ("tobit", "geometric", "midpoint", "lower", "upper"):
                                if m in mm or (m == "midpoint" and "arithmetic" in mm):
                                    return m
            except Exception:
                pass
    return None


# ==============================================================================
# 4. Statistical Distributions, Tobit MLE, Range Strategies, ROS, Transforms
# ==============================================================================

def normal_cdf(x):
    """Standard normal cumulative distribution function Phi(x)."""
    if x < -8.0:
        return 1e-15
    if x > 8.0:
        return 1.0 - 1e-15
    return max(1e-15, min(1.0 - 1e-15, 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))))


def normal_pdf(x):
    """Standard normal probability density function phi(x)."""
    if abs(x) > 8.0:
        return 1e-15
    return (1.0 / math.sqrt(2.0 * math.pi)) * math.exp(-0.5 * x * x)


_RE_CENS = re.compile(r"^(<|<=|>|>=)\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?:\s*[a-zA-Z/%_°][a-zA-Z0-9/%_°.\-]*)?$")
_RE_RANGE = re.compile(r"^([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*(?:-|–|—|\bto\b)\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(?:\s*[a-zA-Z/%_°][a-zA-Z0-9/%_°.\-]*)?$", re.IGNORECASE)
_RE_EMBEDDED = re.compile(r"^([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*[a-zA-Z/%_]+$")
_MISSING_TOKENS = frozenset(("none", "null", "na", "n/a", "nan", "no_water_available", "-", ".", ""))


def parse_raw_measurement_token(raw_val):
    """Categorizes a measurement string into 'point', 'interval', 'left', 'right', or None."""
    if raw_val is None:
        return None
    s = str(raw_val).strip()
    if not s or s.lower() in _MISSING_TOKENS:
        return None

    # Fast path for plain numeric strings (avoids running regexes on the majority of cells)
    if s[0] in "0123456789+-.":
        try:
            val = float(s)
            if not math.isnan(val) and not math.isinf(val):
                return {"type": "point", "a": val, "b": val, "val": val, "sign": "=", "raw": s}
        except ValueError:
            pass

    m_cens = _RE_CENS.match(s)
    if m_cens:
        op, num = m_cens.groups()
        try:
            val = float(num)
            if "<" in op:
                return {"type": "left", "a": 0.0, "b": val, "val": val, "sign": "<", "raw": s}
            return {"type": "right", "a": val, "b": None, "val": val, "sign": ">", "raw": s}
        except ValueError:
            return None

    m_range = _RE_RANGE.match(s)
    if m_range:
        try:
            a, b = float(m_range.group(1)), float(m_range.group(2))
            if a > b:
                a, b = b, a
            return {"type": "interval", "a": a, "b": b, "val": 0.5 * (a + b), "sign": "=", "raw": s}
        except ValueError:
            return None

    m_embedded = _RE_EMBEDDED.match(s)
    if m_embedded:
        try:
            val = float(m_embedded.group(1))
            return {"type": "point", "a": val, "b": val, "val": val, "sign": "=", "raw": s}
        except ValueError:
            pass
    return None


class IntervalTobitModel:
    """Interval-censored lognormal MLE for exact points, bounded intervals, and censored values."""

    def __init__(self, name="Determinand"):
        self.name = name
        self.mu = 0.0
        self.sigma = 1.0
        self.fitted = False
        self.n_total = 0
        self.n_points = 0
        self.n_intervals = 0
        self.n_left = 0
        self.n_right = 0
        self.cache = {}

    def fit_from_token_counts(self, token_counts):
        if not token_counts:
            return
        sample_log_vals = []
        for (rec_type, a, b), cnt in token_counts.items():
            self.n_total += cnt
            if rec_type == "point":
                self.n_points += cnt
                if a and a > 0:
                    sample_log_vals.extend([math.log(a)] * min(cnt, 200))
            elif rec_type == "interval":
                self.n_intervals += cnt
                if a is not None and a <= 0 and b and b > 0:
                    sample_log_vals.extend([math.log(b) - 0.5] * min(cnt, 200))
                elif a and b and a > 0 and b > 0:
                    sample_log_vals.extend([0.5 * (math.log(a) + math.log(b))] * min(cnt, 200))
            elif rec_type == "left":
                self.n_left += cnt
                if b and b > 0:
                    sample_log_vals.extend([math.log(b) - 0.5] * min(cnt, 200))
            elif rec_type == "right":
                self.n_right += cnt
                if a and a > 0:
                    sample_log_vals.extend([math.log(a) + 0.5] * min(cnt, 200))

        if not sample_log_vals:
            self.mu, self.sigma, self.fitted = 0.0, 1.0, True
            return

        sample_log_vals.sort()
        n = len(sample_log_vals)
        init_m = sample_log_vals[n // 2]
        q25 = sample_log_vals[int(0.25 * n)]
        q75 = sample_log_vals[int(0.75 * n)]
        init_s = max(0.15, (q75 - q25) / 1.349) if q75 > q25 else 1.0

        half_log_2pi = 0.5 * math.log(2.0 * math.pi)
        # Precompute log-bounds once outside the ~700-iteration grid/pattern search
        precomp_terms = []
        for (rec_type, a, b), cnt in token_counts.items():
            if rec_type == "point" and a and a > 0:
                precomp_terms.append((0, math.log(a), 0.0, cnt))
            elif (rec_type == "left" or (rec_type == "interval" and (a is None or a <= 0))) and b and b > 0:
                precomp_terms.append((1, 0.0, math.log(b), cnt))
            elif rec_type == "right" and a and a > 0:
                precomp_terms.append((2, math.log(a), 0.0, cnt))
            elif rec_type == "interval" and a and b and a > 0 and b > 0:
                precomp_terms.append((3, math.log(a), math.log(b), cnt))

        def log_likelihood(m, s):
            s = max(0.05, s)
            inv_s = 1.0 / s
            const_pt = -half_log_2pi - math.log(s)
            total = 0.0
            for kind, la, lb, cnt in precomp_terms:
                if kind == 0:
                    z = (la - m) * inv_s
                    total += cnt * (const_pt - la - 0.5 * z * z)
                elif kind == 1:
                    total += cnt * math.log(max(1e-12, normal_cdf((lb - m) * inv_s)))
                elif kind == 2:
                    total += cnt * math.log(max(1e-12, 1.0 - normal_cdf((la - m) * inv_s)))
                else:
                    total += cnt * math.log(max(1e-12, normal_cdf((lb - m) * inv_s) - normal_cdf((la - m) * inv_s)))
            return total

        best_ll = -1e20
        best_m, best_s = init_m, init_s
        m_span = max(1.2, 2.2 * init_s)
        m_cands = [init_m + x * (m_span / 15.0) for x in range(-15, 16)]
        s_cands = [max(0.1, init_s * math.exp(x * 0.15)) for x in range(-10, 11)]
        for mc in m_cands:
            for sc in s_cands:
                ll = log_likelihood(mc, sc)
                if ll > best_ll:
                    best_ll, best_m, best_s = ll, mc, sc

        step_m = (m_cands[1] - m_cands[0]) * 0.5
        step_s = max(0.02, best_s * 0.075)
        for _ in range(12):
            improved = False
            for dm in (-step_m, step_m):
                ll = log_likelihood(best_m + dm, best_s)
                if ll > best_ll:
                    best_ll, best_m, improved = ll, best_m + dm, True
            for ds in (-step_s, step_s):
                cand_s = max(0.05, best_s + ds)
                ll = log_likelihood(best_m, cand_s)
                if ll > best_ll:
                    best_ll, best_s, improved = ll, cand_s, True
            if not improved:
                step_m *= 0.5
                step_s *= 0.5

        self.mu, self.sigma, self.fitted = best_m, best_s, True

    def predict(self, rec_type, a, b, val=None):
        cache_key = (rec_type, a, b, val)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if rec_type == "point":
            res = (val, "=", val, val, 0.0, "exact_point")
            self.cache[cache_key] = res
            return res
        m, s = self.mu, self.sigma
        if rec_type == "left":
            b_pos = max(1e-6, b) if b is not None else 1e-6
            beta = (math.log(b_pos) - m) / s
            exp_log_y = m - s * (normal_pdf(beta) / max(1e-12, normal_cdf(beta)))
            res = (round(math.exp(exp_log_y), 5), "<", 0.0, b, round(s / math.sqrt(12.0), 4), "tobit_interval")
        elif rec_type == "right":
            a_pos = max(1e-6, a) if a is not None else 1e-6
            alpha = (math.log(a_pos) - m) / s
            exp_log_y = m + s * (normal_pdf(alpha) / max(1e-12, 1.0 - normal_cdf(alpha)))
            res = (round(math.exp(exp_log_y), 5), ">", a, None, round(s / math.sqrt(12.0), 4), "tobit_interval")
        elif rec_type == "interval":
            if a is None or a <= 0:
                b_pos = max(1e-6, b) if b is not None else 1e-6
                beta = (math.log(b_pos) - m) / s
                exp_log_y = m - s * (normal_pdf(beta) / max(1e-12, normal_cdf(beta)))
                res = (round(math.exp(exp_log_y), 5), "=", a if a is not None else 0.0, b, round(s / math.sqrt(12.0), 4), "tobit_interval")
            else:
                a_pos = max(1e-6, a)
                b_pos = max(a_pos + 1e-9, b) if b is not None else a_pos + 1e-6
                alpha = (math.log(a_pos) - m) / s
                beta = (math.log(b_pos) - m) / s
                denom = normal_cdf(beta) - normal_cdf(alpha)
                unc = (math.log(b_pos) - math.log(a_pos)) / math.sqrt(12.0)
                if denom < 1e-10:
                    res = (round(math.sqrt(a_pos * b_pos), 5), "=", a, b, round(unc, 4), "geometric_midpoint")
                else:
                    exp_log_y = m + s * ((normal_pdf(alpha) - normal_pdf(beta)) / denom)
                    res = (round(math.exp(exp_log_y), 5), "=", a, b, round(unc, 4), "tobit_interval")
        else:
            res = (val or 0.0, "=", a or 0.0, b or 0.0, 0.0, "unknown")
        self.cache[cache_key] = res
        return res


def compute_alternative_range_result(token, range_mode="geometric"):
    """Computes representative value & uncertainty for geometric/midpoint/lower/upper."""
    rec_type, a, b, val = token["type"], token.get("a"), token.get("b"), token.get("val")
    if rec_type == "point":
        return val, "=", val, val, 0.0, "exact_point"
    if rec_type == "left":
        return b, "<", 0.0, b, round(0.5 * max(0.0, b) / math.sqrt(12.0), 4), f"{range_mode}_left"
    if rec_type == "right":
        return a, ">", a, None, round(0.5 * max(0.0, a) / math.sqrt(12.0), 4), f"{range_mode}_right"
    if rec_type == "interval":
        if range_mode == "geometric":
            a_pos = max(1e-6, a) if a is not None else 1e-6
            b_pos = max(a_pos + 1e-9, b) if b is not None else a_pos + 1e-6
            return (round(math.sqrt(a_pos * b_pos), 5), "=", a, b,
                    round((math.log(b_pos) - math.log(a_pos)) / math.sqrt(12.0), 4), "geometric_midpoint")
        unc = round((b - a) / math.sqrt(12.0), 4)
        if range_mode == "midpoint":
            return round(0.5 * (a + b), 5), "=", a, b, unc, "arithmetic_midpoint"
        if range_mode == "lower":
            return a, "=", a, b, unc, "lower_bound"
        if range_mode == "upper":
            return b, "=", a, b, unc, "upper_bound"
    return val or 0.0, "=", a or 0.0, b or 0.0, 0.0, "unknown"


def extract_unit_from_string(s):
    m = re.search(r"\(([^)]+)\)|\[([^\]]+)\]", str(s))
    return (m.group(1) or m.group(2)).strip() if m else None


def estimate_boxcox_lambda(vals, lower=-2.0, upper=2.0):
    pos = [v for v in vals if v is not None and v > 0]
    n = len(pos)
    if n < 10:
        return 0.5
    log_vals = [math.log(v) for v in pos]
    sum_log = sum(log_vals)

    def ll(lam):
        if abs(lam) < 1e-4:
            z = log_vals
        else:
            inv_lam = 1.0 / lam
            z = [(math.pow(v, lam) - 1.0) * inv_lam for v in pos]
        mz = sum(z) / n
        vz = sum((x - mz) ** 2 for x in z) / n
        return -1e12 if vz < 1e-12 else (-0.5 * n * math.log(vz) + (lam - 1.0) * sum_log)

    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lower, upper
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = ll(c), ll(d)
    while abs(b - a) > 1e-4:
        if fc > fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = ll(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = ll(d)
    return round((b + a) / 2.0, 4)


def estimate_yeojohnson_lambda(vals, lower=-2.0, upper=2.0):
    clean = [v for v in vals if v is not None and not math.isnan(v)]
    n = len(clean)
    if n < 10:
        return 0.5
    pos_vals = [v + 1.0 for v in clean if v >= 0]
    neg_vals = [-v + 1.0 for v in clean if v < 0]
    sum_log_pos = sum(math.log(v) for v in pos_vals)
    sum_log_neg = sum(math.log(v) for v in neg_vals)

    def ll(lam):
        z = []
        if abs(lam) < 1e-4:
            z.extend(math.log(v) for v in pos_vals)
        else:
            inv_lam = 1.0 / lam
            z.extend((math.pow(v, lam) - 1.0) * inv_lam for v in pos_vals)
        lam2 = 2.0 - lam
        if abs(lam2) < 1e-4:
            z.extend(-math.log(v) for v in neg_vals)
        else:
            inv_lam2 = 1.0 / lam2
            z.extend(-((math.pow(v, lam2) - 1.0) * inv_lam2) for v in neg_vals)
        jac = (lam - 1.0) * sum_log_pos + (1.0 - lam) * sum_log_neg
        mz = sum(z) / n
        vz = sum((x - mz) ** 2 for x in z) / n
        return -1e12 if vz < 1e-12 else (-0.5 * n * math.log(vz) + jac)

    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lower, upper
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = ll(c), ll(d)
    while abs(b - a) > 1e-4:
        if fc > fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = ll(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = ll(d)
    return round((b + a) / 2.0, 4)


def apply_forward_transform(val, ttype="log", lam=None, scale=1.0):
    if val is None:
        return None
    try:
        val = float(val)
    except (ValueError, TypeError):
        return None
    if math.isnan(val):
        return None
    t = (ttype or "log").lower().strip()
    if t in ("log", "ln"):
        return math.log(max(1e-9, val)) if val > 0 else None
    if t in ("log10", "log_10"):
        return math.log10(max(1e-9, val)) if val > 0 else None
    if t in ("log1p", "log_1p"):
        return math.log1p(max(0.0, val)) if val >= 0 else None
    if t in ("sqrt", "square_root"):
        return math.sqrt(max(0.0, val)) if val >= 0 else None
    if t in ("cbrt", "cube_root"):
        return math.copysign(math.pow(abs(val), 1.0 / 3.0), val)
    if t in ("fourthrt", "qtrt"):
        return math.pow(max(0.0, val), 0.25) if val >= 0 else None
    if t in ("reciprocal", "inv", "inverse"):
        return -1.0 / max(1e-6, val) if val > 0 else None
    if t in ("boxcox", "box_cox", "bc"):
        if val <= 0:
            return None
        lv = lam if lam is not None else 0.5
        return math.log(val) if abs(lv) < 1e-4 else (math.pow(val, lv) - 1.0) / lv
    if t in ("yeojohnson", "yeo_johnson", "yj"):
        lv = lam if lam is not None else 0.5
        if val >= 0:
            return math.log(val + 1.0) if abs(lv) < 1e-4 else (math.pow(val + 1.0, lv) - 1.0) / lv
        return -math.log(-val + 1.0) if abs(2.0 - lv) < 1e-4 else -((math.pow(-val + 1.0, 2.0 - lv) - 1.0) / (2.0 - lv))
    if t in ("asin", "arcsin", "arcsine", "asin_sqrt"):
        return math.asin(math.sqrt(min(1.0, max(0.0, val / scale))))
    if t == "logit":
        p = min(1.0 - 1e-4, max(1e-4, val / scale))
        return math.log(p / (1.0 - p))
    return float(val)


def print_available_transforms():
    print("\nAvailable Data Transformations (--transform=<type> [--transform-lambda=<val>]):\n")
    print(f"{'Type':<12} | {'Short Aliases':<16} | {'Formula z = f(y)':<18} | {'Valid Domain':<16} | {'Best Applied When'}")
    print("-" * 125)
    for t_type, alias, formula, domain, note in TRANSFORM_CATALOG:
        print(f"{t_type:<12} | {alias:<16} | {formula:<18} | {domain:<16} | {note}")
    print()


_STD_NORM = statistics.NormalDist()


def compute_ros_plotting_positions(records):
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

    J = len(c_limits)
    A = [[] for _ in range(J + 1)]
    C = [[] for _ in range(J + 1)]
    limit_to_j = {lim: idx + 1 for idx, lim in enumerate(c_limits)}

    for idx, (v, is_c) in enumerate(records):
        if is_c:
            C[limit_to_j[v]].append(idx)
        else:
            # Find largest j in 1..J with c_limits[j - 1] <= v (or 0 if v < c_limits[0])
            j = bisect.bisect_right(c_limits, v)
            A[j].append(idx)

    for j in range(J + 1):
        A[j].sort(key=lambda i: records[i][0])

    # Exact O(J) computation of B_j = #{censored <= d_j} + #{uncensored < d_j}
    # Note: A[k] for k < j contains all uncensored v < d_j; C[k] for k <= j contains all censored v <= d_j.
    pref_A = [0] * (J + 2)
    pref_C = [0] * (J + 2)
    for k in range(J + 1):
        pref_A[k + 1] = pref_A[k] + len(A[k])
        pref_C[k + 1] = pref_C[k] + len(C[k])

    p = [0.0] * (J + 2)
    for j in range(J, 0, -1):
        A_j = len(A[j])
        B_j = pref_A[j] + pref_C[j + 1]
        p[j] = p[j + 1] + (1.0 - p[j + 1]) * (A_j / (A_j + B_j)) if (A_j + B_j > 0) else p[j + 1]
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
    mz, ml = statistics.mean(zs), statistics.mean(ly)
    vz = sum((z - mz) ** 2 for z in zs)
    if vz < 1e-9:
        return None
    slope = max(0.05, sum((z - mz) * (l - ml) for z, l in zip(zs, ly)) / vz)
    return (ml - slope * mz, slope)


def apply_ros_to_records(records, fallback_params=None):
    if not any(is_c for _, is_c in records):
        return [v for v, _ in records]
    pp = compute_ros_plotting_positions(records)
    n_unc = sum(1 for v, is_c in records if not is_c and v > 0)
    params = fit_ros_ols(records, pp) if (n_unc >= 3 and (n_unc / len(records)) >= 0.20) else None
    if params is None:
        params = fallback_params
    if params is None:
        return [0.5 * v if is_c else v for v, is_c in records]
    mu_hat, sigma_hat = params
    return [
        min(v, max(1e-6, math.exp(mu_hat + sigma_hat * _STD_NORM.inv_cdf(pp[i])))) if is_c else v
        for i, (v, is_c) in enumerate(records)
    ]


def impute_censored_by_site(raw_by_site, method="half_lod"):
    if method == "half_lod":
        return {nid: [0.5 * v if is_c else v for v, is_c in recs] for nid, recs in raw_by_site.items()}
    flat_nids, flat_recs = [], []
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
    if method == "ros_site":
        return {nid: apply_ros_to_records(recs, fallback_params=glob_params) for nid, recs in raw_by_site.items()}
    raise ValueError(f"Unknown censor_method: {method}")


# ==============================================================================
# 5. GeoPackage / WKB Geometry & Envelope Readers + Planar Metrics
# ==============================================================================

def _coerce_gpkg_blob(blob):
    if not blob:
        return None
    if isinstance(blob, str):
        s = blob.strip()
        if not s:
            return None
        try:
            return bytes.fromhex(s[2:] if s.lower().startswith("0x") else s)
        except ValueError:
            return None
    if isinstance(blob, memoryview):
        return blob.tobytes()
    return blob


def parse_gpkg_envelope(blob):
    """
    Fast O(1) extractor for the GeoPackage Binary (GPB) 2D envelope (minx, maxx, miny, maxy).
    Returns None if no GPB envelope header is present (e.g., raw WKB/GeoParquet).
    """
    blob = _coerce_gpkg_blob(blob)
    if not blob or len(blob) < 40:
        return None
    if blob[0] == 0x47 and blob[1] == 0x50:
        flags = blob[3]
        env_ind = (flags >> 1) & 0x07
        if env_ind in (1, 2, 3, 4):
            bo = "<" if (flags & 0x01) else ">"
            return struct.unpack_from(bo + "4d", blob, 8)
    return None


def parse_gpkg_geometry(blob):
    """Universal GPB / ISO WKB / EWKB parser returning (pts_xy, zs) for LineString & MultiLineString."""
    blob = _coerce_gpkg_blob(blob)
    if not blob or len(blob) < 9:
        return [], []
    if blob[0] == 0x47 and blob[1] == 0x50:
        env_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
        env_ind = (blob[3] >> 1) & 0x07
        offset = 8 + env_sizes.get(env_ind, 0)
    else:
        offset = 0
    if offset + 9 > len(blob):
        return [], []

    bo = "<" if blob[offset] == 1 else ">"
    wkb_type = struct.unpack_from(bo + "I", blob, offset + 1)[0]
    has_srid = (wkb_type & 0x20000000) != 0
    low16 = wkb_type & 0xFFFF
    has_z = (wkb_type & 0x80000000) != 0 or (1000 <= low16 < 4000 and (low16 // 1000) in (1, 3))
    has_m = (wkb_type & 0x40000000) != 0 or (1000 <= low16 < 4000 and (low16 // 1000) in (2, 3))
    ndim = 2 + (1 if has_z else 0) + (1 if has_m else 0)
    base_type = low16 % 1000
    hdr_len = 9 if has_srid else 5

    pts, zs = [], []
    if base_type == 5:
        num_geoms = struct.unpack_from(bo + "I", blob, offset + hdr_len)[0]
        offset += hdr_len + 4
        for _ in range(num_geoms):
            g_bo = "<" if blob[offset] == 1 else ">"
            g_type = struct.unpack_from(g_bo + "I", blob, offset + 1)[0]
            g_srid = (g_type & 0x20000000) != 0
            g_low16 = g_type & 0xFFFF
            g_z = (g_type & 0x80000000) != 0 or (1000 <= g_low16 < 4000 and (g_low16 // 1000) in (1, 3))
            g_m = (g_type & 0x40000000) != 0 or (1000 <= g_low16 < 4000 and (g_low16 // 1000) in (2, 3))
            g_dim = 2 + (1 if g_z else 0) + (1 if g_m else 0)
            g_hdr = 9 if g_srid else 5
            npts = struct.unpack_from(g_bo + "I", blob, offset + g_hdr)[0]
            offset += g_hdr + 4
            coords = struct.unpack_from(g_bo + f"{g_dim * npts}d", blob, offset)
            offset += 8 * g_dim * npts
            if g_dim == 2:
                pts.extend(zip(coords[0::2], coords[1::2]))
            else:
                for i in range(0, g_dim * npts, g_dim):
                    pts.append((coords[i], coords[i + 1]))
                    if g_z:
                        zv = coords[i + 2]
                        if not math.isnan(zv) and -500.0 < zv < 9000.0:
                            zs.append(zv)
    elif base_type == 2:
        npts = struct.unpack_from(bo + "I", blob, offset + hdr_len)[0]
        offset += hdr_len + 4
        coords = struct.unpack_from(bo + f"{ndim * npts}d", blob, offset)
        if ndim == 2:
            pts.extend(zip(coords[0::2], coords[1::2]))
        else:
            for i in range(0, ndim * npts, ndim):
                pts.append((coords[i], coords[i + 1]))
                if has_z:
                    zv = coords[i + 2]
                    if not math.isnan(zv) and -500.0 < zv < 9000.0:
                        zs.append(zv)
    return pts, zs


def parse_gpkg_multilinestring(blob):
    pts, _ = parse_gpkg_geometry(blob)
    return pts


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
        t = 0.0 if seg_len2 == 0.0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / seg_len2))
        projx, projy = x1 + t * dx, y1 + t * dy
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
    n_seg = len(poly) - 1
    seg_lens = [0.0] * n_seg
    tot = 0.0
    for i in range(n_seg):
        p1 = poly[i]
        p2 = poly[i + 1]
        d = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
        seg_lens[i] = d
        tot += d
    half, cur = tot * 0.5, 0.0
    mx, my = poly[0]
    for i, d in enumerate(seg_lens):
        if cur + d >= half and d > 0:
            t = (half - cur) / d
            p1, p2 = poly[i], poly[i + 1]
            mx, my = p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1])
            break
        cur += d
    straight = math.hypot(poly[-1][0] - poly[0][0], poly[-1][1] - poly[0][1])
    sinuosity = max(1.0, tot / straight) if straight > 1e-3 else 1.0
    return tot, sinuosity, (poly[0][0], poly[0][1], mx, my, poly[-1][0], poly[-1][1])
