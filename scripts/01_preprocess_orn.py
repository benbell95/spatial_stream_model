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
import shutil
import sqlite3
import statistics
import struct
import subprocess
import time
import zipfile
from collections import Counter, defaultdict

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPT_DIR)
try:
    BASE_DIR = "." if os.path.samefile(os.getcwd(), _REPO_ROOT) else os.path.relpath(_REPO_ROOT, os.getcwd())
except OSError:
    BASE_DIR = "."


def _rel_join(base, *parts):
    p = os.path.join(base, *parts)
    if p.startswith("./") and len(p) > 2:
        p = p[2:]
    return p


def find_directory_config_file(base_dir=BASE_DIR):
    """Finds existing directory configuration file in configs/ or root, or returns canonical path."""
    candidates = [
        _rel_join(base_dir, "configs", "directories.json"),
        _rel_join(base_dir, "configs", "paths.json"),
        _rel_join(base_dir, "config.json"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return candidates[0]


def ensure_directory_config(base_dir=BASE_DIR, default_data="data", default_gis="gis", default_out="output"):
    """
    Creates configs/directories.json automatically if no directory config file exists yet.
    Returns (cfg_file_path, was_created).
    """
    cfg_file = find_directory_config_file(base_dir)
    created = False
    if not os.path.exists(cfg_file):
        os.makedirs(os.path.dirname(os.path.abspath(cfg_file)), exist_ok=True)
        default_cfg = {
            "data_dir": default_data,
            "gis_dir": default_gis,
            "output_dir": default_out
        }
        with open(cfg_file, "w", encoding="utf-8") as f:
            json.dump(default_cfg, f, indent=2)
            f.write("\n")
        created = True
        print(f"   [Config] Created directory configuration file: {cfg_file}")
    return cfg_file, created


DEFAULT_NETWORK_CONFIGS = json.loads(r"""{
  "orn": {
    "name": "orn",
    "abbr": "orn",
    "label": "OS Open Rivers Network (ORN v2)",
    "short_label": "ORN v2",
    "gpkg_paths": [
      "gis/ORN_v2.gpkg"
    ],
    "zip_fallbacks": [
      "gis/ORN_v2_GeoPackage.zip"
    ],
    "layers": [
      "ORN"
    ],
    "geom_col": "Shape",
    "crs_epsg": 27700,
    "topology": {
      "mode": "nodes",
      "id_col": "RivID",
      "from_node_col": "Fnode",
      "to_node_col": "Tnode",
      "next_down_col": null,
      "catch_id_col": "CatchID",
      "rbd_id_col": "RBD_ID"
    },
    "attributes": {
      "dist_to_mouth_col": "Dist2Mth",
      "dist_to_mouth_unit": "m",
      "afv_weight_col": "US_Accum",
      "afv_weight_unit": "m",
      "shreve_col": "Shreve",
      "strahler_col": "Strahler",
      "slope_col": null,
      "elevation_col": null,
      "river_name_cols": [
        "name1",
        "name2"
      ]
    },
    "covariates": {
      "mode": "orn_100m_table",
      "sampled_table": "ORN_Sampled_100m",
      "bedrock_lookup": "tblBedrockGeologyLookUp",
      "superficial_lookup": "tblSuperficialGeologyLookUp",
      "rbd_lookup": "tblRiverBasinLookUp"
    },
    "core_basins": {
      "humber": 80,
      "thames": 3055,
      "severn": 2786
    },
    "max_snap_dist_m": 250.0
  }
}""")

OPTIONAL_NETWORK_CONFIGS = json.loads(r"""{
  "euhydro": {
    "name": "euhydro",
    "abbr": "euhydro",
    "label": "EU-Hydro v1.3 (Thames + Tweed)",
    "short_label": "EU-Hydro v1.3",
    "gpkg_paths": [
      "gis/EU_hydro_gpkg_eu/euhydro_thames_v013_GPKG/euhydro_thames_v013.gpkg",
      "gis/EU_hydro_gpkg_eu/euhydro_tweed_v013_GPKG/euhydro_tweed_v013.gpkg"
    ],
    "zip_fallbacks": [
      "gis/EU_hydro_gpkg_eu/euhydro_thames_v013_GPKG.zip",
      "gis/EU_hydro_gpkg_eu/euhydro_tweed_v013_GPKG.zip"
    ],
    "layers": [
      "River_Net_l",
      "Canals_l",
      "Ditches_l"
    ],
    "geom_col": "Shape",
    "crs_epsg": 3035,
    "topology": {
      "mode": "pointer",
      "id_col": "OBJECT_ID",
      "from_node_col": "FNODE",
      "to_node_col": "TNODE",
      "next_down_col": "NEXTDOWNID",
      "catch_id_col": "CatchID",
      "rbd_id_col": null
    },
    "attributes": {
      "dist_to_mouth_col": null,
      "dist_to_mouth_unit": "m",
      "afv_weight_col": "CUM_LEN",
      "afv_weight_unit": "m",
      "shreve_col": null,
      "strahler_col": "STRAHLER",
      "slope_col": "PENTE",
      "elevation_col": null,
      "river_name_cols": [
        "nameText",
        "nameTxtInt"
      ],
      "extra_numeric_cols": [
        "LONGPATH"
      ]
    },
    "nodes_table": {
      "table": "Nodes",
      "id_col": "OBJECT_ID",
      "elev_col": "ELEV",
      "dist_to_mouth_col": "LEN_TOM"
    },
    "covariates": {
      "mode": "geometry_and_spatial_join",
      "orn_geology_gpkg": "gis/ORN_v2.gpkg"
    },
    "core_basins": {
      "humber": 75,
      "thames": 306,
      "severn": 243
    },
    "max_snap_dist_m": 250.0
  },
  "euhydro_thames": {
    "name": "euhydro_thames",
    "abbr": "euhydro_thames",
    "label": "EU-Hydro v1.3 (Thames Network)",
    "short_label": "EU-Hydro (Thames)",
    "gpkg_paths": [
      "gis/EU_hydro_gpkg_eu/euhydro_thames_v013_GPKG/euhydro_thames_v013.gpkg"
    ],
    "zip_fallbacks": [
      "gis/EU_hydro_gpkg_eu/euhydro_thames_v013_GPKG.zip"
    ],
    "layers": [
      "River_Net_l",
      "Canals_l",
      "Ditches_l"
    ],
    "geom_col": "Shape",
    "crs_epsg": 3035,
    "topology": {
      "mode": "pointer",
      "id_col": "OBJECT_ID",
      "from_node_col": "FNODE",
      "to_node_col": "TNODE",
      "next_down_col": "NEXTDOWNID",
      "catch_id_col": "CatchID",
      "rbd_id_col": null
    },
    "attributes": {
      "dist_to_mouth_col": null,
      "dist_to_mouth_unit": "m",
      "afv_weight_col": "CUM_LEN",
      "afv_weight_unit": "m",
      "shreve_col": null,
      "strahler_col": "STRAHLER",
      "slope_col": "PENTE",
      "elevation_col": null,
      "river_name_cols": [
        "nameText",
        "nameTxtInt"
      ],
      "extra_numeric_cols": [
        "LONGPATH"
      ]
    },
    "nodes_table": {
      "table": "Nodes",
      "id_col": "OBJECT_ID",
      "elev_col": "ELEV",
      "dist_to_mouth_col": "LEN_TOM"
    },
    "covariates": {
      "mode": "geometry_and_spatial_join",
      "orn_geology_gpkg": "gis/ORN_v2.gpkg"
    },
    "core_basins": {
      "humber": 75,
      "thames": 306,
      "severn": 243
    },
    "max_snap_dist_m": 250.0
  },
  "euhydro_tweed": {
    "name": "euhydro_tweed",
    "abbr": "euhydro_tweed",
    "label": "EU-Hydro v1.3 (Tweed Network)",
    "short_label": "EU-Hydro (Tweed)",
    "gpkg_paths": [
      "gis/EU_hydro_gpkg_eu/euhydro_tweed_v013_GPKG/euhydro_tweed_v013.gpkg"
    ],
    "zip_fallbacks": [
      "gis/EU_hydro_gpkg_eu/euhydro_tweed_v013_GPKG.zip"
    ],
    "layers": [
      "River_Net_l",
      "Canals_l",
      "Ditches_l"
    ],
    "geom_col": "Shape",
    "crs_epsg": 3035,
    "topology": {
      "mode": "pointer",
      "id_col": "OBJECT_ID",
      "from_node_col": "FNODE",
      "to_node_col": "TNODE",
      "next_down_col": "NEXTDOWNID",
      "catch_id_col": "CatchID",
      "rbd_id_col": null
    },
    "attributes": {
      "dist_to_mouth_col": null,
      "dist_to_mouth_unit": "m",
      "afv_weight_col": "CUM_LEN",
      "afv_weight_unit": "m",
      "shreve_col": null,
      "strahler_col": "STRAHLER",
      "slope_col": "PENTE",
      "elevation_col": null,
      "river_name_cols": [
        "nameText",
        "nameTxtInt"
      ],
      "extra_numeric_cols": [
        "LONGPATH"
      ]
    },
    "nodes_table": {
      "table": "Nodes",
      "id_col": "OBJECT_ID",
      "elev_col": "ELEV",
      "dist_to_mouth_col": "LEN_TOM"
    },
    "covariates": {
      "mode": "geometry_and_spatial_join",
      "orn_geology_gpkg": "gis/ORN_v2.gpkg"
    },
    "core_basins": {
      "humber": null,
      "thames": null,
      "severn": null
    },
    "max_snap_dist_m": 250.0
  }
}""")

DEFAULT_SOURCE_CONFIGS = json.loads(r"""{
  "rsn": {
    "name": "rsn",
    "abbr": "RSN",
    "label": "River Surveillance Network (RSN)",
    "description": "Environment Agency River Surveillance Network (1,382 baseline stations across England's mainstems)",
    "data_dir": "data/NCEA_River_Surveillance_Network_ARD_Open",
    "sites_file": "RSN Water Quality Sample Points.csv",
    "eco_sites_file": "RSN Ecology Sample Points.csv",
    "id_col": "network_id",
    "easting_col": "ngr_easting",
    "northing_col": "ngr_northing",
    "crs_epsg": 27700,
    "datasets": {
      "wq": {
        "file": "RSN Water Quality Samples.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime",
        "det_col": "det_code",
        "det_desc_col": "det_desc",
        "val_col": "meas_result",
        "cens_col": "meas_sign",
        "cens_flag": "<",
        "unit_col": "unit_short_desc"
      },
      "lcms": {
        "file": "RSN LCMS GCMS Samples.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime",
        "det_col": "cas_number",
        "det_desc_col": "compound_name",
        "val_col": "meas_result",
        "cens_col": "less_than_flag",
        "cens_flag": "yes",
        "unit_col": "unit_short_desc"
      },
      "inverts": {
        "file": "RSN Macroinvertebrates.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime"
      },
      "diatoms": {
        "file": "RSN Diatoms.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime"
      },
      "macrophytes": {
        "file": "RSN Macrophytes.csv",
        "id_col": "network_id",
        "datetime_col": "survey_date"
      },
      "rhs": {
        "file": "RSN RHS.csv",
        "id_col": "network_id",
        "datetime_col": "survey_date"
      }
    }
  },
  "ssn": {
    "name": "ssn",
    "abbr": "SSN",
    "label": "Small Streams Network (SSN)",
    "description": "Environment Agency Small Streams Network (313 dedicated headwater stations across England)",
    "data_dir": "data/NCEA_Small_Streams_Network_ARD_Open",
    "sites_file": "SSN Water Quality Sample Points.csv",
    "eco_sites_file": "SSN Ecology Sample Points.csv",
    "id_col": "network_id",
    "easting_col": "ngr_easting",
    "northing_col": "ngr_northing",
    "crs_epsg": 27700,
    "datasets": {
      "wq": {
        "file": "SSN Water Quality Samples.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime",
        "det_col": "det_code",
        "det_desc_col": "det_desc",
        "val_col": "meas_result",
        "cens_col": "meas_sign",
        "cens_flag": "<",
        "unit_col": "unit_short_desc"
      },
      "inverts": {
        "file": "SSN Macroinvertebrates.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime"
      },
      "diatoms": {
        "file": "SSN Diatoms.csv",
        "id_col": "network_id",
        "datetime_col": "sample_datetime"
      },
      "macrophytes": {
        "file": "SSN Macrophytes.csv",
        "id_col": "network_id",
        "datetime_col": "survey_date"
      }
    }
  },
  "composites": {
    "composites": {
      "rss": {
        "name": "rss",
        "abbr": "RSS",
        "label": "Combined River Surveillance & Small Streams (RSS)",
        "description": "Integrated River Surveillance Network and Small Streams Network",
        "sources": [
          "rsn",
          "ssn"
        ]
      },
      "all": {
        "name": "all",
        "abbr": "ALL",
        "label": "All Registered Monitoring Networks",
        "description": "Union of all active monitoring networks",
        "sources": [
          "rsn",
          "ssn"
        ]
      }
    }
  }
}""")

OPTIONAL_SOURCE_CONFIGS = json.loads(r"""{
  "waterblitz": {
    "name": "waterblitz",
    "abbr": "WATERB",
    "label": "Waterblitz",
    "range_mode": "tobit",
    "alias_target": "waterblitz_tobit",
    "description": "Custom monitoring dataset waterblitz (defaulting to waterblitz_tobit)",
    "data_dir": "data/custom/waterblitz_tobit",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
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
  },
  "waterblitz_tobit": {
    "name": "waterblitz_tobit",
    "abbr": "WATERB_TOB",
    "label": "Waterblitz [Tobit]",
    "range_mode": "tobit",
    "description": "Custom monitoring dataset formatted by format_custom_data.py (TOBIT mode)",
    "data_dir": "data/custom/waterblitz_tobit",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
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
  },
  "waterblitz_midpoint": {
    "name": "waterblitz_midpoint",
    "abbr": "WATERB_MID",
    "label": "Waterblitz [Midpoint]",
    "range_mode": "midpoint",
    "description": "Custom monitoring dataset formatted by format_custom_data.py (MIDPOINT mode)",
    "data_dir": "data/custom/waterblitz_midpoint",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
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
  },
  "river_guardians": {
    "name": "river_guardians",
    "abbr": "RGCS",
    "label": "River Guardians Citizen Science Network",
    "description": "Volunteer water quality monitoring network collecting monthly nutrients and field chemistry",
    "data_dir": "data/custom/river_guardians",
    "sites_file": "volunteer_survey_points.csv",
    "id_col": "survey_site_id",
    "easting_col": "lon",
    "northing_col": "lat",
    "crs_epsg": 4326,
    "datasets": {
      "wq": {
        "file": "field_observations.csv",
        "id_col": "survey_site_id",
        "datetime_col": "observed_date",
        "det_col": "chemical_parameter",
        "det_desc_col": "parameter_name",
        "val_col": "concentration_ppm",
        "cens_col": "detection_flag",
        "cens_flag": "<",
        "unit_col": "units"
      }
    }
  },
  "rsn_gam": {
    "name": "rsn_gam",
    "abbr": "RSN_GAM",
    "label": "River Surveillance Network GAM Model Predictions",
    "short_label": "RSN GAM",
    "description": "Environment Agency River Surveillance Network GAM Model Estimates with Measurement Error",
    "data_dir": "data/custom/rsn_gam",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
    "meas_method": "gam",
    "datasets": {
      "wq": {
        "file": "samples.csv",
        "id_col": "site_id",
        "datetime_col": "sample_datetime",
        "val_col": "meas_result",
        "unc_col": "meas_uncertainty",
        "meth_col": "meas_method"
      }
    }
  },
  "rsn_gam_mean": {
    "name": "rsn_gam_mean",
    "abbr": "RSN_GAM_MEAN",
    "label": "River Surveillance Network Station Mean GAM Predictions",
    "short_label": "RSN GAM Mean",
    "description": "Environment Agency River Surveillance Network Station Mean GAM Model Estimates with Measurement Error",
    "data_dir": "data/custom/rsn_gam_mean",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
    "meas_method": "gam",
    "datasets": {
      "wq": {
        "file": "samples.csv",
        "id_col": "site_id",
        "datetime_col": "sample_datetime",
        "val_col": "meas_result",
        "unc_col": "meas_uncertainty",
        "meth_col": "meas_method"
      }
    }
  },
  "custom_example": {
    "name": "custom_example",
    "abbr": "CUSTOM",
    "label": "Custom River Monitoring Study",
    "description": "Example configuration template for registering a custom river monitoring dataset",
    "data_dir": "data/custom/example_study",
    "sites_file": "sites.csv",
    "id_col": "site_id",
    "easting_col": "easting",
    "northing_col": "northing",
    "crs_epsg": 27700,
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
}""")


def _network_files_exist(net_cfg, base_dir=BASE_DIR, gis_dir=None):
    if gis_dir is None:
        gis_dir = GIS_DIR
    gpkg_paths = net_cfg.get("gpkg_paths", [])
    parquet_paths = net_cfg.get("parquet_paths", [])
    network_paths = net_cfg.get("network_paths", [])
    zip_fallbacks = net_cfg.get("zip_fallbacks", [])

    all_targets = gpkg_paths + parquet_paths + network_paths + zip_fallbacks
    if not all_targets:
        return False

    for rel_p in all_targets:
        p1 = os.path.join(base_dir, rel_p)
        if os.path.exists(p1):
            return True
        rel_sub = rel_p[4:] if rel_p.startswith("gis/") else rel_p
        p2 = os.path.join(gis_dir, rel_sub)
        if os.path.exists(p2):
            return True
        p3 = os.path.join(gis_dir, rel_p)
        if os.path.exists(p3):
            return True
    return False


def _source_files_exist(src_cfg, base_dir=BASE_DIR, data_dir=None):
    if data_dir is None:
        data_dir = DATA_DIR
    raw_data_dir = src_cfg.get("data_dir")
    if not raw_data_dir:
        return False
    cands = [
        os.path.join(base_dir, raw_data_dir),
        os.path.join(data_dir, raw_data_dir),
        os.path.join(data_dir, raw_data_dir[5:]) if raw_data_dir.startswith("data/") else None,
        os.path.join(data_dir, "custom", os.path.basename(raw_data_dir)),
    ]
    for c in cands:
        if c and os.path.isdir(c):
            try:
                if len(os.listdir(c)) > 0:
                    return True
            except OSError:
                pass
    return False


def ensure_network_configs(base_dir=BASE_DIR):
    """
    Ensures that configs/networks directory exists and contains default network JSON configuration (ORN).
    Optional network configurations (e.g. EU-Hydro) are only created if their underlying GIS files exist in gis/.
    Returns list of paths created.
    """
    net_dir = _rel_join(base_dir, "configs", "networks")
    os.makedirs(net_dir, exist_ok=True)
    created = []

    # 1. Unconditional Default: ORN
    for net_name, net_cfg in DEFAULT_NETWORK_CONFIGS.items():
        dst = os.path.join(net_dir, f"{net_name}.json")
        if not os.path.exists(dst):
            try:
                with open(dst, "w", encoding="utf-8") as f:
                    json.dump(net_cfg, f, indent=2)
                    f.write("\n")
                created.append(dst)
            except Exception as e:
                print(f"   [Warning] Failed to write network config {dst}: {e}")

    # 2. Optional Networks (e.g. EU-Hydro): only create if underlying GIS files exist in gis/
    for net_name, net_cfg in OPTIONAL_NETWORK_CONFIGS.items():
        dst = os.path.join(net_dir, f"{net_name}.json")
        if _network_files_exist(net_cfg, base_dir, GIS_DIR):
            if not os.path.exists(dst):
                try:
                    with open(dst, "w", encoding="utf-8") as f:
                        json.dump(net_cfg, f, indent=2)
                        f.write("\n")
                    created.append(dst)
                except Exception as e:
                    print(f"   [Warning] Failed to write network config {dst}: {e}")
        else:
            if os.path.exists(dst):
                try:
                    os.remove(dst)
                except Exception:
                    pass

    if created:
        print(f"   [Config] Created network configuration(s) in {net_dir}/ ({len(created)} profile(s))")
    return created


def ensure_source_configs(base_dir=BASE_DIR):
    """
    Ensures that configs/sources directory exists and contains default source JSON configurations (RSN, SSN, RSS).
    Optional source configurations (e.g. waterblitz, river_guardians, rsn_gam) are only created if their
    corresponding directories exist in data/custom/.
    Returns list of paths created.
    """
    src_dir = _rel_join(base_dir, "configs", "sources")
    os.makedirs(src_dir, exist_ok=True)
    created = []

    # 1. Unconditional Defaults: RSN, SSN, RSS (composites)
    for src_name, src_cfg in DEFAULT_SOURCE_CONFIGS.items():
        dst = os.path.join(src_dir, f"{src_name}.json")
        if not os.path.exists(dst):
            try:
                with open(dst, "w", encoding="utf-8") as f:
                    json.dump(src_cfg, f, indent=2)
                    f.write("\n")
                created.append(dst)
            except Exception as e:
                print(f"   [Warning] Failed to write source config {dst}: {e}")

    # 2. Optional Sources: only create if underlying directories exist in data/
    for src_name, src_cfg in OPTIONAL_SOURCE_CONFIGS.items():
        dst = os.path.join(src_dir, f"{src_name}.json")
        if _source_files_exist(src_cfg, base_dir, DATA_DIR):
            if not os.path.exists(dst):
                try:
                    with open(dst, "w", encoding="utf-8") as f:
                        json.dump(src_cfg, f, indent=2)
                        f.write("\n")
                    created.append(dst)
                except Exception as e:
                    print(f"   [Warning] Failed to write source config {dst}: {e}")
        else:
            if os.path.exists(dst):
                try:
                    os.remove(dst)
                except Exception:
                    pass

    if created:
        print(f"   [Config] Created source configuration(s) in {src_dir}/ ({len(created)} profile(s))")
    return created


def ensure_all_configs(base_dir=BASE_DIR):
    """Ensures directory config, default network configs, and default source configs all exist."""
    ensure_directory_config(base_dir)
    ensure_network_configs(base_dir)
    ensure_source_configs(base_dir)


def load_directory_config(base_dir=BASE_DIR, auto_create=False):
    """
    Loads directory configuration from configs/directories.json (or paths.json / config.json).
    If auto_create=True, creates the file with default locations if not found.
    Resolves relative paths against base_dir.
    """
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

    data_dir = _resolve(raw_cfg.get("data_dir") or raw_cfg.get("data") or raw_cfg.get("DATA_DIR"), "data")
    gis_dir = _resolve(raw_cfg.get("gis_dir") or raw_cfg.get("gis") or raw_cfg.get("GIS_DIR"), "gis")
    output_dir = _resolve(
        raw_cfg.get("output_dir") or raw_cfg.get("out_dir") or raw_cfg.get("output") or raw_cfg.get("out") or raw_cfg.get("OUTPUT_DIR"),
        "output"
    )
    proc_val = raw_cfg.get("processed_dir") or raw_cfg.get("processed") or raw_cfg.get("PROCESSED_DIR")
    proc_dir = _resolve(proc_val, _rel_join(data_dir, "processed")) if proc_val else _rel_join(data_dir, "processed")

    return {
        "config_file": cfg_file if os.path.exists(cfg_file) else None,
        "data_dir": data_dir,
        "gis_dir": gis_dir,
        "output_dir": output_dir,
        "processed_dir": proc_dir,
    }


DIR_CONFIG = load_directory_config(BASE_DIR, auto_create=False)
GIS_DIR = os.environ.get("S3N_GIS_DIR") or DIR_CONFIG["gis_dir"]
DATA_DIR = os.environ.get("S3N_DATA_DIR") or DIR_CONFIG["data_dir"]
cand_rsn = _rel_join(DATA_DIR, "NCEA_River_Surveillance_Network_ARD_Open")
RSN_DIR = cand_rsn if os.path.isdir(cand_rsn) else (DATA_DIR if os.path.exists(os.path.join(DATA_DIR, "RSN Water Quality Samples.csv")) else cand_rsn)
cand_ssn = _rel_join(DATA_DIR, "NCEA_Small_Streams_Network_ARD_Open")
SSN_DIR = cand_ssn if os.path.isdir(cand_ssn) else (DATA_DIR if os.path.exists(os.path.join(DATA_DIR, "SSN Water Quality Samples.csv")) else cand_ssn)
WQ_DIR = RSN_DIR
CONFIG_DIR = _rel_join(BASE_DIR, "configs", "networks")
SOURCE_CONFIG_DIR = _rel_join(BASE_DIR, "configs", "sources")
OUT_DIR = os.environ.get("S3N_PROCESSED_DIR") or DIR_CONFIG["processed_dir"]


def wgs84_to_bng(lat, lon):
    """
    Converts WGS84 decimal latitude and longitude into British National Grid
    (OSGB36 / EPSG:27700) Eastings and Northings in meters using the standard
    Airy 1830 ellipsoid and Transverse Mercator projection (pure Python standard library).
    """
    a, b = 6377563.396, 6356256.909
    F0 = 0.9996012717
    lat0 = math.radians(49.0)
    lon0 = math.radians(-2.0)
    N0, E0 = -100000.0, 400000.0
    e2 = 1.0 - (b * b) / (a * a)
    n = (a - b) / (a + b)

    phi = math.radians(lat)
    lam = math.radians(lon)

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


def load_source_registry(base_dir=BASE_DIR):
    """
    Loads source configurations from `configs/sources/*.json` and auto-discovers
    subdirectories in `data/custom/`.
    """
    ensure_source_configs(base_dir)
    cfg_dir = os.path.join(base_dir, "configs", "sources")
    sources = {}
    composites = {}

    sources["rsn"] = {
        "name": "rsn",
        "abbr": "RSN",
        "label": "River Surveillance Network (RSN)",
        "description": "Environment Agency River Surveillance Network (1,382 stations)",
        "data_dir": "data/NCEA_River_Surveillance_Network_ARD_Open",
        "sites_file": "RSN Water Quality Sample Points.csv",
        "id_col": "network_id",
        "easting_col": "ngr_easting",
        "northing_col": "ngr_northing",
        "crs_epsg": 27700,
    }
    sources["ssn"] = {
        "name": "ssn",
        "abbr": "SSN",
        "label": "Small Streams Network (SSN)",
        "description": "Environment Agency Small Streams Network (313 headwater stations)",
        "data_dir": "data/NCEA_Small_Streams_Network_ARD_Open",
        "sites_file": "SSN Water Quality Sample Points.csv",
        "id_col": "network_id",
        "easting_col": "ngr_easting",
        "northing_col": "ngr_northing",
        "crs_epsg": 27700,
    }

    composites["rss"] = {
        "name": "rss",
        "abbr": "RSS",
        "label": "Combined River Surveillance & Small Streams (RSS)",
        "description": "Integrated River Surveillance Network and Small Streams Network",
        "sources": ["rsn", "ssn"],
    }
    composites["all"] = {
        "name": "all",
        "abbr": "ALL",
        "label": "All Registered Monitoring Networks",
        "description": "Union of all active monitoring networks",
        "sources": ["rsn", "ssn"],
    }

    if os.path.isdir(cfg_dir):
        for f in glob.glob(os.path.join(cfg_dir, "*.json")):
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

    # Auto-detect data/custom/* subdirectories
    cust_dir = os.path.join(base_dir, "data", "custom")
    if os.path.isdir(cust_dir):
        for sd in os.listdir(cust_dir):
            full_sd = os.path.join(cust_dir, sd)
            if os.path.isdir(full_sd):
                s_key = sd.lower()
                if s_key not in sources:
                    sources[s_key] = {
                        "name": s_key,
                        "abbr": s_key.upper(),
                        "label": f"Custom Dataset: {sd}",
                        "description": f"User custom dataset in data/custom/{sd}",
                        "data_dir": full_sd,
                        "sites_file": "sites.csv",
                        "id_col": "site_id",
                        "easting_col": "easting",
                        "northing_col": "northing",
                        "crs_epsg": 27700,
                    }

    return sources, composites


def resolve_source_tag(source="rsn", registry=None):
    """
    Resolves source tags, aliases, composite groups, and delimiter combinations
    (e.g., 'rsn+ssn', 'rsn,custom'). Returns (source_key, constituent_list, source_abbr, source_label).
    """
    if not source:
        return "rsn", ["rsn"], "RSN", "River Surveillance Network (RSN)"

    if registry is None:
        sources_reg, composites_reg = load_source_registry()
    else:
        sources_reg, composites_reg = registry

    raw_s = str(source).lower().strip()
    raw_s = re.sub(r"\btobitz\b", "tobit", raw_s)
    raw_s = re.sub(r"_tobitz\b", "_tobit", raw_s)
    tokens = [t.strip().lower() for t in re.split(r"[,+ ]+", raw_s) if t.strip()]
    if not tokens:
        tokens = ["rsn"]

    norm_tokens = []
    for tk in tokens:
        if tk.endswith("_tobitz"):
            tk = tk[:-1]
        elif tk == "tobitz":
            tk = "tobit"
        norm_tokens.append(tk)
    tokens = norm_tokens

    if len(tokens) > 1:
        if set(tokens) == {"rsn", "ssn"}:
            return "rss", ["rsn", "ssn"], "RSS", "Combined River Surveillance & Small Streams (RSS)"
        # Expand composite tokens and check aliases
        expanded = []
        for tk in tokens:
            if tk in composites_reg:
                expanded.extend(composites_reg[tk].get("sources", [tk]))
            elif tk in sources_reg and sources_reg[tk].get("alias_target"):
                expanded.append(sources_reg[tk]["alias_target"])
            elif f"{tk}_tobit" in sources_reg:
                expanded.append(f"{tk}_tobit")
            else:
                expanded.append(tk)
        unique_tokens = []
        for t in expanded:
            if t not in unique_tokens:
                unique_tokens.append(t)
        combo_key = "_".join(unique_tokens)
        combo_abbr = "+".join(sources_reg[t].get("abbr", t.upper()) if t in sources_reg else t.upper() for t in unique_tokens)
        combo_label = f"Combined Monitoring Sources ({combo_abbr})"
        return combo_key, unique_tokens, combo_abbr, combo_label

    s = tokens[0]
    if s in ("ssn", "small_stream", "small_streams", "small-streams", "small-stream"):
        return "ssn", ["ssn"], "SSN", "Small Streams Network (SSN)"
    if s in ("rss", "combined", "both", "all", "rsn_ssn", "ssn_rsn", "rsn+ssn"):
        return "rss", ["rsn", "ssn"], "RSS", "Combined River Surveillance & Small Streams (RSS)"
    if s in composites_reg:
        c_info = composites_reg[s]
        c_srcs = c_info.get("sources", ["rsn", "ssn"])
        c_abbr = c_info.get("abbr", s.upper())
        c_label = c_info.get("label", f"Composite: {s}")
        return s, c_srcs, c_abbr, c_label
    if s in sources_reg:
        s_info = sources_reg[s]
        if s_info.get("alias_target") and s_info["alias_target"] in sources_reg:
            t_key = s_info["alias_target"]
            t_info = sources_reg[t_key]
            return t_key, [t_key], t_info.get("abbr", t_key.upper()), t_info.get("label", f"Monitoring Source: {t_key}")
        s_abbr = s_info.get("abbr", s.upper())
        s_label = s_info.get("label", f"Monitoring Source: {s}")
        return s, [s], s_abbr, s_label

    # Check if a range-tagged profile exists (e.g. waterblitz -> waterblitz_tobit)
    for rm in ("tobit", "geometric", "midpoint", "lower", "upper"):
        cand_key = f"{s}_{rm}"
        if cand_key in sources_reg:
            t_info = sources_reg[cand_key]
            return cand_key, [cand_key], t_info.get("abbr", cand_key.upper()), t_info.get("label", f"Monitoring Source: {cand_key}")

    # Check if s has a range mode suffix
    for rm in ("tobit", "geometric", "midpoint", "lower", "upper"):
        if s.endswith(f"_{rm}"):
            base_s = s[:-len(rm)-1]
            if base_s in sources_reg:
                b_info = sources_reg[base_s]
                b_label = b_info.get("label", base_s.title())
                if f"[{rm.title()}]" not in b_label:
                    b_label = f"{b_label} [{rm.title()}]"
                b_abbr = f"{b_info.get('abbr', base_s.upper())}_{rm.upper()[:3]}"
                return s, [s], b_abbr, b_label

    # Check if decomposing an underscore-joined string matches multiple registered sources
    if "_" in s and s not in sources_reg and s not in ("rsn", "ssn", "rss"):
        reg_names = sorted(list(sources_reg.keys()) + ["rsn", "ssn"], key=len, reverse=True)
        matched = []
        rem = s
        while rem:
            found = False
            for sn in reg_names:
                if rem == sn or rem.startswith(f"{sn}_"):
                    matched.append(sn)
                    rem = rem[len(sn):].lstrip("_")
                    found = True
                    break
            if not found:
                break
        if matched and not rem and len(matched) > 1:
            unique_tokens = []
            for t in matched:
                if t in sources_reg and sources_reg[t].get("alias_target"):
                    t = sources_reg[t]["alias_target"]
                if t not in unique_tokens:
                    unique_tokens.append(t)
            combo_key = "_".join(unique_tokens)
            combo_abbr = "+".join(sources_reg[t].get("abbr", t.upper()) if t in sources_reg else t.upper() for t in unique_tokens)
            combo_label = f"Combined Monitoring Sources ({combo_abbr})"
            return combo_key, unique_tokens, combo_abbr, combo_label

    # Generic custom source
    clean_tag = re.sub(r"[^a-z0-9_]+", "_", s)
    return clean_tag, [clean_tag], clean_tag.upper(), f"Custom Monitoring Source ({clean_tag.upper()})"


def detect_range_mode(source_key, sources_reg, custom_samples=None, custom_dir=None, cli_mode="auto"):
    """
    Detects whether range conversion was used and returns the range mode tag
    ('tobit', 'geometric', 'midpoint', 'lower', 'upper', or None).
    """
    if cli_mode and cli_mode not in ("auto", "none"):
        return cli_mode
    if cli_mode == "none":
        return None

    for m in ("tobit", "geometric", "midpoint", "lower", "upper"):
        if source_key.endswith(f"_{m}"):
            return m

    if source_key in sources_reg:
        rm = sources_reg[source_key].get("range_mode")
        if rm:
            return rm

    cand_samples = custom_samples
    if not cand_samples and custom_dir:
        cand_samples = os.path.join(custom_dir, "samples.csv")
    if not cand_samples and source_key in sources_reg:
        s_dir = sources_reg[source_key].get("data_dir", "")
        if s_dir:
            cand = os.path.join(s_dir, "samples.csv")
            if os.path.exists(cand):
                cand_samples = cand
            else:
                cand = os.path.join(BASE_DIR, s_dir, "samples.csv")
                if os.path.exists(cand):
                    cand_samples = cand

    if cand_samples and os.path.exists(cand_samples):
        try:
            with open(cand_samples, "r", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                if reader.fieldnames and "meas_method" in reader.fieldnames:
                    for i, r in enumerate(reader):
                        if i > 500:
                            break
                        mm = (r.get("meas_method") or "").strip().lower()
                        if "tobit" in mm:
                            return "tobit"
                        if "geometric" in mm:
                            return "geometric"
                        if "midpoint" in mm or "arithmetic" in mm:
                            return "midpoint"
                        if "lower" in mm:
                            return "lower"
                        if "upper" in mm:
                            return "upper"
        except Exception:
            pass

    return None


def print_available_sources(base_dir=BASE_DIR):
    sources_reg, composites_reg = load_source_registry(base_dir=base_dir)
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
    """Updates global GIS_DIR, DATA_DIR, RSN_DIR, SSN_DIR, WQ_DIR, and OUT_DIR from CLI or environment overrides."""
    global GIS_DIR, DATA_DIR, RSN_DIR, SSN_DIR, WQ_DIR, OUT_DIR, DIR_CONFIG
    DIR_CONFIG = load_directory_config(BASE_DIR, auto_create=False)
    if gis_dir:
        GIS_DIR = gis_dir
    else:
        GIS_DIR = os.environ.get("S3N_GIS_DIR") or DIR_CONFIG["gis_dir"]

    if data_dir:
        DATA_DIR = data_dir
    else:
        DATA_DIR = os.environ.get("S3N_DATA_DIR") or DIR_CONFIG["data_dir"]

    cand_rsn = _rel_join(DATA_DIR, "NCEA_River_Surveillance_Network_ARD_Open")
    RSN_DIR = cand_rsn if os.path.isdir(cand_rsn) else (DATA_DIR if os.path.exists(os.path.join(DATA_DIR, "RSN Water Quality Samples.csv")) else cand_rsn)
    cand_ssn = _rel_join(DATA_DIR, "NCEA_Small_Streams_Network_ARD_Open")
    SSN_DIR = cand_ssn if os.path.isdir(cand_ssn) else (DATA_DIR if os.path.exists(os.path.join(DATA_DIR, "SSN Water Quality Samples.csv")) else cand_ssn)
    WQ_DIR = RSN_DIR

    if out_dir:
        OUT_DIR = out_dir
    else:
        OUT_DIR = os.environ.get("S3N_PROCESSED_DIR") or DIR_CONFIG["processed_dir"]
    os.makedirs(OUT_DIR, exist_ok=True)


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
    sources_reg, composites_reg = load_source_registry(base_dir=BASE_DIR)
    source_key, constituents, source_abbr, source_label = resolve_source_tag(source, (sources_reg, composites_reg))

    cand_rsn = _rel_join(data_dir, "NCEA_River_Surveillance_Network_ARD_Open")
    rsn_dir = cand_rsn if os.path.isdir(cand_rsn) else (data_dir if os.path.exists(os.path.join(data_dir, "RSN Water Quality Samples.csv")) else cand_rsn)
    cand_ssn = _rel_join(data_dir, "NCEA_Small_Streams_Network_ARD_Open")
    ssn_dir = cand_ssn if os.path.isdir(cand_ssn) else (data_dir if os.path.exists(os.path.join(data_dir, "SSN Water Quality Samples.csv")) else cand_ssn)

    pts = {}
    pts_source = {}
    core_wq_sites = set()
    all_records = []

    for src_item in constituents:
        src_label = src_item.upper()
        # 1. Custom CLI override
        if src_item in ("custom", custom_name.lower()) or custom_sites or custom_dir:
            sites_p = custom_sites or (os.path.join(custom_dir, "sites.csv") if custom_dir else None)
            samples_p = custom_samples or (os.path.join(custom_dir, "samples.csv") if custom_dir else None)
            if not sites_p or not os.path.exists(sites_p):
                # Check data/custom/<src_item>/sites.csv
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

                    # Auto-convert lat/lon in decimal degrees to British National Grid (EPSG:27700)
                    if (abs(px) <= 180.0 and abs(py) <= 90.0) or (custom_crs == 4326):
                        px, py = wgs84_to_bng(py, px)

                    r_compat = {"wims_id": raw_id, "wims_region_code": "CUSTOM"}
                    pts[raw_id] = (px, py, r_compat)
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

                    for r in s_reader:
                        s_nid = r.get(s_id_col, "").strip()
                        if not s_nid:
                            continue
                        if det_col_cand and det_code:
                            row_det = r.get(det_col_cand, "").strip()
                            if row_det and row_det != det_code and (det_code.lower() not in row_det.lower()):
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
                        dt_str = r.get(dt_col, "").strip() if dt_col else "2023-06-01"
                        yr = int(dt_str[0:4]) if len(dt_str) >= 4 and dt_str[0:4].isdigit() else 2023
                        mo = int(dt_str[5:7]) if len(dt_str) >= 7 and dt_str[5:7].isdigit() else 6
                        seas = MONTH_TO_SEASON.get(mo, "summer")
                        unc_val = None
                        if unc_col:
                            raw_unc = r.get(unc_col, "").strip()
                            if raw_unc:
                                try:
                                    unc_val = float(raw_unc)
                                except ValueError:
                                    unc_val = None
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

                    r_compat = {"wims_id": raw_id, "wims_region_code": cfg_s.get("abbr", "SRC")}
                    pts[raw_id] = (px, py, r_compat)
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

                    for r in s_reader:
                        s_nid = r.get(s_id_col, "").strip()
                        if not s_nid:
                            continue
                        if det_col_c and det_code:
                            row_det = r.get(det_col_c, "").strip()
                            if row_det and row_det != det_code and (det_code.lower() not in row_det.lower()):
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
                        dt_str = r.get(dt_col, "").strip() if dt_col else "2023-06-01"
                        yr = int(dt_str[0:4]) if len(dt_str) >= 4 and dt_str[0:4].isdigit() else 2023
                        mo = int(dt_str[5:7]) if len(dt_str) >= 7 and dt_str[5:7].isdigit() else 6
                        seas = MONTH_TO_SEASON.get(mo, "summer")
                        unc_val = None
                        if unc_col:
                            raw_unc = r.get(unc_col, "").strip()
                            if raw_unc:
                                try:
                                    unc_val = float(raw_unc)
                                except ValueError:
                                    unc_val = None
                        meth_val = r.get(meth_col, "").strip().lower() if meth_col else str(cfg_s.get("meas_method", "")).lower()
                        all_records.append((s_nid, val, is_cens, yr, mo, seas, unc_val, meth_val))
            continue

        # 3. Standard RSN / SSN official datasets or auto-discovered custom folder
        if src_item not in ("rsn", "ssn"):
            # Check typo normalization (e.g. tobitz -> tobit)
            norm_item = src_item.replace("tobitz", "tobit")
            if norm_item in sources_reg:
                src_item = norm_item
            cand_cs = os.path.join(data_dir, "custom", src_item)
            if not os.path.isdir(cand_cs):
                cand_cs = os.path.join(BASE_DIR, "data", "custom", src_item)
            if not os.path.isdir(cand_cs) and src_item in sources_reg:
                s_d = sources_reg[src_item].get("data_dir", "")
                if s_d and os.path.isdir(s_d):
                    cand_cs = s_d
            if not os.path.isdir(cand_cs):
                valid_srcs = sorted(list(sources_reg.keys()) + ["rsn", "ssn", "rss"])
                hint = f" Did you mean '{norm_item}'?" if norm_item in sources_reg else ""
                raise ValueError(f"Unknown monitoring data source '{src_item}'.{hint} Available sources: {valid_srcs}")

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
                        r_compat = {"wims_id": nid, "wims_region_code": ""}
                        pts[nid] = (float(r["ngr_easting"]), float(r["ngr_northing"]), r_compat)
                        pts_source[nid] = prefix

        wq_samples_path = os.path.join(s_dir, f"{prefix} Water Quality Samples.csv")
        if os.path.exists(wq_samples_path):
            with open(wq_samples_path, "r", encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    if r.get("det_code") == "0117":
                        core_wq_sites.add(r["network_id"])

        if dataset == "wq":
            if os.path.exists(wq_samples_path):
                with open(wq_samples_path, "r", encoding="utf-8-sig") as f:
                    for r in csv.DictReader(f):
                        if r.get("det_code") == det_code:
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
        elif dataset == "lcms":
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
            z = [(math.pow(v, lam) - 1.0) / lam for v in pos]
        mz = statistics.mean(z)
        vz = statistics.mean([(x - mz) ** 2 for x in z])
        if vz < 1e-12:
            return -1e12
        return -0.5 * n * math.log(vz) + (lam - 1.0) * sum_log

    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lower, upper
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    while abs(b - a) > 1e-4:
        if ll(c) > ll(d):
            b = d
        else:
            a = c
        c = b - gr * (b - a)
        d = a + gr * (b - a)
    return round((b + a) / 2.0, 4)


def estimate_yeojohnson_lambda(vals, lower=-2.0, upper=2.0):
    clean = [v for v in vals if v is not None and not math.isnan(v)]
    n = len(clean)
    if n < 10:
        return 0.5

    def ll(lam):
        z = []
        jac = 0.0
        for y in clean:
            if y >= 0:
                if abs(lam) < 1e-4:
                    zv = math.log(y + 1.0)
                else:
                    zv = (math.pow(y + 1.0, lam) - 1.0) / lam
                jac += (lam - 1.0) * math.log(y + 1.0)
            else:
                if abs(2.0 - lam) < 1e-4:
                    zv = -math.log(-y + 1.0)
                else:
                    zv = -((math.pow(-y + 1.0, 2.0 - lam) - 1.0) / (2.0 - lam))
                jac += (1.0 - lam) * -math.log(-y + 1.0)
            z.append(zv)
        mz = statistics.mean(z)
        vz = statistics.mean([(x - mz) ** 2 for x in z])
        if vz < 1e-12:
            return -1e12
        return -0.5 * n * math.log(vz) + jac

    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lower, upper
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    while abs(b - a) > 1e-4:
        if ll(c) > ll(d):
            b = d
        else:
            a = c
        c = b - gr * (b - a)
        d = a + gr * (b - a)
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
    elif t == "log10":
        return math.log10(max(1e-9, val)) if val > 0 else None
    elif t == "log1p":
        return math.log1p(max(0.0, val)) if val >= 0 else None
    elif t in ("sqrt", "square_root"):
        return math.sqrt(max(0.0, val)) if val >= 0 else None
    elif t in ("cbrt", "cube_root"):
        return math.copysign(math.pow(abs(val), 1.0 / 3.0), val)
    elif t in ("fourthrt", "qtrt"):
        return math.pow(max(0.0, val), 0.25) if val >= 0 else None
    elif t in ("reciprocal", "inv"):
        return -1.0 / max(1e-6, val) if val > 0 else None
    elif t in ("boxcox", "bc"):
        if val <= 0:
            return None
        l_val = lam if lam is not None else 0.5
        if abs(l_val) < 1e-4:
            return math.log(val)
        return (math.pow(val, l_val) - 1.0) / l_val
    elif t in ("yeojohnson", "yj"):
        l_val = lam if lam is not None else 0.5
        if val >= 0:
            return math.log(val + 1.0) if abs(l_val) < 1e-4 else (math.pow(val + 1.0, l_val) - 1.0) / l_val
        else:
            return -math.log(-val + 1.0) if abs(2.0 - l_val) < 1e-4 else -((math.pow(-val + 1.0, 2.0 - l_val) - 1.0) / (2.0 - l_val))
    elif t in ("asin", "arcsin", "arcsine", "asin_sqrt"):
        p = min(1.0, max(0.0, val / scale))
        return math.asin(math.sqrt(p))
    elif t == "logit":
        p = min(1.0 - 1e-4, max(1e-4, val / scale))
        return math.log(p / (1.0 - p))
    else:
        return float(val)


def print_available_transforms():
    rows = [
        ("log", "ln", "ln(y)", "y > 0", "Default for water chemistry & LCMS. Multiplicative concentrations spanning orders of magnitude."),
        ("sqrt", "square_root", "y^(1/2)", "y >= 0", "Strongly recommended for moderately skewed concentrations (e.g. UK Nitrate, skewness drop +1.00 -> +0.09) and taxa counts."),
        ("cbrt", "cube_root", "sign(y)*|y|^(1/3)", "All real y", "Wilson-Hilferty transformation for Gamma-distributed hydrological variables. Compresses high values without crushing near-zero observations."),
        ("log1p", "log_1p", "ln(1 + y)", "y >= 0", "Accommodates true zeros (y=0) without arbitrary LOD substitution. Gentle log compression for small values."),
        ("log10", "log_10", "log10(y)", "y > 0", "Order-of-magnitude chemical concentrations (e.g. pathogen counts, electrical conductivity)."),
        ("fourthrt", "qtrt", "y^(1/4)", "y >= 0", "Intermediate between square root and log; classic ecological standard for heavy-tailed abundance and biomass."),
        ("boxcox", "bc, box_cox", "(y^lambda-1)/lam", "y > 0", "Analytically optimizes residual normality via Profile Maximum Likelihood lambda in [-2, 2] (or custom --lambda)."),
        ("yeojohnson", "yj", "Yeo-Johnson(lam)", "All real y", "Extension of Box-Cox supporting zero and negative values. Automatic Profile MLE estimation of lambda."),
        ("reciprocal", "inv, inverse", "-1 / y", "y > 0", "Monotonic rate-dilution processes. Preserves positive spatial gradient direction (higher y -> higher z)."),
        ("asin", "arcsin, arcsine", "asin(sqrt(y/sc))", "0 <= y <= scale", "Angular transformation for bounded proportions and percentages (e.g. silt cover, diatom valences)."),
        ("logit", "logit", "ln(p / (1 - p))", "0 < y < scale", "S-shaped logistic transformation for bounded proportions and sediment percentages."),
        ("none", "raw, identity", "y", "All real y", "Default for ecological indices (WHPT ASPT, LIFE index, RHS scores) that are already normally distributed.")
    ]
    print("\nAvailable Data Transformations (--transform=<type> [--transform-lambda=<val>]):\n")
    print(f"{'Type':<12} | {'Short Aliases':<16} | {'Formula z = f(y)':<18} | {'Valid Domain':<16} | {'Best Applied When'}")
    print("-" * 125)
    for t_type, alias, formula, domain, note in rows:
        print(f"{t_type:<12} | {alias:<16} | {formula:<18} | {domain:<16} | {note}")
    print()


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
    Universal parser for GeoPackage Binary (GPB) and standard WKB / ISO WKB / EWKB
    (as used in GeoParquet) LineString & MultiLineString geometries.
    Supports 2D (XY), 3D (XYZ / XYM), and 4D (XYZM) in both ISO WKB (1000/2000/3000)
    and OGC EWKB (0x80000000 / 0x40000000 / 0x20000000), returning `(pts_xy, zs)` where:
      - `pts_xy` : list of (x, y) tuples
      - `zs`     : list of valid vertex Z elevations (if 3D/4D Z is present, else [])
    """
    if not blob:
        return [], []
    if isinstance(blob, str):
        s = blob.strip()
        if not s:
            return [], []
        try:
            blob = bytes.fromhex(s[2:] if s.lower().startswith("0x") else s)
        except ValueError:
            return [], []
    elif isinstance(blob, memoryview):
        blob = blob.tobytes()

    if len(blob) < 9:
        return [], []

    # Detect OGC GeoPackage Binary header ("GP" = 0x47, 0x50) vs raw WKB / ISO WKB / EWKB (0x00 or 0x01)
    if blob[0] == 0x47 and blob[1] == 0x50:
        env_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
        flags = blob[3]
        env_ind = (flags >> 1) & 0x07
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

    pts = []
    zs = []
    if base_type == 5:  # MultiLineString
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
            for i in range(0, g_dim * npts, g_dim):
                pts.append((coords[i], coords[i + 1]))
                if g_z:
                    zv = coords[i + 2]
                    if not math.isnan(zv) and -500.0 < zv < 9000.0:
                        zs.append(zv)
    elif base_type == 2:  # LineString
        npts = struct.unpack_from(bo + "I", blob, offset + hdr_len)[0]
        offset += hdr_len + 4
        coords = struct.unpack_from(bo + f"{ndim * npts}d", blob, offset)
        for i in range(0, ndim * npts, ndim):
            pts.append((coords[i], coords[i + 1]))
            if has_z:
                zv = coords[i + 2]
                if not math.isnan(zv) and -500.0 < zv < 9000.0:
                    zs.append(zv)
    return pts, zs


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


def parse_gpkg_multilinestring(blob):
    """Backward-compatible helper returning only 2D (x, y) vertices."""
    pts, _ = parse_gpkg_geometry(blob)
    return pts


def transform_coords_batch(xy_list, from_epsg, to_epsg):
    """
    Transforms a list of (x, y) coordinates from `from_epsg` to `to_epsg` using
    `cs2cs` with traditional GIS easting/northing (`+init=epsg:...`) order.
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
        target_path = os.path.join(CONFIG_DIR, f"{mapped}.json")
        if os.path.exists(target_path):
            with open(target_path, "r", encoding="utf-8") as f:
                return json.load(f)
        if mapped in DEFAULT_NETWORK_CONFIGS:
            return DEFAULT_NETWORK_CONFIGS[mapped]
        if mapped in OPTIONAL_NETWORK_CONFIGS and _network_files_exist(OPTIONAL_NETWORK_CONFIGS[mapped], BASE_DIR, GIS_DIR):
            return OPTIONAL_NETWORK_CONFIGS[mapped]

    if net_clean.lower() in DEFAULT_NETWORK_CONFIGS:
        return DEFAULT_NETWORK_CONFIGS[net_clean.lower()]

    if net_clean.lower() in OPTIONAL_NETWORK_CONFIGS and _network_files_exist(OPTIONAL_NETWORK_CONFIGS[net_clean.lower()], BASE_DIR, GIS_DIR):
        return OPTIONAL_NETWORK_CONFIGS[net_clean.lower()]

    if net_clean.endswith(".json") and os.path.exists(net_clean):
        with open(net_clean, "r", encoding="utf-8") as f:
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


def write_outputs_to_dirs(net_abbr, files_dict, source="rsn"):
    """
    Writes generated output files into `data/processed/<net_abbr>_<source>/` (e.g. `orn_rsn`, `orn_ssn`)
    AND copies them into `data/processed/<net_abbr>/` and `data/processed/` root when source is 'rsn'
    so downstream scripts work seamlessly with or without `--network`.
    """
    source_key, _, _, _ = resolve_source_tag(source)
    tagged_name = f"{net_abbr}_{source_key}"
    net_dir = os.path.join(DATA_DIR, "processed", tagged_name)
    os.makedirs(net_dir, exist_ok=True)
    for fname, src_path in files_dict.items():
        dst_net = os.path.join(net_dir, fname)
        if os.path.abspath(src_path) != os.path.abspath(dst_net):
            shutil.copy2(src_path, dst_net)

    # For default RSN, also ensure outputs are available at bare net_abbr (orn) and proc_root fallback
    if source_key == "rsn":
        bare_dir = os.path.join(DATA_DIR, "processed", net_abbr)
        os.makedirs(bare_dir, exist_ok=True)
        proc_root = os.path.join(DATA_DIR, "processed")
        for fname, src_path in files_dict.items():
            dst_bare = os.path.join(bare_dir, fname)
            if os.path.abspath(src_path) != os.path.abspath(dst_bare):
                shutil.copy2(src_path, dst_bare)
            dst_root = os.path.join(proc_root, fname)
            if os.path.abspath(src_path) != os.path.abspath(dst_root):
                shutil.copy2(src_path, dst_root)


def load_orn_spatial_reference():
    """
    Loads the pre-extracted ORN England reaches (`data/processed/orn_rsn/england_rsn_reaches.csv`)
    into a fast 2 km spatial hash grid in EPSG:27700 so non-ORN networks in England
    can inherit BGS bedrock/superficial geology (`Is_Chalk_Limestone`, `Is_Sandstone`,
    `Is_Aquifer`, `Is_Peat`), `RBD_ID`, and fallback `Elevation`/`Slope` in < 0.2s.
    """
    orn_csv = os.path.join(DATA_DIR, "processed", "orn_rsn", "england_rsn_reaches.csv")
    if not os.path.exists(orn_csv):
        orn_csv = os.path.join(DATA_DIR, "processed", "orn", "england_rsn_reaches.csv")
    if not os.path.exists(orn_csv):
        orn_csv = os.path.join(OUT_DIR, "england_rsn_reaches.csv")
    if not os.path.exists(orn_csv):
        orn_csv = os.path.join(DATA_DIR, "processed", "england_rsn_reaches.csv")
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


def preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, nitrate_by_site, period_label, t0,
                               trans_type="log", eff_lambda=None, pts_source=None, source_key="rsn", range_mode=None):
    """
    Network-agnostic preprocessing pipeline for any configured GeoPackage network
    (such as EU-Hydro Thames, EU-Hydro Tweed, combined EU-Hydro, OS Water Network, etc.).
    """
    if pts_source is None:
        pts_source = {}
    source_key, constituents, source_abbr, source_label = resolve_source_tag(source_key)
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
        poly, zs = parse_gpkg_geometry(blob)
        if len(poly) < 2:
            return
        raw_id_val = rdict.get(id_col) if id_col else None
        if raw_id_val is None:
            raw_id_val = rdict.get("FID", len(raw_reaches) + 1)
        orig_id = str(raw_id_val)
        raw_cid = rdict.get(catch_col) if catch_col else None
        cid_key = (file_idx, _parse_cid_val(raw_cid))

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

    # Assign unique positive integer RivIDs
    key_to_rivid = {}
    global_orig_to_rivid = {}
    reaches_by_id = {}
    for idx, (file_idx, orig_id, cid_key, poly, zs, rdict) in enumerate(filtered_reaches, 1):
        if len(gpkg_paths) == 1 and len(layers) == 1 and str(orig_id).isdigit():
            riv_id = int(orig_id)
        else:
            riv_id = idx
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
            fn_key = str(rdict.get(fnode_col))
            tn_key = str(rdict.get(tnode_col))
            fn_el = (node_lookup.get((f_idx, fn_key)) or node_lookup.get(("*", fn_key)) or (None, None))[0]
            tn_el = (node_lookup.get((f_idx, tn_key)) or node_lookup.get(("*", tn_key)) or (None, None))[0]
            vals_el = [float(x) for x in (fn_el, tn_el) if x is not None]
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

    site_by_reach = defaultdict(list)
    sites_csv_path = os.path.join(OUT_DIR, "rsn_sites.csv")
    with open(sites_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "site_id", "wims_id", "wims_region_code",
            "easting", "northing", "snap_dist_m", "RivID", "CatchID",
            "obs_up_dist_km", "n_samples", "nitrate_mean", "nitrate_median",
            "nitrate_sd", "log_nitrate_mean", "mean_log_nitrate",
            "meas_uncertainty", "meas_method", "included_250m",
            "source"
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
            src_str = pts_source.get(nid, "RSN")
            vals = nitrate_by_site.get(nid, [])
            u_vals = unc_by_site.get(nid, [])
            m_set = meth_by_site.get(nid, set())
            site_u = round(math.sqrt(sum(u * u for u in u_vals) / len(u_vals)), 6) if u_vals else ""
            site_m = sorted(m_set)[0] if m_set else ""
            if vals:
                m_val = statistics.mean(vals)
                med_val = statistics.median(vals)
                sd_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                trans_m_val = apply_forward_transform(m_val, ttype=trans_type, lam=eff_lambda)
                log_m = round(trans_m_val, 5) if trans_m_val is not None else ""
                trans_clean = [
                    tv for tv in (apply_forward_transform(v, ttype=trans_type, lam=eff_lambda) for v in vals)
                    if tv is not None
                ]
                m_log = round(statistics.mean(trans_clean), 5) if trans_clean else ""
                if inc:
                    site_by_reach[riv_id].append((nid, d, along_m, vals, u_vals, m_set))
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), len(vals), round(m_val, 4), round(med_val, 4),
                    round(sd_val, 4), log_m, m_log, site_u, site_m, inc, src_str
                ])
            else:
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), 0, "", "", "", "", "", "", "", inc, src_str
                ])

    reach_period_vals = defaultdict(list)
    reach_period_unc = defaultdict(list)
    reach_period_meth = defaultdict(set)
    reach_along_map = defaultdict(list)
    for rec in imputed_records:
        nid, imp_val, yr, mo, seas = rec[0], rec[1], rec[2], rec[3], rec[4]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
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
            if unc_val is not None:
                reach_period_unc[(riv_id, "all", "all")].append(unc_val)
                reach_period_unc[(riv_id, "year", str(yr))].append(unc_val)
                reach_period_unc[(riv_id, "season", seas)].append(unc_val)
                reach_period_unc[(riv_id, "month", f"{mo:02d}")].append(unc_val)
                reach_period_unc[(riv_id, "year_season", f"{yr}_{seas}")].append(unc_val)
            if meth_val:
                reach_period_meth[(riv_id, "all", "all")].add(meth_val)
                reach_period_meth[(riv_id, "year", str(yr))].add(meth_val)
                reach_period_meth[(riv_id, "season", seas)].add(meth_val)
                reach_period_meth[(riv_id, "month", f"{mo:02d}")].add(meth_val)
                reach_period_meth[(riv_id, "year_season", f"{yr}_{seas}")].add(meth_val)

    periods_csv_path = os.path.join(OUT_DIR, "rsn_reach_periods.csv")
    with open(periods_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "RivID", "CatchID", "period_type", "period_value",
            "obs_up_dist_km", "n_samples", "nitrate_mean", "nitrate_median", "log_nitrate_mean",
            "meas_uncertainty", "meas_method"
        ])
        for (riv_id, p_type, p_val), vals in sorted(reach_period_vals.items()):
            u_info = reaches_by_id[riv_id]
            d2m_km = u_info["d2m_km"]
            along_m = statistics.mean(reach_along_map[riv_id])
            obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
            m_v = statistics.mean(vals)
            med_v = statistics.median(vals)
            trans_mv_val = apply_forward_transform(m_v, ttype=trans_type, lam=eff_lambda)
            log_mv = round(trans_mv_val, 5) if trans_mv_val is not None else ""
            u_list = reach_period_unc.get((riv_id, p_type, p_val), [])
            m_set = reach_period_meth.get((riv_id, p_type, p_val), set())
            p_unc = round(math.sqrt(sum(u * u for u in u_list) / len(u_list)), 6) if u_list else ""
            p_meth = sorted(m_set)[0] if m_set else ""
            w.writerow([
                riv_id, u_info["catch_id"], p_type, p_val,
                round(obs_up_km, 5), len(vals), round(m_v, 4), round(med_v, 4), log_mv,
                p_unc, p_meth
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
        "nitrate_mean", "nitrate_median", "log_nitrate_mean",
        "meas_uncertainty", "meas_method", "n_samples"
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
                all_uncs = [u for s in s_list for u in s[4]]
                all_meths = {m for s in s_list for m in s[5] if m}
                along_m = statistics.mean(s[2] for s in s_list)
                snap_d = min(s[1] for s in s_list)
                obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
                m_no3 = statistics.mean(all_vals)
                med_no3 = statistics.median(all_vals)
                trans_no3_val = apply_forward_transform(m_no3, ttype=trans_type, lam=eff_lambda)
                log_no3 = round(trans_no3_val, 5) if trans_no3_val is not None else ""
                r_unc = round(math.sqrt(sum(u * u for u in all_uncs) / len(all_uncs)), 6) if all_uncs else ""
                r_meth = sorted(all_meths)[0] if all_meths else ""
                row_out = [
                    riv_id, u_info["catch_id"], u_info["rbd_id"], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(u_info["us_accum_km"], 5),
                    round(log_afv[riv_id], 6), u_info["shreve"], u_info["strahler"],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(u_info["elevation"], 3), round(u_info["slope"], 4), round(u_info["sinuosity"], 4), round(u_info["src_elevation"], 3),
                    u_info["is_chalk_lime"], u_info["is_sand"], u_info["is_aq"], u_info["is_peat"],
                    1, s_list[0][0], round(snap_d, 2), round(obs_up_km, 5),
                    round(m_no3, 4), round(med_no3, 4), log_no3,
                    r_unc, r_meth, len(all_vals)
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
                    "", "", "", "", "", 0
                ]
            w.writerow(row_out + extra_vals)

    active_meas_method = "gam" if any("gam" in s[5] for s_list in site_by_reach.values() for s in s_list) else None
    meta_path = os.path.join(OUT_DIR, "network_metadata.json")
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
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_obj, f, indent=2)

    write_outputs_to_dirs(net_abbr, {
        "england_rsn_reaches.csv": out_reaches_path,
        "rsn_sites.csv": sites_csv_path,
        "rsn_reach_periods.csv": periods_csv_path,
        "rsn_catchments_lookup.csv": catch_lu_path,
        "network_metadata.json": meta_path,
    }, source=source_key)
    print(f"Done in {time.time() - t0:.2f} seconds (saved to {OUT_DIR}/).")


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess stream network GeoPackages (ORN v2, EU-Hydro, etc.) and EA RSN samples for S3N modeling."
    )
    parser.add_argument(
        "--gis-dir",
        "--gis",
        dest="gis_dir",
        default=None,
        help="Custom GIS directory containing water network GeoPackages (default: 'gis').",
    )
    parser.add_argument(
        "--data-dir",
        "--data",
        dest="data_dir",
        default=None,
        help="Custom data directory containing RSN sample folders and processed outputs (default: 'data').",
    )
    parser.add_argument(
        "--out-dir",
        "--output-dir",
        "--out",
        "--processed-dir",
        dest="out_dir",
        default=None,
        help="Custom output directory for processed network CSVs (default: '<data-dir>/processed').",
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
        "--source",
        "--data-source",
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
    args = parser.parse_args()
    ensure_all_configs(BASE_DIR)
    if args.list_transforms:
        print_available_transforms()
        return
    if args.list_sources:
        print_available_sources()
        return

    # Auto-route source to custom-name if custom paths supplied
    if (args.custom_sites or args.custom_dir) and args.source == "rsn":
        args.source = args.custom_name

    sources_reg, composites_reg = load_source_registry(base_dir=BASE_DIR)
    source_key, constituents, source_abbr, source_label = resolve_source_tag(args.source, (sources_reg, composites_reg))

    active_range_mode = detect_range_mode(
        source_key=source_key,
        sources_reg=sources_reg,
        custom_samples=args.custom_samples,
        custom_dir=args.custom_dir,
        cli_mode=args.range_mode
    )

    if active_range_mode:
        if source_key != "rsn" and not source_key.endswith(f"_{active_range_mode}"):
            source_key = f"{source_key}_{active_range_mode}"
            if not source_abbr.endswith(active_range_mode.upper()[:3]):
                source_abbr = f"{source_abbr}_{active_range_mode.upper()[:3]}"
            if f"[{active_range_mode.title()}]" not in source_label:
                source_label = f"{source_label} [{active_range_mode.title()}]"
        print(f"   [Range Handling] Detected range mode: {active_range_mode.upper()} (Processed folder tagged with _{active_range_mode})")

    configure_runtime_dirs(gis_dir=args.gis_dir, data_dir=args.data_dir, out_dir=args.out_dir)

    if args.list_networks:
        list_available_networks()
        return

    cfg = load_network_config(args.network, args.network_config)
    net_abbr = cfg["abbr"]

    if not args.out_dir:
        sub_name = f"{net_abbr}_{source_key}"
        OUT_DIR = _rel_join(DATA_DIR, "processed", sub_name)
    os.makedirs(OUT_DIR, exist_ok=True)

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
        source=args.source,
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
        raw_all_by_site[nid].append((val, is_cens))
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

    nitrate_by_site = defaultdict(list)
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
        nitrate_by_site[nid].append(imp_val)
        if unc_val is not None:
            unc_by_site[nid].append(unc_val)
        if meth_val:
            meth_by_site[nid].add(meth_val)
        n_slice_samples += 1
        if is_cens:
            n_slice_cens += 1

    print(
        f"   Found {len(nitrate_by_site)} active stations in period '{period_label}' "
        f"({n_slice_samples} samples; {n_slice_cens} censored '<' "
        f"[{100.0 * n_slice_cens / max(1, n_slice_samples):.2f}%] imputed via '{args.censor_method}')."
    )

    # Resolve transformation type and Profile Maximum Likelihood lambda
    trans_type = (args.transform or "log").lower().strip()
    all_sample_vals = [v for s_vals in nitrate_by_site.values() for v in s_vals if v is not None]
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
        preprocess_generic_network(cfg, args, pts, core_wq_sites, imputed_records, nitrate_by_site, period_label, t0,
                                   trans_type=trans_type, eff_lambda=eff_lambda,
                                   pts_source=pts_source, source_key=source_key,
                                   range_mode=active_range_mode)
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

    active_monitored_sites = set(nitrate_by_site.keys()) | core_wq_sites
    monitored_catchments = {
        best_snap[nid][2]
        for nid in active_monitored_sites
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
            "nitrate_sd", "log_nitrate_mean", "mean_log_nitrate",
            "meas_uncertainty", "meas_method", "included_250m",
            "source"
        ])
        for nid in all_station_ids:
            px, py, meta = pts[nid]
            d, riv_id, catch_id, along, d2m = best_snap[nid]
            if catch_id not in monitored_catchments and nid not in core_wq_sites:
                continue
            obs_up_km = max(0.0, (d2m - along) / 1000.0)
            inc = 1 if (d <= max_snap_m and catch_id in monitored_catchments) else 0
            src_str = pts_source.get(nid, "RSN")
            vals = nitrate_by_site.get(nid, [])
            u_vals = unc_by_site.get(nid, [])
            m_set = meth_by_site.get(nid, set())
            site_u = round(math.sqrt(sum(u * u for u in u_vals) / len(u_vals)), 6) if u_vals else ""
            site_m = sorted(m_set)[0] if m_set else ""
            if vals:
                m_val = statistics.mean(vals)
                med_val = statistics.median(vals)
                sd_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                trans_m_val = apply_forward_transform(m_val, ttype=trans_type, lam=eff_lambda)
                log_m = round(trans_m_val, 5) if trans_m_val is not None else ""
                trans_clean = [
                    tv for tv in (apply_forward_transform(v, ttype=trans_type, lam=eff_lambda) for v in vals)
                    if tv is not None
                ]
                m_log = round(statistics.mean(trans_clean), 5) if trans_clean else ""
                if inc:
                    site_by_reach[riv_id].append((nid, d, along, vals, u_vals, m_set))
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), len(vals), round(m_val, 4), round(med_val, 4),
                    round(sd_val, 4), log_m, m_log, site_u, site_m, inc, src_str
                ])
            else:
                w.writerow([
                    nid, meta["wims_id"], meta["wims_region_code"],
                    px, py, round(d, 3), riv_id, catch_id,
                    round(obs_up_km, 5), 0, "", "", "", "", "", "", "", inc, src_str
                ])

    reach_period_vals = defaultdict(list)
    reach_period_unc = defaultdict(list)
    reach_period_meth = defaultdict(set)
    reach_along_map = defaultdict(list)
    reach_catch_map = {}
    reach_d2m_map = {}
    for rec in imputed_records:
        nid, imp_val, yr, mo, seas = rec[0], rec[1], rec[2], rec[3], rec[4]
        unc_val = rec[6] if len(rec) > 6 else None
        meth_val = rec[7] if len(rec) > 7 else ""
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
            if unc_val is not None:
                reach_period_unc[(riv_id, "all", "all")].append(unc_val)
                reach_period_unc[(riv_id, "year", str(yr))].append(unc_val)
                reach_period_unc[(riv_id, "season", seas)].append(unc_val)
                reach_period_unc[(riv_id, "month", f"{mo:02d}")].append(unc_val)
                reach_period_unc[(riv_id, "year_season", f"{yr}_{seas}")].append(unc_val)
            if meth_val:
                reach_period_meth[(riv_id, "all", "all")].add(meth_val)
                reach_period_meth[(riv_id, "year", str(yr))].add(meth_val)
                reach_period_meth[(riv_id, "season", seas)].add(meth_val)
                reach_period_meth[(riv_id, "month", f"{mo:02d}")].add(meth_val)
                reach_period_meth[(riv_id, "year_season", f"{yr}_{seas}")].add(meth_val)

    periods_csv_path = os.path.join(OUT_DIR, "rsn_reach_periods.csv")
    with open(periods_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "RivID", "CatchID", "period_type", "period_value",
            "obs_up_dist_km", "n_samples", "nitrate_mean", "nitrate_median", "log_nitrate_mean",
            "meas_uncertainty", "meas_method"
        ])
        for (riv_id, p_type, p_val), vals in sorted(reach_period_vals.items()):
            d2m_km = reach_d2m_map[riv_id] / 1000.0
            along_m = statistics.mean(reach_along_map[riv_id])
            obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
            m_v = statistics.mean(vals)
            med_v = statistics.median(vals)
            trans_mv_val = apply_forward_transform(m_v, ttype=trans_type, lam=eff_lambda)
            log_mv = round(trans_mv_val, 5) if trans_mv_val is not None else ""
            u_list = reach_period_unc.get((riv_id, p_type, p_val), [])
            m_set = reach_period_meth.get((riv_id, p_type, p_val), set())
            p_unc = round(math.sqrt(sum(u * u for u in u_list) / len(u_list)), 6) if u_list else ""
            p_meth = sorted(m_set)[0] if m_set else ""
            w.writerow([
                riv_id, reach_catch_map[riv_id], p_type, p_val,
                round(obs_up_km, 5), len(vals), round(m_v, 4), round(med_v, 4), log_mv,
                p_unc, p_meth
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
            "nitrate_mean", "nitrate_median", "log_nitrate_mean",
            "meas_uncertainty", "meas_method", "n_samples"
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
                all_uncs = [u for s in s_list for u in s[4]]
                all_meths = {m for s in s_list for m in s[5] if m}
                along_m = statistics.mean(s[2] for s in s_list)
                snap_d = min(s[1] for s in s_list)
                obs_up_km = max(0.0, d2m_km - along_m / 1000.0)
                m_no3 = statistics.mean(all_vals)
                med_no3 = statistics.median(all_vals)
                trans_no3_val = apply_forward_transform(m_no3, ttype=trans_type, lam=eff_lambda)
                log_no3 = round(trans_no3_val, 5) if trans_no3_val is not None else ""
                r_unc = round(math.sqrt(sum(u * u for u in all_uncs) / len(all_uncs)), 6) if all_uncs else ""
                r_meth = sorted(all_meths)[0] if all_meths else ""
                w.writerow([
                    riv_id, r[3], r[4], parent_dn.get(riv_id, 0),
                    round(d2m_km, 5), round(up_mid_km, 5), round(l_km, 5), round(r[6] / 1000.0, 5),
                    round(log_afv[riv_id], 6), r[7], r[8],
                    round(sx, 1), round(sy, 1), round(mx, 1), round(my, 1), round(ex, 1), round(ey, 1),
                    round(rc[0] or 0.0, 3), round(rc[1] or 0.0, 4), round(rc[2] or 1.0, 4), round(rc[5] or rc[0] or 0.0, 3),
                    is_chalk_lime, is_sand, is_aq, is_peat,
                    1, s_list[0][0], round(snap_d, 2), round(obs_up_km, 5),
                    round(m_no3, 4), round(med_no3, 4), log_no3,
                    r_unc, r_meth, len(all_vals)
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
                    "", "", "", "", "", 0
                ])

    active_meas_method = "gam" if any("gam" in s[5] for s_list in site_by_reach.values() for s in s_list) else None
    meta_path = os.path.join(OUT_DIR, "network_metadata.json")
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
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta_obj, f, indent=2)

    write_outputs_to_dirs(net_abbr, {
        "england_rsn_reaches.csv": out_reaches_path,
        "rsn_sites.csv": sites_csv_path,
        "rsn_reach_periods.csv": periods_csv_path,
        "rsn_catchments_lookup.csv": catch_lu_path,
        "network_metadata.json": meta_path,
    }, source=source_key)
    print(f"Done in {time.time() - t0:.2f} seconds (saved to {OUT_DIR}/).")


if __name__ == "__main__":
    main()
