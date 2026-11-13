# Setup, Architecture, Configuration & Testing

This guide is the single reference for repository layout, environment dependencies, compiling the C++ OpenMP extension, JSON/environment configuration, and running automated tests.

---

## Table of Contents

- [1. Repository & Output Layout](#1-repository--output-layout)
- [2. Software & Hardware Requirements](#2-software--hardware-requirements)
  - [2.1 Python Environment (`01_preprocess_orn.py`, `format_custom_data.py`)](#21-python-environment-01_preprocess_ornpy-format_custom_datapy)
  - [2.2 R Environment](#22-r-environment)
  - [2.3 Compiled C++ / OpenMP Library (`src/s3n_fast_neighbors.cpp`)](#23-compiled-c--openmp-library-srcs3n_fast_neighborscpp)
  - [2.4 Linux System Dependencies](#24-linux-system-dependencies)
  - [2.5 Hardware Profile](#25-hardware-profile)
- [3. Configuration Files & Environment Variables](#3-configuration-files--environment-variables)
  - [3.1 Directory Paths (`configs/directories.json`)](#31-directory-paths-configsdirectoriesjson)
  - [3.2 Determinand Metadata Catalog (`configs/determinands.json`)](#32-determinand-metadata-catalog-configsdeterminandsjson)
- [4. Running the Automated Test Suite](#4-running-the-automated-test-suite)
- [5. License](#5-license)

---

## 1. Repository & Output Layout

```text
spatial_stream_model/
├── configs/
│   ├── directories.json              # Default data, gis, processed, and output roots
│   ├── determinands.json             # (Optional user file) Determinand labels & units
│   ├── networks/
│   │   └── orn.json                  # Built-in ORN v2 network profile
│   └── sources/
│       ├── rsn.json                  # Environment Agency River Surveillance Network
│       ├── ssn.json                  # Environment Agency Small Streams Network
│       └── composites.json           # Built-in composite sources (rss, all)
├── data/
│   ├── NCEA_River_Surveillance_Network_ARD_Open/
│   ├── NCEA_Small_Streams_Network_ARD_Open/
│   └── processed/
│       └── <net>_<source>/           # Standardized 5-file output from 01_preprocess_orn.py
├── gis/
│   └── ORN_v2.gpkg                   # National river network GeoPackage (~10.4 GB)
├── output/
│   ├── results/
│   │   └── det_<code>_<net>_<src>.../# Model predictions, CV tables, parameters, fig1–fig5
│   │       └── descriptives/         # Single-run report & figures from 05_descriptive_stats.R
│   ├── station_snapping/
│   │   └── <net>_<source>/           # Snapping tables & maps from 06_station_snapping_diagnostics.R
│   └── comparisons/
│       └── <comparison_slug>/        # Multi-run reports & figures from 04_compare_model_runs.R
├── scripts/
│   ├── 01_preprocess_orn.py          # Network extraction, ROS imputation, station snapping
│   ├── format_custom_data.py         # Helper to convert user CSVs into source format + JSON
│   ├── run_s3n_determinand.R         # Main end-to-end R runner (fit + CV + kriging + plots)
│   ├── 02_fit_s3n_models.R           # Standalone model fitting step (saves .rds checkpoint)
│   ├── 03_generate_plots.R           # Standalone plotting step (reads .rds checkpoint)
│   ├── 04_compare_model_runs.R       # Multi-run & cross-determinand comparison
│   ├── 05_descriptive_stats.R        # Single-run descriptive stats & gradient tests
│   ├── 06_station_snapping_diagnostics.R
│   ├── s3n_py_utils.py               # Shared Python utilities (WKB, ROS, config, snapping)
│   └── s3n_*.R                       # Shared R modules (config, covariates, transforms, model, plot, stats)
├── src/
│   ├── Makevars                      # OpenMP compiler flags for R CMD SHLIB
│   └── s3n_fast_neighbors.cpp        # OpenMP C++ Euler-tour & NNGP neighbor engine
├── tests/
│   ├── test_python_utils.py          # Python unit tests
│   └── test_r_and_cpp.R              # R modules + C++ shared library unit tests
└── docs/                             # Documentation (see docs/README.md)
```

---

## 2. Software & Hardware Requirements

The pipeline minimizes external dependencies so it runs reliably in headless Linux environments without heavy GIS or plotting stacks (`GDAL`, `geopandas`, `sf`, `SSN2`, or `ggplot2`).

### 2.1 Python Environment (`01_preprocess_orn.py`, `format_custom_data.py`)
- **Python Version**: **Python ≥ 3.8** (uses `statistics.NormalDist()` for ROS normal quantiles).
- **Third-Party Packages**: **None (Zero `pip` or `conda` dependencies)** — uses only Python standard library modules (`sqlite3`, `struct`, `statistics`, `math`, `collections`, `csv`, `json`, `argparse`, `os`, `time`).

### 2.2 R Environment
- **R Version**: **R ≥ 3.5.0**.
- **Required R Packages**:
  | Package | Source | Purpose |
  | :--- | :--- | :--- |
  | **`S3N`** | GitHub (`remotes::install_github("jpierkunke/S3N")`) | Scalable Spatial Stream Network NNGP estimation (`BRISC_estimation_stream`), kriging (`BRISC_prediction_stream`), and parametric bootstrap (`BRISC_bootstrap`). |
  | `RANN`, `rdist`, `matrixStats`, `pbapply`, `dplyr` | CRAN | Runtime dependencies of `jpierkunke/S3N`. |
  | `jsonlite` *(optional / recommended)* | CRAN | Fast JSON parsing for `configs/*.json` and `network_metadata.json` (with built-in base-R fallback in [`scripts/s3n_config.R`](../scripts/s3n_config.R)). |
  | `stats`, `graphics`, `grDevices`, `utils`, `parallel` | Base R | Linear algebra, statistics, and headless PNG plotting (`type = "cairo"`). |

### 2.3 Compiled C++ / OpenMP Library (`src/s3n_fast_neighbors.cpp`)
R scripts dynamically load `src/s3n_fast_neighbors.so`. [`src/Makevars`](../src/Makevars) already configures `$(SHLIB_OPENMP_CXXFLAGS)`, so compile from the repository root with:

```bash
R CMD SHLIB src/s3n_fast_neighbors.cpp
```

### 2.4 Linux System Dependencies
| Dependency | Fedora / RHEL | Ubuntu / Debian | Purpose |
| :--- | :--- | :--- | :--- |
| **C/C++ Compiler & OpenMP** | `gcc`, `gcc-c++`, `libgomp` | `build-essential`, `libomp-dev` | Compiling `src/s3n_fast_neighbors.cpp` and `S3N`. |
| **R Dev Headers & BLAS/LAPACK** | `R-core-devel`, `openblas-devel`, `lapack-devel` | `r-base-dev`, `libopenblas-dev`, `liblapack-dev` | `<R.h>`, LAPACK, and BLAS headers. |
| **Cairo Graphics** | `cairo`, `cairo-devel` | `libcairo2`, `libcairo2-dev` | Headless anti-aliased PNG rendering (`type = "cairo"`). |
| **SQLite3** | `sqlite-libs` | `libsqlite3-0` | Underlying C library used by Python `sqlite3` to read `.gpkg` files. |

### 2.5 Hardware Profile
- **RAM**: ~1.5 GB peak memory for a full national run.
- **CPU**: Defaults to `4` OpenMP threads (`--omp=4`); a full national + 3-basin ORN run with 5-fold CV and 4 figures completes in ~18–25 seconds.

---

## 3. Configuration Files & Environment Variables

### 3.1 Directory Paths (`configs/directories.json`)
All Python and R scripts resolve root paths in this precedence order (highest to lowest):
1. **CLI flags**: `--gis-dir=`, `--data-dir=`, `--processed-dir=`, `--output-dir=`
2. **Environment variables**: `S3N_GIS_DIR`, `S3N_DATA_DIR`, `S3N_PROCESSED_DIR`, `S3N_OUTPUT_DIR`
3. **`configs/directories.json`**:
   ```json
   {
     "gis_dir": "gis",
     "data_dir": "data",
     "processed_dir": "data/processed",
     "output_dir": "output"
   }
   ```

### 3.2 Determinand Metadata Catalog (`configs/determinands.json`)
`load_determinand_catalog()` in [`scripts/s3n_config.R`](../scripts/s3n_config.R) and [`scripts/s3n_py_utils.py`](../scripts/s3n_py_utils.py) reads `configs/determinands.json` to map determinand codes (e.g., `"0117"`, `"whpt_aspt"`) to human-readable plot titles, short labels, and units.

On first run, [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) automatically generates `configs/determinands.json` by scanning available RSN/SSN Water Quality and LCMS/GCMS CSVs alongside the built-in ecological index definitions (preserving any manual edits on later runs; pass `--update-determinands` to rescan for new codes):

```json
{
  "0117": {
    "short": "Nitrate",
    "name": "Nitrate as N",
    "unit": "mg/L"
  },
  "0180": {
    "short": "Orthophosphate",
    "name": "Orthophosphate as P",
    "unit": "mg/L"
  },
  "whpt_aspt": {
    "short": "WHPT ASPT",
    "name": "WHPT Average Score Per Taxon",
    "unit": "ASPT"
  }
}
```

---

## 4. Running the Automated Test Suite

The `tests/` directory contains fast unit tests covering Python preprocessing/formatting utilities, R configuration/covariate/transform helpers, and the compiled C++ OpenMP neighbor engine:

```bash
# 1. Python unit tests (WKB decoding, ROS imputation, source tags, custom data formatting)
python3 -m unittest tests/test_python_utils.py

# 2. R & C++ unit tests (requires src/s3n_fast_neighbors.so and configs/determinands.json)
Rscript tests/test_r_and_cpp.R
```

> [!NOTE]
> Both test suites check determinand metadata resolution (e.g., `"0117"` → `"Nitrate"`, `"whpt_aspt"` → `"ASPT"`), so on a fresh clone run `python3 scripts/01_preprocess_orn.py --list-sources` once if `configs/determinands.json` has not been generated yet.

---

## 5. License

This repository's source code is licensed under the **GNU General Public License v3.0 (`GPL-3.0`)** — see [`LICENSE`](../LICENSE). External river networks and monitoring datasets placed under `gis/` and `data/` are downloaded separately from their upstream providers and remain governed by their own licences and attribution terms.
