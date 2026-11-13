# 🌊 Spatial Stream Network (`S3N`) Water Quality & Ecological Modelling

[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B%20(stdlib%20only)-3776AB?logo=python&logoColor=white)](docs/setup_and_architecture.md)
[![R 3.5+](https://img.shields.io/badge/R-3.5%2B-276DC3?logo=r&logoColor=white)](docs/setup_and_architecture.md)
[![C++ OpenMP](https://img.shields.io/badge/C%2B%2B-OpenMP%20Engine-00599C?logo=cplusplus&logoColor=white)](src/s3n_fast_neighbors.cpp)
[![Base R Graphics](https://img.shields.io/badge/Plots-Pure%20Base%20R-16A085)](docs/usage_guide.md)

A high-speed, headless pipeline for **Spatial Stream Network (`S3N` / NNGP)** regression, 5-fold cross-validation, universal stream kriging, and diagnostic reporting across national and regional river networks (defaulting to the Ordnance Survey Open Rivers Network v2, `ORN v2`, across England).

Designed to run cleanly on Linux servers, laptops, and Databricks clusters **without heavy GIS/plotting dependencies** (`GDAL`, `geopandas`, `sf`, `SSN2`, or `ggplot2`).

---

## ✨ Key Features

- ⚡ **Sub-Second C++/OpenMP Stream-Tree & Joint NNGP Engine** ([`src/s3n_fast_neighbors.cpp`](src/s3n_fast_neighbors.cpp)): Replaces slow R tree walks with an Euler-tour LCA & Nearest-Neighbor Gaussian Process (NNGP) builder, multi-kernel profile MLE solver (`exponential`, `spherical`, `mariah`, `linear`, `auto`), physical confluence `AFV` weighting (`length`, `shreve`, `hydro_yield`), and **Joint Additive `Tail-Up + Tail-Down + Euclidean` Universal Kriging** (`--combine=joint`), processing ~91,000 river reaches in under a second.
- 🐍 **Zero-Dependency Python GIS Extractor** ([`scripts/01_preprocess_orn.py`](scripts/01_preprocess_orn.py)): Reads multi-gigabyte OGC GeoPackages (`.gpkg`) via Python's built-in `sqlite3` and binary WKB decoder, builds directed river trees, computes Strahler/Shreve orders and Additive Function Values (`log_AFV`), snaps stations (`≤ 250 m`), and auto-generates `configs/determinands.json`.
- 🧪 **Censored & Interval-Range Imputation**: Supports Robust Lognormal Regression on Order Statistics (`ROS`) for `<LOD` values and multiple range-handling strategies (`tobit`, `geometric`, `midpoint`, `lower`, `upper`, or `ros_site`).
- 📐 **12 Response Transformations & 4 Normalisation Methods**: Built-in `log`, `log10`, `log1p`, `sqrt`, `cbrt`, `fourthrt`, `reciprocal`, `boxcox`, `yeojohnson`, `asin`, `logit`, and `none` transforms (with Taylor/Duan back-transformation), plus automatic `zscore`, `robust_zscore`, `minmax`, and `percentile` outputs.
- 📉 **Two-Stage GAM & Measurement-Error Kriging**: Propagate site-level temporal/GAM standard errors (`meas_uncertainty`) directly into spatial variance estimation and kriging uncertainty via `--meas-error` / `--gam`.
- 🌲 **Non-Linear Spline (`GAM + S3N`) & Tree-Ensemble (`RF / GBM + S3N`) Regression Kriging**: Capture non-linear environmental thresholds and high-order catchment interactions using REML-penalized cubic regression splines (`--trend=gam`), GCV natural cubic splines (`--trend=spline`), Out-of-Bag Random Forests (`--trend=rf`), or cross-fitted Gradient Boosted Trees (`--trend=gbm`) combined with S3N residual stream kriging.
- 🗺️ **Multi-Source, Multi-Network & Temporal Batching**: Combine built-in (`rsn`, `ssn`, `rss`) and custom CSV data sources (`+` / `,` delimiters via [`scripts/format_custom_data.py`](scripts/format_custom_data.py)), bring custom river networks (`BYOWN`), customize `15` covariates with `--vars=`, and batch-fit by `--by=year|season|month`.
- 📊 **Automated Diagnostics, Statistics & Maps**: Exports publication-grade Base-R maps (`fig1`–`fig5`), station snapping diagnostics ([`06_station_snapping_diagnostics.R`](scripts/06_station_snapping_diagnostics.R)), single-run environmental gradient reports ([`05_descriptive_stats.R`](scripts/05_descriptive_stats.R)), and multi-run / `N:P` stoichiometry comparisons ([`04_compare_model_runs.R`](scripts/04_compare_model_runs.R)).
- 🔗 **Multivariate Stream Co-Kriging** ([`scripts/run_s3n_cokriging.R`](scripts/run_s3n_cokriging.R)): Jointly models $P \ge 2$ heterotopic chemical and ecological targets via a dedicated OpenMP **Factor-LMC Joint NNGP** engine ([`docs/multivariate_cokriging_guide.md`](docs/multivariate_cokriging_guide.md)), reporting both All-Holdout and Auxiliary-Aided CV gains.

---

## 🔄 Pipeline Architecture

```mermaid
flowchart LR
    Raw["Raw Inputs: gis/*.gpkg + data/*.csv"] -->|"Optional CSV formatting"| Custom["format_custom_data.py: sites, samples + JSON"]
    Raw --> Pre["01_preprocess_orn.py: Topology, ROS, Snapping"]
    Custom --> Pre
    Pre --> Proc["data/processed/net_source/ (5 standardized files)"]
    Proc --> Snap["06_station_snapping_diagnostics.R -> output/station_snapping/"]
    Proc --> Run["run_s3n_determinand.R (02_fit + 03_plots + C++ OpenMP)"]
    Run --> Res["output/results/det_code_net_source/ (predictions, CV, params, fig1-fig5)"]
    Res --> Desc["05_descriptive_stats.R -> run_dir/descriptives/"]
    Res --> Comp["04_compare_model_runs.R -> output/comparisons/"]
```

---

## 🛠️ Setup & Installation

### 1. System Prerequisites

**Common (🐧 Linux & 🪟 Windows):**
- **Python ≥ 3.8** — uses the standard library only (**zero `pip` packages required**).
- **R ≥ 3.5.0** (on Windows, ensure `Rscript.exe` is on `PATH`).

**🐧 Linux:**
- Development headers, OpenMP C/C++ compiler, BLAS/LAPACK, Cairo, and SQLite3:
  - **Fedora / RHEL**: `sudo dnf install gcc gcc-c++ libgomp R-core-devel openblas-devel lapack-devel cairo-devel sqlite-libs`
  - **Ubuntu / Debian**: `sudo apt install build-essential libomp-dev r-base-dev libopenblas-dev liblapack-dev libcairo2-dev libsqlite3-0`

**🪟 Windows:**
- Install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)** matching your R version (provides the MinGW/GCC C++ compiler and OpenMP needed to install `S3N` from GitHub and compile [`src/s3n_fast_neighbors.cpp`](src/s3n_fast_neighbors.cpp)).
- *(Alternative)* Use **WSL2** (Windows Subsystem for Linux), which runs the 🐧 Linux setup—including the precompiled `src/s3n_fast_neighbors.so`—natively without Rtools.

### 2. Install R Packages
```r
install.packages(c("RANN", "rdist", "matrixStats", "pbapply", "dplyr", "jsonlite", "remotes"))
remotes::install_github("jpierkunke/S3N")
```

### 3. C++ OpenMP Neighbor Engine (`src/s3n_fast_neighbors.cpp`)
- **🐧 Linux (`x86_64`)**: A precompiled binary ([`src/s3n_fast_neighbors.so`](src/s3n_fast_neighbors.so)) is already included in the repository. You only need to recompile if you edit the C++ source or hit an ABI mismatch on a different Linux distro/CPU architecture:
  ```bash
  R CMD SHLIB src/s3n_fast_neighbors.cpp
  ```
- **🪟 Windows / macOS / Other Architectures**: A Linux `.so` binary cannot be loaded on Windows, so compile a native `.dll` using Rtools from the repository root (`src/Makevars.win` applies OpenMP flags automatically):
  ```powershell
  R CMD SHLIB src/s3n_fast_neighbors.cpp
  ```
  R automatically loads `s3n_fast_neighbors.dll` on Windows (`.Platform$dynlib.ext`) once compiled, or you can run inside WSL2 using the precompiled 🐧 Linux `.so`.

### 4. Configure Data & Output Paths *(Optional)*
By default, scripts read from `gis/`, `data/`, `data/processed/`, and `output/` as defined in [`configs/directories.json`](configs/directories.json). You can override these via `configs/directories.json`, environment variables (`S3N_GIS_DIR`, `S3N_DATA_DIR`, `S3N_PROCESSED_DIR`, `S3N_OUTPUT_DIR`), or CLI flags (`--gis-dir=`, `--data-dir=`, `--processed-dir=`, `--output-dir=`).

---

## 🚀 Quickstart

### Step 1 — Preprocess the River Network & Monitoring Data
```bash
# Default: ORN v2 network + EA River Surveillance Network (auto-creates configs/determinands.json)
python3 scripts/01_preprocess_orn.py --network=orn --source=rsn

# Or combine RSN + Small Streams Network (SSN)
python3 scripts/01_preprocess_orn.py --network=orn --source=rss
```

### Step 2 — Inspect Station-to-Reach Snapping *(Optional)*
```bash
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Thames
```

### Step 3 — Fit S3N Models, Krige & Render Maps
```bash
# Nitrate (0117) with square-root transform and robust Z-score maps
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --norm=robust_zscore

# Auto-select best stream covariance kernel (Exponential, Spherical, Mariah, Linear) + Hydrological-Yield AFV
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --kernel=auto --afv=hydro_yield --combine=joint

# Penalized GAM Smooths + S3N or Tree-Ensemble Regression Kriging (Random Forest / Gradient Boosting + S3N)
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=gam
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=rf --trees=150

# Orthophosphate (0180) with default log transform
Rscript scripts/run_s3n_determinand.R 0180

# Temporal batch mode (fits all 4 meteorological seasons + combined temporal grid)
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --by=season
```

### Step 4 — Generate Statistical Reports & Multi-Run Comparisons
```bash
# Single-run descriptive stats & landscape/stream-order gradient tests
Rscript scripts/05_descriptive_stats.R output/results/det_0117_orn_rsn_trans_sqrt

# Compare Nitrate vs. Orthophosphate (includes N:P stoichiometry & Wald Z-tests)
Rscript scripts/04_compare_model_runs.R \
  output/results/det_0117_orn_rsn_trans_sqrt \
  output/results/det_0180_orn_rsn \
  --labels="Nitrate (0117),Orthophosphate (0180)"
```

> 💡 **Split workflow tip**: You can run model fitting ([`scripts/02_fit_s3n_models.R`](scripts/02_fit_s3n_models.R)) and plotting ([`scripts/03_generate_plots.R`](scripts/03_generate_plots.R)) separately so you can re-style colorbars, basins (`--basin=Eden`), or `--norm=` views without re-fitting models.

---

## 📦 Key Outputs

| Stage / Script | Output Directory | Key Files Produced |
| :--- | :--- | :--- |
| **Preprocessing** (`01`) | `data/processed/<net>_<source>/` | `reaches.csv`, `rsn_sites.csv`, `rsn_site_reach_mapping.csv`, `catchment_lookup.csv`, `network_metadata.json` |
| **Snapping Diagnostics** (`06`) | `output/station_snapping/<net>_<source>/` | `station_snapping_all.csv`, `catchment_snapping_summary.csv`, `station_snapping_report.md`, `fig_snap_*.png` |
| **S3N Fit & Plots** (`run` / `02`+`03`) | `output/results/det_<code>_<net>_<src>.../` | `predictions_*.csv` (+ normalised columns), `cv_comparison_*.csv`, `parameters_*.csv`, `fig1`–`fig5_*.png` |
| **Descriptive Stats** (`05`) | `<run_dir>/descriptives/` | `descriptive_report.md`, `descriptive_stats_summary.csv`, `stratified_results_summary.csv`, `within_run_statistical_tests.csv`, `fig_desc1`–`2.png` |
| **Run Comparison** (`04`) | `output/comparisons/<comparison_slug>/` | `comparison_report.md`, `cv_metrics_comparison.csv`, `statistical_tests_summary.csv`, `reach_relationship_summary.csv`, `parameter_comparison.csv`, `fig_comp1`–`5.png` |

---

## 📂 Repository Structure

<details>
<summary>Click to expand directory layout</summary>

```text
spatial_stream_model/
├── configs/
│   ├── directories.json              # Default gis, data, processed, and output roots
│   ├── determinands.json             # Auto-generated determinand names, short labels & units
│   ├── networks/orn.json             # Built-in ORN v2 network profile
│   └── sources/                      # Monitoring source configs (rsn, ssn, composites)
├── data/                             # Raw EA monitoring archives & processed/<net>_<source>/
├── docs/                             # Full guides & theoretical documentation (see docs/README.md)
├── gis/                              # Raw river network GeoPackages (e.g. ORN_v2.gpkg)
├── output/                           # Generated results/, station_snapping/, and comparisons/
├── scripts/
│   ├── 01_preprocess_orn.py          # GeoPackage extractor, ROS imputation & station snapping
│   ├── format_custom_data.py         # Formatter for user-supplied CSV datasets
│   ├── run_s3n_determinand.R         # End-to-end S3N fit, 5-fold CV, kriging & plotting
│   ├── run_s3n_cokriging.R           # Multivariate Factor-LMC Joint NNGP Co-Kriging runner
│   ├── 02_fit_s3n_models.R           # Standalone model fitter (exports .rds checkpoint)
│   ├── 03_generate_plots.R           # Standalone map & diagnostic plotter
│   ├── 04_compare_model_runs.R       # Multi-run & cross-determinand statistical comparison
│   ├── 05_descriptive_stats.R        # Single-run distributional & gradient diagnostics
│   ├── 06_station_snapping_diagnostics.R
│   └── s3n_*.py / s3n_*.R            # Modular config, covariate, transform, model, plot & stats helpers
├── src/
│   ├── Makevars                      # OpenMP compiler flags
│   └── s3n_fast_neighbors.cpp        # OpenMP C++ Euler-tour & NNGP neighbor engine
└── tests/
    ├── test_python_utils.py          # Python unit tests
    └── test_r_and_cpp.R              # R modules & C++ shared library unit tests
```

</details>

---

## 🧪 Running Automated Tests

Ensure `src/s3n_fast_neighbors.so` is compiled and `configs/determinands.json` exists (created automatically by any `01_preprocess_orn.py` invocation, e.g. `--list-sources`):

```bash
# Python utilities, WKB decoding, ROS imputation & custom formatter tests
python3 -m unittest tests/test_python_utils.py

# R module helpers (config, transforms, covariates, stats) & C++ OpenMP engine tests
Rscript tests/test_r_and_cpp.R
```

---

## 📚 Documentation

See the **[Documentation Hub](docs/README.md)** or jump straight to a guide:

| Category | Guides |
| :--- | :--- |
| **Getting Started & CLI** | [Setup, Architecture & Testing](docs/setup_and_architecture.md) · [End-to-End Usage Guide](docs/usage_guide.md) · [Running on Databricks](docs/databricks_guide.md) |
| **Specialised Workflows** | [Regulatory Exceedance & WFD Risk Mapping](docs/regulatory_exceedance_mapping_guide.md) · [Spatio-Temporal Trends & Borrowing](docs/spatiotemporal_trends_and_borrowing_guide.md) · [Multivariate Stream Co-Kriging](docs/multivariate_cokriging_guide.md) · [Non-Linear (`GAM + S3N`) & Tree-Ensemble (`RF/GBM`) Guide](docs/nonlinear_and_tree_ensemble_trend_guide.md) |
| **Customising the Pipeline** | [Custom Monitoring Sources](docs/custom_data_sources_guide.md) · [Custom Water Networks (`BYOWN`)](docs/custom_water_network_guide.md) · [ORN Covariates Reference (`--vars`)](docs/orn_variables_reference.md) |
| **Diagnostics & Interpretation** | [Station Snapping Guide](docs/station_snapping_guide.md) · [Descriptive Statistics Guide](docs/descriptive_stats_interpretation_guide.md) · [Model Comparison Guide](docs/model_comparison_interpretation_guide.md) |
| **Methodology & Theory** | [S3N Modelling & Kriging Guide](docs/s3n_modelling_and_kriging_guide.md) · [Two-Stage GAM & Measurement-Error Kriging](docs/two_stage_gam_measurement_error_kriging_guide.md) · [Preprocessing & C++ Engine](docs/preprocessing_pipeline.md) |
