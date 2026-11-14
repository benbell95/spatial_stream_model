# Scalable Spatial Stream Network Modelling

[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-3776AB?logo=python&logoColor=white)](docs/setup_and_architecture.md)
[![R 3.5+](https://img.shields.io/badge/R-3.5%2B-276DC3?logo=r&logoColor=white)](docs/setup_and_architecture.md)
[![C++ OpenMP](https://img.shields.io/badge/C%2B%2B-OpenMP%20Engine-00599C?logo=cplusplus&logoColor=white)](src/s3n_fast_neighbors.cpp)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

This project implements a headless pipeline for **Scalable Spatial Stream Network (`S3N` / NNGP)** regression, 5-fold cross-validation, universal stream kriging, and diagnostic reporting, using open-access UK Environment Agency (EA) national river network data, and the Ordnance Survey Open Rivers Network, `ORN v2` riverline dataset across England.

These scripts are fast, and are designed to run on any system without heavy GIS/plotting dependencies (no `GDAL`, `geopandas`, `sf`, `SSN2`, or `ggplot2` packages are required).

Also supports using other environmental (river) datasets including citizen science data, and the use of other riverline datasets, and is not limited to the UK.

> **Attribution:** This pipeline builds on and adapts the **Scalable Spatial Stream Network (`S3N`)** framework ([`jpierkunke/S3N`](https://github.com/jpierkunke/S3N)), which brings linear Nearest-Neighbour Gaussian Process (`NNGP`) scaling to river networks. Full credit belongs to the authors of these packages and theory on which this repository builds:
> - **`S3N` R Package** ([https://github.com/jpierkunke/S3N](https://github.com/jpierkunke/S3N)) — **Jessica P. Kunke, Julian D. Olden, and Tyler H. McCormick**:
>   Kunke, J. P., Olden, J. D., & McCormick, T. H. (2025). *Scalable Spatial Stream Network (S3N) Models*. **arXiv:2512.12398 [stat.ME]**. [https://arxiv.org/abs/2512.12398](https://arxiv.org/abs/2512.12398)
> - **`BRISC` R Package** — **Arkajyoti Saha and Abhirup Datta**:
>   Saha, A., & Datta, A. (2018). *BRISC: Bootstrap for Rapid Inference on Spatial Covariances*. *Stat*, 7(1), e184.
>
> **If you use this pipeline in published work or downstream analysis, please cite both the original `S3N` paper (Kunke et al., 2025, `arXiv:2512.12398`) and the `BRISC` paper (Saha & Datta, 2018, *Stat*, 7(1), e184).**
>
> See [Where This Pipeline Differs from `S3N`](docs/where_this_differs_from_s3n.md) for a full breakdown of what is shared, replaced, or added.

<img alt="model_results" src="https://github.com/user-attachments/assets/fa6210b8-aeb3-47d1-89a8-38fd6ec1a756" />

## 📑 Table of Contents

- [Key Features](#-key-features)
- [Setup & Installation](#-setup--installation)
  - [1. System Prerequisites](#1-system-prerequisites)
  - [2. Quick Setup (R Packages, C++ Library & Paths)](#2-quick-setup-r-packages-c-library--paths)
  - [3. Required Data Downloads](#3-required-data-downloads)
- [Quickstart](#-quickstart)
  - [Step 1: Preprocess the River Network & Monitoring Data](#step-1-preprocess-the-river-network--monitoring-data)
  - [Step 2: Inspect Station-to-Reach Snapping](#step-2-inspect-station-to-reach-snapping-optional)
  - [Step 3: Fit S3N Models, Krige & Render Maps](#step-3-fit-s3n-models-krige--render-maps)
  - [Step 4: Generate Statistical Reports & Multi-Run Comparisons](#step-4-generate-statistical-reports--multi-run-comparisons)
- [Full Documentation](#-full-documentation)
- [AI Disclosure](#-ai-disclosure)
- [License](#-license)


## ✨ Key Features

- ⚡ Sub-Second C++/OpenMP Stream-Tree & Joint NNGP Engine
- 🐍 Low-Dependency Python GIS Extractor 
- 🧪 Censored & Interval-Range Imputation 
- 📐 12 Response Transformations & 4 Normalisation Methods
- 🔗 Univariate Kriging and Multivariate Co-Kriging
- 📉 Non-Linear Spline & Tree-Ensemble Regression Kriging, and Two-Stage GAM & Measurement-Error Kriging
- 🗺️ Multi-Source, Multi-Network & Temporal Batching (Combine built-in (`rsn`, `ssn`) and custom data sources)
- 📊 Automated Diagnostics, Statistics & Maps


## 🛠️ Setup & Installation

### 1. System Prerequisites

- **Python ≥ 3.8** — uses Python standard library only for core `.gpkg`/CSV preprocessing (optional `pyproj` for custom CRS reprojection, `pyarrow` for GeoParquet).
- **R ≥ 3.5.0** with `S3N`, its dependencies, and optional `data.table` (strongly recommended for 5–15× faster CSV I/O). 

On Linux, install C++/OpenMP, BLAS/LAPACK, Cairo, and SQLite dev headers (if needed); on Windows, install **[Rtools](https://cran.r-project.org/bin/windows/Rtools/)** (or use **WSL2**).

### 2. Quick Setup (R Packages, C++ Library & Paths)
```r
# Install required/recommended R dependencies
install.packages(c("RANN", "rdist", "matrixStats", "pbapply", "dplyr", "jsonlite", "remotes", "data.table"))
# Install S3N package from GitHub
remotes::install_github("jpierkunke/S3N")
```
⚠️ Windows and MacOS users must compile the OpenMP C++ engine (`src/s3n_fast_neighbors.cpp`) in order to run these scripts. (This comes precompiled on Linux). More detailed information available in the [setup guide](docs/setup_and_architecture.md#13-compiled-c--openmp-library-srcs3n_fast_neighborscpp). 

Run these commands from your terminal in the repository root directory.

```bash
# Compile OpenMP C++ engine if needed
# macOS (--preclean ensures any prebuilt Linux .so is rebuilt cleanly)
R CMD SHLIB --preclean src/s3n_fast_neighbors.cpp
```

```powershell
# Windows PowerShell (omit full path if R's bin\x64 directory is on PATH; requires Rtools)
& "C:\Program Files\R\R-<ver>\bin\x64\R.exe" CMD SHLIB src/s3n_fast_neighbors.cpp
```

Default workspace folders (`gis/`, `data/`, `output/`) are configured in [`configs/directories.json`](configs/directories.json) and can be overridden via environment variables or CLI flags (`--gis-dir=`, `--data-dir=`, `--processed-dir=`, `--output-dir=`). See [Setup, Architecture & Testing](docs/setup_and_architecture.md) for full OS package lists and configuration details.

### 3. Required Data Downloads

💾 Download and unzip the three external datasets into `gis/` and `data/` before running the default pipeline (or pass `--source=rsn` if omitting the optional `SSN` dataset):

1. **OS Open Rivers Network (`ORN` v2.0 — Required)**: Download **version 2.0** in **GeoPackage** format (`ORN_v2_GeoPackage.zip`) from [openrivers.net/download/](https://openrivers.net/download/) and extract `ORN_v2.gpkg` into `gis/ORN_v2.gpkg`.
2. **River Surveillance Network (`RSN` — Required)**: Download `NCEA_River_Surveillance_Network_ARD_Open.zip` from [Environment Agency NCEA — RSN Analysis Ready Data](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/21763290-f395-4ca1-bcf6-884c2e7548d4?q=analysis+ready+data#access) and extract its CSV files directly into `data/NCEA_River_Surveillance_Network_ARD_Open/`.
3. **Small Streams Network (`SSN` — Optional, recommended)**: Download `NCEA_Small_Streams_Network_ARD_Open.zip` from [Environment Agency NCEA — SSN Analysis Ready Data](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/98aca14e-3736-4329-b4b6-ae7c36fa5db1?q=small+streams+network&jry=qs&pg=1&rpp=20&srt=most_relevant#access) and extract its CSV files directly into `data/NCEA_Small_Streams_Network_ARD_Open/`.

See [Section 2 of Setup, Architecture & Testing](docs/setup_and_architecture.md#2-obtaining--placing-the-required-data) for further details.


## 🚀 Quickstart

ℹ️ These scripts are designed to run "headless" and are called from the terminal, rather than from inside R or Python environments. Run these commands from your terminal. Windows users will need to ensure R is on PATH to use these commands. Alternatively, you may run from inside R/Python if you prefer, see the [Running on Databricks](docs/databricks_guide.md) guide for details (this can be adapted for your personal machine/setup).

> [!TIP]
> Open [`cli_workbench/index.html`](cli_workbench/index.html) for a visual command builder directly in your web browser. Pick your analysis options from grouped dropdowns and copy ready-to-run commands formatted for **Bash**, **Windows PowerShell**, or **Windows Command Prompt**.

<img alt="pipeline_command" src="https://github.com/user-attachments/assets/c2c0b5ce-2741-4338-a890-acda8a3e3ba1" />

### Step 1: Preprocess the River Network & Monitoring Data

You must first pre-process the data. This will happen automatically if running `scripts/run_s3n_determinand.R`, but running the pre-processing scripts allows you to customise the network.

```bash
# Default: ORN v2 network + EA River Surveillance Network (auto-creates configs/determinands.json)
python3 scripts/01_preprocess_orn.py --network=orn --source=rsn

# Or combine RSN + Small Streams Network (SSN)
python3 scripts/01_preprocess_orn.py --network=orn --source=rss
```

### Step 2: Inspect Station-to-Reach Snapping *(Optional)*

This allows you to check the station snapping to the riverline network (choose which catchments to look at).

```bash
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Thames
```

### Step 3: Fit S3N Models, Krige & Render Maps

This is the main script runner to run the models. Use the [`cli_workbench/index.html`](cli_workbench/index.html) for a visual command builder.

```bash
# Nitrate (0117) with square-root transform
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt

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

### Step 4: Generate Statistical Reports & Multi-Run Comparisons

These scripts generate statistical analysis and comparison reports and plots.

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


## 📚 Full Documentation

See the **[Documentation Hub](docs/README.md)** or jump straight to a guide:

| Category | Guides |
| :--- | :--- |
| **Getting Started & CLI** | [Interactive CLI Workbench](cli_workbench/index.html) · [Setup, Architecture & Testing](docs/setup_and_architecture.md) · [End-to-End Usage Guide](docs/usage_guide.md) · [Running on Databricks](docs/databricks_guide.md) |
| **Specialised Workflows** | [Regulatory Exceedance & Threshold Risk Mapping](docs/regulatory_exceedance_mapping_guide.md) · [Optimal Monitoring Network Design & Gap Mapping](docs/optimal_monitoring_design_guide.md) · [Spatio-Temporal Trends & Borrowing](docs/spatiotemporal_trends_and_borrowing_guide.md) · [Multivariate Stream Co-Kriging](docs/multivariate_cokriging_guide.md) · [Multi-Source & Citizen-Science Fusion](docs/multi_source_and_citizen_science_guide.md) · [Non-Linear (`GAM + S3N`) & Tree-Ensemble (`RF/GBM`) Guide](docs/nonlinear_and_tree_ensemble_trend_guide.md) |
| **Customising the Pipeline** | [Custom Monitoring Sources](docs/custom_data_sources_guide.md) · [Custom Water Networks (`BYOWN`)](docs/custom_water_network_guide.md) · [ORN Covariates Reference (`--vars`)](docs/orn_variables_reference.md) |
| **Diagnostics & Interpretation** | [Station Snapping Guide](docs/station_snapping_guide.md) · [Descriptive Statistics Guide](docs/descriptive_stats_interpretation_guide.md) · [Model Comparison Guide](docs/model_comparison_interpretation_guide.md) |
| **Methodology & Theory** | [S3N Modelling & Kriging Guide](docs/s3n_modelling_and_kriging_guide.md) · [Two-Stage GAM & Measurement-Error Kriging](docs/two_stage_gam_measurement_error_kriging_guide.md) · [Preprocessing & C++ Engine](docs/preprocessing_pipeline.md) · [Where This Differs from `S3N`](docs/where_this_differs_from_s3n.md) |


## 🤖 AI Disclosure

This project was created with the assistance of artificial intelligence.

- **Coding Assistant:** Google DeepMind Antigravity
- **Models Used:** Gemini 4 Argon, Gemini 4 Hydrogen
- **Human Oversight:** All code, statistical methods, documentation, and generated outputs were reviewed and verified by a human maintainer. Users remain responsible for validating model assumptions, outputs, and suitability for any downstream scientific or regulatory use.


## 📄 License

This codebase is licensed under the **GNU General Public License v3.0 (`GPL-3.0`)** — see the [LICENSE](LICENSE) file for the full text, consistent with the `GPL-3` / `GPL-2 | GPL-3` licences of `S3N`, and `BRISC`.

External river-network geometries and monitoring datasets (`gis/` and `data/`) are **not** included in this repository; users download them directly from their respective providers and remain responsible for complying with those datasets' terms and attributions (such as the UK Open Government Licence v3.0 for OS Open Rivers / Environment Agency ARD data).
