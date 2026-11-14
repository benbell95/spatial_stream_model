# Scalable Spatial Stream Network (S3N) Pipeline — Usage & CLI Reference Guide

This guide documents all scripts, command-line options, output directory conventions, and practical examples for running the Scalable Spatial Stream Network (S3N) water quality pipeline across the Ordnance Survey Open Rivers Network ([`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg), 90,876 reaches across 141 monitored catchments), the Environment Agency River Surveillance Network ([`RSN`](../data/NCEA_River_Surveillance_Network_ARD_Open/)), the Small Streams Network ([`SSN`](../data/NCEA_Small_Streams_Network_ARD_Open/)), and the combined River Surveillance + Small Streams ([`RSS`](../data/processed/orn_rss/)) monitoring network.

> **Mathematical & Algorithmic Deep-Dive**: For a step-by-step theoretical and mathematical explanation of how the `S3N` package, stream-network covariance functions (`Tail-Up` & `Tail-Down`), Nearest-Neighbour Gaussian Process (`NNGP`) estimation, and Universal Stream Kriging work, see [`./s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md).
> 
> **Two-Stage GAM Modelling & Measurement Error Kriging**: For a dedicated guide on using first-stage Generalised Additive Models (GAMs) with S3N to account for heteroskedastic estimation uncertainty via Approach A Measurement Error Kriging, see [`./two_stage_gam_measurement_error_kriging_guide.md`](./two_stage_gam_measurement_error_kriging_guide.md).
>
> **Interactive CLI Workbench**: If you prefer a friendly visual builder instead of assembling flags by hand, open [`cli_workbench/index.html`](../cli_workbench/index.html) directly in your browser (works completely offline). It groups all workflow flags into clean dropdowns with `i` hover tooltips and switches cleanly between **Bash**, **Windows PowerShell**, and **Windows Command Prompt** (with an optional full `Rscript.exe` path & R-version picker).

> [!NOTE]
> ℹ️ **Run from your system terminal, not inside an interactive R or Python prompt:** All `Rscript` and `python3` commands in this guide are **command-line (CLI) scripts** meant to be run from your OS terminal (Linux/macOS Bash, Windows PowerShell / Command Prompt, or the RStudio **Terminal** tab) in the repository root — not pasted into an interactive `R>` console or `>>>` Python REPL.
>
> If you prefer to launch a run from inside an R/RStudio script or console, invoke it via `system2()`, for example:
> `system2(file.path(R.home("bin"), "Rscript"), c("scripts/run_s3n_determinand.R", "0117", "--transform=sqrt"))`
> *(On Windows, if `Rscript.exe` is not on `PATH`, either use its full path or generate ready-to-paste commands in [`cli_workbench/index.html`](../cli_workbench/index.html).)*

---

## Table of Contents

- [1. Overview of Pipeline Scripts](#1-overview-of-pipeline-scripts)
- [2. Primary Runner: `scripts/run_s3n_determinand.R`](#2-primary-runner-scriptsrun_s3n_determinandr)
  - [2.1 Synopsis](#21-synopsis)
  - [2.2 Complete Option Reference](#22-complete-option-reference)
  - [2.3 Common Determinand Codes (`det_code`)](#23-common-determinand-codes-det_code)
  - [2.4 Ecological & Biological Datasets (`--dataset=inverts|diatoms|macrophytes|rhs|lcms`)](#24-ecological--biological-datasets---datasetinvertsdiatomsmacrophytesrhslcms)
  - [2.5 Quick Decision Cheat Sheet ("Which Flag Should I Choose?")](#25-quick-decision-cheat-sheet-which-flag-should-i-choose)
  - [2.6 Universal Data Transformation Engine (`--transform=<type>`)](#26-universal-data-transformation-engine---transformtype)
  - [2.7 Multi-Network Monitoring Data Sources & Custom Datasets (`--source=<name>`)](#27-multi-network-monitoring-data-sources--custom-datasets---sourcename)
  - [2.8 Practical Examples for `scripts/run_s3n_determinand.R`](#28-practical-examples-for-scriptsrun_s3n_determinandr)
- [3. GIS & Network Preprocessor: `scripts/01_preprocess_orn.py`](#3-gis--network-preprocessor-scripts01_preprocess_ornpy)
  - [3.1 Synopsis](#31-synopsis)
  - [3.2 Complete Option Reference](#32-complete-option-reference)
  - [3.3 Station-to-Reach Snapping Tolerance (`--max-snap-dist`)](#33-station-to-reach-snapping-tolerance---max-snap-dist)
  - [3.4 Practical Examples for `scripts/01_preprocess_orn.py`](#34-practical-examples-for-scripts01_preprocess_ornpy)
- [4. Output Files & Directory Structure](#4-output-files--directory-structure)
  - [4.1 Preprocessed Data Directory (`data/processed/<net>_<source>/`)](#41-preprocessed-data-directory-dataprocessednet_source)
  - [4.2 Output Directory Layout (`output/`) & Model Run Results (`output/results/`)](#42-output-directory-layout-output--model-run-results-outputresults)
- [5. Software Dependencies, C++ Compilation & Automated Tests](#5-software-dependencies-c-compilation--automated-tests)
- [6. Single-Output Descriptive Statistics & Diagnostics: `scripts/05_descriptive_stats.R`](#6-single-output-descriptive-statistics--diagnostics-scripts05_descriptive_statsr)
- [7. Comparing 2 or More Model Runs & Cross-Determinand Relationships: `scripts/04_compare_model_runs.R`](#7-comparing-2-or-more-model-runs--cross-determinand-relationships-scripts04_compare_model_runsr)
- [8. Station-to-Reach Snapping Diagnostics & Catchment Plots: `scripts/06_station_snapping_diagnostics.R`](#8-station-to-reach-snapping-diagnostics--catchment-plots-scripts06_station_snapping_diagnosticsr)
- [9. Co-Kriging Comparisons, Transform Benchmarks & Dashboard Builder (`07`, `08`, `build_showcase.py`)](#9-co-kriging-comparisons-transform-benchmarks--dashboard-builder-07-08-build_showcasepy)
- [10. Regenerating Plots & Standalone Fitters (`scripts/03_generate_plots.R` & `scripts/02_fit_s3n_models.R`)](#10-regenerating-plots--standalone-fitters-scripts03_generate_plotsr--scripts02_fit_s3n_modelsr)

---

## 1. Overview of Pipeline Scripts

| Script | Language | Purpose |
| :--- | :--- | :--- |
| [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) | R + C++ | **Primary all-in-one runner**: Extracts any determinand (`det_code`), applies ROS/Tobit/0.5 × LOD censoring, filters or loops over time periods (`all`, `year`, `season`, `month`, `--st-borrow`), fits S3N models (`England`, `Humber`, `Thames`, `Severn`, plus `--basin`), computes normalised predictions, exceedance probabilities, optimal design sites, 5-fold CV, and renders figures (`fig1`–`fig7`). |
| [`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R) | R + C++ | **Multivariate Co-Kriging runner**: Jointly fits $P \ge 2$ water-quality or ecological targets across heterotopic stations via Scalable Factor-LMC Joint NNGP (`--targets=`, `--factors=`, `--method=lmc\|colocated`) and reports Dual CV (`All-Holdout` & `Aux-Aided`). See [`docs/multivariate_cokriging_guide.md`](./multivariate_cokriging_guide.md). |
| [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) | Python 3 | **Network & GIS preprocessor**: Parses [`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg) (or custom `.gpkg`/`.parquet`), snaps stations to reaches (`--max-snap-dist`), builds directed river tree topology (`parent_dn`, `Dist2Mth_km`, `log_AFV`), and writes canonical CSV tables + `network_metadata.json` under `data/processed/<net>_<source>/`. |
| [`scripts/format_custom_data.py`](../scripts/format_custom_data.py) | Python 3 | **Custom dataset formatter**: Standardises raw monitoring CSVs (coordinates/EPSG, interval-censored `[lo, hi]` Tobit modelling, schema mapping, and optional `configs/sources/<name>.json` generation). |
| [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) | R + C++ | **Standalone model fitter**: Fits S3N models directly on the active preprocessed `england_rsn_reaches.csv` and writes flat CSV/RDS outputs to `output/results/<net>_<source>/`. |
| [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) | R | **Standalone figure generator**: Renders `fig1`–`fig7` (and `--basin` or `--norm=<method>` maps) from saved CSV outputs in any model run directory (`--dir=output/results/...`). |
| [`scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R) | R | **Model run comparator**: Compares 2+ model output directories across networks, sources, transforms, covariates, or temporal slices, writes `comparison_report.md` + CSV summaries under `output/comparisons/`, and renders `fig_comp1`–`fig_comp5`. |
| [`scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R) | R | **Descriptive & stratified analyser**: Computes distributional metrics, representativeness tests, and stratified environmental breakdowns for any single run directory (`<run_dir>/descriptives/`). |
| [`scripts/06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R) | R | **Snapping diagnostics & plotter**: Generates station-to-reach snapping tables (`station_snapping_*.csv`), reports, and zoomed catchment plots under `output/station_snapping/<net>_<source>/`. |
| [`scripts/07_compare_cokriging_runs.R`](../scripts/07_compare_cokriging_runs.R) | R | **Co-kriging run comparator**: Compares two co-kriging runs (e.g., equal vs. source-weighted runs) with CV comparisons, difference maps, and vector SVGs. |
| [`scripts/08_compare_transforms.R`](../scripts/08_compare_transforms.R) | R | **Transform benchmark script**: Evaluates Nitrate (`0117`) across `none`, `log`, `sqrt`, and `boxcox` ($\lambda=0.50$) runs in raw and model scales (`output/comparisons/cs7_transforms_0117/`). |
| [`scripts/build_showcase.py`](../scripts/build_showcase.py) | Python 3 | **Dashboard bundle builder**: Converts rendered `output/` figures to WebP/minified SVG assets and compiles CSV/Markdown outputs into `modelling_dashboard/data.js`. |
| Shared R & Python Modules | R / Python | Modular core shared across scripts: [`s3n_config.R`](../scripts/s3n_config.R), [`s3n_covariates.R`](../scripts/s3n_covariates.R), [`s3n_transforms.R`](../scripts/s3n_transforms.R), [`s3n_model_functions.R`](../scripts/s3n_model_functions.R), [`s3n_multivariate_functions.R`](../scripts/s3n_multivariate_functions.R), [`s3n_stats_functions.R`](../scripts/s3n_stats_functions.R), [`s3n_plot_functions.R`](../scripts/s3n_plot_functions.R), and [`s3n_py_utils.py`](../scripts/s3n_py_utils.py). |

> [!TIP]
> Once [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) has been run once to build `data/processed/orn_rsn/` (`england_rsn_reaches.csv`, `rsn_sites.csv`, `rsn_catchments_lookup.csv`, `network_metadata.json`), you can do all subsequent modelling, temporal slicing, custom basin analysis, and plotting directly with [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) without re-running the GIS preprocessor.

---

## 2. Primary Runner: [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R)

### 2.1 Synopsis
```bash
Rscript scripts/run_s3n_determinand.R --list [--dataset=wq|inverts|diatoms|macrophytes|rhs|lcms]
Rscript scripts/run_s3n_determinand.R --list-basins
Rscript scripts/run_s3n_determinand.R <code_or_var> [OPTIONS]
```

### 2.2 Complete Option Reference

| Option / Argument | Default | Description |
| :--- | :--- | :--- |
| `<code_or_var>` | *(Required)* | Either a 4-digit EA RSN determinand code (e.g., `0117` for Nitrate, `0180` for Orthophosphate, `0061` for pH) **or** an ecological/habitat metric column name (e.g., `whpt_aspt`, `whpt_ntaxa`, `life_family_index`, `psi_family_score`, `rmni`, `hqa_score`). |
| `--source=<name>` | `rsn` | Selects the monitoring station network data source:<br>• `rsn`: **River Surveillance Network** ([`RSN`](../data/NCEA_River_Surveillance_Network_ARD_Open/)), the default nationwide surveillance network (1,382 monitoring stations across wider rivers; saves outputs to `output/results/det_<code>_<net>_rsn/`)<br>• `ssn`: **Small Streams Network** ([`SSN`](../data/NCEA_Small_Streams_Network_ARD_Open/)), dedicated to headwaters and smaller tributaries (313 monitoring stations; saves outputs to `output/results/det_<code>_<net>_ssn/`)<br>• `rss`: **Combined River Surveillance & Small Streams**, merging RSN + SSN into a single integrated dataset (1,695 stations, 1,384 snapped observed reaches on ORN v2; saves outputs to `output/results/det_<code>_<net>_rss/`)<br>• **Custom or registered source**: e.g. `--source=my_study` (from `configs/sources/<name>.json` or `data/custom/<name>/`)<br>• **Multi-source delimiter combination**: e.g. `--source=rsn+my_study` or `--source=rsn,ssn,my_study` (normalised internally with `+`, preserving hyphens inside source names). See [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md) and [`docs/multi_source_and_citizen_science_guide.md`](./multi_source_and_citizen_science_guide.md). |
| `--prefer-source=<s1[,s2]>` | *(None)* | Reach-level source priority order when fusing networks (e.g. `--prefer-source=rsn` or `--prefer-source=rsn,ssn,campaign`): on co-located reaches (`RivID`), excludes lower-priority sources while retaining them on unmonitored reaches. |
| `--source-weight=<s:w,...>` | `1` | Continuous reliability weights in `(0, 1]` per `data_source` (e.g. `rsn:1,waterblitz:0.25`; automatically normalised if values exceed `1`). Controls source-grouped reach aggregation and analytically down-weights secondary sources in NNGP/S3N fitting & CV. |
| `--source-se=<s:sd,...>` | *(None)* | Explicit per-source physical observation standard error in **raw units** (e.g. `rsn:0.02,waterblitz:0.5`). Automatically enables measurement-error propagation (`meas_uncertainty`) through the Delta Method and kriging SE maps. |
| `--range-mode=<mode>` | `auto` | Controls range-interval suffix matching for custom datasets (`auto`, `tobit`, `geometric`, `midpoint`, `lower`, `upper`, `none`). |
| `--custom-samples=<path>` | *(None)* | Direct path to a custom sample observations CSV file (e.g., containing `site_id`, `sample_datetime`, `meas_result`). Ingests custom measurements into the modelling pipeline. |
| `--custom-sites=<path>` | *(None)* | Direct path to a custom station points CSV file (e.g., containing `site_id`, `easting`, `northing`). Automatically snaps custom points to reaches. |
| `--custom-dir=<path>` | *(None)* | Path to a custom directory containing `sites.csv` and/or `samples.csv`. |
| `--custom-name=<name>` | `custom` | Custom dataset tag/name used in directory paths and plot titles (default: `custom`). |
| `--list-sources` | — | Prints all supported and configured monitoring network data sources (`rsn`, `ssn`, `rss`, `all`, custom profiles), abbreviations, full descriptive labels, and resolved data directories. |
| `--dataset=<name>` | Auto-detected (`wq` for 4-digit codes) | Selects which dataset to model:<br>• `wq`: Water Quality Samples (`Water Quality Samples.csv`)<br>• `inverts`: Macroinvertebrate Samples (`Macroinvertebrate Samples.csv`)<br>• `diatoms`: Diatom Samples (`Diatom Samples.csv`)<br>• `macrophytes`: Macrophyte Samples (`Macrophyte Samples.csv`)<br>• `rhs`: River Habitat Surveys (`River Habitat Surveys.csv`, RSN only)<br>• `lcms`: LCMS/GCMS Organic Pollutants (`LCMS GCMS Samples.csv`, RSN only)<br>*(Note: If omitted, the script automatically detects `inverts`, `diatoms`, `macrophytes`, or `rhs` when you pass a column name like `whpt_aspt` or `rmni`)*. |
| `--list`, `-h`, `--help` | — | Lists available variables, sample counts, descriptions, and units for the selected `--dataset` (defaults to top 25 water quality determinands; combine with `--dataset=inverts`, `--dataset=diatoms`, `--dataset=macrophytes`, `--dataset=rhs`, or `--dataset=lcms`). |
| `--list-basins` | — | Prints all 11 Environment Agency River Basin Districts (`RBD_ID`, `RBD_Name`) and the top 30 individual catchments (`CatchID`, `Main_River`, `RBD_Name`, total reaches, and monitored reaches) from `data/processed/<net>_<source>/rsn_catchments_lookup.csv`. |
| `--list-networks` | — | Lists all configured water network profiles (ships with `orn` in [`configs/networks/orn.json`](../configs/networks/orn.json)) and user-added JSON profiles in [`configs/networks/`](../configs/networks). |
| `--network=<name>` | `orn` | Selects the river network dataset (`orn`, or a custom GeoPackage `.gpkg`, GeoParquet `.parquet`/`.geoparquet`, or JSON profile under `configs/networks/`). Automatically preprocesses the network if not yet cached in `data/processed/<name>_<source>/`, writes results under `output/results/det_<code>_<net>_<source>/`, and displays the network name on all figures. See [`docs/custom_water_network_guide.md`](./custom_water_network_guide.md). |
| `--list-vars` (or `--list-covariates`) | — | Lists all 15 available fixed-effect regression covariates, their short CLI aliases, transformations, and default inclusion status. |
| `--vars=<spec>` (or `--covariates=<spec>`) | `default` | Selects which fixed-effect predictors enter the regression design matrix `X` (without needing to re-run [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py)):<br>• `default`: `Elevation_z`, `Slope_z`, `log10_US_Accum_z`, `Is_Chalk_Limestone`, `Is_Sandstone`, `Is_Peat` (or `Is_Aquifer` for single basins), `Easting_z`, `Northing_z` (saves to `output/results/det_<code>_<net>_<src>/`)<br>• **Explicit subset**: e.g. `--vars=elevation,slope,aquifer` (saves to `output/results/det_<code>_<net>_<src>_vars_elev_slope_aquifer/` and adds covariate names to figure subtitles)<br>• **Relative/mixed add & remove (`+var`/`-var`)**: e.g. `--vars=-easting,-northing` (`_vars_no_east_no_north/`), `--vars=+sinuosity,+strahler`, or `--vars=-chalk,-sand,+aquifer`<br>• `none` (or `intercept`): Intercept-only model (`X = 1`, `_vars_none/`)<br>• `all`: All candidate covariates (`_vars_all/`). |
| `--list-transforms` | — | Prints a comprehensive reference table of all available data transformation algorithms, mathematical formulas $z = f(y)$, valid domains, and scientific recommendations across [`run_s3n_determinand.R`](../scripts/run_s3n_determinand.R), [`02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R), and [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py). |
| `--transform=<type>` (or `--trans=`, `--transformation=`) | `log` (`wq`/`lcms`), `none` (ecology) | Selects the mathematical response transformation applied prior to spatial covariance estimation and universal stream kriging. Supported options: `log`, `sqrt`, `cbrt`, `log1p`, `log10`, `fourthrt`, `boxcox`, `yeojohnson`, `reciprocal`, `asin`, `logit`, `none` (or `raw`; `--raw` / `--untransformed` aliases). Non-default transforms tag the output folder automatically (e.g. `_trans_sqrt/`). |
| `--lambda=<val>` (or `--transform-lambda=`) | Auto (Profile MLE) | Power parameter $\lambda$ for Box-Cox and Yeo-Johnson transformations. If omitted, the optimal parameter $\hat{\lambda}$ is automatically estimated from observed samples via Profile Maximum Likelihood over $[-2, 2]$. |
| `--scale=<val>` (or `--transform-scale=`) | Auto (`1.0` for proportions, `100.0` for `%`) | Upper bound scale parameter for bounded angular (`asin`) and logistic (`logit`) transformations. |
| `--list-kernels` | — | Lists all supported stream-network covariance kernels and exits. |
| `--kernel=<type>` | `exponential` | Selects the spatial stream-network covariance function (`exponential`, `spherical`, `mariah`, `linear`, or `auto` to select best kernel by NNGP AIC). |
| `--combine=<mode>` | `joint` | How Tail-Up + Tail-Down components are combined: `joint` (single unified `TU + TD + Euclidean` C++ NNGP Cholesky solve) or `precision` (separate TU & TD fits combined via precision weights). |
| `--afv=<scheme>` | `length` | Additive Function Value confluence weighting for Tail-Up covariance: `length` (upstream channel length), `shreve` (link magnitude), or `hydro_yield` (elevation + aquifer yield proxy). |
| `--list-trends` | — | Lists all supported fixed-effect mean-trend engines (`--trend=`) and exits. |
| `--trend=<engine>` (or `--gam-s3n` for `gam`) | `linear` | Mean-trend regression engine prior to S3N residual kriging: `linear` (OLS), `spline` (natural cubic splines via `splines::ns`), `gam` (REML-penalised splines via `mgcv::gam`), `rf` (OOB Random Forest), or `gbm` (5-fold cross-fitted Gradient Boosted Trees). |
| `--spline-df=<N>` / `--trees=<N>` | `3` / `100` | Degrees of freedom per continuous predictor for `--trend=spline`, and ensemble tree count for `--trend=rf` / `--trend=gbm`. |
| `--reverse-color` (aliases `--reverse-colour`, `--rev-color`, `--rev-colour`, `--reverse-scale`) | `FALSE` (blue = low, red = high) | Reverses the map colour palette (`val_cols_rev`) so that **low values are warm red/orange and high values are cool blue/teal** — useful for ecological/habitat quality metrics where higher scores indicate better condition (e.g., `hqa_score`, `whpt_aspt`, `life_family_index`, `psi_family_score`). Also supported in [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) (alongside `--normal-color` / `--no-reverse-color`). |
| `--norm=<method>` (or `--normalise=<method>` / `--normalize=<method>`) | `none` | Plots **normalised stream predictions** (computed automatically on every run using England-wide reference statistics so basins share a common scale) and appends `_<method>` to map filenames (e.g., `fig1_england_stream_0117_map_zscore.png`). Supported options:<br>• `none`: Standard model/raw concentration units<br>• `zscore` (or `z`): Standardised Z-scores `(x - mean) / sd` on the model scale, rendered with a symmetric diverging blue–white–red palette centred at `0.0`<br>• `robust_zscore` (or `robust`, `mad`): Median/MAD standardised Z-scores resistant to extreme outliers, centred at `0.0`<br>• `minmax` (or `min-max`): Linearly rescaled to `[0, 1]`<br>• `percentile` (or `pct`, `rank`): Scale-free empirical percentile rank `[0%, 100%]`.<br>Also supported in [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) (recalculates normalised columns on the fly if missing). |
| `--vmin=<val>` / `--vmax=<val>` | Auto-shared across periods | Optional explicit minimum and maximum scale bar bounds (in raw measurement units, or in normalised units when `--norm` is active, e.g., `--vmin=-2.5 --vmax=2.5`). Note: When running seasonal, yearly, or monthly plots (`--by=season|year|month` or individual period slices), scale bars are **automatically synchronised** across all periods even without `--vmin`/`--vmax`. |
| `--censor=<method>` | `ros` | Method for imputing left-censored (`<`) observations below the limit of detection:<br>• `ros`: **Robust Lognormal Regression on Order Statistics** (Helsel & Cohn, 1988).<br>• `em` (or `tobit`): **Spatial Tobit / Iterative EM Censored Kriging** — conditions censored reaches on their spatial stream-network neighbours via truncated-normal E-steps and NNGP M-steps (recommended for 30–80% nondetects; tags output folder `_tobit`).<br>• `half_lod`: Simple substitution of `< L` with `0.5 * L`. |
| `--meas-error` (or `--gam`) | Auto-detected | Activates **Approach A: Measurement Error Kriging** (heteroskedastic observation error). Accounts for `meas_uncertainty` via Delta Method variance conversion and microscale/estimation nugget decomposition in C++. Auto-activated whenever `meas_method` is `gam`, `tobit_interval`, `tobit_raw_sd`, `geometric_midpoint`, or `source_se` (or via `--source-se`). |
| `--boot=<N>` | `50` | Number of parametric bootstrap iterations (`BRISC_bootstrap`) used to compute 95% confidence intervals on fixed-effect and covariance parameters. Use `--boot=15` for faster exploratory runs or `--boot=100` for final tables. |
| `--list-validation` / `--list-cv` | — | Lists all supported cross-validation (`--cv-mode=`) strategies and exits. |
| `--cv-mode=<mode>` | `random` | Cross-validation splitting strategy: `random` (5-fold), `catchment` (5-fold Leave-Catchment-Out by `CatchID`), `strahler` (train on mainstems `Strahler >= 3`, test on headwaters `Strahler <= 2`), or `spatial_block` (`spatial` alias; 5-fold k-means on `(mid_X, mid_Y)`). See [`docs/s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md). |
| `--list-thresholds` | — | Lists built-in UK regulatory exceedance thresholds (`--threshold=uk`) and directions, then exits. |
| `--threshold=<val\|uk>` | *(None)* | Raw-scale threshold (or `uk` for built-in indicative UK presets) for Gaussian exceedance probability `exceed_prob`, 5-tier `compliance_tier` classification (`exceedance_summary_<det>.csv`), and management classification maps. See [`docs/regulatory_exceedance_mapping_guide.md`](./regulatory_exceedance_mapping_guide.md). |
| `--threshold-dir=<above\|below>` | `above` (or preset dir) | Exceedance probability tail direction: `above` ($P(Y > T)$, e.g. Nitrate/Phosphorus) or `below` ($P(Y < T)$, e.g. Dissolved Oxygen/minimum pH). |
| `--design-sites=<K>` | `0` | Runs the OpenMP C++ greedy variance-reduction solver to rank the top `K` unmonitored reaches for new station placement (`optimal_monitoring_design_<det>.csv`). See [`docs/s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md). |
| `--omp=<N>` | `1` (recommend `4`) | Number of OpenMP threads used by [`src/s3n_fast_neighbors.so`](../src/s3n_fast_neighbors.so) and `S3N` C++ kernels (defaults to `1`; pass `--omp=4` or your core count for multi-threaded NNGP solves). |
| `--basin=<spec>` | *(None)* | Additionally models, cross-validates, and plots one or more custom river basins alongside `All England`, `Humber`, `Thames`, and `Severn`. Accepts:<br>• **Mainstem river name** (case-insensitive substring, e.g., `--basin=Eden`, `--basin="Great Ouse"`, `--basin=Tyne`, `--basin=Ribble`)<br>• **Numeric `CatchID`** (e.g., `--basin=3125`, `--basin=27`, `--basin=2828`)<br>• **River Basin District name or ID** (e.g., `--basin=Anglian`, `--basin="South West"`, `--basin=rbd:5`)<br>• **Multiple basins**: Comma-separated (`--basin=Eden,Tyne`) or repeated `--basin=` flags. |
| `--year=<YYYY>` | *(All years)* | Filters water quality samples to a single calendar year (e.g., `--year=2021`, `--year=2022`, `--year=2023`, `--year=2024`). Saves outputs in `output/results/det_<code>_<net>_<src>/year_<YYYY>/`. |
| `--season=<season>` | *(All seasons)* | Filters water quality samples to a single meteorological season:<br>• `winter` (Dec, Jan, Feb)<br>• `spring` (Mar, Apr, May)<br>• `summer` (Jun, Jul, Aug)<br>• `autumn` or `fall` (Sep, Oct, Nov)<br>Saves outputs in `output/results/det_<code>_<net>_<src>/season_<season>/`. |
| `--month=<1..12>` | *(All months)* | Filters water quality samples to a single calendar month (`1` to `12`). Saves outputs in `output/results/det_<code>_<net>_<src>/month_<MM>/`. Can be combined with `--year=<YYYY>` (e.g., `--year=2023 --month=7`). |
| `--by=<mode>` | *(None)* | **Automated Temporal Batch Mode**: Automatically loops over all periods in `<mode>`, fitting the full S3N suite and generating all figures in per-period subdirectories, plus combined multi-period summary CSVs (`cv_summary_by_<mode>_<code>.csv`, `parameters_by_<mode>_<code>.csv`, `spatiotemporal_trends_by_<mode>_<code>.csv`) and `fig6_spatiotemporal_trends_by_<mode>_<code>.png`:<br>• `--by=year`: Fits every calendar year (`year_2021/`, `year_2022/`, `year_2023/`, `year_2024/`)<br>• `--by=season`: Fits all 4 meteorological seasons (`season_winter/`, `season_spring/`, `season_summer/`, `season_autumn/`)<br>• `--by=month`: Fits all 12 calendar months (`month_01/` .. `month_12/`). |
| `--st-borrow` | `FALSE` | Enables **Dynamic Spatio-Temporal Borrowing** across periods when `--by=` is active (precision-weighted exponential smoothing using estimated lag-1 temporal correlation $\hat{\rho}_t$; tags output folder `_st`). See [`docs/spatiotemporal_trends_and_borrowing_guide.md`](./spatiotemporal_trends_and_borrowing_guide.md). |
| `--gis-dir=<path>` (or `--gis=<path>`) | `configs/directories.json` (`gis`) | Custom path (relative or absolute) to the directory containing river network GeoPackages (`ORN_v2.gpkg`, etc.). Overrides `configs/directories.json` / `S3N_GIS_DIR`. |
| `--data-dir=<path>` (or `--data=<path>`) | `configs/directories.json` (`data`) | Custom path (relative or absolute) to the root data directory containing monitoring source folders and `processed/`. Overrides `configs/directories.json` / `S3N_DATA_DIR`. |
| `--processed-dir=<path>` (or `--processed=<path>`) | `<data_dir>/processed` (or `S3N_PROCESSED_DIR`) | Custom path to the processed networks root folder (containing `<net>_<source>/`). Overrides `<data_dir>/processed` and `S3N_PROCESSED_DIR`. |
| `--output-dir=<path>` (or `--out-dir=<path>`) | `configs/directories.json` (`output`) | Custom path (relative or absolute) to the root output directory where `results/`, `comparisons/`, and `station_snapping/` folders are written. Overrides `configs/directories.json` / `S3N_OUTPUT_DIR`. |
| `--fig-format=<png\|pdf\|svg>` | `png` | Output figure file format across all plotting scripts (`png`, vector `pdf`, or vector `svg`). |
| `--dpi=<N>` | `200` (per figure) | PNG raster resolution in DPI (36–2400). Pixel dimensions scale proportionally so text and margins keep identical physical layout. |
| `--fig-size=<pct\|scale>` | `100%` (`1.0`) | Scale figure physical dimensions by percentage or multiplier (e.g. `--fig-size=150%`, `--fig-size=+25%`, `--fig-size=-20%`, or `--fig-scale=0.8`). |

---

### 2.3 Common Determinand Codes (`det_code`)

| `det_code` | Determinand Description | Units | Recommended Scale Flag |
| :--- | :--- | :--- | :--- |
| `0117` | Nitrate as N | `mg/l` | Default (`--transform=log`; or `--transform=sqrt`) |
| `0180` | Orthophosphate, reactive as P | `mg/l` | Default (`--transform=log`) |
| `0111` | Ammoniacal Nitrogen as N | `mg/l` | Default (`--transform=log`) |
| `0116` | Nitrogen, Total Oxidised as N (TON) | `mg/l` | Default (`--transform=log`) |
| `0118` | Nitrite as N | `mg/l` | Default (`--transform=log`) |
| `0301` | Carbon, Organic, Dissolved as C (DOC) | `mg/l` | Default (`--transform=log`) |
| `0077` | Conductivity at 25 C | `uS/cm` | Default (`--transform=log`) |
| `0162` | Alkalinity to pH 4.5 as CaCO3 | `mg/l` | Default (`--transform=log`) |
| `0253` | Magnesium, Dissolved | `mg/l` | Default (`--transform=log`) |
| `0241` | Calcium, Dissolved | `mg/l` | Default (`--transform=log`) |
| `6396` | Turbidity : In Situ | `FTU` | Default (`--transform=log`) |
| `0061` | pH : In Situ | `pH units` | `--transform=none` (untransformed) |
| `0076` | Temperature of Water | `CEL` | `--transform=none` (untransformed) |
| `0085` | Oxygen, Dissolved, % Saturation | `%` | `--transform=none` (untransformed) |
| `9924` | Oxygen, Dissolved as O2 : In Situ | `mg/l` | `--transform=none` (untransformed) |

### 2.4 Ecological & Biological Datasets (`--dataset=inverts|diatoms|macrophytes|rhs|lcms`)

| Dataset (`--dataset=...`) | Variable Name | Description | Units | Default Scale |
| :--- | :--- | :--- | :--- | :--- |
| `inverts` | `whpt_aspt` | WHPT Average Score Per Taxon (macroinvertebrate metric) | `ASPT` | Raw (`--transform=none`) |
| `inverts` | `whpt_ntaxa` | WHPT Taxon Richness (N-Taxa) | `taxa` | Raw (`--transform=none`) |
| `inverts` | `life_family_index` | LIFE Flow Velocity Stress Index (Family) | `LIFE` | Raw (`--transform=none`) |
| `inverts` | `psi_family_score` | PSI Fine Sediment Sensitivity Score | `%` | Raw (`--transform=none`) |
| `inverts` | `epsi_family_score` | Empirically Weighted PSI Sediment Score | `%` | Raw (`--transform=none`) |
| `inverts` | `cci` | Community Conservation Index | `CCI` | Raw (`--transform=none`) |
| `inverts` | `dehli` | DEHLI Drought Stress Index | `DEHLI` | Raw (`--transform=none`) |
| `inverts` | `armi_score` | Anglers' Riverfly Monitoring Score | `score` | Raw (`--transform=none`) |
| `diatoms` | `percent_organic_tolerant` | Organic-Tolerant Diatom Relative Abundance | `%` | Raw (`--transform=none`) |
| `diatoms` | `percent_motile` | Motile Diatom Guild Relative Abundance | `%` | Raw (`--transform=none`) |
| `macrophytes` | `rmni` | River Macrophyte Nutrient Index (RMNI) | `RMNI` | Raw (`--transform=none`) |
| `macrophytes` | `rn_a_taxa` | Macrophyte Aquatic Taxon Richness | `taxa` | Raw (`--transform=none`) |
| `rhs` | `hqa_score` | River Habitat Survey Habitat Quality (HQA) | `HQA` | Raw (`--transform=none`) |
| `rhs` | `hms_score` | River Habitat Survey Modification Score (HMS) | `HMS` | Raw (or `--transform=log`) |

### 2.5 Quick Decision Cheat Sheet ("Which Flag Should I Choose?")

| Need | Recommended Flag / Value | Why |
| :--- | :--- | :--- |
| **Strongly right-skewed nutrients / trace metals** (`0180`, `0111`, `0301`) | `--transform=log` | Handles orders-of-magnitude concentration ranges cleanly. |
| **Moderately skewed nutrients** (`0117` Nitrate) | `--transform=sqrt` | Reduces right-skew without stretching clean headwaters toward `-∞`. |
| **Symmetric ecological indices / temperature / pH** (`whpt_aspt`, `hqa_score`, `0061`) | `--transform=none` | Preserves native linear scale; no back-transformation needed. |
| **Data containing exact zeros** | `--transform=log1p` or `--transform=yeojohnson` | Avoids `log(0)` issues without arbitrary positive offsets. |
| **Normalised maps / CSVs across skewed basins** | `--norm=robust_zscore` (or `percentile`) | Median/MAD scaling resists extreme point-source hotspots; `percentile` gives `0–100%` rank maps. |
| **Standard cross-determinand comparison** | `--norm=zscore` (or `minmax`) | Mean/SD standardises to zero mean and unit variance (`Z`). |
| **Combining RSN + SSN monitoring networks** | `--source=rss` (or `--source=rsn+ssn`) | Includes both mainstem RSN and headwater SSN stations in one model. |
| **Left-censored lab data** (`<LOD` flags) | `--censor=ros` (default) or `--censor=em` | `ros` applies Lognormal ROS; `em` runs Spatial Tobit kriging across stream neighbours (`_tobit` folder). |
| **Range-binned custom data** (`[lo, hi]` intervals) | `--range-mode=tobit` (default `auto`) or `midpoint` | `tobit` models interval uncertainty via `IntervalTobitModel` (`meas_method="tobit_interval"`); `midpoint` uses `(lo+hi)/2`. |
| **Testing predictor sensitivity** | `--vars=-chalk,-sand,+aquifer` | Modifies the default 6-covariate set without listing every variable manually. |
| **Selecting stream covariance kernel** | `--kernel=auto` (or `exponential`, `spherical`, `mariah`, `linear`) | Fits compact-support (`spherical`, `linear`) or heavy-tailed (`mariah`) stream kernels, or selects the best kernel by NNGP AIC (`--list-kernels`). |
| **Joint vs. precision-weighted `TU + TD` combination** | `--combine=joint` (default) or `--combine=precision` | `joint` solves a single unified `Tail-Up + Tail-Down + Euclidean` NNGP Cholesky system; `precision` uses legacy spatial information-gain averaging. |
| **Physical confluence weighting (`AFV`)** | `--afv=hydro_yield` (or `length`, `shreve`) | Weights `Tail-Up` confluence splits by channel length (`length`), Shreve link magnitude (`shreve`), or elevation + aquifer hydrological yield (`hydro_yield`). |
| **Non-linear / spline fixed effects (`GAM + S3N`)** | `--trend=gam` (or `--gam-s3n`, `--trend=spline`) | Fits REML-penalised cubic regression splines (`mgcv::gam`) or GCV natural splines (`splines::ns`) before S3N stream kriging (`--list-trends`). |
| **Tree-ensemble regression kriging (`RF / GBM + S3N`)** | `--trend=rf` or `--trend=gbm` (`--trees=100`) | Fits Out-of-Bag Random Forest (`rf`) or 5-fold cross-fitted Gradient Boosted Trees (`gbm`) for high-order covariate interactions + S3N residual kriging. |

---

### 2.6 Universal Data Transformation Engine (`--transform=<type>`)

Geostatistical spatial stream network kriging assumes that the stochastic error process follows a multivariate Gaussian distribution. The pipeline provides 12 analytical transformations with exact inverse back-transformation (reach aggregation computes `obs_y = mean(fwd(y))` before analytical inverse mapping), automatic Profile Maximum Likelihood parameter estimation ($\hat{\lambda}$), and matching colour scales.

```bash
# Print the full transformation catalog from the command line
Rscript scripts/run_s3n_determinand.R --list-transforms
```

#### Supported Transformations Catalog

| Transform (`--transform=`) | Aliases | Forward Formula $z = f(y)$ | Inverse Formula $y = f^{-1}(z)$ | Domain | Scientific Application & Suitability |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `log` | `ln` | $\ln(y)$ | $\exp(z)$ | $y > 0$ | **Default for chemical & LCMS determinands**. Best for multiplicative concentrations spanning $2+$ orders of magnitude (e.g. Orthophosphate `0180`, Dissolved Organic Carbon `0301`, Conductivity `0077`). |
| `sqrt` | `square_root` | $\sqrt{y}$ | $z^2$ | $y \ge 0$ | **Recommended for moderately skewed concentrations and taxa counts**. Ideal for UK Nitrate `0117` (reduces skewness from `+1.00` to `+0.09`), Alkalinity `0162`, and `whpt_ntaxa`. |
| `cbrt` | `cube_root` | $\text{sign}(y) \cdot \lvert y \rvert^{1/3}$ | $z^3$ | All real $y$ | **Wilson-Hilferty transformation** for Gamma-distributed hydrological variables; preserves low values while compressing peaks. |
| `log1p` | `log_1p` | $\ln(1 + y)$ | $\exp(z) - 1$ | $y \ge 0$ | **Zero-tolerant log transform** for low taxa counts or trace occurrences with true zeros ($y = 0$). |
| `log10` | `log_10` | $\log_{10}(y)$ | $10^z$ | $y > 0$ | **Base-10 order-of-magnitude scaling** for multi-order variables (pathogens, conductivity). |
| `fourthrt` | `qtrt` | $y^{1/4}$ | $z^4$ | $y \ge 0$ | **Quartic root transform** for heavy-tailed abundance and biomass data. |
| `boxcox` | `bc`, `box_cox` | $(y^\lambda - 1)/\lambda$ | $(\lambda z + 1)^{1/\lambda}$ | $y > 0$ | **Profile MLE power optimisation** over $\hat{\lambda} \in [-2, 2]$. |
| `yeojohnson` | `yj` | $\text{Yeo-Johnson}(\lambda)$ | $\text{Inverse YJ}(z, \lambda)$ | All real $y$ | **Zero- and negative-tolerant Box-Cox extension** with Profile MLE $\hat{\lambda}$. |
| `reciprocal` | `inv` | $-1 / y$ | $-1 / z$ | $y > 0$ | **Rate-dilution transform** preserving monotonic spatial ordering. |
| `asin` | `arcsin` | $\arcsin(\sqrt{y / \text{scale}})$ | $\text{scale} \cdot \sin^2(z)$ | $0 \le y \le \text{scale}$ | **Angular transform** for bounded proportions and percentages (`percent_silt`, diatom abundances). |
| `logit` | `logistic` | $\ln(p / (1 - p))$ | $\text{scale} / (1 + \exp(-z))$ | $0 < y < \text{scale}$ | **Logistic transform** for bounded proportions (`0 < y < scale`). |
| `none` | `raw`, `identity` | $y$ | $z$ | All real $y$ | **Default for multimetric ecological indices** (`whpt_aspt`, `life_family_index`, `hqa_score`). |

> [!TIP]
> **Why `sqrt` improves UK Nitrate (`0117`) over `log`:**
> Raw UK Nitrate in [`england_rsn_reaches.csv`](../data/processed/orn_rsn/england_rsn_reaches.csv) has moderate right skewness (`+1.001`). A log transform over-compresses low concentrations (skewness `-1.128`), whereas `--transform=sqrt` produces near-Gaussian symmetry (`+0.092`).

---

### 2.7 Multi-Network Monitoring Data Sources & Custom Datasets (`--source=<name>`)

The pipeline supports official statutory networks, custom field campaigns, and arbitrary multi-network combinations via `--source=<name>`:

| Source Flag | Abbreviation | Monitoring Network | Focus / Character | Total Stations | Snapped Reaches (ORN v2) | Total Reaches | Output Folder Tag |
| :---: | :---: | :--- | :--- | ---: | ---: | ---: | :--- |
| `--source=rsn` *(default)* | `RSN` | **River Surveillance Network** | Nationwide network across wider rivers and mainstems | 1,382 | 1,210 | 90,876 | `output/results/det_<code>_<net>_rsn/` |
| `--source=ssn` | `SSN` | **Small Streams Network** | Dedicated headwater streams (< 20 km² drainage) | 313 | 177 | 76,242 | `output/results/det_<code>_<net>_ssn/` |
| `--source=rss` | `RSS` | **Combined Surveillance & Small Streams** | Integrated nationwide surveillance + headwater monitoring | 1,695 | 1,384 | 91,118 | `output/results/det_<code>_<net>_rss/` |
| `--source=<custom>` | User-defined | **Custom User / Field Dataset** | Field study or local monitoring campaign | Variable | Variable | Monitored catchments | `output/results/det_<code>_<net>_<custom>/` |
| `--source=<A>+<B>` | `<A>+<B>` | **Multi-Source Combination** | Merged union of multiple registered or custom sources | Variable | Variable | Union of catchments | `output/results/det_<code>_<net>_<A>+<B>/` |

> [!TIP]
> **Complete Custom Data Guide & Automated Formatter**: For instructions on bringing your own point observation datasets, creating `sites.csv` and `samples.csv`, registering profiles in `configs/sources/`, and handling lat/lon coordinates, see [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md).
> To automatically prepare an unstandardised raw CSV dataset (including interval-censored Tobit modelling for concentration ranges, analytical bin uncertainty, coordinate projection, and column remapping), run:
> ```bash
> python3 scripts/format_custom_data.py path/to/raw_data.csv --source-name=my_project --generate-config
> ```
> Output folders across formatting (`data/custom/`), network preprocessing (`data/processed/`), and modelling (`output/results/`) are tagged with the active range mode (e.g. `_tobit`, `_midpoint`).

> [!NOTE]
> **Dataset Availability Across Sources**:
> - **Water Quality (`--dataset=wq`)**, **Macroinvertebrates (`--dataset=inverts`)**, **Diatoms (`--dataset=diatoms`)**, and **Macrophytes (`--dataset=macrophytes`)** are monitored on **both RSN and SSN**. They can be modelled independently (`--source=rsn`, `--source=ssn`) or integrated directly (`--source=rss`).
> - **LCMS/GCMS Organic Pollutants (`--dataset=lcms`)** and **River Habitat Surveys (`--dataset=rhs`)** are monitored exclusively on RSN. If requested with `--source=ssn`, the pipeline halts with an informative message; if requested with `--source=rss`, it uses available RSN observations.

---

### 2.8 Practical Examples for [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R)

#### A. Discovery Commands
```bash
# List supported monitoring network sources and data paths
Rscript scripts/run_s3n_determinand.R --list-sources
python3 scripts/01_preprocess_orn.py --list-sources

# List the top 25 water quality determinands across monitoring sources
Rscript scripts/run_s3n_determinand.R --list                   # RSN (default)
Rscript scripts/run_s3n_determinand.R --list --source=ssn      # SSN
Rscript scripts/run_s3n_determinand.R --list --source=rss      # Combined RSS

# List available variables in the Macroinvertebrate, Diatom, Macrophyte, or RHS datasets
Rscript scripts/run_s3n_determinand.R --list --dataset=inverts
Rscript scripts/run_s3n_determinand.R --list --dataset=diatoms
Rscript scripts/run_s3n_determinand.R --list --dataset=macrophytes
Rscript scripts/run_s3n_determinand.R --list --dataset=rhs

# List all 11 River Basin Districts and top 30 individual catchments
Rscript scripts/run_s3n_determinand.R --list-basins
Rscript scripts/run_s3n_determinand.R --list-basins --source=ssn
Rscript scripts/run_s3n_determinand.R --list-basins --source=rss
```

#### B. Running Water Chemistry & Ecological / Macroinvertebrate Variables
```bash
# 1. Nitrate as N (0117) on RSN (default), SSN, and Combined RSS
Rscript scripts/run_s3n_determinand.R 0117                   # RSN (default nationwide network)
Rscript scripts/run_s3n_determinand.R 0117 --source=ssn      # Small Streams Network (SSN headwaters)
Rscript scripts/run_s3n_determinand.R 0117 --source=rss      # Combined RSN + SSN (RSS integrated network)

# 2. Macroinvertebrate WHPT ASPT across headwater streams (SSN)
Rscript scripts/run_s3n_determinand.R whpt_aspt --source=ssn --reverse-color

# 3. Macroinvertebrate WHPT ASPT across combined RSS surveillance network
Rscript scripts/run_s3n_determinand.R whpt_aspt --source=rss --reverse-color

# 4. Compare RSN vs. SSN vs. RSS models directly in scripts/04_compare_model_runs.R
Rscript scripts/04_compare_model_runs.R \
  output/results/det_0117_orn_rsn \
  output/results/det_0117_orn_ssn \
  output/results/det_0117_orn_rss \
  --out=output/comparisons/nitrate_rsn_vs_ssn_vs_rss

# 5. Macroinvertebrate WHPT ASPT, LIFE Flow Index, or PSI Sediment Score (auto-detects --dataset=inverts)
Rscript scripts/run_s3n_determinand.R whpt_aspt
Rscript scripts/run_s3n_determinand.R life_family_index
Rscript scripts/run_s3n_determinand.R psi_family_score

# 6. River Macrophyte Nutrient Index (RMNI) and River Habitat Survey Quality (HQA)
Rscript scripts/run_s3n_determinand.R rmni --dataset=macrophytes
Rscript scripts/run_s3n_determinand.R hqa_score --dataset=rhs

# 7. Orthophosphate as P (0180) comparing ROS vs. 0.5 * LOD substitution
Rscript scripts/run_s3n_determinand.R 0180 --censor=ros
Rscript scripts/run_s3n_determinand.R 0180 --censor=half_lod

# 8. In-situ pH (0061) modeled on its natural (raw) scale
Rscript scripts/run_s3n_determinand.R 0061 --transform=none
```

#### C. Two-Stage Modelling & Measurement Error Kriging (Approach A, GAM Estimates)
Once you have created a first-stage GAM prediction dataset in `data/custom/rsn_gam_mean/` (with `meas_uncertainty` and `meas_method = "gam"`):
```bash
# 1. Preprocess reach network using station-mean GAM predictions with uncertainty
python3 scripts/01_preprocess_orn.py --source=rsn_gam_mean

# 2. Fit S3N model using Approach A Measurement Error Kriging (auto-activated via meas_method="gam")
Rscript scripts/run_s3n_determinand.R 0117 --source=rsn_gam_mean --basin=Thames

# 3. Explicit CLI flag on any custom dataset with meas_uncertainty
Rscript scripts/run_s3n_determinand.R 0117 --source=my_gam_data --meas-error
```

#### C2. Non-Linear / Spline Fixed Effects (`GAM + S3N`) & Tree-Ensemble Regression Kriging (`--trend=...`)
```bash
# 1. List all available fixed-effect / mean-trend engines and tuning flags
Rscript scripts/run_s3n_determinand.R --list-trends

# 2. Penalised GAM Smooths + S3N (mgcv::gam REML cubic regression splines + S3N kriging)
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=gam

# 3. Natural Cubic Regression Splines + S3N with 4 degrees of freedom per continuous covariate
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=spline --spline-df=4

# 4. Random Forest Regression Kriging (Out-of-Bag RF trend + S3N stream kriging)
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=rf --trees=150

# 5. Gradient Boosted Regression Trees + S3N (Cross-fitted GBM trend + S3N stream kriging)
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=gbm --trees=100
```

> **Full Guide**: See [`./nonlinear_and_tree_ensemble_trend_guide.md`](./nonlinear_and_tree_ensemble_trend_guide.md) for details on using `--trend=` on standard raw data vs. `--meas-error`, package requirements (base/recommended R only), and interpreting `pred_trend_sc`, `GAM_EDF`, and `VarImp(%)`.

#### D. Modelling & Plotting Additional Custom Basins (`--basin=...`)
```bash
# 1. Add the River Eden catchment (CatchID=3125, Solway Tweed RBD) by river name
Rscript scripts/run_s3n_determinand.R 0117 --basin=Eden

# 2. Add the River Great Ouse catchment (CatchID=27, Anglian RBD) by quoted name or CatchID
Rscript scripts/run_s3n_determinand.R 0117 --basin="Great Ouse"
Rscript scripts/run_s3n_determinand.R 0117 --basin=27

# 3. Add multiple custom basins simultaneously (produces fig5 and fig6 maps + combined Fig 4)
Rscript scripts/run_s3n_determinand.R 0117 --basin=Eden,Tyne,Ribble

# 4. Model an entire Environment Agency River Basin District (e.g. Anglian RBD = 10,555 reaches)
Rscript scripts/run_s3n_determinand.R 0117 --basin=Anglian
Rscript scripts/run_s3n_determinand.R 0117 --basin="South West"
```

#### E. Temporal Filtering & Automated Multi-Period Batch Runs
```bash
# 1. Single calendar year (saved to output/results/det_0117_orn_rsn/year_2023/)
Rscript scripts/run_s3n_determinand.R 0117 --year=2023

# 2. Single meteorological season (saved to output/results/det_0117_orn_rsn/season_summer/)
Rscript scripts/run_s3n_determinand.R 0117 --season=summer

# 3. Single calendar month (saved to output/results/det_0117_orn_rsn/month_07/)
Rscript scripts/run_s3n_determinand.R 0117 --month=7

# 4. Specific year + season combination (saved to output/results/det_0117_orn_rsn/year_2023_season_winter/)
Rscript scripts/run_s3n_determinand.R 0117 --year=2023 --season=winter

# 5. Automated batch run across all 4 meteorological seasons (winter, spring, summer, autumn)
Rscript scripts/run_s3n_determinand.R 0117 --by=season

# 6. Automated batch run across all calendar years (2021, 2022, 2023, 2024) with custom basin
Rscript scripts/run_s3n_determinand.R 0117 --by=year --basin=Eden
```

#### F. Overriding Default Directories (`--gis-dir`, `--data-dir`, `--processed-dir`, `--output-dir`)

Default paths are resolved centrally from [`configs/directories.json`](../configs/directories.json) and `S3N_*` environment variables (see [`./setup_and_architecture.md`](./setup_and_architecture.md)). Any script accepts CLI directory overrides when targeting alternate locations:

```bash
# 1. Preprocess with custom GIS and data directories
python3 scripts/01_preprocess_orn.py --network=orn --gis-dir=custom_gis --data-dir=custom_data

# 2. Fit & plot with custom GIS, data, processed, and output directories
Rscript scripts/run_s3n_determinand.R 0117 \
  --gis-dir=custom_gis --data-dir=custom_data --processed-dir=custom_proc --output-dir=custom_output

# 3. Run single-output stats or multi-run comparisons against a custom output root
Rscript scripts/05_descriptive_stats.R --list --output-dir=custom_output
Rscript scripts/04_compare_model_runs.R det_0117_orn_rsn det_0180_orn_rsn --output-dir=custom_output
```

---

## 3. GIS & Network Preprocessor: [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py)

[`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) builds the stream network topology from [`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg) (or any custom network profile), snaps monitoring stations to their nearest reach polyline within a configurable distance tolerance (`--max-snap-dist`, default `250 m`), computes Shreve/accumulation-based additive function values (`log_AFV`), and exports standardised tables into `data/processed/<net>_<source>/` (default [`data/processed/orn_rsn/`](../data/processed/orn_rsn/)).

### 3.1 Synopsis
```bash
python3 scripts/01_preprocess_orn.py [OPTIONS]
```

### 3.2 Complete Option Reference

| Option / Argument | Default | Description |
| :--- | :--- | :--- |
| `--source <NAME>` | `rsn` | Monitoring network data source:<br>• `rsn`: River Surveillance Network (`RSN`, default nationwide surveillance dataset)<br>• `ssn`: Small Streams Network (`SSN`, dedicated headwater and tributary dataset)<br>• `rss`: Combined River Surveillance & Small Streams (`RSS`, integrated nationwide dataset)<br>• `<custom>`: Any configured profile in `configs/sources/<name>.json` or custom directory in `data/custom/<name>/`. Automatically preserves and propagates measurement uncertainty (`meas_uncertainty`) and estimation method (`meas_method`) to reach tables and metadata.<br>• Multi-source combination: e.g. `--source=rsn+custom` or `--source=rsn,ssn,custom` (normalised internally with `+`, preserving hyphens in source names). Preprocessed outputs are saved to `data/processed/<net>_<source>/` (e.g. `data/processed/orn_rsn/`). See [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md). |
| `--range-mode <MODE>` | `auto` | Controls range-interval strategy suffixing (`auto`, `tobit`, `geometric`, `midpoint`, `lower`, `upper`, `none`). |
| `--custom-sites <PATH>` | `None` | Direct path to a custom monitoring station points CSV file (e.g. with `site_id`, `easting`, `northing`). Automatically snaps custom points to reaches. |
| `--custom-samples <PATH>` | `None` | Direct path to a custom sample observations CSV file (e.g. with `site_id`, `sample_datetime`, `meas_result`). |
| `--custom-dir <DIR>` | `None` | Path to a custom directory containing `sites.csv` and/or `samples.csv`. |
| `--custom-name <NAME>` | `custom` | Short dataset tag/name used in directory paths and metadata (default: `custom`). |
| `--custom-crs <EPSG>` | *(Auto-detected)* | EPSG coordinate system code override for custom points (`27700` for British National Grid metres; `4326` for WGS84 decimal degrees, which are automatically projected to BNG metres). |
| `--list-sources` | — | Prints supported and configured monitoring network data sources, descriptions, and data paths, then exits. |
| `--max-snap-dist <METRES>` | `250.0` | **Maximum snapping distance tolerance in metres** from monitoring stations to river reach polylines. Stations with Euclidean snapping distance `d ≤ max_snap_dist` are retained as monitored reaches (`has_obs = 1`); stations farther than this threshold are marked `included_250m = 0` and excluded from reach observation targets to prevent spurious topological associations. |
| `--dataset <NAME>` | `wq` | Monitoring dataset to extract:<br>• `wq`: Water Quality Samples ([`RSN Water Quality Samples.csv`](../data/NCEA_River_Surveillance_Network_ARD_Open/RSN%20Water%20Quality%20Samples.csv))<br>• `inverts`: Macroinvertebrate Samples<br>• `diatoms`: Diatom Samples<br>• `macrophytes`: Macrophyte Samples<br>• `rhs`: River Habitat Surveys<br>• `lcms`: LCMS/GCMS Organic Pollutant Samples. |
| `--det-code <CODE>` | `0117` | 4-digit EA RSN determinand code (for `--dataset wq`, e.g. `0117` Nitrate, `0180` Orthophosphate) or column metric name (for ecological datasets, e.g. `whpt_aspt`, `rmni`, `hqa_score`). |
| `--transform <TYPE>` | `log` | Data transformation algorithm applied to determinand observations (`trans_obs_mean`, `mean_trans_obs`): `log`, `sqrt`, `cbrt`, `log1p`, `log10`, `fourthrt`, `boxcox`, `yeojohnson`, `reciprocal`, `asin`, `logit`, `none`/`raw`. |
| `--transform-lambda <VAL>` | Auto (Profile MLE) | Power parameter $\lambda$ for Box-Cox or Yeo-Johnson transformations. If omitted, $\hat{\lambda}$ is automatically estimated from observed station samples. |
| `--list-transforms` | — | Prints the full catalog of supported data transformations, formulas, and recommendations, then exits. |
| `--censor-method <METHOD>` | `ros` | Left-censored (`<`) value handling:<br>• `ros`: Pooled Robust Lognormal Regression on Order Statistics across all stations.<br>• `ros_site`: Within-station Lognormal ROS for stations with ≥ 3 uncensored detections and < 80% censoring (falling back to global ROS slope/intercept otherwise).<br>• `half_lod`: Substitution with 0.5 × LOD. |
| `--year <YYYY>` | `None` (all years) | Filters the primary reach table (`england_rsn_reaches.csv`) to samples from calendar year `<YYYY>`. Note: `rsn_reach_periods.csv` is always generated with all years, seasons, and months regardless of this filter. |
| `--season <SEASON>` | `None` (all seasons) | Filters the primary reach table (`england_rsn_reaches.csv`) to `winter`, `spring`, `summer`, or `autumn`/`fall`. |
| `--month <MONTH>` | `None` (all months) | Filters the primary reach table (`england_rsn_reaches.csv`) to a calendar month (`1`–`12` or `Jan`–`Dec`). |
| `--network <NAME\|PATH>` | `orn` | Water network profile name (`orn` or a user-added JSON config in `configs/networks/`) or direct path to a GeoPackage (`.gpkg`), GeoParquet (`.parquet`/`.geoparquet`/`.pq`), or `.json` network config. |
| `--network-config <PATH>` | `None` | Optional path to a custom network JSON configuration file. |
| `--topology-mode <MODE>` | `auto` | Override how river reach connectivity is constructed: `auto` (uses profile), `nodes` (`Fnode`/`Tnode` pairs), `pointer` (`NEXTDOWNID`), or `geometry` (coordinate snapping). |
| `--afv-col <COL>` | `None` | Override attribute column used for Additive Function Value (AFV) confluence weights (e.g. `UP_CELLS`, `CATCHMENT_AREA`). |
| `--afv-scheme <SCHEME>` | `length` | Physical confluence weighting scheme for `log_AFV`: `length` (cumulative upstream channel length `US_Accum_km`), `shreve` (Shreve link magnitude), or `hydro_yield` (elevation/orographic + Principal Aquifer baseflow hydrological yield proxy). |
| `--list-networks` | — | Lists all configured water network profiles and exits. |
| `--update-determinands` | `FALSE` | Scans all source monitoring files and regenerates [`configs/determinands.json`](../configs/determinands.json) with discovered codes, names, and units. |
| `--convert-to-geoparquet <PATH>` | `None` | Converts the selected `--network` GeoPackage layer(s) into an OGC GeoParquet (`.parquet`) file and exits. |
| `--gis-dir <DIR>` (or `--gis <DIR>`) | `gis` (or `$S3N_GIS_DIR`) | Custom directory containing input river network GeoPackage (`.gpkg`), GeoParquet (`.parquet`/`.geoparquet`), or `.zip` files. |
| `--data-dir <DIR>` (or `--data <DIR>`) | `data` (or `$S3N_DATA_DIR`) | Custom root data directory containing monitoring datasets and `processed/`. |
| `--out-dir <DIR>` (or `--output-dir <DIR>`) | `<data-dir>/processed/<net>_<source>` | Custom output directory for preprocessed network CSV files. |

---

### 3.3 Station-to-Reach Snapping Tolerance (`--max-snap-dist`)

Stations are snapped to river reach centrelines using a `2 km × 2 km` spatial hash grid and sub-metre orthogonal polyline projection (`included_250m = 1` when distance $d \le \text{max\_snap\_dist}$).

| Setting | Value | Recommended Use Cases |
| :--- | :---: | :--- |
| **Default Threshold** | `--max-snap-dist 250` (250 m) | **Standard for high-resolution networks (ORN v2, 1:25,000)** — retains 97.6% of RSN stations (median 8.0 m) while excluding terrestrial/borehole sites. |
| **Increased Threshold** | `--max-snap-dist 500` or `1000` | **Coarse networks (e.g. EU-Hydro 1:100k+) or wide floodplains/estuaries** where sites sit >250 m from the mapped centreline. |
| **Reduced Threshold** | `--max-snap-dist 50` or `100` | **Dense parallel drainage ditches/canals** (e.g. Fenlands, Somerset Levels) to avoid snapping across channel divides. |

To inspect exact snapping coordinates (`snap_dist_m`) and zoomed catchment maps, see Section 8 (`scripts/06_station_snapping_diagnostics.R`) and [`./station_snapping_guide.md`](./station_snapping_guide.md).

---

### 3.4 Practical Examples for [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py)

```bash
# 1. Default preprocessing (Nitrate 0117, ORN v2 network, 250 m snapping tolerance -> data/processed/orn_rsn/)
python3 scripts/01_preprocess_orn.py

# 2. Increase snapping tolerance to 500 m to capture edge-case stations on wide floodplains
python3 scripts/01_preprocess_orn.py --max-snap-dist 500

# 3. Restrict snapping tolerance to 100 m for dense canal / ditch networks
python3 scripts/01_preprocess_orn.py --max-snap-dist 100

# 4. Preprocess with Square Root transformation for Nitrate
python3 scripts/01_preprocess_orn.py --transform sqrt

# 5. Preprocess Orthophosphate (0180) for Summer 2023 with custom 300 m snapping
python3 scripts/01_preprocess_orn.py --det-code 0180 --year 2023 --season summer --max-snap-dist 300

# 6. Preprocess a custom GIS file or GeoPackage directly with 500 m snapping
python3 scripts/01_preprocess_orn.py --network gis/my_network.gpkg --max-snap-dist 500
```

---

## 4. Output Files & Directory Structure

### 4.1 Preprocessed Data Directory (`data/processed/<net>_<source>/`)
For the default `orn` + `rsn` run ([`data/processed/orn_rsn/`](../data/processed/orn_rsn/)):

| File | Rows | Description |
| :--- | ---: | :--- |
| [`england_rsn_reaches.csv`](../data/processed/orn_rsn/england_rsn_reaches.csv) | 90,876 | Master stream network table containing reach coordinates, topology (`parent_dn`, `CatchID`), distances (`Dist2Mth_km`, `up_dist_mid_km`, `obs_up_dist_km`), flow weights (`Shreve`, `log_AFV`), geology/elevation covariates, and observed water quality (`has_obs`, `obs_mean`, `obs_median`, `trans_obs_mean`, `meas_uncertainty`, `meas_method`, `n_samples`). |
| [`rsn_sites.csv`](../data/processed/orn_rsn/rsn_sites.csv) | 1,382 | Station-level lookup mapping each monitoring station (`site_id`) to its snapped ORN reach (`RivID`), snapping distance (`snap_dist_m`), and station summary statistics (`obs_mean`, `obs_median`, `obs_sd`, `trans_obs_mean`, `mean_trans_obs`). |
| [`rsn_catchments_lookup.csv`](../data/processed/orn_rsn/rsn_catchments_lookup.csv) | 141 | Lookup table of all 141 monitored catchments with `CatchID`, `Main_River`, `RBD_ID`, `RBD_Name`, `n_reaches`, and `n_obs_reaches`. |
| [`rsn_reach_periods.csv`](../data/processed/orn_rsn/rsn_reach_periods.csv) | 22,566 | Long-format reach summaries broken down by `period_type` (`all`, `year`, `season`, `month`, `year_season`) and `period_value`. |
| [`network_metadata.json`](../data/processed/orn_rsn/network_metadata.json) | — | Active water network & monitoring source metadata used automatically by R modelling and plotting scripts. |

### 4.2 Output Directory Layout (`output/`) & Model Run Results (`output/results/`)

All pipeline outputs are organised into dedicated subfolders under `output/`:
- `output/results/` — Model fitting outputs, predictions, CV metrics, parameter estimates, and maps from [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) (`output/results/det_<code>_<net>_<source>[_trans_<t>][_vars_<v>]/`) and [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) (`output/results/<net>_<source>[_trans_<t>][_vars_<v>]/`).
- `output/station_snapping/` — Station-to-reach snapping tables, summaries, and catchment maps from [`scripts/06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R) (`output/station_snapping/<net>_<source>/`).
- `output/comparisons/` — Cross-run comparison reports, statistical tests, and figures from [`scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R).

Each run of [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) creates a clean folder hierarchy under `output/results/det_<det_code>_<net>_<source>/` (or a temporal subfolder such as `output/results/det_0117_orn_rsn/season_summer/`):

| Output File | Description |
| :--- | :--- |
| `predictions_<code>.csv` | Predictions (`pred_tu`, `pred_td`, `pred_combined`), model-scale predictions/SEs (`pred_model_sc`, `se_model_sc`), and **normalised results** (`pred_zscore`, `pred_robust_zscore`, `pred_minmax`, `pred_percentile`, plus `_tu_`, `_td_`, and `obs_` equivalents computed from England-wide reference statistics) across all reaches in England. |
| `predictions_humber_<code>.csv`, `predictions_thames_<code>.csv`, `predictions_severn_<code>.csv` | Core basin reach predictions, normalised results, and kriging standard errors from each basin-specific S3N fit. |
| `predictions_basin_<slug>_<code>.csv` | *(When `--basin` is used)* Reach-level predictions, normalised results, and standard errors for each custom basin from its basin-specific S3N fit. |
| `cv_comparison_<code>.csv` | 5-fold cross-validation metrics (`RMSE`, `MAE`, `R2`, plus transform metadata `Transform_Type`, `Transform_Lambda`, `Transform_Label`, `Transform_Scale`) comparing OLS, S3N Tail-Down, S3N Tail-Up, and S3N Combined across `All England`, `Humber`, `Thames`, `Severn`, and any custom `--basin`. |
| `cv_predictions_england_<code>.csv`, `cv_predictions_humber_<code>.csv` | Held-out 5-fold cross-validation station/reach predictions (`obs_raw`, `obs_model_sc`, OLS/TU/TD/Combined CV predictions). |
| `parameters_<code>.csv` | Fixed-effect regression coefficients (`beta`) and spatial covariance parameters (`sigma.sq`, `tau.sq`, `phi`, and `tau.sq (microscale)`/`mean_meas_variance` under Measurement Error Kriging) with 95% parametric bootstrap confidence intervals. |
| `fig1_england_stream_<code>_map[_<norm>].png` | **Figure 1**: 2-panel national map across all reaches (Panel A: S3N Predicted Concentration or Normalised Score; Panel B: S3N Kriging SE + monitoring stations). |
| `fig2_humber_basin_<code>_map[_<norm>].png` | **Figure 2**: 2-panel high-resolution map of the Humber Basin (`CatchID=80`). |
| `fig3_thames_severn_<code>_maps[_<norm>].png` | **Figure 3**: Side-by-side S3N maps of the Thames Basin (`CatchID=3055`) and Severn Basin (`CatchID=2786`). |
| `fig4_cv_and_diagnostics_<code>.png` | **Figure 4**: 4-panel diagnostic suite — (A) Out-of-sample 5-fold CV scatter plot, (B) 5-fold CV R² barplot across all modelled basins, (C) 95% bootstrap CI forest plot of environmental drivers, and (D) Boxplot of predicted concentrations across all reaches (overall) vs. upland/lowland/aquifer regimes. |
| `fig5_<slug>_basin_<code>_map[_<norm>].png` | *(When `--basin` is used)* **Figure 5+**: Dedicated 2-panel high-resolution S3N prediction and Kriging uncertainty map for each user-specified custom basin (numbered `fig5`, `fig6`, ... when multiple basins are given). |
| `fig5_england_by_<mode>_<code>_map[_<norm>].png` | *(When `--by=year\|season\|month` is used)* **Temporal Grid Figure**: Multi-panel All-England prediction grid with shared scale bar saved in the run root alongside `cv_summary_by_<mode>_<code>.csv` and `parameters_by_<mode>_<code>.csv`. |
| `spatiotemporal_trends_by_<mode>_<code>.csv` & `fig6_spatiotemporal_trends_by_<mode>_<code>.png` | *(When `--by=` is used)* Reach-level Sen/Theil-Sen temporal slope estimates, borrowing $\hat{\rho}_t$, and multi-period trend maps/profiles. |
| `spatiotemporal_transitions_by_<mode>_<code>.csv` & `fig7_spatiotemporal_transitions_by_<mode>_<code>.png` | *(When `--by=` is used)* Consecutive period-to-period flux/hotspot transitions and change diagnostics (`fig7`). |
| `exceedance_summary_<code>.csv` & `fig_exceedance_*` | *(When `--threshold=` is used)* Gaussian exceedance probability `exceed_prob`, 5-tier `compliance_tier`, and regulatory compliance maps. |
| `optimal_monitoring_design_<code>.csv` & `fig_optimal_design_*` | *(When `--design-sites=<K>` is used)* Top `K` unmonitored reaches selected by the OpenMP C++ greedy variance-reduction solver and site placement maps. |

---

## 5. Software Dependencies, C++ Compilation & Automated Tests

The pipeline runs with **Python ≥ 3.8** (core preprocessing uses the Python standard library; optional `pyarrow`/`pyproj`/`Pillow` unlock GeoParquet, CRS reprojection, and `build_showcase.py` WebP asset bundling) and **R ≥ 3.5.0** (`S3N`, CRAN dependencies, optional `data.table` for fast CSV I/O, and [`src/s3n_fast_neighbors.so`](../src/s3n_fast_neighbors.so) / `.dll`).

- Full dependency tables, `configs/determinands.json` schema, and `R CMD SHLIB` build details: see [`./setup_and_architecture.md`](./setup_and_architecture.md).
- Run automated tests: `python3 -m unittest tests/test_python_utils.py` and `Rscript tests/test_r_and_cpp.R`.

---

## 6. Single-Output Descriptive Statistics & Diagnostics: [`scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R)

> **Detailed Interpretation Guide**: See [`./descriptive_stats_interpretation_guide.md`](./descriptive_stats_interpretation_guide.md) for definitions of all metrics, representativeness tests, stratified gradients, and figures (`fig_desc1_*`, `fig_desc2_*`).

Writes `descriptive_report.md`, `descriptive_stats_summary.csv`, `monitoring_representativeness.csv`, `stratified_results_summary.csv`, `within_run_statistical_tests.csv`, and `fig_desc1_*`–`fig_desc2_*` into `<run_dir>/descriptives/` (or `--out=<dir>`).

```bash
# List available model output directories under output/results/
Rscript scripts/05_descriptive_stats.R --list [--output-dir=<dir>]

# Run descriptive statistics on a single model output directory
Rscript scripts/05_descriptive_stats.R output/results/det_0117_orn_rsn
Rscript scripts/05_descriptive_stats.R output/results/det_0180_orn_rsn --label="Orthophosphate (ORN v2)" --out=output/custom_descriptives --fig-format=pdf --dpi=300
```

---

## 7. Comparing 2 or More Model Runs & Cross-Determinand Relationships: [`scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R)

> **Detailed Interpretation Guide**: See [`./model_comparison_interpretation_guide.md`](./model_comparison_interpretation_guide.md) for details on global vs. co-located CV, hypothesis tests, reach co-occurrence, Jaccard hotspots, `N:P` stoichiometry, Wald `Z`-tests, and `fig_comp1`–`fig_comp5`.

Compares two or more model output directories (same determinand across networks/sources/transforms/covariates or cross-determinand, or temporal subdirectories via `--by=season|year|month`) and writes `comparison_report.md`, CSV summaries, and `fig_comp1`–`fig_comp5` into `output/comparisons/` (or `--out=<dir>`).

```bash
# List available model runs
Rscript scripts/04_compare_model_runs.R --list [--output-dir=<dir>]

# Compare two or more runs with custom labels, output folder, and figure options
Rscript scripts/04_compare_model_runs.R \
  output/results/det_0117_orn_rsn \
  output/results/det_0117_orn_rss \
  --labels="RSN (Default),RSS (Combined)" \
  --out=output/comparisons/nitrate_rsn_vs_rss \
  --fig-format=png --dpi=300 --fig-size=100%

# Automatically compare all seasonal, yearly, or monthly subdirectories inside a parent run
Rscript scripts/04_compare_model_runs.R output/results/det_0117_orn_rsn --by=season
```

---

## 8. Station-to-Reach Snapping Diagnostics & Catchment Plots: [`scripts/06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R)

> **Detailed Snapping Guide**: See [`./station_snapping_guide.md`](./station_snapping_guide.md) for complete column definitions, catchment-plot layout, and diagnostics interpretation.

Writes `station_snapping_all.csv`, `catchment_snapping_summary.csv`, `station_snapping_report.md`, `fig_snap_overview_<network>.png`, and `fig_snap_catchment_<CatchID>_<name>.png` into `output/station_snapping/<net>_<source>/` (or `--out=<dir>`).

```bash
# List catchments in the network with station counts and snapping distances
Rscript scripts/06_station_snapping_diagnostics.R --list-catchments [--source=rss]

# Generate snapping diagnostics for specific catchments (or exclude >250 m sites with --only-included)
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Humber,Thames,Severn,Eden
Rscript scripts/06_station_snapping_diagnostics.R --only-included --catchment=Thames --fig-format=pdf --dpi=300
```

Can also be sourced directly in R (`res <- analyze_station_snapping(catchment = c("Thames", "Eden"), network = "orn")`).

---

## 9. Co-Kriging Comparisons, Transform Benchmarks & Dashboard Builder (`07`, `08`, `build_showcase.py`)

- **[`scripts/07_compare_cokriging_runs.R`](../scripts/07_compare_cokriging_runs.R)** — Compares two multivariate co-kriging runs (e.g., equal vs. source-weighted RSN+SSN+WaterBlitz runs), generating `cs6_cv_comparison.csv`, `cs6_cross_covariance_comparison.csv`, `cs6_reach_comparison_summary.csv`, `fig_cs6_diff_maps_<det>.png`, `fig_cs6_equal_vs_favoured_comparison.svg`, and per-run `fig_cokrig_cross_correlation.svg`:
  ```bash
  Rscript scripts/07_compare_cokriging_runs.R \
    --eq-dir=output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit \
    --wt-dir=output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit \
    --sites-csv=data/processed/orn_rsn+ssn+waterblitz_tobit/rsn_sites.csv \
    --out-dir=output/comparisons/cs6_equal_vs_favoured
  ```
- **[`scripts/08_compare_transforms.R`](../scripts/08_compare_transforms.R)** — Compares Nitrate (`0117`) runs across `none`, `log`, `sqrt`, and `boxcox` ($\lambda=0.50$) in `output/comparisons/cs7_transforms_0117/` (`cs7_transform_summary.csv`, `fig_cs7_transform_comparison.svg`, `fig_cs7_transform_diff_0117.png`):
  ```bash
  Rscript scripts/08_compare_transforms.R --dpi=300
  ```
- **[`scripts/build_showcase.py`](../scripts/build_showcase.py)** — Converts rendered `output/` PNGs to WebP and compiles CSV/Markdown outputs into `modelling_dashboard/data.js` for zero-server viewing in `modelling_dashboard/index.html` (requires `Pillow`):
  ```bash
  python3 scripts/build_showcase.py
  ```

---

## 10. Regenerating Plots & Standalone Fitters ([`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) & [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R))

Because every model run saves its complete reach predictions, cross-validation tables, and parameter estimates as CSVs (`predictions_*.csv`, `cv_comparison_*.csv`, `parameters_*.csv`), you **never need to re-fit S3N models** just to change figure formats, resolution, physical dimensions, colour palettes, normalisations, or basin zooms. [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) reads any existing results directory in seconds and re-renders `fig1`–`fig4` (plus `fig5+` custom basins, `fig_exceedance_*`, and temporal grid/trend figures if present).

### Key Plot-Customisation Flags

| Flag | Default | Purpose |
| :--- | :--- | :--- |
| `--dir=<path>` *(or `--out-dir=` / `--output-dir=`)* | *(Required)* | Existing results folder (e.g. `output/results/det_0117_orn_rsn`, or a temporal batch folder like `output/results/det_0117_orn_rsn_st`). |
| `--fig-format=<png\|pdf\|svg>` | `png` | Output file format: Cairo `png`, vector `cairo_pdf` (`pdf`), or vector `svg`. |
| `--dpi=<N>` | `200` | Raster resolution in DPI for PNGs (`36`–`2400`). Pixel dimensions scale proportionally so fonts and margins keep exact proportions. |
| `--fig-size=<pct\|scale>` | `100%` | Scale figure physical width & height (`150%`, `+25%`, `-20%`, or multiplier `1.25`). |
| `--norm=<mode>` | `none` | Re-plot maps using `zscore`, `robust_zscore`, `minmax`, or `percentile` normalisation (computed on the fly if not already in CSV). |
| `--reverse-color` / `--normal-color` | `FALSE` | Invert palette (`blue = high, red = low` for ecological scores) or restore standard (`blue = low, red = high`). |
| `--vmin=<val>` / `--vmax=<val>` | Auto | Lock colorbar bounds (in raw units, or normalised units when `--norm=` is active). |
| `--basin=<CatchID\|Name>` | *(None)* | Render an additional high-resolution basin map (`fig5+_<basin>_map.*`) from the saved national predictions. |

### Practical Examples

```bash
# 1. Re-export an existing run as publication-ready vector PDFs or SVGs
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --fig-format=pdf
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --fig-format=svg --fig-size=+20%

# 2. High-resolution 300 DPI PNGs at 125% size, or Z-score normalised maps with fixed bounds
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --dpi=300 --fig-size=125%
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --norm=zscore --vmin=-2.5 --vmax=2.5

# 3. Ecological score with reversed palette, or extra basin zoom without re-fitting
Rscript scripts/03_generate_plots.R --dir=output/results/inverts_whpt_aspt_orn_rsn --reverse-color --fig-format=pdf
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --basin=Eden,Ribble --dpi=300
```

> [!TIP]
> `--fig-format=`, `--dpi=`, and `--fig-size=` work identically in [`run_s3n_determinand.R`](../scripts/run_s3n_determinand.R), [`run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R), [`04_compare_model_runs.R`](../scripts/04_compare_model_runs.R), [`05_descriptive_stats.R`](../scripts/05_descriptive_stats.R), and [`06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R).

### Standalone Script Flag Summary ([`02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) & [`03_generate_plots.R`](../scripts/03_generate_plots.R))
- **[`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R)**: Accepts `--source=`, `--range-mode=`, `--network=`, `--vars=`/`--covariates=`, `--transform=`, `--lambda=`, `--scale=`, `--kernel=`, `--combine=`, `--afv=`, `--trend=` (`--gam-s3n`), `--spline-df=`, `--trees=`, `--cv-mode=`, `--threshold=`, `--threshold-dir=`, `--design-sites=`, `--boot=`/`--n-boot=`, `--omp=`, `--meas-error`/`--gam`, `--gis-dir=`, `--data-dir=`, `--processed-dir=`, `--output-dir=`, plus `--list-sources`, `--list-networks`, `--list-vars`, `--list-transforms`, `--list-kernels`, `--list-trends`, `--list-validation`/`--list-cv`, `--list-thresholds`. Saves flat results + `s3n_model_results.rds` to `output/results/<net>_<source>[_trans_<t>][_vars_<v>]/`.
- **[`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R)**: Accepts `--dir=` (`--out-dir=`/`--output-dir=`/`--out=`/`--output=`), `--code=`, `--name=`, `--unit=`, `--transform=`, `--lambda=`, `--scale=`, `--norm=`, `--vmin=`, `--vmax=`, `--reverse-color`, `--normal-color` / `--no-reverse-color`, `--basin=`/`--catchid=`, `--fig-format=`, `--dpi=`, `--fig-size=`, `--data-dir=`, `--processed-dir=`, `--gis-dir=`.
