# S3N Documentation Index

Use this index to find the right guide for your task.

---

## Table of Contents

- [1. Getting Started & CLI Reference](#1-getting-started--cli-reference)
- [2. Multivariate & Multi-Source Co-Kriging](#2-multivariate--multi-source-co-kriging)
- [3. Specialised Modelling Workflows](#3-specialised-modelling-workflows)
- [4. Customising Data, Networks & Covariates](#4-customising-data-networks--covariates)
- [5. Diagnostics & Output Interpretation](#5-diagnostics--output-interpretation)
- [6. Methodology & Statistical Theory](#6-methodology--statistical-theory)

---

## 1. Getting Started & CLI Reference

| Guide | Purpose |
| :--- | :--- |
| [**CLI Workbench (`cli_workbench/index.html`)**](../cli_workbench/index.html) | Friendly standalone offline browser tool for building ready-to-run commands with grouped dropdowns, `i` option tooltips, and shell formatting for **Bash**, **Windows PowerShell**, and **Windows Command Prompt** (including optional full `Rscript.exe` path & R-version selection). |
| [**Setup, Architecture & Testing**](./setup_and_architecture.md) | Repository layout, obtaining & unzipping required datasets (`ORN` v2.0 GeoPackage, `RSN`, optional `SSN`), directory & determinand configuration (`configs/directories.json`, `configs/determinands.json`), Python/R/C++ dependencies, compiling `src/s3n_fast_neighbors.so` / `.dll`, and running `tests/`. |
| [**End-to-End Usage Guide**](./usage_guide.md) | Complete CLI flag reference across all scripts (`01`–`08`, `run_s3n_determinand.R`, `run_s3n_cokriging.R`, `format_custom_data.py`, `build_showcase.py`), decision cheat sheet, and output schemas. |
| [**Running on Databricks**](./databricks_guide.md) | Cell-by-cell notebook pattern for running the hybrid Python/R pipeline on Databricks Git Folders and Unity Catalog Volumes. |

---

## 2. Multivariate & Multi-Source Co-Kriging

| Guide | Purpose |
| :--- | :--- |
| [**Multivariate Stream Co-Kriging (Factor-LMC)**](./multivariate_cokriging_guide.md) | Jointly modelling $P \ge 2$ chemical and ecological targets across heterotopic networks via `scripts/run_s3n_cokriging.R`, Factor-LMC NNGP (`--factors`), and Dual CV benchmarks. |
| [**Multi-Source Fusion & Citizen-Science Co-Kriging**](./multi_source_and_citizen_science_guide.md) | Supplementing statutory networks (`RSN` / `SSN`) with citizen science or campaign data via Factor-LMC / Colocated co-kriging, reach priority (`--prefer-source`), reliability weighting (`--source-weight`), and source SE (`--source-se`). |

---

## 3. Specialised Modelling Workflows

| Guide | Purpose |
| :--- | :--- |
| [**Regulatory Exceedance & Threshold Risk Mapping**](./regulatory_exceedance_mapping_guide.md) | Reach-level Gaussian exceedance probabilities (`--threshold=`, `--threshold-dir=`), built-in indicative UK presets (`--threshold=uk`), 5-tier compliance classification, and dedicated risk maps. |
| [**Optimal Monitoring Network Design & Gap Mapping**](./optimal_monitoring_design_guide.md) | OpenMP C++ active-learning placement (`--design-sites=<K>`), Schur-complement variance reduction, ranked CSV exports, and `fig_optimal_design_*` monitoring gap & gain curve maps. |
| [**Spatio-Temporal Trends & Dynamic Borrowing**](./spatiotemporal_trends_and_borrowing_guide.md) | Multi-year/seasonal/monthly batch mapping (`--by=`), precision-weighted dynamic temporal borrowing (`--st-borrow`), and reach-level WLS trend significance (`fig6`). |
| [**Non-Linear Spline (`GAM + S3N`) & Tree-Ensemble (`RF / GBM + S3N`) Guide**](./nonlinear_and_tree_ensemble_trend_guide.md) | Using `--trend=spline\|gam\|rf\|gbm` on standard raw or GAM-weighted data, zero-external-package architecture, orthogonal OOB/OOF residual kriging, and interpreting `GAM_EDF` / `VarImp(%)`. |

---

## 4. Customising Data, Networks & Covariates

| Guide | Purpose |
| :--- | :--- |
| [**Custom Monitoring Data Sources**](./custom_data_sources_guide.md) | Adding your own monitoring CSVs via `format_custom_data.py`, `configs/sources/*.json`, censored/range handling (`--range-mode`), and composite source combinations (`+`). |
| [**Custom Water Networks (`BYOWN`)**](./custom_water_network_guide.md) | Bringing your own river network GeoPackage/Shapefile via `configs/networks/<id>.json` and understanding the standardised 5-file processed schema. |
| [**ORN Variables & Covariates Reference**](./orn_variables_reference.md) | All `15` covariates in `COVARIATE_CATALOG`, default vs. optional predictors, and mixed `+/-` `--vars=` selection syntax. |

---

## 5. Diagnostics & Output Interpretation

| Guide | Purpose |
| :--- | :--- |
| [**Station Snapping Diagnostics**](./station_snapping_guide.md) | Running `06_station_snapping_diagnostics.R`, inspecting snapping distances/coordinates, and interpreting catchment snapping maps. |
| [**Descriptive Statistics Guide**](./descriptive_stats_interpretation_guide.md) | Interpreting single-run outputs from `05_descriptive_stats.R` (distributional metrics, monitoring bias tests, landscape/hydrogeology/stream-order gradients, within-run hypothesis tests). |
| [**Model Comparison Guide**](./model_comparison_interpretation_guide.md) | Interpreting multi-run outputs from `04_compare_model_runs.R`, equal vs. weighted co-kriging comparisons (`07_compare_cokriging_runs.R`), transform benchmarks (`08_compare_transforms.R`), `N:P` stoichiometry, `ΔZ` divergence, and Wald `Z`-tests. |

---

## 6. Methodology & Statistical Theory

| Guide | Purpose |
| :--- | :--- |
| [**S3N Modelling & Kriging Theory**](./s3n_modelling_and_kriging_guide.md) | Mathematical foundations of ROS imputation, 12 response transformations, Tail-Up / Tail-Down / Joint Additive stream covariance, multi-kernel NNGP estimation, hydrological CV modes (`--cv-mode=`), universal stream kriging, and C++ optimal monitoring design (`--design-sites=`). |
| [**Where This Pipeline Differs from `S3N`**](./where_this_differs_from_s3n.md) | Attribution to `S3N` (`arXiv:2512.12398`) and `BRISC`, which `S3N` functions remain in use, and what the C++/R/Python platform replaces or adds. |
| [**Preprocessing & C++ Neighbour Engine**](./preprocessing_pipeline.md) | Walkthrough of `01_preprocess_orn.py` (GeoPackage extraction, ROS censored imputation, station snapping, AFV/Strahler computation) and the OpenMP C++ Euler-tour engine (`src/s3n_fast_neighbors.cpp`). |
| [**Two-Stage GAM & Measurement-Error Kriging**](./two_stage_gam_measurement_error_kriging_guide.md) | Propagating stage-1 temporal/GAM standard errors (`meas_uncertainty`) into stage-2 S3N kriging (`delta_se` / `total_tau.sq`). |
