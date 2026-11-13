# S3N Documentation Index

Use this index to find the right guide for your task.

---

## 1. Getting Started & CLI Reference

| Guide | Purpose |
| :--- | :--- |
| [**Setup, Architecture & Testing**](./setup_and_architecture.md) | Repository layout, directory & determinand configuration (`configs/directories.json`, `configs/determinands.json`), Python/R/C++ dependencies, compiling `src/s3n_fast_neighbors.so`, and running the test suite. |
| [**End-to-End Usage Guide**](./usage_guide.md) | Complete CLI flag reference, flag/method decision cheat sheet, workflow examples (`01_preprocess_orn.py`, `run_s3n_determinand.R`, `02`/`03` split execution), and output file descriptions. |
| [**Running on Databricks**](./databricks_guide.md) | Cell-by-cell notebook pattern for running the hybrid Python/R pipeline on Databricks Git Folders and Unity Catalog Volumes. |

---

## 2. Specialised Modelling Workflows

| Guide | Purpose |
| :--- | :--- |
| [**Regulatory Exceedance & WFD Risk Mapping**](./regulatory_exceedance_mapping_guide.md) | Reach-level Gaussian exceedance probabilities (`--threshold=`, `--threshold-dir=`), built-in UK/WFD presets (`--threshold=uk`), 5-tier compliance classification, and dedicated risk maps. |
| [**Spatio-Temporal Trends & Dynamic Borrowing**](./spatiotemporal_trends_and_borrowing_guide.md) | Multi-year/seasonal/monthly batch mapping (`--by=`), precision-weighted dynamic temporal borrowing (`--st-borrow`), and reach-level WLS trend significance (`fig6`). |
| [**Multivariate Stream Co-Kriging (Factor-LMC)**](./multivariate_cokriging_guide.md) | Jointly modelling $P \ge 2$ chemical and ecological targets across heterotopic networks via `scripts/run_s3n_cokriging.R`, Factor-LMC NNGP (`--factors`), and Dual CV benchmarks. |
| [**Non-Linear Spline (`GAM + S3N`) & Tree-Ensemble (`RF / GBM + S3N`) Guide**](./nonlinear_and_tree_ensemble_trend_guide.md) | Using `--trend=spline\|gam\|rf\|gbm` on standard raw or GAM-weighted data, zero-external-package architecture, orthogonal OOB/OOF residual kriging, and interpreting `GAM_EDF` / `VarImp(%)`. |

---

## 3. Customising Data, Networks & Covariates

| Guide | Purpose |
| :--- | :--- |
| [**Custom Monitoring Data Sources**](./custom_data_sources_guide.md) | Adding your own monitoring CSVs via `format_custom_data.py`, `configs/sources/*.json`, censored/range handling (`--range-mode`), and composite source combinations (`+`). |
| [**Custom Water Networks (`BYOWN`)**](./custom_water_network_guide.md) | Bringing your own river network GeoPackage/Shapefile via `configs/networks/<id>.json` and understanding the standardized 5-file processed schema. |
| [**ORN Variables & Covariates Reference**](./orn_variables_reference.md) | All `15` covariates in `COVARIATE_CATALOG`, default vs. optional predictors, and mixed `+/-` `--vars=` selection syntax. |

---

## 4. Diagnostics & Output Interpretation

| Guide | Purpose |
| :--- | :--- |
| [**Station Snapping Diagnostics**](./station_snapping_guide.md) | Running `06_station_snapping_diagnostics.R`, inspecting snapping distances/coordinates, and interpreting catchment snapping maps. |
| [**Descriptive Statistics Guide**](./descriptive_stats_interpretation_guide.md) | Interpreting single-run outputs from `05_descriptive_stats.R` (distributional metrics, monitoring bias tests, landscape/hydrogeology/stream-order gradients, within-run hypothesis tests). |
| [**Model Comparison Guide**](./model_comparison_interpretation_guide.md) | Interpreting multi-run outputs from `04_compare_model_runs.R` (cross-network, covariate, transform, temporal, and cross-determinand comparisons, `N:P` stoichiometry, `ΔZ` divergence, and Wald `Z`-tests). |

---

## 5. Methodology & Statistical Theory

| Guide | Purpose |
| :--- | :--- |
| [**S3N Modelling & Kriging Theory**](./s3n_modelling_and_kriging_guide.md) | Mathematical foundations of ROS imputation, 12 response transformations, Tail-Up / Tail-Down / Joint Additive stream covariance, multi-kernel NNGP estimation, hydrological CV modes (`--cv-mode=`), universal stream kriging, and C++ optimal monitoring design (`--design-sites=`). |
| [**Two-Stage GAM & Measurement-Error Kriging**](./two_stage_gam_measurement_error_kriging_guide.md) | Propagating stage-1 temporal/GAM standard errors (`meas_uncertainty`) into stage-2 S3N kriging (`delta_se` / `total_tau.sq`). |
| [**Preprocessing & C++ Neighbor Engine**](./preprocessing_pipeline.md) | Detailed walkthrough of `01_preprocess_orn.py` (GeoPackage extraction, ROS censored imputation, station snapping, AFV/Strahler computation) and the OpenMP C++ Euler-tour engine (`src/s3n_fast_neighbors.cpp`). |
| [**Where This Pipeline Differs from `S3N`**](./where_this_differs_from_s3n.md) | Attribution to `S3N` (`arXiv:2512.12398`) and `BRISC`, which `S3N` functions are still used, and what our C++/R/Python platform replaces or adds. |
