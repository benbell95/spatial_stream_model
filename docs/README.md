# Documentation

Comprehensive guides covering setup, modelling workflows, custom data and river networks, diagnostic interpretation, and statistical theory.

- **Back to Project Home:** [Repository README](../README.md) · [Live Modelling Dashboard](https://benbell95.github.io/spatial_stream_model/) · [Offline Dashboard](../modelling_dashboard/index.html)

> [!TIP]
> **Suggested Reading Paths**
> - **First-time users:** [Setup, Architecture & Testing](./setup_and_architecture.md) &rarr; [End-to-End Usage Guide](./usage_guide.md)
> - **Bringing your own data:** [Custom Monitoring Data Sources](./custom_data_sources_guide.md) &rarr; [Custom Water Networks](./custom_water_network_guide.md)
> - **Statisticians & researchers:** [Spatial Stream Modelling & Kriging Theory](./s3n_modelling_and_kriging_guide.md) &rarr; [Where This Pipeline Differs from `S3N`](./where_this_differs_from_s3n.md)

## Quick Reference: Common Tasks

| I want to… | Read this |
| :--- | :--- |
| Install dependencies, download EA/ORN data, and run tests | [Setup, Architecture & Testing](./setup_and_architecture.md) |
| Map probabilities of exceeding water-quality thresholds | [Regulatory Exceedance & Threshold Risk Mapping](./regulatory_exceedance_mapping_guide.md) |
| Blend citizen-science (`WaterBlitz`) or external data with `RSN` | [Multi-Source Fusion & Citizen-Science Co-Kriging](./multi_source_and_citizen_science_guide.md) |
| Use my own monitoring CSVs or non-UK river network | [Custom Monitoring Data Sources](./custom_data_sources_guide.md) · [Custom Water Networks](./custom_water_network_guide.md) |
| Choose the best locations for new monitoring gauges | [Optimal Monitoring Network Design & Gap Mapping](./optimal_monitoring_design_guide.md) |

## Table of Contents

- [1. Getting Started](#1-getting-started)
- [2. Multivariate & Multi-Source Co-Kriging](#2-multivariate--multi-source-co-kriging)
- [3. Specialised Modelling Workflows](#3-specialised-modelling-workflows)
- [4. Customising Data, Networks & Covariates](#4-customising-data-networks--covariates)
- [5. Diagnostics & Output Interpretation](#5-diagnostics--output-interpretation)
- [6. Methodology & Statistical Theory](#6-methodology--statistical-theory)

---

## 1. Getting Started

| Guide | Purpose |
| :--- | :--- |
| [**Setup, Architecture & Testing**](./setup_and_architecture.md) | Repository layout, obtaining required datasets (`ORN` v2.0, `RSN`, optional `SSN`), directory & determinand configuration, Python/R/C++ dependencies, compiling the OpenMP C++ neighbour library, and running the test suite. |
| [**End-to-End Usage Guide**](./usage_guide.md) | Complete CLI flag reference across preprocessing, modelling, diagnostics, and reporting workflows, decision cheat sheet, and output schemas. |
| [**Running on Databricks**](./databricks_guide.md) | Cell-by-cell notebook pattern for running the hybrid Python/R pipeline on Databricks Git Folders and Unity Catalog Volumes. |

---

## 2. Multivariate & Multi-Source Co-Kriging

| Guide | Purpose |
| :--- | :--- |
| [**Multivariate Stream Co-Kriging (Factor-LMC)**](./multivariate_cokriging_guide.md) | Jointly modelling multiple correlated chemical and ecological targets across partially overlapping monitoring networks using Factor-LMC NNGP and Dual CV benchmarks. |
| [**Multi-Source Fusion & Citizen-Science Co-Kriging**](./multi_source_and_citizen_science_guide.md) | Supplementing EA surveillance networks (`RSN` / `SSN`) with citizen science or campaign data using reach-priority rules, reliability weighting, and source-specific measurement error. |

---

## 3. Specialised Modelling Workflows

| Guide | Purpose |
| :--- | :--- |
| [**Regulatory Exceedance & Threshold Risk Mapping**](./regulatory_exceedance_mapping_guide.md) | Reach-level Gaussian exceedance probabilities, built-in indicative UK thresholds, 5-tier compliance classification, and dedicated risk maps. |
| [**Optimal Monitoring Network Design & Gap Mapping**](./optimal_monitoring_design_guide.md) | Active-learning gauge placement via Schur-complement variance reduction, ranked candidate reach exports, and monitoring gap & gain curve maps. |
| [**Spatio-Temporal Trends & Dynamic Borrowing**](./spatiotemporal_trends_and_borrowing_guide.md) | Multi-year, seasonal, and monthly batch mapping, precision-weighted dynamic temporal borrowing, and reach-level trend significance maps. |
| [**Non-Linear Spline (`GAM + S3N`) & Tree-Ensemble (`RF / GBM + S3N`) Guide**](./nonlinear_and_tree_ensemble_trend_guide.md) | Fitting non-linear spline and tree-ensemble catchment trends without heavy external packages, orthogonal residual kriging, and interpreting effective degrees of freedom and variable importance. |

---

## 4. Customising Data, Networks & Covariates

| Guide | Purpose |
| :--- | :--- |
| [**Custom Monitoring Data Sources**](./custom_data_sources_guide.md) | Formatting and registering your own monitoring CSVs, censored/range-bracket handling, and combining multiple data sources. |
| [**Custom Water Networks**](./custom_water_network_guide.md) | Bringing your own river network GeoPackage or Shapefile and understanding the standardised 5-file processed schema. |
| [**ORN Variables & Covariates Reference**](./orn_variables_reference.md) | All `15` catchment predictors in the covariate catalogue, default vs. optional predictors, and how to include or exclude covariates. |

---

## 5. Diagnostics & Output Interpretation

| Guide | Purpose |
| :--- | :--- |
| [**Station Snapping Diagnostics**](./station_snapping_guide.md) | Auditing station-to-reach snapping, inspecting snapping distances and coordinates, and interpreting catchment snapping maps. |
| [**Descriptive Statistics Guide**](./descriptive_stats_interpretation_guide.md) | Interpreting single-run descriptive statistics, including distributional metrics, monitoring bias tests, landscape/stream-order gradients, and within-run hypothesis tests. |
| [**Model Comparison Guide**](./model_comparison_interpretation_guide.md) | Interpreting multi-run comparisons, equal vs. weighted co-kriging benchmarks, transform comparisons, `N:P` stoichiometry, `ΔZ` divergence, and Wald `Z`-tests. |

---

## 6. Methodology & Statistical Theory

| Guide | Purpose |
| :--- | :--- |
| [**Spatial Stream Modelling & Kriging Theory**](./s3n_modelling_and_kriging_guide.md) | Mathematical foundations of ROS imputation, 12 response transformations, Tail-Up / Tail-Down / Joint Additive stream covariance, multi-kernel NNGP estimation, hydrological cross-validation modes, universal stream kriging, and optimal monitoring design. |
| [**Where This Pipeline Differs from `S3N`**](./where_this_differs_from_s3n.md) | Attribution to `S3N` (`arXiv:2512.12398`) and `BRISC`, which `S3N` functions remain in use, and what this C++/R/Python platform replaces or adds. |
| [**Preprocessing & C++ Neighbour Engine**](./preprocessing_pipeline.md) | Walkthrough of river-network GeoPackage extraction, ROS censored imputation, station snapping, AFV/Strahler computation, and the OpenMP C++ Euler-tour neighbour engine. |
| [**Two-Stage GAM & Measurement-Error Kriging**](./two_stage_gam_measurement_error_kriging_guide.md) | Propagating stage-1 temporal/GAM standard errors into stage-2 spatial kriging variance. |
