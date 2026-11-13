# Where This Pipeline Differs from the S3N Package

## Acknowledgements, Foundations & Attribution

This pipeline is built upon and inspired by the **Scalable Spatial Stream Network (`S3N`)** statistical framework and its reference R implementation, [`jpierkunke/S3N`](https://github.com/jpierkunke/S3N). The key conceptual breakthrough—embedding dendritic river-network covariance structures inside a **Nearest-Neighbor Gaussian Process (`NNGP`)** so that likelihood evaluation and kriging scale linearly $\mathcal{O}(n)$ rather than cubically $\mathcal{O}(n^3)$—originates directly with `S3N`.

Full credit belongs to the authors of the foundational packages and theory on which this repository stands:

- **`S3N` R Package** ([https://github.com/jpierkunke/S3N](https://github.com/jpierkunke/S3N)) — **Jessica P. Kunke, Julian D. Olden, and Tyler H. McCormick**:
  > Kunke, J. P., Olden, J. D., & McCormick, T. H. (2025). *Scalable Spatial Stream Network (S3N) Models*. **arXiv:2512.12398 [stat.ME]**. [https://arxiv.org/abs/2512.12398](https://arxiv.org/abs/2512.12398)
- **`BRISC` R Package** (*Bootstrap for Rapid Inference on Spatial Covariances*) — **Arkajyoti Saha and Abhirup Datta**, whose sparse NNGP Cholesky conditioning, profile likelihood solver, and rapid parametric bootstrap algorithm provide the computational core extended by `S3N`:
  > Saha, A., & Datta, A. (2018). *BRISC: Bootstrap for Rapid Inference on Spatial Covariances*. *Stat*, 7(1), e184.
- **Foundational NNGP & Stream Geostatistics**:
  - Datta, A., Banerjee, S., Finley, A. O., & Gelfand, A. E. (2016). *Hierarchical Nearest-Neighbor Gaussian Process Models for Large Geostatistical Datasets*. *Journal of the American Statistical Association*, 111(514), 800–812.
  - Ver Hoef, J. M., Peterson, E. E., & Theobald, D. (2006). *Spatial Statistical Models that Use Flow and Stream Distance*. *Environmental and Ecological Statistics*, 13(4), 449–464.
  - Ver Hoef, J. M., & Peterson, E. E. (2010). *A Moving Average Approach for Spatial Statistical Models of Stream Networks*. *Journal of the American Statistical Association*, 105(489), 6–18.

In this repository, `S3N` serves as an optional fitting/kriging backend for standard univariate exponential models, while the surrounding codebase wraps and extends it into an end-to-end national preprocessing, modelling, validation, multivariate co-kriging, risk-mapping, and reporting platform.

---

## 1. High-Level Overview: `S3N` vs. This Pipeline

| Dimension | Upstream `S3N` R Package | This Pipeline (`spatial_stream_model`) |
| :--- | :--- | :--- |
| **Scope** | Univariate spatial estimation, bootstrap, and kriging functions in R (`S3N`). | End-to-end headless CLI pipeline (Python + R + OpenMP C++) from raw GIS/CSVs to maps and reports. |
| **GIS & Topology Prep** | Expects pre-processed, pre-snapped stream network topology/order tables prepared externally (e.g. via GIS/`SSN2`). | Zero-GIS Python extractor (`01_preprocess_orn.py`) reading `.gpkg`/`.parquet`, breaking braids/canal loops, and orthogonally snapping stations. |
| **Neighbor Construction** | R-level `order_stream()` and `get_neighbor_stream()` tree walks. | Sub-second OpenMP C++ Euler-tour LCA & `std::partial_sort` neighbor engine (`src/s3n_fast_neighbors.cpp`). |
| **Covariance Structure** | Single-component `Tail-Up` or `Tail-Down` exponential NNGP at a time. | 4 stream kernels (`exponential`, `spherical`, `mariah`, `linear`, plus AIC `auto`), plus **Joint Additive `Tail-Up + Tail-Down + Euclidean` NNGP**. |
| **Multivariate Targets** | Not supported (univariate only). | Dedicated **Factor-LMC Joint NNGP Co-Kriging** (`scripts/run_s3n_cokriging.R`) for $P \ge 2$ heterotopic variables. |
| **Uncertainty & Censoring** | Assumes uncensored, homoskedastic responses. | Multi-threshold Lognormal **ROS**, **Spatial Tobit EM** kriging, and **Two-Stage GAM Measurement-Error Kriging**. |
| **Mean Trend (`Xβ`)** | Parametric linear fixed effects. | Linear GLS plus orthogonalized **GCV Splines**, **REML GAMs**, **OOB Random Forests**, and **OOF Gradient Boosted Trees**. |
| **Decision Support** | Raw point predictions & prediction intervals. | Temporal batching + dynamic borrowing (`--by=`, `--st-borrow`), WFD exceedance probabilities (`--threshold=uk`), and C++ optimal station design (`--design-sites=`). |

---

## 2. Where This Pipeline Uses the `S3N` Package

Inside [`scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R), `load_s3n_backend()` checks whether the `S3N` R package is installed (`.s3n_env$has_s3n_pkg`). When **all three conditions hold**:
1. `S3N` is installed,
2. `--kernel=exponential` is selected, and
3. Heteroskedastic measurement error (`--meas-error` / `--gam`) is **not** active,

the pipeline delegates individual `Tail-Up` and `Tail-Down` component fitting and prediction to three `S3N` functions:

| `S3N` Function | Role in This Pipeline | Where Called |
| :--- | :--- | :--- |
| `BRISC_estimation_stream()` | Estimates GLS fixed-effect slopes `Beta` and spatial parameters `Theta = (sigma.sq, tau.sq, phi)` via sparse NNGP profile MLE (`L-BFGS-B`) on `nb_tu$obs_neighbors` and `nb_td$obs_neighbors`. | [`scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R) (`estimate_component()`) |
| `BRISC_bootstrap()` | Generates 95% parametric bootstrap confidence intervals for `Beta` and `Theta` when `engine == "brisc"`. | [`scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R) (`bootstrap_component()`) |
| `BRISC_prediction_stream()` | Computes universal stream kriging point predictions (`$prediction`) and 95% prediction intervals (`$prediction.ci`) across all `q_tot` reaches for `Tail-Up` and `Tail-Down`. | [`scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R) (`predict_component()`) |

> [!NOTE]
> **Automatic C++ Fallback**: If `S3N` is not installed in the R environment—or whenever non-exponential kernels, measurement-error kriging, Joint Additive mode, or Multivariate Co-Kriging are used—the pipeline runs cleanly on our custom OpenMP C++ NNGP solvers (`fit_nngp_cpp()`, `bootstrap_nngp_cpp()`, `C_s3n_predict_meas_error()`, `C_s3n_joint_nngp()`, `C_s3n_multivariate_nngp()`) in [`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp).

---

## 3. Where This Pipeline Replaces `S3N`

### 3.1 Zero-Dependency GIS Ingestion, Cycle Resolution & Station Snapping
Upstream `S3N` operates on already-cleaned stream topology and distance tables. [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) eliminates heavy GIS dependencies (`GDAL`, `geopandas`, `sf`, `SSN2`) by parsing multi-gigabyte OGC GeoPackages (`.gpkg`) and GeoParquet networks using Python's standard library (`sqlite3` + binary WKB decoding), enforcing strict directed rooted trees (resolving braided channels and canal loops via shortest-path-to-mouth and BFS cycle breaking), computing Strahler/Shreve orders, and snapping monitoring stations orthogonally (`<= 250 m`) to exact along-reach river mouth distances (`obs_up_dist_km`).

### 3.2 Sub-Second C++/OpenMP Euler-Tour Neighbor Engine
Instead of calling `S3N`'s R-level `order_stream()` and `get_neighbor_stream()`, [`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp) builds the exact `obs_neighbors`, `pred_neighbors`, `obsobs_dist`, and `obsobs_wt` structures in compiled C++:
- **`C_build_tree_tin_tout()`**: Runs a single $\mathcal{O}(q_{\text{tot}})$ DFS Euler tour storing entry/exit timestamps (`tin`, `tout`), reducing every flow-connectedness test `FC(a, b)` to an $\mathcal{O}(1)$ integer interval check (`tin[a] <= tin[b] && tout[b] <= tout[a]`).
- **`C_compute_s3n_neighbors()`**: Uses OpenMP threads and `std::partial_sort` to construct `m = 15` nearest stream-network neighborhoods across ~91,000 reaches in under a second (`100×–300×` faster than R tree walks).

### 3.3 Numerically Stable `log_AFV` & Physical Confluence Weighting (`--afv=`)
In large national catchments (e.g., Humber or Thames), raw multiplicative Additive Function Values (`AFV`) can underflow floating-point precision after `40–70` headwater bifurcations. [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) and `C_compute_log_afv()` in [`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp) accumulate `log_AFV` additively in double precision and allow switching dynamically among three physical confluence-mixing schemes:
- `length`: cumulative upstream channel length (`US_Accum_km`)
- `shreve`: Shreve link magnitude
- `hydro_yield`: orographic elevation + aquifer baseflow yield proxy.

### 3.4 Combining `Tail-Up` & `Tail-Down`
Upstream `S3N` fits one stream covariance model (`Tail-Up` **or** `Tail-Down`) at a time. This pipeline fits both and combines them either via a unified **Joint Additive NNGP** (`--combine=joint`, default) or via reach-level **Local Spatial Information-Gain Weighting** (`--combine=precision`: $w(s_0) = \max(10^{-4}, 1 - \text{SE}(s_0)^2 / (\hat{\sigma}^2 + \hat{\tau}^2))$), which shifts unmonitored headwaters toward `Tail-Down` while blending both along monitored mainstems.

---

## 4. Additional Features Available in This Pipeline

Alongside the core `S3N` workflow, this codebase provides the following modelling, validation, and operational features:

1. **Joint Additive (`Tail-Up + Tail-Down + Euclidean`) NNGP (`--combine=joint`)**:
   Simultaneously estimates variance partitioning across downstream advective flow ($\sigma^2_{\text{TU}}$), along-channel dispersal ($\sigma^2_{\text{TD}}$), 2D Euclidean overland/atmospheric proximity ($\sigma^2_{\text{Euc}}$), and microscale nugget ($\tau^2$) inside one $m \times m$ Cholesky system (`C_s3n_joint_nngp()`).
2. **Multi-Kernel Stream Covariances & Automatic AIC Selection (`--kernel=`)**:
   Supports `exponential`, compact-support `spherical` and `linear`, heavy-tailed `mariah`, and data-driven NNGP AIC selection (`--kernel=auto`) via `C_s3n_nngp_loglik()`.
3. **Multivariate Stream Co-Kriging (`Factor-LMC`)** ([`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R)):
   Jointly models $P \ge 2$ heterotopic chemical and ecological targets via `C_s3n_multivariate_nngp()` (and `--method=colocated`) with row-by-row `NA` masking, Dual CV (`All-Holdout` vs. `Aux-Aided`), cross-correlation matrices (`Corr_TU`, `Corr_TD`, `Corr_Nugget`), and an automatic [`cokriging_evaluation_report.md`](./multivariate_cokriging_guide.md).
4. **Two-Stage GAM & Heteroskedastic Measurement-Error Kriging (`--meas-error` / `--gam`)**:
   Converts first-stage station-level standard errors (`meas_uncertainty`) onto the transformed model scale via the Delta Method, decomposes total nugget variance into $\tau^2_{\text{micro}} + \bar{v}$, and solves the heteroskedastic NNGP kriging system in `C_s3n_predict_meas_error()` (and `C_s3n_joint_nngp()`).
5. **Non-Linear Spline & Tree-Ensemble Regression Kriging (`--trend=spline|gam|rf|gbm`)**:
   Replaces strictly linear covariates with orthogonalized GCV natural cubic splines, REML-penalized GAMs, Out-of-Bag Random Forests, or 5-fold cross-fitted Gradient Boosted Trees prior to S3N residual kriging, preserving spatial GLS inference without overfitting residuals.
6. **Left-Censored & Interval-Range Imputation (`--censor=ros|em|tobit|half_lod`, `--range-mode=`)**:
   Handles `<LOD` nondetects via multi-threshold Robust Lognormal ROS or iterative **Spatial Tobit EM Kriging** (truncated-normal E-steps conditioned on stream neighbors), plus custom interval-range modes (`tobit`, `geometric`, `midpoint`, `ros_site`, etc.).
7. **Spatio-Temporal Batching, Dynamic Borrowing & Trend Mapping (`--by=year|season|month`, `--st-borrow`)**:
   Fits ordered time slices, borrows precision-weighted spatial strength across sparsely monitored periods using estimated lag-1 temporal correlation $\hat{\rho}_t$, and computes reach-level WLS temporal trend significance maps (`fig6` / `fig7`).
8. **Regulatory Exceedance & WFD Risk Mapping (`--threshold=uk|<val>`)**:
   Evaluates exact reach-level Gaussian exceedance probabilities (`exceed_prob`) on the model scale and assigns 5-tier management compliance classifications (`Confidently Compliant` through `Confidently Exceeding`) with dedicated risk maps.
9. **Optimal Monitoring Network Design (`--design-sites=K`)**:
   Runs an OpenMP greedy Schur-complement variance-reduction solver (`C_s3n_optimal_design()`) to rank the unmonitored reaches where new monitoring stations most reduce network-wide prediction uncertainty.
10. **Full Operational Workflow, Hydrological Validation & Diagnostics**:
    - **4 Hydrological CV Modes (`--cv-mode=`)**: `random`, `catchment` (Leave-Catchment-Out), `strahler` (mainstem-to-headwater extrapolation), and `spatial_block`.
    - **12 Response Transforms & 4 Normalisations**: Profile-MLE `boxcox`/`yeojohnson`, `sqrt`, `log`, `log1p`, `asin`, `logit`, etc., plus `zscore`, `robust_zscore`, `minmax`, and `percentile` outputs.
    - **Custom Data & Networks**: Combine built-in (`rsn`, `ssn`, `rss`) and user CSV sources (`+`), bring custom river networks (`BYOWN`), and customize 15 covariates via `--vars=`.
    - **Diagnostics & Re-Plottable Multi-Format Figures**: Automated station snapping diagnostics (`06_station_snapping_diagnostics.R`), descriptive environmental gradient reports (`05_descriptive_stats.R`), multi-run / `N:P` stoichiometry comparisons (`04_compare_model_runs.R`), and standalone re-plotting (`03_generate_plots.R`) in `png`, Cairo `pdf`, and `svg` with `--dpi` and `--fig-size`.

---

## 5. License

Consistent with `S3N` (`GPL-3`) and `BRISC` (`GPL-2 | GPL-3`), all code in this repository is licensed under the **GNU General Public License v3.0 (`GPL-3.0`)** — see [`LICENSE`](../LICENSE). Raw GIS networks and monitoring datasets (`gis/`, `data/`) are not distributed with the repository and are governed by their respective upstream provider licences.

---

## 6. Related Documentation

- [**S3N Modelling & Kriging Theory**](./s3n_modelling_and_kriging_guide.md) — Mathematical derivations of Tail-Up, Tail-Down, Joint NNGP, ROS, hydrological CV, and optimal monitoring design.
- [**Preprocessing & C++ Neighbor Engine**](./preprocessing_pipeline.md) — Deep dive into `01_preprocess_orn.py`, Euler-tour DFS timestamps, and `src/s3n_fast_neighbors.cpp`.
- [**Multivariate Stream Co-Kriging (Factor-LMC)**](./multivariate_cokriging_guide.md) — Joint $P \ge 2$ modelling via `scripts/run_s3n_cokriging.R`.
- [**End-to-End Usage Guide**](./usage_guide.md) — Complete CLI flag reference and workflow examples.
