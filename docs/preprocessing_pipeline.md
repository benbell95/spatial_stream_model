# Technical Documentation: Scalable Spatial Stream Network (S3N) Pipeline

This document provides complete technical documentation for the three-stage **Scalable Spatial Stream Network (S3N)** codebase:
1. **Part I — GIS & Water Quality Preprocessing** ([`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py))
2. **Part II — OpenMP C++ Stream-Tree Neighbor Engine** ([`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp))
3. **Part III — R S3N Estimation, Prediction & Cross-Validation** ([`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) & [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R))

> **Full Mathematical & Kriging Guide**: For the complete step-by-step mathematical derivation of the `S3N` covariance functions (`Tail-Up` & `Tail-Down`), Nearest-Neighbor Gaussian Process (`NNGP`) estimation (`BRISC_estimation_stream`), Universal Stream Kriging (`BRISC_prediction_stream`), and spatial information-gain weighting, see [`./s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md).

---

## Table of Contents

- [Part I: Python GIS & Monitoring Preprocessing (`01_preprocess_orn.py`)](#part-i-python-gis--monitoring-preprocessing-01_preprocess_ornpy)
  - [1. Multi-Network Monitoring Data Ingestion (`--source=<name>`)](#1-multi-network-monitoring-data-ingestion---sourcename)
  - [2. Zero-Dependency GeoPackage WKB Geometry Decoding (`parse_gpkg_multilinestring()`)](#2-zero-dependency-geopackage-wkb-geometry-decoding-parse_gpkg_multilinestring)
  - [3. Spatial Grid Snapping & Exact Along-Reach Projection (`point_to_polyline()`)](#3-spatial-grid-snapping--exact-along-reach-projection-point_to_polyline)
  - [4. Subsetting to Monitored Catchments & Aggregating 100 m Covariates](#4-subsetting-to-monitored-catchments--aggregating-100-m-covariates)
  - [5. Directed River Tree Construction & Additive Function Values](#5-directed-river-tree-construction--additive-function-values)
- [Part II: OpenMP C++ Stream-Tree Neighbor & NNGP Engine (`src/s3n_fast_neighbors.cpp`)](#part-ii-openmp-c-stream-tree-neighbor--nngp-engine-srcs3n_fast_neighborscpp)
  - [1. O(1) Flow-Connectedness via Euler-Tour DFS Timestamps (`C_build_tree_tin_tout()` L27–70)](#1-o1-flow-connectedness-via-euler-tour-dfs-timestamps-c_build_tree_tin_tout-l2770)
  - [2. Pairwise Stream Distance & Weight Kernel (`calc_pair` L99–140)](#2-pairwise-stream-distance--weight-kernel-calc_pair-l99140)
  - [3. Compressed Sparse Row (CSR) Neighbor Buffers (`C_compute_s3n_neighbors()` L142–302)](#3-compressed-sparse-row-csr-neighbor-buffers-c_compute_s3n_neighbors-l142302)
- [Part III: R S3N Estimation, Prediction & Cross-Validation (`02_fit_s3n_models.R` & `run_s3n_determinand.R`)](#part-iii-r-s3n-estimation-prediction--cross-validation-02_fit_s3n_modelsr--run_s3n_determinandr)
- [Part IV: System Requirements, Environment & Dependencies](#part-iv-system-requirements-environment--dependencies)

---

## Part I: Python GIS & Monitoring Preprocessing (`01_preprocess_orn.py`)

[`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) transforms the raw **10.4 GB Ordnance Survey Open Rivers Network (`ORN_v2.gpkg`)** and **Environment Agency Monitoring Network** tables into a topologically ordered, model-ready spatial stream network dataset in **~13–15 seconds** using pure Python standard library modules (`sqlite3`, `struct`, `math`, `csv`).

### 1. Multi-Network Monitoring Data Ingestion (`--source=<name>`)
The preprocessor supports official statutory networks, custom field campaigns, and multi-network combinations via `--source=<name>`:
- `rsn` (*default*, `configs/sources/rsn.json`): Environment Agency **River Surveillance Network** (`data/NCEA_River_Surveillance_Network_ARD_Open/`) monitoring 1,382 baseline stations across England's mainstems and major tributaries.
- `ssn` (`configs/sources/ssn.json`): Environment Agency **Small Streams Network** (`data/NCEA_Small_Streams_Network_ARD_Open/`) monitoring 313 dedicated headwater stream stations (< 20 km² catchment area).
- `rss` (`configs/sources/composites.json`): **Combined River Surveillance + Small Streams Network** integrating both datasets (1,695 unique monitoring stations) for complete headwater-to-mouth stream network coverage (`all` also resolves to `rsn+ssn` by default and includes any user-registered sources in `configs/sources/`).
- `<custom>`: Any custom monitoring dataset registered in `configs/sources/<name>.json`, placed in `data/custom/<name>/`, or supplied ad-hoc via `--custom-sites` and `--custom-samples`.
- `<A>+<B>`: On-the-fly combination merging multiple sources (e.g. `--source=rsn+my_study` or `--source=rsn,ssn,custom`, normalised internally with `+` so hyphens inside source names such as `waterblitz-tobit` are preserved).

> **Complete Custom Data Guide**: For complete instructions on bringing your own point observation datasets, creating `sites.csv` and `samples.csv`, registering profiles in `configs/sources/`, and handling lat/lon coordinates, see [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md).

```bash
# List supported monitoring sources and data directories
python3 scripts/01_preprocess_orn.py --list-sources

# Run preprocessing for RSN (default), SSN, combined RSS, or custom sources
python3 scripts/01_preprocess_orn.py                       # data/processed/orn_rsn/ (default)
python3 scripts/01_preprocess_orn.py --source=ssn          # data/processed/orn_ssn/
python3 scripts/01_preprocess_orn.py --source=rss          # data/processed/orn_rss/
python3 scripts/01_preprocess_orn.py --source=rsn+my_study # data/processed/orn_rsn+my_study/
python3 scripts/01_preprocess_orn.py --custom-sites=path/to/sites.csv --custom-samples=path/to/samples.csv
```

Custom base paths (`data_dir`, `gis_dir`, `output_dir`) are resolved via `configs/directories.json`, environment variables (`S3N_DATA_DIR`, `S3N_GIS_DIR`, `S3N_OUTPUT_DIR`, `S3N_PROCESSED_DIR`), or `--data-dir` / `--gis-dir` / `--out-dir`.

- **Station Lookup (`load_monitoring_stations_and_records`)**: Reads `*Water Quality Sample Points.csv` (or custom `sites.csv`) and `*Water Quality Samples.csv` (or custom `samples.csv`) from active source directories, mapping each `network_id` to British National Grid (`EPSG:27700`) coordinates in meters, tracking origin network (`source`), and preserving interval/uncertainty metadata (`meas_uncertainty`, `meas_method`) when available (`--range-mode {auto,tobit,geometric,midpoint,lower,upper,none}`).
- **Determinand, Transformation & Temporal Filtering (`--det-code`, `--transform`, `--transform-lambda`, `--year`, `--season`, `--month`)**:
  - Parses ISO-8601 `sample_datetime` (`YYYY-MM-DDTHH:MM:SS`) into calendar `year` (`2021`–`2024`), `month` (`1..12`), and meteorological `season` (`winter` = Dec–Feb, `spring` = Mar–May, `summer` = Jun–Aug, `autumn` = Sep–Nov).
  - Supports running on all data (default), by `--year`, by `--season`, or by `--month`, applies the requested data transformation (`log` by default for WQ, configurable via `--list-transforms`) to compute `trans_obs_mean` and `mean_trans_obs`, and exports `rsn_reach_periods.csv` with reach-level aggregations for all periods (`all`, `year`, `season`, `month`, `year_season`) plus `network_metadata.json` for downstream scripts.
- **Configurable Left-Censored (`<`) Handling (`--censor-method {ros,ros_site,half_lod}`)**:
  - `ros` (default): Pooled **Robust Lognormal Regression on Order Statistics** (Hirsch & Stedinger 1987; Helsel & Cohn 1988) across multiple detection limits (`compute_ros_plotting_positions` and `fit_ros_ols`).
  - `ros_site`: Station-level Robust Lognormal ROS within each station (when ≥ 3 uncensored detections exist and < 80% censoring), with automatic fallback to the global ROS regression parameters.
  - `half_lod`: Standard 0.5 × LOD substitution (`0.5 * meas_result`).

### 2. Zero-Dependency GeoPackage WKB Geometry Decoding (`parse_gpkg_multilinestring()`)
- Unpacks the OGC GeoPackage Binary (`GPB`) header (`flags = blob[3]`, `env_ind = (flags >> 1) & 0x07`) to read the 32-byte bounding box `[minx, maxx, miny, maxy]` at byte offset `8` (`struct.unpack_from("<4d", blob, 8)`) in ~1 microsecond without decoding vertex arrays unless a monitoring station is nearby.
- Decodes little-endian `LineString` (`wkb_type = 2`) and `MultiLineString` (`wkb_type = 5`) double-precision vertex arrays directly via `struct.unpack_from`.

### 3. Spatial Grid Snapping & Exact Along-Reach Projection (`point_to_polyline()`)
- **2 km Spatial Hash Grid (`site_grid`)**: Buckets all candidate stations (1,382 for RSN; 313 for SSN; 1,695 for combined RSS) into 2 km × 2 km grid cells `(int(px // 2000), int(py // 2000))`.
- **Orthogonal Segment Projection**: Projects candidate stations onto each segment `p_i -> p_{i+1}` (digitized upstream-to-downstream from `Fnode` to `Tnode`), computing Euclidean snapping distance `d` and cumulative distance `along` from the upstream `Fnode`.
- **Exact Station Distance to River Mouth (`obs_up_dist_km`)**: Subtracting `along` from `Dist2Mth` yields the exact hydrologic distance from the snapped station to the sea/estuary mouth: `max(0.0, (Dist2Mth - along) / 1000.0)`.
- **Configurable Snapping Threshold (`--max-snap-dist <meters>`, default 250 m)**: Retains stations within `d <= max_snap_dist` of a reach:
  - **RSN** (*default*): 1,349 of 1,382 stations snapped (median distance 8.0 m), spanning 1,210 unique observed reaches across 141 catchments (90,876 reaches total).
  - **SSN**: 197 of 313 stations snapped (median distance 10.2 m), spanning 177 unique reaches across 58 catchments (76,242 reaches total).
  - **Combined RSS**: 1,546 of 1,695 stations snapped (median distance 8.2 m), spanning 1,384 unique reaches across 150 catchments (91,118 reaches total).
  - Threshold can be relaxed (e.g. `--max-snap-dist 500`) for coarse networks or edge-case meander stations, or tightened (e.g. `--max-snap-dist 100`) for dense parallel ditch systems. Detailed snapping coordinates, catchment summaries, and diagnostic maps can be generated via `scripts/06_station_snapping_diagnostics.R` (saved under `output/station_snapping/<net>_<source>/`).

### 4. Subsetting to Monitored Catchments & Aggregating 100 m Covariates
- Aggregates the 5.98M-row `ORN_Sampled_100m` table via SQL `GROUP BY PolylineID` across monitored catchments (141 catchments / 90,876 reaches for RSN; 58 catchments / 76,242 reaches for SSN; 150 catchments / 91,118 reaches for RSS) to extract mean channel elevation (`Elevation`), gradient (`Slope`), sinuosity (`Sinuosity`), source elevation (`SrcElevation`), and BGS bedrock/superficial geology indicators (`Is_Chalk_Limestone`, `Is_Sandstone`, `Is_Aquifer`, `Is_Peat`).
- Walks halfway (`0.5 * reach_len`) along each reach polyline to interpolate the exact geometric midpoint `(mid_X, mid_Y)`.

### 5. Directed River Tree Construction & Additive Function Values
- Links each reach's downstream node `(CatchID, Tnode)` to the upstream node `(CatchID, Fnode)` of its immediate downstream reach (`parent_dn`).
- Traverses each catchment tree from river mouth (`log_AFV = 0.0`) to headwaters via depth-first search (`topo_order`), accumulating log confluence weights `log_afv[p] = log_afv[u] + log(M[p] / sum(M[q] for q in parents[u]))` according to `--afv-scheme=length|shreve|hydro_yield` (default `length` = `US_Accum_km`; `shreve` = Shreve link magnitude; `hydro_yield` = upstream length weighted by elevation/orographic rainfall and Principal Aquifer baseflow). R modeling scripts can also recompute `log_AFV` on-the-fly in C++ (`C_compute_log_afv`) via `--afv=length|shreve|hydro_yield`.
- Runs 3 forward/backward tree-neighbor imputation passes to copy covariates for short reaches (`< 100 m`) that lacked a 100 m sample point, leaving 0 missing values.
- Writes canonical outputs into `data/processed/<net>_<source>/`: `england_rsn_reaches.csv` (using canonical `obs_mean`, `obs_median`, `trans_obs_mean`, `meas_uncertainty`, `meas_method`), `rsn_sites.csv` (also includes `obs_sd` and `mean_trans_obs`), `rsn_reach_periods.csv`, `rsn_catchments_lookup.csv`, and `network_metadata.json`.

---

## Part II: OpenMP C++ Stream-Tree Neighbor & NNGP Engine (`src/s3n_fast_neighbors.cpp`)

[`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp) replaces `SSN2::createDistMat()` and `S3N::s3n_neighbors()`, computing all exact stream distances, flow-connectedness flags, NNGP neighbor structures, multi-kernel NNGP likelihoods (`C_s3n_nngp_loglik`), and Joint Additive `TU + TD + Euclidean` kriging (`C_s3n_joint_nngp`) across all 90,876 reaches in memory in **~0.88 seconds** on 4 OpenMP threads.

### 1. O(1) Flow-Connectedness via Euler-Tour DFS Timestamps (`C_build_tree_tin_tout()` L27–70)
- Runs a linear `O(q)` depth-first search from every river mouth root (`parent_dn < 0`), recording entry time `tin[u]` and exit time `tout[u]` for every reach `u`.
- Reach `b` lies in the upstream watershed subtree of reach `a` if and only if `tin[a] <= tin[b] && tout[b] <= tout[a]`, reducing flow-connectedness checks to two integer comparisons (`O(1)` time).

### 2. Pairwise Stream Distance & Weight Kernel (`calc_pair` L99–140)
- **Same Catchment & Flow-Connected**:
  - Hydrologic distance: `dist = max(1e-4, abs(up_dist_km[a] - up_dist_km[b]))`.
  - Weight: `wt = exp(0.5 * (min(log_afv[a], log_afv[b]) - max(log_afv[a], log_afv[b])))` for Tail-Up (`model_type == 0`), or `wt = 1.0` for Tail-Down (`model_type == 1`).
- **Same Catchment & Flow-Unconnected (Bifurcated Tributaries)**:
  - Walks `anc = parent_dn[anc]` downstream from `a` until `tin[anc] <= tin[b] && tout[b] <= tout[anc]` to locate the Lowest Common Ancestor confluence reach `anc`.
  - Hydrologic distance across the confluence junction: `dist = max(1e-4, up_dist_km[a] + up_dist_km[b] - 2.0 * d2m_km[anc])`.
  - Tail-Up sets `fc = false, wt = 0.0`; Tail-Down sets `fc = true, wt = 1.0`.
- **Different Catchments**:
  - Sets `fc = false, wt = 0.0` and `dist = 1e6 + euclidean_dist_km` for Vecchia zero-weight padding.

### 3. Compressed Sparse Row (CSR) Neighbor Buffers (`C_compute_s3n_neighbors()` L142–302)
- Partially sorts candidate neighbors (`std::partial_sort`) using `struct Cand` (prioritizing `fc == true`, then smallest `dist`) to populate `obs_neighbors` (`d`, `nnWt`, `nnIndx`, `nnIndxLU`, `D`, `Dwt`, `CIndx`, `Length.D`) and `pred_neighbors` (`nnIndx`, `nnDist`, `nnWght`).

---

## Part III: R S3N Estimation, Prediction & Cross-Validation (`02_fit_s3n_models.R` & `run_s3n_determinand.R`)

1. **Downstream-to-Upstream Vecchia Ordering**: Orders observed reaches by `CatchID` and increasing `obs_up_dist_km` (`1..n_obs`), followed by prediction reaches (`n_obs+1..q_tot`).
2. **Multi-Kernel Maximum Likelihood Estimation (`--kernel=exponential|spherical|mariah|linear|auto`)**: Fits both S3N Tail-Up (`model_type = 0`) and S3N Tail-Down (`model_type = 1`) via sparse NNGP profile maximum likelihood (`BRISC_estimation_stream` or `C_s3n_nngp_loglik`), supporting Exponential, Spherical, Mariah, Linear-with-Sill, and automatic AIC selection (`--kernel=auto`), plus 95% confidence intervals.
3. **Joint Additive (`TU + TD + Euclidean`) & Variance-Weighted Stream Kriging (`--combine=joint|precision`)**:
   - In **`joint` mode** (*default*, `--combine=joint`), estimates joint variance components $(\sigma^2_{\text{TU}}, \sigma^2_{\text{TD}}, \sigma^2_{\text{Euc}}, \tau^2)$ and solves a single unified Cholesky kriging system (`C_s3n_joint_nngp`) combining flow-connected `Tail-Up`, hydrologic-distance `Tail-Down`, and 2D overland `Euclidean` covariance across all 90,876 reaches.
   - In **`precision` mode** (`--combine=precision`), combines separate `Tail-Up` and `Tail-Down` predictions using local spatial information-gain weights `info = pmax(1e-4, 1 - se^2 / (sigma.sq + tau.sq))`.
4. **5-Fold Out-of-Sample Cross-Validation**: Reorders each fold's training reaches (`1..n_tr`) ahead of held-out test reaches (`n_tr+1..n_obs`) to evaluate out-of-sample `RMSE`, `MAE`, and `R^2` against non-spatial OLS.

---

## Part IV: System Requirements, Environment & Dependencies

- **Software & Compilation**: `01_preprocess_orn.py` requires only **Python ≥ 3.8** (standard library only, zero `pip` packages). R scripts require **R ≥ 3.5.0**, `S3N` (plus its CRAN dependencies), optional `jsonlite`, and [`src/s3n_fast_neighbors.so`](../src/s3n_fast_neighbors.so) (`R CMD SHLIB src/s3n_fast_neighbors.cpp`). Shared R helper modules live in [`scripts/s3n_config.R`](../scripts/s3n_config.R), [`scripts/s3n_covariates.R`](../scripts/s3n_covariates.R), [`scripts/s3n_transforms.R`](../scripts/s3n_transforms.R), [`scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R), and [`scripts/s3n_plot_functions.R`](../scripts/s3n_plot_functions.R). See [`./setup_and_architecture.md`](./setup_and_architecture.md) for full OS package tables, `configs/determinands.json` schema, and test commands.
- **Required Data Files**:
  1. [`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg) (~10.4 GB Ordnance Survey Open Rivers Network v2 GeoPackage — only required when running [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py)).
  2. [`data/NCEA_River_Surveillance_Network_ARD_Open/`](../data/NCEA_River_Surveillance_Network_ARD_Open) (~235 MB across EA RSN Water Quality, Macroinvertebrate, Diatom, Macrophyte, River Habitat Survey, and LCMS/GCMS CSV files).
  3. [`data/NCEA_Small_Streams_Network_ARD_Open/`](../data/NCEA_Small_Streams_Network_ARD_Open) (~48 MB across EA SSN Small Streams Network CSV files — required when using `--source=ssn` or `--source=rss`).
- **Memory & CPU**:
  - **RAM**: ~1.5 GB peak memory (fits comfortably on any laptop or standard VM).
  - **CPU**: Uses 4 OpenMP threads by default (`n_omp = 4L` / `--omp=4`); a full national + 3-basin S3N run with 5-fold cross-validation and 4 figures completes in **~18–25 seconds**.

