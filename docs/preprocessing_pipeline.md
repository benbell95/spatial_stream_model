# Technical Documentation: Scalable Spatial Stream Network (S3N) Pipeline

This document provides complete technical documentation for the three-stage **Scalable Spatial Stream Network (S3N)** codebase:
1. **Part I — GIS & Water Quality Preprocessing** ([`scripts/01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py))
2. **Part II — OpenMP C++ Stream-Tree Neighbor Engine** ([`src/s3n_fast_neighbors.cpp`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.cpp))
3. **Part III — R S3N Estimation, Prediction & Cross-Validation** ([`scripts/02_fit_s3n_models.R`](file:///home/ben/code/r/spatial_stream_model/scripts/02_fit_s3n_models.R) & [`scripts/run_s3n_determinand.R`](file:///home/ben/code/r/spatial_stream_model/scripts/run_s3n_determinand.R))

---

## Part I: Python GIS & Water Quality Preprocessing (`01_preprocess_orn.py`)

[`scripts/01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py) transforms the raw **10.4 GB Ordnance Survey Open Rivers Network (`ORN_v2.gpkg`)** and **Environment Agency River Surveillance Network (`RSN`)** tables into a topologically ordered, model-ready spatial stream network dataset in **~13 seconds** using pure Python standard library modules (`sqlite3`, `struct`, `math`, `csv`).

### 1. Water Quality Ingestion & Censored-Value Handling (`ROS` vs. `0.5 * LOD`)
- **Station Lookup (`pts`)**: Reads `RSN Water Quality Sample Points.csv` and maps each `network_id` to its British National Grid (`EPSG:27700`) coordinates `(ngr_easting, ngr_northing)` in meters.
- **Determinand & Temporal Filtering (`--det-code`, `--year`, `--season`, `--month`)**:
  - Parses ISO-8601 `sample_datetime` (`YYYY-MM-DDTHH:MM:SS`) into calendar `year` (`2021`–`2024`), `month` (`1..12`), and meteorological `season` (`winter` = Dec–Feb, `spring` = Mar–May, `summer` = Jun–Aug, `autumn` = Sep–Nov).
  - Supports running on all data (default), by `--year`, by `--season`, or by `--month`, and exports `data/processed/rsn_reach_periods.csv` with reach-level aggregations for all periods (`all`, `year`, `season`, `month`, `year_season`).
- **Configurable Left-Censored (`<`) Handling (`--censor-method {ros,ros_site,half_lod}`)**:
  - `ros` (default): Pooled **Robust Lognormal Regression on Order Statistics** (Hirsch & Stedinger 1987; Helsel & Cohn 1988) across multiple detection limits (`compute_ros_plotting_positions` and `fit_ros_ols`).
  - `ros_site`: Station-level Robust Lognormal ROS within each station (when $\ge 3$ uncensored detections exist and $< 80\%$ censoring), with automatic fallback to the global ROS regression parameters.
  - `half_lod`: Standard $0.5 \times \text{LOD}$ substitution (`0.5 * meas_result`).

### 2. Zero-Dependency GeoPackage WKB Geometry Decoding (`parse_gpkg_multilinestring()` L35–L59)
- Unpacks the OGC GeoPackage Binary (`GPB`) header (`flags = blob[3]`, `env_ind = (flags >> 1) & 0x07`) to read the 32-byte bounding box `[minx, maxx, miny, maxy]` at byte offset `8` (`struct.unpack_from("<4d", blob, 8)`) in 1 microsecond without decoding vertex arrays unless a monitoring station is nearby.
- Decodes little-endian `LineString` (`wkb_type = 2`) and `MultiLineString` (`wkb_type = 5`) double-precision vertex arrays directly via `struct.unpack_from`.

### 3. Spatial Grid Snapping & Exact Along-Reach Projection (`point_to_polyline()` L62–L84, `main()` L106–L180)
- **2 km Spatial Hash Grid (`site_grid`)**: Buckets all 1,382 stations into 2 km × 2 km grid cells `(int(px // 2000), int(py // 2000))`.
- **Orthogonal Segment Projection**: Projects candidate stations onto each segment `p_i -> p_{i+1}` (digitized upstream-to-downstream from `Fnode` to `Tnode`), computing Euclidean snapping distance `d` and cumulative distance `along` from the upstream `Fnode`.
- **Exact Station Distance to River Mouth (`obs_up_dist_km`)**: Subtracting `along` from `Dist2Mth` yields the exact hydrologic distance from the snapped station to the sea/estuary mouth: `max(0.0, (Dist2Mth - along) / 1000.0)`.
- **250 m Snapping Threshold**: Retains stations within `<= 250 m` of a reach (1,349 of 1,382 stations, median snap distance 8.0 m, spanning 1,210 unique observed reaches across 141 catchments).

### 4. Subsetting to Monitored Catchments & Aggregating 100 m Covariates (`main()` L181–L236)
- Aggregates the 5.98M-row `ORN_Sampled_100m` table via SQL `GROUP BY PolylineID` across the 141 monitored catchments (90,876 reaches) to extract mean channel elevation (`Elevation`), gradient (`Slope`), sinuosity (`Sinuosity`), source elevation (`SrcElevation`), and BGS bedrock/superficial geology indicators (`Is_Chalk_Limestone`, `Is_Sandstone`, `Is_Aquifer`, `Is_Peat`).
- Walks halfway (`0.5 * reach_len`) along each reach polyline to interpolate the exact geometric midpoint `(mid_X, mid_Y)`.

### 5. Directed River Tree Construction & Additive Function Values (`main()` L238–L281)
- Links each reach's downstream node `(CatchID, Tnode)` to the upstream node `(CatchID, Fnode)` of its immediate downstream reach (`parent_dn`).
- Traverses each catchment tree from river mouth (`log_AFV = 0.0`) to headwaters via depth-first search (`topo_order`), accumulating log confluence weights proportional to upstream accumulated channel length (`US_Accum`):
  `log_afv[p] = log_afv[u] + log(US_Accum[p] / sum(US_Accum[q] for q in parents[u]))`
- Runs 3 forward/backward tree-neighbor imputation passes to copy covariates for short reaches (`< 100 m`) that lacked a 100 m sample point, leaving 0 missing values.

---

## Part II: OpenMP C++ Stream-Tree Neighbor Engine (`src/s3n_fast_neighbors.cpp`)

[`src/s3n_fast_neighbors.cpp`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.cpp) replaces `SSN2::createDistMat()` and `S3N::s3n_neighbors()`, computing all exact stream distances, flow-connectedness flags, and NNGP neighbor structures across all 90,876 reaches in memory in **~0.88 seconds** on 4 OpenMP threads.

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
2. **Maximum Likelihood Estimation (`BRISC_estimation_stream`)**: Fits both S3N Exponential Tail-Up (`model_type = 0`) and S3N Exponential Tail-Down (`model_type = 1`) via sparse NNGP L-BFGS-B optimization, plus parametric bootstrap 95% CIs (`BRISC_bootstrap`).
3. **Variance-Weighted Combined Stream Kriging (`BRISC_prediction_stream`)**: Combines Tail-Up and Tail-Down predictions across all 90,876 reaches using spatial information-gain weights `info = pmax(1e-4, 1 - se^2 / (sigma.sq + tau.sq))`, seamlessly blending mainstem flow-connected kriging (`TU`) with cross-tributary headwater borrowing (`TD`).
4. **5-Fold Out-of-Sample Cross-Validation**: Reorders each fold's training reaches (`1..n_tr`) ahead of held-out test reaches (`n_tr+1..n_obs`) to evaluate out-of-sample `RMSE`, `MAE`, and `R^2` against non-spatial OLS.

---

## Part IV: System Requirements, Environment & Dependencies

The entire S3N pipeline was deliberately engineered to minimize external software dependencies so that it runs reliably in headless Linux environments without heavy GIS or plotting stacks (`GDAL`, `geopandas`, `sf`, `SSN2`, or `ggplot2`).

### 1. Python Environment & Packages ([`scripts/01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py))
- **Python Version**: **Python $\ge 3.8$** (tested on `Python 3.14.8`). Note: `statistics.NormalDist()` used in the Robust Lognormal ROS censored-value imputation was introduced in Python 3.8.
- **Third-Party Python Packages**: **None (Zero `pip` or `conda` dependencies)**.
- **Standard Library Modules Used**:
  | Module | Role in [`01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py) |
  | :--- | :--- |
  | `sqlite3` | Queries the 10.4 GB [`gis/ORN_v2.gpkg`](file:///home/ben/code/r/spatial_stream_model/gis/ORN_v2.gpkg) OGC GeoPackage directly via SQL without `GDAL`/`fiona`/`geopandas`. |
  | `struct` | Unpacks little-endian binary GeoPackage headers and WKB `LineString`/`MultiLineString` coordinate arrays. |
  | `statistics` | Computes means, medians, standard deviations, and normal quantiles (`NormalDist().inv_cdf`) for ROS imputation. |
  | `math`, `collections`, `csv`, `argparse`, `os`, `time` | Geometry projection, tree traversal, CSV I/O, and CLI argument parsing. |

---

### 2. R Environment & Packages ([`scripts/run_s3n_determinand.R`](file:///home/ben/code/r/spatial_stream_model/scripts/run_s3n_determinand.R), [`scripts/02_fit_s3n_models.R`](file:///home/ben/code/r/spatial_stream_model/scripts/02_fit_s3n_models.R), [`scripts/03_generate_plots.R`](file:///home/ben/code/r/spatial_stream_model/scripts/03_generate_plots.R))
- **R Version**: **R $\ge 3.5.0$** (tested on `R 4.5.3` on `x86_64-redhat-linux-gnu`).
- **Required R Packages** (installed in the default R package library):
  | Package | Version / Source | Location | Purpose |
  | :--- | :--- | :--- | :--- |
  | **`S3N`** | `1.0.5` (`remotes::install_github("jpierkunke/S3N")`) | Standard R Library (`.libPaths()`) | Scalable Spatial Stream Network NNGP estimation (`BRISC_estimation_stream`), kriging prediction (`BRISC_prediction_stream`), and parametric bootstrap (`BRISC_bootstrap`) (Kunke et al., 2024; `https://github.com/jpierkunke/S3N`). |
  | `RANN`, `rdist`, `matrixStats`, `pbapply`, `dplyr` | CRAN | Standard R Library (`.libPaths()`) | Package dependencies required by the upstream `jpierkunke/S3N` package namespace. |
  | `stats`, `graphics`, `grDevices`, `utils`, `parallel` | Base R | Bundled with R | All linear algebra (`qr.solve`), statistical functions, and publication-quality base-R PNG graphics (`png(..., type = "cairo")`, `segments`, `points`, `barplot`, `boxplot`). |
- **No `ggplot2` / `sf` / `SSN2` Required**: All maps and diagnostic plots use pure base R graphics, and all stream-tree neighbor calculations use the compiled C++ OpenMP engine ([`src/s3n_fast_neighbors.cpp`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.cpp)).

---

### 3. Compiled C++ / OpenMP Shared Library ([`src/s3n_fast_neighbors.cpp`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.cpp))
- **Compiled Binary**: [`src/s3n_fast_neighbors.so`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.so), loaded dynamically by the R scripts via `dyn.load(file.path(BASE_DIR, "src/s3n_fast_neighbors.so"))`.
- **Recompilation Command** (if migrating to a new CPU architecture or OS):
  ```bash
  PKG_CXXFLAGS="-fopenmp -O3" PKG_LIBS="-fopenmp" R CMD SHLIB src/s3n_fast_neighbors.cpp
  ```

---

### 4. System / OS Dependencies (Linux)
| Dependency | Fedora / RHEL Package | Ubuntu / Debian Package | Purpose |
| :--- | :--- | :--- | :--- |
| **C/C++ Compiler & OpenMP** | `gcc`, `gcc-c++`, `libgomp` | `build-essential`, `libomp-dev` | Compiling [`src/s3n_fast_neighbors.cpp`](file:///home/ben/code/r/spatial_stream_model/src/s3n_fast_neighbors.cpp) and the `S3N` C/C++ BRISC backend with multi-threaded OpenMP (`-fopenmp`). |
| **R Development Headers & BLAS/LAPACK** | `R-core-devel`, `openblas-devel`, `lapack-devel` | `r-base-dev`, `libopenblas-dev`, `liblapack-dev` | Provides `<R.h>`, `<Rinternals.h>`, `<R_ext/Lapack.h>`, and `<R_ext/BLAS.h>` for `S3N` and `s3n_fast_neighbors.so`. |
| **Cairo Graphics** | `cairo`, `cairo-devel` | `libcairo2`, `libcairo2-dev` | Enables anti-aliased headless PNG rendering (`type = "cairo"`) in R without an X11 display server. |
| **SQLite3** | `sqlite-libs` | `libsqlite3-0` | Underlying C library used by Python's built-in `sqlite3` module to read `ORN_v2.gpkg`. |

---

### 5. Input Data & Hardware Requirements
- **Required Data Files**:
  1. [`gis/ORN_v2.gpkg`](file:///home/ben/code/r/spatial_stream_model/gis/ORN_v2.gpkg) (~10.4 GB Ordnance Survey Open Rivers Network v2 GeoPackage — only required when running [`scripts/01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py)).
  2. [`data/NCEA_River_Surveillance_Network_ARD_Open/`](file:///home/ben/code/r/spatial_stream_model/data/NCEA_River_Surveillance_Network_ARD_Open) (~235 MB total across the EA RSN Water Quality, Macroinvertebrate, Diatom, Macrophyte, River Habitat Survey, and LCMS/GCMS CSV files).
- **Memory & CPU**:
  - **RAM**: ~1.5 GB peak memory (fits comfortably on any laptop or standard VM).
  - **CPU**: Uses 4 OpenMP threads by default (`n_omp = 4L`); a full national + 3-basin S3N run with 5-fold cross-validation and 4 figures completes in **~18–25 seconds**.

