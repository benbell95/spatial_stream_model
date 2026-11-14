# Open Rivers Network (`ORN_v2.gpkg`) Variables Reference for S3N Modelling

This document details all variables extracted from the **Ordnance Survey Open Rivers Network v2** ([`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg)) during preprocessing ([`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py)) and explains how each variable is used across the **Scalable Spatial Stream Network (S3N)** modelling, spatial covariance construction, and plotting stages.

---

## Table of Contents

- [1. Source Tables in `ORN_v2.gpkg`](#1-source-tables-in-orn_v2gpkg)
- [2. Spatial Stream Network Topology & Covariance Variables](#2-spatial-stream-network-topology--covariance-variables)
- [3. Fixed-Effect Regression Covariates (`X`)](#3-fixed-effect-regression-covariates-x)
- [4. Variables Extracted for Visualisation, Basin Selection, Optional Predictors & Reference](#4-variables-extracted-for-visualisation-basin-selection-optional-predictors--reference)
- [5. Selecting & Comparing Fixed-Effect Covariates (`--vars=...` & `--list-vars`)](#5-selecting--comparing-fixed-effect-covariates---vars----list-vars)
  - [5.1 Available Fixed-Effect Terms (`COVARIATE_CATALOG`)](#51-available-fixed-effect-terms-covariate_catalog)
  - [5.2 Example Commands for Comparing Fixed-Effect Specifications](#52-example-commands-for-comparing-fixed-effect-specifications)
  - [5.3 Plot Subtitles & Differentiated Output Folders When `--vars` Is Used](#53-plot-subtitles--differentiated-output-folders-when---vars-is-used)

---

## 1. Source Tables in `ORN_v2.gpkg`

[`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) queries five tables inside [`gis/ORN_v2.gpkg`](../gis/ORN_v2.gpkg) (configured in [`configs/networks/orn.json`](../configs/networks/orn.json)):

| GeoPackage Table | Rows (National) | Role in [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) |
| :--- | ---: | :--- |
| **`ORN`** | 183,322 reaches (90,876 in the 141 monitored RSN catchments) | Primary river network layer providing reach geometry (`Shape`), topological nodes (`Fnode`, `Tnode`), catchment/basin identifiers (`CatchID`, `RBD_ID`), hydrologic distances (`Dist2Mth`, `US_Accum`), stream orders (`Shreve`, `Strahler`, `Hack`), and watercourse names (`name1`). |
| **`ORN_Sampled_100m`** | ~5,980,000 points | High-resolution 100 m sample points along reaches, aggregated via `GROUP BY PolylineID` to obtain reach-level elevation, channel slope, sinuosity, source elevation, and BGS bedrock/superficial geology codes. |
| **`tblBedrockGeologyLookUp`** | Lookup table | Maps integer `LEX_RCS_BEDROCK` codes from `ORN_Sampled_100m` to British Geological Survey (BGS) bedrock lithology descriptions (`CHALK`, `LIMESTONE`, `DOLOSTONE`, `SANDSTONE`, etc.). |
| **`tblSuperficialGeologyLookUp`** | Lookup table | Maps integer `LEX_RCS_SUPERFICIAL` codes from `ORN_Sampled_100m` to BGS superficial deposit descriptions (`PEAT`, etc.). |
| **`tblRiverBasinLookUp`** | 11 River Basin Districts | Maps `RBD_ID` to Environment Agency River Basin District names (`Humber`, `Thames`, `Severn`, `Anglian`, `North West`, etc.). |

---

## 2. Spatial Stream Network Topology & Covariance Variables

These variables are passed to the compiled C++/OpenMP neighbour engine ([`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp)) to construct the directed river tree, test flow-connectedness, compute pairwise hydrologic stream distances, and evaluate Tail-Up and Tail-Down spatial covariance weights:

| Exported Variable | Source Column(s) in `ORN_v2.gpkg` | Derivation in [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) | Role in S3N Modelling ([`s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp), [`s3n_model_functions.R`](../scripts/s3n_model_functions.R)) |
| :--- | :--- | :--- | :--- |
| **`RivID`** | `ORN.RivID` (`ORN_Sampled_100m.PolylineID`) | Unique integer identifier for each river reach. | Primary key linking snapped monitoring stations (`data/processed/orn_<source>/rsn_sites.csv`) to stream reaches and indexing parent–child links. |
| **`CatchID`** | `ORN.CatchID` | Unique integer identifier for each independent river catchment draining to a sea/estuary mouth. | Defines connected river trees (`net_id[a] == net_id[b]` in `C_compute_s3n_neighbors()`); reaches in different `CatchID` trees have zero spatial covariance (`fc = false, wt = 0.0`). Also used to order reaches by catchment and filter individual basins (e.g. Humber `CatchID=80`, Thames `CatchID=3055`, Severn `CatchID=2786`). |
| **`parent_dn`** | `ORN.Fnode`, `ORN.Tnode`, `ORN.CatchID` | Matches each reach's downstream node `(CatchID, Tnode)` to the upstream node `(CatchID, Fnode)` of its immediate downstream reach; set to `0` (`-1` in 0-based C indexing) for river mouth roots. | Used by `C_build_tree_tin_tout()` to compute Euler-tour DFS entry/exit timestamps (`tin`, `tout`) for `O(1)` flow-connectedness checks (`tin[a] <= tin[b] && tout[b] <= tout[a]`), and walked downstream to locate the Lowest Common Ancestor (LCA) confluence `anc` for flow-unconnected tributaries. |
| **`Dist2Mth_km`** | `ORN.Dist2Mth` | Hydrologic distance (converted from m to km) from the reach's upstream node (`Fnode`) to the sea/estuary river mouth. | Provides the exact distance to the mouth (`ncj_up = d2m_km[anc]`) at the confluence junction `anc` where two bifurcated tributaries meet, enabling exact across-confluence stream distance calculation: `d_stream(a, b) = d_up(a) + d_up(b) - 2 * Dist2Mth_km(anc)`. Also selectable as an optional covariate (`--vars=+dist2mth`). |
| **`obs_up_dist_km`** | `ORN.Dist2Mth`, `ORN.Shape` | Exact hydrologic distance (km) from a snapped monitoring station to the river mouth: `max(0, (Dist2Mth - s_along) / 1000)`, where `s_along` is the distance from the upstream `Fnode` to the station's orthogonal projection along the reach polyline. | 1. Orders observed reaches downstream-to-upstream (`order(CatchID, obs_up_dist_km)`) for Vecchia NNGP conditioning sets.<br>2. Provides `up_dist` for observed reaches (`has_obs == 1`) when computing pairwise hydrologic distances `abs(d_up(a) - d_up(b))`. |
| **`up_dist_mid_km`** | `ORN.Dist2Mth`, `ORN.Shape` | Hydrologic distance (km) from the geometric midpoint of the reach to the river mouth: `max(0, Dist2Mth_km - 0.5 * length_km)`. | Provides `up_dist` for unobserved prediction reaches (`has_obs == 0`) when computing stream distances between prediction reaches and their `m = 15` nearest observed neighbours. |
| **`log_AFV`** | `ORN.US_Accum`, `ORN.Fnode`, `ORN.Tnode` | Natural logarithm of the **Additive Function Value (AFV)**. Initialised at `log_AFV(mouth) = 0` and accumulated upstream at every confluence in proportion to each tributary's upstream accumulated network length `A_p = max(1e-4, US_Accum_p)`:<br>`log_AFV(p) = log_AFV(u) + log(A_p / sum(A_q for q in parents(u)))` | Used in the **S3N Tail-Up (`model_type = 0`)** covariance kernel (`s3n_fast_neighbors.cpp`) to compute Ver Hoef & Peterson (2010) flow-connected confluence weights:<br>`wt(a, b) = sqrt(AFV_up / AFV_dn) = exp(0.5 * (min(log_AFV[a], log_AFV[b]) - max(log_AFV[a], log_AFV[b])))` |

---

## 3. Fixed-Effect Regression Covariates (`X`)

These variables are assembled into the standardised design matrix `X` by [`build_design_matrix()`](../scripts/s3n_covariates.R) (defined in [`scripts/s3n_covariates.R`](../scripts/s3n_covariates.R) and shared across [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) and [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R)) for Non-Spatial OLS, `BRISC_estimation_stream`, and `BRISC_prediction_stream`:

| Predictor in Design Matrix `X` | Exported Column in `england_rsn_reaches.csv` | Source in `ORN_v2.gpkg` | Derivation & Transformation |
| :--- | :--- | :--- | :--- |
| **`Elevation_z`** | `Elevation` | `ORN_Sampled_100m.OS_T50_Ele` | Mean OS Terrain 50 elevation (m above ordnance datum) across 100 m sample points on the reach (`AVG(OS_T50_Ele)`), z-score standardised (`scale(Elevation)`). Also used in Figure 4D to classify **Upland (≥ 150 m)** vs. **Lowland (< 150 m)** reaches. |
| **`Slope_z`** | `Slope` | `ORN_Sampled_100m.Slope` | Mean channel gradient (%) along the reach (`AVG(Slope)`), z-score standardised (`scale(Slope)`). |
| **`log10_US_Accum_z`** | `US_Accum_km` | `ORN.US_Accum` | Total upstream accumulated channel length (converted from m to km), transformed via `log10(pmax(1e-4, US_Accum_km))` and z-score standardised as a proxy for catchment drainage area and river discharge. |
| **`Is_Chalk_Limestone`** *(All England default)* | `Is_Chalk_Limestone` | `ORN_Sampled_100m.LEX_RCS_BEDROCK` + `tblBedrockGeologyLookUp` | Binary indicator (`1` or `0`) equal to `1` if the reach's dominant BGS bedrock description contains `CHALK`, `LIMESTONE`, or `DOLOSTONE`. Captures carbonate/chalk groundwater baseflow enrichment. |
| **`Is_Sandstone`** *(All England default)* | `Is_Sandstone` | `ORN_Sampled_100m.LEX_RCS_BEDROCK` + `tblBedrockGeologyLookUp` | Binary indicator (`1` or `0`) equal to `1` if the reach's dominant BGS bedrock description contains `SANDSTONE` and does not contain `MUDSTONE` or `SILTSTONE`. Captures Permo-Triassic and other principal sandstone aquifers. |
| **`Is_Peat`** *(All England default)* | `Is_Peat` | `ORN_Sampled_100m.LEX_RCS_SUPERFICIAL` + `tblSuperficialGeologyLookUp` | Binary indicator (`1` or `0`) equal to `1` if the reach's dominant BGS superficial deposit description contains `PEAT`. Captures upland blanket bog / organic soil regimes. |
| **`Is_Aquifer`** *(Single-basin default)* | `Is_Aquifer` | Derived from `Is_Chalk_Limestone` and `Is_Sandstone` | Composite binary indicator (`1` if `Is_Chalk_Limestone == 1` or `Is_Sandstone == 1`, else `0`). Used in regional/single-basin fits where individual geological classes may have zero variance, and in Figure 4D to stratify **Lowland Non-Aquifer** vs. **Lowland Aquifer** reaches. Also addable manually via `--vars=+aquifer`. |
| **`Easting_z`**, **`Northing_z`** | `mid_X`, `mid_Y` | `ORN.Shape` (WKB `LineString` / `MultiLineString`) | British National Grid (`EPSG:27700`) coordinates (m) interpolated at the exact half-length (`0.5 * L`) of the reach polyline, z-score standardised (`scale(mid_X)`, `scale(mid_Y)`) to capture large-scale climatic and agricultural gradients. Also used as Euclidean tie-breaker distances (`1e6 + d_eucl`) when padding zero-weight neighbours across separate catchments. |

> [!NOTE]
> **Automatic Collinearity / Zero-Variance Guard:** [`build_design_matrix()`](../scripts/s3n_covariates.R) automatically checks `apply(X[obs_idx, -1, drop = FALSE], 2, sd) > 0` and drops any indicator column that is constant across the observed reaches of a sub-basin.
>
> **Short-Reach Covariate Imputation:** For short reaches (< 100 m) that did not contain a point in `ORN_Sampled_100m`, [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) runs 3 forward/backward passes along the river tree (`topo_order`) to inherit `Elevation`, `Slope`, `Sinuosity`, `SrcElevation`, and geology codes from immediate downstream/upstream connected reaches, leaving 0 missing values across all 90,876 RSN reaches.

---

## 4. Variables Extracted for Visualisation, Basin Selection, Optional Predictors & Reference

The remaining ORN variables exported into `data/processed/orn_<source>/england_rsn_reaches.csv` by [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) support map rendering, custom `--basin` filtering, optional fixed-effect terms (`--vars=+sinuosity,+strahler`, etc.), or auxiliary reference:

| Exported Variable | Source in `ORN_v2.gpkg` | Where Used |
| :--- | :--- | :--- |
| **`start_X`, `start_Y`, `end_X`, `end_Y`** | `ORN.Shape` (`poly[0]`, `poly[-1]`) | Upstream (`Fnode`) and downstream (`Tnode`) reach endpoint coordinates in `EPSG:27700` (m). Used by [`plot_stream_network()`](../scripts/s3n_plot_functions.R) to draw the two half-segments `(start_X, start_Y) -> (mid_X, mid_Y)` and `(mid_X, mid_Y) -> (end_X, end_Y)` for every reach. |
| **`Shreve`** | `ORN.Shreve` | Shreve stream magnitude (number of upstream headwater sources). Used in [`plot_stream_network()`](../scripts/s3n_plot_functions.R) to scale plotted river line widths (`lwd = pmin(max_lwd, base_lwd + scale_lwd * log10(pmax(1, Shreve)))`), in [`01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) to identify the mainstem river name (`Main_River`), and as an optional covariate (`--vars=+shreve`). |
| **`Strahler`** | `ORN.Strahler` | Strahler stream order (`1..7`). Exported as a hydromorphological attribute and available as an optional covariate (`--vars=+strahler`). |
| **`RBD_ID`** & **`RBD_Name`** | `ORN.RBD_ID` + `tblRiverBasinLookUp` | Environment Agency River Basin District identifier and name. Exported to `england_rsn_reaches.csv` and `rsn_catchments_lookup.csv`; used by `resolve_basin_spec()` when running `--basin=<RBD_Name>` or `--basin=rbd:<RBD_ID>`. |
| **`name1` (`Main_River`)** | `ORN.name1` | Watercourse name of the highest-Shreve/accumulation reach in each catchment. Saved as `Main_River` in `rsn_catchments_lookup.csv` to support `--basin=<RiverName>` (e.g. `--basin=Eden`, `--basin="Great Ouse"`). |
| **`length_km`** | `ORN.Shape` | Euclidean polyline length of the reach (km), computed by summing vertex segment lengths. Used during preprocessing to interpolate `(mid_X, mid_Y)` and `up_dist_mid_km`, and available as an optional covariate (`--vars=+length`). |
| **`Sinuosity`** | `ORN_Sampled_100m.Sinuosity` | Mean channel sinuosity ratio (`AVG(Sinuosity)`). Exported as an auxiliary attribute and available as an optional covariate (`--vars=+sinuosity`). |
| **`SrcElevation`** | `ORN_Sampled_100m.OS_T50_SrcEle` | Mean source elevation (`AVG(OS_T50_SrcEle)`). Exported as an auxiliary attribute and available as an optional covariate (`--vars=+src_elevation`). |

---

## 5. Selecting & Comparing Fixed-Effect Covariates (`--vars=...` & `--list-vars`)

Because [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) already exports all 15 candidate topographic, hydromorphological, geological, and coordinate variables into `data/processed/orn_<source>/england_rsn_reaches.csv`, **you do not need to re-run preprocessing** to change which predictors enter the fixed-effect regression design matrix `X`.

Both [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) and [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) accept `--vars=<list>` (or `--covariates=<list>`) and `--list-vars` (powered by [`scripts/s3n_covariates.R`](../scripts/s3n_covariates.R)):

### 5.1 Available Fixed-Effect Terms (`COVARIATE_CATALOG`)
| Model Term | Short CLI Aliases (`--vars=...`) | Default Status | Description |
| :--- | :--- | :--- | :--- |
| `Elevation_z` | `elevation`, `elev` | **Default (All domains)** | Mean channel elevation (`Elevation`, z-score) |
| `Slope_z` | `slope` | **Default (All domains)** | Mean channel gradient (`Slope`, z-score) |
| `log10_US_Accum_z` | `us_accum`, `accum` | **Default (All domains)** | `log10(US_Accum_km)` (z-score) |
| `Is_Chalk_Limestone` | `chalk_limestone`, `chalk` | **Default (All England)** | Binary indicator (`0`/`1`) for Chalk / Limestone / Dolostone bedrock |
| `Is_Sandstone` | `sandstone`, `sand` | **Default (All England)** | Binary indicator (`0`/`1`) for Sandstone aquifer bedrock |
| `Is_Peat` | `peat` | **Default (All England)** | Binary indicator (`0`/`1`) for Peat superficial deposits |
| `Is_Aquifer` | `aquifer` | **Default (Single-basin fits)** | Composite binary indicator (`0`/`1`) for Chalk/Limestone or Sandstone |
| `Easting_z` | `easting`, `mid_x` | **Default (All domains)** | Reach midpoint BNG Easting (`mid_X`, z-score) |
| `Northing_z` | `northing`, `mid_y` | **Default (All domains)** | Reach midpoint BNG Northing (`mid_Y`, z-score) |
| `Sinuosity_z` | `sinuosity` | Optional (`--vars=+sinuosity`) | Mean channel sinuosity ratio (`Sinuosity`, z-score) |
| `SrcElevation_z` | `src_elevation`, `src_elev` | Optional (`--vars=+src_elevation`) | Mean headwater source elevation (`SrcElevation`, z-score) |
| `Dist2Mth_z` | `dist2mth` | Optional (`--vars=+dist2mth`) | Hydrologic distance to river mouth (`Dist2Mth_km`, z-score) |
| `Length_z` | `length` | Optional (`--vars=+length`) | Reach polyline length (`length_km`, z-score) |
| `log10_Shreve_z` | `shreve` | Optional (`--vars=+shreve`) | `log10(Shreve)` stream magnitude (z-score) |
| `Strahler_z` | `strahler` | Optional (`--vars=+strahler`) | Strahler stream order (`Strahler`, z-score) |

### 5.2 Example Commands for Comparing Fixed-Effect Specifications
```bash
# 1. List all available fixed-effect covariates and their aliases
Rscript scripts/run_s3n_determinand.R --list-vars

# 2. Default covariates (Elevation, Slope, US_Accum, Geology, Easting, Northing)
Rscript scripts/run_s3n_determinand.R 0117

# 3. Pick a custom subset of covariates (e.g. topography + composite aquifer, no coordinates)
Rscript scripts/run_s3n_determinand.R 0117 --vars=elevation,slope,aquifer

# 4. Relative removal (-var): start from the default set and drop Easting & Northing
Rscript scripts/run_s3n_determinand.R 0117 --vars=-easting,-northing

# 5. Relative addition (+var) or mixed (+/-): add Sinuosity & Strahler order (or swap in +aquifer)
Rscript scripts/run_s3n_determinand.R 0117 --vars=+sinuosity,+strahler
Rscript scripts/run_s3n_determinand.R 0117 --vars=-chalk,-sand,+aquifer

# 6. Intercept-only model (pure spatial stream network kriging with no environmental covariates)
Rscript scripts/run_s3n_determinand.R 0117 --vars=none

# 7. All candidate ORN covariates simultaneously (omitting Is_Aquifer to avoid collinearity with Chalk+Sandstone)
Rscript scripts/run_s3n_determinand.R 0117 --vars=all
```

### 5.3 Plot Subtitles & Differentiated Output Folders When `--vars` Is Used
When `--vars=<spec>` is set to anything other than `default`:
1. **Plot Subtitles on Figures 1–3 (and Custom Basin / Temporal Grid Maps)**:
   - Figures 1, 2, and 3 (as well as any custom `--basin` maps and temporal comparison grid figures) automatically display the active fixed-effect covariates in both the panel subtitle (`Covariates: Elevation, Slope, Aquifer`) and a figure-level top header (`Fixed-Effect Regression Covariates: ...`).
   - [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R) also detects non-default covariate sets directly from `parameters_<code>.csv` when re-rendering figures from a saved folder.
2. **Differentiated Output Folder Names (`output/results/...`)**:
   - Instead of overwriting the default output directory (`output/results/det_0117_orn_rsn/`), [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) appends a descriptive `_vars_<tag>` suffix to the run folder name:
     - `--vars=elevation,slope,aquifer` → `output/results/det_0117_orn_rsn_vars_elev_slope_aquifer/`
     - `--vars=-easting,-northing` → `output/results/det_0117_orn_rsn_vars_no_east_no_north/`
     - `--vars=+sinuosity,+strahler` → `output/results/det_0117_orn_rsn_vars_add_sinuosity_add_strahler/`
     - `--vars=none` → `output/results/det_0117_orn_rsn_vars_none/`
     - `--vars=all` → `output/results/det_0117_orn_rsn_vars_all/`
   - Similarly, [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) writes to `output/results/orn_rsn_vars_<tag>/` when `--vars` is passed.



