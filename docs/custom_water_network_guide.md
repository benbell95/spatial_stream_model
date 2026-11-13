# Bring Your Own Water Network Guide

This guide explains what is required in a river/stream network dataset to work with the Scalable Spatial Stream Network (`S3N`) pipeline, how the preprocessor automatically derives missing hydrographic and topological attributes, and step-by-step instructions for adding and running a new water network.

---

## Table of Contents

- [1. Overview & Architecture](#1-overview--architecture)
- [2. Dataset Requirements](#2-dataset-requirements)
  - [2.1 File Format](#21-file-format)
  - [2.2 Geometry Requirements](#22-geometry-requirements)
  - [2.3 Coordinate Reference System (CRS)](#23-coordinate-reference-system-crs)
  - [2.4 Flow Routing & Topology Requirements](#24-flow-routing--topology-requirements)
  - [2.5 Required vs. Optional vs. Auto-Computed Attributes](#25-required-vs-optional-vs-auto-computed-attributes)
- [3. Step-by-Step Setup Guide](#3-step-by-step-setup-guide)
  - [Step 1: Place Your GeoPackage or GeoParquet File in `gis/` and Inspect Its Schema](#step-1-place-your-geopackage-or-geoparquet-file-in-gis-and-inspect-its-schema)
  - [Step 2: Create a Network Profile in `configs/networks/<abbr>.json`](#step-2-create-a-network-profile-in-configsnetworksabbrjson)
  - [Step 3: Run the Preprocessor on Your Network](#step-3-run-the-preprocessor-on-your-network)
  - [Step 4: Run S3N Models with Your Water Network](#step-4-run-s3n-models-with-your-water-network)
  - [Step 5: Re-Generate or Customize Plots](#step-5-re-generate-or-customize-plots)
- [4. Ad-Hoc Usage (Running Directly on a `.parquet` or `.gpkg` File Without a Config File)](#4-ad-hoc-usage-running-directly-on-a-parquet-or-gpkg-file-without-a-config-file)
- [5. Troubleshooting & Validation Checklist](#5-troubleshooting--validation-checklist)

---

## 1. Overview & Architecture

The `S3N` R package models spatial autocorrelation along river networks using two moving-average covariance structures:
- **Tail-Up (TU)**: Flow-connected covariance that requires a **strictly converging directed tree** per catchment (every reach has at most **1 downstream child**, and no loops/cycles exist) plus **additive function values (`afv`, πᵢ)** that split weights proportionally at confluences.
- **Tail-Down (TD)**: Flow-connected and flow-unconnected covariance that requires the **hydrologic distance to the catchment outlet (`Dist2Mth_km`)** for every reach.

To decouple the `S3N` modelling scripts from any single GIS dataset (such as OS Open Rivers Network `ORN_v2.gpkg` or Copernicus `EU-Hydro`), the pipeline uses:
1. **Declarative Network Profiles (`configs/networks/<abbr>.json`)**: Maps a dataset's native GeoPackage (`.gpkg`) or GeoParquet (`.parquet` / `.geoparquet`) file paths, layer names, CRS, topology columns, and attribute names to the standardized `S3N` schema (`configs/networks/orn.json` is included out of the box; custom profiles can be added alongside it).
2. **Network-Agnostic Preprocessor (`scripts/01_preprocess_orn.py`)**: Reads the GeoPackage or GeoParquet file(s), repairs topology into a valid converging tree, projects coordinates, computes missing hydrological attributes, snaps EA monitoring stations (or custom sources), and writes standardized CSV tables to `data/processed/<abbr>_<source>/`.
3. **Network-Aware Modelling & Plotting (`scripts/run_s3n_determinand.R`, `scripts/02_fit_s3n_models.R`, `scripts/03_generate_plots.R`)**: Reads `data/processed/<abbr>_<source>/`, writes run outputs to `output/results/<run_tag>/` (e.g. `output/results/det_0117_my_net_rsn/`), and stamps all figures with the network's human-readable label.

---

## 2. Dataset Requirements

### 2.1 File Format
- **Supported formats**:
  - **OGC GeoPackage (`.gpkg`)**
  - **OGC GeoParquet (`.parquet`, `.geoparquet`, `.pq`)** — including single GeoParquet files, multiple layer `.parquet` files, or partitioned directories of `.parquet` files
  - **`.zip` archives** containing `.gpkg`, `.parquet`, or `.geoparquet` files.
- **How GeoParquet Files Are Read**:
  - `scripts/01_preprocess_orn.py` inspects the OGC GeoParquet 1.0/1.1 `"geo"` JSON metadata footer (`primary_column`, `encoding`, PROJJSON `crs`) directly from the Parquet footer and streams features and WKB geometries using either `pyarrow.parquet` (if installed) or the system `libgdal` C OGR `Parquet` driver (`ctypes`).
- **Converting from Shapefile / GeoJSON / FileGDB (or GeoPackage to GeoParquet)**:
  - Using `ogr2ogr` (to GeoPackage or GeoParquet):
    ```bash
    # Convert to GeoPackage (.gpkg)
    ogr2ogr -f GPKG gis/my_network.gpkg path/to/input_rivers.shp -nln rivers

    # Convert to GeoParquet (.parquet)
    ogr2ogr -f Parquet gis/my_network.parquet path/to/input_rivers.shp
    ```
  - Or convert any configured GeoPackage network directly to GeoParquet using `scripts/01_preprocess_orn.py`:
    ```bash
    python3 scripts/01_preprocess_orn.py --network=euhydro_thames --convert-to-geoparquet=gis/euhydro_thames.parquet
    ```

### 2.2 Geometry Requirements
- **Geometry Type**: `LineString` or `MultiLineString` in standard OGC WKB or ISO WKB:
  - **2D (`XY`)**: Supported. Reach lengths, sinuosity, and endpoint coordinates are computed directly from the 2D vertices.
  - **3D (`XYZ` / `XYM`) & 4D (`XYZM`)**: Supported (e.g., EU-Hydro uses ISO `MultiLineStringZM`, WKB type `3005`). When `Z` coordinates are present on vertices, the preprocessor automatically extracts mean reach `Elevation` (m) and channel `Slope` (%).
- **Digitized Flow Direction**: Polylines should ideally be digitized **upstream → downstream** (first vertex = upstream start node, last vertex = downstream end node).
  - *Automatic Orientation Repair*: If your dataset has any reaches digitized backwards relative to its topology pointers/nodes, `scripts/01_preprocess_orn.py` automatically detects when the start vertex is closer to the downstream reach than the end vertex and flips the vertex order.

### 2.3 Coordinate Reference System (CRS)
- **Supported CRS**: Any projected planar CRS (in meters) with an EPSG code recognized by `PROJ` / `/usr/bin/cs2cs` (for example, `EPSG:27700` British National Grid or `EPSG:3035` ETRS89-LAEA Europe).
- **How CRS Transformation Works**:
  1. EA RSN monitoring station coordinates (`Eastings`, `Northings` in `EPSG:27700`) are transformed into the water network's native CRS (`crs_epsg`) so station-to-reach perpendicular distance snapping is performed in the network's native planar meters.
  2. Reach coordinates (`start_X`, `start_Y`, `mid_X`, `mid_Y`, `end_X`, `end_Y`) are transformed into `EPSG:27700` so all England and basin maps share a consistent projection.
- *Note*: If your raw network is in geographic coordinates (`EPSG:4326` lat/lon degrees), reproject it to a metric CRS (`EPSG:27700` or `EPSG:3035`) with `ogr2ogr -t_srs EPSG:27700` before running the preprocessor so reach lengths and snap distances are in meters.

### 2.4 Flow Routing & Topology Requirements
Your dataset must provide **at least one** of the following three ways to determine which reach flows into which downstream reach (`topology.mode`):

| Topology Mode | Required Columns in GIS Layer | How Routing is Built | Example Dataset |
| :--- | :--- | :--- | :--- |
| **`"nodes"`** | `from_node_col` (upstream junction ID) and `to_node_col` (downstream junction ID) | Reach A flows into reach B when `to_node(A) == from_node(B)`. | OS Open Rivers Network (`Fnode` → `Tnode`) |
| **`"pointer"`** | `id_col` (unique reach ID) and `next_down_col` (downstream reach ID) | Reach A flows directly into the reach whose `id_col` equals `next_down_col(A)`. | Copernicus EU-Hydro (`OBJECT_ID` → `NEXTDOWNID`) |
| **`"geometry"`** | *None* (uses polyline vertices only) | Snaps each reach's downstream endpoint `(end_X, end_Y)` to the nearest reach's upstream startpoint `(start_X, start_Y)` within tolerance (`snap_tol_m`, default `5.0m`). | Custom digitised polylines without node/pointer attributes |
| **`"auto"`** | Automatically selects `"pointer"`, `"nodes"`, or `"geometry"` based on which columns exist in the layer. | — | Ad-hoc `.gpkg` files |

> **Multi-Layer Networks (e.g., Rivers + Canals + Ditches)**:
> Some datasets split connected waterways across multiple layers (for example, in `euhydro_thames_v013.gpkg`, 147 river segments in `River_Net_l` have `NEXTDOWNID` pointing to canal segments in `Canals_l`). You can list multiple layers in `"layers": ["River_Net_l", "Canals_l", "Ditches_l"]` and the preprocessor will merge them into a single unified routing graph before building trees.

> **Automatic Braided-Channel & Cycle Repair**:
> Natural river datasets often contain braided splits/divergences (where one node has multiple outgoing reaches) or minor topological cycles. The preprocessor automatically resolves divergences (retaining the primary downstream branch with the largest upstream accumulation/length) and breaks cycles so the output graph is guaranteed to be a valid **directed acyclic forest** for `S3N`.

---

### 2.5 Required vs. Optional vs. Auto-Computed Attributes

You do **not** need pre-computed Strahler orders, Shreve orders, distances to mouth, or catchment IDs in your GIS dataset. Anything missing is automatically computed from the geometry and tree topology:

| Standardized Output Column | Used By | Strictly Required in GIS? | Native Column Mapping (`attributes`) | Fallback / Automatic Computation if `null` or Missing |
| :--- | :--- | :---: | :--- | :--- |
| **`RivID`** (`1..N`) | S3N graph & joins | No | `topology.id_col` | Assigned sequential integer IDs `1..N` |
| **`parent_dn`** | S3N adjacency & tree | **Yes** (via `nodes`, `pointer`, or `geometry`) | `from_node_col`/`to_node_col` or `next_down_col` | Built from node matching, pointer lookup, or endpoint coordinate snapping |
| **`CatchID`** | Catchment grouping | No | `topology.catch_id_col` | Weakly connected tree components are automatically assigned integer `CatchID`s |
| **`length_km`** | Routing & `US_Accum_km` | No | Computed from WKB polyline | Exact 2D/3D polyline arc length converted to km |
| **`Dist2Mth_km`** | **Tail-Down** distances | No | `attributes.dist_to_mouth_col` | Computed via root-to-leaf tree traversal: `Dist2Mth(u) = Dist2Mth(down(u)) + Length(u)` |
| **`US_Accum_km`** | **Tail-Up** `log_AFV` & OLS | No | `attributes.afv_weight_col` | Computed via leaf-to-root tree traversal: `US_Accum(u) = Length(u) + sum(US_Accum(p) for p in up(u))` |
| **`Shreve`** | Alt `AFV` weight & width | No | `attributes.shreve_col` | Computed via leaf-to-root tree traversal (`1` at headwaters, `sum(Shreve(p) for p in up(u))` downstream) |
| **`Strahler`** | Optional OLS covariate & map width | No | `attributes.strahler_col` | Computed via leaf-to-root Strahler stream ordering rules |
| **`log_AFV`** | **Tail-Up** weight ($\log \pi_i$) | No | Derived from `US_Accum_km` | Accumulated root-to-leaf log proportion: `log_AFV(u) = log_AFV(down(u)) + log(w(u) / sum(w(s) for s in up(down(u))))` |
| **`Sinuosity`** | Optional OLS covariate | No | Computed from WKB polyline | `Length_m / max(StraightLineDist_m, 1.0) >= 1.0` |
| **`Elevation`** | OLS covariate | No | `attributes.elevation_col` or `nodes_table` | Extracted from 3D WKB vertex `Z` coordinates, or spatially transferred from nearest ORN reach |
| **`Slope`** | OLS covariate | No | `attributes.slope_col` | Extracted from 3D WKB vertex `dZ / Length_m * 100`, or spatially transferred |
| **`Is_Chalk_Limestone` .. `Is_Peat`** | OLS geology covariates | No | `covariates.mode` | Transferred via fast spatial-hash nearest-reach lookup (`2 km` grid) from BGS catchment geology |

---

## 3. Step-by-Step Setup Guide

### Step 1: Place Your GeoPackage or GeoParquet File in `gis/` and Inspect Its Schema

Place your `.gpkg`, `.parquet` / `.geoparquet`, or `.zip` file inside the `gis/` directory (or a custom `--gis-dir`).

- **If using a GeoPackage (`.gpkg`)**, you can inspect its layer names, geometry column, CRS, and attribute columns using `sqlite3`:
  ```bash
  # 1. List feature layers, geometry column name, and EPSG code
  sqlite3 gis/my_network.gpkg "SELECT table_name, column_name, geometry_type_name, srs_id FROM gpkg_geometry_columns;"

  # 2. List all columns in the river network layer
  sqlite3 gis/my_network.gpkg "PRAGMA table_info(my_river_layer);"

  # 3. Preview the first 5 rows of non-geometry attributes
  sqlite3 -header -column gis/my_network.gpkg "SELECT FID, OBJECT_ID, FNODE, TNODE, NEXTDOWNID FROM my_river_layer LIMIT 5;"
  ```
- **If using a GeoParquet (`.parquet` / `.geoparquet`) file**, `scripts/01_preprocess_orn.py` can inspect the GeoParquet `"geo"` metadata footer and column schema automatically:
  ```bash
  python3 -c 'import sys; sys.path.insert(0, "scripts"); import importlib; m = importlib.import_module("01_preprocess_orn"); print(m.inspect_geoparquet_schema("gis/my_network.parquet"))'
  ```

---

### Step 2: Create a Network Profile in `configs/networks/<abbr>.json`

Create a new JSON configuration file in [`configs/networks/`](../configs/networks) (for example, `configs/networks/my_net.json`). You can specify either `"parquet_paths"` (for `.parquet` / `.geoparquet` files), `"gpkg_paths"` (for `.gpkg` or `.parquet` files), or `"network_paths"`. Below is a complete template with annotations explaining every option:

```json
{
  "name": "my_net",
  "abbr": "my_net",
  "label": "My Custom River Network v1.0",
  "short_label": "MyNet v1",
  "parquet_paths": [
    "gis/my_network.parquet"
  ],
  "gpkg_paths": [
    "gis/my_network.parquet"
  ],
  "zip_fallbacks": [
    "gis/my_network.zip"
  ],
  "layers": [
    "my_river_layer"
  ],
  "geom_col": "Shape",
  "crs_epsg": 27700,
  "topology": {
    "mode": "pointer",
    "id_col": "OBJECT_ID",
    "from_node_col": "FNODE",
    "to_node_col": "TNODE",
    "next_down_col": "NEXTDOWNID",
    "catch_id_col": "CatchID",
    "rbd_id_col": null,
    "snap_tol_m": 5.0
  },
  "attributes": {
    "dist_to_mouth_col": null,
    "dist_to_mouth_unit": "m",
    "afv_weight_col": "CUM_LEN",
    "afv_weight_unit": "m",
    "shreve_col": null,
    "strahler_col": "STRAHLER",
    "slope_col": null,
    "elevation_col": null,
    "river_name_cols": ["RiverName"]
  },
  "nodes_table": null,
  "covariates": {
    "mode": "geometry_and_spatial_join",
    "orn_geology_gpkg": "gis/ORN_v2.gpkg"
  },
  "core_basins": {
    "humber": null,
    "thames": null,
    "severn": null
  },
  "max_snap_dist_m": 250.0
}
```

#### Field Reference for `configs/networks/<abbr>.json`
- **`name` / `abbr`**: Short lowercase identifier (e.g. `"my_net"` or `"euhydro_thames"`). Used in CLI `--network=<abbr>`, processed directory `data/processed/<abbr>_<source>/`, and output directory `output/results/det_0117_<abbr>_<source>/`.
- **`label` / `short_label`**: Human-readable network name displayed on figures (`Water Network: <short_label>`) and saved in `network_metadata.json`.
- **`parquet_paths` / `gpkg_paths` / `network_paths`**: List of one or more relative `.parquet` / `.geoparquet` or `.gpkg` paths under `gis/` (or the project root). Multiple files (for example, merging regional sub-basin files or separate GeoParquet layer files) are automatically combined.
- **`zip_fallbacks`**: Optional `.zip` archive paths; if a file in `parquet_paths` / `gpkg_paths` does not exist yet, the preprocessor automatically extracts it from the corresponding `.zip`.
- **`layers`**: List of feature table names to load from each GeoPackage (ignored for single-layer `.parquet` files).
- **`geom_col`**: Name of the WKB geometry column (typically `"Shape"`, `"geometry"`, or `"geom"`; auto-detected from the GeoParquet `"geo"` metadata footer if omitted).
- **`crs_epsg`**: Integer EPSG code of the layer's native CRS (e.g., `27700` for British National Grid, `3035` for ETRS89-LAEA Europe; auto-detected from GeoParquet PROJJSON metadata when running ad-hoc).
- **`topology`**:
  - `mode`: `"nodes"` (matches `to_node_col` → `from_node_col`), `"pointer"` (matches `next_down_col` → `id_col`), `"geometry"` (snaps reach endpoints within `snap_tol_m` meters), or `"auto"`.
  - `id_col`, `from_node_col`, `to_node_col`, `next_down_col`, `catch_id_col`: Column names in the GIS table (set any unused column to `null`).
- **`attributes`**:
  - Set any attribute column that exists in your GIS layer, or set to `null` to have the preprocessor compute it automatically from the network tree and 3D geometry.
  - `dist_to_mouth_unit` / `afv_weight_unit`: `"m"` (automatically converted to `"km"` in `Dist2Mth_km` and `US_Accum_km`) or `"km"`.
- **`nodes_table`**: Optional lookup table containing node elevations (`elev_col`) or distances to mouth (`dist_to_mouth_col`). Can be either a table inside a `.gpkg` (`{"table": "Nodes", ...}`) or a standalone `.parquet` / `.csv` file (`{"path": "gis/my_nodes.parquet", "id_col": "OBJECT_ID", "elev_col": "ELEV", "dist_to_mouth_col": "LEN_TOM"}`). Set to `null` if not applicable.
- **`core_basins`**: Optional integer `CatchID` overrides for the Humber, Thames, and Severn basins used in Figures 2 and 3. If set to `null`, the R scripts automatically identify the Humber, Thames, and Severn catchments from `rsn_catchments_lookup.csv` by river name (`Main_River`) or river basin district.
- **`max_snap_dist_m`**: Maximum distance in meters (default `250.0`) for snapping EA RSN water quality and biological monitoring stations to their nearest river reach.

---

### Step 3: Run the Preprocessor on Your Network

 first verify that your profile appears in `--list-networks`, then run `scripts/01_preprocess_orn.py`:

```bash
# Verify the network profile is registered
python3 scripts/01_preprocess_orn.py --list-networks

# Run preprocessing for your network
python3 scripts/01_preprocess_orn.py --network=my_net
```

#### What `scripts/01_preprocess_orn.py --network=my_net` Produces
All processed files are saved to `data/processed/my_net_rsn/` (or `data/processed/<abbr>_<source>/` for any monitoring source):
1. **`england_rsn_reaches.csv`**: All reaches in monitored catchments with standardized topology (`RivID`, `CatchID`, `parent_dn`, `Dist2Mth_km`, `up_dist_mid_km`, `US_Accum_km`, `log_AFV`, `Shreve`, `Strahler`), coordinates (`start_X/Y`, `mid_X/Y`, `end_X/Y` in `EPSG:27700`), physical/geological covariates (`Elevation`, `Slope`, `Sinuosity`, `SrcElevation`, `Is_Chalk_Limestone`, `Is_Sandstone`, `Is_Aquifer`, `Is_Peat`), and canonical baseline observation columns (`has_obs`, `site_id`, `snap_dist_m`, `obs_up_dist_km`, `obs_mean`, `obs_median`, `trans_obs_mean`, `meas_uncertainty`, `meas_method`, `n_samples`).
2. **`rsn_sites.csv`**: Station-level snapping coordinates and observation summaries (`site_id`, `easting`, `northing`, `snap_dist_m`, `RivID`, `CatchID`, `obs_up_dist_km`, `n_samples`, `obs_mean`, `obs_median`, `obs_sd`, `trans_obs_mean`, `mean_trans_obs`, `meas_uncertainty`, `meas_method`, `included_250m`, `source`), used by `scripts/run_s3n_determinand.R` and `scripts/06_station_snapping_diagnostics.R`.
3. **`rsn_reach_periods.csv`**: Reach-level observations aggregated across temporal slices (`all`, `year`, `season`, `month`, `year_season`).
4. **`rsn_catchments_lookup.csv`**: Summary table of every monitored catchment (`CatchID`, `Main_River`, `RBD_ID`, `RBD_Name`, `n_reaches`, `n_obs_reaches`).
5. **`network_metadata.json`**: Metadata record storing the network's `abbr`, `label`, `short_label`, `source`, reach/station counts, transform metadata, and resolved `core_basins` (`humber`, `thames`, `severn`).

*(The C++ neighbor engine `src/s3n_fast_neighbors.cpp` builds the directed tree directly in memory from `parent_dn` without separate edge-list CSVs.)*

---

### Step 4: Run S3N Models with Your Water Network

Pass `--network=my_net` to [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) or [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R):

```bash
# Fit S3N models for Nitrate (0117) on your network
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net

# Fit S3N models for Orthophosphate (0180) on your network with custom covariates
Rscript scripts/run_s3n_determinand.R 0180 --network=my_net --vars=elevation,us_accum,strahler

# Fit seasonal S3N models on your network (includes All-England seasonal comparison grid)
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net --by=season

# Include an additional custom 2-panel basin map (e.g. Thames or a specific CatchID)
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net --basin=Thames
```

> **Automatic Preprocessing Trigger**: If `data/processed/my_net_rsn/` does not exist yet when you run `Rscript scripts/run_s3n_determinand.R 0117 --network=my_net`, the R script automatically calls `python3 scripts/01_preprocess_orn.py --network=my_net` first.

---

### Step 5: Re-Generate or Customize Plots

All output files are saved in a network- and source-specific folder under `output/results/` so runs never overwrite other configurations:
- Standard covariates: `output/results/det_0117_<abbr>_<source>/` (e.g. `output/results/det_0117_my_net_rsn/`)
- Custom covariates: `output/results/det_0117_<abbr>_<source>_vars_<tag>/` (e.g. `output/results/det_0117_my_net_rsn_vars_elev_accum_strahler/`)

You can re-render Figures 1–4 (and custom `--basin` or `--norm` maps) from any saved output directory using [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R):

```bash
# Re-render Figures 1-4 from saved predictions
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_my_net_rsn

# Also render a dedicated 2-panel map for the Thames basin or any CatchID
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_my_net_rsn --basin=Thames
```

---

## 4. Ad-Hoc Usage (Running Directly on a `.parquet` or `.gpkg` File Without a Config File)

If you want to run a `.parquet` / `.geoparquet` or `.gpkg` file directly without creating a JSON file in `configs/networks/`, pass the file path directly to `--network=`. Both [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) and [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) will automatically inspect the GeoParquet or GeoPackage schema (`geom_col`, `crs_epsg`, topology columns, and hydrological attributes):

```bash
# 1. Preprocess a GeoParquet file directly using automatic schema & topology detection
python3 scripts/01_preprocess_orn.py \
  --network=gis/custom_rivers.parquet \
  --topology-mode=auto \
  --max-snap-dist=300

# 2. Or run S3N modelling & plotting directly on a GeoParquet file (auto-triggers preprocessing)
Rscript scripts/run_s3n_determinand.R 0117 --network=gis/custom_rivers.parquet

# 3. Preprocess a GeoPackage file directly
python3 scripts/01_preprocess_orn.py \
  --network=gis/custom_rivers.gpkg \
  --topology-mode=auto \
  --max-snap-dist=300
```

Available CLI overrides in [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py):
- `--network=<name|path>`: Profile name in `configs/networks/` or direct path to a `.parquet`, `.geoparquet`, `.pq`, or `.gpkg` file.
- `--network-config=<path.json>`: Path to a custom JSON configuration file outside `configs/networks/`.
- `--convert-to-geoparquet=<out.parquet>`: Converts the selected `--network` GeoPackage layer(s) into a single OGC GeoParquet file.
- `--topology-mode=auto|nodes|pointer|geometry`: Override how downstream reach connectivity is constructed.
- `--afv-col=<col>`: Override the attribute column used to compute Tail-Up additive function weights (`afv`).
- `--max-snap-dist=<meters>`: Override the maximum station-to-reach snapping distance (default `250`).
- `--gis-dir=<dir>`: Custom directory containing `.gpkg`, `.parquet` / `.geoparquet`, or `.zip` GIS files (defaults to `configs/directories.json` or `gis`).
- `--data-dir=<dir>`: Custom root data directory containing RSN inputs and `processed/` outputs (defaults to `configs/directories.json` or `data`).
- `--out-dir=<dir>` (or `--output-dir=<dir>`): Custom output directory for preprocessed network CSV tables (defaults to `<data-dir>/processed`).

---

## 5. Troubleshooting & Validation Checklist

1. **Too Few Observed RSN Stations Snapped (`n_obs_reaches` is low)**:
   - Verify that `"crs_epsg"` in your JSON config matches the actual coordinate reference system of the geometries in `gpkg_geometry_columns`.
   - If your river network is generalized/coarse (e.g., 1:250,000 scale), increase `"max_snap_dist_m"` from `250.0` to `500.0` (or pass `--max-snap-dist=500`) so monitoring stations along small tributaries or wide floodplains snap to the nearest mapped mainstem.
2. **Broken Trees or Unresolved `NEXTDOWNID` Pointers**:
   - Check the preprocessor console log for `Loaded ... total reaches | ... valid downstream links`.
   - If your network splits rivers and artificial connectors across separate layers (like EU-Hydro's `River_Net_l`, `Canals_l`, and `Ditches_l`), include all connected layers in the `"layers"` array of your JSON config so downstream pointers across layer boundaries resolve cleanly.
   - If your dataset lacks node/pointer columns or has corrupted IDs, set `"mode": "geometry"` in `"topology"` to build downstream links purely from polyline endpoint proximity.
3. **Inspecting Available Catchments (`CatchID`) on a New Network**:
   - Because each water network defines its own `CatchID` numbering, run:
     ```bash
     Rscript scripts/run_s3n_determinand.R --network=my_net --list-basins
     ```
     to view the top monitored river basins, their `CatchID`, main river name, total reaches, and number of observed RSN stations on your network.

