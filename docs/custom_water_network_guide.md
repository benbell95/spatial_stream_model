# Bring Your Own Water Network Guide

This guide explains what is required in a river/stream network dataset to work with the Scalable Spatial Stream Network (`S3N`) pipeline, how the preprocessor automatically derives missing hydrographic and topological attributes, and step-by-step instructions for adding and running a new water network.

---

## 1. Overview & Architecture

The `S3N` R package models spatial autocorrelation along river networks using two moving-average covariance structures:
- **Tail-Up (TU)**: Flow-connected covariance that requires a **strictly converging directed tree** per catchment (every reach has at most **1 downstream child**, and no loops/cycles exist) plus **additive function values (`afv`, $\pi_i$)** that split weights proportionally at confluences.
- **Tail-Down (TD)**: Flow-connected and flow-unconnected covariance that requires the **hydrologic distance to the catchment outlet (`Dist2Mth_km`)** for every reach.

To decouple the `S3N` modelling scripts from any single GIS dataset (such as OS Open Rivers Network `ORN_v2.gpkg` or Copernicus `EU-Hydro`), the pipeline uses:
1. **Declarative Network Profiles (`configs/networks/<abbr>.json`)**: Maps a dataset's native layer names, CRS, topology columns, and attribute names to the standardized `S3N` schema.
2. **Network-Agnostic Preprocessor (`scripts/01_preprocess_orn.py`)**: Reads the GeoPackage, repairs topology into a valid converging tree, projects coordinates, computes missing hydrological attributes, snaps EA River Surveillance Network (RSN) monitoring stations, and writes standardized CSV tables to `data/processed/<abbr>/`.
3. **Network-Aware Modelling & Plotting (`scripts/run_s3n_determinand.R`, `scripts/02_fit_s3n_models.R`, `scripts/03_generate_plots.R`)**: Reads `data/processed/<abbr>/`, tags output folders with `<abbr>` (e.g. `output/det_0117_euhydro_thames/`), and stamps all figures (`fig1`–`fig5`) with the network's human-readable label.

---

## 2. Dataset Requirements

### 2.1 File Format
- **Supported format**: **OGC GeoPackage (`.gpkg`)** or a `.zip` archive containing a `.gpkg` file.
- **Converting from Shapefile / GeoJSON / FileGDB**: If your river network is in another vector format (`.shp`, `.geojson`, `.gdb`), convert it to a GeoPackage using GDAL's `ogr2ogr` before running the preprocessor:
  ```bash
  ogr2ogr -f GPKG gis/my_network.gpkg path/to/input_rivers.shp -nln rivers
  ```

### 2.2 Geometry Requirements
- **Geometry Type**: `LineString` or `MultiLineString` in standard OGC WKB or ISO WKB:
  - **2D (`XY`)**: Supported. Reach lengths, sinuosity, and endpoint coordinates are computed directly from the 2D vertices.
  - **3D (`XYZ` / `XYM`) & 4D (`XYZM`)**: Supported (e.g., EU-Hydro uses ISO `MultiLineStringZM`, WKB type `3005`). When `Z` coordinates are present on vertices, the preprocessor automatically extracts mean reach `Elevation` (m) and channel `Slope` (%).
- **Digitized Flow Direction**: Polylines should ideally be digitized **upstream $\to$ downstream** (first vertex = upstream start node, last vertex = downstream end node).
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
| **`"nodes"`** | `from_node_col` (upstream junction ID) and `to_node_col` (downstream junction ID) | Reach $A$ flows into reach $B$ when $\text{to\_node}(A) = \text{from\_node}(B)$. | OS Open Rivers Network (`Fnode` $\to$ `Tnode`) |
| **`"pointer"`** | `id_col` (unique reach ID) and `next_down_col` (downstream reach ID) | Reach $A$ flows directly into the reach whose `id_col` equals $\text{next\_down\_col}(A)$. | Copernicus EU-Hydro (`OBJECT_ID` $\to$ `NEXTDOWNID`) |
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
| **`RivID`** (`1..N`) | S3N graph & joins | No | `topology.id_col` (stored as `OrigID`) | Assigned sequential integer IDs `1..N` (`OrigID` preserves original string/int ID) |
| **`NextDownID`** | S3N adjacency & tree | **Yes** (via `nodes`, `pointer`, or `geometry`) | `from_node_col`/`to_node_col` or `next_down_col` | Built from node matching, pointer lookup, or endpoint coordinate snapping |
| **`CatchID`** | Catchment grouping | No | `topology.catch_id_col` | Weakly connected tree components are automatically assigned integer `CatchID`s |
| **`Length_m`** | Routing & `US_Accum_km` | No | Computed from WKB polyline | Exact 2D/3D polyline arc length in meters |
| **`Dist2Mth_km`** | **Tail-Down** distances | No | `attributes.dist_to_mouth_col` | Computed via root-to-leaf tree traversal: $\text{Dist2Mth}(u) = \text{Dist2Mth}(\text{down}(u)) + \text{Length}(u)$ |
| **`US_Accum_km`** | **Tail-Up** `afv` & OLS | No | `attributes.afv_weight_col` | Computed via leaf-to-root tree traversal: $\text{US\_Accum}(u) = \text{Length}(u) + \sum_{p \in \text{up}(u)} \text{US\_Accum}(p)$ |
| **`Shreve`** | Alt `afv` weight | No | `attributes.shreve_col` | Computed via leaf-to-root tree traversal ($1$ at headwaters, $\sum_{p \in \text{up}(u)} \text{Shreve}(p)$ downstream) |
| **`Strahler`** | OLS covariate & map width | No | `attributes.strahler_col` | Computed via leaf-to-root Strahler stream ordering rules |
| **`afv`** | **Tail-Up** weight ($\pi_i$) | No | Derived from `US_Accum_km` | Multiplicative root-to-leaf proportion $\pi_u = \pi_{\text{down}(u)} \cdot \frac{w_u}{\sum_{s \in \text{up}(\text{down}(u))} w_s}$ |
| **`Sinuosity`** | OLS covariate | No | Computed from WKB polyline | $\text{Length\_m} / \max(\text{StraightLineDist\_m}, 1.0) \ge 1.0$ |
| **`Elevation`** | OLS covariate | No | `attributes.elevation_col` or `nodes_table` | Extracted from 3D WKB vertex `Z` coordinates, or spatially interpolated from nearest ORN reach |
| **`Slope`** | OLS covariate | No | `attributes.slope_col` | Extracted from 3D WKB vertex $\Delta Z / \text{Length\_m} \times 100$, or spatially interpolated |
| **`Pct_Chalk` .. `Pct_Peat`** | OLS geology covariates | No | `covariates.mode` | Transferred via fast spatial-hash nearest-reach lookup (`2 km` grid) from BGS catchment geology |

---

## 3. Step-by-Step Setup Guide

### Step 1: Place Your GeoPackage in `gis/` and Inspect Its Schema

Place your `.gpkg` (or `.zip`) file inside the `gis/` directory, then inspect its layer names, geometry column, CRS, and attribute columns using `sqlite3`:

```bash
# 1. List feature layers, geometry column name, and EPSG code
sqlite3 gis/my_network.gpkg "SELECT table_name, column_name, geometry_type_name, srs_id FROM gpkg_geometry_columns;"

# 2. List all columns in the river network layer
sqlite3 gis/my_network.gpkg "PRAGMA table_info(my_river_layer);"

# 3. Preview the first 5 rows of non-geometry attributes
sqlite3 -header -column gis/my_network.gpkg "SELECT FID, OBJECT_ID, FNODE, TNODE, NEXTDOWNID FROM my_river_layer LIMIT 5;"
```

---

### Step 2: Create a Network Profile in `configs/networks/<abbr>.json`

Create a new JSON configuration file in [`configs/networks/`](file:///home/ben/code/r/spatial_stream_model/configs/networks) (for example, `configs/networks/my_net.json`). Below is a complete template with annotations explaining every option:

```json
{
  "name": "my_net",
  "abbr": "my_net",
  "label": "My Custom River Network v1.0",
  "short_label": "MyNet v1",
  "gpkg_paths": [
    "gis/my_network.gpkg"
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
- **`name` / `abbr`**: Short lowercase identifier (e.g. `"euhydro_thames"`). Used in CLI `--network=<abbr>`, cache folder `data/processed/<abbr>/`, and output directory `output/det_0117_<abbr>/`.
- **`label` / `short_label`**: Human-readable network name displayed on Figures 1–5 (`Water Network: <short_label>`) and stored in `Network_Label` inside CSV tables.
- **`gpkg_paths`**: List of one or more relative `.gpkg` paths under the project root. Multiple GeoPackages (such as `euhydro_thames` + `euhydro_tweed` in [`configs/networks/euhydro.json`](file:///home/ben/code/r/spatial_stream_model/configs/networks/euhydro.json)) are automatically merged.
- **`zip_fallbacks`**: Optional `.zip` archive paths; if a `.gpkg` in `gpkg_paths` does not exist yet, the preprocessor automatically extracts it from the corresponding `.zip`.
- **`layers`**: List of feature table names to load from each GeoPackage.
- **`geom_col`**: Name of the GeoPackage WKB geometry column (typically `"Shape"` or `"geom"`).
- **`crs_epsg`**: Integer EPSG code of the layer's native CRS (e.g., `27700` for British National Grid, `3035` for ETRS89-LAEA Europe).
- **`topology`**:
  - `mode`: `"nodes"` (matches `to_node_col` $\to$ `from_node_col`), `"pointer"` (matches `next_down_col` $\to$ `id_col`), `"geometry"` (snaps reach endpoints within `snap_tol_m` meters), or `"auto"`.
  - `id_col`, `from_node_col`, `to_node_col`, `next_down_col`, `catch_id_col`: Column names in the GIS table (set any unused column to `null`).
- **`attributes`**:
  - Set any attribute column that exists in your GIS layer, or set to `null` to have the preprocessor compute it automatically from the network tree and 3D geometry.
  - `dist_to_mouth_unit` / `afv_weight_unit`: `"m"` (automatically converted to `"km"` in `Dist2Mth_km` and `US_Accum_km`) or `"km"`.
- **`nodes_table`**: Optional lookup table in the same GeoPackage (such as EU-Hydro's `"Nodes"` layer) containing node elevations (`elev_col`) or distances to mouth (`dist_to_mouth_col`). Set to `null` if not applicable.
- **`core_basins`**: Optional integer `CatchID` overrides for the Humber, Thames, and Severn basins used in Figures 2 and 3. If set to `null`, the R scripts automatically identify the Humber, Thames, and Severn catchments from `catchment_lookup.csv` by river name (`Main_River`) or river basin district.
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
All processed files are saved to `data/processed/my_net/` (and mirrored to `data/processed/`):
1. **`england_rsn_reaches.csv`**: All reaches in monitored catchments with standardized topology (`RivID`, `NextDownID`, `CatchID`, `afv`), coordinates (`start_X/Y`, `mid_X/Y`, `end_X/Y` in `EPSG:27700`), physical/geological covariates, and baseline nitrate observations (`has_obs`, `nitrate_mean`, `log_nitrate`).
2. **`england_rsn_edges.csv`**: Directed parent-to-child edge list (`from_RivID`, `to_RivID`, `CatchID`) for building the `S3N` adjacency matrix.
3. **`humber_reaches.csv` & `humber_edges.csv`**: Re-indexed `1..N_h` subset for the Humber basin.
4. **`rsn_site_reach_mapping.csv`**: Lookup mapping every snapped EA RSN monitoring site (`Monitoring_Point`, `Eastings`, `Northings`, `snap_dist_m`) to its `RivID` and `CatchID` on your network, used by `scripts/run_s3n_determinand.R` to aggregate any chemical or ecological determinand dynamically.
5. **`catchment_lookup.csv`**: Summary table of every monitored catchment (`CatchID`, `Main_River`, `n_reaches`, `n_obs_reaches`, `max_order`, `outlet_X`, `outlet_Y`).
6. **`network_metadata.json`**: Metadata record storing the network's `abbr`, `label`, `short_label`, reach/station counts, and resolved `core_basins` (`humber`, `thames`, `severn`).

---

### Step 4: Run S3N Models with Your Water Network

Pass `--network=my_net` to [`scripts/run_s3n_determinand.R`](file:///home/ben/code/r/spatial_stream_model/scripts/run_s3n_determinand.R) or [`scripts/02_fit_s3n_models.R`](file:///home/ben/code/r/spatial_stream_model/scripts/02_fit_s3n_models.R):

```bash
# Fit S3N models for Nitrate (0117) on your network
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net

# Fit S3N models for Orthophosphate (0180) on your network with custom covariates
Rscript scripts/run_s3n_determinand.R 0180 --network=my_net --vars=Elevation,US_Accum_km,Strahler

# Fit seasonal S3N models on your network (includes Figure 5 All-England seasonal comparison grid)
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net --by-season

# Include an additional custom 2-panel basin map (e.g. Thames or a specific CatchID)
Rscript scripts/run_s3n_determinand.R 0117 --network=my_net --basin=Thames
```

> **Automatic Preprocessing Trigger**: If `data/processed/my_net/` does not exist yet when you run `Rscript scripts/run_s3n_determinand.R 0117 --network=my_net`, the R script automatically calls `python3 scripts/01_preprocess_orn.py --network=my_net` first.

---

### Step 5: Re-Generate or Customize Plots

All output files are saved in a network-specific folder so they never overwrite runs from other networks:
- Standard covariates: `output/det_0117_<abbr>/` (e.g., [`output/det_0117_euhydro_thames/`](file:///home/ben/code/r/spatial_stream_model/output/det_0117_euhydro_thames))
- Custom covariates: `output/det_0117_<abbr>_vars_<tag>/` (e.g., `output/det_0117_euhydro_thames_vars_3vars_elev_usac_stra/`)

You can re-render Figures 1–5 from any saved output directory using [`scripts/03_generate_plots.R`](file:///home/ben/code/r/spatial_stream_model/scripts/03_generate_plots.R):

```bash
# Re-render Figures 1-4 from saved predictions
Rscript scripts/03_generate_plots.R --dir=output/det_0117_my_net

# Also render a dedicated 2-panel map (Figure 5) for the Thames basin or any CatchID
Rscript scripts/03_generate_plots.R --dir=output/det_0117_my_net --basin=Thames
```

---

## 4. Ad-Hoc Usage (Running Without a Config File)

If you want to test a `.gpkg` file quickly without creating a JSON file in `configs/networks/`, you can pass the `.gpkg` path directly to `--network=` and override topology or weighting options on the command line:

```bash
# Run directly on a GeoPackage file using automatic topology detection
python3 scripts/01_preprocess_orn.py \
  --network=gis/custom_rivers.gpkg \
  --topology-mode=auto \
  --max-snap-dist=300
```

Available CLI overrides in [`scripts/01_preprocess_orn.py`](file:///home/ben/code/r/spatial_stream_model/scripts/01_preprocess_orn.py):
- `--network=<name|path>`: Profile name in `configs/networks/` or direct path to a `.gpkg` file.
- `--network-config=<path.json>`: Path to a custom JSON configuration file outside `configs/networks/`.
- `--topology-mode=auto|nodes|pointer|geometry`: Override how downstream reach connectivity is constructed.
- `--afv-col=<col>`: Override the attribute column used to compute Tail-Up additive function weights (`afv`).
- `--max-snap-dist=<meters>`: Override the maximum station-to-reach snapping distance (default `250`).

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

