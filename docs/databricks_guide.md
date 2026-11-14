# Running the Pipeline on Databricks

Databricks notebooks are well suited to this hybrid pipeline because they natively support **multi-language execution (`%python`, `%r`, `%sh`) in the same notebook**.

---

## Table of Contents

- [1. Repository Setup: Databricks Git Folders](#1-repository-setup-databricks-git-folders)
- [2. Notebook Architecture (Cell-by-Cell Flow)](#2-notebook-architecture-cell-by-cell-flow)
  - [Cell 1: Compile C++ Engine & Check Dependencies (`%sh` and `%r`)](#cell-1-compile-c-engine--check-dependencies-sh-and-r)
  - [Cell 2: Configure Large Data Storage (Unity Catalog Volumes)](#cell-2-configure-large-data-storage-unity-catalog-volumes)
  - [Cell 3: GIS & Data Preprocessing (`%python`)](#cell-3-gis--data-preprocessing-python)
  - [Cell 4: Stream Network Modelling & Interactive Inspection (`%r`)](#cell-4-stream-network-modelling--interactive-inspection-r)
  - [Cell 5: Model Comparison & Diagnostics (`%r`)](#cell-5-model-comparison--diagnostics-r)
  - [Cell 6: Inline Visualisation (`%python`)](#cell-6-inline-visualisation-python)
- [3. Databricks-Specific Advantages](#3-databricks-specific-advantages)

---

## 1. Repository Setup: Databricks Git Folders

1. Clone this repository into **Databricks Git Folders**:
   - In your Databricks workspace, go to **Workspace** → **Repos / Git Folders** → **Add** → Git URL.
2. When you open a notebook inside the Git Folder, the working directory is automatically set to the repo root (`/Workspace/Users/.../spatial_stream_model`). All relative paths (`scripts/`, `configs/`, `data/`) resolve naturally.

---

## 2. Notebook Architecture (Cell-by-Cell Flow)

Set the notebook default language to **R** (or **Python**) and transition between steps using language magics.

### Cell 1: Compile C++ Engine & Check Dependencies (`%sh` and `%r`)

`scripts/01_preprocess_orn.py` uses only the Python standard library for `.gpkg`/CSV preprocessing (optional `pyarrow`/`pyproj` for GeoParquet or CRS reprojection). On the R side, install `S3N`, its CRAN dependencies, recommended `data.table` (fast CSV I/O), and `mgcv`/`rpart` (if using `--trend=gam|rf|gbm`), and compile `src/s3n_fast_neighbors.so` on the driver node:

```sh
%sh
R CMD SHLIB src/s3n_fast_neighbors.cpp
```

```r
required_pkgs <- c("RANN", "rdist", "matrixStats", "pbapply", "dplyr",
                   "jsonlite", "data.table", "remotes", "mgcv", "rpart")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
if (!requireNamespace("S3N", quietly = TRUE)) remotes::install_github("jpierkunke/S3N")
```

*(See [`./setup_and_architecture.md`](./setup_and_architecture.md) for full Linux OS and `jpierkunke/S3N` installation details.)*

---

### Cell 2: Configure Large Data Storage (Unity Catalog Volumes)

On Databricks, store large GIS GeoPackages and monitoring archives in **Unity Catalog Volumes** rather than inside the Git repository:

```python
%python
import json

paths = {
    "data_dir": "/Volumes/main_catalog/water_monitoring/data",
    "gis_dir":  "/Volumes/main_catalog/water_monitoring/gis",
    "output_dir": "/Volumes/main_catalog/water_monitoring/output"
}

with open("configs/directories.json", "w") as f:
    json.dump(paths, f, indent=2)
```

*(Alternatively, set `S3N_DATA_DIR`, `S3N_GIS_DIR`, `S3N_PROCESSED_DIR`, and `S3N_OUTPUT_DIR` as environment variables without modifying `configs/directories.json`.)*

---

### Cell 3: GIS & Data Preprocessing (`%python`)

```python
%python
import subprocess

res = subprocess.run(
    ["python3", "scripts/01_preprocess_orn.py", "--network=orn", "--source=rsn"],
    capture_output=True,
    text=True
)
print(res.stdout)
if res.stderr:
    print("Errors/Warnings:", res.stderr)
```

---

### Cell 4: Stream Network Modelling & Interactive Inspection (`%r`)

Fit a model via `scripts/run_s3n_determinand.R` (or source `scripts/s3n_plot_functions.R` / `scripts/s3n_stats_functions.R`) and inspect results using Databricks' native `display()`:

```sh
%sh
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt
```

```r
source("scripts/s3n_plot_functions.R")
source("scripts/s3n_stats_functions.R")

# List available model runs under output_dir/results/
list_available_model_runs()

# Load a fitted run bundle into memory
bundle <- load_model_run_bundle("output/results/det_0117_orn_rsn_trans_sqrt")

# Databricks display() provides interactive filtering, sorting, and charts
display(bundle$eng_pred)
```

---

### Cell 5: Model Comparison & Diagnostics (`%r`)

```r
source("scripts/04_compare_model_runs.R")

comp <- compare_s3n_runs(
  run_dirs = c(
    "output/results/det_0117_orn_rsn_trans_sqrt",
    "output/results/det_0180_orn_rsn"
  ),
  labels = c("Nitrate [0117]", "Orthophosphate [0180]")
)

display(comp$stat_tests)
```

---

### Cell 6: Inline Visualisation (`%python`)

```python
%python
from IPython.display import Image, display

display(Image(filename="/Volumes/main_catalog/water_monitoring/output/comparisons/det_0117_orn_rsn_trans_sqrt_vs_det_0180_orn_rsn/fig_comp1_cv_performance.png"))
display(Image(filename="/Volumes/main_catalog/water_monitoring/output/comparisons/det_0117_orn_rsn_trans_sqrt_vs_det_0180_orn_rsn/fig_comp4_side_by_side_maps.png"))
```

---

## 3. Databricks-Specific Advantages

1. **`display()` Support**: Wrapping data frames in `display(df)` in `%r` cells provides interactive search, column sorting, pagination, and quick charts without extra plotting code.
2. **Unity Catalog Governance**: Pointing `"output_dir"` to a Unity Catalog Volume (`/Volumes/...`) stores all predictions (`predictions_*.csv`), parameter tables, and comparison reports as governed cloud assets.
3. **Cluster Scaling**: For multi-period batch runs (`--by=year`, `--by=month`) or 100k+ reach networks, compute-optimised multi-core driver nodes can be used with `--omp=<cores>` while keeping the same notebook code.
