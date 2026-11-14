# Setup, Architecture, Configuration & Testing

This guide is the single reference for repository layout, environment dependencies, compiling the C++ OpenMP extension, JSON/environment configuration, and running automated tests.

---

## Table of Contents

- [1. Repository & Output Layout](#1-repository--output-layout)
- [2. Software & Hardware Requirements](#2-software--hardware-requirements)
  - [2.1 Python Environment (`01_preprocess_orn.py`, `format_custom_data.py`, `build_showcase.py`)](#21-python-environment-01_preprocess_ornpy-format_custom_datapy-build_showcasepy)
  - [2.2 R Environment](#22-r-environment)
  - [2.3 Compiled C++ / OpenMP Library (`src/s3n_fast_neighbors.cpp`)](#23-compiled-c--openmp-library-srcs3n_fast_neighborscpp)
  - [2.4 System Dependencies (Linux, macOS & Windows)](#24-system-dependencies-linux-macos--windows)
  - [2.5 Hardware Profile](#25-hardware-profile)
- [3. Obtaining & Placing the Required Data](#3-obtaining--placing-the-required-data)
  - [3.1 Dataset Downloads & Locations](#31-dataset-downloads--locations)
  - [3.2 Unzipping & Expected Directory Layout](#32-unzipping--expected-directory-layout)
  - [3.3 Verifying the Layout & Custom Locations](#33-verifying-the-layout--custom-locations)
  - [3.4 Additional River Networks & Citizen Science Data Sources](#34-additional-river-networks--citizen-science-data-sources)
- [4. Configuration Files & Environment Variables](#4-configuration-files--environment-variables)
  - [4.1 Directory Paths (`configs/directories.json`)](#41-directory-paths-configsdirectoriesjson)
  - [4.2 Determinand Metadata Catalog (`configs/determinands.json`)](#42-determinand-metadata-catalog-configsdeterminandsjson)
- [5. Running the Automated Test Suite](#5-running-the-automated-test-suite)
- [6. License](#6-license)

---

## 1. Repository & Output Layout

```text
spatial_stream_model/
├── cli_workbench/                    # Standalone offline CLI command builder (open index.html)
├── configs/
│   ├── directories.json              # Default data, gis, processed, and output roots
│   ├── determinands.json             # (Optional user file) Determinand labels & units
│   ├── networks/
│   │   └── orn.json                  # Built-in ORN v2 network profile
│   └── sources/
│       ├── rsn.json                  # Environment Agency River Surveillance Network
│       ├── ssn.json                  # Environment Agency Small Streams Network
│       └── composites.json           # Built-in composite sources (rss, all)
├── data/
│   ├── NCEA_River_Surveillance_Network_ARD_Open/
│   ├── NCEA_Small_Streams_Network_ARD_Open/
│   └── processed/
│       └── <net>_<source>/           # Standardised 5-file output from 01_preprocess_orn.py
├── gis/
│   └── ORN_v2.gpkg                   # National river network GeoPackage (~10.4 GB)
├── modelling_dashboard/              # Interactive offline case-study explorer & results dashboard
├── output/
│   ├── results/
│   │   └── det_<code>_<net>_<src>.../# Model predictions, CV tables, parameters, fig1–fig5
│   │       └── descriptives/         # Single-run report & figures from 05_descriptive_stats.R
│   ├── station_snapping/
│   │   └── <net>_<source>/           # Snapping tables & maps from 06_station_snapping_diagnostics.R
│   └── comparisons/
│       └── <comparison_slug>/        # Multi-run reports & figures from 04_compare_model_runs.R
├── scripts/
│   ├── 01_preprocess_orn.py          # Network extraction, ROS imputation, station snapping
│   ├── format_custom_data.py         # Helper to convert user CSVs into source format + JSON
│   ├── build_showcase.py             # Stages WebP/SVG figures & data.js for modelling_dashboard/
│   ├── run_s3n_determinand.R         # Main end-to-end R runner (fit + CV + kriging + plots)
│   ├── run_s3n_cokriging.R           # Multivariate Factor-LMC / colocated co-kriging runner
│   ├── 02_fit_s3n_models.R           # Standalone model fitting step (saves .rds checkpoint)
│   ├── 03_generate_plots.R           # Standalone plotting step (reads .rds checkpoint)
│   ├── 04_compare_model_runs.R       # Multi-run & cross-determinand comparison
│   ├── 05_descriptive_stats.R        # Single-run descriptive stats & gradient tests
│   ├── 06_station_snapping_diagnostics.R
│   ├── 07_compare_cokriging_runs.R   # Equal vs. source-weighted co-kriging comparison
│   ├── 08_compare_transforms.R       # Multi-transform (none/log/sqrt/boxcox) benchmark
│   ├── s3n_py_utils.py               # Shared Python utilities (WKB, ROS, config, snapping)
│   └── s3n_*.R                       # Shared R modules (config, covariates, transforms, model, multivariate, plot, stats)
├── src/
│   ├── Makevars                      # OpenMP compiler flags for Linux/macOS R CMD SHLIB
│   ├── Makevars.win                  # OpenMP compiler flags for Windows Rtools R CMD SHLIB
│   └── s3n_fast_neighbors.cpp        # OpenMP C++ Euler-tour & NNGP neighbor engine
├── tests/
│   ├── test_python_utils.py          # Python unit tests
│   └── test_r_and_cpp.R              # R modules + C++ shared library unit tests
└── docs/                             # Documentation (see docs/README.md)
```

> [!TIP]
> **Generating Commands Visually:** Open [`cli_workbench/index.html`](../cli_workbench/index.html) locally in any browser for a standalone point-and-click command builder supporting **Bash**, **Windows PowerShell**, and **Windows Command Prompt** (including an optional full `Rscript.exe` path & R-version selector when R is not on Windows `PATH`).

---

## 2. Software & Hardware Requirements

The pipeline minimises external dependencies so it runs reliably in headless Linux environments without heavy GIS or plotting stacks (`GDAL`, `geopandas`, `sf`, `SSN2`, or `ggplot2`).

### 2.1 Python Environment (`01_preprocess_orn.py`, `format_custom_data.py`, `build_showcase.py`)
- **Python Version**: **Python ≥ 3.8** (uses `statistics.NormalDist()` for ROS normal quantiles).
- **Third-Party Packages**: **None required** for core GeoPackage (`.gpkg`) and CSV preprocessing (uses standard library modules `sqlite3`, `struct`, `statistics`, `math`, `csv`, `json`, `argparse`).
  - Optional packages: `pyarrow` (or system `libgdal`) for GeoParquet `.parquet` networks, `pyproj` (or `osgeo.osr`) for CRS auto-reprojection in `format_custom_data.py`, and `Pillow` for WebP asset generation in `build_showcase.py`.

### 2.2 R Environment
- **R Version**: **R ≥ 3.5.0**.
- **Required R Packages**:
  | Package | Source | Purpose |
  | :--- | :--- | :--- |
  | **`S3N`** | GitHub (`remotes::install_github("jpierkunke/S3N")`) | Scalable Spatial Stream Network NNGP estimation (`BRISC_estimation_stream`), kriging (`BRISC_prediction_stream`), and parametric bootstrap (`BRISC_bootstrap`). *(Recommended for standard exponential fits; C++ fallback handles custom kernels/ME/missing `S3N`.)* |
  | `RANN`, `rdist`, `matrixStats`, `pbapply`, `dplyr` | CRAN | Runtime dependencies of `jpierkunke/S3N`. |
  | `jsonlite` *(optional / recommended)* | CRAN | Fast JSON parsing for `configs/*.json` and `network_metadata.json` (with built-in base-R fallback in [`scripts/s3n_config.R`](../scripts/s3n_config.R)). |
  | `data.table` *(optional / strongly recommended)* | CRAN | Multi-threaded `fread()` / `fwrite()` CSV I/O (`5–15x` faster on large sample files and ~91k-reach prediction outputs, with automatic base-R `read.csv()` / `write.csv()` fallback via `fast_read_csv()` / `fast_write_csv()` in [`scripts/s3n_config.R`](../scripts/s3n_config.R)). |
  | `splines`, `mgcv`, `rpart` | Bundled with standard R | Used when `--trend=spline\|gam\|rf\|gbm` is enabled (automatic fallback to `spline`/`linear` if unavailable). |
  | `stats`, `graphics`, `grDevices`, `utils`, `parallel` | Base R | Linear algebra, statistics, and headless PNG/PDF/SVG plotting (`type = "cairo"`). |

### 2.3 Compiled C++ / OpenMP Library (`src/s3n_fast_neighbors.cpp`)
R scripts dynamically load `src/s3n_fast_neighbors.so` (Linux/macOS) or `src/s3n_fast_neighbors.dll` (Windows via Rtools). [`src/Makevars`](../src/Makevars) and [`src/Makevars.win`](../src/Makevars.win) already configure `$(SHLIB_OPENMP_CXXFLAGS)`, so compile from the repository root with:

```bash
R CMD SHLIB src/s3n_fast_neighbors.cpp
```
*(A precompiled Linux `x86_64` `.so` is included; recompile after modifying C++ code, on Windows/macOS, or if hitting a distro ABI mismatch.)*

### 2.4 System Dependencies (Linux, macOS & Windows)
| Dependency | Fedora / RHEL | Ubuntu / Debian | Purpose |
| :--- | :--- | :--- | :--- |
| **C/C++ Compiler & OpenMP** | `gcc`, `gcc-c++`, `libgomp` | `build-essential`, `libomp-dev` | Compiling `src/s3n_fast_neighbors.cpp` and `S3N`. |
| **R Dev Headers & BLAS/LAPACK** | `R-core-devel`, `openblas-devel`, `lapack-devel` | `r-base-dev`, `libopenblas-dev`, `liblapack-dev` | `<R.h>`, LAPACK, and BLAS headers. |
| **Cairo Graphics** | `cairo`, `cairo-devel` | `libcairo2`, `libcairo2-dev` | Headless anti-aliased PNG/SVG rendering (`type = "cairo"`). |
| **SQLite3** | `sqlite-libs` | `libsqlite3-0` | Underlying C library used by Python `sqlite3` to read `.gpkg` files. |

#### macOS Notes (Intel & Apple Silicon)
- **Compiler & C++ Recompile**: Install Xcode Command Line Tools (`xcode-select --install`) and run `R CMD SHLIB src/s3n_fast_neighbors.cpp` from the repository root (replacing the precompiled Linux binary).
- **OpenMP (`--omp=N`)**: Stock Apple Clang + CRAN R compiles cleanly without OpenMP (running single-threaded). To enable multi-threaded C++ and `S3N`, install `brew install libomp` and add the corresponding `-Xpreprocessor -fopenmp` / `-lomp` flags in `~/.R/Makevars`.
- **Graphics Devices**: `PNG` (default) and `PDF` fall back automatically to native macOS Quartz/PDF devices when Cairo is unavailable. `SVG` output (`--fig-format=svg` or [`scripts/07_compare_cokriging_runs.R`](../scripts/07_compare_cokriging_runs.R)) requires XQuartz (`capabilities("cairo") == TRUE`).
- **GeoParquet & Arbitrary CRS Transforms**: Standard `.gpkg` and CSV preprocessing need no extra packages. For `.parquet`/`.geoparquet` files, install `pyarrow` (`pip install pyarrow`) or Homebrew `gdal` (`brew install gdal`, searched via `libgdal.dylib` / `/opt/homebrew/lib/libgdal.dylib`); for non-BNG EPSG CRS transforms, install `pyproj` or `brew install proj` (`cs2cs`).

#### Windows Notes
- **Rtools & C++ Recompile**: Install **Rtools** matching your R version so `R CMD SHLIB src/s3n_fast_neighbors.cpp` (and building `S3N` from GitHub) compiles `src/s3n_fast_neighbors.dll` with OpenMP support via [`src/Makevars.win`](../src/Makevars.win).
- **`Rscript.exe` on `PATH`**: CRAN's Windows installer does not add `Rscript.exe` to `PATH` automatically. Either add `C:\Program Files\R\R-<ver>\bin\x64\` to `PATH`, invoke the full `Rscript.exe` path, or use [`cli_workbench/index.html`](../cli_workbench/index.html) (which emits PowerShell and Command Prompt commands with an R-version path selector).
- **Python Auto-Preprocessing**: Install Python 3 from `python.org` (or Conda); [`resolve_python_bin()`](../scripts/s3n_config.R#L64-L78) automatically skips Windows' `WindowsApps` Store stub and falls back to `py -3`.
- **GeoParquet & Arbitrary CRS Transforms**: Standard `.gpkg` and CSV preprocessing work out of the box (SQLite3 is bundled with Python for Windows). For `.parquet`/`.geoparquet` files, install `pyarrow` (`pip install pyarrow`) or ensure OSGeo4W/Conda `gdal*.dll` is on `PATH`; for non-BNG EPSG CRS transforms, install `pyproj` (or ensure `cs2cs` is on `PATH`).

### 2.5 Hardware Profile
- **RAM**: ~1.5 GB peak memory for a full national run.
- **CPU**: Defaults to **`1` OpenMP thread** (`--omp=1` via `parse_n_omp_arg(..., default = 1L)`) for reproducible single-node runs; pass `--omp=4` (or higher) for multi-core speedups. With `--omp=4`, a full national + 3-basin ORN run with 5-fold CV and 4 figures completes in ~18–25 seconds.

---

## 3. Obtaining & Placing the Required Data

External river-network and monitoring datasets are not tracked in the Git repository and must be downloaded separately before running preprocessing (`scripts/01_preprocess_orn.py`) or modelling (`scripts/run_s3n_determinand.R` / `scripts/run_s3n_cokriging.R`).

### 3.1 Dataset Downloads & Locations

| Dataset | Status | Download Source & Required Format / Archive | Target Extracted Path | Config Definition |
| :--- | :--- | :--- | :--- | :--- |
| **OS Open Rivers Network (`ORN` v2.0)** | **Required** | [**openrivers.net/download/**](https://openrivers.net/download/) — select **version 2.0** in **GeoPackage** format (`ORN_v2_GeoPackage.zip`) | `gis/ORN_v2.gpkg` | [`configs/networks/orn.json`](../configs/networks/orn.json) |
| **River Surveillance Network (`RSN`)** | **Required** | [**Environment Agency NCEA — RSN Analysis Ready Data**](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/21763290-f395-4ca1-bcf6-884c2e7548d4?q=analysis+ready+data#access) — download `NCEA_River_Surveillance_Network_ARD_Open.zip` | `data/NCEA_River_Surveillance_Network_ARD_Open/` | [`configs/sources/rsn.json`](../configs/sources/rsn.json) |
| **Small Streams Network (`SSN`)** | **Optional (Recommended)** | [**Environment Agency NCEA — SSN Analysis Ready Data**](https://environment.data.gov.uk/natural-capital-ecosystem-assessment/search/98aca14e-3736-4329-b4b6-ae7c36fa5db1?q=small+streams+network&jry=qs&pg=1&rpp=20&srt=most_relevant#access) — download `NCEA_Small_Streams_Network_ARD_Open.zip` | `data/NCEA_Small_Streams_Network_ARD_Open/` | [`configs/sources/ssn.json`](../configs/sources/ssn.json) |

- **ORN v2.0 GeoPackage**: Make sure to select **version 2.0** in **GeoPackage** format (`ORN_v2_GeoPackage.zip`, ~623 MB compressed, extracting to a ~10.4 GB `ORN_v2.gpkg` file containing layer `ORN` in British National Grid `EPSG:27700`). Note that `configs/networks/orn.json` also lists `gis/ORN_v2_GeoPackage.zip` under `zip_fallbacks`, so `01_preprocess_orn.py` can read directly from the zip if needed, although extracting `gis/ORN_v2.gpkg` is faster across repeated runs.
- **RSN (`NCEA_River_Surveillance_Network_ARD_Open.zip`)**: Required by default for water-quality (`RSN Water Quality Samples.csv`, `RSN Water Quality Sample Points.csv`), LCMS/GCMS (`RSN LCMS GCMS Samples.csv`), and ecological survey tables (`RSN Macroinvertebrate Samples.csv`, `RSN Diatom Samples.csv`, `RSN Macrophyte Samples.csv`, `RSN River Habitat Surveys.csv`).
- **SSN (`NCEA_Small_Streams_Network_ARD_Open.zip`)**: Optional but recommended because it adds headwater stream observations (`--source=rsn+ssn`, the default combined run). If only `RSN` is downloaded, pass `--source=rsn` when running preprocessing and models.

### 3.2 Unzipping & Expected Directory Layout

Create `gis/` and the two dataset subdirectories inside `data/` at the repository root, then extract each archive so that `ORN_v2.gpkg` and the `.csv` files sit directly in their target folders (without an extra nested subfolder of the same name).

**Linux / macOS (Bash):**
```bash
mkdir -p gis data/NCEA_River_Surveillance_Network_ARD_Open data/NCEA_Small_Streams_Network_ARD_Open

# 1. Extract ORN v2.0 GeoPackage -> gis/ORN_v2.gpkg
unzip -o ~/Downloads/ORN_v2_GeoPackage.zip -d gis/

# 2. Extract RSN CSVs -> data/NCEA_River_Surveillance_Network_ARD_Open/*.csv
unzip -j -o ~/Downloads/NCEA_River_Surveillance_Network_ARD_Open.zip -d data/NCEA_River_Surveillance_Network_ARD_Open/

# 3. (Optional, recommended) Extract SSN CSVs -> data/NCEA_Small_Streams_Network_ARD_Open/*.csv
unzip -j -o ~/Downloads/NCEA_Small_Streams_Network_ARD_Open.zip -d data/NCEA_Small_Streams_Network_ARD_Open/
```
*(In the commands above, `-j` ignores any internal wrapper directory inside the NCEA zip archives so the CSV files land directly inside the target directory.)*

**Windows (PowerShell):**
```powershell
New-Item -ItemType Directory -Force -Path gis, data\NCEA_River_Surveillance_Network_ARD_Open, data\NCEA_Small_Streams_Network_ARD_Open | Out-Null

# 1. Extract ORN v2.0 GeoPackage -> gis\ORN_v2.gpkg
Expand-Archive -Force -Path "$HOME\Downloads\ORN_v2_GeoPackage.zip" -DestinationPath "gis"

# 2. Extract RSN archive
Expand-Archive -Force -Path "$HOME\Downloads\NCEA_River_Surveillance_Network_ARD_Open.zip" -DestinationPath "data\NCEA_River_Surveillance_Network_ARD_Open"

# 3. (Optional, recommended) Extract SSN archive
Expand-Archive -Force -Path "$HOME\Downloads\NCEA_Small_Streams_Network_ARD_Open.zip" -DestinationPath "data\NCEA_Small_Streams_Network_ARD_Open"
```

> [!IMPORTANT]
> If Windows Explorer ("Extract All...") or macOS Archive Utility creates a double-nested folder such as `data/NCEA_River_Surveillance_Network_ARD_Open/NCEA_River_Surveillance_Network_ARD_Open/*.csv`, move the `.csv` files up one level so they sit directly under `data/NCEA_River_Surveillance_Network_ARD_Open/` (and likewise for `data/NCEA_Small_Streams_Network_ARD_Open/`).

Expected layout after extraction:

```text
spatial_stream_model/
├── gis/
│   └── ORN_v2.gpkg
└── data/
    ├── NCEA_River_Surveillance_Network_ARD_Open/
    │   ├── RSN Water Quality Sample Points.csv
    │   ├── RSN Water Quality Samples.csv
    │   ├── RSN Ecology Sample Points.csv
    │   ├── RSN Macroinvertebrate Samples.csv
    │   ├── RSN Diatom Samples.csv
    │   ├── RSN Macrophyte Samples.csv
    │   ├── RSN LCMS GCMS Samples.csv
    │   └── ...
    └── NCEA_Small_Streams_Network_ARD_Open/        # optional, recommended
        ├── SSN Water Quality Sample Points.csv
        ├── SSN Water Quality Samples.csv
        ├── SSN Ecology Sample Points.csv
        ├── SSN Macroinvertebrate Samples.csv
        ├── SSN Diatom Samples.csv
        ├── SSN Macrophyte Samples.csv
        └── ...
```

### 3.3 Verifying the Layout & Custom Locations

Once extracted, verify that Python discovers the expected network and monitoring sources before running a full model:

```bash
python3 scripts/01_preprocess_orn.py --list-sources
```

- **Alternate folders**: If `gis/` or `data/` live outside the repository root (for example, on a shared drive or cloud volume), update [`configs/directories.json`](../configs/directories.json), set `S3N_GIS_DIR` / `S3N_DATA_DIR`, or pass `--gis-dir=` / `--data-dir=` (see [Section 4.1](#41-directory-paths-configsdirectoriesjson)).
- **Custom monitoring or river networks**: To ingest non-EA monitoring tables or a non-UK hydrography network instead of `ORN`/`RSN`/`SSN`, see [`./custom_data_sources_guide.md`](./custom_data_sources_guide.md) and [`./custom_water_network_guide.md`](./custom_water_network_guide.md).

### 3.4 Additional River Networks & Citizen Science Data Sources

Beyond the core `ORN` + `RSN`/`SSN` setup, the repository supports additional external river hydrography and citizen science monitoring datasets through custom network profiles (`configs/networks/*.json`), `scripts/format_custom_data.py`, and multi-source co-kriging:

| Dataset | Type | Source URL | Integration Path & Notes |
| :--- | :--- | :--- | :--- |
| **Copernicus EU-Hydro River Network Database** | Alternative / Pan-European River Network | [**land.copernicus.eu — EU-Hydro River Network Database**](https://land.copernicus.eu/en/products/eu-hydro/eu-hydro-river-network-database#download) | GeoPackage (`.gpkg`, e.g. `gis/EU_hydro_gpkg_eu/`) using `"topology": "pointer"` (`OBJECT_ID` → `NEXTDOWNID`) across connected layers (`River_Net_l`, `Canals_l`, `Ditches_l`) with `MultiLineStringZM` elevation/slope extraction. Configure via `configs/networks/<id>.json` and `--network=<id>` (increase `--max-snap-dist` to `500`–`1000` m if needed for coarser 1:100k+ centrelines; see [`./custom_water_network_guide.md`](./custom_water_network_guide.md)). |
| **Earthwatch FreshWater Watch / WaterBlitz** | Citizen Science Water Quality | [**freshwaterwatch.org/pages/explore-our-data**](https://www.freshwaterwatch.org/pages/explore-our-data) | Public nitrate and phosphate citizen science surveys, often reported as range/interval brackets (e.g. `0.2 - 0.5`, `<0.2`, `>10`). Format exported CSVs with `scripts/format_custom_data.py --range-mode=tobit` (producing e.g. `data/custom/waterblitz_tobit/` and [`configs/sources/waterblitz_tobit.json`](../configs/sources/waterblitz_tobit.json)), then fuse with statutory `RSN`/`SSN` runs via `--source=rsn+waterblitz_tobit`, `--source-weight=`, `--source-se=`, or Factor-LMC co-kriging (see [`./multi_source_and_citizen_science_guide.md`](./multi_source_and_citizen_science_guide.md)). |
| **CaSTCo Citizen Science Portal** | UK Catchment Citizen Science Data | [**data.castco.org**](https://data.castco.org/) | Catchment Systems Thinking Cooperative (CaSTCo) community water-quality and ecological monitoring datasets. Export site coordinates and sample tables to CSV, standardise into `data/custom/<source>/` and `configs/sources/<source>.json` using `scripts/format_custom_data.py`, and combine with `RSN`/`SSN` via `--source=rsn+<source>` (see [`./custom_data_sources_guide.md`](./custom_data_sources_guide.md) and [`./multi_source_and_citizen_science_guide.md`](./multi_source_and_citizen_science_guide.md)). |

---

## 4. Configuration Files & Environment Variables

### 4.1 Directory Paths (`configs/directories.json`)
All Python and R scripts resolve root paths in this precedence order (highest to lowest):
1. **CLI flags**: `--gis-dir=`, `--data-dir=`, `--processed-dir=`, `--output-dir=`
2. **Environment variables**: `S3N_GIS_DIR`, `S3N_DATA_DIR`, `S3N_PROCESSED_DIR`, `S3N_OUTPUT_DIR`
3. **`configs/directories.json`**:
   ```json
   {
     "gis_dir": "gis",
     "data_dir": "data",
     "processed_dir": "data/processed",
     "output_dir": "output"
   }
   ```

### 4.2 Determinand Metadata Catalog (`configs/determinands.json`)
`load_determinand_catalog()` in [`scripts/s3n_config.R`](../scripts/s3n_config.R) and [`scripts/s3n_py_utils.py`](../scripts/s3n_py_utils.py) reads `configs/determinands.json` to map determinand codes (e.g., `"0117"`, `"whpt_aspt"`) to human-readable plot titles, short labels, and units.

On first run, [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) automatically generates `configs/determinands.json` by scanning available RSN/SSN Water Quality and LCMS/GCMS CSVs alongside the built-in ecological index definitions (preserving any manual edits on later runs; pass `--update-determinands` to rescan for new codes):

```json
{
  "0117": {
    "short": "Nitrate",
    "name": "Nitrate as N",
    "unit": "mg/L"
  },
  "0180": {
    "short": "Orthophosphate",
    "name": "Orthophosphate as P",
    "unit": "mg/L"
  },
  "whpt_aspt": {
    "short": "WHPT ASPT",
    "name": "WHPT Average Score Per Taxon",
    "unit": "ASPT"
  }
}
```

---

## 5. Running the Automated Test Suite

The `tests/` directory contains fast unit tests covering Python preprocessing/formatting utilities, R configuration/covariate/transform helpers, and the compiled C++ OpenMP neighbour engine:

```bash
# 1. Python unit tests (WKB decoding, ROS imputation, source tags, custom data formatting)
python3 -m unittest tests/test_python_utils.py

# 2. R & C++ unit tests (requires src/s3n_fast_neighbors.so and configs/determinands.json)
Rscript tests/test_r_and_cpp.R
```

> [!NOTE]
> Both test suites check determinand metadata resolution (e.g., `"0117"` → `"Nitrate"`, `"whpt_aspt"` → `"ASPT"`), so on a fresh clone run `python3 scripts/01_preprocess_orn.py --list-sources` once if `configs/determinands.json` has not been generated yet.

---

## 6. License

This repository's source code is licensed under the **GNU General Public License v3.0 (`GPL-3.0`)** — see [`LICENSE`](../LICENSE). External river networks and monitoring datasets placed under `gis/` and `data/` are downloaded separately from their upstream providers and remain governed by their own licences and attribution terms.
