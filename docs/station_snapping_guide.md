# Station Snapping Diagnostics Guide

[`scripts/06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R) inspects how monitoring stations snap onto river network reaches after running [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py). It exports exact snapping coordinates, catchment-level distance summaries, and diagnostic maps so you can verify spatial alignment before fitting S3N models.

---

## Table of Contents

- [1. CLI & R Usage](#1-cli--r-usage)
  - [Command Line](#command-line)
  - [Supported CLI Flags](#supported-cli-flags)
  - [Calling from R](#calling-from-r)
- [2. Generated Outputs](#2-generated-outputs)
- [3. Key Columns in the Station Snapping Output](#3-key-columns-in-the-station-snapping-output)
- [4. Interpreting the Diagnostic Figures](#4-interpreting-the-diagnostic-figures)
- [5. How 2 km Spatial Hash Grid Snapping Handles Bucket Borders](#5-how-2-km-spatial-hash-grid-snapping-handles-bucket-borders)

---

## 1. CLI & R Usage

### Command Line
```bash
# List available catchments / sources
Rscript scripts/06_station_snapping_diagnostics.R --list-catchments
Rscript scripts/06_station_snapping_diagnostics.R --list-sources

# Run diagnostics for the default network (orn) and default source (rsn)
Rscript scripts/06_station_snapping_diagnostics.R

# Target specific catchments by name or CatchID (comma-separated)
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Thames
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Humber,Thames,Severn,Eden

# Target a specific network, monitoring source, or hide > 250 m excluded stations
Rscript scripts/06_station_snapping_diagnostics.R \
  --network=orn --source=rsn --catchment=Thames --only-included
```

### Supported CLI Flags
| Flag | Default | Description |
| :--- | :--- | :--- |
| `--network=<id>` | `orn` | Network ID from `configs/networks/<id>.json`. |
| `--source=<id>` | `rsn` | Monitoring source ID or composite (`rsn`, `ssn`, `rss`, `a+b`). |
| `--range-mode=<mode>` | `auto` | Resolve interval-tagged custom source (`auto`, `tobit`, `geometric`, `midpoint`, `none`). |
| `--catchment=<list>` | Top catchments | Comma-separated catchment names, `CatchID`s, or `'all'` (aliases `--catchments=`, `--basin=`, `--basins=`). |
| `--list-catchments` | `FALSE` | Print available catchments and station counts, then exit (aliases `--list`, `--list-basins`). |
| `--list-sources` | `FALSE` | Print registered monitoring sources and exit. |
| `--only-included` | `FALSE` | Exclude stations snapped beyond `250 m` (`included_250m == 0`) from diagnostic tables and plots (by default, excluded stations are shown). |
| `--out=<dir>` | Auto | Override output folder (default: `output/station_snapping/<net>_<source>/`). |
| `--fig-format=`, `--dpi=`, `--fig-size=` | `png` / `150–200` / `100` | Control exported plot format (`png`, `pdf`, `svg`), resolution, and percentage scaling. |
| `--processed-dir=`, `--data-dir=`, `--gis-dir=`, `--output-dir=` | From `configs/directories.json` | Override workspace root directories. |

### Calling from R
```r
source("scripts/06_station_snapping_diagnostics.R")
res <- analyze_station_snapping(catchment = "Thames", network = "orn", source = "rsn")
```

---

## 2. Generated Outputs

| File | Description |
| :--- | :--- |
| `station_snapping_all.csv` | Full station-to-reach snapping table across the entire network. |
| `catchment_snapping_summary.csv` | Per-catchment station counts (`Stations_Total`, `Stations_Included`) and `snap_dist_m` distribution (`Min`, `Q25`, `Median`, `Mean`, `Q75`, `P95`, `Max`). |
| `station_snapping_catchment_<CatchID>_<slug>.csv` | Station snapping rows filtered to each analysed catchment. |
| `station_snapping_report.md` | Automated Markdown summary of network-wide and catchment-level snapping statistics. |
| `fig_snap_overview_<net>.png` | 4-panel network-wide snapping diagnostic figure. |
| `fig_snap_catchment_<CatchID>_<slug>.png` | Per-catchment map, snapping distance histogram, and 4 metre-scale zoomed insets. |

---

## 3. Key Columns in the Station Snapping Output

| Column | Meaning |
| :--- | :--- |
| `site_id`, `wims_id`, `wims_region_code` | Monitoring station identifiers. |
| `snapped_RivID`, `CatchID`, `Catchment_Name` | Assigned river reach ID and hydrological catchment. |
| `station_X`, `station_Y` | Original monitoring station coordinates (British National Grid, `EPSG:27700` m). |
| `snap_reach_X`, `snap_reach_Y` | Orthogonal projected snap point on the target reach polyline. |
| `reach_mid_X`, `reach_mid_Y` | Midpoint coordinates of the snapped reach. |
| `reach_start_X/Y`, `reach_end_X/Y` | Upstream and downstream endpoint coordinates of the snapped reach. |
| `snap_dist_m` | Euclidean distance (metres) from the station to the snapped reach. |
| `included_250m` | `1` if `snap_dist_m <= 250 m` (used in spatial modelling); `0` if excluded. |

---

## 4. Interpreting the Diagnostic Figures

- **`fig_snap_overview_<net>.png`**: Summarises network-wide snapping distances (histogram + CDF against the `250 m` cutoff), spatial distribution of snapped vs. excluded stations, and top catchments by station count / median snap distance.
- **`fig_snap_catchment_<CatchID>_<slug>.png`**:
  - **Panel A**: Full catchment river network highlighting snapped reaches and station positions.
  - **Panel B**: Catchment snapping distance histogram (`snap_dist_m`) with median and `250 m` threshold markers.
  - **Panel C**: Four metre-scale zoomed insets showing the exact monitoring point (`station_X/Y`), the reach polyline (`start -> mid -> end`), and the orthogonal snapping vector to `snap_reach_X/Y`—useful for verifying whether large snap distances are caused by canal/tributary ambiguity or GPS coordinate offsets.

---

## 5. How 2 km Spatial Hash Grid Snapping Handles Bucket Borders

Station-to-reach snapping in [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py) uses a `2 km × 2 km` spatial hash grid (`site_grid`) strictly as a coarse spatial filter, **not** as a single-bucket lookup. Stations located directly on or near a grid cell boundary are handled safely by expanding every reach's query window across neighbouring cells:

1. **Indexing stations into `2 km × 2 km` cells (`site_grid`)**:
   Each monitoring station `(px, py)` (transformed into the network CRS) is placed into a grid bucket `(int(px // 2000), int(py // 2000))`.
2. **Padding each reach's bounding box before querying**:
   When scanning a river reach with bounding box `[minx, maxx] × [miny, maxy]`, the script expands the box outward on all four sides by:
   ```python
   snap_pad = max(1500.0, max_snap_m + 500.0)
   gx0, gx1 = int((minx - snap_pad) // 2000), int((maxx + snap_pad) // 2000)
   gy0, gy1 = int((miny - snap_pad) // 2000), int((maxy + snap_pad) // 2000)
   ```
   It then collects candidate stations from **every** grid cell `(gx, gy)` in `[gx0..gx1] × [gy0..gy1]` and computes the exact orthogonal Euclidean distance to the reach polyline via `point_to_polyline(px, py, poly)`.
3. **Why border stations are never missed**:
   Because `snap_pad` is at least `1,500 m` (and always strictly greater than `max_snap_m`, which defaults to `250 m`), any reach within `max_snap_m` of a station will have its padded bounding box cross into that station's `2 km` cell—even if the station and reach lie on opposite sides of a grid border. The grid only avoids testing every station in England against every reach; the final snapping decision and `included_250m` filter (`d <= max_snap_m`) use exact Euclidean polyline distance.
4. **Secondary `2 km` grid lookup (`query_orn_spatial_grid`)**:
   For nearest-reach covariate/RBD lookup (`grid[(int(mx // 2000), int(my // 2000))]`), the script similarly avoids single-cell boundary artifacts by searching a `3 × 3` cell neighbourhood (`radius = 1`, guaranteeing at least `2 km` clearance from any cell edge) first, and expanding to a `7 × 7` neighbourhood (`radius = 3`, covering at least `6 km`) if no reach is found.
