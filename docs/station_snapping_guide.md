# Station Snapping Diagnostics Guide

[`scripts/06_station_snapping_diagnostics.R`](../scripts/06_station_snapping_diagnostics.R) inspects how monitoring stations snap onto river network reaches after running [`scripts/01_preprocess_orn.py`](../scripts/01_preprocess_orn.py). It exports exact snapping coordinates, catchment-level distance summaries, and diagnostic maps so you can verify spatial alignment before fitting S3N models.

---

## 1. CLI & R Usage

### Command Line
```bash
# List available catchments in the preprocessed dataset
Rscript scripts/06_station_snapping_diagnostics.R --list-catchments

# Run diagnostics for the default network (orn) and default source (rsn)
Rscript scripts/06_station_snapping_diagnostics.R

# Target specific catchments by name or CatchID (comma-separated)
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Thames
Rscript scripts/06_station_snapping_diagnostics.R --catchment=Humber,Thames,Severn,Eden

# Target a specific network, monitoring source, or snapping threshold
Rscript scripts/06_station_snapping_diagnostics.R \
  --network=orn --source=rsn --catchment=Thames --max-plots=8

# Include stations exceeding the 250 m threshold in plots/tables
Rscript scripts/06_station_snapping_diagnostics.R --include-excluded
```

### Supported CLI Flags
| Flag | Default | Description |
| :--- | :--- | :--- |
| `--network=<id>` | `orn` | Network ID from `configs/networks/<id>.json`. |
| `--source=<id>` | `rsn` | Monitoring source ID or composite (`rsn`, `ssn`, `rss`, `a+b`). |
| `--catchment=<list>` | Top catchments | Comma-separated catchment names or `CatchID`s to plot in detail. |
| `--list-catchments` | `FALSE` | Print available catchments and station counts, then exit. |
| `--max-plots=<N>` | Auto | Maximum number of per-catchment diagnostic plots to generate when `--catchment` is omitted. |
| `--include-excluded` | `FALSE` | Include stations snapped beyond `250 m` (`included_250m == 0`) in outputs. |
| `--out=<dir>` | Auto | Override output folder (default: `output/station_snapping/<net>_<source>/`). |
| `--processed-dir=`, `--data-dir=`, `--gis-dir=`, `--output-dir=` | From `configs/directories.json` | Override workspace root directories. |

### Calling from R
```r
source("scripts/06_station_snapping_diagnostics.R")
res <- analyze_station_snapping(catchment = "Thames", network = "orn", source = "rsn")
```

---

## 2. Generated Outputs (`output/station_snapping/<net>_<source>/`)

| File | Description |
| :--- | :--- |
| `station_snapping_all.csv` | Full station-to-reach snapping table across the entire network. |
| `catchment_snapping_summary.csv` | Per-catchment station counts (`Stations_Total`, `Stations_Included`) and `snap_dist_m` distribution (`Min`, `Q25`, `Median`, `Mean`, `Q75`, `P95`, `Max`). |
| `station_snapping_catchment_<CatchID>_<slug>.csv` | Station snapping rows filtered to each analyzed catchment. |
| `station_snapping_report.md` | Automated Markdown summary of network-wide and catchment-level snapping statistics. |
| `fig_snap_overview_<net>.png` | 4-panel network-wide snapping diagnostic figure. |
| `fig_snap_catchment_<CatchID>_<slug>.png` | Per-catchment map, snapping distance histogram, and 4 metre-scale zoomed insets. |

---

## 3. Key Columns in `station_snapping_all.csv`

| Column | Meaning |
| :--- | :--- |
| `site_id`, `wims_id`, `wims_region_code` | Monitoring station identifiers. |
| `snapped_RivID`, `CatchID`, `Catchment_Name` | Assigned river reach ID and hydrological catchment. |
| `station_X`, `station_Y` | Original monitoring station coordinates (British National Grid, `EPSG:27700` m). |
| `snap_reach_X`, `snap_reach_Y` | Orthogonal projected snap point on the target reach polyline. |
| `reach_mid_X`, `reach_mid_Y` | Midpoint coordinates of the snapped reach. |
| `reach_start_X/Y`, `reach_end_X/Y` | Upstream and downstream endpoint coordinates of the snapped reach. |
| `snap_dist_m` | Euclidean distance (metres) from the station to the snapped reach. |
| `included_250m` | `1` if `snap_dist_m <= 250 m` (used in S3N modelling); `0` if excluded. |

---

## 4. Interpreting the Diagnostic Figures

- **`fig_snap_overview_<net>.png`**: Summarizes network-wide snapping distances (histogram + CDF against the `250 m` cutoff), spatial distribution of snapped vs. excluded stations, and top catchments by station count / median snap distance.
- **`fig_snap_catchment_<CatchID>_<slug>.png`**:
  - **Panel A**: Full catchment river network highlighting snapped reaches and station positions.
  - **Panel B**: Catchment snapping distance histogram (`snap_dist_m`) with median and `250 m` threshold markers.
  - **Panel C**: Four metre-scale zoomed insets showing the exact monitoring point (`station_X/Y`), the reach polyline (`start -> mid -> end`), and the orthogonal snapping vector to `snap_reach_X/Y`—useful for verifying whether large snap distances are caused by canal/tributary ambiguity or GPS coordinate offsets.
