# Optimal Monitoring Network Design & Gap Mapping Guide

When expanding water quality monitoring, simply placing new stations where Kriging Standard Error (`SE`) is highest can waste visits on short, isolated headwater twigs. Passing `--design-sites=<K>` (alias `--optimal-design=<K>`) invokes the OpenMP C++ greedy active-learning solver (`C_s3n_optimal_design`) to select the $K$ unmonitored reaches that yield the largest network-wide reduction in prediction uncertainty, exporting both a ranked candidate CSV and multi-panel **Monitoring Gap & Priority Placement Maps** (`plot_fig_optimal_design()`).

---

## Table of Contents

- [1. CLI Flags & Pipeline Integration](#1-cli-flags--pipeline-integration)
- [2. Mathematical Formulation (`C_s3n_optimal_design`)](#2-mathematical-formulation-c_s3n_optimal_design)
- [3. Visualisation Function (`plot_fig_optimal_design`) & Map Panels](#3-visualisation-function-plot_fig_optimal_design--map-panels)
- [4. Output CSV Schema (`optimal_monitoring_design_<det>.csv`)](#4-output-csv-schema-optimal_monitoring_design_detcsv)
- [5. Practical Interpretation & Combining with Threshold Exceedance](#5-practical-interpretation--combining-with-threshold-exceedance)
- [6. Worked CLI Examples](#6-worked-cli-examples)

---

## 1. CLI Flags & Pipeline Integration

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--design-sites=<K>` | `0` *(Off)* | Number of new unmonitored reaches $K$ to greedily select and rank per fitted spatial domain (`All England`, `Humber`, `Thames`, `Severn`, plus any `--basin=` targets). |
| `--optimal-design=<K>` | `0` *(Off)* | Alias for `--design-sites=<K>`. |
| `--list-validation` | — | Prints all available CV modes, threshold options, and monitoring design flags. |
| `--dpi=<N>` / `--fig-format=<fmt>` | `200` / `png` | Applies to generated `fig_optimal_design_*` maps (`png`, `pdf`, or `svg`). |

### How Each Entrypoint Uses It
- [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R): Fits each spatial domain (`fit_s3n_domain()`), calls `compute_optimal_monitoring_design()`, writes `optimal_monitoring_design_<det>.csv`, and automatically renders `fig_optimal_design_<domain>_<det>.png` via `plot_fig_optimal_design()`.
- [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R): Runs `compute_optimal_monitoring_design()` when `--design-sites=<K>` is passed and exports `optimal_monitoring_design_0117.csv`.
- [`scripts/03_generate_plots.R`](../scripts/03_generate_plots.R): Automatically detects `optimal_monitoring_design_<det>.csv` in the target output directory and re-renders `fig_optimal_design_*` figures (supports standalone styling / SVG export without re-fitting).

---

## 2. Mathematical Formulation (`C_s3n_optimal_design`)

> **In plain English:** At each step, the solver tests every unmonitored reach and scores how much adding a gauge there would lower uncertainty across all connected downstream and neighbouring reaches combined. After picking the best site (`#1`), it updates the remaining uncertainty map before picking `#2`, preventing picks from bunching together on the same stream.

After universal S3N kriging produces the baseline combined prediction standard error $\text{SE}_{\text{curr}}(s_j)$ at every reach $s_j$, the solver ([`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp) / [`compute_optimal_monitoring_design()`](../scripts/s3n_model_functions.R)) models the variance reduction achieved across evaluation reaches $\mathcal{E}$ by placing a station with nugget variance $\tau^2$ at an unmonitored candidate reach $c \in \mathcal{C}$ (`has_obs == 0`).

1. **Candidate & Evaluation Sets**:
   - $\mathcal{C} = \{i : \text{has\_obs}_i = 0\}$ is the pool of unmonitored reaches.
   - To guarantee sub-second interactive turnaround on national networks ($> 25{,}000$–$90{,}000+$ reaches), when $|\mathcal{C}| > 4{,}000$ or $|\mathcal{E}| > 4{,}000$ the pools are deterministically thinned to $4{,}000$ evenly spaced river-order reaches.
2. **Stream + Euclidean Spatial Correlation**:
   For any pair of reaches $(a, b)$, `pair_corr(a, b)` evaluates the fitted hybrid covariance kernel using Euler-tour tree intervals (`tin`, `tout`) and lowest-common-ancestor traversal (`d2m_km`) in $O(1)$–$O(\text{depth})$ time:
   $$
   r(a, b) = w_{\text{stream}} \, k_\phi\big(d_{\text{stream}}(a, b)\big)\,\mathbb{I}(\text{same basin}) + w_{\text{euc}} \, k_\phi\big(d_{\text{euc}}(a, b)\big)
   $$
3. **Greedy Schur-Complement Step ($t = 1, \dots, K$)**:
   At iteration $t$, for every unselected candidate $c \in \mathcal{C}$ with residual variance $v_c = \sigma^2_{\text{curr}}(c)$:
   $$
   \Delta V(c) = \sum_{j \in \mathcal{E}} \min\!\left(v_j,\; \frac{\big(\sqrt{v_j \, v_c} \cdot r(s_j, c)\big)^2}{v_c + \tau^2}\right)
   $$
   The solver selects $c^* = \arg\max_c \Delta V(c)$ in parallel (`#pragma omp parallel for`), records the domain-wide mean SE drop, and **immediately updates both $v_j$ (`v_eval`) and $v_c$ (`v_cand`) in place**:
   - Updating `v_eval` & `v_cand` after each pick prevents redundant clustering — subsequent ranks automatically move to *other* unmonitored catchments or tributaries rather than picking adjacent reaches next to rank `#1`.

---

## 3. Visualisation Function (`plot_fig_optimal_design`) & Map Panels

[`plot_fig_optimal_design()`](../scripts/s3n_plot_functions.R) renders a 3-panel monitoring gap & active-learning layout whenever design results exist:

- **Files Written**:
  - `fig_optimal_design_england_<det>.png` — All England monitoring gaps & top-$K$ placements
  - `fig_optimal_design_humber_<det>.png`, `fig_optimal_design_thames_<det>.png`, `fig_optimal_design_severn_<det>.png` — Basin-scale active-learning designs (each core basin optimises its own local $K$ picks)
  - `fig_optimal_design_basin_<slug>_<det>.png` — Any custom catchment requested with `--basin=`

### Panel Anatomy
| Panel | Visual Elements | What It Shows |
| :--- | :--- | :--- |
| **A. Monitoring Gaps Map** *(left column)* | Stream network shaded by baseline **Kriging SE** (`uncert_cols`: deep navy = high uncertainty gap); existing monitoring stations shown as subtle teal points; **Top 1–10 priority reaches** plotted in vivid crimson with `#1..#10` callout labels; **Ranks 11–K** plotted in amber. | Exactly *where* current network coverage leaves high residual uncertainty and *which* unmonitored tributaries bridge the largest connected gaps. |
| **B. Variance Reduction Curve** *(top right)* | Crimson curve plotting **Cumulative % Variance Reduced** (`Cum_Pct_Reduced`), pale blue bars showing per-site marginal `%` gain (`Pct_Var_Reduced`), and a secondary dashed teal curve tracing **Mean Kriging SE** (`Mean_SE_Before` $\to$ `Mean_SE_After`). | Diminishing returns "elbow" — shows how many extra stations are cost-effective before marginal uncertainty reduction plateaus. |
| **C. Top-12 Priority Reaches** *(bottom right)* | Horizontal bar chart ordered `#1` to `#12` displaying `RivID`, `CatchID`, Strahler stream order (`O`), local baseline Kriging SE, and individual `% dV` reduction. | Actionable candidate table ready for field logistics and catchment-officer site reconnaissance. |

---

## 4. Output CSV Schema (`optimal_monitoring_design_<det>.csv`)

Each fitted domain contributes `K` rows ordered by `Rank` (`1..K`):

| Column | Type | Meaning |
| :--- | :---: | :--- |
| `Domain` | string | Spatial domain label (`All England Stream Network`, `Humber Basin`, etc.) |
| `Rank` | int | Sequential greedy priority (`1` = highest variance reduction first) |
| `RivID` | int | ORN v2 Reach ID of the candidate reach |
| `CatchID` | int | Hydrometric Catchment ID |
| `Strahler` | int | Strahler stream order of the candidate reach |
| `Dist2Mth_km` | float | Distance along river network to tidal mouth (km) |
| `US_Accum_km` | float | Total upstream drainage network length (km) |
| `mid_X`, `mid_Y` | float | British National Grid Easting / Northing of reach midpoint (m) |
| `Elevation` | float | Mean reach elevation (m AOD) |
| `Current_SE_sc` | float | Baseline S3N Kriging SE on the model transform scale before placement |
| `Var_Reduced` | float | Sum of squared-SE reductions $\Delta V(c^*)$ across evaluation reaches |
| `Pct_Var_Reduced` | float | Marginal `%` of initial total domain variance eliminated by this station |
| `Cum_Pct_Reduced` | float | Cumulative `%` of initial total domain variance eliminated by ranks `1..Rank` |
| `Mean_SE_Before` | float | Mean domain Kriging SE across evaluation reaches immediately before step `Rank` |
| `Mean_SE_After` | float | Mean domain Kriging SE across evaluation reaches immediately after step `Rank` |

---

## 5. Practical Interpretation & Combining with Threshold Exceedance

1. **Statistical Variance Design vs. Compliance Risk Gaps**:
   - `--design-sites=<K>` minimises **model-wide epistemic uncertainty** (`Kriging SE`), naturally favouring long, well-connected unmonitored mainstem/tributary corridors in poorly sampled catchments.
   - Threshold exceedance (`--threshold=uk`) identifies **decision uncertainty** (`Uncertain 20–80%`, `0.20 <= P <= 0.80`) where predicted levels sit close to a standard. Cross-referencing both (`--threshold=uk --design-sites=50`) highlights high-SE candidate reaches that *also* lie inside yellow `Uncertain` management zones — delivering maximum regulatory value per extra water sample.
2. **Reading the Marginal Gain Curve**:
   - Large `Pct_Var_Reduced` on early ranks (`1..10`) signals clear structural holes in current monitoring coverage; when `Pct_Var_Reduced` flattens toward `0.05%–0.1%` per site, extra stations mainly reduce local reach nugget error rather than regional interpolation error.
3. **National vs. Catchment Scope**:
   - In an `All England` fit, `K = 50` ranks national structural gaps across all `267` catchments; if your operational budget targets one catchment (e.g. River Eden or Humber), inspect the basin-specific `optimal_monitoring_design_<det>.csv` rows / `fig_optimal_design_<basin>_<det>.png` map or pass `--basin=<RiverName>`.

---

## 6. Worked CLI Examples

```bash
# 1. Rank top 50 new monitoring reaches nationally + per core basin on Nitrate (0117)
Rscript scripts/run_s3n_determinand.R 0117 --design-sites=50 --threshold=uk --dpi=300 --omp=4

# 2. Focus on a single river catchment (--basin=Eden) and rank top 15 local additions
Rscript scripts/run_s3n_determinand.R 0117 --basin=Eden --design-sites=15 --threshold=uk --dpi=300

# 3. Include environmental covariates + Leave-Catchment-Out CV + optimal design
Rscript scripts/run_s3n_determinand.R 0180 --cv-mode=catchment --threshold=uk --design-sites=30 --omp=4

# 4. Re-render an existing run's optimal design map cleanly as scalable SVG
Rscript scripts/03_generate_plots.R --dir=output/results/det_0117_orn_rsn --fig-format=svg
```
