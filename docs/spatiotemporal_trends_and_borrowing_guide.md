# Spatio-Temporal Trend Mapping & Dynamic Borrowing Across Time

This guide explains how to run temporal batch analyses (`--by=year|season|month`), how **Dynamic Temporal Borrowing (`--st-borrow`)** shares statistical strength across time periods, and how **Reach-Level Trend Analysis** identifies statistically significant water-quality improvements or degradation across the entire river network.

---

## 1. Quick Start & CLI Flags

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--by=<mode>` | *(Off)* | Splits observations into ordered time slices and fits the full S3N model per slice:<br>• `year`: Every calendar year present in the dataset (`year_2021/`, `year_2022/`, ...)<br>• `season`: All 4 meteorological seasons (`winter`, `spring`, `summer`, `autumn`)<br>• `month`: Calendar months `01`–`12` (`month_01/` .. `month_12/`). |
| `--st-borrow` *(alias `--temporal-borrow`)* | `FALSE` | Enables **Dynamic Spatio-Temporal Borrowing** across slices when `--by=` is active. Appends `_st` to the output folder name, updates `predictions_*.csv` (and `exceed_prob` / `compliance_tier` if `--threshold=` is set), and reduces kriging SEs at sparsely sampled reaches. |
| `--year=<YYYY>` / `--season=<s>` / `--month=<M>` | *(All)* | Optional pre-filter before `--by=` slicing (e.g., `--year=2023 --by=month` analyses all 12 months of 2023). |

---

## 2. Practical Usage Examples

### 2.1 Multi-Year Trend Mapping (Standard vs. Dynamic Borrowing)
```bash
# Fit every year independently and generate Figure 5 (multi-year grid) + Figure 6 (trend map)
Rscript scripts/run_s3n_determinand.R 0117 --by=year

# Same multi-year run WITH dynamic borrowing across years (recommended when station coverage varies year-to-year)
Rscript scripts/run_s3n_determinand.R 0117 --by=year --st-borrow
```

### 2.2 Combining with Regulatory Exceedance Mapping (`--threshold=uk`)
When `--st-borrow` and `--threshold=` are used together, exceedance probabilities (`exceed_prob`) and 5-tier management classifications (`compliance_tier`) are automatically recomputed from the temporally borrowed mean and reduced kriging SE:
```bash
Rscript scripts/run_s3n_determinand.R 0117 --by=year --st-borrow --threshold=uk
```

### 2.3 Seasonal & Monthly Dynamics
```bash
# Compare all 4 meteorological seasons across England with unified scale bars
Rscript scripts/run_s3n_determinand.R 0180 --by=season --st-borrow

# Compare all 12 months within a single year (2023)
Rscript scripts/run_s3n_determinand.R 0117 --year=2023 --by=month
```

### 2.4 Combining with Custom Basins, Non-Linear Trends, or Normalisation
```bash
# Multi-year River Eden + All England Nitrate trends using GAM + S3N and Z-score normalisation
Rscript scripts/run_s3n_determinand.R 0117 --by=year --st-borrow --trend=gam --basin=Eden --norm=zscore

# Heavily censored Dissolved Lead trends using Spatial Tobit EM kriging
Rscript scripts/run_s3n_determinand.R 0052 --by=year --st-borrow --censor=em
```

---

## 3. How Dynamic Borrowing Across Time Works (`--st-borrow`)

### 3.1 Why Borrow Strength Across Time?
In real monitoring networks, station coverage is rarely identical every year:
- A tributary might be sampled in 2021 and 2023, but **skipped in 2022** due to drought, access issues, or rotating WFD surveillance schedules.
- If 2022 is fitted completely independently, that unmonitored tributary falls back toward the fixed-effect regression mean with a **large kriging standard error (`se_model_sc`)**.
- Yet we know from 2021 and 2023 what that tributary's spatial residual looks like.

### 3.2 Two-Step Empirical Bayes Temporal Smoother
After all $T$ slices ($t = 1, \dots, T$) are fitted by S3N:

1. **Data-Driven Temporal Correlation ($\hat{\rho}_t$)**:
   The pipeline identifies reaches that were directly monitored in consecutive slices $(t, t+1)$ and computes the average Pearson correlation of `obs_model_sc` across adjacent periods, clamped to $\hat{\rho}_t \in [0.15, 0.85]$ (defaulting to `0.60` if fewer than 10 overlapping stations exist).
   - High $\hat{\rho}_t$ (e.g. `0.80`) means water quality is persistent year-to-year, so neighboring years carry strong information.
   - Low $\hat{\rho}_t$ (e.g. `0.25`) means conditions fluctuate rapidly, keeping borrowing conservative.

2. **Precision-Weighted Exponential Smoothing**:
   Let $\hat{z}_{i,s}$ and $v_{i,s} = \text{SE}_{i,s}^2$ be the independent S3N kriging prediction and variance at reach $i$ in period $s$. For target period $t$, each period $s$ receives temporal weight $w_{t,s} = \hat{\rho}_t^{|t - s|}$ multiplied by its local spatial precision $1 / v_{i,s}$:

   $$
   \hat{z}_{i,t}^{\text{ST}} = \frac{\sum_{s=1}^T \hat{\rho}_t^{|t - s|} \, \hat{z}_{i,s} / v_{i,s}}{\sum_{s=1}^T \hat{\rho}_t^{|t - s|} / v_{i,s}}, \qquad \text{SE}_{i,t}^{\text{ST}} = \left(\sum_{s=1}^T \frac{\hat{\rho}_t^{|t - s|}}{v_{i,s}}\right)^{-1/2}
   $$

   **Key intuition**:
   - **Densely monitored reach in period $t$**: $v_{i,t}$ is small, so own-period precision $1/v_{i,t}$ dominates and $\hat{z}_{i,t}^{\text{ST}} \approx \hat{z}_{i,t}$.
   - **Unmonitored reach in period $t$ (but monitored in $t \pm 1$)**: $v_{i,t}$ is large while $v_{i,t\pm 1}$ is small, so $\hat{z}_{i,t}^{\text{ST}}$ borrows smoothly from adjacent periods and $\text{SE}_{i,t}^{\text{ST}} < \text{SE}_{i,t}$.

3. **Rewriting Outputs**:
   Updated `pred_model_sc`, `se_model_sc`, `pred_combined` (back-transformed to raw units), and (if `--threshold=` is active) `exceed_prob` and `compliance_tier` are saved back into each slice's `predictions_*.csv` before generating `fig1`–`fig6`.

---

## 4. How Reach-Level Trend Analysis Works

Whenever `--by=` produces at least $T \ge 2$ valid slices, `compute_spatiotemporal_trends()` automatically computes a reach-by-reach temporal trend across the entire network—whether `--st-borrow` is enabled or not.

### 4.1 Precision-Weighted Least Squares (WLS) Slope
For each reach $i$, we regress the predicted model-scale values $\hat{z}_{i,1}, \dots, \hat{z}_{i,T}$ against the integer period index $t = 1, \dots, T$ using inverse-variance weights $w_{i,t} = 1 / \max(\text{SE}_{i,t}^2, 10^{-8})$:

$$
\hat{\beta}_{\text{time}}(i) = \frac{\sum_{t=1}^T w_{i,t} (t - \bar{t}_i)(\hat{z}_{i,t} - \bar{z}_i)}{\sum_{t=1}^T w_{i,t} (t - \bar{t}_i)^2}, \qquad \text{SE}\big(\hat{\beta}_{\text{time}}(i)\big) = \left(\sum_{t=1}^T w_{i,t} (t - \bar{t}_i)^2\right)^{-1/2}
$$

Why precision weighting? Reaches far from monitoring stations have large $\text{SE}_{i,t}$ (high uncertainty), so they receive lower weight and **smaller $Z$-scores**—preventing noisy extrapolation in unmonitored headwaters from falsely lighting up as "significant trends".

### 4.2 Wald $Z$-Test & 5-Tier Significance Classification
From the slope and its standard error:

$$
Z_i = \frac{\hat{\beta}_{\text{time}}(i)}{\text{SE}\big(\hat{\beta}_{\text{time}}(i)\big)}, \qquad p_i = 2\,\Phi(-|Z_i|)
$$

Every reach is classified into one of five `trend_tier` categories:

| `trend_tier` | Condition | Map Colour (`fig6`) |
| :--- | :--- | :--- |
| `Significant Decrease (p < 0.05)` | $\hat{\beta}_{\text{time}} < 0$ and $p < 0.05$ | Dark Blue (`#2166ac`) |
| `Moderate Decrease (p < 0.20)` | $\hat{\beta}_{\text{time}} < 0$ and $0.05 \le p < 0.20$ | Light Blue (`#67a9cf`) |
| `Stable / No Clear Trend` | $p \ge 0.20$ | Light Grey (`#cfd4da`) |
| `Moderate Increase (p < 0.20)` | $\hat{\beta}_{\text{time}} > 0$ and $0.05 \le p < 0.20$ | Light Orange (`#ef8a62`) |
| `Significant Increase (p < 0.05)` | $\hat{\beta}_{\text{time}} > 0$ and $p < 0.05$ | Dark Red (`#b2182b`) |

The pipeline also reports `delta_raw` $= g^{-1}(\hat{z}_{i,T}) - g^{-1}(\hat{z}_{i,1})$, giving the net estimated change in physical units (e.g. `mg/L`) between the first and last period.

---

## 5. Output Files, Columns & Figures

Under `output/results/det_<code>_<net>_<src>[_st]/`:

| File | Description |
| :--- | :--- |
| `<period>/predictions_*.csv` | Per-period reach predictions (updated with borrowed estimates when `--st-borrow` is active). |
| `cv_summary_by_<by>_<code>.csv` | Combined cross-validation metrics across all periods. |
| `parameters_by_<by>_<code>.csv` | Combined fixed-effect and spatial covariance parameters across all periods. |
| **`spatiotemporal_trends_by_<by>_<code>.csv`** | Reach-level overall trend table with `RivID`, `CatchID`, `Shreve`, coordinates, `trend_slope_sc`, `trend_se_sc`, `trend_z`, `trend_p`, `delta_raw`, and `trend_tier`. |
| **`spatiotemporal_transitions_by_<by>_<code>.csv`** | Reach-level pairwise period-to-period transitions (`Year 2021 -> Year 2022`, `Winter -> Spring -> Summer -> Autumn -> Winter`, `Month 01 -> Month 02`, ...) with `delta_sc`, `delta_se_sc`, `trend_z`, `trend_p`, `delta_raw`, and `trend_tier`. |
| **`fig5_england_by_<by>_<code>_map.png`** | Multi-panel All England comparison grid sharing a unified colour scale across all periods. |
| **`fig6_spatiotemporal_trends_by_<by>_<code>.png`** | Two-panel All England **overall trend** map for `--by=year`, `--by=season`, and `--by=month` (**Panel A**: Temporal slope $\hat{\beta}_{\text{time}}$; **Panel B**: 5-tier Wald $Z$-test significance map). |
| **`fig7_spatiotemporal_transitions_by_<by>_<code>.png`** | Multi-panel **Period-to-Period Transition Grid** (year-to-year, season-to-season, or month-to-month): **Top Row** of each block shows model-scale step change $\Delta z$; **Bottom Row** shows 5-tier Wald $Z$-test significance. |

---

## 6. How to Interpret Figure 6 & Figure 7

### Figure 6 (`fig6_spatiotemporal_trends_by_<by>_<code>.png`) — Overall Trend
`fig6` always shows the **overall linear trend** across all periods in your `--by=` run at every river reach (both monitored and unmonitored).

- **Panel A (Trend Slope $\hat{\beta}_{\text{time}}$)**: The precision-weighted linear slope on the **model scale** (e.g. `log` scale for most water quality determinands, or raw units when `--transform=none`). **Blue** = decreasing, **White** = near-zero change, **Red** = increasing. On a `log` scale, `+0.05` per period is roughly **+5% per period**. Check `delta_raw` in `spatiotemporal_trends_by_<by>_<code>.csv` for exact physical-unit net change.
- **Panel B (Trend Direction & Significance — Wald $Z$-Test)**: Filters the slope through each reach's kriging uncertainty ($\text{SE}_{i,t}$). **Dark Blue / Dark Red (`p < 0.05`)** highlights statistically confident decreases or increases; **Light Blue / Light Orange (`0.05 <= p < 0.20`)** indicates moderate evidence; **Grey (`p >= 0.20`)** means either truly stable or unmonitored reaches with wide kriging SEs.

### Figure 7 (`fig7_spatiotemporal_transitions_by_<by>_<code>.png`) — Period-to-Period Transitions
While `fig6` summarises the overall net trajectory, **`fig7`** breaks down every consecutive step (**year-to-year**, **season-to-season**, or **month-to-month**):
- **For `--by=year`**: Consecutive transitions (`Year 2021 -> Year 2022`, `Year 2022 -> Year 2023`, ...), showing exactly *which* years drove a sudden jump or improvement.
- **For `--by=season`**: Cyclic transitions (`Winter -> Spring`, `Spring -> Summer`, `Summer -> Autumn`, `Autumn -> Winter`), capturing spring flushing, summer uptake/dilution, and autumn re-wetting.
- **For `--by=month`**: Month-to-month transitions (`Month 01 -> Month 02`, ..., wrapping `Month 12 -> Month 01` when all 12 months are present), arranged in blocks of up to 4 transitions per row pair so every panel stays legible.
- Each transition has two stacked panels: **Step Change ($\Delta z = z_{\text{to}} - z_{\text{from}}$)** on a shared diverging blue–white–red scale above its **5-Tier Wald $Z$-Test Significance** map.
