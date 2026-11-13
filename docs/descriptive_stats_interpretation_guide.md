# S3N Descriptive Statistics & Environmental Gradient Analysis — Interpretation Guide

## Table of Contents

- [1. Overview & Purpose](#1-overview--purpose)
- [2. Quick-Start CLI Reference & Generated Outputs](#2-quick-start-cli-reference--generated-outputs)
  - [2.1 Running the Script](#21-running-the-script)
  - [2.2 Output Files Written to `<run_dir>/descriptives/`](#22-output-files-written-to-run_dirdescriptives)
- [3. Section-by-Section Guide to the Report & Metrics](#3-section-by-section-guide-to-the-report--metrics)
  - [3.1 Section 1: Comprehensive Distributional Statistics (`descriptive_stats_summary.csv`)](#31-section-1-comprehensive-distributional-statistics-descriptive_stats_summarycsv)
  - [3.2 Section 2: Monitoring Network Representativeness & Bias Analysis (`monitoring_representativeness.csv`)](#32-section-2-monitoring-network-representativeness--bias-analysis-monitoring_representativenesscsv)
  - [3.3 Section 3: Stratified Environmental, Topological & Basin Analysis (`stratified_results_summary.csv`)](#33-section-3-stratified-environmental-topological--basin-analysis-stratified_results_summarycsv)
  - [3.4 Section 4: Within-Run Statistical Hypothesis Tests (`within_run_statistical_tests.csv`)](#34-section-4-within-run-statistical-hypothesis-tests-within_run_statistical_testscsv)
  - [3.5 Section 5: Cross-Validation Performance Summary (`cv_comparison_*.csv`)](#35-section-5-cross-validation-performance-summary-cv_comparison_csv)
- [4. Visual Diagnostic Guide: Interpreting the Figures](#4-visual-diagnostic-guide-interpreting-the-figures)
  - [4.1 Figure 1: `fig_desc1_distributions_and_representativeness.png`](#41-figure-1-fig_desc1_distributions_and_representativenesspng)
  - [4.2 Figure 2: `fig_desc2_stratified_gradients.png`](#42-figure-2-fig_desc2_stratified_gradientspng)
- [5. Worked Empirical Summary: Nitrate (`0117`) on ORN v2 (`output/results/det_0117_orn_rsn/descriptives/`)](#5-worked-empirical-summary-nitrate-0117-on-orn-v2-outputresultsdet_0117_orn_rsndescriptives)
- [6. Related Documentation](#6-related-documentation)

---

## 1. Overview & Purpose

The **Single-Output Descriptive Statistics & Environmental Gradient Analysis** script ([`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R)) evaluates the full reach-level prediction table, station observations, and 5-fold cross-validation results produced by a single Spatial Stream Network (`S3N`) model run.

While standard model-fitting output focuses on parameter estimates and station cross-validation accuracy, river management and scientific interpretation require understanding what the model actually predicts across the **entire river network**—including the **75%–95% of stream reaches that have no monitoring station**. Specifically, [`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R) answers five core questions:

1. **Distributional Structure**: What is the central tendency, spread, skewness, and tail behavior of predicted concentrations or ecological scores across all river reaches?
2. **Monitoring Network Representativeness**: Are regulatory monitoring stations (such as the Environment Agency's River Surveillance Network) statistically representative of the river network as a whole, or are they biased toward specific landscape positions (e.g., accessible lowland mainstems)?
3. **Environmental & Topological Stratification**: How do predicted concentrations and prediction uncertainties vary across **topography** (Upland vs. Lowland), **hydrogeology** (Principal Aquifers vs. Non-Aquifers), **Strahler stream order** (headwaters to mainstems), and **major river basins**?
4. **Spatial Covariance & Model Component Behavior**: Does incorporating branching river-network topology (`S3N Combined`) statistically outperform standard non-spatial regression (`OLS`) or single-component (`Tail-Up` / `Tail-Down`) models at held-out stations, and how closely do flow-connected (`Tail-Up`) and hydrologic-distance (`Tail-Down`) predictions agree?
5. **Kriging Uncertainty Drivers**: Where in the river network is prediction uncertainty (`Kriging SE`) lowest and highest, and how strongly do monitoring station density and stream order control uncertainty?

---

## 2. Quick-Start CLI Reference & Generated Outputs

### 2.1 Running the Script

```bash
# List available completed model runs in output/results/
Rscript scripts/05_descriptive_stats.R --list

# Analyze a single output directory directly (full path or folder name under output/results/)
Rscript scripts/05_descriptive_stats.R output/results/det_0117_orn_rsn

# Analyze multiple outputs sequentially (each receives its own report in <run_dir>/descriptives/)
Rscript scripts/05_descriptive_stats.R \
  output/results/det_0117_orn_rsn \
  output/results/det_0180_orn_rsn \
  output/results/det_whpt_aspt_orn_rsn

# Custom workspace directories and custom output subfolder
Rscript scripts/05_descriptive_stats.R output/results/det_0117_orn_rsn \
  --label="Nitrate (ORN + RSN)" \
  --gis-dir=gis --data-dir=data --processed-dir=data/processed --output-dir=output \
  --out=output/results/det_0117_orn_rsn/descriptives
```

### 2.2 Output Files Written to `<run_dir>/descriptives/`

| Output File | Format | Description |
| :--- | :--- | :--- |
| `descriptive_report.md` | Markdown | Automated narrative report containing all summary tables, statistical test interpretations, and embedded links to figures. |
| `descriptive_stats_summary.csv` | CSV | 18-column distributional summary across 7 variable groups (observations, combined predictions, unmonitored reaches, Tail-Up, Tail-Down, model-scale predictions, and Kriging SE). |
| `monitoring_representativeness.csv` | CSV | Statistical comparison of observed monitoring stations vs. full-network reach predictions (`Pred / Obs Median Ratio`, Mann–Whitney U test, Two-Sample KS test). |
| `stratified_results_summary.csv` | CSV | Stratified descriptive statistics and environmental gradient ratios across Topography, Hydrogeology, Landscape Regime, Stream Order, and River Basin. |
| `within_run_statistical_tests.csv` | CSV | Formal within-run hypothesis tests (paired CV model accuracy tests, Tail-Up vs. Tail-Down correlation, Kriging uncertainty tests, and environmental gradient tests). |
| `fig_desc1_distributions_and_representativeness.png` | PNG (2400×1800) | 4-panel diagnostic figure covering kernel density distributions, empirical CDFs, Tail-Up vs. Tail-Down agreement, and Kriging SE by monitoring status and Strahler order. |
| `fig_desc2_stratified_gradients.png` | PNG (2400×1800) | 4-panel environmental gradient figure covering landscape/hydrogeology boxplots, Strahler order longitudinal profile, major basin distributions, and altitudinal decile profile. |

---

## 3. Section-by-Section Guide to the Report & Metrics

### 3.1 Section 1: Comprehensive Distributional Statistics (`descriptive_stats_summary.csv`)

Section 1 summarizes seven distinct subsets/components of the model output so you can compare observations, predictions, spatial components, and uncertainty side by side.

#### A. The 7 Variable Groups Explained

| Variable Group | Source Column in `predictions_*.csv` | Filter / Subset | What It Represents & Why It Matters |
| :--- | :--- | :--- | :--- |
| **`Observed Stations (Raw)`** | `obs_value` | `has_obs == 1` | True empirical values measured at monitoring stations on the natural physical scale (e.g., `mg/L` or ecological score). Serves as the empirical ground-truth benchmark. |
| **`S3N Combined Predictions (Raw)`** | `pred_combined` | All reaches (`N` total) | The primary variance-weighted `S3N` prediction across every reach in the network, back-transformed (with Taylor / Duan smearing adjustment according to `--transform`) to the natural scale. |
| **`Unmonitored Reach Predictions (Raw)`** | `pred_combined` | `has_obs == 0` | Predictions restricted strictly to reaches **without** a monitoring station. Comparing this against `Observed Stations (Raw)` reveals how unmonitored reaches (mostly headwaters) differ from monitored reaches. |
| **`S3N Tail-Up Predictions (Raw)`** | `pred_tu` | All reaches | Predictions generated solely by the **flow-connected (`Tail-Up`)** moving-average spatial model, where spatial autocorrelation only propagates along flow-connected pathways and splits at confluences using stream-order / AFV weights. |
| **`S3N Tail-Down Predictions (Raw)`** | `pred_td` | All reaches | Predictions generated solely by the **hydrologic-distance (`Tail-Down`)** spatial model, which allows spatial covariance between both flow-connected and flow-unconnected reaches sharing a common downstream junction. |
| **`S3N Combined (Model Scale)`** | `pred_model_sc` | All reaches | Predictions on the internal regression scale (set by `--transform=...`, e.g., `log(Y + offset)` or untransformed when `--transform=none`). Useful for checking Gaussianity before back-transformation. |
| **`Kriging Standard Error (SE)`** | `se_model_sc` | All reaches | Reach-level universal kriging standard error on the model scale (`√Var(ŷ(s₀) - y(s₀))`). Lower values indicate higher local confidence (near observed stations); higher values indicate extrapolation into unmonitored sub-catchments. |

#### B. How to Interpret Each Distributional Metric

Water quality concentrations (such as Nitrate and Orthophosphate) are typically right-skewed (`Log-Normal`), whereas ecological indices (`WHPT ASPT`, `RHS HQA`) are closer to symmetric (`Gaussian`). Consequently, both parametric (mean/SD) and non-parametric (median/IQR/MAD) statistics are reported:

| Metric | Mathematical Definition | How to Interpret It |
| :--- | :--- | :--- |
| **`N`** | Count of valid (non-`NA`) reaches or stations | Shows sample size. The ratio `N(Observed) / N(All Reaches)` gives the fraction of reaches directly monitored. |
| **`Mean ± SD`** | Arithmetic mean `μ` and standard deviation `σ` | Best for symmetric variables (`WHPT`, `HQA`). For right-skewed nutrients, `Mean` is pulled upward by extreme point-source or agricultural hotspots (`Mean > Median`). |
| **`CV (%)`** | Coefficient of Variation: `100 × (SD / |Mean|)` | Dimensionless measure of relative dispersion:<ul><li>`< 25%`: Low spatial variability (homogeneous across the network, e.g., `WHPT ASPT` or `HQA`).</li><li>`25% – 75%`: Moderate spatial variability.</li><li>`> 75%`: High spatial variability driven by strong landscape contrasts or point-source hotspots (typical of `Nitrate` and `Orthophosphate`).</li></ul> |
| **`Geom Mean`** | Geometric Mean: `exp(mean(log(x[x > 0])))` | Natural central tendency for log-normally distributed water quality variables. For a pure log-normal distribution, `Geom Mean ≈ Median`. |
| **`Skewness`** | Third standardized central moment: `E[(X - μ)³] / σ³` | Measures asymmetry:<ul><li>`≈ 0` (`-0.5` to `+0.5`): Approximately symmetric.</li><li>`> +1.0`: Strongly right-skewed (long tail of high-concentration reaches).</li><li>`< -1.0`: Strongly left-skewed (tail of degraded low-score reaches).</li></ul> |
| **`Kurtosis`** | Excess Kurtosis: `(E[(X - μ)⁴] / σ⁴) - 3` | Measures tail heaviness relative to a normal distribution (`0`):<ul><li>`≈ 0`: Mesokurtic (normal tails).</li><li>`> +1.0`: Leptokurtic (heavy tails / frequent extreme outlier hotspots).</li><li>`< -0.5`: Platykurtic (light tails / bounded uniform-like spread).</li></ul> |
| **`Min`, `P05`, `Q25`, `Median`, `Q75`, `P95`, `Max`** | Empirical quantiles (`0%`, `5%`, `25%`, `50%`, `75%`, `95%`, `100%`) | `Median` (`50%`) is the most robust measure of typical river condition. `P05–P95` defines the central `90%` environmental range across the network while excluding the top/bottom `5%` extremes. |
| **`IQR`** | Interquartile Range: `Q75 - Q25` | Width of the middle `50%` of reaches. Insensitive to extreme outliers. |
| **`MAD`** | Median Absolute Deviation: `1.4826 × median(|x - median(x)|)` | Robust standard-deviation equivalent. If `MAD` is much smaller than `SD` (`MAD << SD`), extreme high-leverage reaches are inflating `SD`. |

> **Why are `S3N Combined Predictions` less dispersed than `Observed Stations` (`SD(Pred) < SD(Obs)` and `Max(Pred) < Max(Obs)`)?**
> Universal kriging and regression models are **conditional expectation estimators** (`E[Y(s₀) | Data]`). By mathematical properties of best linear unbiased prediction (BLUP), kriging smooths out micro-scale white noise (`nugget variance τ²`) that occurs at individual sampling times/sites. Therefore, reach predictions always exhibit slightly lower `CV (%)`, lower `Skewness`, and fewer extreme spikes than raw station measurements.

---

### 3.2 Section 2: Monitoring Network Representativeness & Bias Analysis (`monitoring_representativeness.csv`)

Regulatory monitoring networks are rarely random spatial samples of a river network. Bridges, gauging stations, and statutory compliance points are disproportionately located on **accessible lowland, higher-order (`Strahler ≥ 3`) reaches**, whereas `75%–80%` of total stream segments in topological river networks are **first- and second-order headwaters** (many in upland moorland or forested headwaters).

Section 2 quantifies this **Monitoring Network Selection Bias** by comparing `Observed Stations (Raw)` against `S3N Combined Predictions (Raw)` across the complete network:

| Metric / Test | What It Measures | Interpretation Thresholds |
| :--- | :--- | :--- |
| **`Pred / Obs Median Ratio`** | `Median(All Reaches) / Median(Observed Stations)` | <ul><li>**`< 0.85x`**: **High-Concentration Monitoring Bias**. Monitoring stations are preferentially located in higher-concentration reaches (e.g., lowland arable/urban catchments), so full-network predictions (which include clean upland headwaters) have a substantially lower median.</li><li>**`0.85x – 1.15x`**: **Representative Median Coverage**. Monitoring stations capture the typical network median well.</li><li>**`> 1.15x`**: **Low-Concentration Monitoring Bias**. Unmonitored reaches have higher predicted values than monitored reaches.</li></ul> |
| **`Mann–Whitney U p-val`** | Two-sided non-parametric Wilcoxon rank-sum test comparing the median/rank distribution of all reaches vs. observed stations. | <ul><li>`p < 0.05` (`*`, `**`, `***`): Statistically significant location shift between the monitored station network and the full river network.</li><li>`p ≥ 0.05` (`ns`): No statistically significant difference in rank distribution.</li></ul> |
| **`KS D-Stat (p-val)`** | Two-Sample Kolmogorov–Smirnov statistic: `D = sup_x |F_pred(x) - F_obs(x)|`, measuring the maximum vertical distance between the two cumulative distribution functions (CDFs). | <ul><li>`D < 0.05`: Negligible distributional divergence.</li><li>`0.05 ≤ D < 0.15`: Moderate distributional shift (typical when headwaters differ from mainstems).</li><li>`D ≥ 0.15`: Large structural divergence between what is monitored and the full river network.</li></ul> |

> **Practical Takeaway for Reporting**: If a national regulatory agency reports the "national median river concentration" using only raw monitoring station averages (`Obs Median`), and `Pred / Obs Median Ratio = 0.70x` (as seen for Nitrate across England), the raw station summary **overestimates** the median concentration of all river reaches by `~43%` (`3.35 mg/L` vs. `2.34 mg/L`) because low-nitrate upland headwaters are under-sampled relative to their topological abundance.

---

### 3.3 Section 3: Stratified Environmental, Topological & Basin Analysis (`stratified_results_summary.csv`)

Section 3 partitions every reach in the network across **five environmental and topological stratification dimensions** to reveal where high concentrations, ecological degradation, and prediction uncertainties occur.

#### A. Definitions of the 5 Stratification Dimensions

1. **`1. Topography` (`elevation` threshold at `150 m AOD`)**:
   - **`Lowland (<150m)`**: Reaches below `150 m` elevation—typically arable/intensive pastoral agriculture, higher population density, and lower runoff dilution.
   - **`Upland (>=150m)`**: Reaches at or above `150 m` elevation—typically upland moorland, peatland, rough grazing, and high-precipitation headwaters.
2. **`2. Hydrogeology` (`aquifer` binary flag)**:
   - **`Non-Aquifer`**: Surface-runoff-dominated catchments over low-permeability bedrock or clay.
   - **`Aquifer (Chalk/Limestone/Sandstone)`**: Reaches intersecting BGS 625k Principal Aquifers (e.g., Chalk, Jurassic Limestone, Permo-Triassic Sandstone), where baseflow groundwater pathways dominate hydrochemistry and store legacy solutes (such as multi-decadal nitrate in Chalk aquifers).
3. **`3. Landscape Regime` (Combined 3-Class Typology)**:
   - **`Upland (>=150m)`**: High-elevation source waters.
   - **`Lowland Non-Aquifer`**: Low-elevation surface-water/clay catchments (often flashy runoff and point-source sensitivity).
   - **`Lowland Aquifer`**: Low-elevation groundwater-fed catchments (often highest sustained nitrate concentrations and stable thermal/flow regimes).
4. **`4. Stream Order` (Strahler Topological Hierarchy)**:
   - **`Headwaters (Order 1-2)`**: Small source streams (`~75%–80%` of all network segments).
   - **`Medium Tributaries (Order 3-4)`**: Mid-basin collector rivers (`~17%–20%` of segments).
   - **`Mainstem Rivers (Order 5+)`**: Major lowland rivers (`~3%–5%` of segments) that integrate upstream sub-catchments and major wastewater treatment works (`STWs`).
5. **`5. River Basin` (`HydroID` Topological Sub-Networks)**:
   - Compares **`All England`** against the three major benchmark catchments: **`Humber`** (mixed upland Pennines + industrial/arable lowlands), **`Thames`** (predominantly lowland Chalk/Oolite aquifer + high urban density), and **`Severn`** (Welsh/Midlands upland-to-lowland gradient).

#### B. How to Interpret Stratified Columns

| Column | Formula / Definition | Interpretation Guide |
| :--- | :--- | :--- |
| **`Reaches (N)` & `Share (%)`** | Reach count and percentage of the entire river network in that stratum. | Shows topological dominance (e.g., `Headwaters (Order 1-2)` account for `78.4%` of ORN reaches). |
| **`Obs Sites` & `Monitored (%)`** | Count of stations (`has_obs == 1`) and `100 × Obs Sites / Reaches (N)`. | Reveals **sampling density disparities**: e.g., on ORN v2 for Nitrate, only `12.6%` of `Headwaters (Order 1-2)` are monitored, compared with `57.9%` of `Medium Tributaries (Order 3-4)` and `60.2%` of `Mainstem Rivers (Order 5+)`. |
| **`Pred Mean ± SD` & `Pred Median (IQR)`** | Parametric and non-parametric central tendency and spread within the stratum. | Compare across strata to see how physical geography shifts concentrations. |
| **`Pred P95`** | 95th percentile within the stratum. | Identifies upper-tail hotspot intensity inside that landscape type. |
| **`Mean Z-Score`** | Mean of global standardized scores: `Z_i = (Pred_i - μ_global) / σ_global`. | Measures how many global standard deviations a stratum sits above (`+`) or below (`-`) the network-wide mean:<ul><li>`+0.50` to `+1.00`: Strong positive enrichment (e.g., `Lowland Aquifer` has `Mean Z = +0.58` for Nitrate).</li><li>`-0.50` to `-1.00`: Strong negative depletion (e.g., `Upland (>=150m)` has `Mean Z = -0.82` for Nitrate).</li></ul> |
| **`Kriging SE Med`** | Median universal kriging standard error (`se_model_sc`) in that stratum. | Quantifies prediction uncertainty in each stratum. Strata with low `Monitored (%)` (such as `Upland` and `Headwaters`) systematically exhibit higher `Kriging SE Med`. |
| **`Fold Ratio vs Ref`** | `Median(Stratum) / Median(Reference Stratum)` (where the first stratum in each dimension has `1.00x (ref)`). | Gives an intuitive multiplicative effect size: e.g., in `Landscape Regime` for Nitrate, using `Upland (>=150m)` (`0.843 mg/L`) as `1.00x (ref)`, `Lowland Non-Aquifer` is **`2.79x`** higher (`2.353 mg/L`) and `Lowland Aquifer` is **`4.30x`** higher (`3.621 mg/L`). |

Directly below the table in Section 3, the script prints automatic **Within-Run Environmental Gradient Tests**:
- **2-Group Strata (`Topography`, `Hydrogeology`)**: Tested via two-sided **Mann–Whitney U test** (`p`-value and median fold ratio).
- **Multi-Group Strata (`Landscape Regime`, `Stream Order`, `River Basin`)**: Tested via **Kruskal–Wallis One-Way Analysis of Variance by Ranks** (`χ²` statistic, degrees of freedom `df`, and `p`-value).

---

### 3.4 Section 4: Within-Run Statistical Hypothesis Tests (`within_run_statistical_tests.csv`)

Section 4 conducts up to **10 formal statistical hypothesis tests** across four analytical categories to validate internal model performance and spatial behavior.

| Category | Test Name | Null Hypothesis (`H₀`) | Statistic Reported | How to Interpret the Result |
| :--- | :--- | :--- | :--- | :--- |
| **`Cross-Validation`** | **Paired Student's t-test on Squared CV Error (`S3N Combined` vs `OLS` / `Tail-Up` / `Tail-Down`)** | Both models have equal Mean Squared Error (`MSE`) across held-out 5-fold CV stations (`E[e₁² - e₂²] = 0`). Equivalent to a Diebold–Mariano test for nested/competing spatial predictors. | `t = X` (`% RMSE Reduction`) | A negative `t`-statistic with `p < 0.05` (`*`, `**`, `***`) proves that `S3N Combined` achieves a **statistically significant reduction in squared prediction error** at the exact same monitoring stations compared to `OLS` or single-component models. |
| **`Cross-Validation`** | **Paired Wilcoxon Signed-Rank Test on Absolute CV Error (`|e|`)** | Median absolute CV error (`MAE`) is identical between `S3N Combined` and the competing model. | `V = X` (`ΔMAE`) | Non-parametric counterpart to the paired `t`-test on `MSE`. Because squared errors can be dominated by 1–2 extreme outliers, a significant Wilcoxon `p < 0.05` confirms that `S3N Combined` improves accuracy across the **majority of typical stations**, not just extreme sites. |
| **`Spatial Components`** | **Tail-Up vs. Tail-Down Reach Correlation** | Flow-connected (`Tail-Up`) and hydrologic-distance (`Tail-Down`) reach predictions are uncorrelated across the network. | `r` (Pearson), `ρ` (Spearman) | <ul><li>`r > 0.90`: Very strong agreement—the fixed-effect covariates and dense station network constrain both spatial structures similarly.</li><li>`0.70 ≤ r ≤ 0.90`: Moderate divergence—flow-connected weighting at confluences (`Tail-Up`) produces noticeably different predictions in unmonitored tributaries compared to radial hydrologic distance (`Tail-Down`).</li></ul> |
| **`Kriging Uncertainty`** | **Mann–Whitney U Test (`Unmonitored` vs. `Monitored Reaches`)** | Kriging standard error (`se_model_sc`) has the same median on unmonitored reaches (`has_obs == 0`) and monitored reaches (`has_obs == 1`). | `W = X` (`Ratio = SE_unmon / SE_mon`) | Always expected to be highly significant (`p < 1e-16`, `Ratio > 1.0x`, e.g., `1.13x`). Confirms that universal kriging properly honors local station conditioning: uncertainty shrinks near observed stations and expands in unmonitored headwaters. |
| **`Environmental Gradient`** | **Lowland vs. Upland (`Mann–Whitney U`)**, **Aquifer vs. Non-Aquifer (`Mann–Whitney U`)**, **Strahler Stream Orders (`Kruskal–Wallis`)** | Reach predictions have identical medians across environmental or topological classes. | `W` or `ChiSq (df)` | Quantifies which physical driver exerts the strongest control on the determinand (compare the `Fold Ratio` and `ChiSq` magnitudes). |

---

### 3.5 Section 5: Cross-Validation Performance Summary (`cv_comparison_*.csv`)

Section 5 tabulates the 5-fold out-of-sample cross-validation metrics across each spatial domain (`All England`, `Humber`, `Thames`, `Severn`):
- **`CV R²`**: Out-of-sample coefficient of determination (`1 - SSE / SST`) on the model scale. Values `> 0.50` indicate strong predictive skill; comparing `S3N Combined` vs. `OLS` shows the exact fraction of variance explained by river-network spatial autocorrelation (`ΔR² = R²_S3N - R²_OLS`).
- **`CV RMSE` & `CV MAE`**: Root Mean Squared Error and Mean Absolute Error on the model scale. `RMSE` penalizes large errors heavily; `MAE` reflects typical station error.

---

## 4. Visual Diagnostic Guide: Interpreting the Figures

### 4.1 Figure 1: `fig_desc1_distributions_and_representativeness.png`

This 4-panel figure (`2×2` grid) diagnoses distributional shape, monitoring bias, component agreement, and uncertainty structure:

1. **Panel A — Reach Prediction vs. Observed Station Density**:
   - Compares kernel density curves for **All Reaches (`S3N`)** (solid blue with shaded fill), **Unmonitored Reaches** (dashed green), and **Observed Stations** (dotted red), with vertical dashed lines marking the **All-Reach Median** and **Observed Median**.
   - *What to look for*: If the red dotted curve (`Observed Stations`) is shifted to the right of the blue/green curves, monitoring stations are biased toward higher-concentration lowland/mainstem reaches.
2. **Panel B — Empirical Cumulative Distribution Functions (ECDF)**:
   - Plots `F(x)` from `0.0` to `1.0` for `S3N Combined`, `S3N Tail-Up`, `S3N Tail-Down`, and `Observed Stations`, annotated with the **Two-Sample KS `D`-statistic** and **`Pred / Obs Median Ratio`**.
   - *What to look for*: Horizontal separation between the blue (`S3N Combined`) and red (`Observed`) curves at `F(x) = 0.50` directly visualizes the median monitoring bias. Separation between `Tail-Up` (green) and `Tail-Down` (purple) shows whether flow-connected vs. tail-down covariance alters the overall network distribution.
3. **Panel C — Tail-Up vs. Tail-Down Component Agreement**:
   - Bivariate scatterplot of `pred_td_sc` vs. `pred_tu_sc` (or `pred_td` vs. `pred_tu`) across up to `8,500` randomly subsampled reaches, with the `1:1` identity line (dashed red) and fitted OLS regression line (solid black), annotated with Pearson `r` and Spearman `ρ`.
   - *What to look for*: Tight clustering along the `1:1` line indicates high robustness to spatial covariance specification. Scatter away from the `1:1` line occurs primarily in unmonitored tributaries joining large mainstems (where `Tail-Down` borrows strength from the mainstem, whereas `Tail-Up` isolates the unmonitored tributary).
4. **Panel D — Kriging Standard Error (`SE`) by Monitoring Status & Strahler Order**:
   - Horizontal boxplots comparing `se_model_sc` across `Monitored`, `Unmonitored`, and Strahler Orders `Ord 1` through `Ord 6+`.
   - *What to look for*: A systematic downward staircase in `Kriging SE` from `Ord 1` (highest uncertainty) down to `Ord 5–6+` (lowest uncertainty), demonstrating how higher monitoring density on mainstem rivers constrains kriging variance.

### 4.2 Figure 2: `fig_desc2_stratified_gradients.png`

This 4-panel figure (`2×2` grid) visualizes how physical geography and river topology shape the determinand across the network:

1. **Panel A — Predictions by Landscape & Hydrogeology Regime**:
   - Boxplots across `Upland (>=150m)`, `Lowland Non-Aq`, and `Lowland Aquifer`, with the global network median shown as a horizontal dashed red line and Kruskal–Wallis `χ²` in the header.
   - *What to look for*: Separation of interquartile boxes relative to the red global median line shows whether topography (`Upland` vs. `Lowland`) or bedrock permeability (`Non-Aq` vs. `Aquifer`) is the dominant control.
2. **Panel B — Longitudinal Profile Across Strahler Stream Orders**:
   - Plots the **Stratum Mean** (blue circles, solid line), **Stratum Median** (red squares, dashed line), and **`Q25–Q75` Interquartile Band** (shaded blue polygon) from `Order 1` headwaters to `Order ≥ 6` mainstems, with reach sample sizes (`n = ...`) along the bottom.
   - *What to look for*: An upward trajectory from `Order 1` to `Order 6+` indicates cumulative downstream enrichment (typical of agricultural and point-source nutrients). A flat or downward trajectory indicates headwater dominance or downstream dilution/degradation.
3. **Panel C — Distribution Across Major River Basins**:
   - Boxplots comparing `All England`, `Humber`, `Thames`, and `Severn`.
   - *What to look for*: Contrasts regional hydrochemical regimes (e.g., groundwater-dominated `Thames` vs. upland-influenced `Severn` and `Humber`).
4. **Panel D — Altitudinal Gradient Profile (10 Elevation Deciles)**:
   - Divides all network reaches into 10 equal-sized elevation deciles and plots the **Decile Mean** (green circles), **Decile Median** (dark green triangles), and **`Q25–Q75` IQR polygon** against mean decile elevation (`m AOD`).
   - *What to look for*: Identifies non-linear elevation thresholds (for instance, a sharp drop in Nitrate above `120–180 m AOD` marking the transition from lowland arable land to upland pastoral/moorland cover).

---

## 5. Worked Empirical Summary: Nitrate (`0117`) on ORN v2 (`output/results/det_0117_orn_rsn/descriptives/`)

Below is a concise reference summary illustrating how to synthesize the outputs of [`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R) for a real model run (**Nitrate as N, `0117`, across `37,396` reaches and `8,401` observed stations on ORN v2**):

1. **Distributional Structure (Section 1)**:
   - Across all `37,396` reaches in England, predicted Nitrate (`S3N Combined`) has a **Mean ± SD of `2.707 ± 1.981 mg/L`** (`CV = 73.2%`), a **Median of `2.344 mg/L` (`IQR = 2.671 mg/L`)**, and a **Geometric Mean of `1.881 mg/L`**.
   - The distribution is moderately right-skewed (`Skewness = 1.24`, `Kurtosis = 3.04`), spanning from **`0.068 mg/L` (`P05 = 0.369 mg/L`)** in pristine upland headwaters to **`25.43 mg/L` (`P95 = 6.235 mg/L`)** in intensive lowland catchments.
   - On the log-transformed model scale (`S3N Combined (Model Scale)`), skewness is negative (`-0.84`), and kriging uncertainty has a **Median `SE` of `0.506` (`IQR = 0.089`)**.
2. **Monitoring Network Representativeness (Section 2)**:
   - Observed EA RSN stations (`n = 8,401`, covering `22.5%` of reaches) have a **Median of `3.349 mg/L`**, whereas all network reaches have a **Median of `2.344 mg/L`** (and unmonitored reaches have a **Median of `2.078 mg/L`**).
   - This yields a **`Pred / Obs Median Ratio = 0.700x`** (`Mann–Whitney U p < 1e-16`, `Two-Sample KS D = 0.143`, `p < 1e-16`).
   - **Interpretation**: EA monitoring stations are strongly skewed toward accessible higher-order lowland streams (`58%–60%` monitored in Orders `3–5+` vs. `12.6%` in Orders `1–2`). Using raw station medians overstates whole-network median river nitrate by `42.9%`.
3. **Environmental & Topological Gradients (Section 3)**:
   - **Topography**: Lowland reaches (`< 150 m`, `80.6%` of network) have a median Nitrate of **`2.790 mg/L` (`Mean Z = +0.20`)**, which is **`3.31x` higher** (`p < 1e-16`) than Upland reaches (`≥ 150 m`, `19.4%` of network, median **`0.843 mg/L`**, `Mean Z = -0.82`).
   - **Hydrogeology**: Principal Aquifer reaches (`22.4%` of network) have a median of **`3.448 mg/L` (`Mean Z = +0.44`)**, **`1.66x` higher** (`p < 1e-16`) than Non-Aquifer reaches (**`2.075 mg/L`**), reflecting long-residence-time nitrate accumulation in Chalk, Limestone, and Sandstone groundwater.
   - **Combined Landscape Regime**: Moving from `Upland (>=150m)` (`0.843 mg/L`, `1.00x ref`) to `Lowland Non-Aquifer` (`2.353 mg/L`, **`2.79x`**) and `Lowland Aquifer` (`3.621 mg/L`, **`4.30x`**) explains massive hydrochemical variance (`Kruskal–Wallis χ² = 8,925.5`, `p < 1e-16`).
   - **Stream Order & Uncertainty**: Median Nitrate rises steadily from **`2.129 mg/L` in Headwaters (`Order 1–2`)** to **`3.112 mg/L` (`1.46x`) in Medium Tributaries (`Order 3–4`)** and **`3.591 mg/L` (`1.69x`) in Mainstem Rivers (`Order 5+`)**. Simultaneously, median `Kriging SE` drops from **`0.517` in Headwaters** to **`0.464` in Mainstems** because station density increases nearly five-fold downstream (`12.6%` to `60.2%`).
4. **Within-Run Hypothesis Tests & Spatial Components (Sections 4 & 5)**:
   - In 5-fold cross-validation (`All England`), `S3N Combined` achieves **`CV R² = 0.5658`** (`RMSE = 0.5067`), outperforming non-spatial `OLS` (`CV R² = 0.4252`, `RMSE = 0.5830`) by **`ΔR² = +0.1406`**.
   - Paired hypothesis tests confirm that `S3N Combined` significantly reduces squared CV error compared to `OLS` (`t = -15.32`, **`-13.1%` RMSE reduction**, `p < 1e-16`), `Tail-Up` alone (`t = -12.60`, **`-3.4%` RMSE reduction**, `p < 1e-16`), and `Tail-Down` alone (`t = -3.49`, **`-0.5%` RMSE reduction**, `p = 4.91e-04`).
   - Across all `37,396` reaches, `Tail-Up` and `Tail-Down` predictions are highly correlated (`Pearson r = 0.9462`, `Spearman ρ = 0.9628`), confirming stable spatial interpolation across the network.

---

## 6. Related Documentation

- **S3N Modelling & Universal Stream Kriging Theoretical Guide**: [`./s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md)
- **Multi-Run & Cross-Determinand Comparison Interpretation Guide**: [`./model_comparison_interpretation_guide.md`](./model_comparison_interpretation_guide.md)
- **End-to-End Pipeline Usage Guide**: [`./usage_guide.md`](./usage_guide.md)
- **Custom Water Network Guide (`BYOWN`)**: [`./custom_water_network_guide.md`](./custom_water_network_guide.md)
- **Preprocessing & Graph Topology Guide**: [`./preprocessing_pipeline.md`](./preprocessing_pipeline.md)
- **ORN Variables & Covariates Reference**: [`./orn_variables_reference.md`](./orn_variables_reference.md)
