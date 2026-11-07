# S3N Multi-Run & Cross-Determinand Model Comparison — Interpretation Guide

## 1. Overview & Supported Comparison Modes

The **S3N Model Run Comparison & Statistical Evaluation** script ([`../scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R)) performs automated statistical hypothesis testing, reach-by-reach spatial co-location analysis, environmental gradient comparison, parameter stability testing, and base-R visualization across **two or more completed `S3N` model output directories**. Internally, it also calls [`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R) to compute unified descriptive and stratified statistics for every model run being compared.

[`../scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R) is designed to support **four primary scientific and operational comparison workflows**:

| Comparison Mode | Typical Example | Primary Analytical Goal |
| :--- | :--- | :--- |
| **1. Cross-Network Comparison** | `ORN v2` (`output/det_0117`) vs. `EU-Hydro Thames` (`output/det_0117_euhydro_thames`) | Determine which digital river network topology produces higher out-of-sample cross-validation accuracy (`CV R²`, `CV RMSE`), lower kriging uncertainty (`Kriging SE`), and whether reach predictions are consistent across co-located segments (`≤ 250 m`). |
| **2. Covariate Specification Comparison** | Default 8 covariates vs. reduced model (`--vars=elevation,slope,aquifer`) | Test whether adding or removing catchment land-cover/geology predictors significantly improves station cross-validation accuracy (via paired Diebold–Mariano / `t`-tests and Wilcoxon tests) or alters spatial variance components (`σ²`, `τ²`, `φ`). |
| **3. Temporal Slice Comparison** | Seasonal (`--by=season`: Winter, Spring, Summer, Autumn), Annual (`--by=year`), or Monthly (`--by=month`) | Quantify seasonal flushing vs. summer point-source concentration regimes, multi-year water quality trends, and how spatial autocorrelation range (`φ`) changes across flow regimes. |
| **4. Cross-Determinand Comparison** | `Nitrate (0117)` vs. `Orthophosphate (0180)`, or `Nitrate (0117)` vs. `RHS Habitat Quality (hqa_score)` | Analyze reach-by-reach co-occurrence (`4-Quadrant Concordance`), multi-pollutant hotspot overlap (`Jaccard > P75 / > P90`), power-law scaling (`β`), **freshwater nutrient stoichiometry (`N:P` mass & molar ratios and Redfield nutrient limitation regimes)**, or multi-stressor `ΔZ` divergence. |

---

## 2. Quick-Start CLI Reference & Generated Outputs

### 2.1 Running the Script

```bash
# 1. Compare two or more specific output folders directly
Rscript scripts/04_compare_model_runs.R output/det_0117 output/det_0117_euhydro_thames

# 2. Cross-determinand comparison (e.g., Nitrate vs. Orthophosphate on ORN v2)
Rscript scripts/04_compare_model_runs.R \
  output/det_0117 output/det_0180 \
  --labels="Nitrate (0117),Orthophosphate (0180)" \
  --out=output/comparisons/det_0117_vs_det_0180_orn

# 3. Auto-fit and compare temporal slices (season, year, or month)
Rscript scripts/04_compare_model_runs.R --det=0117 --by=season
Rscript scripts/04_compare_model_runs.R --det=0117 --by=year --slices=2021,2022,2023,2024

# 4. List all available completed model runs in output/
Rscript scripts/04_compare_model_runs.R --list

# 5. Specify custom workspace directories
Rscript scripts/04_compare_model_runs.R output/det_0117 output/det_0180 \
  --gis-dir=gis --data-dir=data --output-dir=output
```

### 2.2 Output Files Written to `output/comparisons/<comparison_slug>/`

| Output File | Format | Description |
| :--- | :--- | :--- |
| `comparison_report.md` | Markdown | Standalone narrative report containing executive winner determination, all comparison tables, statistical test interpretations, reach relationship synthesis, and embedded figure links. |
| `cv_metrics_comparison.csv` | CSV | 5-fold cross-validation metrics (`CV R²`, `CV RMSE`, `CV MAE`, `OLS` baselines, `ΔR² vs OLS`, `% RMSE Reduction`) across all runs and spatial domains. |
| `statistical_tests_summary.csv` | CSV | Up to 12 formal parametric and non-parametric hypothesis tests comparing CV errors, reach distributions, kriging uncertainty, co-located bias, and 4-quadrant co-occurrence. |
| `descriptive_stats_summary.csv` | CSV | Full 18-column distributional summary across all 7 variable groups for every compared run (generated via [`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R)). |
| `stratified_results_summary.csv` | CSV | Multi-run stratified environmental, topological, and basin statistics across all 5 stratification dimensions. |
| `reach_relationship_summary.csv` | CSV | Pairwise reach-level relationship metrics: correlations (`r`, `ρ`), OLS scaling slope, 4-quadrant co-occurrence, Odds Ratio, Jaccard hotspot overlap (`> P75`, `> P90`), and `N:P` stoichiometry or `ΔZ` divergence. |
| `parameter_comparison.csv` | CSV | Side-by-side fixed-effect (`β`) and spatial covariance (`sigma.sq`, `tau.sq`, `phi`) bootstrap estimates with Two-Sample Wald `Z`-test `p`-values. |
| `fig_comp1_cv_performance.png` | PNG (2400×1700) | 4-panel cross-validation comparison (`CV R²`, `CV RMSE`, `CV MAE` vs. `OLS`, and `% RMSE Reduction`). |
| `fig_comp2_reach_distributions.png` | PNG (2400×1700) | 4-panel distributional comparison (kernel densities, ECDFs, Kriging `SE` boxplots, and co-located reach scatter or `Q25/Med/Q75` quantile profiles). |
| `fig_comp3_parameter_comparison.png` | PNG (2400×1600) | 2-panel bootstrap parameter comparison (fixed-effect `β` forest plot with 95% CIs and spatial covariance parameters `σ²`, `τ²`, `φ`). |
| `fig_comp4_spatial_difference_map.png` | PNG (2400×1400) | 3-panel BNG spatial map (Run 1 predictions, Run 2 predictions, and co-located reach difference `Run 2 - Run 1` or standardized `ΔZ = Z₂ - Z₁`). |
| `fig_comp5_determinand_relationships.png` | PNG (2400×1800) | 4-panel reach relationship figure (4-quadrant co-occurrence scatter, stratified median comparison, Strahler longitudinal profiles, and `N:P` molar stoichiometry or `ΔZ` histogram). |

---

## 3. Section-by-Section Guide to the Comparison Report

### 3.1 Sections 1–3: Executive Summary, Winner Determination & Cross-Validation Metrics (`cv_metrics_comparison.csv`)

Sections 1 through 3 evaluate out-of-sample predictive skill across held-out 5-fold cross-validation (`CV`) stations.

#### A. Understanding Cross-Validation Metrics

| Metric | Definition | Direction of Better Performance | What It Measures |
| :--- | :--- | :--- | :--- |
| **`CV R² (S3N)`** | `1 - Σ(y_i - ŷ_i)² / Σ(y_i - ȳ)²` on the model scale | **Higher is better** (`0` to `1`) | Fraction of total station variance explained out-of-sample. Sensitive to the total variance (`SST`) of the sampled domain. |
| **`CV RMSE (S3N)`** | `√(mean((y_i - ŷ_i)²))` on the model scale | **Lower is better** | Root Mean Squared Error in model-scale units (`log(Y)` when both runs are log-transformed). Penalizes large outlier errors heavily. |
| **`CV MAE (S3N)`** | `mean(|y_i - ŷ_i|)` on the model scale | **Lower is better** | Mean Absolute Error. Robust measure of typical station prediction error. |
| **`ΔR² (S3N - OLS)`** | `CV R² (S3N) - CV R² (OLS)` | **Higher positive value** | Incremental predictive power gained specifically from the branching river-network spatial covariance (`Tail-Up + Tail-Down`) beyond catchment land-cover/geology fixed effects (`OLS`). |
| **`RMSE Red. vs OLS (%)`** | `100 × (RMSE_OLS - RMSE_S3N) / RMSE_OLS` | **Higher positive `%`** | Percentage reduction in prediction error achieved by spatial kriging relative to non-spatial regression. |

#### B. Why Global `CV R²` and Global `CV RMSE` Can Rank Runs Differently

When comparing two runs that cover **different geographic extents or station sample sizes** (for example, `ORN v2` across all of England with `8,401` stations vs. `EU-Hydro` restricted to the `Thames` basin with `996` stations), you will often observe that:
- The **national run (`ORN v2`)** has a **higher Global `CV R²`** (`0.5658` vs. `0.4533`), whereas
- The **regional run (`EU-Hydro Thames`)** has a **lower Global `CV RMSE`** (`0.4730` vs. `0.5067`).

**Why does this happen mathematically?**
Because `CV R² = 1 - (RMSE² / Var(y))`. A national dataset spanning pristine upland Cumbria (`< 0.2 mg/L` Nitrate) to intensive arable East Anglia (`> 15 mg/L`) has a much larger total variance `Var(y)` than a single lowland basin (`Thames`). Even if prediction error (`RMSE`) is slightly smaller inside the Thames, dividing by a smaller regional `Var(y)` produces a lower `R²`.

#### C. The Gold Standard Tie-Breaker: Head-to-Head Matched Co-Located Station CV

To eliminate domain-size and variance artifacts when comparing two networks (`ORN` vs. `EU-Hydro`) or two specifications on the same determinand, Section 2 automatically performs a **Head-to-Head Matched Station Evaluation**:
1. It spatially pairs held-out CV stations between Run 1 and Run 2 that lie within **`≤ 250 m` Euclidean BNG distance** (or share the exact `RivID` on the same network) and have identical observed values (`|obs₁ - obs₂| < 1e-6`).
2. On this **exact same set of `M` physical monitoring stations** (e.g., `M = 847` co-located Thames stations), total variance `Var(y)` is identical for both runs, so **`Matched R²`, `Matched RMSE`, and `Matched MAE` are 100% directly comparable**.
3. The script then runs a **Paired Student's t-test (Diebold–Mariano test)** on squared errors (`e₁² - e₂²`) and a **Paired Wilcoxon Signed-Rank test** on absolute errors (`|e₁| - |e₂|`) across those `M` stations to declare whether one network/model is **statistically superior** (`p < 0.05`) or **statistically indistinguishable (tied)** (`p ≥ 0.05`).

> **Note on Cross-Determinand Comparisons**: When comparing two *different* water quality variables (e.g., `Nitrate 0117` vs. `Orthophosphate 0180`), `CV RMSE` and `CV MAE` reflect the intrinsic natural variability and measurement noise of each chemical/ecological variable rather than a "winner/loser" competition. For cross-determinand runs, focus on **`CV R²`** (which is dimensionless) and **`ΔR² vs OLS`** to see which determinand is more strongly governed by catchment covariates vs. river-network spatial autocorrelation.

---

### 3.2 Section 4: Formal Statistical Hypothesis Tests (`statistical_tests_summary.csv`)

Section 4 executes up to **12 formal statistical tests** grouped into five analytical categories. Every test reports the test statistic, degrees of freedom (`df` or sample size `n`), `p`-value, significance stars (`***` `< 0.001`, `**` `< 0.01`, `*` `< 0.05`, `.` `< 0.10`, `ns` `≥ 0.10`), and a plain-English interpretation.

| Category | Test Name | When It Runs | Null Hypothesis (`H₀`) | How to Interpret the Statistic & `p`-Value |
| :--- | :--- | :--- | :--- | :--- |
| **`Cross-Validation`** | **1. Fisher's r-to-z Test on `CV R²`** | Always (Pairwise) | Out-of-sample correlation `r = √R²` is equal between Run 1 and Run 2 (`z = (atanh(r₁) - atanh(r₂)) / √(1/(n₁-3) + 1/(n₂-3))`). | `p < 0.05` indicates one model explains a statistically higher proportion of variance across its station set. |
| **`Cross-Validation`** | **2. Unpaired Mann–Whitney U Test on `|CV Error|`** | Always (Pairwise) | Median absolute out-of-sample station error (`|y - ŷ|`) is equal across all stations in Run 1 (`n₁`) and Run 2 (`n₂`). | Non-parametric test across all available stations (even when station counts differ). Negative `ΔMed |e|` means Run 1 has smaller median error. |
| **`Cross-Validation`** | **3. Unpaired Welch's t-Test on Squared CV Error** | Always (Pairwise) | Mean Squared Error (`MSE`) is equal across all stations in Run 1 and Run 2, allowing unequal variances. | Parametric counterpart to Test 2 across all stations. |
| **`Cross-Validation`** | **4. Paired Student's t-Test on Squared CV Error (Co-located Stations)** | When `M ≥ 10` co-located stations share the same target variable | At the exact same `M` physical stations, `E[e₁² - e₂²] = 0` (Diebold–Mariano equal predictive accuracy test). | **The most rigorous test for comparing two networks or covariate sets.** If `p ≥ 0.05` (`ns`), the two networks/models have statistically indistinguishable MSE at matched stations. |
| **`Cross-Validation`** | **5. Paired Wilcoxon Signed-Rank Test on `|CV Error|` (Co-located Stations)** | When `M ≥ 10` co-located stations share the same target variable | Median difference in absolute error `median(|e₁| - |e₂|) = 0` across matched stations. | Robust to 1–2 extreme outlier stations. If Test 4 is `ns` but Test 5 is `*`, one model performs slightly better at the median station while having similar outlier errors. |
| **`Reach Predictions`** | **6. Welch's Two-Sample t-Test on Reach Predictions** | Always (Pairwise) | Mean reach prediction (`pred_combined`) is equal between Run 1 (`N₁`) and Run 2 (`N₂`). | Tests whether the overall network mean differs (due to regional extent, seasonal shift, or determinand scale). |
| **`Reach Predictions`** | **7. Two-Sample Kolmogorov–Smirnov (KS) Test** | Always (Pairwise) | The cumulative distribution functions `F₁(x)` and `F₂(x)` of reach predictions are identical (`D = max |F₁(x) - F₂(x)|`). | Because `N` is large (`19,000–37,000` reaches), focus on the effect size **`D`**: `D < 0.05` is trivial, `0.05–0.15` is moderate, and `D > 0.20` indicates a major distributional shift. |
| **`Reach Predictions`** | **8. One-Way ANOVA (`F`-test) & 9. Kruskal–Wallis (`χ²`)** | When comparing `K ≥ 3` runs (e.g., 4 seasons or multiple years) | All `K` model runs have identical mean (ANOVA) or median (Kruskal–Wallis) reach predictions. | Omnibus multi-group test confirming significant seasonal or inter-annual variation across the river network. |
| **`Kriging Uncertainty`** | **10. Mann–Whitney U Test on `Kriging SE`** | Always (Pairwise) | Median universal kriging standard error (`se_model_sc`) is equal between Run 1 and Run 2. | Identifies which network or model yields tighter local prediction intervals (`Lower Med SE = Higher Confidence`). |
| **`Kriging Uncertainty`** | **11. Fligner–Killeen Test on `Kriging SE` Variance** | Always (Pairwise) | The dispersion/spread of kriging uncertainty across reaches is homogeneous between Run 1 and Run 2 (`med χ²`). | A significant result indicates one network has more uniform uncertainty across reaches, whereas the other has wider variation between well-monitored and unmonitored sub-catchments. |
| **`Spatial Co-location` & `Co-Occurrence`** | **12. Co-located Reach Correlation (`r`, `ρ`), Paired Bias `t`-Test, & 4-Quadrant Chi-Squared (`χ²`) Test** | When `N_matched ≥ 20` co-located reaches exist | Co-located reaches are uncorrelated (`r = 0`), have zero mean difference (`Δ = 0`), and sort independently into High/Low median quadrants (`Odds Ratio = 1.0`). | Summarizes reach-by-reach agreement (expanded in full detail in Section 6). |

---

### 3.3 Section 5: Descriptive Statistics & Stratified Gradient Comparison (`descriptive_stats_summary.csv` & `stratified_results_summary.csv`)

Section 5 embeds the full multi-run outputs from [`../scripts/05_descriptive_stats.R`](../scripts/05_descriptive_stats.R):
- **Section 5.1 (`Full Distributional Summary`)**: Places all 7 variable groups (`Observed Stations`, `S3N Combined`, `Unmonitored Reaches`, `Tail-Up`, `Tail-Down`, `Model Scale`, and `Kriging SE`) side by side across every run.
- **Section 5.2 (`Monitoring Network Representativeness`)**: Compares `Pred / Obs Median Ratio` and `KS D-Stat` across runs.
  - *Cross-Network Insight*: In the `ORN v2` national network, `Nitrate (0117)` has `Pred / Obs Median Ratio = 0.700x` because England contains large unmonitored upland regions (Pennines, Lake District, Dartmoor). In contrast, in `EU-Hydro (Thames)`, `Pred / Obs Median Ratio = 0.996x` (`p = 0.097 ns`) because the Thames basin is almost entirely lowland (`94.7% < 150 m`), so its monitoring stations match the basin-wide median almost perfectly!
- **Section 5.3 (`Stratified Environmental & Topological Gradients`)**: Compares how each model run responds to `Topography`, `Hydrogeology`, `Landscape Regime`, `Stream Order`, and `River Basin`.
  - *Cross-Determinand Insight*: Comparing `Fold Ratio vs Ref` across determinands immediately reveals contrasting environmental controls. For example, in `Hydrogeology`, moving from `Non-Aquifer` (`1.00x ref`) to `Principal Aquifer` increases median **Nitrate (`0117`) by `1.66x`** (`2.075` to `3.448 mg/L`, due to groundwater nitrate storage), whereas median **Orthophosphate (`0180`) increases by only `1.13x`** (`0.056` to `0.064 mg/L`) because phosphorus binds to soil/aquifer matrices and is primarily controlled by surface point-source `STW` discharges and lowland urban/agricultural runoff.

> For a complete definition of every column in Section 5 (`CV (%)`, `Skewness`, `Kurtosis`, `MAD`, `Mean Z-Score`, `Fold Ratio vs Ref`), see the companion guide: [`./descriptive_stats_interpretation_guide.md`](./descriptive_stats_interpretation_guide.md).

---

### 3.4 Section 6: Reach-Level Relationship, Co-Occurrence & Stoichiometric Analysis (`reach_relationship_summary.csv`)

Section 6 is one of the most powerful analytical components of [`../scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R). It pairs individual river reaches between Run 1 and Run 2 and evaluates how their predictions co-vary across space.

#### A. How Reaches Are Paired Between Model Runs (`Match_Method`)

1. **Exact Topological `RivID` Match (`Same Network`)**:
   - When both runs use the exact same digital river network (e.g., comparing `Nitrate 0117` vs. `Orthophosphate 0180` on `ORN v2`, or comparing `Winter` vs. `Summer` on `ORN v2`), every reach is matched 1-to-1 by its exact `RivID` (`N = 37,396` matched pairs, `100%` coverage).
2. **Spatial Nearest-Neighbor Match (`Different Networks, ≤ 250 m BNG`)**:
   - When comparing two different digital river networks (e.g., `ORN v2` with `37,396` reaches vs. `EU-Hydro Thames` with `19,642` reaches), `RivID` numbering differs and spatial discretization (segment splitting at canals/ditches) varies.
   - The script uses a fast **2D British National Grid (`EPSG:27700`) spatial hash grid** to pair each reach centroid with its nearest counterpart within a strict **`≤ 250 m` tolerance** (e.g., matching `N = 5,321` co-located Thames reaches with a median centroid separation of just `41.5 m`).

#### B. Bivariate Correlation & Scaling Equations

| Metric | Definition | How to Interpret It |
| :--- | :--- | :--- |
| **`Pearson r` (Raw & Log)** | Linear correlation of raw predictions (`pred_combined`) and log-transformed predictions (`log10`). | <ul><li>**Same Determinand across Networks**: Expected to be very high (`r ≈ 0.85–0.98`), confirming both networks capture the same spatial water quality pattern.</li><li>**Cross-Determinand (`Nitrate` vs. `Orthophosphate`)**: Moderate positive correlation (`r_raw = 0.453`, `r_log = 0.678`, `p < 1e-16`), showing shared lowland agricultural/urban enrichment alongside distinct local sources.</li></ul> |
| **`Spearman ρ`** | Rank correlation across all paired reaches. | Insensitive to extreme point-source spikes. If `ρ >> r_raw`, extreme outliers in one determinand are dampening the linear correlation while rank order remains strongly coupled. |
| **`Tail-Up r` & `Tail-Down r`** | Component-specific Pearson correlations (`pred_tu` and `pred_td`). | Shows whether flow-connected (`Tail-Up`) or radial hydrologic-distance (`Tail-Down`) components couple more strongly between the two runs. |
| **`Observed Station r`** | Empirical correlation at co-located monitoring stations (`has_obs == 1` in both runs). | Compares the **raw empirical correlation** against the **kriged network correlation**. Smooth universal kriging surfaces typically exhibit stronger correlation (`r_pred > r_obs`) because sampling-day white noise (`τ²`) is filtered out. |
| **`OLS Linear Fit` & `Log-Log Power-Law Scaling`** | Linear: `Pred_2 = α + β × Pred_1`<br>Power-Law: `log10(Pred_2) = α_log + β_log × log10(Pred_1)` (i.e., `Det_2 ∝ Det_1^β_log`) | In cross-determinand analysis, the **power-law exponent `β_log`** reveals relative scaling:<ul><li>`β_log ≈ 1.0`: Isometric scaling (both determinands increase proportionally).</li><li>`β_log < 1.0`: Sub-linear scaling (e.g., for `Orthophosphate` vs. `Nitrate`, `β_log = 0.695`: a `10-fold` increase in Nitrate corresponds to a `10^0.695 ≈ 4.95-fold` increase in Orthophosphate).</li><li>`β_log > 1.0`: Super-linear scaling.</li></ul> |

#### C. 4-Quadrant Co-Occurrence, Concordance & Hotspot Overlap (`Jaccard Index`)

To help catchment managers identify **multi-pollutant hotspots** vs. **single-stressor reaches**, every paired reach is classified relative to each run's global median (`Med₁` and `Med₂`):

| Quadrant | Classification Rule | Environmental Meaning (e.g., Nitrate `Run 1` vs. Phosphate `Run 2`) |
| :--- | :--- | :--- |
| **`Q1: High / High`** | `Run 1 > Med₁` **AND** `Run 2 > Med₂` | **Co-Enriched / Multi-Stressor Reaches**: Elevated in both pollutants (e.g., intensive lowland mixed agriculture + high urban `STW` effluent loading). |
| **`Q2: Low / Low`** | `Run 1 ≤ Med₁` **AND** `Run 2 ≤ Med₂` | **Co-Clean / Reference Reaches**: Below-median in both pollutants (e.g., upland headwaters, forested/moorland catchments). |
| **`Q3: High Run 1 / Low Run 2`** | `Run 1 > Med₁` **AND** `Run 2 ≤ Med₂` | **Run 1-Dominant Reaches**: E.g., rural Chalk/Limestone groundwater catchments with high agricultural Nitrate leaching (`High N`) but low rural population/sewage effluent (`Low P`). |
| **`Q4: Low Run 1 / High Run 2`** | `Run 1 ≤ Med₁` **AND** `Run 2 > Med₂` | **Run 2-Dominant Reaches**: E.g., impermeable clay/upland-margin catchments with low groundwater Nitrate (`Low N`) but local wastewater discharge or livestock soil erosion (`High P`). |

From this `2×2` contingency table, the script computes five key spatial co-occurrence indices:

1. **`Overall Concordance (%) = Q1 (%) + Q2 (%)`**:
   - Under complete spatial independence, random sorting produces **`50.0%` concordance** (`25%` in each quadrant).
   - Values **`> 65%`** indicate strong spatial coupling; values **`> 85%`** indicate near-identical spatial structure (typical when comparing two networks on the same determinand, e.g., `88.7%` for `ORN` vs. `EU-Hydro`).
2. **`Odds Ratio (OR) = (Q1 × Q2) / (Q3 × Q4)` & `Chi-Squared (χ²)`**:
   - Quantifies how many times more likely a reach is to be above-median in Run 2 given that it is above-median in Run 1:
     - `OR = 1.0`: Independent (`p ≥ 0.05`).
     - `OR = 3.0 – 10.0`: Strong positive co-occurrence (e.g., `OR = 7.42`, `χ² = 7,968`, `p < 1e-16` between Nitrate and Orthophosphate).
     - `OR > 25.0`: Extremely tight spatial concordance (e.g., `OR = 62.6`, `χ² = 3,189` between `ORN` and `EU-Hydro` Nitrate).
     - `OR < 1.0`: Inverse relationship (e.g., comparing a pollutant like `Nitrate` against a positive ecological indicator like `WHPT ASPT`).
3. **`Upper-Quartile Hotspot Overlap (Jaccard > P75)`**:
   - Let `A = {reaches > P75 in Run 1}` (top `25%` worst reaches) and `B = {reaches > P75 in Run 2}`.
   - The **Jaccard Similarity Index** is `J_75 = |A ∩ B| / |A ∪ B|`:
     - Under random independence, `J_75 = 0.25² / (0.25 + 0.25 - 0.25²) = 14.3%`.
     - `J_75 = 35% – 55%`: Substantial multi-pollutant hotspot overlap (e.g., `37.6%` between Nitrate and Orthophosphate).
     - `J_75 > 70%`: Near-identical upper-quartile hotspot map (e.g., `77.1%` between `ORN` and `EU-Hydro`).
4. **`Top-Decile Extreme Hotspot Overlap (Jaccard > P90)`**:
   - Same calculation restricted to the top `10%` extreme hotspots (`> P90`). Under random independence, `J_90 = 5.3%`.
5. **`Conditional Hotspot Probability P(R2 > P75 | R1 > P75)`**:
   - Directly answers the management question: *"If a reach is in the top 25% worst hotspots for Run 1, what is the probability (`%`) that it is also in the top 25% worst hotspots for Run 2?"* (Random baseline = `25.0%`; for Nitrate vs. Orthophosphate it is **`54.6%`**—more than double random chance).

#### D. Freshwater Nutrient Stoichiometry (`N:P` Ratios) vs. Standardized `ΔZ` Divergence

Section 6 automatically adapts its final subsection depending on what is being compared:

##### Case 1: Nitrogen vs. Phosphorus Comparison (`N:P Stoichiometric & Nutrient Limitation Analysis`)
When one run is a **Nitrogen** determinand (`0117` Nitrate as N, `0116` TON, `0111` Ammonia) and the other is a **Phosphorus** determinand (`0180` Orthophosphate as P, `0192` Total P), both are reported in `mg/L` of elemental `N` and `P`. The script automatically computes reach-by-reach **Mass `N:P` Ratios** and **Atomic (Molar) `N:P` Ratios**:

- **Conversion Formula**: Because the atomic weight of Nitrogen is `14.0067 g/mol` and Phosphorus is `30.97376 g/mol`, the molar conversion factor is:
  `Molar N:P = (N_mg_L / 14.0067) / (P_mg_L / 30.97376) = 2.21135 × Mass N:P`
- **Ecological Nutrient Limitation Regimes (Redfield `16:1` / Guildford & Hecky 2000 Freshwater Thresholds)**:
  Every reach in the river network is classified into one of three stoichiometric limitation classes based on its predicted **Molar `N:P` ratio**:
  1. **`P-Limited (Molar N:P > 20:1)`**: Dissolved inorganic nitrogen is in abundant excess relative to phosphorus (typical of groundwater-fed Chalk/Limestone rivers and agricultural nitrate catchments). Benthic algae and macrophyte growth are stoichiometrically constrained by **Phosphorus availability**—meaning **`P` stripping at wastewater treatment works (`STWs`) or agricultural `P` control** is required to curb eutrophication.
  2. **`Co-Limited / Balanced (Molar N:P 10:1 – 20:1)`**: Straddles the classic marine/freshwater **Redfield Ratio (`16:1` molar / `7.2:1` mass)**. Algal communities may switch between N and P limitation seasonally; **dual nutrient management (`N + P`)** is required.
  3. **`N-Limited (Molar N:P < 10:1)`**: Phosphorus is in stoichiometric excess relative to nitrogen (typical of low-nitrate upland/clay streams receiving local sewage effluent or livestock runoff). Primary production is constrained by **Nitrogen availability**.

##### Case 2: General Comparisons (`Standardized Z-Score Divergence ΔZ = Z₂ - Z₁`)
For all non-`N:P` pairs (such as cross-network comparisons, seasonal comparisons, or `Nitrate` vs. `WHPT`/`HQA`), the script standardizes each run to zero mean and unit variance (`Z_i = (Pred_i - μ_i) / σ_i`) and computes reach-level divergence **`ΔZ = Z₂ - Z₁`**:
- **`|ΔZ| ≤ 1.0 SD`**: Concordant relative standing between the two runs.
- **`ΔZ > +1.0 SD`**: Disproportionately higher in **Run 2** relative to Run 1.
- **`ΔZ < -1.0 SD`**: Disproportionately higher in **Run 1** relative to Run 2.

---

### 3.5 Section 7: Fixed-Effect (`β`) & Spatial Covariance Parameter Comparison (`parameter_comparison.csv`)

Section 7 compares the fitted regression coefficients and spatial autocovariance parameters from `parameter_estimates_*.csv` across runs using **Two-Sample Wald `Z`-Tests** derived from the 95% parametric bootstrap confidence intervals:

- **Bootstrap Standard Error Recovery**: For each parameter `θ` with 95% bootstrap interval `[CI_Lower, CI_Upper]`, the asymptotic standard error is `SE(θ) = (CI_Upper - CI_Lower) / (2 × 1.96)`.
- **Two-Sample Wald `Z`-Statistic**:
  `Z = (θ₁ - θ₂) / √(SE(θ₁)² + SE(θ₂)²)`, with two-sided `p = 2 × (1 - Φ(|Z|))`.

#### How to Interpret Parameter Differences

| Parameter Group | Parameter Names | What a Significant Wald `Z`-Test (`p < 0.05`) Means |
| :--- | :--- | :--- |
| **Fixed Effects (`β`)** | `Intercept`, `elevation`, `slope`, `arable`, `improved_grass`, `urban`, `woodland`, `aquifer`, `stw_log` | <ul><li>**Across Networks (`ORN` vs. `EU-Hydro`)**: If `p ≥ 0.05` (`ns`) for most slopes (e.g., `elevation`, `slope`, `arable`, `woodland`, `stw_log`), the physical environmental relationships are **structurally invariant** to network choice. A significant shift in `aquifer` or `Intercept` usually reflects regional subsetting (e.g., Thames Chalk vs. national aquifer mix).</li><li>**Across Determinands (`Nitrate` vs. `Orthophosphate`)**: Highlights divergent physical drivers—for instance, `aquifer` has a much stronger positive slope for `Nitrate` (`+0.391` vs. `+0.172`, `p < 1e-16`), whereas `stw_log` (`STW` sewage load) has a stronger slope for `Orthophosphate` (`+0.072` vs. `+0.044`, `p < 1e-16`).</li></ul> |
| **Partial Sill (`σ²`)** | `Tail-Up: sigma.sq`, `Tail-Down: sigma.sq` | Magnitude of spatially structured variance along the river network. Higher `sigma.sq` indicates stronger spatial autocorrelation among neighboring reaches. |
| **Nugget Effect (`τ²`)** | `Tail-Up: tau.sq`, `Tail-Down: tau.sq` | Micro-scale variability + measurement noise at zero distance (`h = 0`). If `tau.sq << sigma.sq` (high **Structured Variance Ratio** `σ² / (σ² + τ²)`), river-network kriging adds substantial predictive precision over `OLS`. |
| **Spatial Range (`φ`)** | `Tail-Up: phi`, `Tail-Down: phi` | Exponential decay rate (`exp(-φ × d)` in `km⁻¹`, where effective correlation range is `3 / φ km`) or spherical range parameter. Smaller `φ` in `exp(-φ × d)` means spatial correlation persists over longer river distances. |

---

## 4. Visual Diagnostic Guide: Interpreting the 5 Comparison Figures

### 4.1 Figure 1: `fig_comp1_cv_performance.png`
- **Panel A (`CV R²`)**: Grouped bar chart comparing `S3N Combined` vs. `OLS`, `Tail-Up`, and `Tail-Down` across runs. Look at the height of the primary bar (`S3N Combined`) and the gap above the `OLS` bar (`ΔR²`).
- **Panel B (`CV RMSE`)**: Grouped bar chart of Root Mean Squared Error (annotated with `LOWER IS BETTER` when runs share the same determinand).
- **Panel C (`CV MAE`)**: Grouped bar chart of Mean Absolute Error (`S3N` vs. `OLS`).
- **Panel D (`% RMSE Reduction vs. OLS`)**: Horizontal bar chart showing the exact percentage of prediction error eliminated by adding spatial stream network covariance (`100 × (RMSE_OLS - RMSE_S3N) / RMSE_OLS`).

### 4.2 Figure 2: `fig_comp2_reach_distributions.png`
- **Panel A (`Kernel Density Comparison`)**: Overlaid density curves of reach predictions (`log10` scale for skewed nutrients, linear scale for ecological scores), annotated with the Two-Sample KS `D`-statistic.
- **Panel B (`Empirical CDFs`)**: Overlaid cumulative distribution curves `F(x)` from `0` to `1`. Horizontal distance between solid lines (`Predictions`) and dashed lines (`Observed Stations`) shows monitoring bias within each run.
- **Panel C (`Kriging SE Distribution`)**: Boxplots of `se_model_sc` across runs, annotated with Mann–Whitney `p`-values. Lower boxes indicate higher spatial prediction certainty.
- **Panel D (`Co-located Reach Scatter` or `Multi-Run Quantile Profile`)**: For 2 runs, plots matched reach predictions (`Run 2` vs. `Run 1`) with the `1:1` identity line and Pearson `r` / Spearman `ρ`. For `K ≥ 3` runs (e.g., 4 seasons), plots the trajectory of `Median`, `Q25`, and `Q75` across runs.

### 4.3 Figure 3: `fig_comp3_parameter_comparison.png`
- **Panel A (`Fixed-Effect Regression Coefficients β`)**: Forest plot of all covariate slopes (`β`) with horizontal `95%` bootstrap confidence interval bars and a vertical dashed zero line (`β = 0`). Non-overlapping error bars between runs correspond to statistically significant Wald `Z`-test differences (`p < 0.05`).
- **Panel B (`Spatial Covariance Parameters`)**: Grouped bar chart with `95%` bootstrap error bars for `Tail-Up` and `Tail-Down` variance components (`sigma.sq`, `tau.sq`, `phi`).

### 4.4 Figure 4: `fig_comp4_spatial_difference_map.png`
- **Panels 1 & 2 (`Run 1` and `Run 2` Spatial Maps)**: BNG (`EPSG:27700`) reach maps colored by predicted value across the shared geographic bounding box.
- **Panel 3 (`Co-located Reach Difference Map`)**:
  - When comparing the **same determinand** (e.g., `ORN` vs. `EU-Hydro`), maps raw difference **`Δ = Run 2 - Run 1`** using a diverging **Blue–Grey–Red** palette (`Blue` = Run 2 lower, `Grey` = near-zero difference, `Red` = Run 2 higher).
  - When comparing **different determinands** (e.g., `Nitrate` vs. `Orthophosphate`), maps standardized Z-score divergence **`ΔZ = Z₂ - Z₁`** so you can immediately see which sub-catchments are disproportionately enriched in Determinand 2 (`Red`) vs. Determinand 1 (`Blue`).

### 4.5 Figure 5: `fig_comp5_determinand_relationships.png`
- **Panel A (`4-Quadrant Co-Occurrence Scatter`)**: Plots paired reaches colored by quadrant (`Red = Q1 High/High`, `Blue = Q2 Low/Low`, `Orange = Q3 High R1/Low R2`, `Purple = Q4 Low R1/High R2`), overlaid with median crosshairs, OLS regression line, Pearson `r`, Spearman `ρ`, and `% Concordance`.
- **Panel B (`Stratified Median Comparison`)**: Grouped bar chart comparing `Upland`, `Lowland Non-Aquifer`, `Lowland Aquifer`, `Headwaters (Ord 1–2)`, `Medium (Ord 3–4)`, and `Mainstem (Ord 5+)`. Dual `Y`-axes are automatically used when comparing determinands with different measurement scales.
- **Panel C (`Longitudinal Strahler Stream Order Profile`)**: Plots median reach predictions from `Order 1` to `Order 6+` for both runs, revealing whether both variables accumulate at the same rate downstream.
- **Panel D (`Molar N:P Stoichiometry` OR `ΔZ Divergence Histogram`)**:
  - For **Nitrogen vs. Phosphorus**, plots the histogram of `log10(Molar N:P)` across all reaches with shaded limitation zones (`Blue < 10:1 N-Limited`, `Green 10:1–20:1 Co-Limited / Redfield 16:1`, `Red/Orange > 20:1 P-Limited`).
  - For **all other comparisons**, plots the distribution of standardized differences `ΔZ = Z₂ - Z₁` with `±1 SD` thresholds.

---

## 5. Worked Empirical Summaries

### 5.1 Worked Example 1: Cross-Network Comparison — Nitrate (`0117`) on `ORN v2` vs. `EU-Hydro (Thames)`
*(Source report: `output/comparisons/det_0117_vs_det_0117_euhydro_thames/comparison_report.md`)*

1. **Out-of-Sample Cross-Validation & Head-to-Head Matched Evaluation (Sections 2 & 3)**:
   - Across all of England (`8,401` stations), **`ORN v2`** achieves **`Global CV R² = 0.5658`** (`RMSE = 0.5067`, `MAE = 0.3519`), whereas **`EU-Hydro (Thames)`** (`996` stations) achieves **`Global CV R² = 0.4533`** (`RMSE = 0.4730`, `MAE = 0.3125`).
   - When evaluated on the **exact same `M = 847` co-located Thames monitoring stations (`≤ 250 m`)**, both networks achieve virtually identical predictive accuracy:
     - **`ORN v2` (Matched Thames)**: `CV R² = 0.4494`, `RMSE = 0.4836`, `MAE = 0.3186`
     - **`EU-Hydro` (Matched Thames)**: `CV R² = 0.4460`, `RMSE = 0.4851`, `MAE = 0.3177`
   - **Formal Paired Tests (Section 4)**: The paired Diebold–Mariano `t`-test on squared errors (`t = -0.25`, `p = 0.805 ns`) and paired Wilcoxon signed-rank test on absolute errors (`V = 171,334`, `p = 0.246 ns`) confirm **no statistically significant difference in station prediction error** between `ORN v2` and `EU-Hydro` inside the Thames basin. However, **`ORN v2` achieves significantly lower Kriging uncertainty (`Med SE = 0.506` vs. `0.556`, Mann–Whitney `p < 1e-16`)** because national parameter estimation leverages `8.4×` more training stations and `EU-Hydro` contains thousands of unmonitored artificial canal/drainage segments (`19,642` reaches in Thames alone vs. `3,670` topological ORN reaches).
2. **Co-Located Reach Agreement (`N = 5,321` Matched Thames Reaches, Section 6)**:
   - Co-located reaches (`median separation = 41.5 m`) exhibit **very strong spatial correlation (`Pearson r = 0.8758`, `Spearman ρ = 0.8705`, `p < 1e-16`)** and a near-unity scaling slope (`EU-Hydro = 0.308 + 0.961 × ORN`).
   - **4-Quadrant Concordance is `88.7%`** (`Odds Ratio = 62.6`, `χ² = 3,189`, `p < 1e-16`), and **Upper-Quartile Hotspot Overlap (`Jaccard > P75`) is `77.1%`** (`P(EU-Hydro > P75 | ORN > P75) = 87.1%`).
3. **Parameter Stability Across Networks (Section 7)**:
   - Two-Sample Wald `Z`-tests show that key physical covariate slopes—including **`elevation` (`p = 0.654 ns`)**, **`slope` (`p = 0.552 ns`)**, **`arable` (`p = 0.720 ns`)**, **`woodland` (`p = 0.391 ns`)**, and **`stw_log` (`p = 0.460 ns`)**—are **statistically invariant** between the national `ORN v2` network and the `EU-Hydro` Thames network.

---

### 5.2 Worked Example 2: Cross-Determinand Comparison — Nitrate (`0117`) vs. Orthophosphate (`0180`) on `ORN v2`
*(Source report: `output/comparisons/det_0117_vs_det_0180_orn/comparison_report.md`)*

1. **Predictive Skill & Spatial Covariance Structure (Sections 2, 3 & 7)**:
   - **`Nitrate (0117)`** (`8,401` stations) achieves **`CV R² = 0.5658`** (`OLS R² = 0.4252`, **`ΔR² = +0.1406`**, `-13.1%` RMSE reduction from network kriging).
   - **`Orthophosphate (0180)`** (`8,358` stations) achieves **`CV R² = 0.4878`** (`OLS R² = 0.3401`, **`ΔR² = +0.1477`**, `-11.9%` RMSE reduction from network kriging).
   - **Why is `OLS R²` higher for Nitrate (`0.425`) than Orthophosphate (`0.340`), while `ΔR²` from spatial kriging is slightly larger for Orthophosphate (`+0.148`)?**
     Diffuse catchment land cover (`arable`, `improved_grass`) and bedrock hydrogeology (`aquifer` slope `+0.391` for Nitrate vs. `+0.172` for Orthophosphate, Wald `Z = 6.95`, `p < 1e-16`) explain a large share of Nitrate variance via `OLS`. Conversely, Orthophosphate is more heavily driven by point-source wastewater discharges (`stw_log` slope `+0.072` vs. `+0.044`, Wald `Z = -8.78`, `p < 1e-16`) and localized sediment interactions (`tau.sq = 0.251` vs. `0.155`), making longitudinal river-network kriging especially critical for capturing downstream sewage plumes.
2. **Distributional & Environmental Gradient Contrasts (Section 5)**:
   - Across all `37,396` ORN reaches, predicted **Nitrate** has a median of **`2.344 mg/L` (`IQR = 2.671 mg/L`, `CV = 73.2%`)**, whereas **Orthophosphate** is an order of magnitude lower in mass concentration (`Median = 0.058 mg/L`, `IQR = 0.068 mg/L`) but exhibits higher relative dispersion (**`CV = 95.8%`**, `Skewness = 3.04`, `Kurtosis = 19.65`) due to acute point-source spikes downstream of major `STWs`.
   - Along the **Strahler Stream Order gradient**, both nutrients enrich downstream from `Headwaters (Ord 1–2)` to `Mainstem Rivers (Ord 5+)`, with Nitrate increasing **`1.69x`** (`2.129` to `3.591 mg/L`) and Orthophosphate increasing **`1.77x`** (`0.052` to `0.092 mg/L`).
   - Across **Landscape Regimes**, however, **Principal Aquifers** amplify Lowland Nitrate strongly (`Lowland Aquifer = 4.30x` vs. `Upland`, compared with `Lowland Non-Aquifer = 2.79x`), whereas **Lowland Orthophosphate** is nearly identical between `Lowland Non-Aquifer` (`2.56x` vs. `Upland`, `0.069 mg/L`) and `Lowland Aquifer` (`2.51x` vs. `Upland`, `0.068 mg/L`).
3. **Reach-by-Reach Co-Occurrence & Hotspot Overlap (`N = 37,396` Reaches, Section 6)**:
   - Log-transformed reach predictions are **strongly positively correlated (`r_log = 0.6780`, `Spearman ρ = 0.6594`, `p < 1e-16`)**, following a sub-linear power-law scaling **`log10(P) = -1.425 + 0.695 × log10(N)`**.
   - **4-Quadrant Concordance is `73.2%`** (`36.6%` `Q1 High N / High P` co-enriched reaches and `36.6%` `Q2 Low N / Low P` co-clean reaches), yielding an **Odds Ratio of `7.42`** (`χ² = 7,968`, `p < 1e-16`).
   - The remaining **`26.8%` discordant reaches** partition cleanly into **`13.4%` High-Nitrate / Low-Phosphate (`Q3`)** reaches (predominantly rural Chalk/Limestone agricultural groundwater streams) and **`13.4%` Low-Nitrate / High-Phosphate (`Q4`)** reaches (clay/upland-margin streams with local sewage or livestock inputs).
   - **Upper-Quartile Multi-Pollutant Hotspot Overlap (`Jaccard > P75`) is `37.6%`**, meaning **`54.6%` of top-quartile Nitrate reaches are simultaneously top-quartile Orthophosphate hotspots**.
4. **Freshwater `N:P` Stoichiometry & Nutrient Limitation Regimes (Section 6)**:
   - Across all `37,396` reaches in England, the predicted **Median Mass `N:P` ratio is `38.6:1` (`IQR = 29.7`)**, corresponding to a **Median Atomic (Molar) `N:P` ratio of `85.5:1` (`IQR = 65.6`, `P05–P95 = 24.4:1 to 217.0:1`)**—more than five times the canonical Redfield ratio (`16:1`).
   - **Stoichiometric Regime Breakdown across England**:
     - **`P-Limited (Molar N:P > 20:1)`**: **`97.0%` of reaches (`36,277 / 37,396`)** have abundant excess dissolved inorganic nitrogen relative to orthophosphate.
     - **`Co-Limited / Balanced (Molar N:P 10:1 – 20:1)`**: **`2.7%` of reaches (`1,021`)**, primarily in upland-margin catchments.
     - **`N-Limited (Molar N:P < 10:1)`**: **`0.3%` of reaches (`98`)**, restricted to low-nitrate headwaters receiving point-source phosphorus.
   - **Management Conclusion**: Because `97.0%` of English river reaches operate in a strongly **Phosphorus-limited stoichiometric regime (`Molar N:P > 20:1`)**, controlling in-stream benthic algal blooms requires **phosphorus stripping at `STWs` and diffuse agricultural `P` mitigation**, whereas nitrate reduction remains essential for **drinking-water compliance (`50 mg/L NO₃` / `11.3 mg/L NO₃-N`) and downstream coastal/estuarine nitrogen sensitivity**.

---

## 6. Related Documentation

- **S3N Modelling & Universal Stream Kriging Theoretical Guide**: [`./s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md)
- **Single-Output Descriptive Statistics Interpretation Guide**: [`./descriptive_stats_interpretation_guide.md`](./descriptive_stats_interpretation_guide.md)
- **End-to-End Pipeline Usage Guide**: [`./usage_guide.md`](./usage_guide.md)
- **Custom Water Network Guide (`BYOWN`)**: [`./custom_water_network_guide.md`](./custom_water_network_guide.md)
- **Preprocessing & Graph Topology Guide**: [`./preprocessing_pipeline.md`](./preprocessing_pipeline.md)
- **ORN Variables & Covariates Reference**: [`./orn_variables_reference.md`](./orn_variables_reference.md)
