# Multivariate Stream Co-Kriging — Core Architecture, Factor-LMC Theory & Usage Guide

**Multivariate Stream Co-Kriging** jointly models $P \ge 2$ water-quality, chemical, or ecological variables across a dendritic river network using a single **Scalable Factor Linear Model of Coregionalization (`Factor-LMC`) Joint Nearest-Neighbor Gaussian Process (`NNGP`)**. It runs through a dedicated CLI entrypoint, [`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R), leaving all single-determinand scripts ([`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R), [`scripts/02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R)) completely unchanged.

---

## 1. Why Multivariate Stream Co-Kriging?

Single-determinand S3N models each variable completely in isolation. Multivariate Stream Co-Kriging fits $P \ge 2$ variables together inside a unified C++/OpenMP solver ([`src/s3n_fast_neighbors.cpp`](../src/s3n_fast_neighbors.cpp), `C_s3n_multivariate_nngp`) and modular R engine ([`scripts/s3n_multivariate_functions.R`](../scripts/s3n_multivariate_functions.R)) so **variables actively borrow spatial strength from one another** while writing standard `predictions_<domain>_<code>.csv` files compatible with downstream tools.

### How Do Variables Help Predict Each Other?

The Factor-LMC NNGP learns three distinct residual correlation pathways across your determinands (after accounting for shared catchment covariates $X$):

1. **Co-Located Local Surrogate (`Corr_Nugget`)**: If Reach $A$ is missing a target variable (for instance, Orthophosphate `0180` was not sampled or failed QA), yet has a co-located observation of Nitrate (`0117`), Conductivity (`0077`), or macroinvertebrate `whpt_aspt`, co-kriging directly adjusts the target prediction at Reach $A$ using the co-located auxiliary residual.
2. **Flow-Connected Downstream Plumes (`Corr_TU`)**: If an upstream tributary station measured an unexplained spike in Nitrate above what elevation, arable %, and geology predicted, `Corr_TU` propagates that residual signal downstream along `Tail-Up` AFV flow paths to refine predictions of *both* Nitrate and coupled Orthophosphate.
3. **Adjacent Watershed / Headwater Regime (`Corr_TD`)**: Across nearby unconnected tributaries (`Tail-Down`), shared local pressures—such as septic tank density, livestock intensity, or small unmapped point discharges—allow nearby stations observing one variable to constrain predictions and tighten standard errors (`se_model_sc`) of the others.

### Example: Why Model Nitrate (`0117`) and Orthophosphate (`0180`) Together Rather Than Separately?

Even though catchment covariates explain part of both Nitrate and Orthophosphate variation, substantial **unmeasured residuals** remain:
- **Unmapped point and diffuse sources**: A leaky slurry store, agricultural field drain, septic tank cluster, or small Sewage Treatment Works (STW) rarely appears directly in national GIS covariates. Modelling Nitrate and Orthophosphate together lets an observed nutrient spike in one solute inform the spatial surface of the other across connected and neighboring reaches.
- **Uneven ("heterotopic") monitoring coverage**: When combining networks (`--source=rss` or `--custom-samples`), slicing by season/year, or dealing with censored/missed lab runs, some reaches have Nitrate but lack Orthophosphate (or vice versa). Co-kriging natively masks `NA`s row-by-row rather than dropping incomplete stations.
- **Lower uncertainty & clearer WFD tiering**: Borrowing strength across correlated determinands shrinks `se_model_sc`, moving borderline reaches out of the `"Uncertain / Priority Monitoring"` tier when combined with `--threshold=uk`.

> [!NOTE]
> **When Co-Kriging Adds Little Over Univariate S3N:** If off-diagonal residual correlations are near zero once covariates are subtracted (`|Corr_TU|`, `|Corr_TD|`, `|Corr_Nugget|` all $< 0.15$), or if every variable was measured at the exact same set of reaches (*isotopic* design), Co-Kriging gives roughly the same $R^2$ as Univariate S3N. Every run automatically generates **`cokriging_evaluation_report.md`** and **`cv_univariate_vs_cokriging.csv`** so you can verify immediately whether joint modelling helped.

### Key Use Cases

| Use Case | Example `--targets=` | Why Joint Modelling Helps |
| :--- | :--- | :--- |
| **Coupled Nutrient & Eutrophication Mapping** | `--targets=0117,0180,0111 --threshold=uk` | Exploits co-movement of Nitrate, Orthophosphate, and Ammonia to track agricultural vs. sewage plumes and sharpen UK EQS classifications. |
| **Dense Sonde / Easy Proxy $\rightarrow$ Sparse Lab Target** | `--targets=0077,6396,0180` or `--targets=wq:0077,wq:0301,lcms:<cas>` | Uses high-coverage Conductivity (`0077`), Turbidity (`6396`), or DOC (`0301`) to interpolate sparser nutrients or trace LCMS organics. |
| **Chemistry + Ecology Cross-Trophic Synthesis** | `--targets=inverts:whpt_aspt,wq:0117,wq:0180,wq:9924` | Biological indices (`whpt_aspt`, `tdi4`, `rmni`) integrate chronic multi-month water quality and oxygenation stress where spot chemistry was missed—and chemistry informs unmonitored biology reaches. |
| **Multi-Network / Citizen-Science Fusion** | `--targets=0117,0180 --source=rsn+my_study` | Fuses statutory EA monitoring (full suites) with local campaigns or citizen-science monitors that sampled only a subset of determinands. |
| **Multi-Stressor Syndrome Screening ($P \ge 5$)** | `--targets=0117,0180,0111,0301,0077,inverts:whpt_aspt --factors=3` | Summarizes 5–10+ correlated stressors via $K=3$ latent Tail-Up/Tail-Down factors without overfitting. |

---

## 2. Mathematical Model: Scalable Factor-LMC Joint Stream NNGP

For $P$ response variables $y(s) = (y_1(s), \dots, y_P(s))^\top$ observed (possibly partially/heterotopically) along river reaches $s \in \mathcal{S}$, the standardized model-scale spatial process is:

$$y(s) = B^\top x(s) + w_{\text{TU}}(s) + w_{\text{TD}}(s) + \varepsilon(s)$$

where:
- $x(s) \in \mathbb{R}^p$ is the reach covariate vector and $B \in \mathbb{R}^{p \times P}$ is the multivariate regression coefficient matrix estimated via **exact Multivariate Spatial Generalized Least Squares (GLS)** inside the C++ NNGP likelihood.
- $w_{\text{TU}}(s)$ and $w_{\text{TD}}(s)$ follow a **Linear Model of Coregionalization (LMC)** with Tail-Up spatial correlation $C_{\text{TU}}(h; \phi_{\text{tu}})$ (weighted by Additive Function Values `AFV`) and Tail-Down stream correlation $C_{\text{TD}}(h; \phi_{\text{td}})$:

$$\text{Cov}\big(y(s_i), y(s_j)\big) = A_{\text{tu}} C_{\text{TU}}(s_i, s_j; \phi_{\text{tu}}) + A_{\text{td}} C_{\text{TD}}(s_i, s_j; \phi_{\text{td}}) + \delta_{ij} \Psi$$

### How $P = 2\text{–}4$ vs. $P = 5+$ Variables Scale (`--factors=K`)

Each $P \times P$ cross-covariance matrix is parameterized via **Factor Cholesky loadings** $A = L L^\top$, where $L \in \mathbb{R}^{P \times K}$ has rank $K \le P$:

| Regime | Latent Factors $K$ | Behavior | Best Use Case |
| :--- | :---: | :--- | :--- |
| **Full-Rank LMC** | $K = P$ | Unconstrained positive-semidefinite cross-covariance across all pairs ($P = 2, 3, 4$) | Core 2–4 variable pairs (e.g. `0117,0180` or `whpt_aspt,0117,0180,9924`) |
| **Latent-Factor LMC** *(default `K = 3`)* | $K < P$ (e.g. $K = 3$ for $P = 6\text{–}10$) | Captures the top $K$ shared hydrochemical/ecological watershed gradients + variable-specific diagonal nugget $\Psi$ | Large 5–10+ multi-pollutant or cross-taxon syntheses without quadratic parameter explosion |
| **Sequential Residual (`--method=colocated`)** | — | Fast two-stage kriging: fits univariate S3N per target, then blends co-located cross-residual OLS predictions across sites | Rapid exploratory screening across many variables ($P \ge 8$) |

All NNGP conditional mean and block-Cholesky operations in `C_s3n_multivariate_nngp` dynamically mask missing entries (`NA`/`NaN`) at heterotopic stations, so a reach sampled for Nitrate & Conductivity but missing Orthophosphate contributes exact partial information to the joint likelihood and receives full cross-aided predictions.

---

## 3. CLI Reference (`scripts/run_s3n_cokriging.R`)

```bash
Rscript scripts/run_s3n_cokriging.R --targets=<t1>,<t2>[,<t3>,...] [OPTIONS]
```

### Target & Transform Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--targets=<t1>,<t2>,...` | *(Required; or 1st positional arg)* | Comma-separated $P \ge 2$ targets. Prefix with dataset when mixing datasets or avoiding ambiguity: `wq:0117`, `inverts:whpt_aspt`, `diatoms:tdi4`, `macrophytes:rmni`, `rhs:hms`, `lcms:<cas>`. |
| `--dataset=<ds>` | Auto-detected | Fallback dataset name (`wq`, `inverts`, `diatoms`, `macrophytes`, `rhs`, `lcms`) for unprefixed target tokens. |
| `--transforms=<t1,...>` | `log` (`wq`/`lcms`), `none` (ecology) | Comma-separated per-target transforms matching length $P$ (or a single transform broadcast to all targets), e.g. `--transforms=none,log,log,none`. |
| `--censor=ros\|half_lod` | `ros` | Left-censored (`<LOD`) imputation applied consistently per target before transformation. |

### Co-Kriging & Spatial Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--method=lmc\|colocated` | `lmc` | `lmc` fits the joint Factor-LMC NNGP in C++/OpenMP; `colocated` runs fast sequential univariate S3N + co-located residual conditioning. |
| `--factors=<K>` | `3` | Maximum latent spatial factors $K \le P$ for Factor-LMC (`K = P` gives exact full-rank LMC when $P \le 3$; set `--factors=4` for exact full-rank 4-variable LMC). |
| `--neighbors=<m>` | `10` | Nearest stream neighbors per Tail-Up / Tail-Down NNGP DAG. |
| `--kernel=<name>` | `exponential` | Spatial stream correlation kernel (`exponential`, `spherical`, `mariah`, `linear`). |
| `--vars=<v1,v2,...>` | `default` | Shared catchment fixed-effect predictors across $X$ (run `--list-vars`). |
| `--cv-mode=<mode>` | `random` | Out-of-sample CV splitting scheme (`random`, `catchment`, `strahler`, `spatial_block`; run `--list-cv`). |
| `--threshold=uk\|<vals>` | *(None)* | Regulatory exceedance evaluation: pass `uk` for automatic UK DWI/WFD thresholds per target (where defined), or comma-separated raw thresholds aligned with `--targets`. Combine with `--threshold-dir=above\|below,...`. |
| `--basin=<spec>` | *(None)* | Custom river basin(s) or CatchID(s) to fit alongside England, Humber, Thames, and Severn. |
| `--source=<src>` | `rsn` | Monitoring source (`rsn`, `ssn`, `rss`, or custom profiles). |
| `--year` / `--season` / `--month` | *(All)* | Restrict all targets to a matching survey year, season, or calendar month. |
| `--omp=<N>` | `1` | OpenMP threads used by `C_s3n_multivariate_nngp`. |
| `--fig-format=<png\|pdf\|svg>` | `png` | Output figure format (`png`, vector `pdf`, or vector `svg`). |
| `--dpi=<N>` | `200` / `220` | Output PNG resolution in DPI (pixel dimensions scale proportionally). |
| `--fig-size=<pct>` | `100%` | Scale figure physical dimensions (e.g. `125%`, `+25%`, `-20%`, or `1.25`). |

---

## 4. Practical Examples

### Example 1: Coupled Nutrient Pair (Nitrate + Orthophosphate)
```bash
Rscript scripts/run_s3n_cokriging.R --targets=0117,0180 --threshold=uk --omp=4
```
Fits a full-rank 2-variable Tail-Up + Tail-Down LMC across England + Humber, Thames, and Severn, outputs cross-correlations (`Corr_TU`, `Corr_TD`, `Corr_Nugget`), and checks exceedance against UK Nitrate (`11.3 mg/L`) and WFD Phosphate benchmarks.

### Example 2: Cross-Dataset Water Quality & Ecology Co-Kriging ($P = 4$)
```bash
Rscript scripts/run_s3n_cokriging.R \
  --targets=inverts:whpt_aspt,wq:0117,wq:0180,wq:9924 \
  --transforms=none,log,log,none \
  --factors=4 --cv-mode=spatial_block --omp=4
```
Jointly models macroinvertebrate water quality sensitivity (`whpt_aspt`), Nitrate (`0117`), Orthophosphate (`0180`), and Dissolved Oxygen (`9924`) across partially overlapping monitoring stations, testing spatial transferability with spatially blocked 5-fold CV.

### Example 3: Scalable 6-Variable Catchment Syndrome Synthesis ($K = 3$ Latent Factors)
```bash
Rscript scripts/run_s3n_cokriging.R \
  --targets=wq:0117,wq:0180,wq:0111,wq:0301,wq:0077,inverts:whpt_aspt \
  --transforms=log,log,log,log,log,none \
  --factors=3 --neighbors=10 --omp=8
```
Uses $K = 3$ latent Tail-Up/Tail-Down factors to summarize nutrient enrichment, organic matter, ionic strength, and biological response across 6 heterotopic datasets.

### Example 4: Fast Exploratory Screening (`--method=colocated`)
```bash
Rscript scripts/run_s3n_cokriging.R \
  --targets=0117,0180,0111,0301,6396 --method=colocated --basin=Eden
```

---

## 5. Dual Cross-Validation & Output Interpretation

Every run evaluates **three comparable out-of-sample CV accuracy benchmarks** in `cv_univariate_vs_cokriging.csv` (and Panel B of `fig_cokrig_cross_correlation.png`):

1. **`Univariate S3N`**: Independent 5-fold CV fit using only that single target determinand and covariates.
2. **`Co-Kriging (All-Holdout)`**: 5-fold CV where **all $P$ variables** at every holdout station are masked simultaneously (`Y_test = NA`). Tests whether joint multivariate covariance learning improves predictions at completely unmonitored stream reaches.
3. **`Co-Kriging (Aux-Aided)`**: Exact Leave-One-Target-Out conditioning where target $p$ is masked at test station $i$, while any co-located auxiliary targets $q \ne p$ observed at station $i$ (plus all $m$ stream-network neighbors) remain observed. Measures the exact accuracy gain from leveraging cheaper/denser co-located proxies to fill gaps in harder-to-measure targets.

### Output Directory & Artifacts (`output/results/cokrig_<targets>_<net>_<src>/`)

| File | Description |
| :--- | :--- |
| `cokriging_evaluation_report.md` | **Automated Value-Add Report**: Explicitly tells you whether Co-Kriging was worth doing (`Strong Benefit`, `Moderate Benefit`, `Marginal/Parity`, or `No Gain`), **where** it helped (best targets, river basins, heterotopic gap-fill counts, and WFD tier clarity), and **how** (`Tail-Up` vs `Tail-Down` vs `Nugget` cross-correlations). |
| `cv_univariate_vs_cokriging.csv` | Paired RMSE, MAE, and $R^2$ comparing `Univariate S3N`, `Co-Kriging (All-Holdout)`, and `Co-Kriging (Aux-Aided)` across every domain and target. |
| `cokriging_cross_covariance.csv` | Decomposed Tail-Up (`Cov_TU`, `Corr_TU`), Tail-Down (`Cov_TD`, `Corr_TD`), Co-Located Nugget (`Cov_Nugget`, `Corr_Nugget`), and Total cross-covariances/correlations on each variable's model scale. |
| `predictions_<domain>_<code>.csv` | Standard reach-level prediction table (`RivID`, `pred_combined`, `pred_model_sc`, `se_model_sc`, optional `exceed_prob` & `compliance_tier`) for each target and domain. |
| `exceedance_summary_cokriging.csv` | Tier counts & percentages across all targets when `--threshold=` is specified. |
| `fig_cokrig_cross_correlation.png` | **Diagnostic Summary**: Top row displays decomposed Tail-Up, Tail-Down, and Nugget cross-correlation heatmaps; bottom panel compares Univariate vs. All-Holdout vs. Aux-Aided CV $R^2$ and RMSE reduction. |
| `fig_cokrig_england_<code>.png` | **Figure 1 (per target)**: National two-panel prediction (`pred_model_sc`) & Co-Kriging SE (`se_model_sc`) map (plus `fig_exceedance_england_<code>.png` when thresholds are active). |

---

## 6. How to Interpret the Co-Kriging Figures

### Interpreting `fig_cokrig_england_<code>.png` ("Figure 1" National Map)

Each target determinand produces a national two-panel Figure 1 (`fig_cokrig_england_<code>.png`) using [`plot_fig1_england_map()`](../scripts/s3n_plot_functions.R):

- **Panel A — Co-Kriged Spatial Surface (`S3N Predicted <Determinand> Across England`)**:
  - **Line color** shows the predicted concentration or ecological score along every river reach (`90,876` reaches on default ORN v2 RSN), scaled using the determinand's model transformation (`log`, `none`, etc.) with raw-unit tick marks on the colorbar.
  - By default, **cool blue = lower values** and **warm orange/red = higher values** (e.g. elevated Nitrate or Orthophosphate). If `--reverse-color` is passed—recommended for ecological quality metrics such as `whpt_aspt` or `9924` Dissolved Oxygen—low/degraded values appear red and high/healthy values appear blue.
  - **Line width** scales with **Shreve stream order**, making mainstem rivers bold and headwaters fine so downstream plume accumulation stands out visually.
  - Compared with a univariate `fig1` map, reaches lacking local observations of that target may show sharper sub-catchment plumes or finer mainstem transitions where auxiliary variables helped infer the unmeasured residual.

- **Panel B — Co-Kriging Uncertainty & Monitoring Stations (`Prediction Uncertainty & Stations`)**:
  - **Line color** displays the **Co-Kriging Standard Error (`se_model_sc`)** on the model scale (`log` units or raw units). **Darker/cooler/lower SE** indicates high prediction confidence; **brighter/higher SE** highlights poorly constrained reaches far from informative stations.
  - **Circular points** mark reaches where **this specific target** has direct observations (`has_obs == 1`), colored by the observed reach mean (`obs_value`) on the same palette as Panel A.
  - **What to look for in Co-Kriging**: Even along reaches **without a station dot for Target $p$**, Panel B's Co-Kriging SE can still be visibly lower than in a univariate run whenever another correlated target (shown on its own `fig_cokrig_england_*.png` map) was observed nearby or co-located on that reach (`Aux-Aided` uncertainty reduction).

### Interpreting `fig_cokrig_cross_correlation.png`

- **Panels A1–A3 (Top Row — Cross-Correlation Matrices)**:
  - Symmetric heatmaps (`-1.0` deep blue through `0.0` neutral white to `+1.0` deep red) showing *how* variables co-vary after fixed catchment effects are removed:
    - **A1 (`Tail-Up`)**: Cross-correlation travelling downstream along flow-connected paths. High magnitude indicates unmeasured tributary plumes moving together.
    - **A2 (`Tail-Down`)**: Cross-correlation among nearby reaches regardless of flow direction. High magnitude indicates shared sub-catchment land-use or hydrogeological regimes.
    - **A3 (`Co-Located Nugget`)**: Instantaneous correlation at the same reach/station. High magnitude means measuring Variable B at a site strongly predicts missing Variable A there (`Aux-Aided` power).
- **Panel B (Bottom Row — Dual Cross-Validation Bar Chart)**:
  - Side-by-side comparison of **`Univariate S3N`** (grey) vs. **`Co-Kriging (All-Holdout)`** (blue) vs. **`Co-Kriging (Aux-Aided)`** (orange) CV $R^2$, annotated with the percentage **Aux-Aided RMSE change** below each target name. A tall orange bar relative to grey confirms strong surrogate/gap-filling value; a taller blue bar confirms better extrapolation even to completely unmonitored reaches.

---

## Related Guides
- [`docs/s3n_modelling_and_kriging_guide.md`](./s3n_modelling_and_kriging_guide.md) — Univariate S3N NNGP theory, multi-kernel options, CV modes, and optimal design.
- [`docs/usage_guide.md`](./usage_guide.md) — Full CLI flag reference for single-determinand modeling.
- [`docs/regulatory_exceedance_mapping_guide.md`](./regulatory_exceedance_mapping_guide.md) — UK/WFD exceedance threshold presets and 5-tier management mapping.
