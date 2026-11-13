# Non-Linear Spline (`GAM + S3N`) & Tree-Ensemble (`RF / GBM + S3N`) Regression Kriging Guide

> **Quick Summary**: The `--trend=linear|spline|gam|rf|gbm` option replaces the linear fixed-effects mean $x(s)^\top\beta$ with penalized splines, Generalized Additive Models, Random Forests, or Gradient Boosted Regression Trees prior to S3N stream-network kriging. It works directly on **standard raw/preprocessed data**, requires **zero external CRAN packages**, and can be freely combined with `--kernel=`, `--combine=`, `--afv=`, and `--meas-error`.

---

## 1. Key Distinction: `--trend=gam` vs. Two-Stage `meas_method="gam"`

Because both features involve Generalized Additive Models, it is important not to confuse them:

| Feature | CLI / Config Trigger | Data Requirement | What It Models | Where It Acts in S3N |
| :--- | :--- | :--- | :--- | :--- |
| **Non-Linear Spatial Trend (`GAM + S3N`, `RF + S3N`, `GBM + S3N`)** | `--trend=gam` (or `spline`, `rf`, `gbm`) | **Standard raw/preprocessed data** (e.g., `orn_rsn`, `orn_rss`, or custom CSVs) | **Non-linear spatial covariate effects & interactions** across the catchment landscape (Elevation, Slope, US_Accum, Geology, Easting/Northing) | Mean trend $\mu(s_i)$ |
| **Two-Stage GAM Measurement-Error Kriging** | `meas_method="gam"` or `--meas-error` | Preprocessed dataset containing site-level `meas_uncertainty` | **Temporal / site-level observation uncertainty** ($v_{\text{meas}, i}$) estimated from stage-1 site time-series | Nugget / covariance diagonal ($\tau^2_{\text{micro}} + v_{\text{meas}, i}$) |

**Can they be used together?**  
Yes. You can run `--trend=gam` (or `rf` / `gbm`) on **standard raw data** without `--meas-error`, **or** you can combine `--trend=gam --meas-error` on a Two-Stage GAM dataset so that stage-1 temporal uncertainty is propagated into the spatial nugget while stage-2 landscape covariates are modeled non-linearly.

See [`./two_stage_gam_measurement_error_kriging_guide.md`](./two_stage_gam_measurement_error_kriging_guide.md) for details on stage-1 measurement-error propagation.

---

## 2. Package Requirements (Zero Extra Installs Needed)

None of the `--trend=` engines require third-party CRAN packages such as `randomForest`, `ranger`, `xgboost`, or `gbm`. Every mode uses packages that ship with a standard R installation:

| Trend Mode | Engine Used | Package Status in R | Automatic Fallback |
| :--- | :--- | :--- | :--- |
| **`linear`** *(default)* | QR / GLS linear regression | **Base R** (`stats`) | — |
| **`spline`** | Orthogonalized natural cubic B-splines (`splines::ns`) + analytic LOO-GCV ridge selection (`chol2inv`) | **Base R** (`splines`) | Falls back to `linear` if `n_obs < 12` |
| **`gam`** | REML-penalized cubic regression splines (`mgcv::gam(..., bs='cr', method='REML')`) | **Recommended** (`mgcv`, bundled with standard R) | Automatically falls back to `spline` if `mgcv` is unavailable or `n_obs` is too small for REML smooths |
| **`rf`** | Bagged random-subspace decision-tree forest (`rpart::rpart`, `mtry = 0.7p`) with Out-of-Bag (`OOB`) evaluation | **Recommended** (`rpart`, bundled with standard R) | Falls back to `linear` if `n_obs < 12` |
| **`gbm`** | Stochastic gradient-boosted trees (`rpart::rpart`, `maxdepth = 3`, $\nu = 0.06$, `0.8` subsampling) with 5-fold cross-fitting | **Recommended** (`rpart`, bundled with standard R) | Falls back to `linear` if `n_obs < 12` |

---

## 3. Supported Trend Modes & Tuning Flags

Run `Rscript scripts/run_s3n_determinand.R --list-trends` at any time to print the built-in catalog:

| Mode (`--trend=`) | Shorthand Flag | Tuning Flags | Mathematical Formulation | Best Use Case |
| :--- | :--- | :--- | :--- | :--- |
| **`linear`** | *(default)* | — | $\mu(s) = x(s)^\top \beta$ | Fast baseline; linear/monotonic relationships on the transformed scale; clean parametric slopes. |
| **`spline`** | `--spline` | `--spline-df=3` *(range `2..10`)* | $\mu(s) = x(s)^\top\beta + \sum_j \tilde{B}_j(x_j)^\top \gamma_j$, penalized by LOO-GCV $\lambda$ | Smooth univariate non-linearities with lightweight matrix algebra and guaranteed stability on small basins. |
| **`gam`** | `--gam-s3n` | `--spline-df=3` *(sets basis $k = \text{df} + 2$)* | $\mu(s) = \beta_0 + \sum_{j \in \text{cont}} s_j(x_j) + \sum_{l \in \text{disc}} \beta_l x_l$ via REML | Data-adaptive smooth curves (e.g., upland-to-lowland thresholds, dilution plateaus) with automatic EDF selection. |
| **`rf`** | `--rf`, `--rf-s3n` | `--trees=100` *(range `10..2000`)* | $\hat{f}_{\text{RF}}(s) = \frac{1}{B}\sum_{b=1}^B T_b(x(s))$ (`mtry = 0.7p`, `maxdepth = 6`) | Strong **covariate interactions** (e.g., Chalk Aquifer × Low Elevation × Eastern Region) and non-smooth thresholds. |
| **`gbm`** | `--gbm`, `--gbm-s3n` | `--trees=100` *(capped at `250` rounds)* | $\hat{f}_{\text{GBM}}(s) = \bar{y} + \sum_{m=1}^{M^*} 0.06 \, T_m(x(s))$ (`maxdepth = 3`) | High predictive accuracy via sequential boosting with automatic 5-fold early stopping ($M^*$). |

---

## 4. How Orthogonal Out-of-Fold Residual Kriging Works

A naive "fit Random Forest/GAM, then krige the in-sample residuals" approach fails in spatial statistics for two reasons:
1. **In-sample variance collapse**: A flexible tree ensemble can fit observed training stations so closely that in-sample residuals $y_i - \hat{f}(x_i)$ approach zero, tricking NNGP into estimating $\hat{\sigma}^2 \approx 0$ and disabling spatial kriging on unmonitored reaches.
2. **Spatial confounding / loss of linear inference**: Downstream plotting ([`../scripts/03_generate_plots.R`](../scripts/03_generate_plots.R)) and cross-run Wald tests ([`../scripts/04_compare_model_runs.R`](../scripts/04_compare_model_runs.R)) expect fixed-effect slopes $\hat{\beta}_{\text{GLS}}$ and 95% confidence intervals.

`fit_nonlinear_trend()` in [`../scripts/s3n_model_functions.R`](../scripts/s3n_model_functions.R) solves both problems cleanly:

1. **Strictly Out-of-Sample Trend Values at Observed Stations ($\hat{f}_{\text{OOF}}$)**:
   - **`rf`**: Each station's trend prediction $\hat{f}_{\text{OOF}}(s_i)$ is averaged **only over bootstrap trees that did not contain station $i$** (Out-of-Bag).
   - **`gbm`**: Uses **5-fold cross-fitting** on the training stations so each station is predicted by a boosted ensemble that never saw it, while also selecting the optimal stopping iteration $M^*$ by minimizing OOF MSE.
   - **`gam` / `spline`**: Applies exact Leave-One-Out hat-matrix adjustment:
     $$\hat{f}_{\text{LOO}}(s_i) = y_i - \frac{y_i - \hat{f}_{\text{in}}(s_i)}{1 - h_{ii}}$$
   - **Range Clamping for `gam`**: Evaluation covariates are clamped to `[min - 0.25*range, max + 0.25*range]` of the training range so cubic polynomials cannot explode at extreme unmonitored headwaters.
2. **Orthogonal Projection Against Linear Design Matrix $X$**:
   - Projects $\hat{f}$ onto the column space of $X_{\text{obs}}$ via $\hat{\beta}_{\text{proj}} = (X_{\text{obs}}^\top X_{\text{obs}})^{-1} X_{\text{obs}}^\top \hat{f}_{\text{OOF}}$ to isolate the **purely non-linear / interaction component**:
     $$\hat{g}_{\perp}(s) = \hat{f}(s) - X(s)\hat{\beta}_{\text{proj}}, \qquad X_{\text{obs}}^\top \hat{g}_{\perp}(s_{\text{obs}}) = 0$$
3. **S3N Fitting & Kriging on $Y^*(s) = Y(s) - \hat{g}_{\perp}(s)$**:
   - Fits S3N `Tail-Up`, `Tail-Down`, and `Joint (TU + TD + Euc)` to $Y^*(s_{\text{obs}})$ with design matrix $X_{\text{obs}}$, then adds $\hat{g}_{\perp}(s_0)$ back to every prediction:
     $$\hat{Y}(s_0) = \underbrace{x(s_0)^\top \hat{\beta}_{\text{GLS}} + \hat{g}_{\perp}(s_0)}_{\text{Non-Linear / Ensemble Trend } \hat{\mu}(s_0)} + \sum_{k=1}^m \lambda_k(s_0) \left( Y(s_{j_k}) - \hat{\mu}_{\text{OOF}}(s_{j_k}) \right)$$
4. **Zero-Leakage 5-Fold Cross-Validation**:
   - Inside each CV fold $k \in \{1,\dots,5\}$, `fit_nonlinear_trend()` is re-fitted **exclusively on the 80% training stations**, predicting $\hat{g}_{\perp}$ onto the 20% held-out test stations before fold-specific S3N kriging.

---

## 5. Practical Usage Examples

```bash
# 1. Inspect available trend modes and tuning options
Rscript scripts/run_s3n_determinand.R --list-trends

# 2. GAM + S3N on standard RSN Nitrate (0117) with sqrt transform
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=gam

# 3. Natural Cubic Splines + S3N with 4 degrees of freedom per continuous covariate
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=spline --spline-df=4

# 4. Random Forest Regression Kriging (150 trees) + S3N
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=rf --trees=150

# 5. Gradient Boosted Trees + S3N + Auto-Kernel + Hydrological-Yield AFV
Rscript scripts/run_s3n_determinand.R 0117 --transform=sqrt --trend=gbm --trees=100 --kernel=auto --afv=hydro_yield

# 6. Combine GAM + S3N with Two-Stage Measurement Error Kriging (optional)
Rscript scripts/run_s3n_determinand.R 0117 --source=rsn_gam_mean --trend=gam --meas-error

# 7. Compare Linear S3N vs. GAM + S3N vs. RF + S3N in one report
Rscript scripts/04_compare_model_runs.R \
  output/results/det_0117_orn_rsn_trans_sqrt \
  output/results/det_0117_orn_rsn_trans_sqrt_trend_gam \
  output/results/det_0117_orn_rsn_trans_sqrt_trend_rf \
  --labels="Linear S3N,GAM + S3N,RF + S3N"
```

---

## 6. Interpreting Outputs & Diagnostics

1. **Output Directory Tagging**:
   Non-default trend runs are automatically saved to `output/results/<run>_trend_<mode>/` (e.g., `det_0117_orn_rsn_trans_sqrt_trend_gam`), keeping linear and non-linear experiments cleanly separated.
2. **Predictions CSV (`predictions_*.csv`)**:
   - **`pred_trend_sc`**: Exports the non-linear/ensemble mean trend $x(s_0)^\top \hat{\beta}_{\text{GLS}} + \hat{g}_{\perp}(s_0)$ on the model scale, allowing you to separate **landscape-driven variation** from **spatial stream-network kriging adjustments** (`pred_model_sc - pred_trend_sc`).
3. **Cross-Validation Summary (`cv_comparison_*.csv`)**:
   - The `Spatial_Method` column is tagged with `[SPLINE]`, `[GAM]`, `[RF]`, or `[GBM]` (e.g., `Standard S3N [GAM]` or `Measurement Error Kriging [RF]`).
   - Compare `RMSE`, `MAE`, and `R2` of `S3N Combined (Tail-Up + Tail-Down)` against `OLS (Non-Spatial)` and against `Trend_CV_R2` to quantify how much accuracy comes from non-linear covariates vs. dendritic spatial correlation.
4. **Parameters CSV (`parameters_*.csv`)**:
   - Retains `S3N Tail-Up (FC + AFV)`, `S3N Tail-Down (Stream)`, and `S3N Joint (TU + TD + Euc)` parameter rows so `03_generate_plots.R`, `04_compare_model_runs.R`, and `05_descriptive_stats.R` work without modification.
   - Appends a `S3N Trend (...)` diagnostic block:
     - **`spline`**: `Spline_NL:<var>` (L2 norm of orthogonal non-linear spline coefficients) and `Spline_Total_EDF`.
     - **`gam`**: `GAM_EDF:<var>` (Effective Degrees of Freedom; `≈ 1.0` indicates a linear relationship, `> 2.0` indicates strong curvature).
     - **`rf` / `gbm`**: `RF_VarImp(%):<var>` or `GBM_VarImp(%):<var>` (normalized variable importance summing to `100%`) plus `GBM_Best_Trees` (optimal early-stopping round $M^*$).
     - **All modes**: `Trend_CV_RMSE` and `Trend_CV_R2` (5-fold out-of-sample skill of the non-spatial trend model alone).
