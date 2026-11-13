# Multi-Source Fusion & Citizen-Science Co-Kriging Guide

> **Supplementing Statutory River Monitoring Networks (`RSN` / `SSN`) with Citizen Science, Telemetry Sensors, and Local Campaigns**

Statutory monitoring networks such as the Environment Agency **River Surveillance Network (`RSN`)** and **Small Streams Network (`SSN`)** provide high-precision, laboratory-accredited analytical measurements, yet they leave vast stretches of tributaries, headwaters, and local catchments unsampled in any given season or year. Conversely, **citizen-science campaigns** (e.g., WaterBlitz, Anglers' Riverfly Monitoring Initiative), catchment-partnership grab samples, and continuous sondes offer dense spatial or temporal coverage, often with coarser precision, colorimetric range kits, or different measurement protocols.

This guide explains **when** to keep networks as separate co-kriged variables versus **when** to pool them with explicit source priority and precision weighting, documents every multi-source CLI option (`--prefer-source`, `--source-weight`, `--source-se`, `--method=lmc|colocated`), and walks through a complete end-to-end workflow.

---

## 1. Why Co-Kriging Beats Naive Data Pooling

Simply appending citizen-science or campaign readings into the same target variable as `RSN` laboratory assays causes four structural problems in spatial stream modelling:

1. **Methodological Bias & Scale Offset**: Colorimetric dip-strips, handheld probes, or semi-quantitative scorecards often carry additive offsets or multiplicative scaling relative to IC/ICP-MS lab assays. Pooling them shifts the spatial mean field $\mathbf{X}\boldsymbol{\beta}$ and inflates the nugget $\tau^2$.
2. **Heteroskedastic Observation Noise**: Lab titration repeatability might have $\mathrm{SE} \approx 0.02\text{ mg/L}$, whereas a citizen colorimetric test strip might effectively have $\mathrm{SE} \approx 0.5\text{ mg/L}$ or censored interval ranges (`1-2 mg/L`). An unweighted kriging solver trusts a noisy reading almost as much as a statutory station, distorting nearby reaches.
3. **Different Analytes or Surrogate Metrics**: Volunteer programmes frequently record ecological or optical proxies (e.g., **Riverfly ARMI**, **Secchi turbidity**, **field conductivity**) rather than regulated chemical determinands (**Nitrate `0117`**, **Orthophosphate `0180`**, **WHPT ASPT**). You cannot concatenate non-identical units into one vector.
4. **Visit-Frequency Imbalance**: A volunteer catchment blitz can upload 30–50 repeat readings on one urban reach while a statutory site is visited once per quarter. Naive arithmetic reach-averaging allows the high-frequency campaign to numerically drown out the accredited measurement.

### How the Pipeline Solves Each Issue

| Fusion Architecture | Script & Flags | When to Use | How Bias, Noise & Imbalance Are Handled |
| :--- | :--- | :--- | :--- |
| **Separate-Target Multivariate Co-Kriging (`Factor-LMC`)** | [`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R)<br>`--targets=0117,wb_no3` `--method=lmc` | **Default recommendation** when citizen data uses a different kit, range scale, proxy metric (e.g. Riverfly vs. `whpt_aspt`), or systematic offset relative to RSN. | Fits separate intercepts, covariate slopes $\boldsymbol{\beta}_p$, transforms, and error variances $\tau_p^2$ per source variable, while borrowing spatial signal through learned **Tail-Up**, **Tail-Down**, and **Co-Located Nugget** cross-correlations. |
| **Colocated Residual Co-Kriging (`RSN`-Anchored)** | [`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R)<br>`--targets=0117,wb_no3` `--method=colocated` | When you want strict statutory primacy: predictions on the primary target (`0117`) are anchored directly to the `RSN` univariate model and adjusted only where correlated auxiliary networks provide verified residual signal. | Fits Target 1 (`RSN`) independently first, then updates predictions via weighted residual conditioning (`lm.wfit` + auxiliary differences), guaranteeing zero likelihood contamination of `RSN` parameters. |
| **Pooled Single-Determinand with Priority & Reliability Weights** | [`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R)<br>`--source=rsn+survey` + `--prefer-source` / `--source-weight` / `--source-se` | When both sources measure the **exact same determinand and lab/calibrated method** (e.g., EA `RSN` + catchment trust lab-accredited Nitrate, or `RSN` + `SSN`), and you want one unified surface while still trusting `RSN` on shared reaches. | Enforces per-reach source priority (`--prefer-source`), aggregates shared reaches **within each `data_source` first** so repeat counts never drown out `RSN`, and down-weights noisy/secondary sources analytically via penalty variance and/or physical measurement SE. |

> [!TIP]
> You can combine both worlds: pass `--source=rsn+waterblitz` into `run_s3n_cokriging.R` with separate targets (`--targets=0117,wb_no3`) **and** still apply `--source-weight=rsn:1,waterblitz:0.4` or `--source-se=rsn:0.02,waterblitz:0.5` so noisy volunteer points are smoothed appropriately while still contributing regional spatial structure!

---

## 2. Core Multi-Source & Network-Preference Options

All three preference/weighting flags work uniformly across both **single-determinand modelling** ([`scripts/run_s3n_determinand.R`](../scripts/run_s3n_determinand.R)) and **multivariate co-kriging** ([`scripts/run_s3n_cokriging.R`](../scripts/run_s3n_cokriging.R)).

### Option Summary

| Flag | Example | Effect |
| :--- | :--- | :--- |
| `--source=<s1>+<s2>` | `--source=rsn+waterblitz` | Merges the preprocessed observations and spatial reach sets of multiple registered or custom sources. |
| `--prefer-source=<s1,s2>` | `--prefer-source=rsn` or `--prefer-source=rsn,ssn,campaign` | **Hard reach-level priority**: on any stream reach (`RivID`) observed by a higher-priority source, samples from lower-priority sources on that exact reach are excluded; lower-priority sources are retained everywhere else. |
| `--source-weight=<s:w,...>` | `--source-weight=rsn:1,waterblitz:0.25` (or `rsn:4,waterblitz:1`) | **Soft continuous down-weighting**: assigns relative reliability weights $w_s \in (0, 1]$ per `data_source` (automatically normalized by $\max_s w_s$ if values exceed `1`). Governs cross-source reach averaging, OLS/GLS trend estimation, spatial covariance inference, kriging predictions, SE, and CV. |
| `--source-se=<s:sd,...>` | `--source-se=rsn:0.02,waterblitz:0.5` | **Explicit physical measurement SE** (in the determinand's **raw measurement units**): attaches `meas_uncertainty` to every sample from source $s$, converts it onto the model transform scale via the Delta Method, activates Approach A measurement-error kriging automatically, and propagates uncertainty into reach predictions and SE maps. |

---

### Mathematical Formulation of `--source-weight` and `--source-se`

#### 1. Source-Grouped Weighted Reach Aggregation

If a stream reach $r$ contains samples from multiple sources $s \in \mathcal{S}_r$, [`aggregate_samples_by_reach_weighted()`](../scripts/s3n_config.R) computes the within-source mean $\bar{y}_{r,s}$ and sample count $n_{r,s}$ **first**, then averages across sources using weights $\tilde{w}_s = w_s / \sum_{k \in \mathcal{S}_r} w_k$:

$$
\bar{y}_r = \sum_{s \in \mathcal{S}_r} \tilde{w}_s \, \bar{y}_{r,s}, \qquad v_{r,\text{phys}} = \sum_{s \in \mathcal{S}_r} \tilde{w}_s^2 \left( \frac{\sum_{i \in r,s} \mathrm{SE}_i^2}{n_{r,s}^2} \right)
$$

Because aggregation happens **by source first**, 50 volunteer dip-strip readings on reach $r$ cannot override 1 statutory `RSN` station by sheer sample count. Meanwhile, the reach's effective observation weight is carried forward as $w_r = \sum_{s \in \mathcal{S}_r} \tilde{w}_s \, w_s$ (which reduces cleanly to $w_s$ on single-source reaches).

#### 2. Observation-Error Variance Decomposition in S3N & Factor-LMC NNGP

For any observed reach $i$ and target $p$, two sources of diagonal observation variance combine without biasing physical nugget estimation:

- **Physical Measurement-Error Variance** ($V_{\text{phys}, ip}$): supplied either by per-sample uncertainty columns (e.g., GAM Stage-1 SE via `--meas-error`) or by `--source-se`, transformed onto the model scale using the exact Delta-Method derivative [`calc_response_meas_variance()`](../scripts/s3n_config.R):

$$
V_{\text{phys}, ip} \approx \frac{\left(g_p'(\mu_{ip}) \, \mathrm{SE}_{ip}\right)^2}{\mathrm{SD}_p^2}
$$

- **Soft Weight-Penalty Variance** ($V_{\text{wt}, ip}$): derived from `--source-weight` $w_{ip} \in (0, 1]$ using the baseline microscale nugget $\tau_{\text{base}, p}^2$:

$$
V_{\text{wt}, ip} = \left(\frac{1}{w_{ip}} - 1\right) \tau_{\text{base}, p}^2
$$

Their sum $V_{\text{eff}, ip} = V_{\text{phys}, ip} + V_{\text{wt}, ip}$ is injected directly into the diagonal of the observation covariance matrix in:

- **Univariate S3N** ([`fit_s3n_domain()`](../scripts/s3n_model_functions.R)): passed as per-site $v_{\text{meas}, i}$ to `C_s3n_predict_meas_error` and `C_s3n_joint_nngp`, while subtracting **only** physical measurement error $\bar{v}_{\text{phys}}$ from the microscale nugget so an artificial weight penalty never collapses $\tau^2_{\text{phys}}$.
- **Multivariate Factor-LMC Co-Kriging** ([`C_s3n_multivariate_nngp_weighted`](../src/s3n_fast_neighbors.cpp)): enters **Mode 0** (joint REML NNGP likelihood & GLS fixed effects), **Mode 1** (full-network co-kriged predictions & prediction SE), **Mode 2** (exact Leave-One-Target-Out Aux-Aided CV), and **5-fold All-Holdout CV**.

Notice the limiting behaviour:
- At $w = 1$ (default), $V_{\text{wt}} = 0$, matching standard S3N / unweighted Factor-LMC NNGP bit-for-bit.
- As $w \to 0^+$, $V_{\text{wt}} \to \infty$: the kriging weights shrink smoothly toward zero on those uncertain points, preventing local spikes while still preserving whatever broader regional spatial covariance structure the network captures.

---

## 3. High-Value Real-World Use Cases

### Use Case 1: National `RSN` + Citizen Science WaterBlitz Colorimetric Kits

**Scenario**: You have national `RSN` Nitrate (`0117`) and Orthophosphate (`0180`) laboratory assays (`~1,210` reaches) plus a spring citizen-science WaterBlitz survey reporting Nitrate and Phosphate in colorimetric brackets (`0-0.2`, `0.2-0.5`, `1-2`, `5-10 mg/L`) across `2,500+` tributary sites.

**Recommended Strategy**:
1. Format the WaterBlitz CSV using `format_custom_data.py` with `--range-mode=tobit` (default), which builds an interval-censored Tobit representation and derives both point estimates and per-interval uncertainty without collapsing intervals heuristically.
2. Assign distinct determinand codes to the citizen kit outputs (e.g., `wb_no3` and `wb_po4`) rather than overwriting `0117` and `0180`, so `RSN` lab measurements and citizen strips keep separate means, slopes, and nuggets.
3. Run 4-target **Factor-LMC co-kriging** (`--targets=0117,0180,wb_no3,wb_po4 --source=rsn+waterblitz --source-weight=rsn:1,waterblitz:0.35 --method=lmc`).
4. **Why it works**:
   - Where `wb_no3` strongly correlates with `0117` in Tail-Up / Tail-Down / Co-located Nugget cross-correlation matrices (`fig_cokrig_cross_correlation.png`), the solver projects calibrated `0117`-equivalent concentrations onto thousands of unmonitored headwaters while automatically learning and discounting the citizen kit's scale offset and extra nugget noise.
   - Check `cokriging_evaluation_report.md`: if **Aux-Aided CV RMSE** and **All-Holdout CV $R^2$** beat Univariate `RSN` S3N, the citizen dataset genuinely improved regulatory prediction accuracy.

---

### Use Case 2: Ecological Statutory `WHPT ASPT` + Volunteer `Riverfly ARMI`

**Scenario**: Statutory Environment Agency macroinvertebrate sampling (`--dataset=inverts`, `whpt_aspt`) covers key surveillance sites, while angler groups run monthly **Riverfly ARMI** counts on smaller tributaries and urban brooks.

**Recommended Strategy**:
- Register Riverfly scores as a custom target (`riverfly_armi`) and co-krige jointly with statutory water quality and invertebrate indices:
  ```bash
  Rscript scripts/run_s3n_cokriging.R \
    --targets=inverts:whpt_aspt,wq:0117,wq:0180,wq:riverfly_armi \
    --source=rsn+riverfly \
    --transforms=none,log,log,none \
    --source-weight=rsn:1,riverfly:0.5 \
    --factors=3 \
    --reverse-color
  ```
- Co-kriging handles the non-identical units (`WHPT ASPT` vs `ARMI` integer scores) natively and shows in `cokriging_cross_covariance.csv` whether Riverfly tracks ecological health (`Corr_Nugget > 0`, `Corr_TU > 0`) and sewage/nutrient pressure (`Corr < 0` with `0180`).

---

### Use Case 3: Accredited `RSN` + Regional Water Company / Catchment Trust Grab Samples (Same Method)

**Scenario**: A Catchment Partnership or Water Company runs accredited laboratory Nitrate (`0117`) analyses matching `RSN` laboratory standards, adding 600 extra monitoring locations in one river basin, with occasional repeat visits on `RSN` reaches.

**Recommended Strategy**:
- Because analytical protocols match, pool into a single-determinand surface while guaranteeing statutory `RSN` primacy on shared reaches and discounting slightly higher field/transport variance if applicable:
  ```bash
  Rscript scripts/run_s3n_determinand.R 0117 \
    --source=rsn+catchment_lab \
    --prefer-source=rsn \
    --source-weight=rsn:1,catchment_lab:0.75 \
    --source-se=rsn:0.02,catchment_lab:0.08
  ```

---

### Use Case 4: High-Frequency Sondes / In-Situ Telemetry + Spot `RSN` Monitoring

**Scenario**: Continuous UV Nitrate or turbidity/DO sensors record thousands of hourly timestamps at 40 strategic gauges, alongside grab-sampled `RSN` sites nationwide.

**Recommended Strategy**:
- Even when pooled via `--source=rsn+sondes --source-weight=rsn:1,sondes:0.6`, the pipeline's source-grouped reach aggregation computes the sonde mean per `RivID` first before cross-source fusion, preventing thousands of sonde timestamps at one gauge from overpowering neighboring grab samples.
- Alternatively, pass sensor drift / fouling uncertainty directly via `--source-se=rsn:0.02,sondes:0.15` (or export per-site temporal variance in `meas_uncertainty` and pass `--meas-error`), ensuring sensor calibration uncertainty is reflected in the final uncertainty map (`fig1` Panel B).

---

## 4. End-to-End Worked Workflow (`RSN` + Citizen Science WaterBlitz)

This complete walk-through assumes you have raw EA `RSN` data under `data/` and a citizen-science CSV `data/raw_waterblitz_2024.csv` containing columns `Latitude`, `Longitude`, `SurveyDate`, `Nitrate_ppm`, and `Phosphate_ppm` (including interval values like `1 - 2`, `<0.2`, or `>10`).

### Step 1: Standardise Citizen Science Data (with Tobit Range & Site Grouping)

First convert the WaterBlitz CSVs into the S3N custom schema (`data/custom/waterblitz/sites.csv`, `samples.csv`) and auto-generate `configs/sources/waterblitz.json`:

```bash
# 1a. Format Nitrate (code: wb_no3)
python3 scripts/format_custom_data.py data/raw_waterblitz_2024.csv \
  --source-name=waterblitz \
  --x-col="Longitude" \
  --y-col="Latitude" \
  --date-col="SurveyDate" \
  --val-col="Nitrate_ppm" \
  --det-code="wb_no3" \
  --det-name="WaterBlitz Nitrate (strip)" \
  --unit="mg/l" \
  --crs=4326 \
  --group-sites-by=coords \
  --site-id-prefix="WB_" \
  --range-mode=tobit \
  --generate-config
```

*(If your CSV is in **long format** with a determinand column and both Nitrate and Phosphate rows, pass `--det-col="Determinand" --val-col="Result" --unit-col="Unit"` once instead. Alternatively, for wide multi-column CSVs like `--val-col="Nitrate_ppm,Phosphate_ppm"`, review [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md) to verify how determinand codes are assigned and rename citizen codes to `wb_no3` / `wb_po4` if separating from `0117` / `0180`.)*

Verify registered sources:

```bash
python3 scripts/01_preprocess_orn.py --list-sources
```

### Step 2: Snap Both Networks to the Stream Network (`ORN v2`)

Run the spatial/topological preprocessor once for the combined composite source `rsn+waterblitz`. This assigns every `RSN` and `WaterBlitz` station to the exact same stream network (`RivID`), computes spatial stream coordinates, and writes `data/processed/orn_rsn+waterblitz/`:

```bash
python3 scripts/01_preprocess_orn.py \
  --network=orn \
  --source=rsn+waterblitz \
  --max-snap-dist=350
```

> [!NOTE]
> Statutory `RSN` coordinates are typically accurate to 10–50 m (`--max-snap-dist=250` default). Smartphone GPS locations logged by volunteers on riverbanks or bridges can be slightly offset; increasing `--max-snap-dist=350` (or `500`) recovers legitimate volunteer readings without snapping sites onto unrelated distant catchments.

### Step 3: Verify Station Snapping Quality

Always inspect snapping distances before fitting spatial models to confirm volunteer sites snapped onto the intended rivers rather than parallel ditches:

```bash
Rscript scripts/06_station_snapping_diagnostics.R \
  --network=orn \
  --source=rsn+waterblitz \
  --catchment=Thames
```

Open `output/station_snapping/orn_rsn+waterblitz/station_snapping_report.md` and `fig_snap_catchment_*.png` to inspect `snap_dist_m` quantiles and zoomed reach insets.

### Step 4: Run Baseline Univariate `RSN` vs. Co-Kriged Multi-Source Model

First inspect seasonal coverage, then run joint 4-target **Factor-LMC Co-Kriging** combining statutory `RSN` (`0117`, `0180`) with citizen `WaterBlitz` (`wb_no3`, `wb_po4`), giving higher reliability weight to `RSN` and realistic physical measurement error:

```bash
Rscript scripts/run_s3n_cokriging.R \
  --targets=0117,0180,wb_no3,wb_po4 \
  --source=rsn+waterblitz \
  --transforms=log,log,log,log \
  --method=lmc \
  --factors=3 \
  --source-weight=rsn:1,waterblitz:0.35 \
  --source-se=rsn:0.05,waterblitz:0.40 \
  --cv-mode=spatial \
  --threshold=uk \
  --omp=4
```

If you want **strict zero-compromise `RSN` primacy** (locking the `0117`/`0180` spatial model strictly to `RSN` and only transferring verified local auxiliary differences), run `--method=colocated`:

```bash
Rscript scripts/run_s3n_cokriging.R \
  --targets=0117,wb_no3 \
  --source=rsn+waterblitz \
  --method=colocated \
  --source-weight=rsn:1,waterblitz:0.35 \
  --source-se=rsn:0.05,waterblitz:0.40 \
  --omp=4
```

### Step 5: Alternative Single-Determinand Weighted Pooling (Matched Methods Only)

If your extra source uses the **same laboratory method** (same code `0117`) and you prefer a fast single-determinand run with hard reach priority on shared `RSN` reaches plus continuous weighting elsewhere:

```bash
Rscript scripts/run_s3n_determinand.R 0117 \
  --source=rsn+partner_lab \
  --prefer-source=rsn \
  --source-weight=rsn:1,partner_lab:0.6 \
  --source-se=rsn:0.02,partner_lab:0.10 \
  --threshold=uk \
  --omp=4
```

### Step 6: Evaluate Added Value & Inspect Diagnostics

Inspect the generated outputs in `output/results/cokrig_0117_0180_wb_no3_wb_po4_orn_rsn+waterblitz/`:

1. **`cokriging_evaluation_report.md`**: Check whether the pipeline verdict rates multi-source co-kriging as `Strong Benefit`, `Moderate Benefit`, or `Marginal/Parity` on out-of-sample `RSN` prediction.
2. **`cv_univariate_vs_cokriging.csv`**: Compare **Univariate S3N** vs. **Co-Kriging (All-Holdout)** vs. **Co-Kriging (Aux-Aided)** RMSE and $R^2$ for `0117` and `0180`.
   - Look for positive **Aux-Aided $R^2$ gain**: this proves that where a reach has WaterBlitz data (`wb_no3`) but lacks `RSN`, knowing the citizen reading reduces prediction error on statutory Nitrate.
   - Look at **All-Holdout $R^2$**: even on completely unmonitored reaches, a gain confirms WaterBlitz sharpened the underlying spatial stream fields (`Tail-Up` / `Tail-Down`).
3. **`cokriging_cross_covariance.csv` & `fig_cokrig_cross_correlation.png`**: Inspect `Corr_TU`, `Corr_TD`, and `Corr_Nugget` between `0117` and `wb_no3`. If cross-correlations are near zero (e.g., an uncalibrated kit dominated by noise or wrong season), co-kriging automatically shrinks cross-variable influence toward zero rather than corrupting `RSN`.
4. **`fig_cokrig_england_0117.png`**: Compare Panel B (Co-Kriging SE) against a standalone `RSN`-only run (`run_s3n_determinand.R 0117 --source=rsn`) to confirm uncertainty shrinks meaningfully in citizen-sampled sub-catchments without introducing salt-and-pepper artefacts.

---

## 5. Practical Caveats & Checklist Before Trusting a Fused Map

- **Include Overlapping Calibration Sites**: For `Factor-LMC` co-kriging to estimate `Corr_Nugget` accurately between `0117` and `wb_no3`, at least some citizen-science readings (typically `>= 20–50` reaches) should fall on or immediately adjacent to `RSN` reaches during comparable hydrological conditions. Having volunteers sample right next to existing EA gauges during a blitz provides high-value empirical calibration!
- **Match Seasons or Flow Regimes**: A single **summer** weekend WaterBlitz has very different dilution and biological uptake than **winter** nitrate peaks. Filter both networks consistently with `--season=spring` / `--season=summer` or `--year=2024` (or check temporal balance in `05_descriptive_stats.R`) so seasonal contrast does not mimic spatial pollution hotspots.
- **Avoid Arbitrary Weights Without Inspection**: Check both `--source-weight` and `--source-se` against out-of-sample CV (`cv_univariate_vs_cokriging.csv` or `cv_diagnostics` in univariate runs). Start with conservative reliability weights (e.g. `0.25–0.5` for uncalibrated field kits, `0.6–0.8` for accredited partner labs) and verify that `RSN` holdout RMSE decreases or holds steady.
- **Keep Distinct Determinand Codes When Scales Differ**: Never merge citizen colorimetric strips or surrogate indices into `0117` with `--source-weight` alone if there is systematic multiplicative or additive bias—keep them as `wb_no3` and let `run_s3n_cokriging.R` absorb offset and scale differences through target-specific $\boldsymbol{\beta}_p$ and $\tau_p^2$.

---

## Related Guides

- [`docs/multivariate_cokriging_guide.md`](./multivariate_cokriging_guide.md) — Factor-LMC NNGP theory, `--method=lmc|colocated`, Dual CV methodology, and figure interpretation.
- [`docs/custom_data_sources_guide.md`](./custom_data_sources_guide.md) — Formatting custom CSVs (`format_custom_data.py`), Tobit interval handling (`--range-mode`), and `configs/sources/*.json` schemas.
- [`docs/two_stage_gam_measurement_error_kriging_guide.md`](./two_stage_gam_measurement_error_kriging_guide.md) — Mathematical details on `meas_uncertainty` and Delta-Method Approach A variance propagation.
- [`docs/station_snapping_guide.md`](./station_snapping_guide.md) — Auditing station snapping distances and reach assignments with `06_station_snapping_diagnostics.R`.
