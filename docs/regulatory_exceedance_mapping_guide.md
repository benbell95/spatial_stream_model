# Regulatory Exceedance Probability & WFD Risk Mapping Guide

Instead of classifying river reaches solely from a point estimate $\hat{y}(s_0)$, S3N universal stream kriging produces a predictive mean $\hat{z}(s_0)$ and kriging standard error $\text{SE}(s_0)$ on the transformed model scale $z = g(y)$. Using the Gaussian predictive distribution $\mathcal{N}(\hat{z}(s_0), \text{SE}(s_0)^2)$ and a raw-scale regulatory threshold $T$, [`run_s3n_determinand.R`](../scripts/run_s3n_determinand.R) and [`02_fit_s3n_models.R`](../scripts/02_fit_s3n_models.R) compute exact reach-level exceedance probabilities and Water Framework Directive (WFD) management risk tiers.

---

## 1. CLI Flags

| Flag | Default | Description |
| :--- | :--- | :--- |
| `--threshold=<val\|uk>` | *(Off)* | Raw-scale regulatory threshold $T$ in the determinand's physical units, **or `uk` (aliases `wfd`, `default`, `ea`)** to automatically load the built-in UK regulatory / WFD Good-Moderate boundary for supported determinands. |
| `--threshold-dir=<dir>` | `above` *(or preset default)* | Direction of regulatory failure:<br>• `above` *(aliases `gt`, `high`, `exceed`)*: $P(Y(s_0) > T)$ — use for pollutants where high values breach standards.<br>• `below` *(aliases `lt`, `low`, `deficit`)*: $P(Y(s_0) < T)$ — use for dissolved oxygen or ecological indices where low values indicate failure. |
| `--list-thresholds` | — | Prints all built-in UK / WFD regulatory threshold presets (`S3N_UK_THRESHOLDS`). |
| `--list-validation` | — | Prints available validation, threshold, and network design flags. |

---

## 2. Built-In UK & WFD Regulatory Thresholds (`--threshold=uk`)

Passing `--threshold=uk` automatically sets both the threshold value and direction (`above` or `below`) from the following UK statutory and indicative Water Framework Directive (WFD) standards:

| Code | Determinand | Default Threshold | Direction | Regulatory Basis & Notes |
| :--- | :--- | :--- | :--- | :--- |
| `0117` | Nitrate as N | `11.3 mg/L` | `above` | UK Drinking Water Inspectorate (DWI) / Nitrates Directive ($11.3\text{ mg/L as N} \equiv 50\text{ mg/L as NO}_3$) |
| `0116` | Nitrogen, Total Oxidised as N | `11.3 mg/L` | `above` | UK Nitrates Directive proxy ($11.3\text{ mg/L as N}$) |
| `0180` | Orthophosphate, reactive as P | `0.069 mg/L` | `above` | UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate SRP boundary) |
| `0111` | Ammoniacal Nitrogen as N | `0.60 mg/L` | `above` | UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate 90-percentile $\text{NH}_4\text{-N}$) |
| `0119` | Ammonia un-ionised as N | `0.021 mg/L` | `above` | UK WFD / Freshwater Fish EQS (annual mean un-ionised $\text{NH}_3\text{-N}$) |
| `0085` | BOD : 5 Day ATU | `5.0 mg/L` | `above` | UK WFD 2015 Directions (indicative lowland high-alkalinity Good/Moderate 90-percentile BOD) |
| `9924` | Oxygen, Dissolved as O2 | `6.0 mg/L` | `below` | UK WFD 2015 Directions (indicative salmonid Good oxygenation threshold) |
| `0061` | pH : In Situ | `6.0 pH` | `below` | UK WFD 2015 Directions (minimum Good pH boundary; use `--threshold=9.0 --threshold-dir=above` for upper bound) |
| `0076` | Temperature of Water | `25.0 °C` | `above` | UK WFD 2015 Directions (cyprinid Good/Moderate 98-percentile temperature) |
| `0052` | Lead, Dissolved | `1.2 µg/L` | `above` | UK WFD Environmental Quality Standard (bioavailable annual average) |
| `0106` | Cadmium, Dissolved | `0.08 µg/L` | `above` | UK WFD Environmental Quality Standard (annual average, low-hardness baseline) |
| `whpt_aspt` | WHPT ASPT (`--dataset=inverts`) | `5.0 ASPT` | `below` | UK WFD / UKTAG (indicative Good/Moderate macroinvertebrate ASPT boundary) |

> **Important Note on River Typology**: Under *The Water Framework Directive (Standards and Classification) Directions (England and Wales) 2015*, statutory WFD boundaries for Reactive Phosphorus, Total Ammonia, and Dissolved Oxygen vary continuously with site altitude and alkalinity (e.g., upland low-alkalinity streams have stricter phosphorus limits than lowland chalk rivers). The `--threshold=uk` presets provide practical indicative lowland/high-alkalinity defaults; you can override any threshold for a specific catchment using `--threshold=<val>`.

### Official UK Regulatory References
- **WFD (Standards and Classification) Directions (England and Wales) 2015** (`legislation.gov.uk`): [https://www.legislation.gov.uk/](https://www.legislation.gov.uk/) (search *"Water Framework Directive (Standards and Classification) Directions (England and Wales) 2015"*)
- **UK Technical Advisory Group (UKTAG) — River & Lake Environmental Standards**: [https://www.wfduk.org/](https://www.wfduk.org/)
- **Drinking Water Inspectorate (DWI) — Nitrate & Drinking Water Standards**: [https://www.dwi.gov.uk/](https://www.dwi.gov.uk/)

---

## 3. Mathematical Formulation

Given the selected response transformation $g(\cdot)$ (e.g. $g(y) = \ln(y)$ for `--transform=log`, $g(y) = y$ for `--transform=none`):

1. Transform the regulatory threshold onto the model scale: $z_T = g(T)$.
2. Compute the standardized normal deviate using the combined S3N prediction $\hat{z}(s_0)$ and kriging standard error $\text{SE}(s_0)$:

   $$
   z_0 = \frac{z_T - \hat{z}(s_0)}{\max(\text{SE}(s_0), 10^{-8})}
   $$

3. Evaluate the exact Gaussian exceedance probability:
   - **Upper threshold (`--threshold-dir=above`)**:

     $$
     P\big(Y(s_0) > T\big) = 1 - \Phi(z_0)
     $$

   - **Lower threshold (`--threshold-dir=below`)**:

     $$
     P\big(Y(s_0) < T\big) = \Phi(z_0)
     $$

Because both $\hat{z}(s_0)$ and $\text{SE}(s_0)$ are evaluated on the transformed scale where S3N normality holds, this approach naturally handles skewed water-quality concentrations (such as log-normal Nitrate or Orthophosphate) without back-transformation bias.

---

## 4. Five-Tier Management Classification

Every reach is assigned a `compliance_tier` based on its exceedance probability $P$:

| Exceedance Probability $P$ | `compliance_tier` | Management Interpretation |
| :--- | :--- | :--- |
| $P < 0.05$ | `Confidently Compliant` | High statistical confidence ($>95\%$) that the reach meets the standard |
| $0.05 \le P < 0.20$ | `Likely Compliant` | Probably compliant, low regulatory risk |
| $0.20 \le P \le 0.80$ | `Uncertain / Priority Monitoring` | Near the threshold or high kriging $\text{SE}(s_0)$; prime candidate for targeted field sampling |
| $0.80 < P \le 0.95$ | `Likely Exceeding` | Probable regulatory failure |
| $P > 0.95$ | `Confidently Exceeding` | High statistical confidence ($>95\%$) that the reach breaches the standard |

---

## 5. Output Files, Plots & Columns

When `--threshold=<val>` is supplied:

1. **`predictions_<domain>_<det>.csv`** gains two columns:
   - `exceed_prob`: Reach-level exceedance probability in `[0, 1]` (rounded to 5 decimal places).
   - `compliance_tier`: One of the 5 management tier labels above.
2. **`exceedance_summary_<det>.csv`** is written in the output directory, summarising across all fitted domains (`All England`, `Humber`, `Thames`, `Severn`, and any `--basin=` domains):
   - `Domain`, `Threshold_Raw`, `Threshold_Sc`, `Direction`, `Tier`, `Reach_Count`, `Reach_Pct`.
3. **Dedicated Management Classification & Exceedance Probability Plots** (generated automatically by `run_s3n_determinand.R` and `03_generate_plots.R`):
   - `fig_exceedance_england_<det>.png`: Two-panel All England map showing **Panel A: Assigned Management Classification** (5-tier green-to-red map with reach counts and percentages in the legend) and **Panel B: Continuous Exceedance Probability** ($0\%$–$100\%$).
   - `fig_exceedance_basins_<det>.png`: Three-panel regional map for **Humber**, **Thames**, and **Severn** showing reaches coloured by their assigned management classification and monitoring site locations.
   - `fig_exceedance_basin_<slug>_<det>.png`: Two-panel management classification + exceedance probability map for each custom `--basin=<name>` domain.

---

## 6. Usage Examples

```bash
# Use built-in UK DWI / Nitrates Directive standard (automatically 11.3 mg/L N, above)
Rscript scripts/run_s3n_determinand.R 0117 --threshold=uk

# Use built-in UK WFD Dissolved Oxygen threshold (automatically 6.0 mg/L, below)
Rscript scripts/run_s3n_determinand.R 9924 --transform=none --threshold=uk

# List all built-in UK / WFD threshold presets
Rscript scripts/run_s3n_determinand.R --list-thresholds

# Override with a custom site-specific threshold (e.g. upland SRP 0.036 mg/L)
Rscript scripts/run_s3n_determinand.R 0180 --threshold=0.036 --threshold-dir=above
```
