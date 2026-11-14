// Curated popular determinands and ecological/habitat metrics from docs/usage_guide.md (Sections 2.3 & 2.4)
window.WORKBENCH_TARGETS = [
  {
    "id": "wq",
    "label": "Water Quality — Common Determinands (RSN + SSN)",
    "dataset": "wq",
    "items": [
      { "value": "0117", "label": "0117 — Nitrate as N (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0180", "label": "0180 — Orthophosphate, reactive as P (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0111", "label": "0111 — Ammoniacal Nitrogen as N (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0116", "label": "0116 — Nitrogen, Total Oxidised as N (TON, mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0118", "label": "0118 — Nitrite as N (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0301", "label": "0301 — Carbon, Organic, Dissolved as C (DOC, mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0077", "label": "0077 — Conductivity at 25 C (uS/cm)", "ds": "wq", "rsnOnly": false },
      { "value": "0162", "label": "0162 — Alkalinity to pH 4.5 as CaCO3 (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0253", "label": "0253 — Magnesium, Dissolved (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "0241", "label": "0241 — Calcium, Dissolved (mg/l)", "ds": "wq", "rsnOnly": false },
      { "value": "6396", "label": "6396 — Turbidity : In Situ (FTU)", "ds": "wq", "rsnOnly": false },
      { "value": "0061", "label": "0061 — pH : In Situ (pH units)", "ds": "wq", "rsnOnly": false },
      { "value": "0076", "label": "0076 — Temperature of Water (CEL)", "ds": "wq", "rsnOnly": false },
      { "value": "0085", "label": "0085 — Oxygen, Dissolved, % Saturation (%)", "ds": "wq", "rsnOnly": false },
      { "value": "9924", "label": "9924 — Oxygen, Dissolved as O2 : In Situ (mg/l)", "ds": "wq", "rsnOnly": false }
    ]
  },
  {
    "id": "inverts",
    "label": "Biology — Macroinvertebrates (--dataset=inverts)",
    "dataset": "inverts",
    "items": [
      { "value": "whpt_aspt", "label": "whpt_aspt — WHPT Average Score Per Taxon (ASPT)", "ds": "inverts", "rsnOnly": false },
      { "value": "whpt_ntaxa", "label": "whpt_ntaxa — WHPT Taxon Richness (N-Taxa)", "ds": "inverts", "rsnOnly": false },
      { "value": "life_family_index", "label": "life_family_index — LIFE Flow Velocity Stress Index (LIFE)", "ds": "inverts", "rsnOnly": false },
      { "value": "psi_family_score", "label": "psi_family_score — PSI Fine Sediment Sensitivity Score (%)", "ds": "inverts", "rsnOnly": false },
      { "value": "epsi_family_score", "label": "epsi_family_score — Empirically Weighted PSI Score (%)", "ds": "inverts", "rsnOnly": false },
      { "value": "cci", "label": "cci — Community Conservation Index (CCI)", "ds": "inverts", "rsnOnly": false },
      { "value": "dehli", "label": "dehli — DEHLI Drought Stress Index (DEHLI)", "ds": "inverts", "rsnOnly": false },
      { "value": "armi_score", "label": "armi_score — Anglers' Riverfly Monitoring Score (score)", "ds": "inverts", "rsnOnly": false }
    ]
  },
  {
    "id": "diatoms",
    "label": "Biology — Diatoms (--dataset=diatoms)",
    "dataset": "diatoms",
    "items": [
      { "value": "percent_organic_tolerant", "label": "percent_organic_tolerant — Organic-Tolerant Diatom Abundance (%)", "ds": "diatoms", "rsnOnly": false },
      { "value": "percent_motile", "label": "percent_motile — Motile Diatom Guild Abundance (%)", "ds": "diatoms", "rsnOnly": false }
    ]
  },
  {
    "id": "macrophytes",
    "label": "Biology — Macrophytes (--dataset=macrophytes)",
    "dataset": "macrophytes",
    "items": [
      { "value": "rmni", "label": "rmni — River Macrophyte Nutrient Index (RMNI)", "ds": "macrophytes", "rsnOnly": false },
      { "value": "rn_a_taxa", "label": "rn_a_taxa — Macrophyte Aquatic Taxon Richness (taxa)", "ds": "macrophytes", "rsnOnly": false }
    ]
  },
  {
    "id": "rhs",
    "label": "Habitat — River Habitat Surveys (--dataset=rhs, RSN only)",
    "dataset": "rhs",
    "items": [
      { "value": "hqa_score", "label": "hqa_score — Habitat Quality Assessment (HQA, RSN only)", "ds": "rhs", "rsnOnly": true },
      { "value": "hms_score", "label": "hms_score — Habitat Modification Score (HMS, RSN only)", "ds": "rhs", "rsnOnly": true }
    ]
  },
  {
    "id": "multi",
    "label": "Multi-Target Suites (Co-Kriging)",
    "dataset": "multi",
    "items": [
      { "value": "0117,0180", "label": "0117,0180 — Nitrate + Orthophosphate", "ds": "multi", "rsnOnly": false },
      { "value": "0117,0180,0111", "label": "0117,0180,0111 — N + P + Ammonia Nutrient Triad", "ds": "multi", "rsnOnly": false },
      { "value": "0117,0180,9924", "label": "0117,0180,9924 — Nutrients + Dissolved Oxygen", "ds": "multi", "rsnOnly": false },
      { "value": "0117,0180,whpt_aspt", "label": "0117,0180,whpt_aspt — Chemistry + WHPT Ecology", "ds": "multi", "rsnOnly": false },
      { "value": "whpt_aspt,psi_family_score,life_family_index", "label": "whpt_aspt,psi_family_score,life_family_index — Macroinvertebrate Triad", "ds": "multi", "rsnOnly": false },
      { "value": "rmni,0180,0117", "label": "rmni,0180,0117 — Macrophyte RMNI + Nutrients", "ds": "multi", "rsnOnly": false }
    ]
  }
];
