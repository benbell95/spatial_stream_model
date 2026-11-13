/* ==========================================================================
   Spatial Stream Network Modelling Dashboard - Zero-Dependency Vanilla JS
   ========================================================================== */
(function () {
  "use strict";

  const DATA = window.MODELLING_DASHBOARD_DATA || { figures: {}, tables: {}, reports: {}, meta: {} };

  const FIGURE_META = {
    no3_england: {
      title: "England Nitrate (0117) — Joint TU+TD+Euc Surface & Uncertainty",
      note: "90,876 reaches, 1,210 RSN stations. Notice sharp east-west arable Nitrate contrast and reduced kriging SE along monitored mainstems."
    },
    no3_gam_england: {
      title: "England Nitrate (0117) — Non-Linear REML GAM + S3N Trend",
      note: "Penalized cubic regression splines capture non-linear elevation/drainage thresholds before NNGP spatial residual kriging."
    },
    no3_humber: {
      title: "Humber Basin (CatchID=80) — High-Resolution Stream Nitrate",
      note: "222 stations across 14,048 reaches. Combined CV R2 = 0.6995 (vs 0.6546 OLS)."
    },
    no3_thames_severn: {
      title: "Thames (3055) & Severn (2786) Basins — Nitrate Predictions & SE",
      note: "Regional sub-models allow aquifer and slope fixed effects to adapt locally per drainage basin."
    },
    no3_diag: {
      title: "4-Panel Nitrate Diagnostics, Cross-Validation & Parameter CIs (SVG)",
      note: "Compares OLS, Tail-Down, Tail-Up, and Combined NNGP CV accuracy alongside bootstrap covariate & variance intervals."
    },
    no3_gam_diag: {
      title: "4-Panel GAM + S3N Diagnostics & Cross-Validation (SVG)",
      note: "Vector SVG diagnostic summary for non-linear --trend=gam Nitrate model (CV R2 = 0.6743)."
    },
    whpt_england: {
      title: "England Macroinvertebrate Health (WHPT ASPT, --reverse-color)",
      note: "Inverted ecological palette: blue = high ecological quality (upland clean streams), red = impacted lowland/urban corridors."
    },
    whpt_thames_severn: {
      title: "Thames & Severn Ecological WHPT ASPT Maps",
      note: "Biological communities integrate chronic organic/nutrient stress; S3N captures strong stream-network spatial continuity."
    },
    whpt_diag: {
      title: "WHPT ASPT Cross-Validation & Diagnostics (SVG)",
      note: "Combined S3N boosts All-England WHPT ASPT CV R2 from 0.4932 (OLS) to 0.6402."
    },
    no3_z_england: {
      title: "England Nitrate Robust Z-Score Map (--norm=robust_zscore)",
      note: "Expresses deviations as (x - median) / (1.4826 * MAD) on the model scale to pinpoint regional pollution anomalies."
    },
    no3_z_humber: {
      title: "Humber Basin Nitrate Robust Z-Score Map (--norm=robust_zscore)",
      note: "High-resolution Humber sub-basin view highlighting localized nitrate anomalies relative to national median/MAD."
    },
    no3_z_thames_severn: {
      title: "Thames & Severn Nitrate Robust Z-Score Maps (--norm=robust_zscore)",
      note: "Contrasts anomalous Chalk/arable Thames reaches (+Z) against mixed Severn gradients on a dimensionless robust scale."
    },
    no3_exceed_england: {
      title: "Nitrate 5-Tier Exceedance Risk (Threshold: 11.3 mg/L as N)",
      note: "Probabilistic Gaussian risk tiers account for both mean prediction and kriging SE; 61.3% of reaches are Confidently Compliant."
    },
    po4_exceed_england: {
      title: "Orthophosphate 5-Tier Exceedance Risk (Threshold: 0.069 mg/L as P)",
      note: "47.5% of England reaches fall in 'Uncertain / Priority Monitoring' (0.20 <= P <= 0.80 around 0.069 mg/L)."
    },
    no3_exceed_basins: {
      title: "Humber, Thames & Severn Nitrate Exceedance Risk Maps",
      note: "Basin-level threshold risk view identifying localized high-risk tributaries where exceedance confidence exceeds 95%."
    },
    opt_design_england: {
      title: "All England Optimal Monitoring Design & Gap Map (--design-sites=50)",
      note: "Top 50 greedy C++ active-learning placements reduce evaluated England variance by 80.74% (Mean Kriging SE 0.7961 -> 0.3228 log)."
    },
    opt_design_humber: {
      title: "Humber Basin Optimal Monitoring Design Map (CatchID=80, K=50)",
      note: "Ranks top 50 Humber gap-bridging tributaries, achieving 32.31% cumulative variance reduction (Mean SE 0.7177 -> 0.5714)."
    },
    opt_design_thames: {
      title: "Thames Basin Optimal Monitoring Design Map (CatchID=3055, K=50)",
      note: "Targets under-monitored Thames headwater and tributary gaps for 43.55% cumulative variance reduction (Mean SE -> 0.5207)."
    },
    opt_design_severn: {
      title: "Severn Basin Optimal Monitoring Design Map (CatchID=2786, K=50)",
      note: "Places 50 priority stations along Severn mid-reach and headwater gaps for 43.01% cumulative variance reduction."
    },
    st_season_grid: {
      title: "All England Nitrate by Meteorological Season (fig5 Grid)",
      note: "Synchronized shared color scale across Winter (DJF), Spring (MAM), Summer (JJA), and Autumn (SON)."
    },
    st_year_grid: {
      title: "All England Nitrate by Calendar Year 2021–2024 (fig5 Grid)",
      note: "Shared-scale annual spatial snapshots enabling direct multi-year visual comparison across England."
    },
    st_year_trends: {
      title: "Reach-Level Annual Nitrate Trend Slopes & Wald Z-Tests (fig6)",
      note: "Highlights reaches with statistically significant multi-year worsening (red) or improving (blue) trajectories."
    },
    st_season_trends: {
      title: "Reach-Level Seasonal Trend Slopes & Significance Tiers (fig6)",
      note: "Dynamic temporal borrowing stabilizes noisy seasons before testing slope significance."
    },
    st_season_trans: {
      title: "Consecutive Season-to-Season Nitrate Transitions (fig7)",
      note: "Visualizes Winter -> Spring -> Summer -> Autumn hydrological flushing and dilution shifts."
    },
    cokrig_xcorr: {
      title: "Factor-LMC Cross-Correlation Matrices & Dual CV Gains (SVG)",
      note: "Decomposes Nitrate-Orthophosphate coupling into Tail-Up (0.362), Tail-Down (0.691), and Co-Located Nugget (-0.734) correlations plus CV gains."
    },
    cokrig_0117: {
      title: "Co-Kriged England Nitrate Surface & Joint SE (--targets=0117,0180)",
      note: "Leverages shared latent spatial factors with Orthophosphate to lower uncertainty across catchment networks."
    },
    cokrig_0180: {
      title: "Co-Kriged England Orthophosphate Surface & Joint SE",
      note: "Simultaneously estimated with Nitrate; Humber & Thames All-Holdout and Aux-Aided accuracy both outperform Univariate S3N."
    },
    rss_england: {
      title: "Multi-Source RSN + SSN Nitrate Surface (--prefer-source=rsn --source-weight=rsn:1,ssn:0.75)",
      note: "Expands observed network to 1,384 reaches across 91,118 stream segments while prioritizing RSN on co-located sites."
    },
    rss_diag: {
      title: "Multi-Source RSN+SSN Diagnostics & Measurement Error CV (SVG)",
      note: "Analytic (1/w - 1)*tau^2 penalty weighting preserves statutory RSN primacy while extracting headwater signal from SSN."
    },
    desc_dist: {
      title: "Nitrate Network Distributions & Sampling Bias Diagnostics",
      note: "Shows why monitored RSN mean (4.09 mg/L) exceeds predicted full-network mean (2.61 mg/L, p = 5.1e-45)."
    },
    desc_grad: {
      title: "Stratified Landscape, Aquifer & Strahler Order Gradients (SVG)",
      note: "Vector SVG breakdown comparing upland vs lowland, Principal Aquifers, Strahler stream orders, and major RBDs."
    },
    comp_cv: {
      title: "Nitrate vs. Orthophosphate Cross-Validation Benchmark (SVG)",
      note: "Side-by-side R2 & RMSE comparison across England, Humber, Thames, and Severn."
    },
    comp_forest: {
      title: "Covariate Effect Forest Plot: Nitrate vs. Orthophosphate (SVG)",
      note: "Two-sample Wald Z-tests highlight how Chalk/Limestone aquifers increase Nitrate (+0.18) while slope dilutes Orthophosphate."
    },
    comp_maps: {
      title: "Side-by-Side National Maps: Nitrate (0117) vs. Orthophosphate (0180)",
      note: "Visualizes contrasting arable groundwater Nitrate hotspots versus urban/point-source Orthophosphate corridors."
    },
    comp_rel: {
      title: "Reach-Level N:P Stoichiometry & 4-Quadrant Co-Occurrence",
      note: "Identifies reaches with coupled nutrient enrichment vs N-dominated or P-dominated catchments."
    },
    snap_overview: {
      title: "National ORN v2 Station-to-Reach Snapping Overview (1,583 Stations)",
      note: "Distribution of snapping chord distances across 141 catchments under the default 250 m inclusion threshold."
    },
    snap_thames: {
      title: "Thames Basin Station Snapping Diagnostics & Metre-Scale Zoom Insets",
      note: "154 Thames stations; 83.6% snap within 10 m (3.13 m median distance) with 4 zoomed station-to-reach vector insets."
    },
    snap_humber: {
      title: "Humber Basin Station Snapping Diagnostics & Zoom Insets",
      note: "294 Humber stations; 88.6% snap within 10 m (2.31 m median distance)."
    },
    wb_diff_0117: {
      title: "Nitrate (0117) Co-Kriging Difference Map: Favoured minus Equal Sources",
      note: "Spatial shift across 95,773 England reaches; co-located statutory reaches rise (+0.6252 median abs shift) as WaterBlitz dilution is removed."
    },
    wb_diff_0180: {
      title: "Orthophosphate (0180) Co-Kriging Difference Map: Favoured minus Equal Sources",
      note: "Highlights threshold tier divergence (30.84% tier match on WaterBlitz-only reaches) when statutory RSN+SSN anchors dominate."
    },
    wb_eq_0117: {
      title: "Nitrate (0117) Equal-Sources Co-Kriging Map (RSN+SSN+WaterBlitz)",
      note: "All 58,765 samples (8,346 reaches) weighted equally; England mean Kriging SE = 0.9306, latent total cross-corr = +0.2989."
    },
    wb_wt_0117: {
      title: "Nitrate (0117) Favoured RSN+SSN Co-Kriging Map (--prefer-source=rsn,ssn)",
      note: "Drops 2,903 citizen samples on 338 co-located reaches and applies 0.25 weight (3*tau^2 variance penalty) elsewhere; mean SE drops -14.53% to 0.7954."
    },
    wb_eq_0180: {
      title: "Orthophosphate (0180) Equal-Sources Co-Kriging Map (RSN+SSN+WaterBlitz)",
      note: "55,356 samples across 7,512 reaches; 48,802 reaches (51.0%) fall in Uncertain Compliance Priority tier (0.20 <= P <= 0.80 at 0.069 mg/L)."
    },
    wb_wt_0180: {
      title: "Orthophosphate (0180) Favoured RSN+SSN Co-Kriging Map (waterblitz:0.25)",
      note: "Drops 2,672 co-located citizen samples; Uncertain Phosphate reaches fall from 51.0% to 41.7% (39,969 reaches) and national SE drops -23.43%."
    },
    wb_snap_overview: {
      title: "Combined RSN + SSN + WaterBlitz ORN v2 Station Snapping Overview",
      note: "38,921 stations plotted across 268 catchments (Thames 18,017; Humber 4,402); 33,013 stations snapped <=350 m onto 8,346 active reaches."
    },
    wb_comp_svg: {
      title: "Equal vs. Favoured RSN+SSN Multi-Stratum Summary Comparison (SVG)",
      note: "Compares median raw predictions, SE inflation/reduction across strata, compliance tier agreement, and latent cross-correlation."
    },
    wb_eq_xcorr: {
      title: "Equal-Sources Latent Cross-Correlation Matrix (0117 vs. 0180, SVG)",
      note: "All-England Corr_TU = +0.5214, Corr_TD = +0.5230, Corr_Nugget = +0.1335, Corr_Total = +0.2989."
    },
    wb_wt_xcorr: {
      title: "Favoured RSN+SSN Latent Cross-Correlation Matrix (0117 vs. 0180, SVG)",
      note: "De-weighting citizen noise lifts England Corr_Nugget (+0.1335 -> +0.7083) and Corr_Total (+0.2989 -> +0.4308)."
    },
    tr_diff_map: {
      title: "Nitrate Transform Divergence Map (Raw vs. Log, Negatives & Sqrt vs. Log)",
      note: "Highlights 6,349 physically impossible negative reaches (min -4.817 mg/L) under Raw (none) vs. clean positivity and convergence under Sqrt / Box-Cox."
    },
    tr_raw_map: {
      title: "England Nitrate — Untransformed Raw (--transform=none)",
      note: "Native mg/L scale (CV R2 = 0.5335) violates physical positivity with 47 negative CV predictions and P05 = -0.288 mg/L."
    },
    tr_sqrt_map: {
      title: "England Nitrate — Fixed Square Root (--transform=sqrt)",
      note: "Model-scale CV R2 = 0.6349, 0 negative reaches, lowest MAE (1.5705 mg/L), and near-Gaussian CV residuals (+0.184)."
    },
    tr_boxcox_map: {
      title: "England Nitrate — Profile-MLE Box-Cox (lambda = 0.50)",
      note: "Automatic MLE converges to lambda = 0.50 (CV R2 = 0.6352, Raw MAE = 1.5704 mg/L, 0 negative reaches)."
    },
    tr_side_by_side: {
      title: "4-Transform National Nitrate Prediction Maps Side-by-Side (0117)",
      note: "Identical 1,210 RSN stations on 90,876 reaches: Raw (none), Log (ln), Sqrt (0.50), and Profile-MLE Box-Cox (0.50)."
    },
    tr_comp_svg: {
      title: "6-Panel Transform Diagnostics Summary (SVG)",
      note: "Compares CV R2 (model vs. raw scale), Raw RMSE/MAE, Skewness & Negatives, River Quantiles, Management Tiers, and Spearman Rank Concordance."
    },
    tr_pred_dens: {
      title: "Back-Transformed Nitrate Density, Quantiles & Positivity Bounds",
      note: "Shows Raw crossing 0 mg/L (6.99% negative reaches), Log median downward shift (-1.0078 mg/L raw bias), and Sqrt/Box-Cox alignment."
    },
    tr_param_forest: {
      title: "Model-Scale Variance Components & Parameter Comparison",
      note: "Tail-Up, Tail-Down, Nugget, and effective range shifts across Raw, Log, Sqrt, and Profile-MLE Box-Cox."
    }
  };

  const FEATURE_CATALOG = [
    { cat: "major", badge: "CS7", title: "Data Transformation Comparison (Raw, Log, Sqrt, Box-Cox)", desc: "Benchmarks 4 Nitrate transformations across physical positivity, Jensen median bias, residual symmetry (0.118 skew), and Profile-MLE Box-Cox (lambda=0.50).", tab: "cs7" },
    { cat: "major", badge: "CS6", title: "Citizen Science WaterBlitz Co-Kriging Fusion", desc: "Ingest 44,982 England WaterBlitz citizen observations via format_custom_data.py and compare Equal vs. Statutory-Favoured (rsn,ssn priority + 0.25 weight) bivariate Co-Kriging.", tab: "cs6" },
    { cat: "major", badge: "Stage 1", title: "Zero-Dependency Python WKB GIS Extractor", desc: "Pure stdlib sqlite3 + binary WKB decoder extracts 90,876 ORN v2 reaches in ~8-15s without GDAL or sf.", tab: "overview" },
    { cat: "major", badge: "Stage 2", title: "C++/OpenMP Euler-Tour Joint NNGP", desc: "Sub-second O(n) Tail-Up + Tail-Down + Euclidean additive kriging via compiled s3n_fast_neighbors.so.", tab: "cs1" },
    { cat: "major", badge: "CS1", title: "Non-Linear GAM, Spline, RF & GBM Trends", desc: "Replace linear X*beta with --trend=gam, spline, rf, or gbm plus orthogonal OOB/OOF spatial residual kriging.", tab: "cs1" },
    { cat: "major", badge: "CS2", title: "5-Tier Threshold Exceedance Risk Mapping", desc: "--threshold=uk applies indicative UK thresholds (11.3 mg/L NO3-N, 0.069 mg/L P) to map 5-tier Gaussian exceedance probability.", tab: "cs2" },
    { cat: "major", badge: "CS2", title: "Ecological Metrics & Inverted Scales", desc: "Supports inverts (whpt_aspt), diatoms, macrophytes, rhs & lcms with --reverse-color palettes.", tab: "cs2" },
    { cat: "major", badge: "CS3", title: "Spatio-Temporal Batching & Dynamic Borrowing", desc: "--by=year|season|month + --st-borrow precision weighting, fig6 trend tests & fig7 period transitions.", tab: "cs3" },
    { cat: "major", badge: "CS4", title: "Multivariate Factor-LMC Joint Co-Kriging", desc: "Models P >= 2 heterotopic targets simultaneously via --method=lmc|colocated with Dual CV benchmarks.", tab: "cs4" },
    { cat: "major", badge: "CS4", title: "Multi-Source Fusion & Network Preference", desc: "Combine statutory RSN + SSN + citizen science using --prefer-source, --source-weight & --source-se.", tab: "cs4" },
    { cat: "major", badge: "CS5", title: "Descriptive Stats, Bias & N:P Stoichiometry", desc: "Scripts 05 & 04 run Mann-Whitney representativeness tests, gradient breakdowns, Fisher z & Wald Z-tests.", tab: "cs5" },
    { cat: "minor", badge: "Kernels", title: "5 Stream Covariance Kernels (--kernel)", desc: "exponential, spherical, mariah, linear, or auto (selects best kernel by Profile MLE likelihood).", tab: "cs1" },
    { cat: "minor", badge: "Confluence", title: "3 Physical AFV Weighting Schemes (--afv)", desc: "Branch weights by upstream length, Shreve link magnitude, or hydro_yield (orographic + aquifer proxy).", tab: "cs1" },
    { cat: "minor", badge: "Covariates", title: "15 Covariates & Mixed +/- Selection (--vars)", desc: "Customize predictors with --vars=+Soil_Idx,-Slope_z; inspect full catalog via --list-vars.", tab: "cs1" },
    { cat: "minor", badge: "Transforms", title: "12 Transforms & 4 Normalisations", desc: "boxcox, yeojohnson, log, sqrt, logit, asin + Taylor/Duan back-transform; robust_zscore, zscore, minmax, percentile.", tab: "cs7" },
    { cat: "minor", badge: "Censoring", title: "ROS, Site-ROS & Spatial Tobit EM (< LOD)", desc: "Handle left-censored '<' detection limits (--censor=ros|ros_site|em) and colorimetric ranges (--range-mode=tobit).", tab: "cs2" },
    { cat: "minor", badge: "Uncertainty", title: "Two-Stage Approach A Meas-Error Kriging", desc: "--meas-error / --gam and --source-se propagate Delta-Method observation SE into prediction SE maps.", tab: "cs2" },
    { cat: "major", badge: "CS8", title: "Optimal Monitoring Network Design & Gap Mapping", desc: "--design-sites=<K> runs OpenMP C++ Schur-complement active learning and renders multi-panel gap maps (80.74% variance drop at K=50).", tab: "cs8" },
    { cat: "minor", badge: "Validation", title: "Hydrological CV Modes (--cv-mode)", desc: "--cv-mode=random|spatial|catchment|downstream|upstream evaluates geographic & downstream extrapolation.", tab: "cs1" },
    { cat: "minor", badge: "Networks", title: "BYOWN Custom GeoPackage & GeoParquet", desc: "Bring your own river shapefile/gpkg via JSON profile (--network) or convert via --convert-to-geoparquet.", tab: "overview" },
    { cat: "minor", badge: "Graphics", title: "Hybrid DPI, Vector SVG & Size Scaling", desc: "--fig-format=png|pdf|svg, --dpi=<N>, and --fig-size=<pct> across standalone & runner scripts.", tab: "workbench" }
  ];

  function esc(val) {
    return String(val == null ? "" : val)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function fmt(val, digits) {
    if (val === null || val === undefined || val === "") return "—";
    if (typeof val === "number") return Number.isInteger(val) ? val.toLocaleString() : val.toFixed(digits != null ? digits : 4);
    return esc(val);
  }

  function switchTab(tabId, updateHash) {
    document.querySelectorAll(".tab-btn").forEach((btn) => {
      btn.classList.toggle("active", btn.dataset.tab === tabId);
    });
    document.querySelectorAll(".tab-panel").forEach((p) => {
      p.classList.toggle("active", p.id === "tab-" + tabId);
    });
    if (updateHash !== false) {
      history.replaceState(null, "", "#" + tabId);
    }
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  function renderFigureInto(slotId, figKey) {
    const slot = document.getElementById(slotId);
    if (!slot) return;
    const fig = (DATA.figures || {})[figKey];
    const meta = FIGURE_META[figKey] || { title: figKey, note: "" };
    if (!fig) {
      slot.innerHTML = `<div class="caveat-box">Figure asset <code>${esc(figKey)}</code> not staged yet. Run <code>python3 scripts/build_showcase.py</code>.</div>`;
      return;
    }
    slot.innerHTML = `
      <div class="figure-frame" data-open-fig="${esc(figKey)}" title="Click to expand high-resolution view">
        <img src="${esc(fig.path)}" alt="${esc(meta.title)}" loading="lazy" />
      </div>
      <div class="figure-caption">
        <span><strong>${esc(meta.title)}</strong> &mdash; ${esc(meta.note)}</span>
        <span class="pill pill-green">${esc(fig.format)} &middot; ${esc(fig.size_kb)} KB</span>
      </div>
    `;
  }

  function initFigureToggles() {
    const defaults = [
      ["fig-slot-cs1-map", "no3_england"],
      ["fig-slot-cs1-diag", "no3_diag"],
      ["fig-slot-cs2-ecol", "whpt_england"],
      ["fig-slot-cs2-risk", "no3_exceed_england"],
      ["fig-slot-cs3-grid", "st_season_grid"],
      ["fig-slot-cs3-trend", "st_year_trends"],
      ["fig-slot-cs4-ck", "cokrig_xcorr"],
      ["fig-slot-cs4-ms", "rss_england"],
      ["fig-slot-cs5-audits", "desc_grad"],
      ["fig-slot-cs5-snap", "snap_overview"],
      ["fig-slot-cs6-maps", "wb_diff_0117"],
      ["fig-slot-cs6-charts", "wb_comp_svg"],
      ["fig-slot-cs7-maps", "tr_diff_map"],
      ["fig-slot-cs7-charts", "tr_comp_svg"],
      ["fig-slot-cs8-design", "opt_design_england"]
    ];
    defaults.forEach(([slotId, figKey]) => renderFigureInto(slotId, figKey));

    document.querySelectorAll(".fig-toggle-bar[data-group]").forEach((bar) => {
      const group = bar.dataset.group;
      const slotId = "fig-slot-" + group;
      bar.querySelectorAll(".fig-toggle-btn").forEach((btn) => {
        btn.addEventListener("click", () => {
          bar.querySelectorAll(".fig-toggle-btn").forEach((b) => b.classList.remove("active"));
          btn.classList.add("active");
          renderFigureInto(slotId, btn.dataset.fig);
        });
      });
    });
  }

  function renderTable(elId, headers, rows, rowClassFn) {
    const el = document.getElementById(elId);
    if (!el) return;
    if (!rows || rows.length === 0) {
      el.innerHTML = `<p style="padding:12px; color:var(--text-muted); margin:0;">No rows found in bundle.</p>`;
      return;
    }
    const ths = headers.map((h) => `<th class="${h.num ? "num" : ""}">${esc(h.label)}</th>`).join("");
    const trs = rows.map((r) => {
      const cls = rowClassFn ? rowClassFn(r) : "";
      const tds = headers.map((h) => {
        const raw = h.fn ? h.fn(r) : r[h.key];
        const txt = typeof raw === "number" ? fmt(raw, h.digits) : (raw == null ? "—" : esc(raw));
        const tdCls = [h.num ? "num" : "", h.clsFn ? h.clsFn(r) : ""].filter(Boolean).join(" ");
        return `<td class="${tdCls}">${txt}</td>`;
      }).join("");
      return `<tr class="${cls}">${tds}</tr>`;
    }).join("");
    el.innerHTML = `<table><thead><tr>${ths}</tr></thead><tbody>${trs}</tbody></table>`;
  }

  function buildInlineBarSvg(title, bars, maxVal, unitSuffix) {
    const h = 36 + bars.length * 32;
    const rows = bars.map((b, i) => {
      const y = 30 + i * 32;
      const w = Math.max(4, Math.min(340, (b.val / maxVal) * 340));
      return `
        <text x="12" y="${y + 15}" fill="var(--text-primary)" font-size="12" font-weight="600">${esc(b.label)}</text>
        <rect x="215" y="${y + 3}" width="340" height="16" rx="4" fill="var(--bg-elevated)" />
        <rect x="215" y="${y + 3}" width="${w.toFixed(1)}" height="16" rx="4" fill="${b.color || "#38bdf8"}" />
        <text x="565" y="${y + 15}" fill="var(--text-primary)" font-size="12" font-weight="700">${b.val.toFixed(4)}${esc(unitSuffix || "")}</text>
      `;
    }).join("");
    return `<svg viewBox="0 0 650 ${h}" width="100%" style="max-width:650px;">
      <text x="12" y="18" fill="var(--text-secondary)" font-size="11" font-weight="700">${esc(title)}</text>
      ${rows}
    </svg>`;
  }

  function populateCaseStudyData() {
    const t = DATA.tables || {};

    // --- CS1: Nitrate + GAM ---
    const eng0117 = (t.cv_0117 || []).filter((r) => r.Domain === "All England RSN Network");
    const engGam = (t.cv_0117_gam || []).filter((r) => r.Domain === "All England RSN Network");
    const combLin = eng0117.find((r) => String(r.Model).includes("Combined")) || {};
    const olsLin = eng0117.find((r) => String(r.Model).includes("OLS")) || {};
    const combGam = engGam.find((r) => String(r.Model).includes("Combined")) || {};

    const kpi1 = document.getElementById("cs1-kpis");
    if (kpi1) {
      kpi1.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Linear S3N CV R²</div><div class="kpi-val pos">${fmt(combLin.R2, 4)}</div><div class="kpi-note">RMSE ${fmt(combLin.RMSE, 4)} (vs ${fmt(olsLin.RMSE, 4)} OLS)</div></div>
        <div class="kpi-card"><div class="kpi-label">GAM + S3N CV R²</div><div class="kpi-val pos">${fmt(combGam.R2, 4)}</div><div class="kpi-note">REML penalized cubic regression splines</div></div>
        <div class="kpi-card"><div class="kpi-label">Spatial Variance Gain</div><div class="kpi-val">+0.2126</div><div class="kpi-note">R² lift from Tail-Up + Tail-Down + Euclidean</div></div>
        <div class="kpi-card"><div class="kpi-label">Observed Sites</div><div class="kpi-val">1,210 reaches</div><div class="kpi-note">Predicted out-of-sample across 90,876 reaches</div></div>
      `;
    }

    const r2Svg = document.getElementById("cs1-r2-svg");
    if (r2Svg && eng0117.length) {
      const bars = eng0117.map((r) => ({
        label: r.Model.replace(" (Flow-Connected + AFV)", "").replace(" (Tail-Up + Tail-Down)", "").replace(" (Stream Distance)", ""),
        val: r.R2 || 0,
        color: r.Model.includes("Combined") ? "#34d399" : (r.Model.includes("OLS") ? "#64748b" : "#38bdf8")
      }));
      if (combGam.R2) bars.push({ label: "GAM + S3N Combined", val: combGam.R2, color: "#a78bfa" });
      r2Svg.innerHTML = buildInlineBarSvg("ALL ENGLAND NITRATE (0117) 5-FOLD CV R2 COMPARISON", bars, 0.85, "");
    }

    renderTable("cs1-cv-table",
      [
        { label: "Domain", key: "Domain" },
        { label: "Model", key: "Model" },
        { label: "n_obs", key: "n_obs", num: true },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "MAE", key: "MAE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      (t.cv_0117 || []).filter((r) => r.Model && (r.Model.includes("OLS") || r.Model.includes("Combined"))),
      (r) => String(r.Model).includes("Combined") ? "highlight-row" : ""
    );

    renderTable("cs1-param-table",
      [
        { label: "Parameter (All England Combined)", key: "Parameter" },
        { label: "Estimate", key: "Estimate", num: true, digits: 4 },
        { label: "CI Lower", key: "CI_Lower", num: true, digits: 4 },
        { label: "CI Upper", key: "CI_Upper", num: true, digits: 4 }
      ],
      t.params_0117 || []
    );

    // --- CS2: Exceedance + WHPT ASPT ---
    const ex0117 = (t.exceed_0117 || []).filter((r) => r.Domain === "All England RSN Network");
    const ex0180 = (t.exceed_0180 || []).filter((r) => r.Domain === "All England RSN Network");
    const exSvg = document.getElementById("cs2-exceed-svg");
    if (exSvg && (ex0117.length || ex0180.length)) {
      const tierColors = {
        "Confidently Compliant": "#0284c7",
        "Likely Compliant": "#38bdf8",
        "Uncertain / Priority Monitoring": "#fbbf24",
        "Likely Exceeding": "#f97316",
        "Confidently Exceeding": "#ef4444"
      };
      const mkStack = (rows, y, label) => {
        let x = 165;
        const rects = rows.map((r) => {
          const w = Math.max(0, (r.Reach_Pct || 0) * 4.5);
          const col = tierColors[r.Tier] || "#94a3b8";
          const s = `<rect x="${x.toFixed(1)}" y="${y}" width="${w.toFixed(1)}" height="22" fill="${col}"><title>${esc(r.Tier)}: ${r.Reach_Pct}%</title></rect>`;
          x += w;
          return s;
        }).join("");
        return `<text x="10" y="${y + 15}" fill="var(--text-primary)" font-size="12" font-weight="600">${esc(label)}</text>${rects}`;
      };
      const legend = Object.entries(tierColors).map(([name, c], i) => {
        const lx = 10 + (i % 3) * 210;
        const ly = 96 + Math.floor(i / 3) * 20;
        return `<rect x="${lx}" y="${ly}" width="11" height="11" rx="2" fill="${c}"/><text x="${lx + 16}" y="${ly + 10}" fill="var(--text-secondary)" font-size="11">${esc(name)}</text>`;
      }).join("");
      exSvg.innerHTML = `<svg viewBox="0 0 650 142" width="100%" style="max-width:650px;">
        <text x="10" y="16" fill="var(--text-secondary)" font-size="11" font-weight="700">ENGLAND 5-TIER COMPLIANCE PROBABILITY SHARE (90,876 REACHES)</text>
        ${mkStack(ex0117, 28, "Nitrate (11.3 mg/L)")}
        ${mkStack(ex0180, 58, "Phosphate (0.069)")}
        ${legend}
      </svg>`;
    }

    renderTable("cs2-exceed-table",
      [
        { label: "Determinand", fn: (r) => r.Determinand_Name || (r.Threshold_Raw === 11.3 ? "Nitrate" : "Orthophosphate") },
        { label: "Tier", key: "Tier" },
        { label: "Reaches", key: "Reach_Count", num: true },
        { label: "Share %", key: "Reach_Pct", num: true, digits: 1 }
      ],
      [...ex0117, ...ex0180]
    );

    renderTable("cs2-cv-table",
      [
        { label: "Determinand", key: "Determinand_Name" },
        { label: "Domain", key: "Domain" },
        { label: "Model", key: "Model" },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      [...(t.cv_whpt || []), ...(t.cv_0180 || [])].filter((r) => String(r.Model).includes("OLS") || String(r.Model).includes("Combined")),
      (r) => String(r.Model).includes("Combined") ? "highlight-row" : ""
    );

    // --- CS3: Spatio-Temporal ---
    renderTable("cs3-cv-table",
      [
        { label: "Period Slice", key: "Period" },
        { label: "n_obs", key: "n_obs", num: true },
        { label: "Model", fn: () => "S3N Combined" },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "MAE", key: "MAE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      [...(t.st_year_cv || []), ...(t.st_season_cv || [])],
      (r) => String(r.Period).includes("winter") ? "highlight-row" : ""
    );

    // --- CS4: Co-Kriging & Multi-Source ---
    const ckEng0117 = (t.cokrig_cv || []).filter((r) => r.Domain === "All England RSN Network" && r.Target_Code === "0117");
    const ckUni = ckEng0117.find((r) => String(r.Model).includes("Univariate")) || {};
    const ckAux = ckEng0117.find((r) => String(r.Model).includes("Aux-Aided")) || {};
    const rssComb = (t.cv_rss || []).find((r) => r.Domain === "All England RSS Combined Network" && String(r.Model).includes("Combined")) || {};

    const kpi4 = document.getElementById("cs4-kpis");
    if (kpi4) {
      kpi4.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Co-Kriging Verdict</div><div class="kpi-val pos">8 / 8 Strong</div><div class="kpi-note">Outperformed Univariate across all basins &amp; targets</div></div>
        <div class="kpi-card"><div class="kpi-label">England 0117 Aux-Aided</div><div class="kpi-val pos">${fmt(ckAux.R2, 4)}</div><div class="kpi-note">Up from ${fmt(ckUni.R2, 4)} Univariate (-22.3% RMSE)</div></div>
        <div class="kpi-card"><div class="kpi-label">Median Aux-Aided Gain</div><div class="kpi-val pos">+0.2062 R²</div><div class="kpi-note">Median -16.6% RMSE across all 8 Domain x Target pairs</div></div>
        <div class="kpi-card"><div class="kpi-label">RSN+SSN Combined Network</div><div class="kpi-val">1,384 reaches</div><div class="kpi-note">CV R² = ${fmt(rssComb.R2, 4)} across 91,118 reaches (prefer RSN)</div></div>
      `;
    }

    const ckSvg = document.getElementById("cs4-cv-svg");
    if (ckSvg && ckEng0117.length) {
      const bars = (t.cokrig_cv || [])
        .filter((r) => r.Domain === "All England RSN Network" || r.Domain === "Severn Basin (CatchID=2786)")
        .map((r) => ({
          label: `${r.Domain.split(" ")[0]} ${r.Target_Code} ${r.Model.replace("Co-Kriging ", "")}`,
          val: r.R2 || 0,
          color: r.Model.includes("Aux-Aided") ? "#34d399" : (r.Model.includes("All-Holdout") ? "#38bdf8" : "#64748b")
        }));
      ckSvg.innerHTML = buildInlineBarSvg("DUAL CV R2 COMPARISON: UNIVARIATE vs ALL-HOLDOUT vs AUX-AIDED", bars, 0.85, "");
    }

    // Weight Penalty Analytical Curve SVG
    const wtSvg = document.getElementById("cs4-weight-svg");
    if (wtSvg) {
      const pts = [];
      for (let w = 0.1; w <= 1.001; w += 0.05) {
        const x = 60 + ((w - 0.1) / 0.9) * 480;
        const v = (1 / w - 1); // normalized units of tau^2
        const y = 135 - Math.min(9, v) / 9 * 105;
        pts.push(`${x.toFixed(1)},${y.toFixed(1)}`);
      }
      wtSvg.innerHTML = `<svg viewBox="0 0 600 165" width="100%" style="max-width:600px;">
        <text x="60" y="18" fill="var(--text-secondary)" font-size="11" font-weight="700">ANALYTIC PENALTY VARIANCE MULTIPLIER V_wt / tau^2 = (1 / w) - 1</text>
        <line x1="60" y1="135" x2="545" y2="135" stroke="var(--border-strong)" stroke-width="1.2"/>
        <line x1="60" y1="28" x2="60" y2="135" stroke="var(--border-strong)" stroke-width="1.2"/>
        <polyline fill="none" stroke="#fbbf24" stroke-width="2.6" points="${pts.join(" ")}" />
        <circle cx="540" cy="135" r="4.5" fill="#34d399"/><text x="480" y="125" fill="#34d399" font-size="11" font-weight="700">RSN w=1.0 (0x)</text>
        <circle cx="193" cy="96" r="4.5" fill="#38bdf8"/><text x="198" y="88" fill="#38bdf8" font-size="11" font-weight="700">WaterBlitz w=0.35 (1.86x)</text>
        <text x="60" y="153" fill="var(--text-muted)" font-size="11">w = 0.10 (9x tau^2 noise)</text>
        <text x="485" y="153" fill="var(--text-muted)" font-size="11">w = 1.00 (exact)</text>
      </svg>`;
    }

    renderTable("cs4-cv-table",
      [
        { label: "Domain", fn: (r) => String(r.Domain).replace(" (CatchID=80)", "").replace(" (CatchID=3055)", "").replace(" (CatchID=2786)", "") },
        { label: "Target", key: "Target" },
        { label: "CV Benchmark", key: "Model" },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      t.cokrig_cv || [],
      (r) => String(r.Model).includes("Aux-Aided") ? "highlight-row" : ""
    );

    renderTable("cs4-cov-table",
      [
        { label: "Domain", fn: (r) => String(r.Domain).split(" ")[0] },
        { label: "Pair", fn: (r) => `${r.Var_A} / ${r.Var_B}` },
        { label: "Corr TU", key: "Corr_TU", num: true, digits: 3 },
        { label: "Corr TD", key: "Corr_TD", num: true, digits: 3 },
        { label: "Corr Nugget", key: "Corr_Nugget", num: true, digits: 3 },
        { label: "Total Corr", key: "Corr_Total", num: true, digits: 3 }
      ],
      (t.cokrig_cov || []).filter((r) => r.Var_A !== r.Var_B)
    );

    // --- CS5: Audits, Bias, N:P & Snapping ---
    renderTable("cs5-bias-table",
      [
        { label: "Network / Code", fn: (r) => `${r.Run_Label} (${r.Determinand_Code})` },
        { label: "Obs Median", key: "Obs_Median_Raw", num: true, digits: 3 },
        { label: "Pred Median", key: "Pred_Median_Raw", num: true, digits: 3 },
        { label: "Ratio", key: "Pred_to_Obs_Med_Ratio", num: true, digits: 3 },
        { label: "MW p-value", fn: (r) => Number(r.MW_P_Value).toExponential(2), num: true }
      ],
      t.desc_bias || []
    );

    renderTable("cs5-tests-table",
      [
        { label: "Category", key: "Category" },
        { label: "Statistical Test", key: "Test_Name" },
        { label: "Metric", key: "Metric" },
        { label: "p-value", fn: (r) => typeof r.P_Value === "number" ? (r.P_Value < 0.001 ? r.P_Value.toExponential(2) : r.P_Value.toFixed(4)) : esc(r.P_Value), num: true },
        { label: "Sig", key: "Significance" }
      ],
      (t.comp_tests || []).slice(0, 10)
    );

    renderTable("cs5-snap-table",
      [
        { label: "Catchment", key: "Main_River" },
        { label: "RBD", key: "RBD_Name" },
        { label: "Included", key: "Stations_Included", num: true },
        { label: "Median (m)", key: "Median_Snap_m", num: true, digits: 2 },
        { label: "P95 (m)", key: "P95_Snap_m", num: true, digits: 2 },
        { label: "<=10m %", key: "Pct_Within_10m", num: true, digits: 1 }
      ],
      t.snapping_top || []
    );

    // --- CS6: WaterBlitz Citizen Science Fusion ---
    const kpi6 = document.getElementById("cs6-kpis");
    if (kpi6) {
      kpi6.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">England Reaches Snap-Fused</div><div class="kpi-val">8,346</div><div class="kpi-note">33,013 RSN+SSN+WB stations &le;350 m across 268 catchments</div></div>
        <div class="kpi-card"><div class="kpi-label">Co-Located Reach Dilution Fix</div><div class="kpi-val pos">2.62 &rarr; 3.43</div><div class="kpi-note">NO3-N mg/L median on 384 shared reaches (-46.7% SE)</div></div>
        <div class="kpi-card"><div class="kpi-label">National Kriging SE Change</div><div class="kpi-val pos">-14.5% / -23.4%</div><div class="kpi-note">0117 (0.931&rarr;0.795) &amp; 0180 (1.102&rarr;0.844) all 95,773 reaches</div></div>
        <div class="kpi-card"><div class="kpi-label">PO4 Uncertain Tier Reaches</div><div class="kpi-val pos">51.0% &rarr; 41.7%</div><div class="kpi-note">48,802 &rarr; 39,969 reaches; N-P Corr_Total +0.299 &rarr; +0.431</div></div>
      `;
    }

    renderTable("cs6-reach-table",
      [
        { label: "Var", key: "Target_Code" },
        { label: "Stratum", key: "Stratum" },
        { label: "n", key: "n_reaches", num: true },
        { label: "Eq Med", key: "Equal_Med_Raw", num: true, digits: 4 },
        { label: "Fav Med", key: "Favoured_Med_Raw", num: true, digits: 4 },
        { label: "|Shift|", key: "Med_Abs_Shift", num: true, digits: 4 },
        { label: "Pearson r", key: "Pearson_r_Scale", num: true, digits: 4 },
        { label: "Eq SE", key: "Equal_Mean_SE", num: true, digits: 4 },
        { label: "Fav SE", key: "Favoured_Mean_SE", num: true, digits: 4 },
        { label: "SE Delta", fn: (r) => `${r.SE_Delta_Pct >= 0 ? "+" : ""}${fmt(r.SE_Delta_Pct, 2)}%` },
        { label: "Tier Match", fn: (r) => `${fmt(r.Tier_Match_Pct, 2)}%` }
      ],
      t.cs6_reach_comp || [],
      (r) => String(r.Stratum).includes("Co-Located") || String(r.Stratum).includes("WaterBlitz-Only") ? "highlight-row" : ""
    );

    renderTable("cs6-cv-table",
      [
        { label: "Domain", key: "Domain_Short" },
        { label: "Var", key: "Target_Code" },
        { label: "Model", fn: (r) => String(r.Model).replace("Co-Kriging ", "") },
        { label: "n", key: "n_obs", num: true },
        { label: "Eq RMSE", key: "RMSE_Equal", num: true, digits: 4 },
        { label: "Fav RMSE", key: "RMSE_Favoured", num: true, digits: 4 },
        { label: "Eq R²", key: "R2_Equal", num: true, digits: 4 },
        { label: "Fav R²", key: "R2_Favoured", num: true, digits: 4 }
      ],
      (t.cs6_cv_comp || []).filter((r) => r.Domain_Short === "All England" || r.Model === "Co-Kriging (Aux-Aided)"),
      (r) => String(r.Model).includes("Aux-Aided") ? "highlight-row" : ""
    );

    renderTable("cs6-cov-table",
      [
        { label: "Domain", key: "Domain_Short" },
        { label: "Eq TU r", key: "Corr_TU_Equal", num: true, digits: 4 },
        { label: "Fav TU r", key: "Corr_TU_Favoured", num: true, digits: 4 },
        { label: "Eq TD r", key: "Corr_TD_Equal", num: true, digits: 4 },
        { label: "Fav TD r", key: "Corr_TD_Favoured", num: true, digits: 4 },
        { label: "Eq Nug r", key: "Corr_Nugget_Equal", num: true, digits: 4 },
        { label: "Fav Nug r", key: "Corr_Nugget_Favoured", num: true, digits: 4 },
        { label: "Eq Total r", key: "Corr_Total_Equal", num: true, digits: 4 },
        { label: "Fav Total r", key: "Corr_Total_Favoured", num: true, digits: 4 }
      ],
      t.cs6_cov_comp || [],
      (r) => String(r.Domain_Short) === "All England" ? "highlight-row" : ""
    );

    // --- CS7: Data Transforms Comparison ---
    const kpi7 = document.getElementById("cs7-kpis");
    if (kpi7) {
      kpi7.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Profile-MLE Box-Cox &lambda;</div><div class="kpi-val pos">0.50 (Sqrt)</div><div class="kpi-note">MLE converges to power 0.50, cutting skew 1.001 &rarr; 0.118</div></div>
        <div class="kpi-card"><div class="kpi-label">Raw Positivity Violations</div><div class="kpi-val warn">6,349 (6.99%)</div><div class="kpi-note">47 negative CV stations; min -4.817 mg/L vs. 0% under Sqrt/Log/Box-Cox</div></div>
        <div class="kpi-card"><div class="kpi-label">Best Model-Scale vs. Raw Error</div><div class="kpi-val pos">0.6784 / 1.57</div><div class="kpi-note">Log tops model R&sup2;; Sqrt &amp; Box-Cox minimize Raw MAE (1.57 mg/L)</div></div>
        <div class="kpi-card"><div class="kpi-label">Log Median Bias vs. Rank</div><div class="kpi-val amber">-1.01 mg/L</div><div class="kpi-note">exp(&mu;) Jensen shift vs. Sqrt/Box-Cox &rho;=0.982 rank agreement</div></div>
      `;
    }

    renderTable("cs7-summary-table",
      [
        { label: "Transform", key: "Transform" },
        { label: "Lambda", key: "Lambda" },
        { label: "Model R²", key: "R2_ModelSc", num: true, digits: 4 },
        { label: "Raw R²", key: "R2_Raw_mgL", num: true, digits: 4 },
        { label: "Raw RMSE", key: "RMSE_Raw_mgL", num: true, digits: 4 },
        { label: "Raw MAE", key: "MAE_Raw_mgL", num: true, digits: 4 },
        { label: "Raw Bias", key: "Bias_Raw_mgL", num: true, digits: 4 },
        { label: "CV Skew", key: "CV_Resid_Skew", num: true, digits: 3 },
        { label: "Neg Reaches", fn: (r) => `${Number(r.Neg_Reaches).toLocaleString()} (${fmt(r.Neg_Reach_Pct, 2)}%)`, num: true },
        { label: "Min mg/L", key: "Eng_Min_mgL", num: true, digits: 3 },
        { label: "Uncertain %", fn: (r) => `${fmt(r.Uncertain_Pct, 2)}%`, num: true },
        { label: "rho vs Log", key: "Spearman_vs_Log", num: true, digits: 4 }
      ],
      t.cs7_transform_summary || [],
      (r) => String(r.Key) === "boxcox" || String(r.Key) === "sqrt" ? "highlight-row" : ""
    );

    renderTable("cs7-cv-table",
      [
        { label: "Transform", key: "Run_Label" },
        { label: "Domain", key: "Canon_Domain" },
        { label: "n Obs", key: "n_obs", num: true },
        { label: "Model", key: "Canon_Model" },
        { label: "CV R²", key: "R2", num: true, digits: 4 },
        { label: "CV RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "CV MAE", key: "MAE", num: true, digits: 4 },
        { label: "Gain vs OLS", fn: (r) => `${r.Delta_R2_vs_OLS >= 0 ? "+" : ""}${fmt(r.Delta_R2_vs_OLS, 4)}`, num: true }
      ],
      (t.cs7_cv_comp || []).filter((r) => String(r.Canon_Model) === "S3N Combined"),
      (r) => String(r.Canon_Domain).includes("All England") ? "highlight-row" : ""
    );

    // --- CS8: Optimal Monitoring Network Design ---
    const kpi8 = document.getElementById("cs8-kpis");
    if (kpi8) {
      kpi8.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">England Variance Drop (+50)</div><div class="kpi-val pos">80.74%</div><div class="kpi-note">Evaluated ORN variance reduced across 90,876 reaches</div></div>
        <div class="kpi-card"><div class="kpi-label">National Mean Kriging SE</div><div class="kpi-val pos">0.796 &rarr; 0.323</div><div class="kpi-note">-59.5% SE reduction on log scale (1,210 &rarr; +50 sites)</div></div>
        <div class="kpi-card"><div class="kpi-label">Top Single Placement (#1)</div><div class="kpi-val">Riv 18175</div><div class="kpi-note">Humber CatchID 80 (Strahler O1): +4.58% variance drop</div></div>
        <div class="kpi-card"><div class="kpi-label">Per-Basin Gains (+50 Each)</div><div class="kpi-val pos">32.3 / 43.6 / 43.0%</div><div class="kpi-note">Humber / Thames / Severn independent variance drops</div></div>
      `;
    }

    renderTable("cs8-design-summary-table",
      [
        { label: "Domain (K=50 Added Sites)", key: "Domain" },
        { label: "Rank-50 Reach", fn: (r) => `Riv ${r.RivID} (C${r.CatchID})` },
        { label: "Order", key: "Strahler", num: true },
        { label: "Cum Var Drop", fn: (r) => `${fmt(r.Cum_Pct_Reduced, 2)}%`, num: true },
        { label: "Mean SE Final", key: "Mean_SE_After", num: true, digits: 4 }
      ],
      t.cs2_design_summary || [],
      (r) => String(r.Domain).includes("All England") ? "highlight-row" : ""
    );

    renderTable("cs8-design-table",
      [
        { label: "Rank", key: "Rank", num: true },
        { label: "RivID", key: "RivID", num: true },
        { label: "CatchID", key: "CatchID", num: true },
        { label: "Order", key: "Strahler", num: true },
        { label: "Dist2Mth (km)", key: "Dist2Mth_km", num: true, digits: 1 },
        { label: "US Accum (km)", key: "US_Accum_km", num: true, digits: 1 },
        { label: "Elev (m)", key: "Elevation", num: true, digits: 1 },
        { label: "Local SE", key: "Current_SE_sc", num: true, digits: 4 },
        { label: "Var Reduced", key: "Var_Reduced", num: true, digits: 2 },
        { label: "+dV %", fn: (r) => `+${fmt(r.Pct_Var_Reduced, 2)}%`, num: true },
        { label: "Cum %", fn: (r) => `${fmt(r.Cum_Pct_Reduced, 2)}%`, num: true },
        { label: "Mean SE", key: "Mean_SE_After", num: true, digits: 4 }
      ],
      t.cs2_design_eng || [],
      (r) => Number(r.Rank) <= 5 ? "highlight-row" : ""
    );
  }

  function initFeatureCatalog() {
    const grid = document.getElementById("featureCatalogGrid");
    const btns = document.querySelectorAll("#featureFilterBtns .fig-toggle-btn");
    if (!grid) return;

    function draw(cat) {
      const items = cat === "all" ? FEATURE_CATALOG : FEATURE_CATALOG.filter((f) => f.cat === cat);
      grid.innerHTML = items.map((f) => `
        <div class="feature-item">
          <span class="pill ${f.cat === "major" ? "pill-green" : "pill-amber"}">${esc(f.badge)}</span>
          <h4>${esc(f.title)}</h4>
          <p>${esc(f.desc)}</p>
          <div style="margin-top:8px;">
            <button class="copy-btn" data-jump-tab="${esc(f.tab)}">Open Case Study &rarr;</button>
          </div>
        </div>
      `).join("");
      grid.querySelectorAll("[data-jump-tab]").forEach((b) => {
        b.addEventListener("click", () => switchTab(b.dataset.jumpTab));
      });
    }

    draw("all");
    btns.forEach((btn) => {
      btn.addEventListener("click", () => {
        btns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        draw(btn.dataset.fcat);
      });
    });
  }

  function initSearch() {
    const input = document.getElementById("featureSearch");
    const drop = document.getElementById("searchDropdown");
    if (!input || !drop) return;

    const index = [
      ...FEATURE_CATALOG.map((f) => ({ title: f.title, sub: f.desc, tab: f.tab })),
      ...Object.entries(FIGURE_META).map(([k, m]) => ({
        title: m.title,
        sub: m.note,
        tab: k === "no3_z_thames_severn" || k.startsWith("cs1") || k.startsWith("no3_gam") || k === "no3_england" || k === "no3_humber" || k === "no3_thames_severn" || k === "no3_diag" ? "cs1"
          : k.startsWith("whpt") || k.includes("exceed") || k.includes("no3_z") ? "cs2"
          : k.startsWith("st_") ? "cs3"
          : k.startsWith("cokrig") || k.startsWith("rss_") ? "cs4"
          : k.startsWith("wb_") ? "cs6"
          : k.startsWith("tr_") ? "cs7"
          : k.startsWith("opt_design_") ? "cs8" : "cs5"
      }))
    ];

    input.addEventListener("input", () => {
      const q = input.value.trim().toLowerCase();
      if (!q) { drop.classList.remove("open"); return; }
      const hits = index.filter((item) => (item.title + " " + item.sub).toLowerCase().includes(q)).slice(0, 8);
      if (!hits.length) {
        drop.innerHTML = `<div class="search-item">No matches for "${esc(q)}"</div>`;
      } else {
        drop.innerHTML = hits.map((h) => `
          <button class="search-item" data-hit-tab="${esc(h.tab)}">
            <strong>${esc(h.title)}</strong>
            <small>${esc(h.sub)}</small>
          </button>
        `).join("");
        drop.querySelectorAll("[data-hit-tab]").forEach((btn) => {
          btn.addEventListener("click", () => {
            switchTab(btn.dataset.hitTab);
            drop.classList.remove("open");
            input.value = "";
          });
        });
      }
      drop.classList.add("open");
    });

    document.addEventListener("click", (e) => {
      if (!e.target.closest(".search-box")) drop.classList.remove("open");
    });
  }

  function initWorkbench() {
    const ids = ["wbMode", "wbTarget", "wbDataset", "wbSource", "wbTrend", "wbKernel", "wbTransform", "wbThreshold", "wbTemporal", "wbPreferSource", "wbSourceWeight", "wbSourceSe", "wbDesignSites", "wbFormat", "wbOmp"];
    const els = {};
    ids.forEach((id) => { els[id] = document.getElementById(id); });
    const outEl = document.getElementById("wbOutputCmd");
    const notesEl = document.getElementById("wbNotesList");
    if (!outEl) return;

    function rebuild() {
      const mode = els.wbMode.value;
      const target = els.wbTarget.value.trim() || "0117";
      const ds = els.wbDataset.value;
      const src = els.wbSource.value.trim() || "rsn";
      const trend = els.wbTrend.value;
      const kernel = els.wbKernel.value;
      const trans = els.wbTransform.value;
      const thr = els.wbThreshold.value;
      const temp = els.wbTemporal.value;
      const pref = els.wbPreferSource.value.trim();
      const wt = els.wbSourceWeight.value.trim();
      const se = els.wbSourceSe.value.trim();
      const des = els.wbDesignSites ? Number(els.wbDesignSites.value || 0) : 0;
      const fmtVal = els.wbFormat.value;
      const omp = els.wbOmp.value;

      const notes = [];
      let cmd = "";

      if (mode === "prep") {
        cmd = `python3 scripts/01_preprocess_orn.py --source=${src}`;
        if (ds !== "auto") cmd += ` --dataset=${ds}`;
        if (trans !== "default") cmd += ` --transform=${trans}`;
        notes.push("Outputs snapped topology, log_AFV, ROS-imputed samples & network metadata to data/processed/.");
      } else if (mode === "cokrig") {
        const ckTargets = target.includes(",") ? target : `${target},0180`;
        const method = trend === "colocated" ? "colocated" : "lmc";
        cmd = `Rscript scripts/run_s3n_cokriging.R --targets=${ckTargets} --source=${src} --method=${method}`;
        if (ds !== "auto") cmd += ` --dataset=${ds}`;
        if (kernel !== "exponential") cmd += ` --kernel=${kernel}`;
        if (thr === "uk") cmd += ` --threshold=uk`;
        if (pref) cmd += ` --prefer-source=${pref}`;
        if (wt) cmd += ` --source-weight=${wt}`;
        if (se) cmd += ` --source-se=${se}`;
        if (des > 0) {
          cmd += ` --design-sites=${des}`;
          notes.push(`Runs OpenMP C++ greedy active learning (C_s3n_optimal_design) to rank top ${des} sites & plot fig_optimal_design_* maps.`);
        }
        if (fmtVal === "png_300") cmd += ` --dpi=300`;
        else cmd += ` --fig-format=${fmtVal}`;
        cmd += ` --omp=${omp}`;
        notes.push("Evaluates Univariate S3N vs Co-Kriging (All-Holdout) vs Co-Kriging (Aux-Aided) CV benchmarks automatically.");
        if (wt || se) notes.push("Routes through OpenMP weighted C++ solver C_s3n_multivariate_nngp_weighted.");
      } else {
        cmd = `Rscript scripts/run_s3n_determinand.R ${target}`;
        if (ds !== "auto") cmd += ` --dataset=${ds}`;
        if (src !== "rsn") cmd += ` --source=${src}`;
        if (["spline", "gam", "rf", "gbm"].includes(trend)) cmd += ` --trend=${trend}`;
        if (kernel !== "exponential") cmd += ` --kernel=${kernel}`;
        if (trans !== "default") cmd += ` --transform=${trans}`;
        if (thr === "uk") cmd += ` --threshold=uk`;
        if (temp === "by_season") cmd += ` --by=season --st-borrow`;
        else if (temp === "by_year") cmd += ` --by=year --st-borrow`;
        else if (temp === "winter") cmd += ` --season=winter`;
        if (target === "whpt_aspt" || ds === "inverts") cmd += ` --reverse-color`;
        if (pref) cmd += ` --prefer-source=${pref}`;
        if (wt) cmd += ` --source-weight=${wt}`;
        if (se) cmd += ` --source-se=${se}`;
        if (des > 0) {
          cmd += ` --design-sites=${des}`;
          notes.push(`Runs OpenMP C++ greedy active learning (C_s3n_optimal_design) to rank top ${des} sites & plot fig_optimal_design_* maps.`);
        }
        if (fmtVal === "png_300") cmd += ` --dpi=300`;
        else cmd += ` --fig-format=${fmtVal}`;
        cmd += ` --omp=${omp}`;
        notes.push("Runs England + Humber + Thames + Severn models with 5-fold hydrological/random CV and maps.");
        if (fmtVal === "svg") notes.push("Tip: Use --dpi=300 PNG for 90k-reach national maps and --fig-format=svg for diagnostic bar/forest charts.");
      }
      outEl.textContent = cmd;
      notesEl.innerHTML = notes.map((n) => `<li>${esc(n)}</li>`).join("");
    }

    Object.values(els).forEach((el) => {
      if (el) el.addEventListener("change", rebuild);
      if (el && el.tagName === "INPUT") el.addEventListener("input", rebuild);
    });

    const presets = {
      nitrate: { wbMode: "det", wbTarget: "0117", wbDataset: "wq", wbSource: "rsn", wbTrend: "linear", wbKernel: "auto", wbTransform: "log", wbThreshold: "uk", wbTemporal: "all", wbPreferSource: "", wbSourceWeight: "", wbSourceSe: "", wbFormat: "png_300" },
      ecology: { wbMode: "det", wbTarget: "whpt_aspt", wbDataset: "inverts", wbSource: "rsn", wbTrend: "gam", wbKernel: "exponential", wbTransform: "none", wbThreshold: "none", wbTemporal: "all", wbPreferSource: "", wbSourceWeight: "", wbSourceSe: "", wbFormat: "png_300" },
      seasonal: { wbMode: "det", wbTarget: "0117", wbDataset: "wq", wbSource: "rsn", wbTrend: "linear", wbKernel: "exponential", wbTransform: "log", wbThreshold: "none", wbTemporal: "by_season", wbPreferSource: "", wbSourceWeight: "", wbSourceSe: "", wbDesignSites: "0", wbFormat: "png_300" },
      waterblitz: { wbMode: "cokrig", wbTarget: "0117,0180", wbDataset: "wq", wbSource: "rsn+ssn+waterblitz_tobit", wbTrend: "lmc", wbKernel: "exponential", wbTransform: "log", wbThreshold: "uk", wbTemporal: "all", wbPreferSource: "rsn,ssn", wbSourceWeight: "rsn:1,ssn:1,waterblitz:0.25", wbSourceSe: "", wbDesignSites: "0", wbFormat: "png_300" },
      design: { wbMode: "det", wbTarget: "0117", wbDataset: "wq", wbSource: "rsn", wbTrend: "linear", wbKernel: "exponential", wbTransform: "log", wbThreshold: "uk", wbTemporal: "all", wbPreferSource: "", wbSourceWeight: "", wbSourceSe: "", wbDesignSites: "50", wbFormat: "png_300" }
    };

    document.querySelectorAll("#wbPresetBtns [data-preset]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const p = presets[btn.dataset.preset];
        if (!p) return;
        Object.entries(p).forEach(([k, v]) => { if (els[k]) els[k].value = v; });
        rebuild();
      });
    });

    rebuild();
  }

  function initCopyButtons() {
    document.querySelectorAll("[data-copy-target]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const target = document.getElementById(btn.dataset.copyTarget);
        if (!target) return;
        const txt = target.textContent;
        const done = () => {
          const orig = btn.textContent;
          btn.textContent = "Copied!";
          setTimeout(() => { btn.textContent = orig; }, 1400);
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(txt).then(done).catch(() => fallbackCopy(txt, done));
        } else {
          fallbackCopy(txt, done);
        }
      });
    });
  }

  function fallbackCopy(txt, cb) {
    const ta = document.createElement("textarea");
    ta.value = txt;
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    try { document.execCommand("copy"); cb(); } catch (e) {}
    document.body.removeChild(ta);
  }

  function initLightbox() {
    const modal = document.getElementById("figModal");
    const img = document.getElementById("modalImg");
    const bodyEl = modal ? modal.querySelector(".modal-body") : null;
    const title = document.getElementById("modalTitle");
    const metaPill = document.getElementById("modalMeta");
    const footer = document.getElementById("modalFooter");
    const zoomBtn = document.getElementById("modalZoomBtn");
    const closeBtn = document.getElementById("modalCloseBtn");
    const prevBtn = document.getElementById("modalPrevBtn");
    const nextBtn = document.getElementById("modalNextBtn");
    if (!modal) return;

    let gallery = [];
    let currentIdx = 0;

    function getActiveTabGallery(triggerEl) {
      const panel = triggerEl ? triggerEl.closest(".tab-panel") : document.querySelector(".tab-panel.active");
      const keys = [];
      const seen = new Set();
      if (panel) {
        panel.querySelectorAll(".fig-toggle-bar .fig-toggle-btn[data-fig], [data-open-fig]").forEach((el) => {
          const k = el.dataset.fig || el.dataset.openFig;
          if (k && !seen.has(k) && (DATA.figures || {})[k]) {
            seen.add(k);
            keys.push(k);
          }
        });
      }
      return keys;
    }

    function showAtIndex(idx) {
      if (!gallery.length) return;
      currentIdx = ((idx % gallery.length) + gallery.length) % gallery.length;
      const key = gallery[currentIdx];
      const fig = (DATA.figures || {})[key];
      const info = FIGURE_META[key] || { title: key, note: "" };
      if (!fig) return;
      img.src = fig.path;
      img.alt = info.title;
      img.classList.remove("zoomed");
      if (bodyEl) {
        bodyEl.scrollTop = 0;
        bodyEl.scrollLeft = 0;
      }
      title.textContent = info.title;
      const counter = gallery.length > 1 ? `${currentIdx + 1} / ${gallery.length} · ` : "";
      metaPill.textContent = `${counter}${fig.format} · ${fig.size_kb} KB`;
      footer.textContent = info.note;
      const single = gallery.length <= 1;
      if (prevBtn) prevBtn.hidden = single;
      if (nextBtn) nextBtn.hidden = single;
    }

    function stepGallery(delta) {
      if (gallery.length <= 1) return;
      showAtIndex(currentIdx + delta);
    }

    document.addEventListener("click", (e) => {
      const frame = e.target.closest("[data-open-fig]");
      if (!frame) return;
      const key = frame.dataset.openFig;
      gallery = getActiveTabGallery(frame);
      if (!gallery.includes(key)) {
        gallery = (DATA.figures || {})[key] ? [key] : [];
      }
      if (!gallery.length) return;
      showAtIndex(gallery.indexOf(key));
      modal.classList.add("open");
    });

    if (prevBtn) prevBtn.addEventListener("click", () => stepGallery(-1));
    if (nextBtn) nextBtn.addEventListener("click", () => stepGallery(1));
    zoomBtn.addEventListener("click", () => img.classList.toggle("zoomed"));
    img.addEventListener("click", () => img.classList.toggle("zoomed"));
    closeBtn.addEventListener("click", () => modal.classList.remove("open"));
    modal.addEventListener("click", (e) => { if (e.target === modal) modal.classList.remove("open"); });
    document.addEventListener("keydown", (e) => {
      if (!modal.classList.contains("open")) return;
      if (e.key === "Escape") modal.classList.remove("open");
      else if (e.key === "ArrowLeft") { e.preventDefault(); stepGallery(-1); }
      else if (e.key === "ArrowRight") { e.preventDefault(); stepGallery(1); }
    });
  }

  function initThemeToggle() {
    const btn = document.getElementById("themeToggleBtn");
    const root = document.documentElement;
    function apply(theme) {
      root.setAttribute("data-theme", theme);
      btn.textContent = theme === "dark" ? "☀️ Light Mode" : "🌙 Dark Mode";
      try { localStorage.setItem("modelling_dashboard_theme", theme); } catch (e) {}
    }
    let saved = "dark";
    try { saved = localStorage.getItem("modelling_dashboard_theme") || "dark"; } catch (e) {}
    apply(saved);
    btn.addEventListener("click", () => apply(root.getAttribute("data-theme") === "dark" ? "light" : "dark"));
  }

  // Boot
  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => switchTab(btn.dataset.tab));
  });
  document.querySelectorAll(".pipe-node[data-jump]").forEach((node) => {
    node.addEventListener("click", () => switchTab(node.dataset.jump));
  });

  initThemeToggle();
  initFigureToggles();
  populateCaseStudyData();
  initFeatureCatalog();
  initSearch();
  initWorkbench();
  initCopyButtons();
  initLightbox();

  const initialHash = window.location.hash.replace("#", "");
  if (initialHash && document.getElementById("tab-" + initialHash)) {
    switchTab(initialHash, false);
  }
})();
