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
      title: "4-Panel Nitrate Diagnostics, Cross-Validation & Parameter CIs",
      note: "Compares OLS, Tail-Down, Tail-Up, and Combined NNGP CV accuracy alongside bootstrap covariate & variance intervals."
    },
    no3_gam_diag: {
      title: "4-Panel GAM + S3N Diagnostics & Cross-Validation",
      note: "Diagnostic summary for non-linear --trend=gam Nitrate model (CV R2 = 0.6743)."
    },
    ntaxa_england: {
      title: "England WHPT N-Taxa (whpt_ntaxa) — Joint TU+TD+Euc Surface & SE (--reverse-color)",
      note: "90,876 reaches, 1,160 RSN stations (3,408 samples). Inverted ecological palette: blue = high taxon richness, red = low richness."
    },
    ntaxa_gam_england: {
      title: "England WHPT N-Taxa — Non-Linear REML GAM + S3N Trend (--reverse-color)",
      note: "Penalized cubic regression splines capture non-linear elevation and catchment-size richness thresholds before NNGP spatial kriging (CV R2 = 0.2664)."
    },
    ntaxa_humber: {
      title: "Humber Basin (CatchID=80) — Stream Macroinvertebrate N-Taxa",
      note: "211 stations across 14,048 reaches. Combined CV R2 = 0.1625 (vs 0.0987 OLS), RMSE = 5.1285 taxa."
    },
    ntaxa_thames_severn: {
      title: "Thames (3055) & Severn (2786) Basins — WHPT N-Taxa & Kriging SE",
      note: "Severn CV R2 rises from 0.0200 (OLS) to 0.1750 (Joint S3N), while Thames exhibits localized heterogeneity at n=107."
    },
    ntaxa_diag: {
      title: "4-Panel WHPT N-Taxa Diagnostics, CV & Parameter CIs",
      note: "Compares OLS (R2 = 0.1371), Tail-Up (0.1834), Tail-Down (0.2391), and Combined Joint NNGP (0.2648) across England."
    },
    ntaxa_gam_diag: {
      title: "4-Panel GAM + S3N WHPT N-Taxa Diagnostics & CV",
      note: "Diagnostic summary for non-linear --trend=gam WHPT N-Taxa model (All-England CV R2 = 0.2664, RMSE = 5.3926 taxa)."
    },
    ntaxa_z_england: {
      title: "England WHPT N-Taxa Robust Z-Score Map (--norm=robust_zscore --reverse-color)",
      note: "Expresses richness deviations as (x - median) / (1.4826 * MAD) to highlight ecological richness anomalies nationally."
    },
    ntaxa_z_humber: {
      title: "Humber Basin WHPT N-Taxa Robust Z-Score Map (--norm=robust_zscore)",
      note: "High-resolution Humber sub-basin view contrasting upland Pennine richness (+Z) against lowland/urban corridors (-Z)."
    },
    ntaxa_z_thames_severn: {
      title: "Thames & Severn WHPT N-Taxa Robust Z-Score Maps (--norm=robust_zscore)",
      note: "Contrasts Thames chalk/urban richness anomalies against Severn upland-to-lowland gradients on a dimensionless robust scale."
    },
    eut_cokrig_rmni: {
      title: "Co-Kriged Macrophyte RMNI (Factor-LMC Joint with P & N)",
      note: "Borrows spatial strength from 1,210 Orthophosphate & Nitrate stations to predict RMNI (All-England Aux-Aided R2 = 0.7580; Thames R2 rescued from -0.2671 to 0.3165)."
    },
    eut_cokrig_0180: {
      title: "Co-Kriged Orthophosphate (0180) in P-N-RMNI LMC",
      note: "Jointly modelled log-Orthophosphate surface (Aux-Aided CV R2 = 0.6560 vs 0.4665 Univariate)."
    },
    eut_cokrig_xcorr: {
      title: "LMC Residual Cross-Correlation Matrix (0180, 0117, RMNI)",
      note: "Strong Tail-Down residual correlations between Orthophosphate-RMNI (+0.752) and Nitrate-RMNI (+0.790) power cross-variable borrowing."
    },
    eut_rmni_linear: {
      title: "England Macrophyte RMNI — Linear Joint S3N (300 DPI)",
      note: "All-England CV R2 = 0.7398 (RMSE = 0.7139) vs 0.6289 OLS; higher RMNI indicates more eutrophic macrophyte communities."
    },
    eut_rmni_gam: {
      title: "England Macrophyte RMNI — Non-Linear GAM + S3N (300 DPI)",
      note: "Highest national RMNI accuracy: CV R2 = 0.7607 (RMSE = 0.6846), capturing non-linear slope and drainage-area thresholds."
    },
    eut_rmni_diag: {
      title: "4-Panel RMNI Linear S3N Diagnostics & Parameter CIs",
      note: "CV and parameter summary for Linear Joint S3N on macrophytes:rmni."
    },
    eut_rmni_gam_diag: {
      title: "4-Panel RMNI GAM + S3N Diagnostics & Parameter CIs",
      note: "Diagnostic summary for --trend=gam RMNI model (All-England CV R2 = 0.7607)."
    },
    eut_do_season_grid: {
      title: "England Dissolved Oxygen % Saturation (9901) Across 4 Seasons (--reverse-color)",
      note: "Shared-scale [69.3%, 109.2%] seasonal DO maps with precision-weighted temporal borrowing (rho_t = 0.150)."
    },
    eut_do_summer: {
      title: "Summer (JJA) Dissolved Oxygen % Saturation & Kriging SE (--reverse-color)",
      note: "Low-flow summer oxygen sag zones (red); instantaneous DO has higher diel variance (CV R2 = 0.1574) than biological RMNI."
    },
    eut_p_rmni_maps: {
      title: "Side-by-Side Orthophosphate (0180) vs. Macrophyte RMNI Predictions",
      note: "90,876 matched reaches show Pearson r = 0.9005 (Spearman rho = 0.8542) and 86.55% 4-quadrant concordance (Odds Ratio = 41.40)."
    },
    eut_p_rmni_rel: {
      title: "Orthophosphate vs. RMNI Reach-Level Relationship & Strata",
      note: "Thames Basin has 81.86% High-P / High-RMNI reaches vs. 96.71% Low-P / Low-RMNI in Uplands (>= 150 m)."
    },
    nw_no3_map: {
      title: "North West RBD (RBD=12) — Nitrate (0117) Predictions & Kriging SE",
      note: "12,465 reaches, 136 RSN stations. Joint S3N CV R2 = 0.6702 (vs 0.6523 OLS), RMSE = 0.8429 log mg/L."
    },
    nw_no3_exceed: {
      title: "North West RBD — 5-Tier Nitrate Exceedance Risk (11.3 mg/L as N)",
      note: "93.8% Confidently Compliant (11,686 reaches), 5.7% Likely Compliant (704), 0.6% Uncertain (75), 0% Exceeding."
    },
    nw_no3_winter: {
      title: "North West RBD — Winter (DJF) Nitrate Map (--by=season --st-borrow)",
      note: "Shared-scale [0.033, 10.000] mg/L N winter surface with temporal borrowing (rho_t = 0.363); CV R2 = 0.6197."
    },
    nw_no3_spring: {
      title: "North West RBD — Spring (MAM) Nitrate Map (--by=season --st-borrow)",
      note: "Shared-scale spring Nitrate surface across North West RBD; Joint S3N CV R2 = 0.6576."
    },
    nw_no3_summer: {
      title: "North West RBD — Summer (JJA) Nitrate Map (--by=season --st-borrow)",
      note: "Shared-scale summer Nitrate surface across North West RBD; Joint S3N CV R2 = 0.6704."
    },
    nw_no3_autumn: {
      title: "North West RBD — Autumn (SON) Nitrate Map (--by=season --st-borrow)",
      note: "Shared-scale autumn Nitrate surface across North West RBD; Joint S3N CV R2 = 0.6799."
    },
    nw_ntaxa_map: {
      title: "North West RBD (RBD=12) — WHPT N-Taxa Richness (--reverse-color)",
      note: "12,465 reaches, 131 RSN stations. Joint S3N lifts CV R2 from 0.1252 (OLS) to 0.1950 (RMSE = 5.1640 taxa)."
    },
    nw_snap_ribble: {
      title: "River Ribble (CatchID=161) — Station-to-Reach Snapping Diagnostics",
      note: "2,625 reaches, 31 included stations (2 excluded > 250 m), median snap 1.88 m (87.1% <= 10 m)."
    },
    nw_snap_mersey: {
      title: "River Mersey (CatchID=177) — Station-to-Reach Snapping Diagnostics",
      note: "1,776 reaches, 23 included stations (0 excluded), median snap 2.15 m (95.7% <= 10 m)."
    },
    nw_snap_weaver: {
      title: "River Weaver (CatchID=173) — Station-to-Reach Snapping Diagnostics",
      note: "856 reaches, 23 included stations (5 excluded > 250 m), median snap 2.18 m (100.0% <= 10 m)."
    },
    nw_snap_lune: {
      title: "River Lune (CatchID=157) — Station-to-Reach Snapping Diagnostics",
      note: "1,655 reaches, 17 included stations (0 excluded), median snap 3.48 m (76.5% <= 10 m)."
    },
    nw_snap_derwent: {
      title: "River Derwent (CatchID=138) — Station-to-Reach Snapping Diagnostics",
      note: "994 reaches in Cumbria/Lake District evaluated for station-to-reach snapping accuracy."
    },
    nw_opt_design: {
      title: "North West RBD (RBD=12) Optimal Monitoring Design & Gap Map (K=50)",
      note: "Top 50 greedy C++ active-learning sites cut North West kriging variance by 45.95% (Mean SE 0.6956 -> 0.4887 log mg/L)."
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
      title: "WHPT ASPT Cross-Validation & Diagnostics",
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
      title: "Factor-LMC Cross-Correlation Matrices & Dual CV Gains",
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
      title: "Multi-Source RSN+SSN Diagnostics & Measurement Error CV",
      note: "Analytic (1/w - 1)*tau^2 penalty weighting preserves surveillance RSN primacy while extracting headwater signal from SSN."
    },
    desc_dist: {
      title: "Nitrate Network Distributions & Sampling Bias Diagnostics",
      note: "Shows why monitored RSN mean (4.09 mg/L) exceeds predicted full-network mean (2.61 mg/L, p = 5.1e-45)."
    },
    desc_grad: {
      title: "Stratified Landscape, Aquifer & Strahler Order Gradients",
      note: "Stratified breakdown comparing upland vs lowland, Principal Aquifers, Strahler stream orders, and major RBDs."
    },
    comp_cv: {
      title: "Nitrate vs. Orthophosphate Cross-Validation Benchmark",
      note: "Side-by-side R2 & RMSE comparison across England, Humber, Thames, and Severn."
    },
    comp_forest: {
      title: "Covariate Effect Forest Plot: Nitrate vs. Orthophosphate",
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
      note: "Spatial shift across 95,773 England reaches; co-located surveillance reaches rise (+0.6252 median abs shift) as WaterBlitz dilution is removed."
    },
    wb_diff_0180: {
      title: "Orthophosphate (0180) Co-Kriging Difference Map: Favoured minus Equal Sources",
      note: "Highlights threshold tier divergence (30.84% tier match on WaterBlitz-only reaches) when surveillance RSN+SSN anchors dominate."
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
      title: "Equal vs. Favoured RSN+SSN Multi-Stratum Summary Comparison",
      note: "Compares median raw predictions, SE inflation/reduction across strata, compliance tier agreement, and latent cross-correlation."
    },
    wb_eq_xcorr: {
      title: "Equal-Sources Latent Cross-Correlation Matrix (0117 vs. 0180)",
      note: "All-England Corr_TU = +0.5214, Corr_TD = +0.5230, Corr_Nugget = +0.1335, Corr_Total = +0.2989."
    },
    wb_wt_xcorr: {
      title: "Favoured RSN+SSN Latent Cross-Correlation Matrix (0117 vs. 0180)",
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
      title: "6-Panel Transform Diagnostics Summary",
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
    { cat: "major", badge: "CS1", title: "Nitrate-N (0117) & Non-Linear GAM, Spline, RF & GBM Trends", desc: "Replace linear X*beta with --trend=gam, spline, rf, or gbm plus orthogonal OOB/OOF spatial residual kriging.", tab: "cs1" },
    { cat: "major", badge: "CS2", title: "Data Transformation Comparison (Raw, Log, Sqrt, Box-Cox)", desc: "Benchmarks 4 Nitrate transformations across physical positivity, Jensen median bias, residual symmetry (0.118 skew), and Profile-MLE Box-Cox (lambda=0.50).", tab: "cs7" },
    { cat: "major", badge: "CS3", title: "Macroinvertebrate Taxon Richness (WHPT N-Taxa) & GAM", desc: "Models 3,408 England inverts samples with --reverse-color, --trend=gam, and Robust Z-Score anomaly mapping (CV R2 0.1371 -> 0.2664).", tab: "cs9" },
    { cat: "major", badge: "CS4", title: "Ecological Health (WHPT ASPT) & Inverted Scales", desc: "Supports inverts (whpt_aspt), diatoms, macrophytes, rhs & lcms with --reverse-color palettes (CV R2 = 0.6402).", tab: "cs2" },
    { cat: "major", badge: "CS5", title: "Factor-LMC Co-Kriging & Multi-Source Fusion", desc: "Models P >= 2 heterotopic targets via --method=lmc|colocated and combines EA RSN + headwater SSN networks.", tab: "cs4" },
    { cat: "major", badge: "CS6", title: "Citizen Science WaterBlitz Co-Kriging Fusion", desc: "Ingest 44,982 England WaterBlitz citizen observations via format_custom_data.py and compare Equal vs. Surveillance-Favoured (rsn,ssn priority + 0.25 weight) bivariate Co-Kriging.", tab: "cs6" },
    { cat: "major", badge: "CS7", title: "Multi-Approach Eutrophication Synthesis (P, N, RMNI & DO)", desc: "Compares Factor-LMC Co-Kriging (0180+0117+rmni), GAM+S3N (R2=0.7607), and Seasonal DO + P-RMNI concordance (r=0.9005).", tab: "cs10" },
    { cat: "major", badge: "CS8", title: "Spatio-Temporal Batching & Dynamic Borrowing", desc: "--by=year|season|month + --st-borrow precision weighting, reach-level trend tests & period transitions.", tab: "cs3" },
    { cat: "major", badge: "CS9", title: "5-Tier Threshold Exceedance Risk Mapping", desc: "--threshold=uk applies indicative UK thresholds (11.3 mg/L NO3-N, 0.069 mg/L P) to map 5-tier Gaussian exceedance probability.", tab: "cs15" },
    { cat: "major", badge: "CS10", title: "Descriptive Stats, Monitoring Bias & N:P Stoichiometry", desc: "Scripts 05 & 04 run Mann-Whitney representativeness tests, gradient breakdowns, Fisher z & Wald Z-tests.", tab: "cs5" },
    { cat: "major", badge: "CS11", title: "Station-to-Reach Snapping Audits (Script 06)", desc: "Generates snapping distance summaries, unclipped maps & off-network outlier tables across ORN v2.", tab: "cs16" },
    { cat: "major", badge: "CS12", title: "Optimal Monitoring Network Design & Gap Mapping", desc: "--design-sites=<K> runs OpenMP C++ Schur-complement active learning and renders multi-panel gap maps (80.74% variance drop at K=50).", tab: "cs8" },
    { cat: "major", badge: "CS13", title: "North West Nitrate, Exceedance & Seasonal Trends (RBD=12)", desc: "Regional North West RBD (12,465 reaches) Nitrate predictions, 11.3 mg/L UK exceedance tiers (99.4% compliant), and 4-season borrowing.", tab: "cs11" },
    { cat: "major", badge: "CS14", title: "North West Inverts N-Taxa (RBD=12)", desc: "Regional North West WHPT N-Taxa richness (131 obs reaches, CV R2 0.1252 -> 0.1950, RMSE = 5.1640 taxa).", tab: "cs12" },
    { cat: "major", badge: "CS15", title: "North West Station-to-Reach Snapping Audits", desc: "Snapping precision QC across Ribble, Mersey, Weaver, Lune, and Derwent catchments (median 1.88-3.48 m).", tab: "cs13" },
    { cat: "major", badge: "CS16", title: "North West Optimal Monitoring Design (K=50)", desc: "OpenMP C++ Schur-complement active learning in North West RBD cuts variance by 45.95% (Mean SE 0.6956 -> 0.4887).", tab: "cs14" },
    { cat: "minor", badge: "Kernels", title: "5 Stream Covariance Kernels (--kernel)", desc: "exponential, spherical, mariah, linear, or auto (selects best kernel by Profile MLE likelihood).", tab: "cs1" },
    { cat: "minor", badge: "Confluence", title: "3 Physical AFV Weighting Schemes (--afv)", desc: "Branch weights by upstream length, Shreve link magnitude, or hydro_yield (orographic + aquifer proxy).", tab: "cs1" },
    { cat: "minor", badge: "Covariates", title: "15 Covariates & Mixed +/- Selection (--vars)", desc: "Customize predictors with --vars=+sinuosity,-slope; inspect full catalog via --list-vars.", tab: "cs1" },
    { cat: "minor", badge: "Transforms", title: "12 Transforms & 4 Normalisations", desc: "boxcox, yeojohnson, log, sqrt, logit, asin + Taylor/Duan back-transform; robust_zscore, zscore, minmax, percentile.", tab: "cs7" },
    { cat: "minor", badge: "Multi-Source", title: "Network Preference & Reliability Weighting", desc: "Combine EA surveillance RSN + headwater SSN networks using --prefer-source, --source-weight & --source-se.", tab: "cs4" },
    { cat: "minor", badge: "Censoring", title: "ROS, Site-ROS & Spatial Tobit EM (< LOD)", desc: "Handle left-censored '<' detection limits (--censor=ros|ros_site|em) and colorimetric ranges (--range-mode=tobit).", tab: "cs15" },
    { cat: "minor", badge: "Uncertainty", title: "Two-Stage Approach A Meas-Error Kriging", desc: "--meas-error / --gam and --source-se propagate Delta-Method observation SE into prediction SE maps.", tab: "cs15" },
    { cat: "minor", badge: "Validation", title: "Hydrological CV Modes (--cv-mode)", desc: "--cv-mode=random|catchment|strahler|spatial_block evaluates geographic & hydrological extrapolation.", tab: "cs1" }
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

  const renderedTabs = new Set();

  function ensureTabRendered(tabId) {
    if (!tabId || renderedTabs.has(tabId)) return;
    renderedTabs.add(tabId);
    const panel = document.getElementById("tab-" + tabId);
    if (!panel) return;
    renderTabDefaultFigures(tabId);
    attachKpiExplainers(panel);
    initWorkbenchLinks(panel);
    initCaseStudyPagerForTab(tabId);
  }

  function switchTab(tabId, updateHash) {
    if (tabId === "workbench") {
      window.location.href = "../cli_workbench/index.html";
      return;
    }
    ensureTabRendered(tabId);
    const targetBtn = document.querySelector(`.tab-btn[data-tab="${tabId}"]`);
    const activeGroup = targetBtn ? targetBtn.dataset.group : (tabId === "how-it-works" ? "how-it-works" : "overview");
    const prevGroupBtn = document.querySelector(".group-btn.active");
    const prevGroup = prevGroupBtn ? prevGroupBtn.dataset.group : "overview";
    const groupChanged = activeGroup !== prevGroup;

    document.querySelectorAll(".group-btn").forEach((btn) => {
      const isGrp = btn.dataset.group === activeGroup;
      btn.classList.toggle("active", isGrp);
      btn.setAttribute("aria-selected", isGrp ? "true" : "false");
    });

    const subtabStrip = document.getElementById("subtabStrip");
    if (subtabStrip) {
      subtabStrip.classList.toggle("hidden", activeGroup === "overview" || activeGroup === "how-it-works");
    }

    document.querySelectorAll(".tab-btn").forEach((btn) => {
      const isActiveTab = btn.dataset.tab === tabId;
      btn.classList.toggle("hidden-in-group", btn.dataset.group !== activeGroup);
      btn.classList.toggle("active", isActiveTab);
      btn.setAttribute("aria-selected", isActiveTab ? "true" : "false");
    });
    document.querySelectorAll(".tab-panel").forEach((p) => {
      p.classList.toggle("active", p.id === "tab-" + tabId);
    });
    if (updateHash !== false && window.location.hash !== "#" + tabId) {
      if (groupChanged) {
        history.pushState({ tab: tabId }, "", "#" + tabId);
      } else {
        history.replaceState({ tab: tabId }, "", "#" + tabId);
      }
    }
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  const CS5_NP_EXPLAINERS = {
    comp_cv: "<strong>What N vs P CV Performance Shows:</strong> Compares 5-fold out-of-sample cross-validation accuracy (<code>R&sup2;</code> and <code>RMSE</code>) across OLS and S3N models for Nitrate-N (<code>0117</code>) vs. Orthophosphate-P (<code>0180</code>), tested formally via Fisher's r-to-z test.",
    comp_forest: "<strong>What the Covariate Forest Plot Shows:</strong> Displays standardised fixed-effect regression coefficients (<code>&beta;</code> with 95% bootstrap CIs) and spatial variance partitioning (<code>&sigma;&sup2;</code> vs. <code>&tau;&sup2;</code>) side-by-side for Nitrate and Orthophosphate, with two-sample Wald Z-test significance flags.",
    comp_maps: "<strong>What the Side-by-Side Maps Show:</strong> Places the full-network England predictions (<code>90,876</code> reaches) for Nitrate-N and Orthophosphate-P next to each other to highlight shared agricultural/urban nutrient hotspots versus distinct geological or point-source corridors.",
    comp_rel: "<strong>What N:P Stoichiometry &amp; Quadrant Show:</strong> Maps reach-by-reach 4-quadrant co-occurrence (<code>High-N/High-P</code>, <code>Low-N/Low-P</code>, and discordant regimes) alongside molar <code>N:P</code> ratio distributions across basins, showing where reaches are P-limited (<code>&gt; 20:1</code>), Co-limited (<code>10:1&ndash;20:1</code>), or N-limited (<code>&lt; 10:1</code>)."
  };

  function renderFigureInto(slotId, figKey) {
    const slot = document.getElementById(slotId);
    if (!slot) return;
    if (slotId === "fig-slot-cs5-np") {
      const expEl = document.getElementById("fig-explainer-cs5-np");
      if (expEl) expEl.innerHTML = CS5_NP_EXPLAINERS[figKey] || "";
    }
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
      </div>
    `;
  }

  const DEFAULT_FIGURES_BY_TAB = {
    cs1: [["fig-slot-cs1-map", "no3_england"], ["fig-slot-cs1-diag", "no3_diag"]],
    cs2: [["fig-slot-cs2-ecol", "whpt_england"], ["fig-slot-cs2-diag", "whpt_diag"]],
    cs15: [["fig-slot-cs15-risk", "no3_exceed_england"]],
    cs3: [["fig-slot-cs3-grid", "st_season_grid"], ["fig-slot-cs3-trend", "st_year_trends"]],
    cs4: [["fig-slot-cs4-ck", "cokrig_xcorr"], ["fig-slot-cs4-ms", "rss_england"]],
    cs5: [["fig-slot-cs5-bias", "desc_grad"], ["fig-slot-cs5-np", "comp_cv"]],
    cs16: [["fig-slot-cs16-snap", "snap_overview"]],
    cs6: [["fig-slot-cs6-maps", "wb_diff_0117"], ["fig-slot-cs6-charts", "wb_comp_svg"]],
    cs7: [["fig-slot-cs7-maps", "tr_diff_map"], ["fig-slot-cs7-charts", "tr_comp_svg"]],
    cs8: [["fig-slot-cs8-design", "opt_design_england"]],
    cs9: [["fig-slot-cs9-map", "ntaxa_england"], ["fig-slot-cs9-diag", "ntaxa_diag"]],
    cs10: [["fig-slot-cs10-map", "eut_cokrig_rmni"], ["fig-slot-cs10-xcorr", "eut_cokrig_xcorr"], ["fig-slot-cs10-do", "eut_do_season_grid"], ["fig-slot-cs10-diag", "eut_rmni_gam_diag"]],
    cs11: [["fig-slot-cs11-map", "nw_no3_map"]],
    cs12: [["fig-slot-cs12-map", "nw_ntaxa_map"]],
    cs13: [["fig-slot-cs13-map", "nw_snap_ribble"]],
    cs14: [["fig-slot-cs14-design", "nw_opt_design"]]
  };

  function renderTabDefaultFigures(tabId) {
    const pairs = DEFAULT_FIGURES_BY_TAB[tabId];
    if (!pairs) return;
    pairs.forEach(([slotId, figKey]) => renderFigureInto(slotId, figKey));
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

    // --- CS9: Macroinvertebrate Taxon Richness (WHPT N-Taxa) + GAM ---
    const cs9KpiEl = document.getElementById("cs9-kpis");
    if (cs9KpiEl) {
      cs9KpiEl.innerHTML = [
        { val: "0.2648", lbl: "All-England Joint S3N R²", sub: "Up from 0.1371 OLS (+0.1277)" },
        { val: "0.2664", lbl: "GAM + S3N Joint R²", sub: "RMSE = 5.3926 taxa (vs 5.8486 OLS)" },
        { val: "0.1750", lbl: "Severn Basin Joint R²", sub: "Up from 0.0200 OLS (0.1790 GAM)" },
        { val: "3,408", lbl: "Inverts Samples (1,313 Sites)", sub: "1,160 England RSN observed reaches" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs9SvgEl = document.getElementById("cs9-r2-svg");
    if (cs9SvgEl) {
      cs9SvgEl.innerHTML = buildInlineBarSvg("All-England Held-Out CV R² (WHPT N-Taxa, Untransformed)", [
        { label: "OLS (Non-Spatial Linear)", val: 0.1371, color: "#64748b" },
        { label: "S3N Tail-Up (FC + AFV)", val: 0.1834, color: "#0284c7" },
        { label: "S3N Tail-Down (Stream)", val: 0.2391, color: "#38bdf8" },
        { label: "S3N Combined (TU+TD+Euc)", val: 0.2648, color: "#34d399" },
        { label: "GAM + S3N Combined", val: 0.2664, color: "#a78bfa" }
      ], 0.35);
    }

    const cvNtaxaAll = [
      ...(t.cv_whpt_ntaxa || []).filter((r) => r.Domain === "All England RSN Network" || (r.Model || "").includes("Combined")),
      ...(t.cv_whpt_ntaxa_gam || []).filter((r) => r.Domain === "All England RSN Network" && (r.Model || "").includes("Combined"))
    ];
    renderTable("cs9-cv-table", [
      { key: "Domain", label: "Domain" },
      { key: "Model", label: "Model" },
      { key: "Spatial_Method", label: "Trend" },
      { key: "n_obs", label: "n_obs", num: true },
      { key: "RMSE", label: "RMSE", num: true, digits: 4 },
      { key: "R2", label: "CV R²", num: true, digits: 4, clsFn: (r) => (r.R2 >= 0.25 ? "highlight-cell" : "") }
    ], cvNtaxaAll);

    const pNtaxa = (t.params_whpt_ntaxa || []).filter((r) => r.Domain === "All England RSN Network" && (r.Model || "").includes("Joint"));
    renderTable("cs9-param-table", [
      { key: "Parameter", label: "Parameter" },
      { key: "Estimate", label: "Estimate", num: true, digits: 4 }
    ], pNtaxa);

    // --- CS10: Eutrophication Modelling Synthesis ---
    const cs10KpiEl = document.getElementById("cs10-kpis");
    if (cs10KpiEl) {
      cs10KpiEl.innerHTML = [
        { val: "0.7607", lbl: "Best National RMNI R² (GAM+S3N)", sub: "RMSE = 0.6846 (vs 0.6289 OLS)" },
        { val: "0.7580", lbl: "LMC Co-Kriging RMNI R² (Aux)", sub: "Thames rescued -0.2671 -> +0.3165" },
        { val: "0.9005", lbl: "P (0180) vs RMNI Reach r", sub: "86.55% 4-quadrant concordance (OR=41.4)" },
        { val: "0.1574", lbl: "Summer DO (9901) Joint R²", sub: "Diel noise limits grab-sample DO CV" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs10SvgEl = document.getElementById("cs10-r2-svg");
    if (cs10SvgEl) {
      cs10SvgEl.innerHTML = buildInlineBarSvg("All-England Held-Out CV R² Across Eutrophication Models", [
        { label: "Summer DO % Sat (9901 S3N)", val: 0.1574, color: "#f97316" },
        { label: "RMNI Linear OLS", val: 0.6289, color: "#64748b" },
        { label: "RMNI Non-Spatial GAM", val: 0.7319, color: "#38bdf8" },
        { label: "RMNI Linear Joint S3N", val: 0.7398, color: "#0284c7" },
        { label: "RMNI LMC Co-Kriging (Aux)", val: 0.7580, color: "#34d399" },
        { label: "RMNI GAM + Joint S3N", val: 0.7607, color: "#a78bfa" }
      ], 0.85);
    }

    const eutCvRows = [
      ...(t.cv_eut_cokrig || []).filter((r) => r.Target_Code === "rmni" && (r.Domain === "All England RSN Network" || r.Domain.includes("Thames") || r.Domain.includes("Severn"))).map((r) => ({
        Approach: "1. LMC Co-Kriging",
        Domain: r.Domain,
        Model: r.Model,
        RMSE: r.RMSE,
        R2: r.R2
      })),
      ...(t.cv_eut_rmni || []).filter((r) => r.Domain === "All England RSN Network" && (r.Model.includes("OLS") || r.Model.includes("Combined"))).map((r) => ({
        Approach: "2. Univariate Linear",
        Domain: r.Domain,
        Model: r.Model,
        RMSE: r.RMSE,
        R2: r.R2
      })),
      ...(t.cv_eut_rmni_gam || []).filter((r) => (r.Domain === "All England RSN Network" || r.Domain.includes("Thames")) && r.Model.includes("Combined")).map((r) => ({
        Approach: "2. Univariate GAM+S3N",
        Domain: r.Domain,
        Model: r.Model + " [GAM]",
        RMSE: r.RMSE,
        R2: r.R2
      })),
      ...(t.cv_eut_do_summer || []).filter((r) => r.Domain === "All England RSN Network" && r.Model.includes("Combined")).map((r) => ({
        Approach: "3. Summer DO (9901)",
        Domain: r.Domain,
        Model: r.Model + " (log %)",
        RMSE: r.RMSE,
        R2: r.R2
      }))
    ];
    renderTable("cs10-cv-table", [
      { key: "Approach", label: "Approach" },
      { key: "Domain", label: "Domain" },
      { key: "Model", label: "Model" },
      { key: "RMSE", label: "RMSE", num: true, digits: 4 },
      { key: "R2", label: "CV R²", num: true, digits: 4, clsFn: (r) => (r.R2 >= 0.73 ? "highlight-cell" : "") }
    ], eutCvRows);

    renderTable("cs10-rel-table", [
      { key: "Stratum", label: "Stratum (0180 vs RMNI)" },
      { key: "n_reaches", label: "Reaches", num: true },
      { key: "Pearson_r", label: "Pearson r", num: true, digits: 4 },
      { key: "Pct_Q1_HighHigh", label: "% High-High", num: true, digits: 2 },
      { key: "Pct_Concordant", label: "% Concordant", num: true, digits: 2, clsFn: (r) => (r.Pct_Concordant >= 85 ? "highlight-cell" : "") }
    ], t.eut_p_rmni_rel || []);

    // --- CS11: North West Nitrate, Exceedance & Seasonal Trends ---
    const cs11KpiEl = document.getElementById("cs11-kpis");
    if (cs11KpiEl) {
      cs11KpiEl.innerHTML = [
        { val: "0.6702", lbl: "North West Joint S3N R²", sub: "12,465 reaches, 136 obs (vs 0.6523 OLS)" },
        { val: "93.8%", lbl: "Confidently Compliant", sub: "11,686 reaches < 11.3 mg/L (0% exceeding)" },
        { val: "0.6799", lbl: "Peak Autumn Seasonal R²", sub: "Summer 0.6704 | Spring 0.6576 | Winter 0.6197" },
        { val: "0.8429", lbl: "Joint S3N RMSE (log mg/L)", sub: "Temporal borrowing rho_t = 0.363" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs11SvgEl = document.getElementById("cs11-r2-svg");
    if (cs11SvgEl) {
      cs11SvgEl.innerHTML = buildInlineBarSvg("North West RBD (RBD=12) Nitrate Held-Out CV R²", [
        { label: "Annual OLS (Non-Spatial)", val: 0.6523, color: "#64748b" },
        { label: "Annual Joint S3N", val: 0.6702, color: "#34d399" },
        { label: "Winter (DJF) Joint S3N", val: 0.6197, color: "#38bdf8" },
        { label: "Spring (MAM) Joint S3N", val: 0.6576, color: "#0284c7" },
        { label: "Summer (JJA) Joint S3N", val: 0.6704, color: "#f97316" },
        { label: "Autumn (SON) Joint S3N", val: 0.6799, color: "#a78bfa" }
      ], 0.80);
    }

    const nw0117Rows = [
      ...(t.cv_nw_0117 || []).filter((r) => (r.Domain || "").includes("North West")),
      ...(t.cv_nw_0117_winter || []).filter((r) => (r.Domain || "").includes("North West") && (r.Model || "").includes("Combined")),
      ...(t.cv_nw_0117_spring || []).filter((r) => (r.Domain || "").includes("North West") && (r.Model || "").includes("Combined")),
      ...(t.cv_nw_0117_summer || []).filter((r) => (r.Domain || "").includes("North West") && (r.Model || "").includes("Combined")),
      ...(t.cv_nw_0117_autumn || []).filter((r) => (r.Domain || "").includes("North West") && (r.Model || "").includes("Combined"))
    ];
    renderTable("cs11-cv-table", [
      { key: "Period", label: "Period" },
      { key: "Model", label: "Model" },
      { key: "n_obs", label: "n_obs", num: true },
      { key: "RMSE", label: "RMSE", num: true, digits: 4 },
      { key: "R2", label: "CV R²", num: true, digits: 4, clsFn: (r) => (r.R2 >= 0.67 ? "highlight-cell" : "") }
    ], nw0117Rows);

    renderTable("cs11-exceed-table", [
      { key: "Tier", label: "Exceedance Tier (11.3 mg/L)" },
      { key: "Reaches", label: "Reaches", num: true },
      { key: "Pct", label: "% of NW", num: true, digits: 1, clsFn: (r) => (r.Tier === "Confidently Compliant" ? "highlight-cell" : "") }
    ], [
      { Tier: "Confidently Compliant", Reaches: 11686, Pct: 93.8 },
      { Tier: "Likely Compliant", Reaches: 704, Pct: 5.7 },
      { Tier: "Uncertain / Priority Monitoring", Reaches: 75, Pct: 0.6 },
      { Tier: "Likely Exceeding", Reaches: 0, Pct: 0.0 },
      { Tier: "Confidently Exceeding", Reaches: 0, Pct: 0.0 }
    ]);

    // --- CS12: North West Inverts N-Taxa ---
    const cs12KpiEl = document.getElementById("cs12-kpis");
    if (cs12KpiEl) {
      cs12KpiEl.innerHTML = [
        { val: "0.1950", lbl: "North West N-Taxa Joint R²", sub: "Up from 0.1252 OLS (+0.0698)" },
        { val: "5.1640", lbl: "Joint S3N RMSE (taxa)", sub: "131 observed reaches across 12,465" },
        { val: "13.92", lbl: "Euclidean Variance (sigma2_Euc)", sub: "Dominant spatial gradient in NW RBD" },
        { val: "7.56", lbl: "Tail-Up Variance (sigma2_TU)", sub: "Flow-connected AFV continuity (tau2=4.30)" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs12SvgEl = document.getElementById("cs12-r2-svg");
    if (cs12SvgEl) {
      cs12SvgEl.innerHTML = buildInlineBarSvg("North West RBD (RBD=12) WHPT N-Taxa Held-Out CV R²", [
        { label: "OLS (Non-Spatial)", val: 0.1252, color: "#64748b" },
        { label: "S3N Tail-Down (Stream)", val: 0.1090, color: "#38bdf8" },
        { label: "S3N Tail-Up (FC + AFV)", val: 0.1229, color: "#0284c7" },
        { label: "S3N Combined (TU+TD+Euc)", val: 0.1950, color: "#34d399" }
      ], 0.25);
    }

    const nwNtaxaRows = (t.cv_nw_ntaxa || []).filter((r) => (r.Domain || "").includes("North West"));
    renderTable("cs12-cv-table", [
      { key: "Domain", label: "Domain" },
      { key: "Model", label: "Model" },
      { key: "n_obs", label: "n_obs", num: true },
      { key: "RMSE", label: "RMSE", num: true, digits: 4 },
      { key: "R2", label: "CV R²", num: true, digits: 4, clsFn: (r) => (r.R2 >= 0.19 ? "highlight-cell" : "") }
    ], nwNtaxaRows);

    // --- CS13: North West Station Snapping Diagnostics ---
    const nwSnapIds = new Set(["161", "177", "173", "157", "138"]);
    const nwSnapRows = (t.nw_snapping || []).filter((r) => nwSnapIds.has(String(r.CatchID)) || (r.RBD_Name || "") === "North West");

    const cs13KpiEl = document.getElementById("cs13-kpis");
    if (cs13KpiEl) {
      cs13KpiEl.innerHTML = [
        { val: "1.88 m", lbl: "Ribble Median Snap Offset", sub: "31 included, 2 excluded (87.1% <= 10 m)" },
        { val: "100.0%", lbl: "Weaver Snapping <= 10 m", sub: "23 included, median 2.18 m" },
        { val: "95.7%", lbl: "Mersey Snapping <= 10 m", sub: "23 included, 0 excluded, median 2.15 m" },
        { val: "76.5%", lbl: "Lune Snapping <= 10 m", sub: "17 included, 0 excluded, median 3.48 m" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs13SvgEl = document.getElementById("cs13-snap-svg");
    if (cs13SvgEl && nwSnapRows.length) {
      cs13SvgEl.innerHTML = buildInlineBarSvg("North West Catchments: % of Included Stations Snapped <= 10 m", nwSnapRows.slice(0, 6).map((r) => ({
        label: `${r.Main_River} (${r.CatchID})`,
        val: Number(r.Pct_Within_10m || 0),
        color: Number(r.Pct_Within_10m) >= 90 ? "#34d399" : "#38bdf8"
      })), 100, "%");
    }

    renderTable("cs13-snap-table", [
      { key: "Main_River", label: "Catchment" },
      { key: "Total_Reaches", label: "Reaches", num: true },
      { key: "Stations_Included", label: "Included", num: true },
      { key: "Stations_Excluded", label: "Excl.", num: true },
      { key: "Median_Snap_m", label: "Median (m)", num: true, digits: 2 },
      { key: "Pct_Within_10m", label: "<= 10 m (%)", num: true, digits: 1, clsFn: (r) => (r.Pct_Within_10m >= 85 ? "highlight-cell" : "") }
    ], nwSnapRows);

    // --- CS14: North West Optimal Monitoring Design ---
    const nwDesignRows = (t.nw_design || []).filter((r) => (r.Domain || "").includes("North West"));
    const cs14KpiEl = document.getElementById("cs14-kpis");
    if (cs14KpiEl) {
      cs14KpiEl.innerHTML = [
        { val: "45.95%", lbl: "NW Variance Reduction (K=50)", sub: "Ranked across 12,465 reaches in 21.33s" },
        { val: "0.4887", lbl: "Mean Kriging SE After K=50", sub: "Down -29.7% from 0.6956 log mg/L" },
        { val: "136", lbl: "Existing NW RSN Stations", sub: "93.8% already confidently compliant" },
        { val: "+2.77%", lbl: "#1 Site Variance Drop", sub: "Riv 38848 (Ribble, Order 6)" }
      ].map((k) => `
        <div class="kpi-card">
          <div class="kpi-label">${esc(k.lbl)}</div>
          <div class="kpi-val">${k.val}</div>
          <div class="kpi-note">${esc(k.sub)}</div>
        </div>
      `).join("");
    }

    const cs14SvgEl = document.getElementById("cs14-design-svg");
    if (cs14SvgEl && nwDesignRows.length) {
      const milestones = [1, 5, 10, 25, 50].map((k) => nwDesignRows.find((r) => Number(r.Rank) === k)).filter(Boolean);
      cs14SvgEl.innerHTML = buildInlineBarSvg("North West Cumulative Variance Reduction by Added Stations (K)", milestones.map((r) => ({
        label: `Top ${r.Rank} Station${r.Rank === 1 ? "" : "s"} (SE ${Number(r.Mean_SE_After).toFixed(3)})`,
        val: Number(r.Cum_Pct_Reduced || 0),
        color: Number(r.Rank) === 50 ? "#34d399" : "#38bdf8"
      })), 50, "%");
    }

    renderTable("cs14-design-table", [
      { key: "Rank", label: "Rank", num: true },
      { key: "RivID", label: "RivID", num: true },
      { key: "CatchID", label: "CatchID", num: true },
      { key: "Strahler", label: "Order", num: true },
      { key: "Elevation", label: "Elev (m)", num: true, digits: 1 },
      { key: "Current_SE_sc", label: "Local SE", num: true, digits: 4 },
      { key: "Pct_Var_Reduced", label: "dV (%)", num: true, digits: 2 },
      { key: "Cum_Pct_Reduced", label: "Cum. dV (%)", num: true, digits: 2, clsFn: (r) => (r.Rank <= 5 ? "highlight-cell" : "") }
    ], nwDesignRows.slice(0, 20));

    // --- CS2: WHPT ASPT Ecological Health ---
    const engWhpt = (t.cv_whpt || []).filter((r) => r.Domain === "All England RSN Network");
    const whptComb = engWhpt.find((r) => String(r.Model).includes("Combined")) || {};
    const whptOls = engWhpt.find((r) => String(r.Model).includes("OLS")) || {};
    const kpi2 = document.getElementById("cs2-kpis");
    if (kpi2) {
      kpi2.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Joint S3N CV R²</div><div class="kpi-val pos">${fmt(whptComb.R2, 4)}</div><div class="kpi-note">RMSE ${fmt(whptComb.RMSE, 4)} (vs ${fmt(whptOls.RMSE, 4)} OLS)</div></div>
        <div class="kpi-card"><div class="kpi-label">Spatial Variance Gain</div><div class="kpi-val pos">+0.1470</div><div class="kpi-note">R² lift from Tail-Up + Tail-Down + Euclidean</div></div>
        <div class="kpi-card"><div class="kpi-label">Humber Basin CV R²</div><div class="kpi-val pos">0.5704</div><div class="kpi-note">Up from 0.5028 OLS (Severn 0.3416, Thames 0.3051)</div></div>
        <div class="kpi-card"><div class="kpi-label">Observed Sites</div><div class="kpi-val">1,160 reaches</div><div class="kpi-note">Predicted out-of-sample across 90,876 reaches</div></div>
      `;
    }
    const whptSvg = document.getElementById("cs2-r2-svg");
    if (whptSvg && engWhpt.length) {
      const bars = engWhpt.map((r) => ({
        label: r.Model.replace(" (Flow-Connected + AFV)", "").replace(" (Tail-Up + Tail-Down)", "").replace(" (Stream Distance)", ""),
        val: r.R2 || 0,
        color: r.Model.includes("Combined") ? "#34d399" : (r.Model.includes("OLS") ? "#64748b" : "#38bdf8")
      }));
      whptSvg.innerHTML = buildInlineBarSvg("ALL ENGLAND WHPT ASPT 5-FOLD CV R2 COMPARISON", bars, 0.75, "");
    }

    renderTable("cs2-cv-table",
      [
        { label: "All England Model", key: "Model" },
        { label: "n_obs", key: "n_obs", num: true },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "MAE", key: "MAE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      engWhpt,
      (r) => String(r.Model).includes("Combined") ? "highlight-row" : ""
    );

    renderTable("cs2-basin-table",
      [
        { label: "Regional Basin", key: "Domain" },
        { label: "Model", key: "Model" },
        { label: "n_obs", key: "n_obs", num: true },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      (t.cv_whpt || []).filter((r) => r.Domain !== "All England RSN Network" && (String(r.Model).includes("OLS") || String(r.Model).includes("Combined"))),
      (r) => String(r.Model).includes("Combined") ? "highlight-row" : ""
    );

    // --- CS15: Threshold Exceedance Risk & Anomalies ---
    const ex0117 = (t.exceed_0117 || []).filter((r) => r.Domain === "All England RSN Network");
    const ex0180 = (t.exceed_0180 || []).filter((r) => r.Domain === "All England RSN Network");
    const eng0180 = (t.cv_0180 || []).filter((r) => r.Domain === "All England RSN Network");
    const pComb = eng0180.find((r) => String(r.Model).includes("Combined")) || {};
    const pOls = eng0180.find((r) => String(r.Model).includes("OLS")) || {};
    const kpi15 = document.getElementById("cs15-kpis");
    if (kpi15) {
      kpi15.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Nitrate Compliant (&le; 11.3)</div><div class="kpi-val pos">90.8%</div><div class="kpi-note">63.6% Confident + 27.3% Likely (0.08% exceeding)</div></div>
        <div class="kpi-card"><div class="kpi-label">PO4 Uncertain Tier</div><div class="kpi-val warn">47.5%</div><div class="kpi-note">43,202 reaches with P(exceed 0.069) in [0.20, 0.80]</div></div>
        <div class="kpi-card"><div class="kpi-label">Orthophosphate CV R²</div><div class="kpi-val pos">${fmt(pComb.R2, 4)}</div><div class="kpi-note">RMSE ${fmt(pComb.RMSE, 4)} vs ${fmt(pOls.RMSE, 4)} OLS (log)</div></div>
        <div class="kpi-card"><div class="kpi-label">National Reach Coverage</div><div class="kpi-val">90,876</div><div class="kpi-note">1,210 RSN WQ stations &middot; 5-tier Gaussian risk</div></div>
      `;
    }
    const exSvg = document.getElementById("cs15-exceed-svg");
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

    const r2Svg15 = document.getElementById("cs15-r2-svg");
    if (r2Svg15 && eng0180.length) {
      const bars = eng0180.map((r) => ({
        label: r.Model.replace(" (Flow-Connected + AFV)", "").replace(" (Tail-Up + Tail-Down)", "").replace(" (Stream Distance)", ""),
        val: r.R2 || 0,
        color: r.Model.includes("Combined") ? "#34d399" : (r.Model.includes("OLS") ? "#64748b" : "#38bdf8")
      }));
      r2Svg15.innerHTML = buildInlineBarSvg("ALL ENGLAND ORTHOPHOSPHATE (0180) 5-FOLD CV R2 COMPARISON", bars, 0.75, "");
    }

    renderTable("cs15-exceed-table",
      [
        { label: "Determinand", fn: (r) => r.Determinand_Name || (r.Threshold_Raw === 11.3 ? "Nitrate" : "Orthophosphate") },
        { label: "Tier", key: "Tier" },
        { label: "Reaches", key: "Reach_Count", num: true },
        { label: "Share %", key: "Reach_Pct", num: true, digits: 1 }
      ],
      [...ex0117, ...ex0180]
    );

    renderTable("cs15-cv-table",
      [
        { label: "Domain", key: "Domain" },
        { label: "Model", key: "Model" },
        { label: "n_obs", key: "n_obs", num: true },
        { label: "RMSE", key: "RMSE", num: true, digits: 4 },
        { label: "MAE", key: "MAE", num: true, digits: 4 },
        { label: "R²", key: "R2", num: true, digits: 4 }
      ],
      (t.cv_0180 || []).filter((r) => String(r.Model).includes("OLS") || String(r.Model).includes("Combined")),
      (r) => String(r.Model).includes("Combined") ? "highlight-row" : ""
    );

    // --- CS3: Spatio-Temporal ---
    const stWin = (t.st_season_cv || []).find((r) => String(r.Period).includes("winter")) || {};
    const stSum = (t.st_season_cv || []).find((r) => String(r.Period).includes("summer")) || {};
    const stYrR2 = (t.st_year_cv || []).map((r) => r.R2).filter((v) => typeof v === "number");
    const stYrMin = stYrR2.length ? Math.min(...stYrR2) : 0.5772;
    const stYrMax = stYrR2.length ? Math.max(...stYrR2) : 0.6218;
    const kpi3 = document.getElementById("cs3-kpis");
    if (kpi3) {
      kpi3.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Winter Peak CV R²</div><div class="kpi-val pos">${fmt(stWin.R2, 4)}</div><div class="kpi-note">RMSE ${fmt(stWin.RMSE, 4)} across ${fmt(stWin.n_obs, 0)} obs reaches</div></div>
        <div class="kpi-card"><div class="kpi-label">Summer Low-Flow CV R²</div><div class="kpi-val warn">${fmt(stSum.R2, 4)}</div><div class="kpi-note">RMSE ${fmt(stSum.RMSE, 4)} across ${fmt(stSum.n_obs, 0)} obs reaches</div></div>
        <div class="kpi-card"><div class="kpi-label">Annual CV R² (2021–2024)</div><div class="kpi-val pos">${fmt(stYrMin, 3)}&ndash;${fmt(stYrMax, 3)}</div><div class="kpi-note">434&ndash;439 obs/yr across 90,876 England reaches</div></div>
        <div class="kpi-card"><div class="kpi-label">Temporal Trend Workflow</div><div class="kpi-val">4 Seasons &middot; 4 Years</div><div class="kpi-note">Precision-borrowed WLS slopes (fig6) &amp; shifts (fig7)</div></div>
      `;
    }

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
    const ckEng0117 = (t.cokrig_cv || []).filter((r) => r.Domain === "All England RSN Network" && (r.Target_Code === "0117" || Number(r.Target_Code) === 117));
    const ckUni = ckEng0117.find((r) => String(r.Model).includes("Univariate")) || {};
    const ckAux = ckEng0117.find((r) => String(r.Model).includes("Aux-Aided")) || {};
    const ckEng0180 = (t.cokrig_cv || []).filter((r) => r.Domain === "All England RSN Network" && (r.Target_Code === "0180" || Number(r.Target_Code) === 180));
    const ck180Uni = ckEng0180.find((r) => String(r.Model).includes("Univariate")) || {};
    const ck180Aux = ckEng0180.find((r) => String(r.Model).includes("Aux-Aided")) || {};
    const ck180Drop = (ck180Uni.RMSE && ck180Aux.RMSE) ? (((1 - ck180Aux.RMSE / ck180Uni.RMSE) * 100).toFixed(1)) : "—";
    const rssComb = (t.cv_rss || []).find((r) => r.Domain === "All England RSS Combined Network" && String(r.Model).includes("Combined")) || {};

    const kpi4 = document.getElementById("cs4-kpis");
    if (kpi4) {
      kpi4.innerHTML = `
        <div class="kpi-card"><div class="kpi-label">Co-Kriging Verdict</div><div class="kpi-val pos">8 / 8 Strong</div><div class="kpi-note">Outperformed Univariate across all basins &amp; targets</div></div>
        <div class="kpi-card"><div class="kpi-label">England 0117 Aux-Aided</div><div class="kpi-val pos">${fmt(ckAux.R2, 4)}</div><div class="kpi-note">Up from ${fmt(ckUni.R2, 4)} Univariate (-22.3% RMSE)</div></div>
        <div class="kpi-card"><div class="kpi-label">England 0180 Aux-Aided</div><div class="kpi-val pos">${fmt(ck180Aux.R2, 4)}</div><div class="kpi-note">Up from ${fmt(ck180Uni.R2, 4)} Univariate (-${ck180Drop}% RMSE)</div></div>
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
        <circle cx="406.7" cy="131.1" r="4.5" fill="#38bdf8"/><text x="325" y="112" fill="#38bdf8" font-size="11" font-weight="700">SSN w=0.75 (0.33x)</text>
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

    renderTable("cs16-snap-table",
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
    const searchInput = document.getElementById("catalogSearch");
    if (!grid) return;

    let activeCat = "all";
    let searchTerm = "";

    function draw() {
      const q = searchTerm.trim().toLowerCase();
      const items = FEATURE_CATALOG.filter((f) => {
        if (activeCat !== "all" && f.cat !== activeCat) return false;
        if (!q) return true;
        return (f.title + " " + f.desc + " " + f.badge).toLowerCase().includes(q);
      });
      if (!items.length) {
        grid.innerHTML = `<div class="feature-item" style="grid-column:1/-1;"><p style="margin:0;">No catalog features match "<code>${esc(searchTerm)}</code>".</p></div>`;
        return;
      }
      grid.innerHTML = items.map((f) => `
        <div class="feature-item">
          <span class="pill ${f.cat === "major" ? "pill-green" : "pill-amber"}">${esc(f.badge)}</span>
          <h4>${esc(f.title)}</h4>
          <p>${esc(f.desc)}</p>
          <div style="margin-top:8px;">
            <button class="copy-btn" data-jump-tab="${esc(f.tab)}">${f.tab === "workbench" ? "Open CLI Workbench &rarr;" : "Open Case Study &rarr;"}</button>
          </div>
        </div>
      `).join("");
    }

    draw();
    btns.forEach((btn) => {
      btn.addEventListener("click", () => {
        btns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        activeCat = btn.dataset.fcat || "all";
        draw();
      });
    });
    if (searchInput) {
      searchInput.addEventListener("input", () => {
        searchTerm = searchInput.value || "";
        draw();
      });
    }
  }

  function initSearch() {
    const input = document.getElementById("featureSearch");
    const drop = document.getElementById("searchDropdown");
    if (!input || !drop) return;

    const index = [
      ...FEATURE_CATALOG.map((f) => ({ title: f.title, sub: f.desc, tab: f.tab })),
      ...Object.entries(FIGURE_META)
        .filter(([k]) => Boolean((DATA.figures || {})[k]))
        .map(([k, m]) => ({
        title: m.title,
        sub: m.note,
        tab: k === "no3_z_england" || k === "no3_z_thames_severn" || k.startsWith("cs1") || k.startsWith("no3_gam") || k === "no3_england" || k === "no3_humber" || k === "no3_thames_severn" || k === "no3_diag" ? "cs1"
          : k.startsWith("whpt") ? "cs2"
          : k.includes("exceed") || k.includes("no3_z") ? "cs15"
          : k.startsWith("st_") ? "cs3"
          : k.startsWith("cokrig") || k.startsWith("rss_") ? "cs4"
          : k.startsWith("wb_") ? "cs6"
          : k.startsWith("tr_") ? "cs7"
          : k.startsWith("opt_design_") ? "cs8"
          : k.startsWith("ntaxa_") ? "cs9"
          : k.startsWith("eut_") ? "cs10"
          : k.startsWith("nw_no3_") ? "cs11"
          : k.startsWith("snap_") ? "cs16"
          : k.startsWith("nw_snap_") ? "cs13"
          : k.startsWith("nw_opt_") ? "cs14"
          : k.startsWith("nw_") ? "cs12" : "cs5"
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
    const downloadBtn = document.getElementById("modalDownloadBtn");
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
      if (downloadBtn) {
        downloadBtn.href = fig.path;
        downloadBtn.download = (fig.path || "figure.png").split("/").pop();
      }
      img.classList.remove("zoomed");
      if (bodyEl) {
        bodyEl.scrollTop = 0;
        bodyEl.scrollLeft = 0;
      }
      title.textContent = info.title;
      if (gallery.length > 1) {
        metaPill.textContent = `${currentIdx + 1} / ${gallery.length}`;
        metaPill.hidden = false;
      } else {
        metaPill.textContent = "";
        metaPill.hidden = true;
      }
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
      btn.textContent = theme === "dark" ? "☀️" : "🌙";
      btn.title = theme === "dark" ? "Switch to Light Mode" : "Switch to Dark Mode";
      try { localStorage.setItem("modelling_dashboard_theme", theme); } catch (e) {}
    }
    let saved = "dark";
    try { saved = localStorage.getItem("modelling_dashboard_theme") || "dark"; } catch (e) {}
    apply(saved);
    btn.addEventListener("click", () => apply(root.getAttribute("data-theme") === "dark" ? "light" : "dark"));
  }

  // Boot
  document.querySelectorAll(".group-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const targetTab = btn.dataset.defaultTab || "overview";
      switchTab(targetTab);
    });
  });
  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => switchTab(btn.dataset.tab));
  });
  const STAGE_EXPLAINERS = {
    "1": {
      title: "Stage 1 · Python Stdlib Preprocessing (GIS, ROS & Snapping)",
      body: "Reads OS Open Rivers Network v2 (ORN v2) GeoPackage topology and binary WKB geometries directly using pure-stdlib Python SQLite—no GDAL or sf required. Builds the directed river graph, snaps monitoring stations to nearest reaches (flagging off-network outliers), and imputes left-censored below-LOD or range-limited values via Lognormal ROS or Tobit preprocessing."
    },
    "2": {
      title: "Stage 2 · C++ / OpenMP Euler-Tour Stream NNGP Engine",
      body: "Computes O(1) lowest-common-ancestor (LCA) hydrological distances via Euler-tour RMQ in compiled C++ with OpenMP parallelism. Constructs sparse Nearest-Neighbor Gaussian Process (NNGP) precision factors combining flow-connected Tail-Up (with AFV confluence weights), flow-unconnected Tail-Down, and Euclidean kernels alongside linear, GAM, RF, or GBM fixed-effect trends."
    },
    "3": {
      title: "Stage 3 · Stream Modelling (Univariate & Co-Kriging)",
      body: "Fits either single-determinand Univariate S3N models (water quality or ecological indices such as WHPT ASPT / N-Taxa / RMNI) or multivariate Factor-LMC / colocated Co-Kriging models across coupled targets (e.g., Nitrate + Orthophosphate + RMNI). Supports spatio-temporal period batching with precision-weighted borrowing and multi-source surveillance/citizen-science weighting."
    },
    "4": {
      title: "Stage 4 · Diagnostics, Thresholds & Quality Control",
      body: "Evaluates out-of-sample 5-fold cross-validation, renders national/basin prediction and kriging SE maps, and computes 5-tier probabilistic threshold exceedance risk. Also runs station-to-reach snapping distance audits, Mann-Whitney network representativeness bias tests, and cross-run N:P or P–RMNI concordance analyses."
    }
  };

  function selectPipelineStage(stageId) {
    const info = STAGE_EXPLAINERS[stageId];
    const box = document.getElementById("pipelineStageExplainer");
    if (!info || !box) return;
    document.querySelectorAll(".pipe-node[data-stage]").forEach((n) => {
      n.classList.toggle("active", n.dataset.stage === stageId);
    });
    box.dataset.stage = stageId;
    box.innerHTML = `<h4 style="margin-top:0;">${esc(info.title)}</h4><p style="margin:0; color:var(--text-secondary); font-size:0.86rem; line-height:1.55;">${esc(info.body)}</p>`;
  }

  document.querySelectorAll(".pipe-node[data-stage]").forEach((node) => {
    const activate = () => selectPipelineStage(node.dataset.stage);
    node.addEventListener("mouseenter", activate);
    node.addEventListener("focus", activate);
    node.addEventListener("click", activate);
  });
  selectPipelineStage("1");

  const KPI_EXPLAINERS = {
    "Linear S3N CV R²": "5-fold cross-validation R² for Nitrate (0117) using linear environmental covariates plus Joint Tail-Up, Tail-Down, and Euclidean NNGP stream covariance.",
    "GAM + S3N CV R²": "5-fold CV R² when linear fixed effects are replaced by REML-penalised cubic regression splines (mgcv::gam) before fitting spatial stream-network residuals.",
    "Joint S3N CV R²": "5-fold out-of-sample CV R² for WHPT ASPT ecological quality across England using Joint Tail-Up + Tail-Down + Euclidean NNGP stream covariance.",
    "Humber Basin CV R²": "5-fold CV R² for WHPT ASPT within the Humber Basin (CatchID=80), improving from 0.5028 (non-spatial OLS) to 0.5704.",
    "Nitrate Compliant (≤ 11.3)": "Share of England reaches (90.8%) classified as Confidently or Likely Compliant below the 11.3 mg/L N Nitrates Directive standard.",
    "PO4 Uncertain Tier": "Percentage of England reaches (47.5%, 43,202 reaches) with exceedance probability P(exceed 0.069 mg/L P) in [0.20, 0.80] needing priority monitoring.",
    "Orthophosphate CV R²": "5-fold CV R² for reactive Orthophosphate (0180) on the natural-log scale under Joint S3N.",
    "National Reach Coverage": "Total OS Open Rivers Network (ORN v2) reaches across 141 monitored England catchments covered by full-network S3N prediction and tier mapping.",
    "Spatial Variance Gain": "Additional out-of-sample variance explained (+0.2126 R²) purely from modelling river-network and Euclidean spatial correlation over non-spatial OLS.",
    "Observed Sites": "1,210 ORN v2 reaches with EA surveillance (RSN) Nitrate samples used to condition predictions across all 90,876 reaches in 141 monitored catchments.",

    "Co-Kriging Verdict": "Factor-LMC co-kriging outperformed single-target univariate S3N in all 8 Domain × Target cross-validation benchmarks (England, Humber, Thames, Severn × Nitrate, Orthophosphate).",
    "England 0117 Aux-Aided": "Cross-validation R² for Nitrate when Nitrate is held out at a reach while co-located Orthophosphate (0180) remains observed, cutting RMSE by 22.3%.",
    "England 0180 Aux-Aided": "Cross-validation R² for Orthophosphate (0180) when Phosphate is held out at a reach while co-located Nitrate (0117) remains observed, borrowing strength through Factor-LMC cross-correlation.",
    "Median Aux-Aided Gain": "Median improvement across all 8 Domain × Target pairs when borrowing strength from an auxiliary co-located determinand at monitoring gaps.",
    "RSN+SSN Combined Network": "Combining River Surveillance Network (RSN) and Small Streams Network (SSN) increases coverage to 1,384 observed reaches across 91,118 reaches, with RSN prioritised on shared reaches.",

    "England Reaches Snap-Fused": "Total active ORN v2 reaches with at least one snapped RSN, SSN, or citizen-science WaterBlitz observation within 350 m across 268 catchments.",
    "Co-Located Reach Dilution Fix": "On 384 reaches sampled by both EA surveillance monitors (RSN/SSN) and WaterBlitz, weighting surveillance sources higher (--prefer-source / --source-weight) prevents noisy citizen strips from diluting RSN/SSN medians.",
    "National Kriging SE Change": "Percentage drop in mean national kriging standard error across all 95,773 reaches for Nitrate (0117) and Orthophosphate (0180) after adding WaterBlitz coverage.",
    "PO4 Uncertain Tier Reaches": "Share of reaches falling in the Uncertain regulatory exceedance band for Orthophosphate drops from 51.0% to 41.7% when WaterBlitz observations reduce spatial uncertainty.",

    "Profile-MLE Box-Cox λ": "Maximum-likelihood Box-Cox power parameter λ = 0.50 matches the square-root transform, reducing Nitrate residual skewness from 1.001 to 0.118.",
    "Raw Positivity Violations": "Number of reaches where untransformed linear kriging predicts physically impossible negative Nitrate concentrations (minimum -4.817 mg/L); Sqrt, Log, and Box-Cox guarantee >= 0.",
    "Best Model-Scale vs. Raw Error": "Log achieves the highest R² on its own transformed scale (0.6784), whereas Sqrt and Box-Cox minimise back-transformed prediction error in raw mg/L (MAE = 1.57 mg/L).",
    "Log Median Bias vs. Rank": "Because exp(E[log Y]) estimates the geometric rather than arithmetic mean (Jensen's inequality), raw-scale median predictions sit -1.01 mg/L lower while preserving 0.982 rank correlation.",

    "England Variance Drop (+50)": "Percentage of total evaluated kriging variance across 90,876 England reaches eliminated by adding 50 optimally placed monitoring sites via the greedy Schur-complement solver.",
    "National Mean Kriging SE": "Reduction in average log-scale kriging standard error across England when expanding from 1,210 existing RSN sites to 1,260 (+50 optimal placements).",
    "Top Single Placement (#1)": "Highest-ranked individual new monitoring reach (RivID 18175 in Humber CatchID 80), delivering a +4.58% national variance reduction on its own.",
    "Per-Basin Gains (+50 Each)": "Independent variance reductions achieved when placing 50 optimal stations within each standalone river basin (Humber 32.3%, Thames 43.6%, Severn 43.0%).",

    "All-England Joint S3N R²": "5-fold CV R² for macroinvertebrate taxon richness (WHPT N-Taxa) across England using Joint Tail-Up + Tail-Down + Euclidean NNGP covariance.",
    "GAM + S3N Joint R²": "5-fold CV R² for WHPT N-Taxa using non-linear GAM covariates plus spatial stream-network residuals (RMSE = 5.3926 taxa vs 5.8486 for non-spatial OLS).",
    "Severn Basin Joint R²": "Regional 5-fold CV R² for WHPT N-Taxa within the Severn basin, lifting explained variance from 0.0200 (OLS) to 0.1750 via stream-network spatial structure.",
    "Inverts Samples (1,313 Sites)": "Total macroinvertebrate survey samples across 1,313 RSN sites mapped onto 1,160 unique England ORN v2 reaches.",

    "Best National RMNI R² (GAM+S3N)": "5-fold CV R² for River Macrophyte Nutrient Index (RMNI) using REML-penalised GAM fixed effects and Joint S3N spatial residuals.",
    "LMC Co-Kriging RMNI R² (Aux)": "Factor-LMC co-kriging R² for RMNI when borrowing strength from co-located Orthophosphate (0180), rescuing Thames basin CV R² from -0.2671 to +0.3165.",
    "P (0180) vs RMNI Reach r": "Reach-level Pearson correlation (r = 0.9005) between predicted Orthophosphate and RMNI across England, with 86.55% agreement across high/low risk quadrants.",
    "Summer DO (9901) Joint R²": "5-fold CV R² for summer dissolved oxygen saturation (9901); lower R² reflects strong sub-daily photosynthesis/respiration cycles in daytime grab samples.",

    "North West Joint S3N R²": "5-fold CV R² for Nitrate (0117) across 12,465 reaches and 136 observed RSN reaches in the North West River Basin District.",
    "Confidently Compliant": "Share of North West reaches (93.8%) whose upper kriging confidence bound remains safely below the 11.3 mg/L Nitrates Directive threshold.",
    "Peak Autumn Seasonal R²": "Highest seasonal Nitrate CV R² (0.6799 in Autumn) from multi-season fitting with precision-weighted dynamic temporal borrowing.",
    "Joint S3N RMSE (log mg/L)": "Out-of-sample root-mean-square error on the log Nitrate scale for the North West RBD fit with temporal borrowing (rho_t = 0.363).",

    "North West N-Taxa Joint R²": "5-fold CV R² for macroinvertebrate WHPT N-Taxa in the North West RBD, improving +0.0698 over non-spatial OLS (0.1252).",
    "Joint S3N RMSE (taxa)": "Out-of-sample cross-validation error in raw taxon count units across 131 observed North West macroinvertebrate reaches.",
    "Euclidean Variance (sigma2_Euc)": "Estimated partial sill for straight-line spatial covariance (13.92), capturing regional geology and climate gradients across adjacent catchments.",
    "Tail-Up Variance (sigma2_TU)": "Estimated partial sill for flow-connected downstream dispersion weighted by Additive Function Values (7.56, with nugget tau2 = 4.30).",

    "Ribble Median Snap Offset": "Median orthogonal distance from raw Ribble monitoring coordinates to the mapped ORN v2 river centreline (87.1% within 10 m).",
    "Weaver Snapping <= 10 m": "Percentage of included Weaver stations snapping within 10 m of an ORN v2 reach (100.0% across 23 sites, median offset 2.18 m).",
    "Mersey Snapping <= 10 m": "Percentage of included Mersey stations snapping within 10 m of an ORN v2 reach (95.7% across 23 sites, 0 excluded).",
    "Lune Snapping <= 10 m": "Percentage of included Lune stations snapping within 10 m of an ORN v2 reach (76.5% across 17 sites, median offset 3.48 m).",

    "NW Variance Reduction (K=50)": "Total evaluated kriging variance reduction across 12,465 North West reaches from adding 50 greedy optimal monitoring stations.",
    "Mean Kriging SE After K=50": "Average log-scale kriging standard error across the North West RBD after adding 50 optimal stations (-29.7% from 0.6956).",
    "Winter Peak CV R²": "Winter (DJF) achieves the highest seasonal Nitrate CV R² (0.7153) because high baseflow and winter runoff connect agricultural soils consistently along the river network.",
    "Summer Low-Flow CV R²": "Summer (JJA) Nitrate CV R² (0.5598) is lower due to reduced hydrological connectivity, uptake/denitrification, and sparser summer flow paths.",
    "Annual CV R² (2021–2024)": "Range of All-England 5-fold CV R² across individual calendar years 2021–2024 (~434–439 observed reaches per year).",
    "Temporal Trend Workflow": "Combines --by=season / --by=year slicing with --st-borrow precision weighting before running reach-level WLS trend tests (fig6) and period transitions (fig7).",
    "Existing NW RSN Stations": "Number of current EA surveillance (RSN) Nitrate monitoring reaches in the North West RBD baseline before adding new design sites.",
    "#1 Site Variance Drop": "Variance reduction delivered by the single highest-priority new North West station (RivID 38848 on the Ribble, Strahler order 6)."
  };

  function attachKpiExplainers(rootEl) {
    (rootEl || document).querySelectorAll(".kpi-card .kpi-label, .kpi-card .kpi-lbl").forEach((lblEl) => {
      if (lblEl.querySelector(".info-wrap")) return;
      const title = lblEl.textContent.trim();
      const text = KPI_EXPLAINERS[title];
      if (!text) return;
      const span = document.createElement("span");
      span.textContent = title;
      const wrap = document.createElement("span");
      wrap.className = "info-wrap";
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "info-tip-btn";
      btn.setAttribute("aria-label", `Info about ${title}`);
      btn.textContent = "i";
      const pop = document.createElement("div");
      pop.className = "info-popover";
      pop.setAttribute("role", "tooltip");
      const strong = document.createElement("strong");
      strong.textContent = title;
      pop.appendChild(strong);
      pop.appendChild(document.createTextNode(text));
      wrap.appendChild(btn);
      wrap.appendChild(pop);
      lblEl.textContent = "";
      lblEl.appendChild(span);
      lblEl.appendChild(wrap);
    });
  }

  function initInfoTooltips() {
    function positionPopover(wrap) {
      const btn = wrap.querySelector(".info-tip-btn");
      const pop = wrap.querySelector(".info-popover");
      if (!btn || !pop) return;
      const margin = 12;
      const gap = 8;
      const rect = btn.getBoundingClientRect();
      const popW = pop.offsetWidth || 280;
      const popH = pop.offsetHeight || 110;
      const vw = window.innerWidth;
      const vh = window.innerHeight;

      // Prefer right-aligned with the (i) button, clamped inside viewport
      let left = rect.right - popW;
      if (left < margin) left = margin;
      if (left + popW > vw - margin) left = Math.max(margin, vw - popW - margin);

      // Prefer below the (i) button; flip above if near bottom edge and room exists above
      let top = rect.bottom + gap;
      if (top + popH > vh - margin && rect.top - gap - popH >= margin) {
        top = rect.top - gap - popH;
      } else if (top + popH > vh - margin) {
        top = Math.max(margin, vh - popH - margin);
      }

      pop.style.left = `${Math.round(left)}px`;
      pop.style.top = `${Math.round(top)}px`;
    }

    function closeAll(exceptWrap) {
      document.querySelectorAll(".info-wrap.open").forEach((wrap) => {
        if (wrap !== exceptWrap) wrap.classList.remove("open");
      });
    }

    document.addEventListener("mouseover", (e) => {
      const wrap = e.target.closest(".info-wrap");
      if (wrap && !wrap.contains(e.relatedTarget)) positionPopover(wrap);
    });

    document.addEventListener("focusin", (e) => {
      const wrap = e.target.closest(".info-wrap");
      if (wrap) positionPopover(wrap);
    });

    document.addEventListener("click", (e) => {
      const btn = e.target.closest(".info-tip-btn");
      if (btn) {
        e.preventDefault();
        e.stopPropagation();
        const wrap = btn.closest(".info-wrap");
        if (!wrap) return;
        const willOpen = !wrap.classList.contains("open");
        closeAll(wrap);
        if (willOpen) positionPopover(wrap);
        wrap.classList.toggle("open", willOpen);
        document.querySelectorAll(".info-tip-btn").forEach((b) => {
          b.setAttribute("aria-expanded", b === btn && willOpen ? "true" : "false");
        });
        return;
      }
      if (!e.target.closest(".info-wrap")) {
        closeAll(null);
        document.querySelectorAll(".info-tip-btn").forEach((b) => b.setAttribute("aria-expanded", "false"));
      }
    });

    window.addEventListener("resize", () => {
      document.querySelectorAll(".info-wrap.open").forEach(positionPopover);
    });

    window.addEventListener("scroll", () => {
      document.querySelectorAll(".info-wrap.open").forEach(positionPopover);
    }, true);

    document.addEventListener("keydown", (e) => {
      if (e.key === "Escape") {
        closeAll(null);
        document.querySelectorAll(".info-tip-btn").forEach((b) => b.setAttribute("aria-expanded", "false"));
      }
    });
  }

  function renderOverviewBenchmark() {
    const el = document.getElementById("overview-benchmark-table");
    if (!el) return;
    const t = DATA.tables || {};
    function pickPair(rows, domainFilter) {
      const sub = (rows || []).filter(domainFilter);
      const ols = sub.find((r) => String(r.Model).includes("OLS")) || {};
      const comb = sub.find((r) => String(r.Model).includes("Combined")) || {};
      return { olsR2: ols.R2, combR2: comb.R2, rmse: comb.RMSE };
    }
    const specs = [
      { cs: "CS1", tab: "cs1", study: "England Nitrate (0117, Linear)", ...pickPair(t.cv_0117, (r) => r.Domain === "All England RSN Network") },
      { cs: "CS1", tab: "cs1", study: "England Nitrate (0117, GAM + S3N)", ...pickPair(t.cv_0117_gam, (r) => r.Domain === "All England RSN Network") },
      { cs: "CS3", tab: "cs9", study: "England Inverts WHPT N-Taxa", ...pickPair(t.cv_whpt_ntaxa, (r) => r.Domain === "All England RSN Network") },
      { cs: "CS4", tab: "cs2", study: "England Ecology WHPT ASPT", ...pickPair(t.cv_whpt, (r) => r.Domain === "All England RSN Network") },
      { cs: "CS7", tab: "cs10", study: "England Macrophyte RMNI", ...pickPair(t.cv_eut_rmni, (r) => r.Domain === "All England RSN Network") }
    ].filter((row) => Number.isFinite(row.olsR2) && Number.isFinite(row.combR2));

    el.innerHTML = `
      <table>
        <thead>
          <tr>
            <th>Case Study</th>
            <th>Model Benchmark (All England RSN)</th>
            <th class="num">OLS CV R&sup2;</th>
            <th class="num">S3N Combined CV R&sup2;</th>
            <th class="num">Spatial R&sup2; Lift</th>
            <th class="num">S3N RMSE</th>
            <th>Action</th>
          </tr>
        </thead>
        <tbody>
          ${specs.map((r) => {
            const lift = r.combR2 - r.olsR2;
            return `
              <tr>
                <td><span class="pill pill-blue">${esc(r.cs)}</span></td>
                <td>${esc(r.study)}</td>
                <td class="num">${fmt(r.olsR2, 4)}</td>
                <td class="num"><strong>${fmt(r.combR2, 4)}</strong></td>
                <td class="num" style="color:var(--accent-secondary); font-weight:600;">+${fmt(lift, 4)}</td>
                <td class="num">${fmt(r.rmse, 4)}</td>
                <td><button class="copy-btn" type="button" data-jump-tab="${esc(r.tab)}">Open &rarr;</button></td>
              </tr>
            `;
          }).join("")}
        </tbody>
      </table>
    `;
  }

  function initCaseStudyPagerForTab(tabId) {
    const btns = Array.from(document.querySelectorAll("#subtabStrip .tab-btn[data-tab]"));
    const idx = btns.findIndex((b) => b.dataset.tab === tabId);
    if (idx === -1) return;
    const panel = document.getElementById("tab-" + tabId);
    if (!panel || panel.querySelector(".cs-pager")) return;
    const prev = idx > 0 ? btns[idx - 1] : null;
    const next = idx < btns.length - 1 ? btns[idx + 1] : null;
    const pager = document.createElement("div");
    pager.className = "cs-pager";
    pager.innerHTML = `
      <div>${prev ? `<button class="cs-pager-btn" type="button" data-pager-tab="${esc(prev.dataset.tab)}">&larr; ${esc(prev.textContent)}</button>` : ""}</div>
      <div>${next ? `<button class="cs-pager-btn" type="button" data-pager-tab="${esc(next.dataset.tab)}">${esc(next.textContent)} &rarr;</button>` : ""}</div>
    `;
    panel.appendChild(pager);
  }

  function initWorkbenchLinks(rootEl) {
    (rootEl || document).querySelectorAll(".code-bar").forEach((bar) => {
      if (bar.querySelector("[data-wb-link]")) return;
      const copyBtn = bar.querySelector("[data-copy-target]");
      if (!copyBtn) return;
      const wrap = document.createElement("div");
      wrap.style.display = "flex";
      wrap.style.alignItems = "center";
      wrap.style.gap = "6px";
      const wbLink = document.createElement("a");
      wbLink.className = "wiki-btn";
      wbLink.href = "../cli_workbench/index.html";
      wbLink.setAttribute("data-wb-link", "1");
      wbLink.innerHTML = "Open in CLI Workbench &rarr;";
      copyBtn.parentNode.insertBefore(wrap, copyBtn);
      wrap.appendChild(wbLink);
      wrap.appendChild(copyBtn);
    });
  }

  function initDelegatedActions() {
    document.addEventListener("click", (e) => {
      const figBtn = e.target.closest(".fig-toggle-bar[data-group] .fig-toggle-btn[data-fig]");
      if (figBtn) {
        const bar = figBtn.closest(".fig-toggle-bar[data-group]");
        const slotId = "fig-slot-" + bar.dataset.group;
        bar.querySelectorAll(".fig-toggle-btn").forEach((b) => b.classList.remove("active"));
        figBtn.classList.add("active");
        renderFigureInto(slotId, figBtn.dataset.fig);
        return;
      }

      const jumpBtn = e.target.closest("[data-jump-tab], [data-pager-tab]");
      if (jumpBtn) {
        const targetTab = jumpBtn.dataset.jumpTab || jumpBtn.dataset.pagerTab;
        if (targetTab) switchTab(targetTab);
        return;
      }

      const copyBtn = e.target.closest("[data-copy-target]");
      if (copyBtn) {
        const target = document.getElementById(copyBtn.dataset.copyTarget);
        if (!target) return;
        const txt = target.textContent;
        const done = () => {
          const orig = copyBtn.textContent;
          copyBtn.textContent = "Copied!";
          setTimeout(() => { copyBtn.textContent = orig; }, 1400);
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(txt).then(done).catch(() => fallbackCopy(txt, done));
        } else {
          fallbackCopy(txt, done);
        }
      }
    });

    document.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") {
        const jumpEl = e.target.closest("[data-jump-tab], [data-pager-tab]");
        if (jumpEl && jumpEl.tagName !== "BUTTON" && jumpEl.tagName !== "A") {
          e.preventDefault();
          const targetTab = jumpEl.dataset.jumpTab || jumpEl.dataset.pagerTab;
          if (targetTab) switchTab(targetTab);
        }
      }
    });
  }

  function initChromeHeightSync() {
    const chrome = document.querySelector(".top-chrome");
    if (!chrome) return;
    const update = () => {
      document.documentElement.style.setProperty("--chrome-h", chrome.offsetHeight + "px");
    };
    update();
    if (typeof ResizeObserver !== "undefined") {
      new ResizeObserver(update).observe(chrome);
    } else {
      window.addEventListener("resize", update);
    }
  }

  initThemeToggle();
  initChromeHeightSync();
  populateCaseStudyData();
  renderOverviewBenchmark();
  initInfoTooltips();
  initFeatureCatalog();
  initSearch();
  initDelegatedActions();
  initLightbox();

  function syncTabFromHash() {
    const hashTab = window.location.hash.replace("#", "") || "overview";
    const target = document.getElementById("tab-" + hashTab) ? hashTab : "overview";
    switchTab(target, false);
  }

  window.addEventListener("popstate", syncTabFromHash);
  window.addEventListener("hashchange", syncTabFromHash);

  const initialHash = window.location.hash.replace("#", "");
  if (initialHash && document.getElementById("tab-" + initialHash)) {
    history.replaceState({ tab: initialHash }, "", "#" + initialHash);
    switchTab(initialHash, false);
  } else {
    ensureTabRendered("overview");
    history.replaceState({ tab: "overview" }, "", window.location.pathname + window.location.search);
  }
})();
