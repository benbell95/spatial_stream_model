#!/usr/bin/env Rscript
# ==============================================================================
# 05_descriptive_stats.R
#
# Standalone script and reusable R functions for computing comprehensive
# descriptive statistics, monitoring network representativeness tests,
# stratified environmental/topological breakdowns, and base-R diagnostic
# plots for a single S3N model output directory.
#
# Also sourced by `scripts/04_compare_model_runs.R` so multi-run comparisons
# call the exact same descriptive statistics functions (`compute_detailed_descriptives`
# and `compute_stratified_descriptives`).
#
# Usage (Single Model Output):
#   Rscript scripts/05_descriptive_stats.R --list
#   Rscript scripts/05_descriptive_stats.R output/det_0117
#   Rscript scripts/05_descriptive_stats.R output/det_0180_orn
#   Rscript scripts/05_descriptive_stats.R output/det_0117_euhydro_thames --out=output/descriptives/thames_euhydro
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
BASE_DIR <- resolve_base_dir()

source(file.path(BASE_DIR, "scripts/s3n_stats_functions.R"))

# ------------------------------------------------------------------------------
# Core Function 3: Within-Run Statistical Hypothesis Tests for a Single Output
# ------------------------------------------------------------------------------
compute_single_run_stat_tests <- function(b1) {
  ep <- b1$eng_pred
  stat_rows <- list()
  add_stat_row <- function(comp_lbl, category, test_name, metric, n1, n2, val1, val2, diff_val,
                           stat_name, stat_val, p_val, interpretation) {
    stat_rows[[length(stat_rows) + 1L]] <<- data.frame(
      Comparison     = comp_lbl,
      Category       = category,
      Test_Name      = test_name,
      Metric         = metric,
      n_1            = n1,
      n_2            = n2,
      Value_Group1   = round(val1, 5),
      Value_Group2   = round(val2, 5),
      Diff_G2_G1     = round(diff_val, 5),
      Statistic_Name = stat_name,
      Statistic_Val  = round(stat_val, 5),
      P_Value        = p_val,
      Significance   = sig_stars(p_val),
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }

  # 1. Cross-Validation Model Comparisons (S3N Combined vs OLS, Tail-Up, Tail-Down)
  if (!is.null(b1$eng_cv) && nrow(b1$eng_cv) >= 5) {
    e_comb <- b1$eng_cv$obs_y - b1$eng_cv$cv_comb
    if ("cv_ols" %in% names(b1$eng_cv)) {
      e_ols  <- b1$eng_cv$obs_y - b1$eng_cv$cv_ols
      tt_ols <- safe_paired_ttest(e_comb^2, e_ols^2)
      wt_ols <- safe_paired_wilcox(abs(e_comb), abs(e_ols))
      add_stat_row("S3N Combined vs. OLS", "Cross-Validation",
                   "Paired t-Test / Diebold-Mariano (Station Squared CV Error)",
                   "Station MSE", length(e_ols), length(e_comb),
                   mean(e_ols^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                   mean(e_comb^2, na.rm = TRUE) - mean(e_ols^2, na.rm = TRUE),
                   "t", as.numeric(tt_ols$statistic), tt_ols$p.value,
                   sprintf("S3N Combined reduces station MSE vs non-spatial OLS (%s)", sig_stars(tt_ols$p.value)))
      add_stat_row("S3N Combined vs. OLS", "Cross-Validation",
                   "Paired Wilcoxon Signed-Rank Test (Station |CV Error|)",
                   "Mean |CV Error|", length(e_ols), length(e_comb),
                   mean(abs(e_ols), na.rm = TRUE), mean(abs(e_comb), na.rm = TRUE),
                   mean(abs(e_comb), na.rm = TRUE) - mean(abs(e_ols), na.rm = TRUE),
                   "V", as.numeric(wt_ols$statistic), wt_ols$p.value,
                   sprintf("S3N Combined reduces absolute station CV error vs OLS (%s)", sig_stars(wt_ols$p.value)))
    }
    if ("cv_tu" %in% names(b1$eng_cv)) {
      e_tu  <- b1$eng_cv$obs_y - b1$eng_cv$cv_tu
      tt_tu <- safe_paired_ttest(e_comb^2, e_tu^2)
      add_stat_row("S3N Combined vs. Tail-Up", "Cross-Validation",
                   "Paired t-Test (Station Squared CV Error: Combined vs Tail-Up)",
                   "Station MSE", length(e_tu), length(e_comb),
                   mean(e_tu^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                   mean(e_comb^2, na.rm = TRUE) - mean(e_tu^2, na.rm = TRUE),
                   "t", as.numeric(tt_tu$statistic), tt_tu$p.value,
                   sprintf("Variance-weighted combination vs pure Tail-Up MSE (%s)", sig_stars(tt_tu$p.value)))
    }
    if ("cv_td" %in% names(b1$eng_cv)) {
      e_td  <- b1$eng_cv$obs_y - b1$eng_cv$cv_td
      tt_td <- safe_paired_ttest(e_comb^2, e_td^2)
      add_stat_row("S3N Combined vs. Tail-Down", "Cross-Validation",
                   "Paired t-Test (Station Squared CV Error: Combined vs Tail-Down)",
                   "Station MSE", length(e_td), length(e_comb),
                   mean(e_td^2, na.rm = TRUE), mean(e_comb^2, na.rm = TRUE),
                   mean(e_comb^2, na.rm = TRUE) - mean(e_td^2, na.rm = TRUE),
                   "t", as.numeric(tt_td$statistic), tt_td$p.value,
                   sprintf("Variance-weighted combination vs pure Tail-Down MSE (%s)", sig_stars(tt_td$p.value)))
    }
  }

  fmt_med_ratio <- function(num, den) {
    if (!is.na(num) && !is.na(den) && abs(den) > 1e-9) sprintf("%.2fx", num / den) else "N/A"
  }
  add_mw_group_test <- function(comp_lbl, category, test_name, metric, v1, v2, agg = median, fmt_interp = NULL) {
    if (length(v1) < 5 || length(v2) < 5) return(invisible(NULL))
    mw <- suppressWarnings(wilcox.test(v2, v1, paired = FALSE))
    s1 <- agg(v1, na.rm = TRUE)
    s2 <- agg(v2, na.rm = TRUE)
    interp <- if (!is.null(fmt_interp)) fmt_interp(s1, s2, mw$p.value) else sprintf("Ratio = %s (%s)", fmt_med_ratio(s2, s1), sig_stars(mw$p.value))
    add_stat_row(comp_lbl, category, test_name, metric, length(v1), length(v2),
                 s1, s2, s2 - s1, "W", as.numeric(mw$statistic), mw$p.value, interp)
  }

  # 2. Monitoring Representativeness: Observed Stations vs. Full Network Predictions
  obs_v <- ep$obs_value[ep$has_obs == 1 & !is.na(ep$obs_value)]
  if (length(obs_v) >= 5) {
    mw_rep <- suppressWarnings(wilcox.test(ep$pred_combined, obs_v, paired = FALSE))
    ks_rep <- suppressWarnings(ks.test(ep$pred_combined, obs_v))
    add_stat_row("Network vs. Monitored Sites", "Monitoring Representativeness",
                 "Mann-Whitney U Test (Full Network Predictions vs Observed Stations)",
                 sprintf("Median (%s)", b1$det_unit), length(obs_v), nrow(ep),
                 median(obs_v, na.rm = TRUE), median(ep$pred_combined, na.rm = TRUE),
                 median(ep$pred_combined, na.rm = TRUE) - median(obs_v, na.rm = TRUE),
                 "W", as.numeric(mw_rep$statistic), mw_rep$p.value,
                 sprintf("Pred/Obs median ratio = %s (%s)",
                         fmt_med_ratio(median(ep$pred_combined, na.rm = TRUE), median(obs_v, na.rm = TRUE)),
                         sig_stars(mw_rep$p.value)))
    add_stat_row("Network vs. Monitored Sites", "Monitoring Representativeness",
                 "Two-Sample Kolmogorov-Smirnov Test (Predictions vs Observed Stations)",
                 "Empirical CDF D-stat", length(obs_v), nrow(ep),
                 mean(obs_v, na.rm = TRUE), mean(ep$pred_combined, na.rm = TRUE),
                 mean(ep$pred_combined, na.rm = TRUE) - mean(obs_v, na.rm = TRUE),
                 "D", as.numeric(ks_rep$statistic), ks_rep$p.value,
                 sprintf("Distributional shift between gauged sites and full network (D = %.3f, %s)",
                         ks_rep$statistic, sig_stars(ks_rep$p.value)))
  }

  # 3. S3N Tail-Up vs. Tail-Down Component Concordance across All Reaches
  if ("pred_tu" %in% names(ep) && "pred_td" %in% names(ep)) {
    tu_v <- if ("pred_tu_sc" %in% names(ep)) ep$pred_tu_sc else ep$pred_tu
    td_v <- if ("pred_td_sc" %in% names(ep)) ep$pred_td_sc else ep$pred_td
    ct_tud <- suppressWarnings(cor.test(tu_v, td_v, method = "pearson"))
    cs_tud <- suppressWarnings(cor.test(tu_v, td_v, method = "spearman"))
    add_stat_row("Tail-Up vs. Tail-Down", "S3N Components",
                 "Pearson Correlation Test (Tail-Up vs Tail-Down Reach Predictions)",
                 "Pearson r (Spearman rho)", nrow(ep), nrow(ep),
                 as.numeric(ct_tud$estimate), as.numeric(cs_tud$estimate),
                 as.numeric(ct_tud$estimate)^2,
                 "t", as.numeric(ct_tud$statistic), ct_tud$p.value,
                 sprintf("r = %.4f (Spearman rho = %.4f, R2 = %.4f) across %s reaches (%s)",
                         ct_tud$estimate, cs_tud$estimate, ct_tud$estimate^2,
                         format(nrow(ep), big.mark = ","), sig_stars(ct_tud$p.value)))
  }

  # 4. Kriging Uncertainty at Monitored vs. Unmonitored Reaches
  se_mon   <- ep$se_model_sc[ep$has_obs == 1]
  se_unmon <- ep$se_model_sc[ep$has_obs == 0]
  add_mw_group_test(
    "Unmonitored vs. Monitored SE", "Kriging Uncertainty",
    "Mann-Whitney U Test (Kriging SE: Unmonitored vs Monitored Reaches)",
    "Mean Kriging SE", se_mon, se_unmon, agg = mean,
    fmt_interp = function(s1, s2, pv) sprintf("Unmonitored reaches have %+.4f higher mean Kriging SE than monitored reaches (%s)", s2 - s1, sig_stars(pv))
  )

  # 5. Environmental Gradient Hypothesis Tests (Topography, Hydrogeology, Stream Order)
  if ("Elevation" %in% names(ep)) {
    v_up  <- ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation >= 150.0]
    v_low <- ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation <  150.0]
    add_mw_group_test(
      "Lowland vs. Upland", "Environmental Gradient",
      "Mann-Whitney U Test (Lowland < 150 m vs Upland >= 150 m)",
      sprintf("Median (%s)", b1$det_unit), v_up, v_low, agg = median,
      fmt_interp = function(s1, s2, pv) sprintf("Lowland/Upland median ratio = %s (%s)", fmt_med_ratio(s2, s1), sig_stars(pv))
    )
  }
  if ("Is_Aquifer" %in% names(ep)) {
    v_noaq <- ep$pred_combined[!is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 0]
    v_aq   <- ep$pred_combined[!is.na(ep$Is_Aquifer) & ep$Is_Aquifer == 1]
    add_mw_group_test(
      "Aquifer vs. Non-Aquifer", "Environmental Gradient",
      "Mann-Whitney U Test (Chalk/Limestone Aquifer vs Non-Aquifer Bedrock)",
      sprintf("Median (%s)", b1$det_unit), v_noaq, v_aq, agg = median,
      fmt_interp = function(s1, s2, pv) sprintf("Aquifer/Non-Aquifer median ratio = %s (%s)", fmt_med_ratio(s2, s1), sig_stars(pv))
    )
  }
  if ("Strahler" %in% names(ep)) {
    ord_grp <- ifelse(ep$Strahler <= 2, "Headwaters", ifelse(ep$Strahler <= 4, "Mid-Order", "Mainstem"))
    v_hw <- ep$pred_combined[ord_grp == "Headwaters"]
    v_ms <- ep$pred_combined[ord_grp == "Mainstem"]
    if (length(unique(ord_grp[!is.na(ord_grp)])) >= 2) {
      kw_ord <- suppressWarnings(kruskal.test(ep$pred_combined ~ as.factor(ord_grp)))
      add_stat_row("Strahler Stream Orders", "Environmental Gradient",
                   "Kruskal-Wallis Test Across Stream Order Classes (1-2, 3-4, >=5)",
                   sprintf("Mainstem vs Headw. Med (%s)", b1$det_unit), length(v_hw), length(v_ms),
                   median(v_hw, na.rm = TRUE), median(v_ms, na.rm = TRUE),
                   median(v_ms, na.rm = TRUE) - median(v_hw, na.rm = TRUE),
                   "ChiSq", as.numeric(kw_ord$statistic), kw_ord$p.value,
                   sprintf("Mainstem/Headwater median ratio = %s (%s)",
                           fmt_med_ratio(median(v_ms, na.rm = TRUE), median(v_hw, na.rm = TRUE)),
                           sig_stars(kw_ord$p.value)))
    }
  }

  do.call(rbind, stat_rows)
}

# ------------------------------------------------------------------------------
# Single-Run Base-R Plotting Functions
# ------------------------------------------------------------------------------
plot_single_run_descriptives <- function(b1, desc_df, rep_df, strat_df, out_dir) {
  ep        <- b1$eng_pred
  det_short <- shorten_det_name(b1$det_name)
  det_unit  <- b1$det_unit
  u_lbl     <- if (b1$trans_type == "log") sprintf("log(%s)", det_unit) else if (b1$trans_type == "none") det_unit else sprintf("%s(%s)", b1$trans_type, det_unit)

  # Figure 1: 4-Panel Distributions, Representativeness, Tail-Up vs Tail-Down & Uncertainty
  fig1_file <- file.path(out_dir, "fig_desc1_distributions_and_representativeness.png")
  open_png_device(fig1_file, width = 2600, height = 2000, res = 200)
  tryCatch({
    par(mfrow = c(2, 2), mar = c(4.8, 4.8, 4.0, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Kernel Density of Full Network Predictions vs. Observed Monitoring Stations
  obs_sc   <- ep$obs_model_sc[ep$has_obs == 1 & !is.na(ep$obs_model_sc)]
  unmon_sc <- ep$pred_model_sc[ep$has_obs == 0]
  d_all    <- density(ep$pred_model_sc, na.rm = TRUE, n = 512)
  d_unmon  <- if (length(unmon_sc) >= 10) density(unmon_sc, na.rm = TRUE, n = 512) else d_all
  d_obs    <- if (length(obs_sc) >= 10)   density(obs_sc,   na.rm = TRUE, n = 512) else d_all

  xr_a <- range(c(d_all$x, d_obs$x), na.rm = TRUE)
  yr_a <- c(0, max(c(d_all$y, d_obs$y), na.rm = TRUE) * 1.22)
  plot(NULL, xlim = xr_a, ylim = yr_a,
       xlab = sprintf("%s [%s]", det_short, u_lbl),
       ylab = "Probability Density",
       main = "",
       cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  polygon(c(d_all$x, rev(d_all$x)), c(d_all$y, rep(0, length(d_all$y))),
          col = adjustcolor("#1f77b4", alpha.f = 0.22), border = NA)
  lines(d_all$x,   d_all$y,   col = "#1f77b4", lwd = 2.6, lty = 1)
  lines(d_unmon$x, d_unmon$y, col = "#2ca02c", lwd = 2.0, lty = 2)
  lines(d_obs$x,   d_obs$y,   col = "#d62728", lwd = 2.6, lty = 1)
  abline(v = median(ep$pred_model_sc, na.rm = TRUE), col = "#1f77b4", lty = 3, lwd = 1.4)
  abline(v = median(obs_sc, na.rm = TRUE),           col = "#d62728", lty = 3, lwd = 1.4)
  leg_pos_a <- if (median(ep$pred_model_sc, na.rm = TRUE) > mean(xr_a)) "topleft" else "topright"
  legend(leg_pos_a,
         legend = c(
           sprintf("All Reaches (Med = %.3g %s, n = %s)", median(ep$pred_combined, na.rm = TRUE), det_unit, format(nrow(ep), big.mark = ",")),
           sprintf("Unmonitored Reaches (Med = %.3g %s)", median(ep$pred_combined[ep$has_obs == 0], na.rm = TRUE), det_unit),
           sprintf("Observed Stations (Med = %.3g %s, n = %d)", median(ep$obs_value[ep$has_obs == 1], na.rm = TRUE), det_unit, length(obs_sc))
         ),
         col = c("#1f77b4", "#2ca02c", "#d62728"), lty = c(1, 2, 1), lwd = c(2.6, 2.0, 2.6),
         bty = "n", cex = 0.78)
  sub_a_txt <- if (!is.null(rep_df) && nrow(rep_df) > 0) {
    sprintf("Pred/Obs Median Ratio = %.2fx | Mann-Whitney p = %s (%s) | KS D = %.3f (p = %s)",
            rep_df$Pred_to_Obs_Med_Ratio[1], fmt_pval(rep_df$MW_P_Value[1]),
            rep_df$Significance[1], rep_df$KS_D_Stat[1], fmt_pval(rep_df$KS_P_Value[1]))
  } else NULL
  draw_plot_header(
    main_title = "A. Network Predictions vs. Observed Monitoring Stations",
    sub_title  = sub_a_txt,
    max_main = 60, max_sub = 100, cex_main = 1.10, cex_sub = 0.76, adj = 0.5
  )

  # Panel B: Empirical CDF Comparing S3N Combined, Tail-Up, Tail-Down & Observed Stations
  tu_sc <- if ("pred_tu_sc" %in% names(ep)) ep$pred_tu_sc else ep$pred_tu
  td_sc <- if ("pred_td_sc" %in% names(ep)) ep$pred_td_sc else ep$pred_td
  xr_b  <- quantile(c(ep$pred_model_sc, obs_sc), probs = c(0.005, 0.995), na.rm = TRUE)
  xs_b  <- seq(xr_b[1], xr_b[2], length.out = 300)
  plot(NULL, xlim = xr_b, ylim = c(0, 1),
       xlab = sprintf("%s [%s]", det_short, u_lbl),
       ylab = "Cumulative Proportion",
       main = wrap_title("B. Empirical CDF: S3N Combined, Tail-Up, Tail-Down & Observed", 60),
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  lines(xs_b, ecdf(ep$pred_model_sc)(xs_b), col = "#1f77b4", lwd = 2.5, lty = 1)
  lines(xs_b, ecdf(tu_sc)(xs_b),            col = "#2ca02c", lwd = 2.0, lty = 2)
  lines(xs_b, ecdf(td_sc)(xs_b),            col = "#9467bd", lwd = 2.0, lty = 4)
  if (length(obs_sc) >= 5) {
    lines(xs_b, ecdf(obs_sc)(xs_b),         col = "#d62728", lwd = 2.4, lty = 1)
  }
  legend("bottomright",
         legend = c(
           sprintf("S3N Combined (Med = %.3g %s)", median(ep$pred_combined, na.rm = TRUE), det_unit),
           sprintf("S3N Tail-Up (Med = %.3g %s)",  median(ep$pred_tu, na.rm = TRUE), det_unit),
           sprintf("S3N Tail-Down (Med = %.3g %s)",median(ep$pred_td, na.rm = TRUE), det_unit),
           sprintf("Observed Stations (Med = %.3g %s)", median(ep$obs_value[ep$has_obs == 1], na.rm = TRUE), det_unit)
         ),
         col = c("#1f77b4", "#2ca02c", "#9467bd", "#d62728"),
         lty = c(1, 2, 4, 1), lwd = c(2.5, 2.0, 2.0, 2.4), bty = "n", cex = 0.78)

  # Panel C: Reach-by-Reach Scatter of S3N Tail-Up vs. S3N Tail-Down Predictions
  set.seed(42)
  idx_pt <- if (nrow(ep) > 8500) sample(nrow(ep), 8500) else seq_len(nrow(ep))
  lim_c  <- quantile(c(tu_sc, td_sc), probs = c(0.005, 0.995), na.rm = TRUE)
  plot(tu_sc[idx_pt], td_sc[idx_pt],
       xlim = lim_c, ylim = lim_c, pch = 16, cex = 0.45, col = rgb(0.16, 0.50, 0.73, 0.25),
       xlab = sprintf("S3N Tail-Up Prediction [%s]", u_lbl),
       ylab = sprintf("S3N Tail-Down Prediction [%s]", u_lbl),
       main = wrap_title(sprintf("C. S3N Tail-Up vs. Tail-Down Component Agreement (n = %s)",
                                format(nrow(ep), big.mark = ",")), 60),
       cex.main = 1.10, cex.lab = 0.95, las = 1)
  grid(col = "#e0e0e0", lty = "dotted")
  abline(0, 1, col = "#c0392b", lwd = 2.0, lty = 2)
  fit_c <- lm(td_sc ~ tu_sc)
  abline(fit_c, col = "#111111", lwd = 2.0, lty = 1)
  r_tud   <- cor(tu_sc, td_sc, use = "complete.obs", method = "pearson")
  rho_tud <- cor(tu_sc, td_sc, use = "complete.obs", method = "spearman")
  rmsd_c  <- sqrt(mean((td_sc - tu_sc)^2, na.rm = TRUE))
  legend("topleft",
         legend = c(
           sprintf("Pearson r = %.4f (Spearman rho = %.4f)", r_tud, rho_tud),
           sprintf("RMSD = %.4f | OLS Slope = %.3f", rmsd_c, coef(fit_c)[2]),
           "OLS Regression Fit",
           "1:1 Agreement Line"
         ),
         col = c(NA, NA, "#111111", "#c0392b"), lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 2),
         bty = "n", cex = 0.80)

  # Panel D: Kriging Standard Error (se_model_sc) by Monitoring Status & Stream Order
  se_groups <- list(
    "Monitored\nReaches"   = ep$se_model_sc[ep$has_obs == 1],
    "Unmonitored\nReaches" = ep$se_model_sc[ep$has_obs == 0],
    "Order 1-2\n(Headw.)"  = if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler <= 2] else numeric(0),
    "Order 3-4\n(Medium)"  = if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler %in% c(3, 4)] else numeric(0),
    "Order >=5\n(Mainstem)"= if ("Strahler" %in% names(ep)) ep$se_model_sc[!is.na(ep$Strahler) & ep$Strahler >= 5] else numeric(0)
  )
  se_groups <- se_groups[sapply(se_groups, length) >= 5]
  cols_d <- c("#27ae60", "#e67e22", "#3498db", "#2980b9", "#1f618d")[seq_along(se_groups)]
  bp_tmp <- boxplot(se_groups, plot = FALSE)
  ylim_d <- c(min(bp_tmp$stats[1, ], na.rm = TRUE) * 0.92,
              max(bp_tmp$stats[5, ], na.rm = TRUE) * 1.14)
  bp_d <- boxplot(se_groups, col = adjustcolor(cols_d, alpha.f = 0.65), border = "#222222",
                  outline = FALSE, las = 1, ylim = ylim_d,
                  ylab = sprintf("Kriging Standard Error (%s)", if (b1$trans_type == "none") det_unit else sprintf("%s scale", b1$trans_type)),
                  main = wrap_title("D. Kriging Prediction Uncertainty (se_model_sc) Breakdown", 60),
                  cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.84)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  se_mns <- sapply(se_groups, mean, na.rm = TRUE)
  points(seq_along(se_groups), se_mns, pch = 18, cex = 1.4, col = "#111111")
  text(seq_along(se_groups), bp_d$stats[5, ], labels = sprintf("Mean=%.3f", se_mns), pos = 3, cex = 0.76, font = 2)

    mtext(wrap_title(sprintf("S3N Single-Run Descriptive Diagnostics: %s (%s — %s)",
                             b1$det_name, b1$label, b1$dir_name), 180),
          side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
  }, finally = {
    if (dev.cur() > 1L) dev.off()
  })

  # Figure 2: 4-Panel Stratified Environmental & Topological Gradients
  fig2_file <- file.path(out_dir, "fig_desc2_stratified_gradients.png")
  open_png_device(fig2_file, width = 2600, height = 2000, res = 200)
  tryCatch({
    par(mfrow = c(2, 2), mar = c(5.0, 4.8, 3.6, 1.8), oma = c(0.5, 0.5, 2.2, 0.5), bg = "white")

  # Panel A: Boxplot by Combined Landscape Regime (Upland, Lowland Non-Aquifer, Lowland Aquifer)
  if (("Elevation" %in% names(ep)) && ("Is_Aquifer" %in% names(ep))) {
    reg_list <- list(
      "Upland\n(>= 150 m)"         = ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation >= 150.0],
      "Lowland\nNon-Aquifer"       = ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 0],
      "Lowland\nAquifer (Chalk/Lm)"= ep$pred_model_sc[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 1]
    )
    raw_meds_a <- c(
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation >= 150.0], na.rm = TRUE),
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 0], na.rm = TRUE),
      median(ep$pred_combined[!is.na(ep$Elevation) & ep$Elevation < 150.0 & ep$Is_Aquifer == 1], na.rm = TRUE)
    )
    bp_t <- boxplot(reg_list, plot = FALSE)
    ylim_a <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.1 * diff(range(bp_t$stats, na.rm = TRUE)),
                max(bp_t$stats[5, ], na.rm = TRUE) + 0.18 * diff(range(bp_t$stats, na.rm = TRUE)))
    bp_a <- boxplot(reg_list, col = adjustcolor(c("#2980b9", "#27ae60", "#e67e22"), alpha.f = 0.65),
                    border = "#222222", outline = FALSE, las = 1, ylim = ylim_a,
                    ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                    main = wrap_title("A. Predictions Across Landscape & Hydrogeological Regimes", 60),
                    cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.88)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    text(1:3, bp_a$stats[5, ], labels = sprintf("Med = %.3g %s\n(n = %s)", raw_meds_a, det_unit,
                                                format(sapply(reg_list, length), big.mark = ",")),
         pos = 3, cex = 0.76, font = 2)
  } else {
    plot.new()
  }

  # Panel B: Boxplot Across Individual Strahler Stream Orders
  if ("Strahler" %in% names(ep)) {
    ord_capped <- pmin(6L, pmax(1L, as.integer(ep$Strahler)))
    ord_levels <- sort(unique(ord_capped[!is.na(ord_capped)]))
    ord_list   <- lapply(ord_levels, function(o) ep$pred_model_sc[ord_capped == o])
    names(ord_list) <- ifelse(ord_levels == 6L, "Order >=6", paste0("Order ", ord_levels))
    raw_meds_b <- sapply(ord_levels, function(o) median(ep$pred_combined[ord_capped == o], na.rm = TRUE))
    bp_t <- boxplot(ord_list, plot = FALSE)
    ylim_b <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.08 * diff(range(bp_t$stats, na.rm = TRUE)),
                max(bp_t$stats[5, ], na.rm = TRUE) + 0.16 * diff(range(bp_t$stats, na.rm = TRUE)))
    pal_b <- colorRampPalette(c("#85c1e9", "#2e86c1", "#1b4f72"))(length(ord_list))
    bp_b <- boxplot(ord_list, col = adjustcolor(pal_b, alpha.f = 0.70), border = "#222222",
                    outline = FALSE, las = 1, ylim = ylim_b,
                    ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                    main = wrap_title("B. Longitudinal Gradient Across Strahler Stream Orders", 60),
                    cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.86)
    grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
    text(seq_along(ord_list), bp_b$stats[5, ], labels = sprintf("%.3g", raw_meds_b),
         pos = 3, cex = 0.76, font = 2)
  } else {
    plot.new()
  }

  # Panel C: Boxplot Across Major River Basins
  bas_masks <- list(
    "All England" = rep(TRUE, nrow(ep)),
    "Humber"      = if (!is.na(b1$hum_cid)) (ep$CatchID == b1$hum_cid) else rep(FALSE, nrow(ep)),
    "Thames"      = if (!is.na(b1$tha_cid)) (ep$CatchID == b1$tha_cid) else rep(FALSE, nrow(ep)),
    "Severn"      = if (!is.na(b1$sev_cid)) (ep$CatchID == b1$sev_cid) else rep(FALSE, nrow(ep))
  )
  bas_masks <- bas_masks[sapply(bas_masks, sum) >= 5]
  bas_list  <- lapply(bas_masks, function(m) ep$pred_model_sc[m])
  raw_meds_c <- sapply(bas_masks, function(m) median(ep$pred_combined[m], na.rm = TRUE))
  bp_t <- boxplot(bas_list, plot = FALSE)
  ylim_c <- c(min(bp_t$stats[1, ], na.rm = TRUE) - 0.08 * diff(range(bp_t$stats, na.rm = TRUE)),
              max(bp_t$stats[5, ], na.rm = TRUE) + 0.16 * diff(range(bp_t$stats, na.rm = TRUE)))
  cols_c <- c("#34495e", "#2980b9", "#d35400", "#27ae60")[seq_along(bas_list)]
  bp_c <- boxplot(bas_list, col = adjustcolor(cols_c, alpha.f = 0.65), border = "#222222",
                  outline = FALSE, las = 1, ylim = ylim_c,
                  ylab = sprintf("Predicted %s [%s]", det_short, u_lbl),
                  main = wrap_title("C. Regional Comparison Across Major River Basins", 60),
                  cex.main = 1.10, cex.lab = 0.95, cex.axis = 0.88)
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  text(seq_along(bas_list), bp_c$stats[5, ],
       labels = sprintf("Med = %.3g %s\n(n = %s)", raw_meds_c, det_unit,
                        format(sapply(bas_list, length), big.mark = ",")),
       pos = 3, cex = 0.75, font = 2)

  # Panel D: Binned Decile Profile of Predicted Value vs. Elevation (m)
  if ("Elevation" %in% names(ep) && sum(!is.na(ep$Elevation)) >= 50) {
    q_el <- unique(quantile(ep$Elevation, probs = seq(0, 1, by = 0.1), na.rm = TRUE))
    if (length(q_el) >= 4) {
      el_bin <- cut(ep$Elevation, breaks = q_el, include.lowest = TRUE)
      el_mid <- tapply(ep$Elevation,     el_bin, median, na.rm = TRUE)
      y_med  <- tapply(ep$pred_combined, el_bin, median, na.rm = TRUE)
      y_q25  <- tapply(ep$pred_combined, el_bin, quantile, probs = 0.25, na.rm = TRUE)
      y_q75  <- tapply(ep$pred_combined, el_bin, quantile, probs = 0.75, na.rm = TRUE)
      ok <- !is.na(el_mid) & !is.na(y_med)
      el_mid <- as.numeric(el_mid[ok]); y_med <- as.numeric(y_med[ok])
      y_q25  <- as.numeric(y_q25[ok]);  y_q75 <- as.numeric(y_q75[ok])

      plot(el_mid, y_med, type = "n",
           xlim = range(el_mid, na.rm = TRUE),
           ylim = c(0, max(y_q75, na.rm = TRUE) * 1.15),
           xlab = "Reach Elevation Decile Midpoint (m AOD)",
           ylab = sprintf("Predicted %s (%s, Median & IQR)", det_short, det_unit),
           main = wrap_title("D. Altitudinal Gradient (Decile Median & IQR Profile)", 60),
           cex.main = 1.10, cex.lab = 0.95, las = 1)
      grid(col = "#e0e0e0", lty = "dotted")
      polygon(c(el_mid, rev(el_mid)), c(y_q25, rev(y_q75)),
              col = adjustcolor("#1f77b4", alpha.f = 0.22), border = NA)
      lines(el_mid, y_med, col = "#1f77b4", lwd = 2.6, type = "b", pch = 16, cex = 1.1)
      abline(v = 150.0, col = "#c0392b", lwd = 1.6, lty = 2)
      legend("topright",
             legend = c("Decile Median Prediction", "Interquartile Range (Q25 - Q75)", "150 m Upland Threshold"),
             col = c("#1f77b4", adjustcolor("#1f77b4", alpha.f = 0.35), "#c0392b"),
             lwd = c(2.6, 8, 1.6), lty = c(1, 1, 2), pch = c(16, NA, NA), bty = "n", cex = 0.80)
    } else {
      plot.new()
    }
  } else {
    plot.new()
  }

    mtext(wrap_title(sprintf("S3N Stratified Environmental & Topological Gradients: %s (%s)",
                             b1$det_name, b1$label), 180),
          side = 3, line = 0.4, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
  }, finally = {
    if (dev.cur() > 1L) dev.off()
  })

  invisible(c(fig1_file, fig2_file))
}

# ------------------------------------------------------------------------------
# Single-Run Markdown Report Writer
# ------------------------------------------------------------------------------
write_descriptive_markdown_report <- function(b1, desc_df, rep_df, strat_df, stat_df, out_dir, report_path) {
  mdw      <- new_md_writer(256L)
  add_line <- mdw$add_line

  add_line("# S3N Single-Run Descriptive Statistics & Diagnostic Report")
  add_line("")
  add_line("- **Determinand**: `%s` — **%s** (`%s`, %s scale)",
           b1$det_code, b1$det_name, b1$det_unit, b1$trans_type)
  add_line("- **Data Transformation**: **%s** (`%s`%s)",
           b1$trans_label, b1$trans_type,
           if (!is.na(b1$trans_lambda)) sprintf(", λ = %.4f", b1$trans_lambda) else "")
  add_line("- **Water Network**: %s (`%s` reaches across `%d` catchments, `%d` observed monitoring sites)",
           b1$net_label, format(nrow(b1$eng_pred), big.mark = ","),
           length(unique(b1$eng_pred$CatchID)), sum(b1$eng_pred$has_obs == 1))
  add_line("- **Model Run Directory**: `%s`", b1$dir_name)
  add_line("- **Fixed-Effect Covariates**: %s", b1$vars_label)
  add_line("- **Output Directory**: `%s`", format_rel_display_path(out_dir, BASE_DIR))
  add_line("")

  add_line("## 1. Full Network & Monitoring Station Descriptive Statistics")
  add_line("")
  add_line("| Variable Group | Unit | N | Mean ± SD | CV (%) | Geom Mean | Skew | Kurt | Min | P05 | Q25 | Median | Q75 | P95 | Max | IQR | MAD |")
  add_line("| :--- | :---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(desc_df))) {
    d <- desc_df[r_i, ]
    gm_str <- if (is.na(d$Geom_Mean)) "—" else sprintf("%.4g", d$Geom_Mean)
    cv_str <- if (is.na(d$CV_Pct))    "—" else sprintf("%.1f%%", d$CV_Pct)
    add_line("| **%s** | `%s` | %s | `%.4g ± %.4g` | `%s` | `%s` | `%+.2f` | `%+.2f` | `%.4g` | `%.4g` | `%.4g` | **`%.4g`** | `%.4g` | `%.4g` | `%.4g` | `%.4g` | `%.4g` |",
             d$Variable_Group, d$Unit, format(d$N, big.mark = ","),
             d$Mean, d$SD, cv_str, gm_str, d$Skewness, d$Excess_Kurtosis,
             d$Min, d$P05, d$Q25, d$Median, d$Q75, d$P95, d$Max, d$IQR, d$MAD)
  }
  add_line("")

  if (!is.null(rep_df) && nrow(rep_df) > 0) {
    add_line("## 2. Monitoring Network Representativeness (Observed Stations vs. Full River Network)")
    add_line("")
    add_line("| Monitored Sites | Network Reaches | Obs Median (Mean) | Network Pred Median (Mean) | Pred / Obs Median Ratio | Mann–Whitney p | KS D-stat (p) | Interpretation |")
    add_line("| ---: | ---: | :--- | :--- | ---: | ---: | ---: | :--- |")
    for (r_i in seq_len(nrow(rep_df))) {
      rp <- rep_df[r_i, ]
      add_line("| %d | %s | `%.4g` (`%.4g`) | `%.4g` (`%.4g`) | **`%.3fx`** | `%s` (`%s`) | `%.3f` (`%s`) | %s |",
               rp$n_Obs_Stations, format(rp$n_Total_Reaches, big.mark = ","),
               rp$Obs_Median_Raw, rp$Obs_Mean_Raw, rp$Pred_Median_Raw, rp$Pred_Mean_Raw,
               rp$Pred_to_Obs_Med_Ratio, fmt_pval(rp$MW_P_Value), rp$Significance,
               rp$KS_D_Stat, fmt_pval(rp$KS_P_Value), rp$Interpretation)
    }
    add_line("")
  }

  if (!is.null(strat_df) && nrow(strat_df) > 0) {
    add_line("## 3. Stratified Model Output Breakdown (Topography, Hydrogeology, Stream Order & Major Basins)")
    add_line("")
    add_line("| Dimension | Stratum | Reaches (% Net) | Obs Sites (Med) | Pred Mean ± SD | Pred Median (IQR) | Pred 5th–95th Pct | Mean Z-Score | Kriging SE Med | Fold Ratio vs Ref | Within-Run Test (p-val) |")
    add_line("| :--- | :--- | ---: | ---: | :--- | :--- | :--- | ---: | ---: | ---: | :--- |")
    for (r_i in seq_len(nrow(strat_df))) {
      st <- strat_df[r_i, ]
      obs_s <- if (is.na(st$Obs_Median_Raw)) sprintf("%d (—)", st$n_obs) else sprintf("%d (`%.3g`)", st$n_obs, st$Obs_Median_Raw)
      rat_s <- if (is.na(st$Ratio_vs_Ref)) "—" else sprintf("%.2fx", st$Ratio_vs_Ref)
      tst_s <- if (is.na(st$Within_Run_P)) st$Within_Run_Test else sprintf("%s: p = %s (%s)", st$Within_Run_Test, fmt_pval(st$Within_Run_P), st$Within_Run_Sig)
      add_line("| %s | **%s** | %s (`%.1f%%`) | %s | `%.4g ± %.4g` | **`%.4g`** (`%.4g`) | `%.4g` – `%.4g` | `%+.2f` | `%.4f` | `%s` | %s |",
               st$Dimension, st$Stratum, format(st$n_reaches, big.mark = ","), st$Pct_Network,
               obs_s, st$Pred_Mean_Raw, st$Pred_SD_Raw, st$Pred_Median_Raw, st$Pred_IQR_Raw,
               st$Pred_P05_Raw, st$Pred_P95_Raw, st$Pred_Mean_Z, st$SE_Median, rat_s, tst_s)
    }
    add_line("")
  }

  if (!is.null(stat_df) && nrow(stat_df) > 0) {
    add_line("## 4. Within-Run Statistical Hypothesis Tests")
    add_line("")
    add_line("| Comparison | Category | Statistical Test | Metric | Group 1 Val | Group 2 Val | Diff / R² | Stat | p-value | Sig | Finding |")
    add_line("| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- | ---: | :---: | :--- |")
    for (r_i in seq_len(nrow(stat_df))) {
      s <- stat_df[r_i, ]
      add_line("| %s | %s | %s | %s | `%.4f` | `%.4f` | `%+.4f` | `%s = %.3f` | `%s` | `%s` | %s |",
               s$Comparison, s$Category, s$Test_Name, s$Metric,
               s$Value_Group1, s$Value_Group2, s$Diff_G2_G1,
               s$Statistic_Name, s$Statistic_Val, fmt_pval(s$P_Value),
               s$Significance, s$Interpretation)
    }
    add_line("")
  }

  add_line("## 5. Cross-Validation Performance Summary")
  add_line("")
  add_line("| Domain | Model | n_obs | n_reaches | CV R² | CV RMSE | CV MAE |")
  add_line("| :--- | :--- | ---: | ---: | ---: | ---: | ---: |")
  for (r_i in seq_len(nrow(b1$cv_tab))) {
    r <- b1$cv_tab[r_i, ]
    add_line("| %s | %s | %d | %s | `%.4f` | `%.4f` | `%.4f` |",
             r$Canon_Domain, r$Canon_Model, r$n_obs,
             format(r$n_reaches, big.mark = ","), r$R2, r$RMSE, r$MAE)
  }
  add_line("")

  add_line("## 6. Generated Files")
  add_line("")
  add_line("- `descriptive_report.md` — Full single-run descriptive & diagnostic Markdown report.")
  add_line("- `descriptive_stats_summary.csv` — 18-column distributional summary across all model output and observation variables.")
  add_line("- `monitoring_representativeness.csv` — Mann–Whitney U and KS test comparing gauged stations against the full stream network.")
  add_line("- `stratified_results_summary.csv` — Stratified descriptive statistics across topography, hydrogeology, stream order, and basins.")
  add_line("- `within_run_statistical_tests.csv` — Formal statistical tests for CV model improvements, component agreement, and environmental gradients.")
  add_line("- `fig_desc1_distributions_and_representativeness.png` — 4-panel base-R plot of distributions, ECDF, Tail-Up vs Tail-Down agreement, and Kriging SE.")
  add_line("- `fig_desc2_stratified_gradients.png` — 4-panel base-R plot of landscape regimes, Strahler stream orders, major river basins, and altitudinal profile.")
  add_line("")

  writeLines(mdw$get_lines(), report_path)
}

# ------------------------------------------------------------------------------
# Main Single-Run Entry Point: run_descriptive_stats()
# ------------------------------------------------------------------------------
run_descriptive_stats <- function(run_dir, label = NULL, out_dir = NULL, out_root = NULL, quiet = FALSE) {
  if (is.null(out_root) || !nzchar(out_root)) {
    out_root <- load_directory_config(BASE_DIR)$output_dir
  }
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(rel_path(out_root, "descriptives"), recursive = TRUE, showWarnings = FALSE)
  b1 <- load_model_run_bundle(run_dir, custom_label = label, out_root = out_root)

  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- rel_path(b1$dir_path, "descriptives")
  } else if (!is_abs_path(out_dir)) {
    out_dir <- rel_path(BASE_DIR, out_dir)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!quiet) {
    cat("==============================================================================\n")
    cat(sprintf("Running Single-Output Descriptive Statistics & Diagnostics:\n"))
    cat(sprintf("  Run Label:    %s (%s)\n", b1$label, b1$dir_name))
    cat(sprintf("  Determinand:  %s [%s] (%s, %s scale)\n",
                b1$det_name, b1$det_code, b1$det_unit, b1$trans_type))
    cat(sprintf("  Network:      %s (%s reaches, %d observed stations)\n",
                b1$net_label, format(nrow(b1$eng_pred), big.mark = ","), sum(b1$eng_pred$has_obs == 1)))
    cat(sprintf("  Output Dir:   %s\n", out_dir))
    cat("==============================================================================\n")
  }

  desc_bundle <- compute_detailed_descriptives(list(b1))
  desc_df     <- desc_bundle$desc_df
  rep_df      <- desc_bundle$rep_df
  strat_df    <- compute_stratified_descriptives(list(b1))
  stat_df     <- compute_single_run_stat_tests(b1)

  desc_csv  <- file.path(out_dir, "descriptive_stats_summary.csv")
  rep_csv   <- file.path(out_dir, "monitoring_representativeness.csv")
  strat_csv <- file.path(out_dir, "stratified_results_summary.csv")
  stat_csv  <- file.path(out_dir, "within_run_statistical_tests.csv")
  rep_md    <- file.path(out_dir, "descriptive_report.md")

  fast_write_csv(desc_df,  desc_csv)
  if (!is.null(rep_df)) fast_write_csv(rep_df, rep_csv)
  fast_write_csv(strat_df, strat_csv)
  fast_write_csv(stat_df,  stat_csv)

  fig_paths <- plot_single_run_descriptives(b1, desc_df, rep_df, strat_df, out_dir)
  write_descriptive_markdown_report(b1, desc_df, rep_df, strat_df, stat_df, out_dir, rep_md)

  if (!quiet) {
    cat("\nDescriptive analysis complete! Generated files:\n")
    cat(sprintf("  - Markdown Report:      %s\n", rep_md))
    cat(sprintf("  - Descriptives CSV:     %s\n", desc_csv))
    if (!is.null(rep_df)) {
      cat(sprintf("  - Representativeness:   %s\n", rep_csv))
    }
    cat(sprintf("  - Stratified CSV:       %s\n", strat_csv))
    cat(sprintf("  - Statistical Tests:    %s\n", stat_csv))
    cat(sprintf("  - Figure 1 (Dist/SE):   %s\n", fig_paths[1]))
    cat(sprintf("  - Figure 2 (Gradients): %s\n", fig_paths[2]))
  }

  invisible(list(
    bundle      = b1,
    out_dir     = out_dir,
    report_path = rep_md,
    desc_df     = desc_df,
    rep_df      = rep_df,
    strat_df    = strat_df,
    stat_df     = stat_df,
    fig_paths   = fig_paths
  ))
}

# ==============================================================================
# Command-Line Interface (only runs when 05_descriptive_stats.R is invoked directly)
# ==============================================================================
if (sys.nframe() == 0L) {
  args    <- commandArgs(trailingOnly = TRUE)
  dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat("Usage:\n")
    cat("  Rscript scripts/05_descriptive_stats.R --list [--output-dir=<dir>]\n")
    cat("  Rscript scripts/05_descriptive_stats.R <output_dir> [--label=\"Custom Label\"] [--out=<dir>] [--output-dir=<dir>]\n")
    cat("Examples:\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0117\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0180_orn\n")
    cat("  Rscript scripts/05_descriptive_stats.R output/det_0117_euhydro_thames\n")
    quit(save = "no", status = 0)
  }

  if (any(args %in% c("--list", "--list-runs"))) {
    list_available_model_runs(out_root = dir_cfg$out_dir)
    quit(save = "no", status = 0)
  }

  pos_dirs  <- character(0)
  label_cli <- NULL
  out_cli   <- NULL

  for (a in args) {
    if (grepl("^--dir=", a)) {
      pos_dirs <- c(pos_dirs, sub("^--dir=", "", a))
    } else if (grepl("^--label=", a)) {
      label_cli <- sub("^--label=", "", a)
    } else if (grepl("^--(out|output)=", a)) {
      out_cli <- sub("^--(out|output)=", "", a)
    } else if (grepl("^--(output-dir|out-dir|data-dir|data|gis-dir|gis)=", a)) {
      # Handled via parse_custom_dir_args()
    } else if (!grepl("^--", a)) {
      pos_dirs <- c(pos_dirs, a)
    }
  }

  if (length(pos_dirs) == 0) {
    stop("Please provide a model output directory (e.g., Rscript scripts/05_descriptive_stats.R output/det_0117).")
  }

  for (d_i in seq_along(pos_dirs)) {
    out_i <- if (!is.null(out_cli) && length(pos_dirs) > 1) {
      rel_path(out_cli, basename(pos_dirs[d_i]))
    } else {
      out_cli
    }
    run_descriptive_stats(run_dir = pos_dirs[d_i], label = label_cli, out_dir = out_i, out_root = dir_cfg$out_dir)
  }
}
