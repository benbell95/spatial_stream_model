#!/usr/bin/env Rscript
# ==============================================================================
# 07_compare_cokriging_runs.R
#
# Statistical & Spatial Comparison of Two Multivariate Co-Kriging Runs
# (Designed for Case Study 6: Equal Source Treatment vs. Favouring RSN + SSN
# over Citizen Science WaterBlitz).
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
source(file.path(BASE_DIR, "scripts/s3n_multivariate_functions.R"))
BASE_DIR <- resolve_base_dir()

get_cli_val <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0L) sub(paste0("^", flag, "="), "", hit[length(hit)]) else default
}

eq_dir    <- get_cli_val("--eq-dir", file.path(BASE_DIR, "output/cs6_equal/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit"))
wt_dir    <- get_cli_val("--wt-dir", file.path(BASE_DIR, "output/cs6_weighted/results/cokrig_0117_0180_orn_rsn+ssn+waterblitz_tobit"))
sites_csv <- get_cli_val("--sites-csv", file.path(BASE_DIR, "data/processed/orn_rsn+ssn+waterblitz_tobit/rsn_sites.csv"))
out_dir   <- get_cli_val("--out-dir", file.path(BASE_DIR, "output/comparisons/cs6_equal_vs_favoured"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Comparing Co-Kriging Runs:\n")
cat("  Run 1 (Equal Sources)     :", eq_dir, "\n")
cat("  Run 2 (Favoured RSN+SSN)  :", wt_dir, "\n")
cat("  Output Directory          :", out_dir, "\n")

# ------------------------------------------------------------------------------
# 1. Reconstruct per-reach source membership from rsn_sites.csv
# ------------------------------------------------------------------------------
sites_df <- read.csv(sites_csv, stringsAsFactors = FALSE)
sites_df <- sites_df[sites_df$included_250m == 1, , drop = FALSE]
sites_df$src_low <- tolower(trimws(as.character(sites_df$source)))
sites_df$is_rsn <- grepl("river surveillance|rsn", sites_df$src_low)
sites_df$is_ssn <- grepl("small streams|ssn", sites_df$src_low)
sites_df$is_stat <- sites_df$is_rsn | sites_df$is_ssn
sites_df$is_wb  <- grepl("waterblitz|waterb", sites_df$src_low)

stat_rivs <- unique(sites_df$RivID[sites_df$is_stat])
wb_rivs   <- unique(sites_df$RivID[sites_df$is_wb])
cat(sprintf("  Reach tiers: %d statutory (RSN/SSN), %d WaterBlitz, %d co-located overlap\n",
            length(stat_rivs), length(wb_rivs), length(intersect(stat_rivs, wb_rivs))))

# ------------------------------------------------------------------------------
# 2. Export vector SVGs of fig_cokrig_cross_correlation for both runs
# ------------------------------------------------------------------------------
export_cokrig_svg <- function(r_dir, label_txt) {
  cv_d <- read.csv(file.path(r_dir, "cv_univariate_vs_cokriging.csv"), stringsAsFactors = FALSE)
  cc_d <- read.csv(file.path(r_dir, "cokriging_cross_covariance.csv"), stringsAsFactors = FALSE)
  cv_e <- cv_d[grepl("^All England", cv_d$Domain), , drop = FALSE]
  cc_e <- cc_d[grepl("^All England", cc_d$Domain), , drop = FALSE]
  configure_plot_options(format = "svg")
  out_svg <- file.path(r_dir, "fig_cokrig_cross_correlation.svg")
  plot_fig_cokrig_diagnostics(
    cross_cov_df  = cc_e,
    cv_compare_df = cv_e,
    var_names     = c("Nitrate as N", "Orthophosphate"),
    out_file      = out_svg,
    network_label = "OS Open Rivers Network (ORN v2)",
    source_label  = label_txt,
    method_label  = "Factor-LMC Joint NNGP (K=2)"
  )
  cat("  + Exported vector SVG:", out_svg, "\n")
}
export_cokrig_svg(eq_dir, "RSN + SSN + WaterBlitz [Equal Sources]")
export_cokrig_svg(wt_dir, "RSN + SSN + WaterBlitz [Prefer RSN,SSN | w_WB=0.25]")

# ------------------------------------------------------------------------------
# 3. CV, Cross-Covariance & Management Tier Comparison Table
# ------------------------------------------------------------------------------
cv_eq <- read.csv(file.path(eq_dir, "cv_univariate_vs_cokriging.csv"), stringsAsFactors = FALSE, colClasses = c(Target_Code = "character"))
cv_wt <- read.csv(file.path(wt_dir, "cv_univariate_vs_cokriging.csv"), stringsAsFactors = FALSE, colClasses = c(Target_Code = "character"))
cc_eq <- read.csv(file.path(eq_dir, "cokriging_cross_covariance.csv"), stringsAsFactors = FALSE)
cc_wt <- read.csv(file.path(wt_dir, "cokriging_cross_covariance.csv"), stringsAsFactors = FALSE)
ex_eq <- read.csv(file.path(eq_dir, "exceedance_summary_cokriging.csv"), stringsAsFactors = FALSE)
ex_wt <- read.csv(file.path(wt_dir, "exceedance_summary_cokriging.csv"), stringsAsFactors = FALSE)

m_cv <- merge(cv_eq, cv_wt, by = c("Domain", "Target", "Target_Code", "Dataset", "CV_Mode", "Method", "n_obs", "Model"),
              suffixes = c("_Equal", "_Favoured"), sort = FALSE)
m_cv$Domain_Short <- sub(" \\(CatchID=.*\\)$", "", sub(" RSN\\+SSN\\+WATERB_TOB Network$", "", m_cv$Domain))
m_cv$R2_Diff   <- round(m_cv$R2_Favoured - m_cv$R2_Equal, 4)
m_cv$RMSE_Diff <- round(m_cv$RMSE_Favoured - m_cv$RMSE_Equal, 4)
write.csv(m_cv, file.path(out_dir, "cs6_cv_comparison.csv"), row.names = FALSE)

cc_pair_eq <- cc_eq[cc_eq$Var_A != cc_eq$Var_B, ]
cc_pair_wt <- cc_wt[cc_wt$Var_A != cc_wt$Var_B, ]
m_cc <- merge(cc_pair_eq, cc_pair_wt, by = c("Domain", "Var_A", "Var_B"), suffixes = c("_Equal", "_Favoured"), sort = FALSE)
m_cc$Domain_Short <- sub(" \\(CatchID=.*\\)$", "", sub(" RSN\\+SSN\\+WATERB_TOB Network$", "", m_cc$Domain))
write.csv(m_cc, file.path(out_dir, "cs6_cross_covariance_comparison.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Reach-by-Reach Stratified Comparison & 300-DPI Signed Difference Maps
# ------------------------------------------------------------------------------
div_pal <- colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#f7f7f7", "#f4a582", "#d6604d", "#67001f"))(255)
se_div  <- colorRampPalette(c("#00441b", "#1b7837", "#a6dba0", "#f7f7f7", "#c2a5cf", "#762a83", "#40004b"))(255)

analyze_det <- function(det_code, det_label, unit_lbl) {
  p1 <- read.csv(file.path(eq_dir, sprintf("predictions_england_%s.csv", det_code)), stringsAsFactors = FALSE)
  p2 <- read.csv(file.path(wt_dir, sprintf("predictions_england_%s.csv", det_code)), stringsAsFactors = FALSE)
  ix <- match(p1$RivID, p2$RivID)
  p2 <- p2[ix, , drop = FALSE]

  p1$diff_raw      <- p2$pred_combined - p1$pred_combined
  p1$diff_model_sc <- p2$pred_model_sc - p1$pred_model_sc
  p1$diff_se_sc    <- p2$se_model_sc   - p1$se_model_sc
  p1$is_stat       <- p1$RivID %in% stat_rivs
  p1$is_wb         <- p1$RivID %in% wb_rivs
  p1$stratum <- ifelse(p1$is_stat & p1$is_wb, "Co-Located Statutory + WaterBlitz",
                ifelse(p1$is_stat, "Statutory RSN/SSN Only",
                ifelse(p1$is_wb, "WaterBlitz-Only Reaches", "Unmonitored Gap Reaches")))

  strata_list <- list(
    list(name = "All England Reaches", mask = rep(TRUE, nrow(p1))),
    list(name = "Co-Located (RSN/SSN + WaterBlitz)", mask = p1$is_stat & p1$is_wb),
    list(name = "Statutory (RSN/SSN Only)", mask = p1$is_stat & !p1$is_wb),
    list(name = "WaterBlitz-Only Citizen Reaches", mask = !p1$is_stat & p1$is_wb),
    list(name = "Unmonitored Gap Reaches", mask = !p1$is_stat & !p1$is_wb)
  )

  rows <- list()
  for (st in strata_list) {
    mk <- st$mask
    n_s <- sum(mk)
    pr_r <- suppressWarnings(cor(p1$pred_model_sc[mk], p2$pred_model_sc[mk], use = "complete.obs"))
    tier_ag <- 100.0 * mean(p1$compliance_tier[mk] == p2$compliance_tier[mk], na.rm = TRUE)
    # Prediction drift from local observation where observed (especially co-located where Equal mixed WB into RSN/SSN)
    obs_ok <- mk & (p1$has_obs == 1) & is.finite(p1$obs_model_sc) & is.finite(p2$obs_model_sc)
    rmse_eq_obs <- if (sum(obs_ok) > 0) sqrt(mean((p1$pred_model_sc[obs_ok] - p1$obs_model_sc[obs_ok])^2)) else NA_real_
    rmse_wt_obs <- if (sum(obs_ok) > 0) sqrt(mean((p2$pred_model_sc[obs_ok] - p2$obs_model_sc[obs_ok])^2)) else NA_real_

    se1_m <- mean(p1$se_model_sc[mk], na.rm = TRUE)
    se2_m <- mean(p2$se_model_sc[mk], na.rm = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      Target_Code      = det_code,
      Determinand      = det_label,
      Stratum          = st$name,
      n_reaches        = n_s,
      Equal_Med_Raw    = round(median(p1$pred_combined[mk], na.rm = TRUE), 4),
      Favoured_Med_Raw = round(median(p2$pred_combined[mk], na.rm = TRUE), 4),
      Med_Abs_Shift    = round(median(abs(p1$diff_raw[mk]), na.rm = TRUE), 4),
      Pearson_r_Scale  = round(pr_r, 4),
      Equal_Mean_SE    = round(se1_m, 4),
      Favoured_Mean_SE = round(se2_m, 4),
      SE_Delta_Pct     = round(100.0 * (se2_m - se1_m) / max(1e-8, se1_m), 2),
      Tier_Match_Pct   = round(tier_ag, 2),
      Anchor_RMSE_Eq   = round(rmse_eq_obs, 4),
      Anchor_RMSE_Fav  = round(rmse_wt_obs, 4),
      stringsAsFactors = FALSE
    )
  }

  # Render 300-DPI PNG difference map (Favoured minus Equal)
  configure_plot_options(format = "png", dpi = 300)
  out_png <- file.path(out_dir, sprintf("fig_cs6_diff_maps_%s.png", det_code))
  open_png_device(out_png, width = 2600, height = 1450, res = 300)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.6, 1.2), oma = c(0.6, 0.6, 2.8, 0.6), bg = "white")

  lim_p <- as.numeric(quantile(abs(p1$diff_model_sc), 0.985, na.rm = TRUE))
  lim_p <- max(0.15, ceiling(lim_p * 10) / 10)
  ticks_p <- seq(-lim_p, lim_p, length.out = 5)
  plot_stream_network(
    p1, "diff_model_sc", vmin = -lim_p, vmax = lim_p, pal = div_pal,
    main_title = sprintf("A. %s Prediction Shift (Favoured RSN+SSN - Equal)", det_label),
    sub_title  = sprintf("Model-scale delta across %s reaches (r = %.4f; Red = Favoured higher)",
                         format(nrow(p1), big.mark = ","), cor(p1$pred_model_sc, p2$pred_model_sc)),
    cbar_ticks = ticks_p, cbar_labels = sprintf("%+.2f", ticks_p),
    cbar_title = sprintf("Delta log(%s)", det_code),
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )

  lim_s <- as.numeric(quantile(abs(p1$diff_se_sc), 0.985, na.rm = TRUE))
  lim_s <- max(0.10, ceiling(lim_s * 10) / 10)
  ticks_s <- seq(-lim_s, lim_s, length.out = 5)
  plot_stream_network(
    p1, "diff_se_sc", vmin = -lim_s, vmax = lim_s, pal = se_div,
    main_title = sprintf("B. %s Kriging Uncertainty Shift (SE_Favoured - SE_Equal)", det_label),
    sub_title  = sprintf("Green = lower uncertainty when prioritizing RSN+SSN (mean SE %+.1f%%)",
                         100 * (mean(p2$se_model_sc) - mean(p1$se_model_sc)) / mean(p1$se_model_sc)),
    cbar_ticks = ticks_s, cbar_labels = sprintf("%+.2f", ticks_s),
    cbar_title = "Delta SE\n(log scale)",
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )
  mtext(sprintf("Case Study 6: Co-Kriging Equal Sources vs. Favouring RSN+SSN over WaterBlitz (%s, ORN v2)", det_label),
        side = 3, line = 1.1, outer = TRUE, cex = 1.05, font = 2, col = "#1b2631")
  dev.off()
  cat("  + Rendered 300-DPI diff map:", out_png, "\n")

  rm(p1, p2); gc(FALSE)
  do.call(rbind, rows)
}

reach_sum <- rbind(
  analyze_det("0117", "Nitrate as N", "mg/L"),
  analyze_det("0180", "Orthophosphate", "mg/L")
)
write.csv(reach_sum, file.path(out_dir, "cs6_reach_comparison_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Render Vector SVG Comparison Chart (CV, Cross-Correlation & Tier Uncertainty)
# ------------------------------------------------------------------------------
out_svg <- file.path(out_dir, "fig_cs6_equal_vs_favoured_comparison.svg")
configure_plot_options(format = "svg")
open_png_device(out_svg, width = 2000, height = 1150, res = 150)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.2, 1.5), oma = c(0.8, 0.8, 2.8, 0.8), bg = "white")

# Panel A: England Mean Kriging SE
se_eq <- c(reach_sum$Equal_Mean_SE[1], reach_sum$Equal_Mean_SE[6])
se_wt <- c(reach_sum$Favoured_Mean_SE[1], reach_sum$Favoured_Mean_SE[6])
mat_se <- rbind(se_eq, se_wt)
bp1 <- barplot(mat_se, beside = TRUE, names.arg = c("Nitrate (0117)", "Orthophosphate (0180)"),
               col = c("#3b82f6", "#10b981"), ylim = c(0, max(mat_se) * 1.25),
               ylab = "Mean Network Kriging SE (log scale)",
               main = "A. National Uncertainty (Mean Kriging SE)", border = NA)
legend("topleft", legend = c("Equal Sources (w_all = 1)", "Favoured RSN+SSN (w_WB = 0.25)"),
       fill = c("#3b82f6", "#10b981"), bty = "n", cex = 0.88)
text(bp1, mat_se + max(mat_se) * 0.04, labels = sprintf("%.3f", mat_se), cex = 0.85, font = 2)

# Panel B: Uncertain Tier Reaches (%) across Domains for Orthophosphate (0180)
ex0_eq <- ex_eq[ex_eq$Target == "Orthophosphate" & ex_eq$Tier == "Uncertain / Priority Monitoring", ]
ex0_wt <- ex_wt[ex_wt$Target == "Orthophosphate" & ex_wt$Tier == "Uncertain / Priority Monitoring", ]
ex0_m  <- merge(ex0_eq, ex0_wt, by = c("Domain", "Target", "Tier"), suffixes = c("_eq", "_wt"), sort = FALSE)
dom_lbl <- sub(" \\(CatchID=.*\\)$", "", sub(" RSN\\+SSN\\+WATERB_TOB Network$", "", ex0_m$Domain))
mat_unc <- rbind(ex0_m$Reach_Pct_eq, ex0_m$Reach_Pct_wt)
bp2 <- barplot(mat_unc, beside = TRUE, names.arg = dom_lbl, col = c("#f59e0b", "#10b981"),
               ylim = c(0, 82), ylab = "Uncertain Tier Reaches (%)",
               main = "B. Orthophosphate (0.069 mg/L) Tier Uncertainty Reduction", border = NA)
legend("topleft", legend = c("Equal Sources", "Favoured RSN+SSN"), fill = c("#f59e0b", "#10b981"), bty = "n", cex = 0.88)
text(bp2, mat_unc + 2.5, labels = sprintf("%.1f%%", mat_unc), cex = 0.82, font = 2)

# Panel C: Cross-Correlation Decomposition (Nitrate <-> Orthophosphate, England)
ce <- cc_pair_eq[grepl("^All England", cc_pair_eq$Domain), ][1, ]
cw <- cc_pair_wt[grepl("^All England", cc_pair_wt$Domain), ][1, ]
mat_cc <- rbind(c(ce$Corr_TU, ce$Corr_TD, ce$Corr_Nugget, ce$Corr_Total),
                c(cw$Corr_TU, cw$Corr_TD, cw$Corr_Nugget, cw$Corr_Total))
bp3 <- barplot(mat_cc, beside = TRUE, names.arg = c("Tail-Up", "Tail-Down", "Nugget", "Total Corr"),
               col = c("#6366f1", "#06b6d4"), ylim = c(0, 0.85),
               ylab = "Cross-Correlation (r)",
               main = "C. Learned NO3-PO4 Cross-Correlation Structure", border = NA)
legend("topleft", legend = c("Equal Sources", "Favoured RSN+SSN"), fill = c("#6366f1", "#06b6d4"), bty = "n", cex = 0.88)
text(bp3, mat_cc + 0.03, labels = sprintf("%+.2f", mat_cc), cex = 0.82, font = 2)

# Panel D: Aux-Aided Co-Kriging CV R2 Across Basins (Nitrate 0117 & Phosphate 0180)
ax_eq <- cv_eq[cv_eq$Model == "Co-Kriging (Aux-Aided)", ]
ax_wt <- cv_wt[cv_wt$Model == "Co-Kriging (Aux-Aided)", ]
ax_m  <- merge(ax_eq, ax_wt, by = c("Domain", "Target_Code"), suffixes = c("_eq", "_wt"), sort = FALSE)
ax_m$Tag <- paste0(sub(" \\(CatchID=.*\\)$", "", sub(" RSN\\+SSN\\+WATERB_TOB Network$", "", ax_m$Domain)), "\n", ax_m$Target_Code)
mat_r2 <- rbind(ax_m$R2_eq, ax_m$R2_wt)
bp4 <- barplot(mat_r2, beside = TRUE, names.arg = ax_m$Tag, col = c("#8b5cf6", "#ec4899"),
               ylim = c(0, max(mat_r2, na.rm = TRUE) * 1.25), cex.names = 0.72,
               ylab = "Aux-Aided CV R2 (Unweighted Holdout)",
               main = "D. Aux-Aided R2 Against Mixed Citizen+Statutory CV Folds", border = NA)
legend("topright", legend = c("Equal Sources", "Favoured RSN+SSN"), fill = c("#8b5cf6", "#ec4899"), bty = "n", cex = 0.85)
text(bp4, mat_r2 + 0.02, labels = sprintf("%.2f", mat_r2), cex = 0.72, font = 2)

mtext("RSN + SSN + WaterBlitz Co-Kriging: Equal Sources vs. Prioritizing Statutory RSN+SSN",
      side = 3, line = 0.8, outer = TRUE, cex = 1.15, font = 2, col = "#1b2631")
dev.off()
cat("  + Exported comparison SVG:", out_svg, "\n")
print(reach_sum)
