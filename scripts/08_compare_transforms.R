#!/usr/bin/env Rscript
# ==============================================================================
# 08_compare_transforms.R
#
# Compares S3N Nitrate (0117) models across 4 response transformations:
#   1. Raw / Identity (--transform=none)
#   2. Natural Logarithm (--transform=log)
#   3. Square Root (--transform=sqrt)
#   4. Profile-MLE Box-Cox (--transform=boxcox, lambda = 0.50)
#
# Evaluates both model-scale and back-transformed raw-scale (mg/L) 5-fold CV
# accuracy, residual normality/skewness, negative-prediction rates across all
# 90,876 ORN v2 England reaches, and exports high-DPI diff maps + vector SVG
# charts for the Modelling Dashboard.
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
dpi_arg <- grep("^--dpi=", args, value = TRUE)
dpi_val <- if (length(dpi_arg) > 0L) as.numeric(sub("^--dpi=", "", dpi_arg[1])) else 300
if (!is.finite(dpi_val) || dpi_val < 72) dpi_val <- 300

comp_dir <- file.path(BASE_DIR, "output/comparisons/cs7_transforms_0117")
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

runs <- list(
  list(key = "none",   label = "Raw (none)",     lambda = NA_real_, dir = file.path(BASE_DIR, "output/results/det_0117_orn_rsn_trans_raw")),
  list(key = "log",    label = "Log (ln)",       lambda = NA_real_, dir = file.path(BASE_DIR, "output/results/det_0117_orn_rsn")),
  list(key = "sqrt",   label = "Sqrt",           lambda = NA_real_, dir = file.path(BASE_DIR, "output/results/det_0117_orn_rsn_trans_sqrt")),
  list(key = "boxcox", label = "Box-Cox (0.50)", lambda = 0.50,     dir = file.path(BASE_DIR, "output/results/det_0117_orn_rsn_trans_boxcox_lam0.50"))
)

calc_skew <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  if (!is.finite(s) || s <= 0) return(0)
  mean(((x - m) / s)^3)
}

rows <- list()
pred_list <- list()
cv_raw_list <- list()
cv_resid_list <- list()

for (i in seq_along(runs)) {
  r <- runs[[i]]
  tr <- build_transformer(r$key, lambda = r$lambda)
  p_df  <- fast_read_csv(file.path(r$dir, "predictions_0117.csv"))
  cv_df <- fast_read_csv(file.path(r$dir, "cv_comparison_0117.csv"))
  cv_p  <- fast_read_csv(file.path(r$dir, "cv_predictions_england_0117.csv"))
  ex_f  <- file.path(r$dir, "exceedance_summary_0117.csv")
  ex_df <- if (file.exists(ex_f)) fast_read_csv(ex_f) else data.frame(Domain = character(0), Tier = character(0), Reach_Pct = numeric(0))

  p_df <- p_df[order(p_df$RivID), ]
  pred_list[[r$key]] <- p_df

  # Merge raw reach-mean observation value (obs_value) on 1210 CV stations
  obs_m <- p_df[p_df$has_obs == 1, c("RivID", "obs_value")]
  m_cv  <- merge(cv_p, obs_m, by = "RivID", sort = FALSE)
  y_obs_raw  <- m_cv$obs_value
  y_pred_raw <- tr$inv(m_cv$cv_comb)
  # Ensure Box-Cox / sqrt never yield NaN; leave Raw untransformed to show negative values honestly
  err_raw <- y_pred_raw - y_obs_raw
  rmse_raw <- sqrt(mean(err_raw^2, na.rm = TRUE))
  mae_raw  <- mean(abs(err_raw), na.rm = TRUE)
  bias_raw <- mean(err_raw, na.rm = TRUE)
  med_ae_raw <- median(abs(err_raw), na.rm = TRUE)
  r2_raw <- 1.0 - sum(err_raw^2, na.rm = TRUE) / sum((y_obs_raw - mean(y_obs_raw, na.rm = TRUE))^2, na.rm = TRUE)
  neg_cv <- sum(y_pred_raw < 0, na.rm = TRUE)

  resid_sc <- m_cv$obs_y - m_cv$cv_comb
  resid_skew <- calc_skew(resid_sc)
  obs_sc_skew <- calc_skew(m_cv$obs_y)

  eng_comb <- cv_df[grepl("^All England", cv_df$Domain) & grepl("Combined", cv_df$Model), ][1, ]
  eng_ols  <- cv_df[grepl("^All England", cv_df$Domain) & grepl("OLS", cv_df$Model), ][1, ]

  ex_eng <- ex_df[grepl("^All England", ex_df$Domain), ]
  get_pct <- function(tiers) sum(ex_eng$Reach_Pct[ex_eng$Tier %in% tiers], na.rm = TRUE)

  neg_n <- sum(p_df$pred_combined < 0, na.rm = TRUE)
  q <- quantile(p_df$pred_combined, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)

  cv_raw_list[[r$key]] <- data.frame(obs_raw = y_obs_raw, pred_raw = y_pred_raw, resid_sc = resid_sc)
  cv_resid_list[[r$key]] <- resid_sc

  rows[[i]] <- data.frame(
    Transform = r$label,
    Key = r$key,
    Lambda = if (is.na(r$lambda)) ifelse(r$key == "sqrt", "0.50 (fixed)", "—") else "0.50 (MLE)",
    R2_ModelSc = round(eng_comb$R2, 4),
    Delta_R2_vs_OLS = round(eng_comb$R2 - eng_ols$R2, 4),
    RMSE_ModelSc = round(eng_comb$RMSE, 4),
    R2_Raw_mgL = round(r2_raw, 4),
    RMSE_Raw_mgL = round(rmse_raw, 4),
    MAE_Raw_mgL = round(mae_raw, 4),
    MedAE_Raw_mgL = round(med_ae_raw, 4),
    Bias_Raw_mgL = round(bias_raw, 4),
    CV_Neg_Count = neg_cv,
    ModelScale_Obs_Skew = round(obs_sc_skew, 3),
    CV_Resid_Skew = round(resid_skew, 3),
    Eng_Min_mgL = round(min(p_df$pred_combined, na.rm = TRUE), 3),
    Eng_P05_mgL = round(as.numeric(q[1]), 3),
    Eng_Med_mgL = round(as.numeric(q[2]), 3),
    Eng_Mean_mgL = round(mean(p_df$pred_combined, na.rm = TRUE), 3),
    Eng_P95_mgL = round(as.numeric(q[3]), 3),
    Eng_Max_mgL = round(max(p_df$pred_combined, na.rm = TRUE), 3),
    Neg_Reaches = neg_n,
    Neg_Reach_Pct = round(100.0 * neg_n / nrow(p_df), 2),
    Mean_Krig_SE_Sc = round(mean(p_df$se_model_sc, na.rm = TRUE), 4),
    Compliant_Pct = round(get_pct(c("Confidently Compliant", "Likely Compliant")), 2),
    Uncertain_Pct = round(get_pct(c("Uncertain / Priority Monitoring")), 2),
    Exceeding_Pct = round(get_pct(c("Likely Exceeding", "Confidently Exceeding")), 2),
    stringsAsFactors = FALSE
  )
}

sum_df <- do.call(rbind, rows)

# Spearman rank correlation across all reaches vs Log (ln) aligned by RivID
log_df <- pred_list[["log"]]
sum_df$Spearman_vs_Log <- sapply(names(pred_list), function(k) {
  pk <- pred_list[[k]]
  ix_log <- match(pk$RivID, log_df$RivID)
  round(cor(pk$pred_combined, log_df$pred_combined[ix_log], method = "spearman", use = "complete.obs"), 4)
})
fast_write_csv(sum_df, file.path(comp_dir, "cs7_transform_summary.csv"))
print(sum_df)

# ------------------------------------------------------------------------------
# Render Vector SVG Comparison Chart (CV R2, Raw RMSE/MAE, Skewness & Dist)
# ------------------------------------------------------------------------------
svg_path <- file.path(comp_dir, "fig_cs7_transform_comparison.svg")
configure_plot_options(format = "svg")
open_png_device(svg_path, width = 2000, height = 1180, res = 150)
tryCatch({
  par(mfrow = c(2, 2), mar = c(4.5, 4.6, 3.2, 1.5), oma = c(0.8, 0.8, 2.8, 0.8), bg = "white")

lbls <- sum_df$Transform
cols4 <- c("#ef4444", "#2563eb", "#10b981", "#8b5cf6")

# Panel A: Model-Scale CV R2 vs Back-Transformed Raw-Scale (mg/L) CV R2
mat_r2 <- rbind(sum_df$R2_ModelSc, sum_df$R2_Raw_mgL)
bp1 <- barplot(mat_r2, beside = TRUE, names.arg = lbls, col = c("#2563eb", "#f59e0b"),
               ylim = c(0, 0.82), ylab = "5-Fold CV Variance Explained (R2)",
               main = "A. Model-Scale R2 vs. Back-Transformed Raw-Scale (mg/L) R2", border = NA)
legend("topright", legend = c("Model Scale (transformed)", "Back-Transformed Raw (mg/L)"),
       fill = c("#2563eb", "#f59e0b"), bty = "n", cex = 0.88)
text(bp1, mat_r2 + 0.025, labels = sprintf("%.3f", mat_r2), cex = 0.82, font = 2)

# Panel B: Back-Transformed Raw-Scale Errors in Comparable mg/L Units
mat_err <- rbind(sum_df$RMSE_Raw_mgL, sum_df$MAE_Raw_mgL, sum_df$MedAE_Raw_mgL)
bp2 <- barplot(mat_err, beside = TRUE, names.arg = lbls, col = c("#ef4444", "#06b6d4", "#10b981"),
               ylim = c(0, max(mat_err) * 1.28), ylab = "Raw Nitrate Error (mg/L as N)",
               main = "B. Out-of-Sample CV Error on Comparable Raw Scale (mg/L)", border = NA)
legend("topright", legend = c("RMSE (mg/L)", "MAE (mg/L)", "Median |Err| (mg/L)"),
       fill = c("#ef4444", "#06b6d4", "#10b981"), bty = "n", cex = 0.85)
text(bp2, mat_err + 0.08, labels = sprintf("%.2f", mat_err), cex = 0.76, font = 2)

# Panel C: Model-Scale Observed Skewness vs CV Residual Skewness
mat_sk <- rbind(sum_df$ModelScale_Obs_Skew, sum_df$CV_Resid_Skew)
bp3 <- barplot(mat_sk, beside = TRUE, names.arg = lbls, col = c("#6366f1", "#ec4899"),
               ylim = c(min(mat_sk, 0) - 0.35, max(mat_sk) + 0.45),
               ylab = "Sample Skewness (0 = Symmetric)",
               main = "C. Distribution & Residual Gaussianity (Skewness)", border = NA)
abline(h = 0, lty = 2, col = "#64748b", lwd = 1.2)
legend("topright", legend = c("Station Response Skewness", "5-Fold CV Residual Skewness"),
       fill = c("#6366f1", "#ec4899"), bty = "n", cex = 0.85)
text(bp3, mat_sk + ifelse(mat_sk >= 0, 0.09, -0.14), labels = sprintf("%+.2f", mat_sk), cex = 0.80, font = 2)

  # Panel D: National P05 / Median / Mean / P95 vs Physical Positivity Bound
  p05_raw_val <- sum_df$Eng_P05_mgL[sum_df$Key == "none"][1]
  ylim_d_lo <- min(-0.5, min(sum_df$Eng_P05_mgL, na.rm = TRUE) * 1.25)
  ylim_d_hi <- max(1.0, max(sum_df$Eng_P95_mgL, na.rm = TRUE) * 1.22)
  bp4 <- barplot(rbind(sum_df$Eng_P05_mgL, sum_df$Eng_Med_mgL, sum_df$Eng_Mean_mgL, sum_df$Eng_P95_mgL),
                 beside = TRUE, names.arg = lbls, col = c("#94a3b8", "#3b82f6", "#0ea5e9", "#f97316"),
                 ylim = c(ylim_d_lo, ylim_d_hi),
                 ylab = "Predicted Nitrate (mg/L as N)",
                 main = "D. National Reach Quantiles & Physical Positivity Bound", border = NA)
  abline(h = 0, lty = 2, col = "#dc2626", lwd = 1.4)
  legend("topleft", legend = c(sprintf("P05 (Raw=%.2f)", p05_raw_val), "Median", "Mean", "P95"),
         fill = c("#94a3b8", "#3b82f6", "#0ea5e9", "#f97316"), bty = "n", cex = 0.84, ncol = 2)
  text(colMeans(bp4), sum_df$Eng_P95_mgL + 0.05 * ylim_d_hi,
       labels = ifelse(sum_df$Neg_Reaches > 0,
                       sprintf("%s neg\n(%.1f%%)", format(sum_df$Neg_Reaches, big.mark=","), sum_df$Neg_Reach_Pct),
                       "0 negative\n(100% >= 0)"),
       cex = 0.78, font = 2, col = ifelse(sum_df$Neg_Reaches > 0, "#dc2626", "#059669"))

  mtext("Nitrate (0117) Transformation Comparison: Untransformed (Raw) vs. Log (ln) vs. Sqrt vs. Box-Cox (MLE = 0.50)",
        outer = TRUE, cex = 1.15, font = 2, line = 0.8)
}, finally = {
  if (dev.cur() > 1L) dev.off()
})
cat("  + Exported transform SVG:", svg_path, "\n")

# ------------------------------------------------------------------------------
# Render 300 DPI Two-Panel National Spatial Map (Raw Negatives & Raw - Log Diff)
# ------------------------------------------------------------------------------
configure_plot_options(format = "png", dpi = dpi_val)
p_raw <- pred_list[["none"]]
p_log <- pred_list[["log"]]

diff_df <- p_raw
ix_log_map <- match(p_raw$RivID, p_log$RivID)
diff_df$raw_minus_log <- p_raw$pred_combined - p_log$pred_combined[ix_log_map]
neg_mask <- p_raw$pred_combined < 0

map_png <- file.path(comp_dir, "fig_cs7_transform_diff_0117.png")
open_png_device(map_png, width = 2600, height = 1400, res = dpi_val)
tryCatch({
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.8, 1.2), oma = c(0.6, 0.6, 2.4, 0.6), bg = "white")

  raw_pal <- colorRampPalette(c("#7e22ce", "#e2e8f0", "#38bdf8", "#facc15", "#ef4444"))(256)
  p_raw$pred_clip <- pmax(-1.6, pmin(10.0, p_raw$pred_combined))
  plot_stream_network(
    p_raw, "pred_clip", vmin = -1.6, vmax = 10.0, pal = raw_pal,
    main_title = "A. Untransformed Raw (none) S3N Nitrate (mg/L)",
    sub_title  = sprintf("%s of %s reaches (%.1f%%) predict < 0 mg/L (purple; min = %.2f)",
                         format(sum(neg_mask, na.rm = TRUE), big.mark=","),
                         format(nrow(p_raw), big.mark=","),
                         100*mean(neg_mask, na.rm = TRUE), min(p_raw$pred_combined, na.rm = TRUE)),
    cbar_ticks  = c(-1.6, 0, 2.5, 5.0, 7.5, 10.0),
    cbar_labels = c("<= -1.6", "0.0", "2.5", "5.0", "7.5", ">= 10"),
    cbar_title  = "Raw NO3-N\n(mg/L)", show_obs = FALSE, base_lwd = 0.24, max_lwd = 2.0
  )

  diff_pal <- colorRampPalette(c("#1d4ed8", "#93c5fd", "#f8fafc", "#fca5a5", "#b91c1c"))(256)
  lim_d <- 3.0
  diff_df$raw_minus_log_clip <- pmax(-lim_d, pmin(lim_d, diff_df$raw_minus_log))
  plot_stream_network(
    diff_df, "raw_minus_log_clip", vmin = -lim_d, vmax = lim_d, pal = diff_pal,
    main_title = "B. Raw (none) minus Log (ln) Predicted Difference",
    sub_title  = sprintf("Red = Raw higher (Chalk/arithmetic mean shift); Blue = Raw lower (median diff %+.2f mg/L)",
                         median(diff_df$raw_minus_log, na.rm = TRUE)),
    cbar_ticks  = c(-3.0, -1.5, 0.0, 1.5, 3.0),
    cbar_labels = c("<= -3.0", "-1.5", "0.0", "+1.5", ">= +3.0"),
    cbar_title  = "Raw - Log\n(mg/L)", show_obs = FALSE, base_lwd = 0.24, max_lwd = 2.0
  )

  mtext(sprintf("Spatial Consequence of Transformation Choice on ORN v2 (%s England Reaches)",
                format(nrow(p_raw), big.mark = ",")),
        outer = TRUE, cex = 1.08, font = 2, line = 0.6)
}, finally = {
  if (dev.cur() > 1L) dev.off()
})
cat("  + Exported 300-DPI diff map:", map_png, "\n")
