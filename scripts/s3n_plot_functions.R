# ==============================================================================
# s3n_plot_functions.R
#
# Modular base-R plotting library for Scalable Spatial Stream Network (S3N)
# models. Sourced by both:
#   - scripts/run_s3n_determinand.R
#   - scripts/03_generate_plots.R
#
# Provides:
#   - Color palettes & scale helpers: val_cols, uncert_cols, map_to_colors(),
#     draw_colorbar(), compute_color_scale(), compute_se_scale()
#   - Stream network renderer: plot_stream_network()
#   - Schema normalizers: normalize_pred_df(), normalize_cv_pts(), normalize_cv_tab()
#   - Publication figure generators:
#       * plot_fig1_england_map()
#       * plot_basin_map()
#       * plot_fig2_humber_map()
#       * plot_fig3_thames_severn_map()
#       * plot_fig4_diagnostics()
# ==============================================================================

val_cols     <- colorRampPalette(c("#1b4f72", "#2874a6", "#1abc9c", "#f4d03f", "#e67e22", "#c0392b", "#641e16"))(256)
val_cols_rev <- rev(val_cols)
uncert_cols  <- colorRampPalette(c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(256)

get_val_pal <- function(reverse_color = FALSE) {
  if (isTRUE(reverse_color)) val_cols_rev else val_cols
}

map_to_colors <- function(vals, vmin, vmax, pal) {
  u <- pmin(1, pmax(0, (vals - vmin) / max(1e-9, vmax - vmin)))
  idx <- pmax(1, pmin(length(pal), round(u * (length(pal) - 1)) + 1))
  pal[idx]
}

draw_colorbar <- function(x0, y0, x1, y1, vmin, vmax, pal, ticks, tick_labels, title_txt) {
  n <- length(pal)
  ys <- seq(y0, y1, length.out = n + 1)
  for (i in seq_len(n)) {
    rect(x0, ys[i], x1, ys[i + 1], col = pal[i], border = NA)
  }
  rect(x0, y0, x1, y1, border = "#333333", lwd = 0.8)
  for (j in seq_along(ticks)) {
    ty <- y0 + (ticks[j] - vmin) / max(1e-9, vmax - vmin) * (y1 - y0)
    if (ty >= y0 - 0.02 * (y1 - y0) && ty <= y1 + 0.02 * (y1 - y0)) {
      segments(x1, ty, x1 + (x1 - x0) * 0.35, ty, col = "#333333", lwd = 0.8)
      text(x1 + (x1 - x0) * 0.45, ty, labels = tick_labels[j], adj = c(0, 0.5), cex = 0.78, col = "#222222")
    }
  }
  text(x0, y1 + (y1 - y0) * 0.08, labels = title_txt, adj = c(0, 0), cex = 0.82, font = 2, col = "#111111")
}

compute_color_scale <- function(vals_model_sc, use_log, vmin_raw = NULL, vmax_raw = NULL) {
  qs <- quantile(vals_model_sc, probs = c(0.01, 0.99), na.rm = TRUE)
  vmin <- as.numeric(qs[1])
  vmax <- as.numeric(qs[2])
  if (!is.null(vmin_raw) && !is.na(vmin_raw)) {
    vmin <- if (use_log) log(max(1e-9, as.numeric(vmin_raw))) else as.numeric(vmin_raw)
  }
  if (!is.null(vmax_raw) && !is.na(vmax_raw)) {
    vmax <- if (use_log) log(max(1e-9, as.numeric(vmax_raw))) else as.numeric(vmax_raw)
  }
  if (vmax - vmin < 1e-4) {
    vmin <- vmin - 0.5
    vmax <- vmax + 0.5
  }
  if (use_log) {
    cands <- c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5,
               1, 2, 2.5, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
    raw_in <- cands[cands >= exp(vmin) * 0.95 & cands <= exp(vmax) * 1.05]
    if (length(raw_in) > 7) {
      raw_in <- raw_in[unique(round(seq(1, length(raw_in), length.out = 6)))]
    }
    if (length(raw_in) >= 3) {
      ticks <- log(raw_in)
      labels <- formatC(raw_in, format = "fg", digits = 2)
      if (is.null(vmin_raw) || is.na(vmin_raw)) vmin <- min(vmin, min(ticks))
      if (is.null(vmax_raw) || is.na(vmax_raw)) vmax <- max(vmax, max(ticks))
    } else {
      ticks <- seq(vmin, vmax, length.out = 6)
      labels <- formatC(exp(ticks), format = "fg", digits = 2)
    }
  } else {
    ticks <- pretty(c(vmin, vmax), n = 6)
    ticks <- ticks[ticks >= vmin & ticks <= vmax]
    if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
    labels <- formatC(ticks, format = "fg", digits = 3)
  }
  list(vmin = vmin, vmax = vmax, ticks = ticks, labels = labels)
}

compute_se_scale <- function(se_vals, vmin_se = NULL, vmax_se = NULL) {
  qs <- quantile(se_vals, probs = c(0.02, 0.98), na.rm = TRUE)
  vmin <- max(0.01, floor(qs[1] * 20) / 20)
  vmax <- max(vmin + 0.15, ceiling(qs[2] * 20) / 20)
  if (!is.null(vmin_se) && !is.na(vmin_se)) vmin <- max(0.001, as.numeric(vmin_se))
  if (!is.null(vmax_se) && !is.na(vmax_se)) vmax <- max(vmin + 0.05, as.numeric(vmax_se))
  ticks <- pretty(c(vmin, vmax), n = 5)
  ticks <- ticks[ticks >= vmin & ticks <= vmax]
  if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
  list(vmin = vmin, vmax = vmax, ticks = ticks, labels = sprintf("%.2f", ticks))
}

# Schema normalizers so both 02_fit_s3n_models.R outputs and run_s3n_determinand.R outputs work seamlessly
normalize_pred_df <- function(df) {
  if (!("pred_model_sc" %in% names(df)) && ("pred_log_combined" %in% names(df))) {
    df$pred_model_sc <- df$pred_log_combined
  }
  if (!("se_model_sc" %in% names(df)) && ("se_log_combined" %in% names(df))) {
    df$se_model_sc <- df$se_log_combined
  }
  if (!("obs_model_sc" %in% names(df)) && ("obs_log_nitrate" %in% names(df))) {
    df$obs_model_sc <- df$obs_log_nitrate
  }
  if (!("pred_tu_sc" %in% names(df)) && ("pred_log_tu" %in% names(df))) {
    df$pred_tu_sc <- df$pred_log_tu
  }
  if (!("pred_td_sc" %in% names(df)) && ("pred_log_td" %in% names(df))) {
    df$pred_td_sc <- df$pred_log_td
  }
  if (!("pred_combined" %in% names(df)) && ("pred_mgL_combined" %in% names(df))) {
    df$pred_combined <- df$pred_mgL_combined
  }
  df
}

normalize_cv_pts <- function(cv_df) {
  if (is.null(cv_df)) return(NULL)
  if (!("obs_y" %in% names(cv_df)) && ("obs_log" %in% names(cv_df))) {
    cv_df$obs_y <- cv_df$obs_log
  }
  cv_df
}

normalize_cv_tab <- function(cv_tab) {
  if (is.null(cv_tab)) return(NULL)
  if (!("R2" %in% names(cv_tab)) && ("R2_log" %in% names(cv_tab))) {
    cv_tab$R2 <- cv_tab$R2_log
  }
  if (!("RMSE" %in% names(cv_tab)) && ("RMSE_log" %in% names(cv_tab))) {
    cv_tab$RMSE <- cv_tab$RMSE_log
  }
  if (!("MAE" %in% names(cv_tab)) && ("MAE_log" %in% names(cv_tab))) {
    cv_tab$MAE <- cv_tab$MAE_log
  }
  cv_tab
}

plot_stream_network <- function(df, val_col, vmin, vmax, pal, main_title, sub_title,
                                cbar_ticks, cbar_labels, cbar_title,
                                show_obs = FALSE, obs_col = NULL, base_lwd = 0.35, max_lwd = 2.6) {
  ord <- order(df$Shreve)
  d   <- df[ord, ]

  xr <- range(c(d$start_X, d$end_X), na.rm = TRUE)
  yr <- range(c(d$start_Y, d$end_Y), na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  plot(NA, xlim = c(xr[1] - 0.03 * dx, xr[2] + 0.16 * dx),
       ylim = c(yr[1] - 0.03 * dy, yr[2] + 0.06 * dy),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#d0d4d9", lwd = 1)

  lwds <- base_lwd + (max_lwd - base_lwd) * pmin(1, log10(pmax(1, d$Shreve)) / log10(max(10, max(d$Shreve))))
  cols <- map_to_colors(d[[val_col]], vmin, vmax, pal)
  segments(d$start_X, d$start_Y, d$end_X, d$end_Y, col = cols, lwd = lwds)

  if (show_obs && !is.null(obs_col)) {
    obs_d <- d[d$has_obs == 1, ]
    obs_c <- map_to_colors(obs_d[[obs_col]], vmin, vmax, pal)
    points(obs_d$mid_X, obs_d$mid_Y, pch = 21, bg = obs_c, col = "#111111", cex = 0.85, lwd = 0.7)
  }

  # Scale bar
  sb_len <- if (dx > 250000) 100000 else if (dx > 80000) 25000 else 15000
  sb_lab <- paste0(sb_len / 1000, " km")
  sb_x0  <- if (dx > 250000) xr[2] - sb_len - 0.02 * dx else xr[1] + 0.02 * dx
  sb_y0  <- yr[1] + 0.015 * dy
  rect(sb_x0 - 0.01 * dx, sb_y0 - 0.015 * dy, sb_x0 + sb_len + 0.01 * dx, sb_y0 + 0.038 * dy,
       col = rgb(1, 1, 1, 0.85), border = NA)
  segments(sb_x0, sb_y0, sb_x0 + sb_len, sb_y0, lwd = 2.5, col = "#222222")
  segments(sb_x0, sb_y0 - 0.008 * dy, sb_x0, sb_y0 + 0.008 * dy, lwd = 1.5, col = "#222222")
  segments(sb_x0 + sb_len, sb_y0 - 0.008 * dy, sb_x0 + sb_len, sb_y0 + 0.008 * dy, lwd = 1.5, col = "#222222")
  text(sb_x0 + 0.5 * sb_len, sb_y0 + 0.020 * dy, labels = sb_lab, cex = 0.78, font = 2, col = "#222222")

  # Colorbar
  cb_x0 <- xr[2] + 0.02 * dx
  cb_x1 <- xr[2] + 0.055 * dx
  cb_y0 <- yr[1] + 0.12 * dy
  cb_y1 <- yr[1] + 0.72 * dy
  draw_colorbar(cb_x0, cb_y0, cb_x1, cb_y1, vmin, vmax, pal, cbar_ticks, cbar_labels, cbar_title)

  mtext(main_title, side = 3, line = 1.2, adj = 0, cex = 0.96, font = 2, col = "#111111")
  mtext(sub_title,  side = 3, line = 0.15, adj = 0, cex = 0.80, col = "#444444")
}

shorten_det_name <- function(det_name) {
  tag <- if (grepl("\\[.*\\]$", det_name)) sub("^.*(\\[.*\\])$", " \\1", det_name) else ""
  base <- sub(" \\[.*\\]$", "", det_name)
  if (grepl("\\{DOC\\}", base)) return(paste0("DOC", tag))
  if (grepl("^Oxygen, Dissolved", base)) return(paste0("Dissolved Oxygen", tag))
  if (grepl("^Nitrogen, Total Oxidised", base)) return(paste0("Total Oxidised N", tag))
  if (grepl("^Nitrogen, Total", base)) return(paste0("Total Nitrogen", tag))
  if (grepl("^Phosphorus, Total", base)) return(paste0("Total Phosphorus", tag))
  s <- sub(",.*| :-.*| to pH.*| at 25.*", "", base)
  paste0(substr(s, 1, 22), tag)
}

# Figure 1: Full England Stream Network Map
plot_fig1_england_map <- function(eng_pred, det_code, det_name, det_unit, use_log, out_file,
                                  reverse_color = FALSE, sc_val = NULL, sc_se = NULL) {
  eng_pred <- normalize_pred_df(eng_pred)
  pal_use  <- get_val_pal(reverse_color)
  png(out_file, width = 2600, height = 1450, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2), bg = "white")

  if (is.null(sc_val)) sc_val <- compute_color_scale(eng_pred$pred_model_sc, use_log)
  if (is.null(sc_se))  sc_se  <- compute_se_scale(eng_pred$se_model_sc)
  n_obs_eng <- sum(eng_pred$has_obs == 1)
  n_tot_eng <- nrow(eng_pred)
  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, det_unit)

  plot_stream_network(
    eng_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. S3N Predicted %s Across England (%s Reaches)", det_short, format(n_tot_eng, big.mark = ",")),
    sub_title  = sprintf("Scalable Spatial Stream Network (S3N) model across %d monitored catchments",
                         length(unique(eng_pred$CatchID))),
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )

  plot_stream_network(
    eng_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. S3N Prediction Uncertainty & RSN Stations (n = %s)", format(n_obs_eng, big.mark = ",")),
    sub_title  = sprintf("Kriging standard error (%s) with observed RSN stations overlaid",
                         if (use_log) "log scale" else det_unit),
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = if (use_log) "Kriging SE\n(log scale)" else sprintf("Kriging SE\n(%s)", det_unit),
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )
  obs_eng <- eng_pred[eng_pred$has_obs == 1, ]
  points(obs_eng$mid_X, obs_eng$mid_Y, pch = 21,
         bg = map_to_colors(obs_eng$obs_model_sc, sc_val$vmin, sc_val$vmax, pal_use),
         col = "#111111", cex = 0.65, lwd = 0.5)
  legend("topleft", legend = sprintf("Observed RSN Station (colored by %s)", det_short),
         pch = 21, pt.bg = "#e67e22", col = "#111111", pt.cex = 1.0, bty = "n", cex = 0.80)
}

# Generic 2-Panel River Basin High-Resolution S3N Map (used for Figure 2 Humber & Figure 5+ Custom Basins)
plot_basin_map <- function(basin_pred, basin_title, domain_label, cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                           reverse_color = FALSE, sc_se = NULL) {
  basin_pred <- normalize_pred_df(basin_pred)
  cv_tab     <- normalize_cv_tab(cv_tab)
  pal_use    <- get_val_pal(reverse_color)
  png(out_file, width = 2600, height = 1350, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2), bg = "white")

  if (is.null(sc_se)) sc_se <- compute_se_scale(basin_pred$se_model_sc)
  n_obs_basin <- sum(basin_pred$has_obs == 1)
  det_short   <- shorten_det_name(det_name)
  cv_row      <- subset(cv_tab, Domain == domain_label & Model == "S3N Combined (Tail-Up + Tail-Down)")
  sub_a <- if (nrow(cv_row) > 0) {
    sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | Circles show %d observed RSN stations",
            cv_row$R2[1], cv_row$RMSE[1], n_obs_basin)
  } else {
    sprintf("England S3N model subset | Circles show %d observed RSN stations in %s", n_obs_basin, basin_title)
  }

  plot_stream_network(
    basin_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. %s: S3N Predicted %s (%s Reaches)",
                         basin_title, det_short, format(nrow(basin_pred), big.mark = ",")),
    sub_title  = sub_a,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels,
    cbar_title = sprintf("%s\n(%s)", sub(" \\[.*\\]$", "", det_short), det_unit),
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.35, max_lwd = 2.8
  )
  legend("topleft", legend = "Observed RSN Station (same color scale)",
         pch = 21, pt.bg = "#f4d03f", col = "#111111", pt.cex = 1.1, bty = "n", cex = 0.82)

  plot_stream_network(
    basin_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. %s: S3N Kriging Standard Error (%s)", basin_title, det_short),
    sub_title  = "Uncertainty is lowest along monitored flow-connected paths and higher in headwaters",
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = if (use_log) "Kriging SE\n(log scale)" else sprintf("Kriging SE\n(%s)", det_unit),
    show_obs = FALSE, base_lwd = 0.35, max_lwd = 2.8
  )
  obs_b <- basin_pred[basin_pred$has_obs == 1, ]
  if (nrow(obs_b) > 0) {
    points(obs_b$mid_X, obs_b$mid_Y, pch = 3, col = "#c0392b", cex = 0.55, lwd = 0.8)
  }
  legend("topleft", legend = sprintf("RSN Monitoring Site Location (n = %d)", n_obs_basin),
         pch = 3, col = "#c0392b", pt.cex = 0.9, bty = "n", cex = 0.82)
}

# Figure 2: Humber River Basin High-Resolution S3N Map
plot_fig2_humber_map <- function(hum_pred, cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                                 reverse_color = FALSE, sc_se = NULL) {
  plot_basin_map(hum_pred, "Humber Basin", "Humber Basin (CatchID=80)",
                 cv_tab, sc_val, det_code, det_name, det_unit, use_log, out_file,
                 reverse_color = reverse_color, sc_se = sc_se)
}

# Figure 3: Thames & Severn River Basins S3N Maps
plot_fig3_thames_severn_map <- function(tha_pred, sev_pred, sc_val, det_code, det_name, det_unit, use_log, out_file, reverse_color = FALSE) {
  tha_pred <- normalize_pred_df(tha_pred)
  sev_pred <- normalize_pred_df(sev_pred)
  pal_use  <- get_val_pal(reverse_color)
  png(out_file, width = 2600, height = 1250, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 3.0, 1.2), bg = "white")

  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, det_unit)

  plot_stream_network(
    tha_pred, "pred_tu_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. Thames Basin: S3N Tail-Up Predicted %s (%s Reaches)",
                         det_short, format(nrow(tha_pred), big.mark = ",")),
    sub_title  = sprintf("n = %d observed RSN stations | Flow-connected Tail-Up stream network predictions",
                         sum(tha_pred$has_obs == 1)),
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.40, max_lwd = 2.8
  )

  plot_stream_network(
    sev_pred, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("B. Severn Basin: S3N Predicted %s (%s Reaches)",
                         det_short, format(nrow(sev_pred), big.mark = ",")),
    sub_title  = sprintf("n = %d observed RSN stations | Combined Tail-Up + Tail-Down S3N predictions",
                         sum(sev_pred$has_obs == 1)),
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = "obs_model_sc", base_lwd = 0.40, max_lwd = 2.8
  )
}

# Figure 4: 4-Panel Model Diagnostics, Cross-Validation & Parameter CIs
plot_fig4_diagnostics <- function(eng_pred, eng_cv, hum_cv, cv_tab, param_tab, det_code, det_name, det_unit, use_log, out_file,
                                  extra_cv_list = list(), lims_a = NULL, ylim_d = NULL) {
  eng_pred <- normalize_pred_df(eng_pred)
  eng_cv   <- normalize_cv_pts(eng_cv)
  hum_cv   <- normalize_cv_pts(hum_cv)
  cv_tab   <- normalize_cv_tab(cv_tab)

  png(out_file, width = 2400, height = 1800, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.5, 5.2, 3.2, 1.8), bg = "white")

  det_short <- shorten_det_name(det_name)
  sc_lbl    <- if (use_log) sprintf("log(%s, %s)", det_short, det_unit) else sprintf("%s (%s)", det_short, det_unit)

  # Panel A: Observed vs Predicted (5-Fold CV) on Humber, Custom Basin(s) & All England
  if (is.null(lims_a)) {
    lims <- range(c(eng_cv$obs_y, eng_cv$cv_comb), na.rm = TRUE)
    pad  <- 0.08 * diff(lims)
    lims <- c(lims[1] - pad, lims[2] + pad)
  } else {
    lims <- lims_a
  }

  plot(eng_cv$obs_y, eng_cv$cv_comb,
       pch = 16, col = rgb(0.16, 0.45, 0.65, 0.38), cex = 0.85,
       xlim = lims, ylim = lims,
       xlab = paste("Observed", sc_lbl),
       ylab = paste("5-Fold CV Predicted", sc_lbl),
       main = "A. Out-of-Sample 5-Fold Cross-Validation Accuracy")
  grid(col = "#e0e0e0", lty = "dotted")
  if (!is.null(hum_cv) && nrow(hum_cv) > 0) {
    points(hum_cv$obs_y, hum_cv$cv_comb, pch = 21, bg = rgb(0.90, 0.49, 0.13, 0.75), col = "#222222", cex = 0.95, lwd = 0.6)
  }

  eng_row <- subset(cv_tab, Domain == "All England RSN Network" & Model == "S3N Combined (Tail-Up + Tail-Down)")
  hum_row <- subset(cv_tab, Domain == "Humber Basin (CatchID=80)" & Model == "S3N Combined (Tail-Up + Tail-Down)")

  leg_txt  <- sprintf("All England (n = %d): R2 = %.3f, RMSE = %.3f", nrow(eng_cv), eng_row$R2[1], eng_row$RMSE[1])
  leg_pch  <- 16
  leg_lty  <- NA
  leg_lwd  <- NA
  leg_ptbg <- NA
  leg_col  <- rgb(0.16, 0.45, 0.65, 0.8)

  if (nrow(hum_row) > 0) {
    leg_txt  <- c(leg_txt, sprintf("Humber Basin (n = %d): R2 = %.3f, RMSE = %.3f", hum_row$n_obs[1], hum_row$R2[1], hum_row$RMSE[1]))
    leg_pch  <- c(leg_pch, 21)
    leg_lty  <- c(leg_lty, NA)
    leg_lwd  <- c(leg_lwd, NA)
    leg_ptbg <- c(leg_ptbg, rgb(0.90, 0.49, 0.13, 0.85))
    leg_col  <- c(leg_col, "#222222")
  }

  extra_pchs <- c(24, 22, 23, 25)
  extra_bgs  <- c(rgb(0.15, 0.68, 0.38, 0.85), rgb(0.56, 0.27, 0.68, 0.85),
                  rgb(0.95, 0.61, 0.07, 0.85), rgb(0.09, 0.63, 0.52, 0.85))
  if (length(extra_cv_list) > 0) {
    for (k in seq_along(extra_cv_list)) {
      ecv <- extra_cv_list[[k]]
      pts_k <- normalize_cv_pts(ecv$cv_pts)
      if (!is.null(pts_k) && nrow(pts_k) > 0) {
        pk <- extra_pchs[(k - 1) %% length(extra_pchs) + 1]
        bk <- extra_bgs[(k - 1) %% length(extra_bgs) + 1]
        points(pts_k$obs_y, pts_k$cv_comb, pch = pk, bg = bk, col = "#111111", cex = 1.0, lwd = 0.7)
        erow <- subset(cv_tab, Domain == ecv$domain_label & Model == "S3N Combined (Tail-Up + Tail-Down)")
        if (nrow(erow) > 0) {
          leg_txt  <- c(leg_txt, sprintf("%s (n = %d): R2 = %.3f, RMSE = %.3f", ecv$short_title, erow$n_obs[1], erow$R2[1], erow$RMSE[1]))
          leg_pch  <- c(leg_pch, pk)
          leg_lty  <- c(leg_lty, NA)
          leg_lwd  <- c(leg_lwd, NA)
          leg_ptbg <- c(leg_ptbg, bk)
          leg_col  <- c(leg_col, "#111111")
        }
      }
    }
  }

  abline(0, 1, col = "#c0392b", lwd = 2, lty = 2)
  leg_txt  <- c(leg_txt, "1:1 Perfect Prediction Line")
  leg_pch  <- c(leg_pch, NA)
  leg_lty  <- c(leg_lty, 2)
  leg_lwd  <- c(leg_lwd, 2)
  leg_ptbg <- c(leg_ptbg, NA)
  leg_col  <- c(leg_col, "#c0392b")

  legend("topleft", legend = leg_txt, pch = leg_pch, lty = leg_lty, lwd = leg_lwd,
         pt.bg = leg_ptbg, col = leg_col, bty = "n", cex = 0.80)

  # Panel B: 5-Fold CV R2 Comparison Across Domains
  avail_doms <- unique(cv_tab$Domain)
  r2_mat <- sapply(avail_doms, function(d) {
    sub <- cv_tab[cv_tab$Domain == d, ]
    pmax(0, c(
      sub$R2[sub$Model == "OLS (Non-Spatial)"],
      sub$R2[sub$Model == "S3N Tail-Down (Stream Distance)"],
      sub$R2[sub$Model == "S3N Tail-Up (Flow-Connected + AFV)"],
      sub$R2[sub$Model == "S3N Combined (Tail-Up + Tail-Down)"]
    ))
  })
  short_map <- c(
    "All England RSN Network"     = sprintf("All England\n(n=%d)", eng_row$n_obs[1]),
    "Humber Basin (CatchID=80)"   = "Humber Basin",
    "Severn Basin (CatchID=2786)" = "Severn Basin",
    "Thames Basin (CatchID=3055)" = "Thames Basin"
  )
  colnames(r2_mat) <- ifelse(avail_doms %in% names(short_map), short_map[avail_doms],
                             sub(" \\((CatchID|RBD)=.*\\)$", "", avail_doms))
  ymax_b <- max(0.88, max(r2_mat, na.rm = TRUE) * 1.18)
  cex_names_b <- if (ncol(r2_mat) > 4) 0.78 else 0.90
  cex_txt_b   <- if (ncol(r2_mat) > 4) 0.64 else 0.72
  bp <- barplot(r2_mat, beside = TRUE, ylim = c(0, ymax_b),
                col = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"),
                border = "#222222", cex.names = cex_names_b, ylab = "5-Fold Cross-Validation R2",
                main = "B. Out-of-Sample Variance Explained (R2) by Model")
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  text(bp, r2_mat + 0.03 * ymax_b, labels = sprintf("%.2f", r2_mat), cex = cex_txt_b, font = 2)
  legend("topright",
         legend = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream)", "S3N Tail-Up (FC + AFV)", "S3N Combined (TU + TD)"),
         fill = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"), border = "#222222", bty = "n", cex = 0.80)

  # Panel C: Forest Plot of Fixed-Effect Estimates & 95% Bootstrap CIs (All England S3N)
  eng_td_p <- subset(param_tab, Domain == "All England RSN Network" & Model == "S3N Tail-Down (Stream)" &
                     !(Parameter %in% c("Intercept", "sigma.sq", "tau.sq", "phi")))
  eng_tu_p <- subset(param_tab, Domain == "All England RSN Network" & Model == "S3N Tail-Up (FC + AFV)" &
                     !(Parameter %in% c("Intercept", "sigma.sq", "tau.sq", "phi")))
  if (nrow(eng_td_p) == 0) {
    eng_td_p <- subset(param_tab, Domain == "All England RSN Network" & Model == "S3N Tail-Down (Stream)" &
                       !(Parameter %in% c("sigma.sq", "tau.sq", "phi")))
    eng_tu_p <- subset(param_tab, Domain == "All England RSN Network" & Model == "S3N Tail-Up (FC + AFV)" &
                       !(Parameter %in% c("sigma.sq", "tau.sq", "phi")))
  }

  nice_labels <- c(
    Intercept          = "Intercept (Spatial Mean)",
    Elevation_z        = "Elevation (std)",
    Slope_z            = "Channel Slope (std)",
    log10_US_Accum_z   = "Upstream Length (std)",
    Sinuosity_z        = "Channel Sinuosity (std)",
    SrcElevation_z     = "Source Elevation (std)",
    Dist2Mth_z         = "Distance to Mouth (std)",
    Length_z           = "Reach Length (std)",
    log10_Shreve_z     = "Shreve Order (log10 std)",
    Strahler_z         = "Strahler Order (std)",
    Is_Chalk_Limestone = "Chalk / Limestone Aquifer",
    Is_Sandstone       = "Sandstone Aquifer",
    Is_Aquifer         = "Aquifer Bedrock",
    Is_Peat            = "Peat Superficial Geology",
    Easting_z          = "Easting (std)",
    Northing_z         = "Northing (std)"
  )
  np <- nrow(eng_td_p)
  y_pos <- rev(seq_len(np))
  par_MAR_old <- par("mar")
  par(mar = c(4.5, 10.2, 3.2, 1.5))
  xr_ci <- range(c(eng_tu_p$CI_Lower, eng_tu_p$CI_Upper, eng_td_p$CI_Lower, eng_td_p$CI_Upper, -0.2, 0.2), na.rm = TRUE)
  xr_pad <- 0.12 * diff(xr_ci)
  plot(NA, xlim = c(xr_ci[1] - xr_pad, xr_ci[2] + xr_pad), ylim = c(0.4, np + 0.6),
       yaxt = "n", xlab = sprintf("Effect Size on %s [95%% Bootstrap CI]", sc_lbl), ylab = "",
       main = sprintf("C. Environmental Drivers of %s (England S3N)", det_short))
  grid(nx = NULL, ny = NA, col = "#e0e0e0", lty = "dotted")
  abline(v = 0, col = "#c0392b", lty = 2, lwd = 1.5)
  lbls_c <- ifelse(eng_td_p$Parameter %in% names(nice_labels), nice_labels[eng_td_p$Parameter], eng_td_p$Parameter)
  axis(2, at = y_pos, labels = lbls_c, las = 1, cex.axis = 0.85)

  segments(eng_tu_p$CI_Lower, y_pos + 0.14, eng_tu_p$CI_Upper, y_pos + 0.14, col = "#27ae60", lwd = 2.2)
  points(eng_tu_p$Estimate, y_pos + 0.14, pch = 16, col = "#27ae60", cex = 1.1)

  segments(eng_td_p$CI_Lower, y_pos - 0.14, eng_td_p$CI_Upper, y_pos - 0.14, col = "#2980b9", lwd = 2.2)
  points(eng_td_p$Estimate, y_pos - 0.14, pch = 17, col = "#2980b9", cex = 1.1)

  legend("bottomright", legend = c("S3N Tail-Up (95% Bootstrap CI)", "S3N Tail-Down (95% Bootstrap CI)"),
         col = c("#27ae60", "#2980b9"), pch = c(16, 17), lwd = 2.2, bty = "n", cex = 0.80)
  par(mar = par_MAR_old)

  # Panel D: Distribution across All Reaches (Overall) and by Hydrogeology & Elevation Regime
  reg_str <- ifelse(eng_pred$Elevation > 150, "Upland\n(>150m)",
             ifelse(eng_pred$Is_Aquifer == 1, "Lowland Aquifer\n(Chalk/Limestone/Sand)", "Lowland\nNon-Aquifer"))
  box_list <- list(
    "All Reaches\n(Overall)"                  = eng_pred$pred_combined,
    "Upland\n(>150m)"                         = eng_pred$pred_combined[reg_str == "Upland\n(>150m)"],
    "Lowland\nNon-Aquifer"                    = eng_pred$pred_combined[reg_str == "Lowland\nNon-Aquifer"],
    "Lowland Aquifer\n(Chalk/Limestone/Sand)" = eng_pred$pred_combined[reg_str == "Lowland Aquifer\n(Chalk/Limestone/Sand)"]
  )

  bp_d <- boxplot(box_list, outline = FALSE, ylim = ylim_d,
                  col = c("#95a5a6", "#2874a6", "#1abc9c", "#f4d03f"), border = "#222222",
                  cex.axis = 0.84,
                  ylab = sprintf("S3N Predicted %s (%s)", det_short, det_unit),
                  xlab = "Overall Network vs. Hydrogeological & Topographic Regime",
                  main = sprintf("D. Predicted %s Across %s Reaches by Regime",
                                 det_short, format(nrow(eng_pred), big.mark = ",")))
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(v = 1.5, col = "#7f8c8d", lty = 2, lwd = 1.3)
  medians <- sapply(box_list, median, na.rm = TRUE)
  yr_d <- par("usr")[4] - par("usr")[3]
  y_txt <- pmin(par("usr")[4] - 0.07 * yr_d, bp_d$stats[4, ] + 0.07 * yr_d)
  text(1:4, y_txt, labels = sprintf("Median:\n%s %s", formatC(medians, format = "fg", digits = 3), det_unit),
       font = 2, cex = 0.78, col = "#111111")
}

# Figure 5 (Temporal Multi-Panel Grid): All England S3N Comparison Across Seasons, Years, or Months
plot_fig5_temporal_england_grid <- function(period_items, by_mode, sc_val, det_code, det_name, det_unit,
                                            use_log, out_file, reverse_color = FALSE) {
  period_items <- Filter(Negate(is.null), period_items)
  n_pan <- length(period_items)
  if (n_pan == 0) return(invisible(NULL))

  # Order seasonal panels in canonical Spring -> Summer -> Autumn -> Winter sequence
  if (tolower(by_mode) == "season") {
    s_order <- c("spring", "summer", "autumn", "winter")
    keys <- sapply(period_items, function(x) {
      lbl <- tolower(x$label)
      hit <- which(sapply(s_order, function(s) grepl(s, lbl)))
      if (length(hit) > 0) hit[1] else 99L
    })
    period_items <- period_items[order(keys)]
  }

  if (n_pan <= 2) {
    nrow_g <- 1L; ncol_g <- 2L; w_px <- 2600L; h_px <- 1450L
  } else if (n_pan <= 4) {
    nrow_g <- 2L; ncol_g <- 2L; w_px <- 2600L; h_px <- 2750L
  } else if (n_pan <= 6) {
    nrow_g <- 2L; ncol_g <- 3L; w_px <- 3700L; h_px <- 2700L
  } else if (n_pan <= 9) {
    nrow_g <- 3L; ncol_g <- 3L; w_px <- 3700L; h_px <- 3800L
  } else {
    nrow_g <- 3L; ncol_g <- 4L; w_px <- 4600L; h_px <- 3700L
  }

  pal_use   <- get_val_pal(reverse_color)
  det_short <- sub(" \\[.*\\]$", "", shorten_det_name(det_name))
  cbar_lbl  <- sprintf("%s\n(%s)", det_short, det_unit)
  mode_cap  <- paste0(toupper(substr(by_mode, 1, 1)), tolower(substr(by_mode, 2, nchar(by_mode))))

  png(out_file, width = w_px, height = h_px, res = 200, type = "cairo")
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nrow_g, ncol_g), mar = c(1.0, 1.0, 2.8, 1.0), oma = c(0.5, 0.5, 3.4, 0.5), bg = "white")

  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  for (i in seq_along(period_items)) {
    item   <- period_items[[i]]
    eng_df <- normalize_pred_df(item$eng_pred)
    n_obs_i <- sum(eng_df$has_obs == 1, na.rm = TRUE)
    med_raw <- median(eng_df$pred_combined, na.rm = TRUE)

    # Format clean panel title
    raw_lbl <- item$label
    if (grepl("^Season:?\\s*", raw_lbl, ignore.case = TRUE)) {
      s_nm <- tolower(trimws(sub("^Season:?\\s*", "", raw_lbl, ignore.case = TRUE)))
      s_pretty <- switch(s_nm,
        spring = "Spring (Mar-May)",
        summer = "Summer (Jun-Aug)",
        autumn = "Autumn (Sep-Nov)",
        winter = "Winter (Dec-Feb)",
        paste0(toupper(substr(s_nm, 1, 1)), substr(s_nm, 2, nchar(s_nm)))
      )
      pan_lbl <- s_pretty
    } else if (grepl("^Month\\s+([0-9]+)$", raw_lbl, ignore.case = TRUE)) {
      m_idx <- as.integer(sub("^Month\\s+([0-9]+)$", "\\1", raw_lbl, ignore.case = TRUE))
      pan_lbl <- if (!is.na(m_idx) && m_idx >= 1 && m_idx <= 12) {
        sprintf("%s (Month %02d)", month_names[m_idx], m_idx)
      } else {
        raw_lbl
      }
    } else {
      pan_lbl <- raw_lbl
    }

    r2_str <- ""
    if (!is.null(item$cv_tab)) {
      cvt <- normalize_cv_tab(item$cv_tab)
      erow <- subset(cvt, Domain == "All England RSN Network" & Model == "S3N Combined (Tail-Up + Tail-Down)")
      if (nrow(erow) > 0 && !is.na(erow$R2[1])) {
        r2_str <- sprintf(" | CV R2 = %.3f", erow$R2[1])
      }
    }

    main_t <- sprintf("%s. %s — Predicted %s", LETTERS[i], pan_lbl, det_short)
    sub_t  <- sprintf("n = %s obs reaches%s | Network Median = %s %s",
                      format(n_obs_i, big.mark = ","), r2_str,
                      formatC(med_raw, format = "fg", digits = 3), det_unit)

    plot_stream_network(
      eng_df, "pred_model_sc", vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
      main_title = main_t,
      sub_title  = sub_t,
      cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
      show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
    )
  }

  n_reaches_tot <- nrow(period_items[[1]]$eng_pred)
  n_catch_tot   <- length(unique(period_items[[1]]$eng_pred$CatchID))
  mtext(sprintf("All England S3N Stream Network Comparison by %s: %s (%s)",
                mode_cap, det_short, det_unit),
        side = 3, line = 1.5, outer = TRUE, cex = 1.25, font = 2, col = "#111111")
  mtext(sprintf("Unified color scale across all %d %s panels (%s reaches across %d monitored catchments)",
                n_pan, tolower(by_mode), format(n_reaches_tot, big.mark = ","), n_catch_tot),
        side = 3, line = 0.2, outer = TRUE, cex = 0.92, col = "#444444")
}

# ==============================================================================
# SHARED FIXED-EFFECT DESIGN MATRIX BUILDER & COVARIATE SELECTION (--vars=...)
# ==============================================================================

COVARIATE_CATALOG <- list(
  Elevation_z = list(
    aliases = c("elevation", "elev", "elevation_z", "os_t50_ele"),
    col     = "Elevation",
    desc    = "Mean channel elevation (m AOD, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$Elevation))
  ),
  Slope_z = list(
    aliases = c("slope", "slope_z"),
    col     = "Slope",
    desc    = "Mean channel gradient (%, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$Slope))
  ),
  log10_US_Accum_z = list(
    aliases = c("us_accum", "accum", "us_accum_km", "log10_us_accum", "log10_us_accum_z"),
    col     = "US_Accum_km",
    desc    = "log10 upstream accumulated channel length (km, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(log10(pmax(1e-4, df$US_Accum_km))))
  ),
  Is_Chalk_Limestone = list(
    aliases = c("chalk_limestone", "chalk", "limestone", "is_chalk_limestone"),
    col     = "Is_Chalk_Limestone",
    desc    = "Binary indicator (0/1) for BGS Chalk / Limestone / Dolostone bedrock",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Chalk_Limestone)
  ),
  Is_Sandstone = list(
    aliases = c("sandstone", "sand", "is_sandstone"),
    col     = "Is_Sandstone",
    desc    = "Binary indicator (0/1) for BGS Sandstone aquifer bedrock",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Sandstone)
  ),
  Is_Peat = list(
    aliases = c("peat", "is_peat"),
    col     = "Is_Peat",
    desc    = "Binary indicator (0/1) for BGS Peat superficial deposits",
    default = "All England (default)",
    eval_fn = function(df) as.numeric(df$Is_Peat)
  ),
  Is_Aquifer = list(
    aliases = c("aquifer", "is_aquifer"),
    col     = "Is_Aquifer",
    desc    = "Composite binary indicator (0/1) for Chalk/Limestone or Sandstone aquifer",
    default = "Single-basin fits (default)",
    eval_fn = function(df) as.numeric(df$Is_Aquifer)
  ),
  Easting_z = list(
    aliases = c("easting", "mid_x", "easting_z"),
    col     = "mid_X",
    desc    = "Reach midpoint BNG Easting coordinate (m, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$mid_X))
  ),
  Northing_z = list(
    aliases = c("northing", "mid_y", "northing_z"),
    col     = "mid_Y",
    desc    = "Reach midpoint BNG Northing coordinate (m, z-score standardized)",
    default = "All domains",
    eval_fn = function(df) as.numeric(scale(df$mid_Y))
  ),
  Sinuosity_z = list(
    aliases = c("sinuosity", "sinuosity_z"),
    col     = "Sinuosity",
    desc    = "Mean channel sinuosity ratio (z-score standardized)",
    default = "Optional (--vars=+sinuosity)",
    eval_fn = function(df) as.numeric(scale(df$Sinuosity))
  ),
  SrcElevation_z = list(
    aliases = c("src_elevation", "srcelevation", "src_elev", "src_elevation_z"),
    col     = "SrcElevation",
    desc    = "Mean headwater source elevation (m AOD, z-score standardized)",
    default = "Optional (--vars=+src_elevation)",
    eval_fn = function(df) as.numeric(scale(df$SrcElevation))
  ),
  Dist2Mth_z = list(
    aliases = c("dist2mth", "dist2mth_km", "dist2mth_z"),
    col     = "Dist2Mth_km",
    desc    = "Hydrologic distance from upstream node to river mouth (km, z-score)",
    default = "Optional (--vars=+dist2mth)",
    eval_fn = function(df) as.numeric(scale(df$Dist2Mth_km))
  ),
  Length_z = list(
    aliases = c("length", "length_km", "length_z"),
    col     = "length_km",
    desc    = "Reach polyline length (km, z-score standardized)",
    default = "Optional (--vars=+length)",
    eval_fn = function(df) as.numeric(scale(df$length_km))
  ),
  log10_Shreve_z = list(
    aliases = c("shreve", "log10_shreve", "log10_shreve_z"),
    col     = "Shreve",
    desc    = "log10 Shreve stream magnitude (z-score standardized)",
    default = "Optional (--vars=+shreve)",
    eval_fn = function(df) as.numeric(scale(log10(pmax(1, df$Shreve))))
  ),
  Strahler_z = list(
    aliases = c("strahler", "strahler_z"),
    col     = "Strahler",
    desc    = "Strahler stream order (1..7, z-score standardized)",
    default = "Optional (--vars=+strahler)",
    eval_fn = function(df) as.numeric(scale(df$Strahler))
  )
)

print_available_covariates <- function() {
  cat("Available Fixed-Effect Regression Covariates (--vars=<list>):\n\n")
  cat(sprintf("%-20s | %-28s | %-18s | %s\n",
              "Model Term", "Short CLI Aliases", "Default Status", "Description"))
  cat(strrep("-", 118), "\n")
  for (nm in names(COVARIATE_CATALOG)) {
    info <- COVARIATE_CATALOG[[nm]]
    al_str <- paste(head(info$aliases, 2), collapse = ", ")
    cat(sprintf("%-20s | %-28s | %-18s | %s\n",
                nm, al_str, substr(info$default, 1, 18), info$desc))
  }
  cat("\nSpecial --vars Modes:\n")
  cat("  --vars=default                     # Default covariates (Elevation, Slope, US_Accum, Geology, Easting, Northing)\n")
  cat("  --vars=elevation,slope,aquifer     # Explicit subset of covariates (Intercept is always included)\n")
  cat("  --vars=-easting,-northing          # Start from default covariates and remove Easting & Northing\n")
  cat("  --vars=+sinuosity,+strahler        # Start from default covariates and add Sinuosity & Strahler order\n")
  cat("  --vars=none                        # Intercept-only model (pure spatial stream network kriging)\n")
  cat("  --vars=all                         # All available ORN topographic, channel, geological, and coordinate covariates\n")
}

match_covariate_token <- function(tok) {
  tok_clean <- tolower(trimws(tok))
  for (nm in names(COVARIATE_CATALOG)) {
    if (tok_clean == tolower(nm) || tok_clean %in% COVARIATE_CATALOG[[nm]]$aliases) {
      return(nm)
    }
  }
  valid_aliases <- sapply(names(COVARIATE_CATALOG), function(nm) COVARIATE_CATALOG[[nm]]$aliases[1])
  stop(sprintf("Unknown fixed-effect covariate '%s'. Run with --list-vars to see valid options (%s).",
               tok, paste(valid_aliases, collapse = ", ")))
}

resolve_covariate_terms <- function(df, include_peat = TRUE, var_specs = NULL) {
  has_split_geol <- include_peat &&
    ("Is_Peat" %in% colnames(df)) && ("Is_Sandstone" %in% colnames(df)) &&
    sd(df$Is_Peat[df$has_obs == 1], na.rm = TRUE) > 0 &&
    sd(df$Is_Sandstone[df$has_obs == 1], na.rm = TRUE) > 0

  default_terms <- if (has_split_geol) {
    c("Elevation_z", "Slope_z", "log10_US_Accum_z",
      "Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat",
      "Easting_z", "Northing_z")
  } else {
    c("Elevation_z", "Slope_z", "log10_US_Accum_z",
      "Is_Aquifer", "Easting_z", "Northing_z")
  }

  if (is.null(var_specs) || length(var_specs) == 0) {
    return(default_terms)
  }
  specs <- trimws(var_specs[nchar(trimws(var_specs)) > 0])
  if (length(specs) == 0 || (length(specs) == 1 && tolower(specs[1]) == "default")) {
    return(default_terms)
  }
  if (length(specs) == 1 && tolower(specs[1]) %in% c("none", "intercept", "null", "0")) {
    return(character(0))
  }
  if (length(specs) == 1 && tolower(specs[1]) == "all") {
    all_t <- c("Elevation_z", "Slope_z", "log10_US_Accum_z",
               "Sinuosity_z", "SrcElevation_z", "Dist2Mth_z", "Length_z",
               "log10_Shreve_z", "Strahler_z",
               if (has_split_geol) c("Is_Chalk_Limestone", "Is_Sandstone", "Is_Peat") else "Is_Aquifer",
               "Easting_z", "Northing_z")
    return(all_t[sapply(all_t, function(nm) COVARIATE_CATALOG[[nm]]$col %in% colnames(df))])
  }

  # Check if relative modifier mode (+var / -var) is used
  is_mod <- grepl("^[+-]", specs)
  if (all(is_mod)) {
    terms <- default_terms
    for (sp in specs) {
      op  <- substr(sp, 1, 1)
      nm  <- match_covariate_token(substr(sp, 2, nchar(sp)))
      if (op == "+") {
        if (!(nm %in% terms)) terms <- c(terms, nm)
      } else {
        # If removing 'aquifer', also remove split geology terms (and vice versa) if user wrote -aquifer or -geology
        if (nm == "Is_Aquifer") {
          terms <- setdiff(terms, c("Is_Aquifer", "Is_Chalk_Limestone", "Is_Sandstone"))
        } else {
          terms <- setdiff(terms, nm)
        }
      }
    }
    return(terms)
  }

  # Explicit list of covariates (also allowing + prefix)
  terms <- character(0)
  for (sp in specs) {
    sp_clean <- sub("^\\+", "", sp)
    if (tolower(sp_clean) %in% c("intercept", "none")) next
    nm <- match_covariate_token(sp_clean)
    if (!(nm %in% terms)) terms <- c(terms, nm)
  }
  terms
}

build_design_matrix <- function(df, include_peat = TRUE, var_specs = NULL) {
  terms <- resolve_covariate_terms(df, include_peat = include_peat, var_specs = var_specs)
  X <- matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, "Intercept"))
  if (length(terms) == 0L) {
    return(X)
  }
  col_list <- list(Intercept = rep(1.0, nrow(df)))
  for (nm in terms) {
    info <- COVARIATE_CATALOG[[nm]]
    if (!(info$col %in% colnames(df))) {
      warning(sprintf("Column '%s' for covariate '%s' not found in dataframe; skipping.", info$col, nm))
      next
    }
    vec <- info$eval_fn(df)
    if (all(!is.na(vec)) && sd(vec[df$has_obs == 1], na.rm = TRUE) > 0) {
      col_list[[nm]] <- vec
    }
  }
  do.call(cbind, col_list)
}

