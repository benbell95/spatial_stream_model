# ==============================================================================
# s3n_plot_functions.R
#
# Modular base-R plotting library for Scalable Spatial Stream Network (S3N)
# models. Automatically sources:
#   - scripts/s3n_config.R      (paths, directories, source/network/det configs)
#   - scripts/s3n_transforms.R  (transforms, delta-method variance, normalisation)
#   - scripts/s3n_covariates.R  (covariate catalog, design matrix, core basins)
# ==============================================================================

.s3n_find_scripts_dir <- function() {
  if (file.exists("scripts/s3n_config.R")) return("scripts")
  if (file.exists("s3n_config.R")) return(".")
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- sys.frame(i)
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      d <- dirname(get("ofile", envir = fr, inherits = FALSE))
      if (file.exists(file.path(d, "s3n_config.R"))) return(d)
    }
  }
  f_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f_arg) > 0) {
    d <- dirname(sub("^--file=", "", f_arg[1]))
    if (file.exists(file.path(d, "s3n_config.R"))) return(d)
  }
  "scripts"
}

.s3n_mod_dir <- .s3n_find_scripts_dir()
for (.mod in c("s3n_config.R", "s3n_transforms.R", "s3n_covariates.R")) {
  .mod_p <- file.path(.s3n_mod_dir, .mod)
  if (file.exists(.mod_p)) source(.mod_p, local = FALSE)
}
rm(.mod, .mod_p, .s3n_mod_dir)


open_png_device <- function(out_file, width, height, res = 200) {
  if (isTRUE(capabilities("cairo"))) {
    png(out_file, width = width, height = height, res = res, type = "cairo")
  } else {
    png(out_file, width = width, height = height, res = res)
  }
}

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
  rect(x0, ys[seq_len(n)], x1, ys[seq_len(n) + 1L], col = pal, border = NA)
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

compute_color_scale <- function(vals_model_sc, trans_obj, vmin_raw = NULL, vmax_raw = NULL) {
  tr <- as_transformer(trans_obj)
  qs <- quantile(vals_model_sc, probs = c(0.01, 0.99), na.rm = TRUE)
  vmin <- as.numeric(qs[1])
  vmax <- as.numeric(qs[2])
  has_vmin <- !is.null(vmin_raw) && length(vmin_raw) > 0 && is.finite(suppressWarnings(as.numeric(vmin_raw[1])))
  has_vmax <- !is.null(vmax_raw) && length(vmax_raw) > 0 && is.finite(suppressWarnings(as.numeric(vmax_raw[1])))
  if (has_vmin) {
    vmin <- tr$fwd(as.numeric(vmin_raw[1]))
  }
  if (has_vmax) {
    vmax <- tr$fwd(as.numeric(vmax_raw[1]))
  }
  if (!is.finite(vmin) || !is.finite(vmax) || vmax - vmin < 1e-4) {
    vmin <- ifelse(is.finite(vmin), vmin, 0.0) - 0.5
    vmax <- ifelse(is.finite(vmax), vmax, 0.0) + 0.5
  }
  if (tr$is_transformed) {
    raw_min <- tr$inv(vmin)
    raw_max <- tr$inv(vmax)
    if (raw_max < raw_min) {
      tmp <- raw_min; raw_min <- raw_max; raw_max <- tmp
    }
    if (tr$type %in% c("log", "log10", "reciprocal") ||
        (tr$type %in% c("boxcox", "yeojohnson") && !is.na(tr$lambda) && tr$lambda <= 0.15)) {
      cands <- c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5,
                 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
    } else {
      p_ticks <- pretty(c(raw_min, raw_max), n = 6)
      log_cands <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 100, 200, 500, 1000)
      cands <- sort(unique(c(p_ticks, log_cands)))
      cands <- cands[tr$valid_fn(cands)]
    }
    lo_b <- if (raw_min > 0) raw_min * 0.95 else raw_min - 0.02 * max(1e-6, raw_max - raw_min)
    hi_b <- if (raw_max > 0) raw_max * 1.05 else raw_max + 0.02 * max(1e-6, raw_max - raw_min)
    raw_in <- cands[cands >= lo_b & cands <= hi_b]
    if (length(raw_in) > 7) {
      raw_in <- raw_in[unique(round(seq(1, length(raw_in), length.out = 6)))]
    }
    if (length(raw_in) >= 3) {
      ticks <- tr$fwd(raw_in)
      labels <- formatC(raw_in, format = "fg", digits = 2)
      if (!has_vmin) vmin <- min(vmin, min(ticks))
      if (!has_vmax) vmax <- max(vmax, max(ticks))
    } else {
      ticks <- seq(vmin, vmax, length.out = 6)
      labels <- formatC(tr$inv(ticks), format = "fg", digits = 2)
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
  if (!is.null(vmin_se) && length(vmin_se) > 0 && is.finite(suppressWarnings(as.numeric(vmin_se[1])))) vmin <- max(0.001, as.numeric(vmin_se[1]))
  if (!is.null(vmax_se) && length(vmax_se) > 0 && is.finite(suppressWarnings(as.numeric(vmax_se[1])))) vmax <- max(vmin + 0.05, as.numeric(vmax_se[1]))
  ticks <- pretty(c(vmin, vmax), n = 5)
  ticks <- ticks[ticks >= vmin & ticks <= vmax]
  if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
  list(vmin = vmin, vmax = vmax, ticks = ticks, labels = sprintf("%.2f", ticks))
}

# ------------------------------------------------------------------------------
# Prediction Normalisation Utilities (Z-score, Robust Z-score, Min-Max, Percentile)
# ------------------------------------------------------------------------------
parse_norm_mode <- function(args = character(0), default = "none") {
  mode <- default
  for (a in args) {
    if (grepl("^--(norm|normalise|normalize)=", a, ignore.case = TRUE)) {
      raw_m <- tolower(trimws(sub("^--(norm|normalise|normalize)=", "", a, ignore.case = TRUE)))
      mode <- switch(raw_m,
        "none"          = "none",
        "raw"           = "none",
        "z"             = "zscore",
        "zscore"        = "zscore",
        "z-score"       = "zscore",
        "z_score"       = "zscore",
        "standard"      = "zscore",
        "robust"        = "robust_zscore",
        "robust_z"      = "robust_zscore",
        "robust-z"      = "robust_zscore",
        "robust_zscore" = "robust_zscore",
        "robust-zscore" = "robust_zscore",
        "mad"           = "robust_zscore",
        "minmax"        = "minmax",
        "min-max"       = "minmax",
        "min_max"       = "minmax",
        "01"            = "minmax",
        "pct"           = "percentile",
        "rank"          = "percentile",
        "percentile"    = "percentile",
        stop(sprintf(
          "Unknown --norm='%s'. Valid options: none, zscore, robust_zscore, minmax, percentile",
          raw_m
        ))
      )
    }
  }
  mode
}

get_norm_meta <- function(norm_mode = "none") {
  m <- tolower(trimws(if (is.null(norm_mode) || !nzchar(norm_mode)) "none" else norm_mode))
  switch(m,
    zscore = list(
      mode         = "zscore",
      pred_col     = "pred_zscore",
      obs_col      = "obs_zscore",
      tu_col       = "pred_tu_zscore",
      td_col       = "pred_td_zscore",
      unit_label   = "Z-score (SD)",
      short_label  = "Z-score",
      file_tag     = "_zscore",
      is_diverging = TRUE
    ),
    robust_zscore = list(
      mode         = "robust_zscore",
      pred_col     = "pred_robust_zscore",
      obs_col      = "obs_robust_zscore",
      tu_col       = "pred_tu_robust_zscore",
      td_col       = "pred_td_robust_zscore",
      unit_label   = "Robust Z (MAD)",
      short_label  = "Robust Z",
      file_tag     = "_robust_zscore",
      is_diverging = TRUE
    ),
    minmax = list(
      mode         = "minmax",
      pred_col     = "pred_minmax",
      obs_col      = "obs_minmax",
      tu_col       = "pred_tu_minmax",
      td_col       = "pred_td_minmax",
      unit_label   = "Min-Max [0, 1]",
      short_label  = "Min-Max",
      file_tag     = "_minmax",
      is_diverging = FALSE
    ),
    percentile = list(
      mode         = "percentile",
      pred_col     = "pred_percentile",
      obs_col      = "obs_percentile",
      tu_col       = "pred_tu_percentile",
      td_col       = "pred_td_percentile",
      unit_label   = "Percentile (%)",
      short_label  = "Percentile",
      file_tag     = "_percentile",
      is_diverging = FALSE
    ),
    list(
      mode         = "none",
      pred_col     = "pred_model_sc",
      obs_col      = "obs_model_sc",
      tu_col       = "pred_tu_sc",
      td_col       = "pred_td_sc",
      unit_label   = NULL,
      short_label  = NULL,
      file_tag     = "",
      is_diverging = FALSE
    )
  )
}

get_norm_pal <- function(norm_mode = "none", reverse_color = FALSE) {
  nm <- get_norm_meta(norm_mode)
  if (isTRUE(nm$is_diverging)) {
    cols <- c("#1f4e78", "#3b82b6", "#8cc0de", "#f2f2ef", "#f4a582", "#d6604d", "#99000d")
    if (isTRUE(reverse_color)) cols <- rev(cols)
    return(colorRampPalette(cols)(256))
  }
  get_val_pal(reverse_color)
}

compute_norm_color_scale <- function(vals, norm_mode = "none", vmin_raw = NULL, vmax_raw = NULL) {
  nm <- get_norm_meta(norm_mode)
  has_vmin <- !is.null(vmin_raw) && length(vmin_raw) > 0 && is.finite(suppressWarnings(as.numeric(vmin_raw[1])))
  has_vmax <- !is.null(vmax_raw) && length(vmax_raw) > 0 && is.finite(suppressWarnings(as.numeric(vmax_raw[1])))
  if (nm$mode %in% c("zscore", "robust_zscore")) {
    qs <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)
    lim <- max(abs(qs), na.rm = TRUE)
    if (!is.finite(lim) || lim < 0.5) lim <- 2.0
    lim <- ceiling(lim * 2) / 2
    vmin <- if (has_vmin) as.numeric(vmin_raw[1]) else -lim
    vmax <- if (has_vmax) as.numeric(vmax_raw[1]) else lim
    if (vmax <= vmin) { vmin <- -2.0; vmax <- 2.0 }
    ticks <- pretty(c(vmin, vmax), n = 6)
    ticks <- ticks[ticks >= vmin - 1e-9 & ticks <= vmax + 1e-9]
    if (!(0 %in% ticks) && vmin < 0 && vmax > 0) ticks <- sort(unique(c(ticks, 0)))
    if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
    labels <- sprintf("%+.1f", ticks)
    labels[abs(ticks) < 1e-9] <- "0.0"
    return(list(vmin = vmin, vmax = vmax, ticks = ticks, labels = labels))
  }
  if (nm$mode == "minmax") {
    vmin <- if (has_vmin) as.numeric(vmin_raw[1]) else 0.0
    vmax <- if (has_vmax) as.numeric(vmax_raw[1]) else 1.0
    if (vmax <= vmin) { vmin <- 0.0; vmax <- 1.0 }
    ticks <- pretty(c(vmin, vmax), n = 5)
    ticks <- ticks[ticks >= vmin - 1e-9 & ticks <= vmax + 1e-9]
    if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
    return(list(vmin = vmin, vmax = vmax, ticks = ticks, labels = sprintf("%.2f", ticks)))
  }
  if (nm$mode == "percentile") {
    vmin <- if (has_vmin) as.numeric(vmin_raw[1]) else 0.0
    vmax <- if (has_vmax) as.numeric(vmax_raw[1]) else 100.0
    if (vmax <= vmin) { vmin <- 0.0; vmax <- 100.0 }
    ticks <- pretty(c(vmin, vmax), n = 5)
    ticks <- ticks[ticks >= vmin - 1e-9 & ticks <= vmax + 1e-9]
    if (length(ticks) < 3) ticks <- seq(vmin, vmax, length.out = 5)
    return(list(vmin = vmin, vmax = vmax, ticks = ticks, labels = sprintf("%.0f%%", ticks)))
  }
  compute_color_scale(vals, build_transformer("none"), vmin_raw = vmin_raw, vmax_raw = vmax_raw)
}


draw_plot_header <- function(main_title = NULL, sub_title = NULL,
                             max_main = 65, max_sub = 90,
                             cex_main = 0.96, cex_sub = 0.80,
                             font_main = 2, font_sub = 1,
                             col_main = "#111111", col_sub = "#444444",
                             adj = 0.5,
                             main_base = NULL, sub_base = 0.18) {
  main_w <- if (!is.null(main_title) && nzchar(main_title)) {
    strsplit(wrap_title(main_title, max_len = max_main), "\n", fixed = TRUE)[[1]]
  } else character(0)
  sub_w  <- if (!is.null(sub_title) && nzchar(sub_title)) {
    strsplit(wrap_title(sub_title, max_len = max_sub), "\n", fixed = TRUE)[[1]]
  } else character(0)

  n_sub <- length(sub_w)
  step_sub  <- cex_sub * 1.10
  step_main <- cex_main * 1.18

  if (n_sub > 0) {
    for (s_i in seq_len(n_sub)) {
      s_line <- sub_base + (n_sub - s_i) * step_sub
      mtext(sub_w[s_i], side = 3, line = s_line, adj = adj, cex = cex_sub, font = font_sub, col = col_sub)
    }
  }

  if (is.null(main_base)) {
    gap <- 0.65 + 0.30 * cex_main
    main_base <- if (n_sub > 0) {
      sub_base + (n_sub - 1) * step_sub + gap
    } else {
      if (adj == 0) 0.35 else 0.70
    }
  }

  n_main <- length(main_w)
  if (n_main > 0) {
    for (m_i in seq_len(n_main)) {
      m_line <- main_base + (n_main - m_i) * step_main
      mtext(main_w[m_i], side = 3, line = m_line, adj = adj, cex = cex_main, font = font_main, col = col_main)
    }
  }
}

plot_stream_network <- function(df, val_col, vmin, vmax, pal, main_title, sub_title,
                                cbar_ticks, cbar_labels, cbar_title,
                                show_obs = FALSE, obs_col = NULL, base_lwd = 0.35, max_lwd = 2.6,
                                max_main_chars = 65, max_sub_chars = 90) {
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

  # Wrap and render multi-line titles without overlapping or clipping
  draw_plot_header(
    main_title = main_title, sub_title = sub_title,
    max_main = max_main_chars, max_sub = max_sub_chars,
    cex_main = 0.96, cex_sub = 0.80, adj = 0
  )
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

format_top_banner <- function(network_label = NULL, vars_label = NULL, source_label = NULL, max_len = 180) {
  parts <- character(0)
  if (!is.null(network_label) && nzchar(network_label)) {
    parts <- c(parts, sprintf("Water Network: %s", network_label))
  }
  if (!is.null(source_label) && nzchar(source_label)) {
    parts <- c(parts, sprintf("Data Source: %s", source_label))
  }
  if (!is.null(vars_label) && nzchar(vars_label)) {
    parts <- c(parts, sprintf("Fixed-Effect Regression Covariates: %s", vars_label))
  }
  raw_banner <- paste(parts, collapse = "   |   ")
  wrap_title(raw_banner, max_len = max_len)
}

get_banner_oma <- function(banner) {
  if (is.null(banner) || !nzchar(banner)) return(c(0, 0, 0, 0))
  n_b <- length(strsplit(banner, "\n", fixed = TRUE)[[1]])
  c(0, 0, if (n_b > 1) 3.2 else 2.2, 0)
}

draw_top_banner <- function(banner, cex = 0.88, col = "#2c3e50") {
  if (is.null(banner) || !nzchar(banner)) return(invisible(NULL))
  b_splits <- strsplit(banner, "\n", fixed = TRUE)[[1]]
  n_b <- length(b_splits)
  if (n_b == 1) {
    mtext(b_splits[1], side = 3, line = 0.6, outer = TRUE, cex = cex, font = 3, col = col)
  } else {
    for (i in seq_len(n_b)) {
      b_line <- 0.35 + (n_b - i) * 0.90
      mtext(b_splits[i], side = 3, line = b_line, outer = TRUE, cex = cex * 0.95, font = 3, col = col)
    }
  }
}

# Figure 1: Full England Stream Network Map
plot_fig1_england_map <- function(eng_pred, det_code, det_name, det_unit, trans_obj, out_file,
                                  reverse_color = FALSE, sc_val = NULL, sc_se = NULL,
                                  vars_label = NULL, network_label = NULL, source_label = NULL,
                                  norm_mode = "none") {
  tr       <- as_transformer(trans_obj)
  nm       <- get_norm_meta(norm_mode)
  if (nm$mode != "none" && !(nm$pred_col %in% names(eng_pred))) {
    eng_pred <- add_normalised_columns(eng_pred)
  }
  pal_use  <- get_norm_pal(nm$mode, reverse_color)
  has_vars <- !is.null(vars_label) && nzchar(vars_label)
  has_net  <- !is.null(network_label) && nzchar(network_label)
  banner   <- format_top_banner(network_label, vars_label, source_label)
  open_png_device(out_file, width = 2600, height = 1450, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  if (is.null(sc_val)) {
    sc_val <- if (nm$mode == "none") {
      compute_color_scale(eng_pred$pred_model_sc, tr)
    } else {
      compute_norm_color_scale(eng_pred[[nm$pred_col]], nm$mode)
    }
  }
  if (is.null(sc_se)) sc_se <- compute_se_scale(eng_pred$se_model_sc)
  n_obs_eng <- sum(eng_pred$has_obs == 1)
  n_tot_eng <- nrow(eng_pred)
  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  unit_disp <- if (nm$mode == "none") det_unit else nm$unit_label
  norm_tag  <- if (nm$mode == "none") "" else sprintf(" [%s]", nm$short_label)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, unit_disp)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag   <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  sub_a <- if (has_vars) {
    sprintf("S3N across %d catchments%s | Covariates: %s",
            length(unique(eng_pred$CatchID)), net_tag, vars_label)
  } else {
    sprintf("Scalable Spatial Stream Network (S3N) across %d catchments%s",
            length(unique(eng_pred$CatchID)), net_tag)
  }

  plot_stream_network(
    eng_pred, nm$pred_col, vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. S3N Predicted %s%s Across England (%s Reaches)", det_short, norm_tag, format(n_tot_eng, big.mark = ",")),
    sub_title  = sub_a,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )

  se_scale_txt <- if (tr$type == "none") det_unit else sprintf("%s scale", tr$short_label)
  se_cbar_txt  <- if (tr$type == "none") sprintf("Kriging SE\n(%s)", det_unit) else sprintf("Kriging SE\n(%s scale)", tr$short_label)

  plot_stream_network(
    eng_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. S3N Prediction Uncertainty & %s Stations (n = %s)", src_tag, format(n_obs_eng, big.mark = ",")),
    sub_title  = sprintf("Kriging standard error (%s) with observed %s stations overlaid%s",
                         se_scale_txt, src_tag, net_tag),
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = se_cbar_txt,
    show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6
  )
  obs_eng <- eng_pred[eng_pred$has_obs == 1, ]
  points(obs_eng$mid_X, obs_eng$mid_Y, pch = 21,
         bg = map_to_colors(obs_eng[[nm$obs_col]], sc_val$vmin, sc_val$vmax, pal_use),
         col = "#111111", cex = 0.65, lwd = 0.5)
  legend("topleft", legend = sprintf("Observed %s Station (colored by %s)", src_tag, det_short),
         pch = 21, pt.bg = "#e67e22", col = "#111111", pt.cex = 1.0, bty = "n", cex = 0.80)
  draw_top_banner(banner)
}

# Generic 2-Panel River Basin High-Resolution S3N Map (used for Figure 2 Humber & Figure 5+ Custom Basins)
plot_basin_map <- function(basin_pred, basin_title, domain_label, cv_tab, sc_val, det_code, det_name, det_unit, trans_obj, out_file,
                           reverse_color = FALSE, sc_se = NULL, vars_label = NULL, network_label = NULL, source_label = NULL,
                           norm_mode = "none") {
  tr         <- as_transformer(trans_obj)
  nm         <- get_norm_meta(norm_mode)
  if (nm$mode != "none" && !(nm$pred_col %in% names(basin_pred))) {
    basin_pred <- add_normalised_columns(basin_pred)
  }
  pal_use    <- get_norm_pal(nm$mode, reverse_color)
  has_vars   <- !is.null(vars_label) && nzchar(vars_label)
  has_net    <- !is.null(network_label) && nzchar(network_label)
  banner     <- format_top_banner(network_label, vars_label, source_label)
  open_png_device(out_file, width = 2600, height = 1350, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  if (is.null(sc_se)) sc_se <- compute_se_scale(basin_pred$se_model_sc)
  n_obs_basin <- sum(basin_pred$has_obs == 1)
  det_short   <- shorten_det_name(det_name)
  unit_disp   <- if (nm$mode == "none") det_unit else nm$unit_label
  norm_tag    <- if (nm$mode == "none") "" else sprintf(" [%s]", nm$short_label)
  net_tag     <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag     <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  cv_row      <- subset(cv_tab, (Domain == domain_label | grepl(paste0("^", basin_title), Domain)) &
                                Model == "S3N Combined (Tail-Up + Tail-Down)")
  sub_a <- if (nrow(cv_row) > 0) {
    if (has_vars) {
      sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | n = %d sites%s | Covariates: %s",
              cv_row$R2[1], cv_row$RMSE[1], n_obs_basin, net_tag, vars_label)
    } else {
      sprintf("5-Fold CV R2 = %.3f (RMSE = %.3f) | %d observed %s stations%s",
              cv_row$R2[1], cv_row$RMSE[1], n_obs_basin, src_tag, net_tag)
    }
  } else {
    if (has_vars) {
      sprintf("England S3N subset (n = %d sites)%s | Covariates: %s", n_obs_basin, net_tag, vars_label)
    } else {
      sprintf("England S3N model subset | %d observed %s stations in %s%s", n_obs_basin, src_tag, basin_title, net_tag)
    }
  }

  plot_stream_network(
    basin_pred, nm$pred_col, vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. %s: S3N Predicted %s%s (%s Reaches)",
                         basin_title, det_short, norm_tag, format(nrow(basin_pred), big.mark = ",")),
    sub_title  = sub_a,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels,
    cbar_title = sprintf("%s\n(%s)", sub(" \\[.*\\]$", "", det_short), unit_disp),
    show_obs = TRUE, obs_col = nm$obs_col, base_lwd = 0.35, max_lwd = 2.8
  )
  legend("topleft", legend = sprintf("Observed %s Station (same color scale)", src_tag),
         pch = 21, pt.bg = "#f4d03f", col = "#111111", pt.cex = 1.1, bty = "n", cex = 0.82)

  se_cbar_txt <- if (tr$type == "none") sprintf("Kriging SE\n(%s)", det_unit) else sprintf("Kriging SE\n(%s scale)", tr$short_label)

  plot_stream_network(
    basin_pred, "se_model_sc", vmin = sc_se$vmin, vmax = sc_se$vmax, pal = uncert_cols,
    main_title = sprintf("B. %s: S3N Kriging Standard Error (%s)", basin_title, det_short),
    sub_title  = sprintf("Lower along monitored flow-connected paths, higher in headwaters%s", net_tag),
    cbar_ticks = sc_se$ticks, cbar_labels = sc_se$labels,
    cbar_title = se_cbar_txt,
    show_obs = FALSE, base_lwd = 0.35, max_lwd = 2.8
  )
  obs_b <- basin_pred[basin_pred$has_obs == 1, ]
  if (nrow(obs_b) > 0) {
    points(obs_b$mid_X, obs_b$mid_Y, pch = 3, col = "#c0392b", cex = 0.55, lwd = 0.8)
  }
  legend("topleft", legend = sprintf("%s Monitoring Site Location (n = %d)", src_tag, n_obs_basin),
         pch = 3, col = "#c0392b", pt.cex = 0.9, bty = "n", cex = 0.82)
  draw_top_banner(banner)
}

# Figure 2: Humber River Basin High-Resolution S3N Map
plot_fig2_humber_map <- function(hum_pred, cv_tab, sc_val, det_code, det_name, det_unit, trans_obj, out_file,
                                 reverse_color = FALSE, sc_se = NULL, vars_label = NULL, network_label = NULL, source_label = NULL,
                                 norm_mode = "none") {
  plot_basin_map(hum_pred, "Humber Basin", "Humber Basin (CatchID=80)",
                 cv_tab, sc_val, det_code, det_name, det_unit, trans_obj, out_file,
                 reverse_color = reverse_color, sc_se = sc_se,
                 vars_label = vars_label, network_label = network_label, source_label = source_label,
                 norm_mode = norm_mode)
}

# Figure 3: Thames & Severn River Basins S3N Maps
plot_fig3_thames_severn_map <- function(tha_pred, sev_pred, sc_val, det_code, det_name, det_unit, trans_obj, out_file,
                                        reverse_color = FALSE, vars_label = NULL, network_label = NULL, source_label = NULL,
                                        norm_mode = "none") {
  nm       <- get_norm_meta(norm_mode)
  if (nm$mode != "none") {
    if (!(nm$pred_col %in% names(tha_pred))) tha_pred <- add_normalised_columns(tha_pred)
    if (!(nm$pred_col %in% names(sev_pred))) sev_pred <- add_normalised_columns(sev_pred)
  }
  pal_use  <- get_norm_pal(nm$mode, reverse_color)
  has_vars <- !is.null(vars_label) && nzchar(vars_label)
  has_net  <- !is.null(network_label) && nzchar(network_label)
  banner   <- format_top_banner(network_label, vars_label, source_label)
  open_png_device(out_file, width = 2600, height = 1250, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  det_short <- shorten_det_name(det_name)
  det_cbar  <- sub(" \\[.*\\]$", "", det_short)
  unit_disp <- if (nm$mode == "none") det_unit else nm$unit_label
  norm_tag  <- if (nm$mode == "none") "" else sprintf(" [%s]", nm$short_label)
  cbar_lbl  <- sprintf("%s\n(%s)", det_cbar, unit_disp)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  src_tag   <- if (!is.null(source_label) && nzchar(source_label)) format_source_abbr(source_label) else "RSN"
  sub_tha <- if (has_vars) {
    sprintf("n = %d %s sites%s | S3N Tail-Up | Covariates: %s", sum(tha_pred$has_obs == 1), src_tag, net_tag, vars_label)
  } else {
    sprintf("n = %d observed %s stations%s | Flow-connected Tail-Up predictions",
            sum(tha_pred$has_obs == 1), src_tag, net_tag)
  }
  sub_sev <- if (has_vars) {
    sprintf("n = %d %s sites%s | S3N Combined | Covariates: %s", sum(sev_pred$has_obs == 1), src_tag, net_tag, vars_label)
  } else {
    sprintf("n = %d observed %s stations%s | Combined Tail-Up + Tail-Down predictions",
            sum(sev_pred$has_obs == 1), src_tag, net_tag)
  }

  plot_stream_network(
    tha_pred, nm$tu_col, vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("A. Thames Basin: S3N Tail-Up Predicted %s%s (%s Reaches)",
                         det_short, norm_tag, format(nrow(tha_pred), big.mark = ",")),
    sub_title  = sub_tha,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = nm$obs_col, base_lwd = 0.40, max_lwd = 2.8
  )

  plot_stream_network(
    sev_pred, nm$pred_col, vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
    main_title = sprintf("B. Severn Basin: S3N Predicted %s%s (%s Reaches)",
                         det_short, norm_tag, format(nrow(sev_pred), big.mark = ",")),
    sub_title  = sub_sev,
    cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
    show_obs = TRUE, obs_col = nm$obs_col, base_lwd = 0.40, max_lwd = 2.8
  )
  draw_top_banner(banner)
}

# Figure 4: 4-Panel Model Diagnostics, Cross-Validation & Parameter CIs
plot_fig4_diagnostics <- function(eng_pred, eng_cv, hum_cv, cv_tab, param_tab, det_code, det_name, det_unit, trans_obj, out_file,
                                  extra_cv_list = list(), lims_a = NULL, ylim_d = NULL,
                                  vars_label = NULL, network_label = NULL, source_label = NULL) {
  tr       <- as_transformer(trans_obj)
  banner   <- format_top_banner(network_label, vars_label, source_label)

  open_png_device(out_file, width = 2400, height = 1800, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.5, 5.2, 3.8, 1.8),
      oma = get_banner_oma(banner), bg = "white")

  det_short <- shorten_det_name(det_name)
  sc_lbl <- if (tr$type == "none") {
    sprintf("%s (%s)", det_short, det_unit)
  } else {
    sprintf("%s(%s, %s)", tr$fn_prefix, det_short, det_unit)
  }

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
       main = wrap_title("A. Out-of-Sample 5-Fold Cross-Validation Accuracy", 60))
  grid(col = "#e0e0e0", lty = "dotted")
  if (!is.null(hum_cv) && nrow(hum_cv) > 0) {
    points(hum_cv$obs_y, hum_cv$cv_comb, pch = 21, bg = rgb(0.90, 0.49, 0.13, 0.75), col = "#222222", cex = 0.95, lwd = 0.6)
  }

  eng_dom <- if (any(grepl("All England", cv_tab$Domain))) {
    unique(cv_tab$Domain[grep("All England", cv_tab$Domain)])[1]
  } else {
    "All England RSN Network"
  }
  eng_row <- subset(cv_tab, Domain == eng_dom & Model == "S3N Combined (Tail-Up + Tail-Down)")
  if (nrow(eng_row) == 0) {
    eng_row <- subset(cv_tab, Model == "S3N Combined (Tail-Up + Tail-Down)")
  }
  hum_row <- subset(cv_tab, grepl("^Humber Basin", Domain) & Model == "S3N Combined (Tail-Up + Tail-Down)")

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
      pts_k <- ecv$cv_pts
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
  model_order <- c("OLS (Non-Spatial)", "S3N Tail-Down (Stream Distance)",
                   "S3N Tail-Up (Flow-Connected + AFV)", "S3N Combined (Tail-Up + Tail-Down)")
  r2_mat <- vapply(avail_doms, function(d) {
    sub <- cv_tab[cv_tab$Domain == d, ]
    vals <- vapply(model_order, function(mn) {
      hit <- sub$R2[sub$Model == mn]
      if (length(hit) > 0 && is.finite(hit[1])) max(0, hit[1]) else NA_real_
    }, numeric(1))
    unname(vals)
  }, numeric(4))
  if (is.null(dim(r2_mat))) r2_mat <- matrix(r2_mat, nrow = 4, ncol = length(avail_doms))
  short_map <- c(
    "Humber Basin (CatchID=80)"   = "Humber Basin",
    "Severn Basin (CatchID=2786)" = "Severn Basin",
    "Thames Basin (CatchID=3055)" = "Thames Basin"
  )
  short_map[eng_dom] <- sprintf("All England\n(n=%d)", eng_row$n_obs[1])
  colnames(r2_mat) <- ifelse(avail_doms %in% names(short_map), short_map[avail_doms],
                             sub(" \\((CatchID|RBD)=.*\\)$", "", avail_doms))
  ymax_b <- max(0.88, max(r2_mat, na.rm = TRUE) * 1.18)
  cex_names_b <- if (ncol(r2_mat) > 4) 0.78 else 0.90
  cex_txt_b   <- if (ncol(r2_mat) > 4) 0.64 else 0.72
  bp <- barplot(r2_mat, beside = TRUE, ylim = c(0, ymax_b),
                col = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"),
                border = "#222222", cex.names = cex_names_b, ylab = "5-Fold Cross-Validation R2",
                main = wrap_title("B. Out-of-Sample Variance Explained (R2) by Model", 60))
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  text(bp, r2_mat + 0.03 * ymax_b, labels = sprintf("%.2f", r2_mat), cex = cex_txt_b, font = 2)
  legend("topright",
         legend = c("OLS (Non-Spatial)", "S3N Tail-Down (Stream)", "S3N Tail-Up (FC + AFV)", "S3N Combined (TU + TD)"),
         fill = c("#95a5a6", "#3498db", "#2ecc71", "#e67e22"), border = "#222222", bty = "n", cex = 0.80)

  # Panel C: Forest Plot of Fixed-Effect Estimates & 95% Bootstrap CIs (All England S3N)
  eng_p_dom <- if (!is.null(param_tab) && any(grepl("All England", param_tab$Domain))) {
    unique(param_tab$Domain[grep("All England", param_tab$Domain)])[1]
  } else {
    "All England RSN Network"
  }
  non_cov_params <- c("Intercept", "sigma.sq", "tau.sq", "phi", "tau.sq (microscale)", "mean_meas_variance")
  eng_td_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Down (Stream)" &
                     !(Parameter %in% non_cov_params))
  eng_tu_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Up (FC + AFV)" &
                     !(Parameter %in% non_cov_params))
  if (nrow(eng_td_p) == 0) {
    eng_td_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Down (Stream)" &
                       !(Parameter %in% setdiff(non_cov_params, "Intercept")))
    eng_tu_p <- subset(param_tab, Domain == eng_p_dom & Model == "S3N Tail-Up (FC + AFV)" &
                       !(Parameter %in% setdiff(non_cov_params, "Intercept")))
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
  par(mar = c(4.5, 10.2, 3.8, 1.5))
  xr_ci <- range(c(eng_tu_p$CI_Lower, eng_tu_p$CI_Upper, eng_td_p$CI_Lower, eng_td_p$CI_Upper, -0.2, 0.2), na.rm = TRUE)
  xr_pad <- 0.12 * diff(xr_ci)
  plot(NA, xlim = c(xr_ci[1] - xr_pad, xr_ci[2] + xr_pad), ylim = c(0.4, np + 0.6),
       yaxt = "n", xlab = sprintf("Effect Size on %s [95%% Bootstrap CI]", sc_lbl), ylab = "",
       main = wrap_title(sprintf("C. Environmental Drivers of %s (England S3N)", det_short), 60))
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
                  main = wrap_title(sprintf("D. Predicted %s Across %s Reaches by Regime",
                                           det_short, format(nrow(eng_pred), big.mark = ",")), 60))
  grid(nx = NA, ny = NULL, col = "#e0e0e0", lty = "dotted")
  abline(v = 1.5, col = "#7f8c8d", lty = 2, lwd = 1.3)
  medians <- sapply(box_list, median, na.rm = TRUE)
  yr_d <- par("usr")[4] - par("usr")[3]
  y_txt <- pmin(par("usr")[4] - 0.07 * yr_d, bp_d$stats[4, ] + 0.07 * yr_d)
  text(1:4, y_txt, labels = sprintf("Median:\n%s %s", formatC(medians, format = "fg", digits = 3), det_unit),
       font = 2, cex = 0.78, col = "#111111")
  draw_top_banner(banner)
}

# Figure 5 (Temporal Multi-Panel Grid): All England S3N Comparison Across Seasons, Years, or Months
plot_fig5_temporal_england_grid <- function(period_items, by_mode, sc_val, det_code, det_name, det_unit,
                                            trans_obj, out_file, reverse_color = FALSE,
                                            vars_label = NULL, network_label = NULL, source_label = NULL,
                                            norm_mode = "none") {
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

  nm        <- get_norm_meta(norm_mode)
  pal_use   <- get_norm_pal(nm$mode, reverse_color)
  det_short <- sub(" \\[.*\\]$", "", shorten_det_name(det_name))
  unit_disp <- if (nm$mode == "none") det_unit else nm$unit_label
  cbar_lbl  <- sprintf("%s\n(%s)", det_short, unit_disp)
  mode_cap  <- paste0(toupper(substr(by_mode, 1, 1)), tolower(substr(by_mode, 2, nchar(by_mode))))

  open_png_device(out_file, width = w_px, height = h_px, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(nrow_g, ncol_g), mar = c(1.0, 1.0, 3.6, 1.0), oma = c(0.5, 0.5, 3.8, 0.5), bg = "white")

  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  for (i in seq_along(period_items)) {
    item   <- period_items[[i]]
    eng_df <- item$eng_pred
    if (nm$mode != "none" && !(nm$pred_col %in% names(eng_df))) {
      eng_df <- add_normalised_columns(eng_df)
    }
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
      cvt <- item$cv_tab
      erow <- subset(cvt, grepl("All England", Domain) & Model == "S3N Combined (Tail-Up + Tail-Down)")
      if (nrow(erow) > 0 && !is.na(erow$R2[1])) {
        r2_str <- sprintf(" | CV R2 = %.3f", erow$R2[1])
      }
    }

    main_t <- sprintf("%s. %s — Predicted %s", LETTERS[i], pan_lbl, det_short)
    sub_t  <- sprintf("n = %s obs reaches%s | Network Median = %s %s",
                      format(n_obs_i, big.mark = ","), r2_str,
                      formatC(med_raw, format = "fg", digits = 3), det_unit)

    max_m_t <- if (ncol_g >= 3) 48 else 65
    max_s_t <- if (ncol_g >= 3) 60 else 90
    plot_stream_network(
      eng_df, nm$pred_col, vmin = sc_val$vmin, vmax = sc_val$vmax, pal = pal_use,
      main_title = main_t,
      sub_title  = sub_t,
      cbar_ticks = sc_val$ticks, cbar_labels = sc_val$labels, cbar_title = cbar_lbl,
      show_obs = FALSE, base_lwd = 0.18, max_lwd = 1.6,
      max_main_chars = max_m_t, max_sub_chars = max_s_t
    )
  }

  n_reaches_tot <- nrow(period_items[[1]]$eng_pred)
  n_catch_tot   <- length(unique(period_items[[1]]$eng_pred$CatchID))
  has_vars      <- !is.null(vars_label) && nzchar(vars_label)
  has_net       <- !is.null(network_label) && nzchar(network_label)
  has_src       <- !is.null(source_label) && nzchar(source_label)
  src_tag       <- if (has_src) format_source_abbr(source_label) else ""
  sub_outer <- sprintf("Unified color scale across all %d %s panels (%s reaches across %d catchments)%s%s%s",
                       n_pan, tolower(by_mode), format(n_reaches_tot, big.mark = ","), n_catch_tot,
                       if (has_src) paste0(" | Source: ", src_tag) else "",
                       if (has_net) paste0(" | Network: ", network_label) else "",
                       if (has_vars) paste0(" | Covariates: ", vars_label) else "")
  mtext(wrap_title(sprintf("All England S3N Stream Network Comparison by %s: %s (%s)",
                           mode_cap, det_short, det_unit), 180),
        side = 3, line = 1.6, outer = TRUE, cex = 1.20, font = 2, col = "#111111")
  mtext(wrap_title(sub_outer, 190), side = 3, line = 0.2, outer = TRUE, cex = 0.90, col = "#444444")
}

# ==============================================================================
# Regulatory Exceedance Probability & 5-Tier WFD Management Classification Plots
# ==============================================================================

COMPLIANCE_TIER_LEVELS <- c(
  "Confidently Compliant",
  "Likely Compliant",
  "Uncertain / Priority Monitoring",
  "Likely Exceeding",
  "Confidently Exceeding"
)

COMPLIANCE_TIER_COLORS <- c(
  "Confidently Compliant"           = "#1a9850",
  "Likely Compliant"                = "#91cf60",
  "Uncertain / Priority Monitoring" = "#f4d03f",
  "Likely Exceeding"                = "#f46d43",
  "Confidently Exceeding"           = "#d73027"
)

EXCEED_PROB_PAL <- colorRampPalette(c("#1a9850", "#91cf60", "#f4d03f", "#f46d43", "#d73027"))(256)

plot_compliance_tier_network <- function(df, main_title, sub_title,
                                         base_lwd = 0.35, max_lwd = 2.6,
                                         show_obs = TRUE, show_legend = TRUE,
                                         max_main_chars = 65, max_sub_chars = 90) {
  ord <- order(df$Shreve)
  d   <- df[ord, ]

  xr <- range(c(d$start_X, d$end_X), na.rm = TRUE)
  yr <- range(c(d$start_Y, d$end_Y), na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)

  right_pad <- if (show_legend) 0.04 * dx else 0.03 * dx
  plot(NA, xlim = c(xr[1] - 0.03 * dx, xr[2] + right_pad),
       ylim = c(yr[1] - 0.03 * dy, yr[2] + 0.06 * dy),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#d0d4d9", lwd = 1)

  lwds <- base_lwd + (max_lwd - base_lwd) * pmin(1, log10(pmax(1, d$Shreve)) / log10(max(10, max(d$Shreve))))
  tier_chr <- as.character(d$compliance_tier)
  cols <- unname(COMPLIANCE_TIER_COLORS[tier_chr])
  cols[is.na(cols)] <- "#cccccc"
  segments(d$start_X, d$start_Y, d$end_X, d$end_Y, col = cols, lwd = lwds)

  if (show_obs) {
    obs_d <- d[d$has_obs == 1, ]
    if (nrow(obs_d) > 0L) {
      points(obs_d$mid_X, obs_d$mid_Y, pch = 3, col = "#222222", cex = 0.45, lwd = 0.7)
    }
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

  if (show_legend) {
    tab <- table(factor(tier_chr, levels = COMPLIANCE_TIER_LEVELS))
    pct <- 100 * as.numeric(tab) / max(1L, nrow(d))
    leg_labels <- sprintf("%s (%s, %.1f%%)",
                          COMPLIANCE_TIER_LEVELS,
                          format(as.integer(tab), big.mark = ","),
                          pct)
    legend("topleft",
           legend = leg_labels,
           col    = unname(COMPLIANCE_TIER_COLORS[COMPLIANCE_TIER_LEVELS]),
           lwd    = 3.5,
           bg     = rgb(1, 1, 1, 0.90),
           box.col = "#cccccc",
           cex    = 0.76,
           title  = "Management Classification",
           title.cex = 0.80)
  }

  draw_plot_header(
    main_title = main_title, sub_title = sub_title,
    max_main = max_main_chars, max_sub = max_sub_chars,
    cex_main = 0.96, cex_sub = 0.80, adj = 0
  )
}

plot_fig_exceedance_england <- function(eng_pred, det_code, det_name, det_unit, out_file,
                                        threshold_val = NA_real_, threshold_dir = "above",
                                        vars_label = NULL, network_label = NULL, source_label = NULL) {
  if (!("compliance_tier" %in% names(eng_pred)) || !("exceed_prob" %in% names(eng_pred))) {
    return(invisible(NULL))
  }
  det_short <- shorten_det_name(det_name)
  has_net   <- !is.null(network_label) && nzchar(network_label)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  banner    <- format_top_banner(network_label, vars_label, source_label)
  dir_sym   <- if (identical(tolower(as.character(threshold_dir)), "below")) "<" else ">"
  thresh_txt <- if (is.finite(threshold_val)) {
    sprintf("Threshold: Y %s %s %s", dir_sym, formatC(threshold_val, format = "fg", digits = 4), det_unit)
  } else {
    "Regulatory Threshold Classification"
  }

  open_png_device(out_file, width = 2600, height = 1350, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  plot_compliance_tier_network(
    eng_pred,
    main_title = sprintf("A. All England: Assigned Management Classification — %s", det_short),
    sub_title  = sprintf("%s | %s reaches%s",
                         thresh_txt, format(nrow(eng_pred), big.mark = ","), net_tag),
    base_lwd = 0.25, max_lwd = 2.0, show_obs = FALSE, show_legend = TRUE
  )

  plot_stream_network(
    eng_pred, "exceed_prob", vmin = 0, vmax = 1, pal = EXCEED_PROB_PAL,
    main_title = sprintf("B. All England: Exceedance Probability P(Y %s T) — %s", dir_sym, det_short),
    sub_title  = sprintf("Gaussian predictive probability from S3N kriging mean & SE%s", net_tag),
    cbar_ticks = c(0, 0.05, 0.20, 0.50, 0.80, 0.95, 1.0),
    cbar_labels = c("0%", "5%", "20%", "50%", "80%", "95%", "100%"),
    cbar_title = sprintf("P(Y %s T)\n(%%)", dir_sym),
    show_obs = FALSE, base_lwd = 0.25, max_lwd = 2.0
  )
  draw_top_banner(banner)
}

plot_fig_exceedance_basins <- function(hum_pred, tha_pred, sev_pred, det_code, det_name, det_unit, out_file,
                                       threshold_val = NA_real_, threshold_dir = "above",
                                       vars_label = NULL, network_label = NULL, source_label = NULL) {
  basins <- list(
    list(df = hum_pred, title = "A. Humber Basin"),
    list(df = tha_pred, title = "B. Thames Basin"),
    list(df = sev_pred, title = "C. Severn Basin")
  )
  basins <- Filter(function(x) !is.null(x$df) && nrow(x$df) > 0L && ("compliance_tier" %in% names(x$df)), basins)
  if (length(basins) == 0L) return(invisible(NULL))

  det_short  <- shorten_det_name(det_name)
  banner     <- format_top_banner(network_label, vars_label, source_label)
  dir_sym    <- if (identical(tolower(as.character(threshold_dir)), "below")) "<" else ">"
  thresh_txt <- if (is.finite(threshold_val)) {
    sprintf("Threshold: Y %s %s %s", dir_sym, formatC(threshold_val, format = "fg", digits = 4), det_unit)
  } else {
    "Management Classification"
  }

  open_png_device(out_file, width = 3000, height = 1150, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, length(basins)), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  for (b in basins) {
    plot_compliance_tier_network(
      b$df,
      main_title = sprintf("%s: Management Classification (%s)", b$title, det_short),
      sub_title  = sprintf("%s | %s reaches (n = %d obs)",
                           thresh_txt, format(nrow(b$df), big.mark = ","), sum(b$df$has_obs == 1)),
      base_lwd = 0.35, max_lwd = 2.4, show_obs = TRUE, show_legend = TRUE,
      max_main_chars = 52, max_sub_chars = 68
    )
  }
  draw_top_banner(banner)
}

plot_fig_exceedance_basin <- function(basin_pred, basin_title, det_code, det_name, det_unit, out_file,
                                      threshold_val = NA_real_, threshold_dir = "above",
                                      vars_label = NULL, network_label = NULL, source_label = NULL) {
  if (is.null(basin_pred) || !("compliance_tier" %in% names(basin_pred)) || !("exceed_prob" %in% names(basin_pred))) {
    return(invisible(NULL))
  }
  det_short  <- shorten_det_name(det_name)
  has_net    <- !is.null(network_label) && nzchar(network_label)
  net_tag    <- if (has_net) sprintf(" [%s]", network_label) else ""
  banner     <- format_top_banner(network_label, vars_label, source_label)
  dir_sym    <- if (identical(tolower(as.character(threshold_dir)), "below")) "<" else ">"
  thresh_txt <- if (is.finite(threshold_val)) {
    sprintf("Threshold: Y %s %s %s", dir_sym, formatC(threshold_val, format = "fg", digits = 4), det_unit)
  } else {
    "Regulatory Threshold Classification"
  }

  open_png_device(out_file, width = 2600, height = 1350, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  plot_compliance_tier_network(
    basin_pred,
    main_title = sprintf("A. %s: Assigned Management Classification (%s)", basin_title, det_short),
    sub_title  = sprintf("%s | %s reaches%s",
                         thresh_txt, format(nrow(basin_pred), big.mark = ","), net_tag),
    base_lwd = 0.35, max_lwd = 2.8, show_obs = TRUE, show_legend = TRUE
  )

  plot_stream_network(
    basin_pred, "exceed_prob", vmin = 0, vmax = 1, pal = EXCEED_PROB_PAL,
    main_title = sprintf("B. %s: Exceedance Probability P(Y %s T) (%s)", basin_title, dir_sym, det_short),
    sub_title  = sprintf("Gaussian predictive probability from S3N kriging mean & SE%s", net_tag),
    cbar_ticks = c(0, 0.05, 0.20, 0.50, 0.80, 0.95, 1.0),
    cbar_labels = c("0%", "5%", "20%", "50%", "80%", "95%", "100%"),
    cbar_title = sprintf("P(Y %s T)\n(%%)", dir_sym),
    show_obs = FALSE, base_lwd = 0.35, max_lwd = 2.8
  )
  obs_b <- basin_pred[basin_pred$has_obs == 1, ]
  if (nrow(obs_b) > 0L) {
    points(obs_b$mid_X, obs_b$mid_Y, pch = 3, col = "#222222", cex = 0.55, lwd = 0.8)
  }
  draw_top_banner(banner)
}

# ==============================================================================
# Figure 6: Spatio-Temporal Trend Mapping (--by=year|season|month)
# ==============================================================================

TREND_TIER_LEVELS <- c(
  "Significant Decrease (p < 0.05)",
  "Moderate Decrease (p < 0.20)",
  "Stable / No Clear Trend",
  "Moderate Increase (p < 0.20)",
  "Significant Increase (p < 0.05)"
)

TREND_TIER_COLORS <- c(
  "Significant Decrease (p < 0.05)" = "#2166ac",
  "Moderate Decrease (p < 0.20)"    = "#67a9cf",
  "Stable / No Clear Trend"         = "#cfd4da",
  "Moderate Increase (p < 0.20)"    = "#ef8a62",
  "Significant Increase (p < 0.05)" = "#b2182b"
)

TREND_DIV_PAL <- colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))(256)

plot_fig6_spatiotemporal_trends <- function(trend_df, by_mode, period_labels, det_code, det_name, det_unit,
                                            trans_obj, out_file, borrowed = FALSE, rho_t = NA_real_,
                                            vars_label = NULL, network_label = NULL, source_label = NULL) {
  if (is.null(trend_df) || nrow(trend_df) == 0L) return(invisible(NULL))
  tr        <- as_transformer(trans_obj)
  det_short <- shorten_det_name(det_name)
  has_net   <- !is.null(network_label) && nzchar(network_label)
  net_tag   <- if (has_net) sprintf(" [%s]", network_label) else ""
  banner    <- format_top_banner(network_label, vars_label, source_label)
  mode_cap  <- paste0(toupper(substr(by_mode, 1, 1)), substr(by_mode, 2, nchar(by_mode)))
  st_tag    <- if (isTRUE(borrowed) && is.finite(rho_t)) {
    sprintf("Dynamic Temporal Borrowing (rho_t = %.2f)", rho_t)
  } else {
    "Precision-Weighted WLS Trend"
  }

  slp <- as.numeric(trend_df$trend_slope_sc)
  lim <- max(1e-4, as.numeric(quantile(abs(slp), probs = 0.985, na.rm = TRUE)))
  ticks_a <- seq(-lim, lim, length.out = 5L)
  labs_a  <- sprintf("%+.3f", ticks_a)

  open_png_device(out_file, width = 2600, height = 1350, res = 200)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1.2, 1.2, 4.0, 1.2),
      oma = get_banner_oma(banner), bg = "white")

  plot_stream_network(
    trend_df, "trend_slope_sc", vmin = -lim, vmax = lim, pal = TREND_DIV_PAL,
    main_title = sprintf("A. All England: %s-to-%s Trend Slope (%s)", mode_cap, mode_cap, det_short),
    sub_title  = sprintf("%s -> %s (%s)%s", period_labels[1], tail(period_labels, 1), st_tag, net_tag),
    cbar_ticks = ticks_a, cbar_labels = labs_a,
    cbar_title = sprintf("Slope / %s\n(%s scale)", tolower(by_mode), tr$short_label),
    show_obs = FALSE, base_lwd = 0.25, max_lwd = 2.0
  )

  # Panel B: Trend Significance Classification
  ord <- order(trend_df$Shreve)
  d   <- trend_df[ord, ]
  xr  <- range(c(d$start_X, d$end_X), na.rm = TRUE)
  yr  <- range(c(d$start_Y, d$end_Y), na.rm = TRUE)
  dx  <- diff(xr)
  dy  <- diff(yr)
  plot(NA, xlim = c(xr[1] - 0.03 * dx, xr[2] + 0.04 * dx),
       ylim = c(yr[1] - 0.03 * dy, yr[2] + 0.06 * dy),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#d0d4d9", lwd = 1)
  lwds <- 0.25 + (2.0 - 0.25) * pmin(1, log10(pmax(1, d$Shreve)) / log10(max(10, max(d$Shreve))))
  tier_chr <- as.character(d$trend_tier)
  cols <- unname(TREND_TIER_COLORS[tier_chr])
  cols[is.na(cols)] <- "#cfd4da"
  segments(d$start_X, d$start_Y, d$end_X, d$end_Y, col = cols, lwd = lwds)

  tab <- table(factor(tier_chr, levels = TREND_TIER_LEVELS))
  pct <- 100 * as.numeric(tab) / max(1L, nrow(d))
  leg_labels <- sprintf("%s (%s, %.1f%%)", TREND_TIER_LEVELS, format(as.integer(tab), big.mark = ","), pct)
  legend("topleft", legend = leg_labels, col = unname(TREND_TIER_COLORS[TREND_TIER_LEVELS]),
         lwd = 3.5, bg = rgb(1, 1, 1, 0.90), box.col = "#cccccc", cex = 0.75,
         title = "Trend Significance (Wald Z-test)", title.cex = 0.80)
  draw_plot_header(
    main_title = sprintf("B. All England: Trend Direction & Significance (%s)", det_short),
    sub_title  = sprintf("Two-sided Wald Z-test across %d %ss (%s reaches)%s",
                         length(period_labels), tolower(by_mode), format(nrow(d), big.mark = ","), net_tag),
    max_main = 65, max_sub = 90, cex_main = 0.96, cex_sub = 0.80, adj = 0
  )
  draw_top_banner(banner)
}

plot_fig7_spatiotemporal_transitions <- function(transitions, by_mode, det_code, det_name, det_unit,
                                                 trans_obj, out_file, borrowed = FALSE, rho_t = NA_real_,
                                                 vars_label = NULL, network_label = NULL, source_label = NULL) {
  if (is.null(transitions) || length(transitions) == 0L) return(invisible(NULL))
  tr        <- as_transformer(trans_obj)
  det_short <- shorten_det_name(det_name)
  banner    <- format_top_banner(network_label, vars_label, source_label)
  mode_cap  <- paste0(toupper(substr(by_mode, 1, 1)), substr(by_mode, 2, nchar(by_mode)))
  st_tag    <- if (isTRUE(borrowed) && is.finite(rho_t)) {
    sprintf("Dynamic Temporal Borrowing (rho_t = %.2f)", rho_t)
  } else {
    "Precision-Weighted Wald Z"
  }

  k_tr    <- length(transitions)
  n_cols  <- min(4L, k_tr)
  n_blk   <- ceiling(k_tr / n_cols)
  all_dz  <- unlist(lapply(transitions, function(x) as.numeric(x$df$delta_sc)))
  lim     <- max(1e-4, as.numeric(quantile(abs(all_dz), probs = 0.985, na.rm = TRUE)))
  ticks_a <- seq(-lim, lim, length.out = 5L)
  labs_a  <- sprintf("%+.2f", ticks_a)

  w_px <- max(2200L, 900L * n_cols)
  h_px <- 1950L * n_blk
  open_png_device(out_file, width = w_px, height = h_px, res = 200)
  on.exit(dev.off(), add = TRUE)
  oma_b    <- get_banner_oma(banner)
  oma_b[3] <- max(oma_b[3], 3.8)
  par(mfrow = c(2L * n_blk, n_cols), mar = c(1.0, 1.0, 3.6, 1.0), oma = oma_b, bg = "white")

  panel_lbl <- function(idx) {
    if (idx <= 26L) LETTERS[idx] else paste0("P", idx)
  }

  for (blk in seq_len(n_blk)) {
    i_start <- (blk - 1L) * n_cols + 1L
    i_end   <- min(k_tr, blk * n_cols)
    idx_vec <- i_start:i_end

    # Row A of block: Model-scale step change (Delta z)
    for (ti in idx_vec) {
      item <- transitions[[ti]]
      plot_stream_network(
        item$df, "delta_sc", vmin = -lim, vmax = lim, pal = TREND_DIV_PAL,
        main_title = sprintf("%s. %s (%s)", panel_lbl(ti), item$transition_label, det_short),
        sub_title  = sprintf("Model-scale change (%s)", tr$short_label),
        cbar_ticks = ticks_a, cbar_labels = labs_a,
        cbar_title = sprintf("Change\n(%s sc)", tr$short_label),
        show_obs = FALSE, base_lwd = 0.20, max_lwd = 1.6,
        max_main_chars = 46, max_sub_chars = 58
      )
    }
    if (length(idx_vec) < n_cols) {
      for (pad in seq_len(n_cols - length(idx_vec))) plot.new()
    }

    # Row B of block: 5-tier Wald Z-test significance classification
    for (ti in idx_vec) {
      item <- transitions[[ti]]
      d  <- item$df[order(item$df$Shreve), ]
      xr <- range(c(d$start_X, d$end_X), na.rm = TRUE)
      yr <- range(c(d$start_Y, d$end_Y), na.rm = TRUE)
      dx <- diff(xr); dy <- diff(yr)
      plot(NA, xlim = c(xr[1] - 0.03 * dx, xr[2] + 0.04 * dx),
           ylim = c(yr[1] - 0.03 * dy, yr[2] + 0.06 * dy),
           asp = 1, axes = FALSE, xlab = "", ylab = "")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#f8f9fa", border = "#d0d4d9", lwd = 1)
      lwds <- 0.20 + (1.6 - 0.20) * pmin(1, log10(pmax(1, d$Shreve)) / log10(max(10, max(d$Shreve))))
      tier_chr <- as.character(d$trend_tier)
      cols <- unname(TREND_TIER_COLORS[tier_chr])
      cols[is.na(cols)] <- "#cfd4da"
      segments(d$start_X, d$start_Y, d$end_X, d$end_Y, col = cols, lwd = lwds)

      tab <- table(factor(tier_chr, levels = TREND_TIER_LEVELS))
      pct <- 100 * as.numeric(tab) / max(1L, nrow(d))
      leg_labels <- sprintf("%s (%.1f%%)", TREND_TIER_LEVELS, pct)
      legend("topleft", legend = leg_labels, col = unname(TREND_TIER_COLORS[TREND_TIER_LEVELS]),
             lwd = 3.0, bg = rgb(1, 1, 1, 0.90), box.col = "#cccccc", cex = 0.68,
             title = "Wald Z Significance", title.cex = 0.72)
      draw_plot_header(
        main_title = sprintf("%s. %s — Significance", panel_lbl(k_tr + ti), item$transition_label),
        sub_title  = "Two-sided Wald Z-test (SE-weighted)",
        max_main = 46, max_sub = 58, cex_main = 0.90, cex_sub = 0.76, adj = 0
      )
    }
    if (length(idx_vec) < n_cols) {
      for (pad in seq_len(n_cols - length(idx_vec))) plot.new()
    }
  }

  mtext(wrap_title(sprintf("All England %s-to-%s Transitions: %s (%s) — %s",
                           mode_cap, mode_cap, det_short, det_unit, st_tag), 180),
        side = 3, line = 1.6, outer = TRUE, cex = 1.15, font = 2, col = "#111111")
  draw_top_banner(banner)
}

#' Plot Multivariate Stream Co-Kriging Cross-Correlation & Dual-CV Diagnostics
#'
#' Renders a two-tier diagnostic PNG: top row shows decomposed Tail-Up, Tail-Down,
#' and Co-Located Nugget cross-correlation heatmaps across modeled targets; bottom
#' panel compares Univariate S3N against Co-Kriging (All-Holdout) and Co-Kriging
#' (Aux-Aided) out-of-sample CV R2 and RMSE across every target determinand.
plot_fig_cokrig_diagnostics <- function(cross_cov_df,
                                        cv_compare_df,
                                        var_names,
                                        out_file,
                                        network_label = NULL,
                                        vars_label = NULL,
                                        source_label = NULL,
                                        method_label = "Factor-LMC Joint NNGP") {
  banner <- format_top_banner(network_label, vars_label, source_label)

  P <- length(var_names)
  short_lbl <- vapply(var_names, function(v) {
    v <- trimws(as.character(v))
    if (nchar(v) > 16L) paste0(substr(v, 1, 14), "..") else v
  }, character(1))

  build_corr_mat <- function(col_nm) {
    M <- diag(1.0, nrow = P, ncol = P)
    for (r in seq_len(nrow(cross_cov_df))) {
      i <- match(cross_cov_df$Var_A[r], var_names)
      j <- match(cross_cov_df$Var_B[r], var_names)
      if (!is.na(i) && !is.na(j)) {
        val <- suppressWarnings(as.numeric(cross_cov_df[[col_nm]][r]))
        if (is.finite(val)) {
          M[i, j] <- max(-1.0, min(1.0, val))
          M[j, i] <- M[i, j]
        }
      }
    }
    M
  }

  div_pal <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
  draw_corr_heatmap <- function(M, title_txt) {
    par(mar = c(5.2, 6.2, 3.0, 1.2))
    # image() expects x along rows; reverse columns vertically so (1,1) is top-left
    Z <- t(M[P:1, , drop = FALSE])
    image(seq_len(P), seq_len(P), Z, zlim = c(-1, 1), col = div_pal,
          axes = FALSE, xlab = "", ylab = "", main = "")
    box(col = "#555555")
    axis(1, at = seq_len(P), labels = short_lbl, las = 2, cex.axis = 0.78, tick = FALSE, line = -0.3)
    axis(2, at = seq_len(P), labels = rev(short_lbl), las = 2, cex.axis = 0.78, tick = FALSE, line = -0.3)
    cex_cell <- if (P <= 4L) 0.92 else if (P <= 6L) 0.76 else 0.62
    for (i in seq_len(P)) {
      for (j in seq_len(P)) {
        v <- M[i, j]
        txt_col <- if (abs(v) > 0.58) "white" else "#111111"
        text(j, P - i + 1L, sprintf("%.2f", v), col = txt_col, cex = cex_cell, font = if (i == j) 2 else 1)
      }
    }
    mtext(title_txt, side = 3, line = 0.8, cex = 0.95, font = 2, col = "#222222")
  }

  png_w <- max(2800L, 1900L + 220L * P)
  open_png_device(out_file, width = png_w, height = 2100L, res = 220)
  on.exit(dev.off(), add = TRUE)
  oma_v <- get_banner_oma(banner)
  oma_v[3] <- oma_v[3] + 1.2
  par(oma = oma_v, bg = "white")
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(1.0, 1.15))

  draw_corr_heatmap(build_corr_mat("Corr_TU"),     "A1. Tail-Up Network Correlation")
  draw_corr_heatmap(build_corr_mat("Corr_TD"),     "A2. Tail-Down Network Correlation")
  draw_corr_heatmap(build_corr_mat("Corr_Nugget"), "A3. Co-Located Nugget Correlation")

  # Panel B: Paired Grouped Bar Chart of Dual CV R2 vs Univariate S3N
  par(mar = c(6.4, 5.2, 3.8, 1.8))
  mod_levels <- c("Univariate S3N", "Co-Kriging (All-Holdout)", "Co-Kriging (Aux-Aided)")
  mod_cols   <- c("#7f8c8d", "#2b8cbe", "#d95f02")
  r2_mat   <- matrix(NA_real_, nrow = 3L, ncol = P)
  rmse_mat <- matrix(NA_real_, nrow = 3L, ncol = P)
  for (mi in seq_along(mod_levels)) {
    for (p in seq_len(P)) {
      idx <- which(cv_compare_df$Target == var_names[p] & cv_compare_df$Model == mod_levels[mi])
      if (length(idx) > 0L) {
        r2_mat[mi, p]   <- suppressWarnings(as.numeric(cv_compare_df$R2[idx[1]]))
        rmse_mat[mi, p] <- suppressWarnings(as.numeric(cv_compare_df$RMSE[idx[1]]))
      }
    }
  }
  r2_min <- min(0.0, min(r2_mat, na.rm = TRUE) - 0.08)
  r2_max <- min(1.12, max(0.35, max(r2_mat, na.rm = TRUE) + 0.16))
  bp <- barplot(
    r2_mat, beside = TRUE, col = mod_cols, border = "#333333",
    ylim = c(r2_min, r2_max), names.arg = rep("", P),
    ylab = expression("Out-of-Sample Cross-Validation " * R^2),
    las = 1, cex.axis = 0.88, cex.lab = 0.95
  )
  abline(h = 0, col = "#444444", lwd = 0.9)
  grid(nx = NA, ny = NULL, col = "#e5e7eb", lty = "dotted")
  # Re-draw bars cleanly over grid lines
  bp <- barplot(
    r2_mat, beside = TRUE, col = mod_cols, border = "#333333",
    ylim = c(r2_min, r2_max), names.arg = rep("", P),
    add = TRUE, axes = FALSE
  )
  for (p in seq_len(P)) {
    for (mi in seq_len(3L)) {
      val <- r2_mat[mi, p]
      if (is.finite(val)) {
        y_txt <- if (val >= 0) val + 0.018 * (r2_max - r2_min) else val - 0.030 * (r2_max - r2_min)
        text(bp[mi, p], y_txt, sprintf("%.3f", val), cex = 0.68, font = 2, col = "#222222")
      }
    }
    r_uni <- rmse_mat[1L, p]
    r_aux <- rmse_mat[3L, p]
    lbl_sub <- var_names[p]
    if (is.finite(r_uni) && is.finite(r_aux) && r_uni > 0) {
      d_pct <- 100 * (r_uni - r_aux) / r_uni
      lbl_sub <- sprintf("%s\n(Aux RMSE %+.1f%%)", var_names[p], -d_pct)
    }
    mtext(lbl_sub, side = 1, at = mean(bp[, p]), line = 1.0, cex = if (P <= 5L) 0.78 else 0.68)
  }
  legend("topleft", legend = mod_levels, fill = mod_cols, border = "#333333",
         horiz = TRUE, bty = "n", cex = 0.85, inset = c(0.01, -0.10), xpd = TRUE)
  mtext("B. Dual Cross-Validation: Univariate S3N vs All-Holdout vs Aux-Aided Co-Kriging",
        side = 3, line = 1.8, adj = 0, cex = 0.98, font = 2, col = "#222222")

  dom_lbl <- if (nrow(cv_compare_df) > 0L) cv_compare_df$Domain[1] else "Domain"
  mtext(wrap_title(sprintf("Multivariate Stream Co-Kriging Diagnostics (%s) — %s", dom_lbl, method_label), 140),
        side = 3, line = 1.3, outer = TRUE, cex = 1.12, font = 2, col = "#111111")
  draw_top_banner(banner)
}
