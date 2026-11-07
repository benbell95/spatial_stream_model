#!/usr/bin/env Rscript
# ==============================================================================
# 06_station_snapping_diagnostics.R
#
# Standalone script and reusable R function (`analyze_station_snapping()`) that
# provides detailed diagnostics, coordinate tables, and catchment-level base-R
# visualizations of monitoring station snapping to river network reaches.
#
# Key Outputs (saved to `output/station_snapping_<network>/` by default):
#   1. `station_snapping_all.csv`:
#        Full station-to-reach snapping table containing:
#          - Station ID (`site_id`, `wims_id`, `wims_region_code`)
#          - Snapped Reach ID (`snapped_RivID`) & Catchment (`CatchID`, `Catchment_Name`)
#          - Station Coordinates (`station_X`, `station_Y` in BNG EPSG:27700 m)
#          - Snapped Reach Coordinates:
#              * Orthogonal projected snap point (`snap_reach_X`, `snap_reach_Y`)
#              * Reach midpoint (`reach_mid_X`, `reach_mid_Y`)
#              * Reach upstream/downstream nodes (`reach_start_X/Y`, `reach_end_X/Y`)
#          - Snapping Distance in metres (`snap_dist_m`) & inclusion status (`included_250m`)
#   2. `catchment_snapping_summary.csv`:
#        Catchment-by-catchment summary of station counts and snapping distance
#        statistics (`Min`, `Q25`, `Median`, `Mean`, `Q75`, `P95`, `Max` in metres).
#   3. `station_snapping_<catchment>.csv`:
#        Per-catchment station snapping table for each analyzed catchment.
#   4. `station_snapping_report.md`:
#        Markdown summary report with network-wide and catchment-level tables.
#   5. `fig_snap_overview_<network>.png`:
#        4-panel network-wide station snapping diagnostic figure.
#   6. `fig_snap_catchment_<CatchID>_<Name>.png`:
#        Catchment-level station & snapped reach map (Panel A: full catchment map
#        highlighting snapped reaches and stations; Panel B: catchment snapping
#        distance histogram; Panel C: 4 metre-scale zoomed insets showing the
#        exact station location, snapped reach geometry, and snapping vector).
#
# Usage (Command Line):
#   Rscript scripts/06_station_snapping_diagnostics.R --list-catchments
#   Rscript scripts/06_station_snapping_diagnostics.R
#   Rscript scripts/06_station_snapping_diagnostics.R --catchment=Thames
#   Rscript scripts/06_station_snapping_diagnostics.R --catchment=Humber,Thames,Severn,Eden
#   Rscript scripts/06_station_snapping_diagnostics.R --network=euhydro_thames --catchment=306
#
# Usage (From R):
#   source("scripts/06_station_snapping_diagnostics.R")
#   res <- analyze_station_snapping(catchment = "Thames", network = "orn")
# ==============================================================================

BASE_DIR <- if (file.exists("scripts/s3n_plot_functions.R")) "." else if (file.exists("s3n_plot_functions.R")) ".." else "."
source(file.path(BASE_DIR, "scripts/s3n_plot_functions.R"))
BASE_DIR <- resolve_base_dir()

# ------------------------------------------------------------------------------
# Helper: Compute exact snapped point (snap_reach_X, snap_reach_Y) on the river
# reach at distance `snap_dist_m` from `(st_x, st_y)` in the orthogonal direction
# of the reach polyline `(start_X, start_Y) -> (mid_X, mid_Y) -> (end_X, end_Y)`.
# ------------------------------------------------------------------------------
compute_snapped_reach_point <- function(st_x, st_y, sx, sy, mx, my, ex, ey, snap_dist_m, along_m = NA_real_, len_m = NA_real_) {
  n <- length(st_x)
  out_x <- numeric(n)
  out_y <- numeric(n)

  proj_seg <- function(px, py, x1, y1, x2, y2) {
    dx <- x2 - x1
    dy <- y2 - y1
    l2 <- dx * dx + dy * dy
    if (l2 <= 1e-12) return(c(x1, y1, (px - x1)^2 + (py - y1)^2))
    t_val <- max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / l2))
    qx <- x1 + t_val * dx
    qy <- y1 + t_val * dy
    c(qx, qy, (px - qx)^2 + (py - qy)^2)
  }

  for (i in seq_len(n)) {
    if (is.na(sx[i]) || is.na(sy[i]) || is.na(mx[i]) || is.na(my[i]) || is.na(ex[i]) || is.na(ey[i])) {
      out_x[i] <- NA_real_
      out_y[i] <- NA_real_
      next
    }
    p1 <- proj_seg(st_x[i], st_y[i], sx[i], sy[i], mx[i], my[i])
    p2 <- proj_seg(st_x[i], st_y[i], mx[i], my[i], ex[i], ey[i])
    qx <- if (p1[3] <= p2[3]) p1[1] else p2[1]
    qy <- if (p1[3] <= p2[3]) p1[2] else p2[2]

    vx <- qx - st_x[i]
    vy <- qy - st_y[i]
    v_norm <- sqrt(vx * vx + vy * vy)
    d_true <- ifelse(is.na(snap_dist_m[i]), v_norm, snap_dist_m[i])

    if (v_norm > 1e-6) {
      out_x[i] <- round(st_x[i] + d_true * (vx / v_norm), 2)
      out_y[i] <- round(st_y[i] + d_true * (vy / v_norm), 2)
    } else {
      # Station sits directly on the chord; offset along normal to reach direction
      rx <- ex[i] - sx[i]
      ry <- ey[i] - sy[i]
      rn <- sqrt(rx * rx + ry * ry)
      if (rn > 1e-6) {
        out_x[i] <- round(st_x[i] + d_true * (-ry / rn), 2)
        out_y[i] <- round(st_y[i] + d_true * (rx / rn), 2)
      } else {
        out_x[i] <- round(st_x[i] + d_true, 2)
        out_y[i] <- round(st_y[i], 2)
      }
    }
  }
  list(x = out_x, y = out_y)
}

# ------------------------------------------------------------------------------
# Helper: Load and merge `rsn_sites.csv`, `england_rsn_reaches.csv`, and
# `rsn_catchments_lookup.csv` into a rich station snapping data frame.
# ------------------------------------------------------------------------------
build_station_snapping_table <- function(proc_dir, include_excluded = TRUE) {
  sites_csv   <- file.path(proc_dir, "rsn_sites.csv")
  reaches_csv <- file.path(proc_dir, "england_rsn_reaches.csv")
  catch_csv   <- file.path(proc_dir, "rsn_catchments_lookup.csv")

  if (!file.exists(sites_csv)) {
    stop(sprintf("Station snapping file not found: %s. Run scripts/01_preprocess_orn.py first.", sites_csv))
  }
  if (!file.exists(reaches_csv)) {
    stop(sprintf("Processed reaches file not found: %s. Run scripts/01_preprocess_orn.py first.", reaches_csv))
  }

  sites_df   <- read.csv(sites_csv, stringsAsFactors = FALSE)
  reaches_df <- read.csv(reaches_csv, stringsAsFactors = FALSE)
  catch_lu   <- if (file.exists(catch_csv)) read.csv(catch_csv, stringsAsFactors = FALSE) else NULL

  # Clean numeric fields in sites_df
  sites_df$easting        <- suppressWarnings(as.numeric(sites_df$easting))
  sites_df$northing       <- suppressWarnings(as.numeric(sites_df$northing))
  sites_df$snap_dist_m    <- suppressWarnings(as.numeric(sites_df$snap_dist_m))
  sites_df$RivID          <- suppressWarnings(as.integer(sites_df$RivID))
  sites_df$CatchID        <- suppressWarnings(as.integer(sites_df$CatchID))
  sites_df$obs_up_dist_km <- suppressWarnings(as.numeric(sites_df$obs_up_dist_km))
  sites_df$n_samples      <- suppressWarnings(as.integer(sites_df$n_samples))
  sites_df$included_250m  <- suppressWarnings(as.integer(sites_df$included_250m))
  sites_df <- sites_df[!is.na(sites_df$easting) & !is.na(sites_df$northing) & !is.na(sites_df$RivID), ]

  if (!include_excluded) {
    sites_df <- sites_df[sites_df$included_250m == 1L, ]
  }

  # Match each station's snapped RivID to `reaches_df`
  r_idx <- match(sites_df$RivID, reaches_df$RivID)

  reach_start_X  <- reaches_df$start_X[r_idx]
  reach_start_Y  <- reaches_df$start_Y[r_idx]
  reach_mid_X    <- reaches_df$mid_X[r_idx]
  reach_mid_Y    <- reaches_df$mid_Y[r_idx]
  reach_end_X    <- reaches_df$end_X[r_idx]
  reach_end_Y    <- reaches_df$end_Y[r_idx]
  reach_len_m    <- round(reaches_df$length_km[r_idx] * 1000.0, 1)
  d2m_km         <- reaches_df$Dist2Mth_km[r_idx]
  along_m        <- round(pmax(0.0, (d2m_km - sites_df$obs_up_dist_km) * 1000.0), 1)
  strahler       <- reaches_df$Strahler[r_idx]
  shreve         <- reaches_df$Shreve[r_idx]
  elevation_m    <- reaches_df$Elevation[r_idx]

  # Compute orthogonal snapped reach point (snap_reach_X, snap_reach_Y)
  snap_pts <- compute_snapped_reach_point(
    sites_df$easting, sites_df$northing,
    reach_start_X, reach_start_Y,
    reach_mid_X, reach_mid_Y,
    reach_end_X, reach_end_Y,
    snap_dist_m = sites_df$snap_dist_m,
    along_m = along_m, len_m = reach_len_m
  )

  # Attach Catchment Name and RBD Name
  catch_name <- rep(NA_character_, nrow(sites_df))
  rbd_name   <- rep(NA_character_, nrow(sites_df))
  if (!is.null(catch_lu) && nrow(catch_lu) > 0) {
    c_idx <- match(sites_df$CatchID, catch_lu$CatchID)
    catch_name <- ifelse(is.na(c_idx), paste0("Catchment ", sites_df$CatchID), catch_lu$Main_River[c_idx])
    rbd_name   <- ifelse(is.na(c_idx), "", catch_lu$RBD_Name[c_idx])
  } else {
    catch_name <- paste0("Catchment ", sites_df$CatchID)
  }

  snap_df <- data.frame(
    site_id          = sites_df$site_id,
    source           = if ("source" %in% colnames(sites_df)) sites_df$source else "RSN",
    wims_id          = sites_df$wims_id,
    wims_region_code = sites_df$wims_region_code,
    CatchID          = sites_df$CatchID,
    Catchment_Name   = catch_name,
    RBD_Name         = rbd_name,
    station_X        = round(sites_df$easting, 1),
    station_Y        = round(sites_df$northing, 1),
    snapped_RivID    = sites_df$RivID,
    snap_reach_X     = snap_pts$x,
    snap_reach_Y     = snap_pts$y,
    reach_mid_X      = round(reach_mid_X, 1),
    reach_mid_Y      = round(reach_mid_Y, 1),
    reach_start_X    = round(reach_start_X, 1),
    reach_start_Y    = round(reach_start_Y, 1),
    reach_end_X      = round(reach_end_X, 1),
    reach_end_Y      = round(reach_end_Y, 1),
    snap_dist_m      = round(sites_df$snap_dist_m, 3),
    included_250m    = sites_df$included_250m,
    snap_status      = ifelse(sites_df$included_250m == 1L, "Included (<=250m)", "Excluded (>250m)"),
    reach_length_m   = reach_len_m,
    along_reach_m    = along_m,
    obs_up_dist_km   = round(sites_df$obs_up_dist_km, 4),
    Strahler         = strahler,
    Shreve           = shreve,
    Elevation_m      = round(elevation_m, 2),
    n_samples        = sites_df$n_samples,
    stringsAsFactors = FALSE
  )

  # Sort by CatchID and snap_dist_m descending (or site_id)
  snap_df <- snap_df[order(snap_df$CatchID, -snap_df$included_250m, snap_df$site_id), ]
  rownames(snap_df) <- NULL

  list(
    snap_df    = snap_df,
    reaches_df = reaches_df,
    catch_lu   = catch_lu
  )
}

# ------------------------------------------------------------------------------
# Helper: Compute catchment-level snapping summary table
# ------------------------------------------------------------------------------
compute_catchment_snapping_summary <- function(snap_df, reaches_df, catch_lu = NULL) {
  cids <- sort(unique(c(reaches_df$CatchID, snap_df$CatchID[!is.na(snap_df$reach_mid_X)])))
  rows <- vector("list", length(cids))

  for (i in seq_along(cids)) {
    cid <- cids[i]
    sub_r <- reaches_df[reaches_df$CatchID == cid, ]
    sub_s <- snap_df[snap_df$CatchID == cid, ]
    sub_inc <- sub_s[sub_s$included_250m == 1L, ]
    sub_exc <- sub_s[sub_s$included_250m == 0L, ]

    c_name <- paste0("Catchment ", cid)
    r_name <- ""
    if (!is.null(catch_lu) && cid %in% catch_lu$CatchID) {
      idx <- match(cid, catch_lu$CatchID)
      c_name <- catch_lu$Main_River[idx]
      r_name <- catch_lu$RBD_Name[idx]
    } else if (nrow(sub_s) > 0 && nzchar(sub_s$Catchment_Name[1])) {
      c_name <- sub_s$Catchment_Name[1]
      r_name <- sub_s$RBD_Name[1]
    }

    d_inc <- sub_inc$snap_dist_m[!is.na(sub_inc$snap_dist_m)]
    d_all <- sub_s$snap_dist_m[!is.na(sub_s$snap_dist_m)]
    d_use <- if (length(d_inc) > 0) d_inc else d_all

    if (length(d_use) > 0) {
      qs <- quantile(d_use, probs = c(0.25, 0.50, 0.75, 0.95), na.rm = TRUE, names = FALSE)
      min_d  <- round(min(d_use, na.rm = TRUE), 2)
      q25_d  <- round(qs[1], 2)
      med_d  <- round(qs[2], 2)
      mean_d <- round(mean(d_use, na.rm = TRUE), 2)
      q75_d  <- round(qs[3], 2)
      p95_d  <- round(qs[4], 2)
      max_d  <- round(max(d_use, na.rm = TRUE), 2)
      pct_10m <- round(100.0 * mean(d_use <= 10.0), 1)
      pct_50m <- round(100.0 * mean(d_use <= 50.0), 1)
    } else {
      min_d <- q25_d <- med_d <- mean_d <- q75_d <- p95_d <- max_d <- pct_10m <- pct_50m <- NA_real_
    }

    rows[[i]] <- data.frame(
      CatchID            = cid,
      Main_River         = c_name,
      RBD_Name           = r_name,
      Total_Reaches      = nrow(sub_r),
      Snapped_Reaches    = length(unique(sub_inc$snapped_RivID)),
      Stations_Included  = nrow(sub_inc),
      Stations_Excluded  = nrow(sub_exc),
      Stations_Total     = nrow(sub_s),
      Min_Snap_m         = min_d,
      Q25_Snap_m         = q25_d,
      Median_Snap_m      = med_d,
      Mean_Snap_m        = mean_d,
      Q75_Snap_m         = q75_d,
      P95_Snap_m         = p95_d,
      Max_Snap_m         = max_d,
      Pct_Within_10m     = pct_10m,
      Pct_Within_50m     = pct_50m,
      stringsAsFactors   = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out <- out[order(-out$Stations_Included, -out$Total_Reaches, out$CatchID), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# Helper: Resolve user-supplied catchment specifications (IDs or river names)
# ------------------------------------------------------------------------------
resolve_snapping_catchments <- function(specs, catch_sum) {
  if (length(specs) == 0 || is.null(specs)) return(integer(0))
  if (length(specs) == 1 && tolower(trimws(specs[1])) == "all") {
    return(catch_sum$CatchID[catch_sum$Stations_Total > 0])
  }
  resolved <- integer(0)
  for (sp in specs) {
    sp_clean <- trimws(sp)
    if (!nzchar(sp_clean)) next
    if (grepl("^[0-9]+$", sp_clean)) {
      cid <- as.integer(sp_clean)
      if (cid %in% catch_sum$CatchID) {
        resolved <- c(resolved, cid)
      } else {
        warning(sprintf("CatchID %d not found in network.", cid))
      }
    } else {
      m_idx <- which(
        grepl(sp_clean, catch_sum$Main_River, ignore.case = TRUE) |
        grepl(sp_clean, catch_sum$RBD_Name, ignore.case = TRUE)
      )
      if (length(m_idx) > 0) {
        # Take the top matching catchment with the most stations
        best_cid <- catch_sum$CatchID[m_idx[1]]
        resolved <- c(resolved, best_cid)
      } else {
        warning(sprintf("No catchment matching '%s' found. Run with --list-catchments to see valid names.", sp_clean))
      }
    }
  }
  unique(resolved)
}

# ------------------------------------------------------------------------------
# Color classification helper for snapping distances (in metres)
# ------------------------------------------------------------------------------
get_snap_dist_colors <- function(dist_m, inc_flag) {
  cols <- character(length(dist_m))
  for (i in seq_along(dist_m)) {
    if (is.na(inc_flag[i]) || inc_flag[i] == 0L || dist_m[i] > 250.0) {
      cols[i] <- "#B2182B" # Crimson Red (> 250 m excluded)
    } else if (dist_m[i] <= 5.0) {
      cols[i] <- "#1B7837" # Dark Green (<= 5 m)
    } else if (dist_m[i] <= 25.0) {
      cols[i] <- "#2166AC" # Royal Blue (5 - 25 m)
    } else if (dist_m[i] <= 100.0) {
      cols[i] <- "#E68600" # Amber Gold (25 - 100 m)
    } else {
      cols[i] <- "#D95F02" # Dark Orange (100 - 250 m)
    }
  }
  cols
}

# ------------------------------------------------------------------------------
# Plot 1: Catchment-Level Station & Snapped Reach Diagnostic Plot
#   - Panel A (Left): Full catchment river network, highlighted snapped reaches,
#                     station locations colored by snapping distance, and vectors.
#   - Panel B (Top-Right): Histogram & quantiles of station snapping distances (m).
#   - Panel C (Bottom-Right 2x2): Metre-scale zoomed insets of 4 stations in the
#                     catchment showing station point, snapped reach polyline,
#                     reach midpoint, orthogonal snap point, and metre offset.
# ------------------------------------------------------------------------------
plot_catchment_snapping <- function(cid, snap_df, reaches_df, catch_sum, net_label, out_png) {
  sub_r <- reaches_df[reaches_df$CatchID == cid, ]
  sub_s <- snap_df[snap_df$CatchID == cid & !is.na(snap_df$reach_mid_X), ]
  if (nrow(sub_r) == 0 || nrow(sub_s) == 0) return(invisible(NULL))

  cs_row <- catch_sum[catch_sum$CatchID == cid, ]
  c_name <- if (nrow(cs_row) > 0) cs_row$Main_River[1] else paste0("Catchment ", cid)
  r_name <- if (nrow(cs_row) > 0 && nzchar(cs_row$RBD_Name[1])) paste0(" (", cs_row$RBD_Name[1], " RBD)") else ""

  png(out_png, width = 2500, height = 1750, res = 150, type = "cairo")
  on.exit(dev.off(), add = TRUE)

  # Layout: Left half = Panel A (1, 3, 5), Top-right = Panel B (2), Bottom-right 2x2 = Insets (4, 6, 7, 8)
  lay_mat <- matrix(c(
    1, 1, 2, 2,
    1, 1, 3, 4,
    1, 1, 5, 6
  ), nrow = 3, byrow = TRUE)
  layout(lay_mat, widths = c(1.3, 1.3, 1.0, 1.0), heights = c(1.05, 1.0, 1.0))

  # ============================================================================
  # PANEL A: Full Catchment River Network, Snapped Reaches & Station Locations
  # ============================================================================
  par(mar = c(4.2, 4.2, 4.2, 1.5), bg = "#FFFFFF")
  all_x <- c(sub_r$start_X, sub_r$end_X, sub_s$station_X)
  all_y <- c(sub_r$start_Y, sub_r$end_Y, sub_s$station_Y)
  xr <- range(all_x, na.rm = TRUE)
  yr <- range(all_y, na.rm = TRUE)
  pad_x <- max(500, diff(xr) * 0.06)
  pad_y <- max(500, diff(yr) * 0.06)

  plot(NA, xlim = c(xr[1] - pad_x, xr[2] + pad_x), ylim = c(yr[1] - pad_y, yr[2] + pad_y),
       asp = 1, xlab = "Easting (m, BNG EPSG:27700)", ylab = "Northing (m, BNG EPSG:27700)",
       main = "", cex.lab = 1.1, cex.axis = 0.95)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F8FAFC", border = "#CBD5E1")
  grid(col = "#E2E8F0", lty = 1, lwd = 0.7)

  title(main = sprintf("A. %s%s [CatchID = %d] — Station Snapping Map", c_name, r_name, cid),
        cex.main = 1.28, font.main = 2, line = 2.2)
  mtext(sprintf("%s | %s total reaches | %d snapped reaches | %d stations (<=250m: %d, >250m: %d)",
                net_label, format(nrow(sub_r), big.mark = ","),
                length(unique(sub_s$snapped_RivID[sub_s$included_250m == 1L])),
                nrow(sub_s), sum(sub_s$included_250m == 1L), sum(sub_s$included_250m == 0L)),
        side = 3, line = 0.6, cex = 0.88, col = "#334155")

  # 1. Draw all unmonitored reaches in light grey (2-segment polylines: start -> mid -> end)
  unmon_r <- sub_r[!(sub_r$RivID %in% sub_s$snapped_RivID), ]
  if (nrow(unmon_r) > 0) {
    lwd_u <- pmax(0.45, pmin(2.2, 0.35 * unmon_r$Strahler))
    segments(unmon_r$start_X, unmon_r$start_Y, unmon_r$mid_X, unmon_r$mid_Y, col = "#CBD5E1", lwd = lwd_u)
    segments(unmon_r$mid_X, unmon_r$mid_Y, unmon_r$end_X, unmon_r$end_Y, col = "#CBD5E1", lwd = lwd_u)
  }

  # 2. Highlight all SNAPPED reaches in bold vibrant cyan/blue
  snapped_r <- sub_r[sub_r$RivID %in% sub_s$snapped_RivID, ]
  if (nrow(snapped_r) > 0) {
    segments(snapped_r$start_X, snapped_r$start_Y, snapped_r$mid_X, snapped_r$mid_Y, col = "#0284C7", lwd = 3.4)
    segments(snapped_r$mid_X, snapped_r$mid_Y, snapped_r$end_X, snapped_r$end_Y, col = "#0284C7", lwd = 3.4)
  }

  # 3. Draw connecting lines between station (station_X, station_Y) and snapped reach (snap_reach_X, snap_reach_Y)
  segments(sub_s$station_X, sub_s$station_Y, sub_s$snap_reach_X, sub_s$snap_reach_Y,
           col = "#DC2626", lwd = 1.6, lty = 2)

  # 4. Plot snapped reach points (blue squares) and station points (colored circles / red crosses)
  points(sub_s$snap_reach_X, sub_s$snap_reach_Y, pch = 15, col = "#0369A1", cex = 0.95)

  pt_cols <- get_snap_dist_colors(sub_s$snap_dist_m, sub_s$included_250m)
  inc_k   <- which(sub_s$included_250m == 1L)
  exc_k   <- which(sub_s$included_250m == 0L)

  if (length(inc_k) > 0) {
    points(sub_s$station_X[inc_k], sub_s$station_Y[inc_k],
           pch = 21, bg = pt_cols[inc_k], col = "#0F172A", cex = 1.45, lwd = 1.1)
  }
  if (length(exc_k) > 0) {
    points(sub_s$station_X[exc_k], sub_s$station_Y[exc_k],
           pch = 4, col = "#B2182B", cex = 1.65, lwd = 2.2)
  }

  # Label the top 6 stations with the largest snapping distance in the catchment
  ord_lbl <- order(-sub_s$snap_dist_m)
  top_lbl <- head(ord_lbl, min(6L, nrow(sub_s)))
  text(sub_s$station_X[top_lbl], sub_s$station_Y[top_lbl],
       labels = sprintf("%s (%.1fm)", sub_s$site_id[top_lbl], sub_s$snap_dist_m[top_lbl]),
       pos = 3, offset = 0.45, cex = 0.72, font = 2, col = "#1E293B")

  # Legend on Panel A
  legend("topleft",
         legend = c(
           "Unmonitored Reach (Strahler scaled)",
           "Snapped River Reach (snapped_RivID)",
           "Snapped Reach Point (snap_reach_X/Y)",
           "Station <= 5 m from reach",
           "Station 5 - 25 m from reach",
           "Station 25 - 100 m from reach",
           "Station 100 - 250 m from reach",
           "Excluded Station > 250 m"
         ),
         col    = c("#94A3B8", "#0284C7", "#0369A1", "#0F172A", "#0F172A", "#0F172A", "#0F172A", "#B2182B"),
         pt.bg  = c(NA, NA, NA, "#1B7837", "#2166AC", "#E68600", "#D95F02", NA),
         pch    = c(NA, NA, 15, 21, 21, 21, 21, 4),
         lty    = c(1, 1, NA, NA, NA, NA, NA, NA),
         lwd    = c(1.5, 3.2, NA, 1.1, 1.1, 1.1, 1.1, 2.0),
         pt.cex = c(NA, NA, 1.0, 1.3, 1.3, 1.3, 1.3, 1.4),
         bg = "#FFFFFFEE", box.col = "#CBD5E1", cex = 0.80, title = "Map Symbol Key")

  # ============================================================================
  # PANEL B: Snapping Distance Distribution (in Metres) for the Catchment
  # ============================================================================
  par(mar = c(4.2, 4.2, 3.8, 1.5))
  d_inc <- sub_s$snap_dist_m[sub_s$included_250m == 1L]
  d_plot <- if (length(d_inc) > 0) d_inc else sub_s$snap_dist_m
  brks <- seq(0, max(25, ceiling(max(d_plot, na.rm = TRUE) * 1.05)), length.out = 22)

  hist(d_plot, breaks = brks, col = "#38BDF8", border = "#0369A1",
       main = sprintf("B. Snapping Distance Distribution (%s, n = %d)", c_name, length(d_plot)),
       xlab = "Station-to-Reach Snapping Distance (metres)", ylab = "Number of Stations",
       cex.main = 1.12, cex.lab = 1.0)
  grid(col = "#E2E8F0", lty = 1)
  rug(d_plot, col = "#0F172A", lwd = 1.2)

  med_d <- median(d_plot, na.rm = TRUE)
  p95_d <- quantile(d_plot, probs = 0.95, na.rm = TRUE, names = FALSE)
  mean_d <- mean(d_plot, na.rm = TRUE)
  abline(v = med_d, col = "#15803D", lwd = 2.4, lty = 2)
  abline(v = p95_d, col = "#D97706", lwd = 2.2, lty = 3)

  legend("topright",
         legend = c(
           sprintf("Median Snap Distance: %.2f m", med_d),
           sprintf("Mean Snap Distance:   %.2f m", mean_d),
           sprintf("95th Percentile (P95): %.2f m", p95_d),
           sprintf("Min / Max Distance:   %.2f / %.2f m", min(d_plot, na.rm = TRUE), max(d_plot, na.rm = TRUE)),
           sprintf("Within <= 10 m:       %.1f%% of stations", 100 * mean(d_plot <= 10.0)),
           sprintf("Within <= 50 m:       %.1f%% of stations", 100 * mean(d_plot <= 50.0))
         ),
         col = c("#15803D", NA, "#D97706", NA, NA, NA),
         lty = c(2, NA, 3, NA, NA, NA),
         lwd = c(2.2, NA, 2.2, NA, NA, NA),
         bg = "#FFFFFFEE", box.col = "#CBD5E1", cex = 0.84, title = "Catchment Snapping Summary")

  # ============================================================================
  # PANEL C (4 Zoomed Sub-Panels): Metre-Scale Local View of Station -> Reach Snap
  # ============================================================================
  # Select up to 4 informative stations in this catchment:
  # Top 3 largest-distance included stations + median station (or largest overall)
  s_inc_idx <- which(sub_s$included_250m == 1L)
  if (length(s_inc_idx) >= 4) {
    ord_inc <- s_inc_idx[order(-sub_s$snap_dist_m[s_inc_idx])]
    med_idx <- ord_inc[max(1L, round(length(ord_inc) * 0.5))]
    sel_idx <- unique(c(ord_inc[1:3], med_idx))
    if (length(sel_idx) < 4) sel_idx <- head(ord_inc, 4)
  } else {
    sel_idx <- head(order(-sub_s$snap_dist_m), min(4L, nrow(sub_s)))
  }

  for (k in 1:4) {
    par(mar = c(3.6, 3.8, 2.9, 1.0), mgp = c(2.0, 0.6, 0))
    if (k > length(sel_idx)) {
      plot.new()
      next
    }
    st <- sub_s[sel_idx[k], ]
    # Focus zoom box around the station and its snapped reach point, scaled appropriately to snap_dist_m
    span <- max(55.0, st$snap_dist_m * 2.15)
    cx <- 0.5 * (st$station_X + st$snap_reach_X)
    cy <- 0.5 * (st$station_Y + st$snap_reach_Y)

    plot(NA, xlim = c(cx - span, cx + span), ylim = c(cy - span, cy + span),
         asp = 1, xlab = "Easting (m)", ylab = "Northing (m)",
         main = sprintf("C%d. Zoom: %s -> Reach %d (d = %.2f m)",
                        k, st$site_id, st$snapped_RivID, st$snap_dist_m),
         cex.main = 0.94, cex.lab = 0.85, cex.axis = 0.78)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F8FAFC", border = "#CBD5E1")
    grid(col = "#E2E8F0", lty = 1)

    # Draw nearby reaches within window
    near_r <- sub_r[
      pmax(sub_r$start_X, sub_r$end_X) >= (cx - span * 2.0) &
      pmin(sub_r$start_X, sub_r$end_X) <= (cx + span * 2.0) &
      pmax(sub_r$start_Y, sub_r$end_Y) >= (cy - span * 2.0) &
      pmin(sub_r$start_Y, sub_r$end_Y) <= (cy + span * 2.0),
    ]
    if (nrow(near_r) > 0) {
      other_r <- near_r[near_r$RivID != st$snapped_RivID, ]
      if (nrow(other_r) > 0) {
        segments(other_r$start_X, other_r$start_Y, other_r$mid_X, other_r$mid_Y, col = "#94A3B8", lwd = 2.0)
        segments(other_r$mid_X, other_r$mid_Y, other_r$end_X, other_r$end_Y, col = "#94A3B8", lwd = 2.0)
      }
    }

    # Draw the target snapped reach passing through start -> snap_reach -> mid -> end (ordered along reach)
    if (!is.na(st$along_reach_m) && !is.na(st$reach_length_m) && st$along_reach_m > 0.5 * st$reach_length_m) {
      rx_pts <- c(st$reach_start_X, st$reach_mid_X, st$snap_reach_X, st$reach_end_X)
      ry_pts <- c(st$reach_start_Y, st$reach_mid_Y, st$snap_reach_Y, st$reach_end_Y)
    } else {
      rx_pts <- c(st$reach_start_X, st$snap_reach_X, st$reach_mid_X, st$reach_end_X)
      ry_pts <- c(st$reach_start_Y, st$snap_reach_Y, st$reach_mid_Y, st$reach_end_Y)
    }
    lines(rx_pts, ry_pts, col = "#0284C7", lwd = 4.2)

    # Mark reach midpoint (hollow diamond) and orthogonal snapped point (filled blue square)
    points(st$reach_mid_X, st$reach_mid_Y, pch = 5, col = "#1D4ED8", cex = 1.4, lwd = 1.8)
    points(st$snap_reach_X, st$snap_reach_Y, pch = 15, col = "#0369A1", cex = 1.55)

    # Draw the red dashed snapping line from Station to Snapped Reach Point
    segments(st$station_X, st$station_Y, st$snap_reach_X, st$snap_reach_Y,
             col = "#DC2626", lwd = 2.5, lty = 2)

    # Plot Station Point
    st_col <- get_snap_dist_colors(st$snap_dist_m, st$included_250m)
    points(st$station_X, st$station_Y, pch = 21, bg = st_col, col = "#0F172A", cex = 1.85, lwd = 1.4)

    # Label station and snapped reach point inside the zoom inset
    text(st$station_X, st$station_Y, labels = st$site_id, pos = 3, offset = 0.55, cex = 0.76, font = 2, col = "#0F172A")
    text(st$snap_reach_X, st$snap_reach_Y, labels = sprintf("RivID %d", st$snapped_RivID), pos = 1, offset = 0.55, cex = 0.74, col = "#0369A1", font = 2)

    # Annotate coordinates and distance
    mtext(sprintf("Stn: (%.1f, %.1f) | SnapPt: (%.1f, %.1f)",
                  st$station_X, st$station_Y, st$snap_reach_X, st$snap_reach_Y),
          side = 3, line = 0.12, cex = 0.66, col = "#475569")
    if (k == 1) {
      legend("bottomleft",
             legend = c("Snapped Reach", "Snap Point", "Reach Midpt", "Station", "Snap Vector"),
             col    = c("#0284C7", "#0369A1", "#1D4ED8", "#0F172A", "#DC2626"),
             pt.bg  = c(NA, NA, NA, st_col, NA),
             pch    = c(NA, 15, 5, 21, NA),
             lty    = c(1, NA, NA, NA, 2),
             lwd    = c(3.5, NA, 1.5, 1.2, 2.0),
             bg = "#FFFFFFDD", box.col = "#CBD5E1", cex = 0.66)
    }
  }

  invisible(out_png)
}

# ------------------------------------------------------------------------------
# Plot 2: Network-Wide Station Snapping Overview Figure (4 Panels)
# ------------------------------------------------------------------------------
plot_network_snapping_overview <- function(snap_df, catch_sum, net_label, out_png) {
  png(out_png, width = 2400, height = 1750, res = 150, type = "cairo")
  on.exit(dev.off(), add = TRUE)

  par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3.8, 1.8), bg = "#FFFFFF")

  # Panel A: Spatial Map of All Stations Colored by Snapping Distance Tier
  xr <- range(snap_df$station_X, na.rm = TRUE)
  yr <- range(snap_df$station_Y, na.rm = TRUE)
  plot(NA, xlim = xr, ylim = yr, asp = 1,
       xlab = "Easting (m, BNG EPSG:27700)", ylab = "Northing (m, BNG EPSG:27700)",
       main = sprintf("A. All RSN Monitoring Stations by Snapping Distance (%s)", net_label),
       cex.main = 1.15, cex.lab = 1.0)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "#F8FAFC", border = "#CBD5E1")
  grid(col = "#E2E8F0", lty = 1)

  pt_cols <- get_snap_dist_colors(snap_df$snap_dist_m, snap_df$included_250m)
  inc_k <- which(snap_df$included_250m == 1L)
  exc_k <- which(snap_df$included_250m == 0L)
  if (length(inc_k) > 0) {
    points(snap_df$station_X[inc_k], snap_df$station_Y[inc_k],
           pch = 21, bg = pt_cols[inc_k], col = "#0F172A", cex = 0.95, lwd = 0.6)
  }
  if (length(exc_k) > 0) {
    points(snap_df$station_X[exc_k], snap_df$station_Y[exc_k],
           pch = 4, col = "#B2182B", cex = 1.05, lwd = 1.4)
  }
  legend("topleft",
         legend = c(
           sprintf("<= 5 m (n = %d)", sum(snap_df$included_250m == 1L & snap_df$snap_dist_m <= 5)),
           sprintf("5 - 25 m (n = %d)", sum(snap_df$included_250m == 1L & snap_df$snap_dist_m > 5 & snap_df$snap_dist_m <= 25)),
           sprintf("25 - 100 m (n = %d)", sum(snap_df$included_250m == 1L & snap_df$snap_dist_m > 25 & snap_df$snap_dist_m <= 100)),
           sprintf("100 - 250 m (n = %d)", sum(snap_df$included_250m == 1L & snap_df$snap_dist_m > 100 & snap_df$snap_dist_m <= 250)),
           sprintf("Excluded > 250 m (n = %d)", length(exc_k))
         ),
         col   = c("#0F172A", "#0F172A", "#0F172A", "#0F172A", "#B2182B"),
         pt.bg = c("#1B7837", "#2166AC", "#E68600", "#D95F02", NA),
         pch   = c(21, 21, 21, 21, 4),
         pt.cex = 1.2, bg = "#FFFFFFEE", box.col = "#CBD5E1", cex = 0.85, title = "Snapping Distance Tier")

  # Panel B: Empirical CDF & Histogram of Included Station Snapping Distances (0 - 250 m)
  d_inc <- snap_df$snap_dist_m[snap_df$included_250m == 1L]
  hist(d_inc, breaks = 30, col = "#38BDF8", border = "#0284C7",
       main = sprintf("B. Snapping Distance Distribution (Included <= 250 m, n = %d)", length(d_inc)),
       xlab = "Snapping Distance to Nearest River Reach (metres)", ylab = "Station Count",
       cex.main = 1.15, cex.lab = 1.0)
  grid(col = "#E2E8F0", lty = 1)
  med_inc <- median(d_inc, na.rm = TRUE)
  p95_inc <- quantile(d_inc, 0.95, na.rm = TRUE, names = FALSE)
  abline(v = med_inc, col = "#15803D", lwd = 2.4, lty = 2)
  abline(v = p95_inc, col = "#D97706", lwd = 2.4, lty = 3)
  legend("topright",
         legend = c(
           sprintf("Median: %.2f m (IQR: %.2f m)", med_inc, IQR(d_inc, na.rm = TRUE)),
           sprintf("Mean:   %.2f m (SD:  %.2f m)", mean(d_inc, na.rm = TRUE), sd(d_inc, na.rm = TRUE)),
           sprintf("95th Percentile (P95): %.2f m", p95_inc),
           sprintf("<= 10 m: %.1f%% | <= 50 m: %.1f%%", 100 * mean(d_inc <= 10), 100 * mean(d_inc <= 50))
         ),
         col = c("#15803D", NA, "#D97706", NA), lty = c(2, NA, 3, NA), lwd = 2.2,
         bg = "#FFFFFFEE", box.col = "#CBD5E1", cex = 0.86)

  # Panel C: Snapping Distance by Strahler Stream Order (Included Stations)
  sub_inc <- snap_df[snap_df$included_250m == 1L & !is.na(snap_df$Strahler), ]
  if (nrow(sub_inc) > 0) {
    ord_grp <- ifelse(sub_inc$Strahler >= 6, "6+", as.character(sub_inc$Strahler))
    ord_lvls <- intersect(c("1", "2", "3", "4", "5", "6+"), unique(ord_grp))
    ord_fac <- factor(ord_grp, levels = ord_lvls)
    boxplot(snap_dist_m ~ ord_fac, data = sub_inc,
            col = c("#E0F2FE", "#BAE6FD", "#7DD3FC", "#38BDF8", "#0284C7", "#0369A1")[seq_along(ord_lvls)],
            border = "#0F172A", outline = FALSE,
            main = "C. Station Snapping Distance by Strahler Stream Order",
            xlab = "Strahler Stream Order of Snapped Reach",
            ylab = "Snapping Distance (metres, outliers hidden)",
            cex.main = 1.15, cex.lab = 1.0)
    grid(col = "#E2E8F0", lty = 1)
    abline(h = med_inc, col = "#DC2626", lty = 2, lwd = 1.8)
    cnts <- table(ord_fac)
    mtext(sprintf("n = %s", paste(paste0("Ord ", names(cnts), ": ", cnts), collapse = " | ")),
          side = 3, line = 0.3, cex = 0.78, col = "#475569")
  } else {
    plot.new()
  }

  # Panel D: Top 12 Catchments by Number of Snapped Stations (with Median Distance Labels)
  top_c <- head(catch_sum[catch_sum$Stations_Included > 0, ], 12)
  if (nrow(top_c) > 0) {
    par(mar = c(4.5, 9.5, 3.8, 2.5))
    top_c <- top_c[rev(seq_len(nrow(top_c))), ]
    lbls <- sprintf("%s (#%d)", substr(top_c$Main_River, 1, 18), top_c$CatchID)
    bp <- barplot(top_c$Stations_Included, names.arg = lbls, horiz = TRUE, las = 1,
                  col = "#0284C7", border = "#0F172A",
                  xlim = c(0, max(top_c$Stations_Included) * 1.28),
                  main = "D. Top 12 Catchments by Snapped Stations (with Median Snap m)",
                  xlab = "Number of Included Monitoring Stations (<= 250 m)",
                  cex.names = 0.84, cex.main = 1.12, cex.lab = 1.0)
    grid(nx = NULL, ny = NA, col = "#E2E8F0", lty = 1)
    text(top_c$Stations_Included, bp,
         labels = sprintf(" %d stns (Med: %.1f m)", top_c$Stations_Included, top_c$Median_Snap_m),
         pos = 4, cex = 0.80, font = 2, col = "#1E293B")
  } else {
    plot.new()
  }

  invisible(out_png)
}

# ------------------------------------------------------------------------------
# Main Exported Function: `analyze_station_snapping()`
# ------------------------------------------------------------------------------
analyze_station_snapping <- function(
  catchment        = NULL,
  network          = "orn",
  source           = "rsn",
  include_excluded = TRUE,
  max_catch_plots  = 4L,
  out_dir          = NULL,
  gis_dir          = NULL,
  data_dir         = NULL,
  output_root      = NULL,
  verbose          = TRUE
) {
  g_dir <- if (!is.null(gis_dir) && nzchar(gis_dir)) gis_dir else rel_path(BASE_DIR, "gis")
  d_dir <- if (!is.null(data_dir) && nzchar(data_dir)) data_dir else rel_path(BASE_DIR, "data")
  o_root <- if (!is.null(output_root) && nzchar(output_root)) output_root else rel_path(BASE_DIR, "output")

  src_tag      <- resolve_source_tag(source)
  source_abbr  <- format_source_abbr(src_tag)
  source_label <- format_source_label(src_tag)
  src_ftag     <- format_source_folder_tag(src_tag)

  proc_dir <- resolve_network_processed_dir(BASE_DIR, network, data_dir = d_dir, gis_dir = g_dir, source = src_tag)
  net_meta <- load_network_metadata(proc_dir)

  if (is.null(out_dir) || !nzchar(out_dir)) {
    folder_name <- paste0("station_snapping_", net_meta$abbr)
    if (nzchar(src_ftag)) {
      folder_name <- paste0(folder_name, "_", src_ftag)
    }
    out_dir <- rel_path(o_root, folder_name)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    cat("======================================================================\n")
    cat(sprintf("  Station-to-Reach Snapping Diagnostics [%s | Source: %s]\n", net_meta$label, source_label))
    cat("======================================================================\n")
  }

  built      <- build_station_snapping_table(proc_dir, include_excluded = include_excluded)
  snap_df    <- built$snap_df
  reaches_df <- built$reaches_df
  catch_lu   <- built$catch_lu
  catch_sum  <- compute_catchment_snapping_summary(snap_df, reaches_df, catch_lu)

  # Write full station snapping table and catchment summary table
  all_csv_path   <- file.path(out_dir, "station_snapping_all.csv")
  catch_csv_path <- file.path(out_dir, "catchment_snapping_summary.csv")
  write.csv(snap_df, all_csv_path, row.names = FALSE)
  write.csv(catch_sum, catch_csv_path, row.names = FALSE)

  # Render network-wide overview figure
  overview_png <- file.path(out_dir, sprintf("fig_snap_overview_%s.png", net_meta$abbr))
  plot_network_snapping_overview(snap_df, catch_sum, paste0(net_meta$short_label, " [", source_abbr, "]"), overview_png)

  # Resolve which catchments to plot & export individual catchment tables for
  target_cids <- resolve_snapping_catchments(catchment, catch_sum)
  if (length(target_cids) == 0) {
    target_cids <- head(catch_sum$CatchID[catch_sum$Stations_Included > 0], max_catch_plots)
  }

  catch_outputs <- list()
  for (cid in target_cids) {
    cs_row <- catch_sum[catch_sum$CatchID == cid, ]
    r_slug <- if (nrow(cs_row) > 0) gsub("[^a-zA-Z0-9]+", "_", tolower(cs_row$Main_River[1])) else paste0("cid_", cid)
    r_slug <- gsub("^_+|_+$", "", r_slug)

    sub_s <- snap_df[snap_df$CatchID == cid, ]
    c_csv <- file.path(out_dir, sprintf("station_snapping_catchment_%d_%s.csv", cid, r_slug))
    c_png <- file.path(out_dir, sprintf("fig_snap_catchment_%d_%s.png", cid, r_slug))

    write.csv(sub_s, c_csv, row.names = FALSE)
    plot_catchment_snapping(cid, snap_df, reaches_df, catch_sum, net_meta$short_label, c_png)

    catch_outputs[[as.character(cid)]] <- list(
      CatchID    = cid,
      Main_River = if (nrow(cs_row) > 0) cs_row$Main_River[1] else paste0("Catchment ", cid),
      csv        = c_csv,
      png        = c_png,
      n_stations = nrow(sub_s)
    )
    if (verbose) {
      cat(sprintf("  + Catchment %d (%s): wrote %d stations to %s & plot %s\n",
                  cid, catch_outputs[[as.character(cid)]]$Main_River, nrow(sub_s),
                  format_rel_display_path(c_csv, BASE_DIR),
                  format_rel_display_path(c_png, BASE_DIR)))
    }
  }

  # Write Markdown summary report
  rep_path <- file.path(out_dir, "station_snapping_report.md")
  d_inc <- snap_df$snap_dist_m[snap_df$included_250m == 1L]
  md_lines <- c(
    sprintf("# Station-to-Reach Snapping Diagnostics Report — %s [%s]", net_meta$label, source_abbr),
    "",
    "## 1. Network-Wide Snapping Summary",
    "",
    sprintf("- **Water Network**: %s (`%s`)", net_meta$label, net_meta$abbr),
    sprintf("- **Data Source**: %s (`%s`)", source_label, source_abbr),
    sprintf("- **Total Network Reaches**: `%s` across `%d` catchments", format(nrow(reaches_df), big.mark = ","), nrow(catch_sum)),
    sprintf("- **Total %s Stations Evaluated**: `%d`", source_abbr, nrow(snap_df)),
    sprintf("- **Stations Snapped Within <= 250 m (`included_250m == 1`)**: `%d` (`%.1f%%`) across `%d` unique reaches",
            length(d_inc), 100 * length(d_inc) / max(1, nrow(snap_df)),
            length(unique(snap_df$snapped_RivID[snap_df$included_250m == 1L]))),
    sprintf("- **Stations Excluded (`> 250 m`)**: `%d` (`%.1f%%`)",
            sum(snap_df$included_250m == 0L), 100 * sum(snap_df$included_250m == 0L) / max(1, nrow(snap_df))),
    sprintf("- **Included Snapping Distance (`m`)**: Min = `%.2f m`, Q25 = `%.2f m`, **Median = `%.2f m`**, Mean = `%.2f m`, Q75 = `%.2f m`, P95 = `%.2f m`, Max = `%.2f m`",
            min(d_inc, na.rm = TRUE), quantile(d_inc, 0.25, na.rm = TRUE),
            median(d_inc, na.rm = TRUE), mean(d_inc, na.rm = TRUE),
            quantile(d_inc, 0.75, na.rm = TRUE), quantile(d_inc, 0.95, na.rm = TRUE),
            max(d_inc, na.rm = TRUE)),
    sprintf("- **Precision Breakdown (Included Stations)**: `%.1f%%` within `<= 5 m`, `%.1f%%` within `<= 10 m`, `%.1f%%` within `<= 25 m`, `%.1f%%` within `<= 50 m`",
            100 * mean(d_inc <= 5), 100 * mean(d_inc <= 10), 100 * mean(d_inc <= 25), 100 * mean(d_inc <= 50)),
    "",
    "## 2. Top Catchments by Snapped Monitoring Stations",
    "",
    "| CatchID | Main River | RBD | Reaches | Included Stns | Excluded Stns | Median Snap (m) | Mean Snap (m) | P95 Snap (m) | Max Snap (m) | <= 10 m (%) |",
    "| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
  )
  top_show <- head(catch_sum[catch_sum$Stations_Total > 0, ], 20)
  for (i in seq_len(nrow(top_show))) {
    r <- top_show[i, ]
    md_lines <- c(md_lines, sprintf(
      "| %d | %s | %s | %s | %d | %d | %.2f | %.2f | %.2f | %.2f | %.1f%% |",
      r$CatchID, r$Main_River, r$RBD_Name, format(r$Total_Reaches, big.mark = ","),
      r$Stations_Included, r$Stations_Excluded, r$Median_Snap_m, r$Mean_Snap_m,
      r$P95_Snap_m, r$Max_Snap_m, r$Pct_Within_10m
    ))
  }

  md_lines <- c(
    md_lines,
    "",
    "## 3. Sample of Station-to-Reach Snapping Table (`station_snapping_all.csv`)",
    "",
    "| Station ID | CatchID | Catchment | Station (X, Y) [m] | Snapped RivID | Snapped Point (X, Y) [m] | Reach Midpoint (X, Y) [m] | Snap Dist (m) | Status |",
    "| :--- | ---: | :--- | :--- | ---: | :--- | :--- | ---: | :--- |"
  )
  # Show sample rows from the first analyzed catchment
  samp_s <- head(snap_df[snap_df$CatchID %in% target_cids, ], 25)
  for (i in seq_len(nrow(samp_s))) {
    s <- samp_s[i, ]
    md_lines <- c(md_lines, sprintf(
      "| `%s` | %d | %s | `(%.1f, %.1f)` | `%d` | `(%.1f, %.1f)` | `(%.1f, %.1f)` | **%.2f m** | %s |",
      s$site_id, s$CatchID, s$Catchment_Name,
      s$station_X, s$station_Y, s$snapped_RivID,
      s$snap_reach_X, s$snap_reach_Y, s$reach_mid_X, s$reach_mid_Y,
      s$snap_dist_m, s$snap_status
    ))
  }

  md_lines <- c(
    md_lines,
    "",
    "## 4. Generated Files & Figures",
    "",
    sprintf("- **Full Station Snapping Table**: `station_snapping_all.csv` (`%d` rows)", nrow(snap_df)),
    sprintf("- **Catchment Summary Table**: `catchment_snapping_summary.csv` (`%d` catchments)", nrow(catch_sum)),
    sprintf("- **Network Overview Figure**: `%s`", basename(overview_png))
  )
  for (co in catch_outputs) {
    md_lines <- c(md_lines, sprintf(
      "- **Catchment %d (%s)**: Table `%s` (`%d` stations) & Figure `%s`",
      co$CatchID, co$Main_River, basename(co$csv), co$n_stations, basename(co$png)
    ))
  }
  writeLines(md_lines, rep_path)

  if (verbose) {
    cat("----------------------------------------------------------------------\n")
    cat(sprintf("  Full Snapping Table : %s (%d stations)\n", format_rel_display_path(all_csv_path, BASE_DIR), nrow(snap_df)))
    cat(sprintf("  Catchment Summary   : %s (%d catchments)\n", format_rel_display_path(catch_csv_path, BASE_DIR), nrow(catch_sum)))
    cat(sprintf("  Overview Plot       : %s\n", format_rel_display_path(overview_png, BASE_DIR)))
    cat(sprintf("  Markdown Report     : %s\n", format_rel_display_path(rep_path, BASE_DIR)))
    cat("======================================================================\n")
  }

  invisible(list(
    snap_table        = snap_df,
    catchment_summary = catch_sum,
    catchment_outputs = catch_outputs,
    overview_png      = overview_png,
    report_md         = rep_path,
    out_dir           = out_dir
  ))
}

# ==============================================================================
# CLI Entrypoint (only runs when invoked directly via `Rscript`)
# ==============================================================================
if (sys.nframe() == 0L) {
  args    <- commandArgs(trailingOnly = TRUE)
  dir_cfg <- parse_custom_dir_args(args, base_dir = BASE_DIR)

  src_cli <- parse_source_arg(args)
  if ("--list-sources" %in% args) {
    print_available_sources()
    quit(save = "no", status = 0)
  }

  net_arg <- grep("^--network=", args, value = TRUE)
  network <- if (length(net_arg) > 0) sub("^--network=", "", net_arg[1]) else "orn"

  if (any(args %in% c("--help", "-h"))) {
    cat("Usage:\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R [OPTIONS] [CATCHMENTS...]\n\n")
    cat("Options:\n")
    cat("  --source=rsn|ssn|rss                      Monitoring network dataset (default: rsn)\n")
    cat("  --list-sources                            List supported station data sources\n")
    cat("  --list, --list-catchments, --list-basins  List all catchments and snapping counts\n")
    cat("  --catchment=<ID_or_Name[,ID2...]>         Catchment ID(s), river name(s), or 'all'\n")
    cat("  --basin=<ID_or_Name[,ID2...]>             Alias for --catchment\n")
    cat("  --network=<name>                          Water network ('orn', 'euhydro_thames', etc.)\n")
    cat("  --only-included                           Exclude stations > 250 m threshold\n")
    cat("  --out=<dir>                               Custom output directory for snapping tables/plots\n")
    cat("  --gis-dir=<dir> --data-dir=<dir> --output-dir=<dir>\n\n")
    cat("Examples:\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R --list-catchments\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R --source=ssn --list-catchments\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R --source=rss --catchment=Thames\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R --catchment=Humber,Thames,Severn,Eden\n")
    cat("  Rscript scripts/06_station_snapping_diagnostics.R --network=euhydro_thames --catchment=306\n")
    quit(save = "no", status = 0)
  }

  proc_dir <- resolve_network_processed_dir(BASE_DIR, network, data_dir = dir_cfg$data_dir, gis_dir = dir_cfg$gis_dir, source = src_cli)
  net_meta <- load_network_metadata(proc_dir)

  if (any(args %in% c("--list", "--list-catchments", "--list-basins"))) {
    built <- build_station_snapping_table(proc_dir, include_excluded = TRUE)
    csum  <- compute_catchment_snapping_summary(built$snap_df, built$reaches_df, built$catch_lu)
    cat(sprintf("\nAvailable Catchments in %s (%d catchments):\n", net_meta$label, nrow(csum)))
    cat(sprintf("%-8s | %-24s | %-16s | %-8s | %-9s | %-9s | %-10s | %-10s\n",
                "CatchID", "Main_River", "RBD_Name", "Reaches", "Inc_Stns", "Exc_Stns", "Med_Snap_m", "Max_Snap_m"))
    cat(paste(rep("-", 108), collapse = ""), "\n")
    for (i in seq_len(min(50L, nrow(csum)))) {
      r <- csum[i, ]
      cat(sprintf("%-8d | %-24s | %-16s | %-8d | %-9d | %-9d | %-10.2f | %-10.2f\n",
                  r$CatchID, substr(r$Main_River, 1, 24), substr(r$RBD_Name, 1, 16),
                  r$Total_Reaches, r$Stations_Included, r$Stations_Excluded,
                  ifelse(is.na(r$Median_Snap_m), 0, r$Median_Snap_m),
                  ifelse(is.na(r$Max_Snap_m), 0, r$Max_Snap_m)))
    }
    quit(save = "no", status = 0)
  }

  catch_args <- grep("^--(catchment|catchments|basin|basins)=", args, value = TRUE)
  catch_specs <- if (length(catch_args) > 0) {
    unlist(strsplit(sub("^--(catchment|catchments|basin|basins)=", "", catch_args), ",", fixed = TRUE))
  } else {
    character(0)
  }
  pos_args <- args[!grepl("^--", args)]
  if (length(pos_args) > 0) {
    catch_specs <- c(catch_specs, unlist(strsplit(pos_args, ",", fixed = TRUE)))
  }

  only_inc <- ("--only-included" %in% args)
  out_arg  <- grep("^--out=", args, value = TRUE)
  custom_out <- if (length(out_arg) > 0) sub("^--out=", "", out_arg[1]) else NULL

  analyze_station_snapping(
    catchment        = if (length(catch_specs) > 0) catch_specs else NULL,
    network          = network,
    source           = src_cli,
    include_excluded = !only_inc,
    out_dir          = custom_out,
    gis_dir          = dir_cfg$gis_dir,
    data_dir         = dir_cfg$data_dir,
    output_root      = dir_cfg$out_dir,
    verbose          = TRUE
  )
}
