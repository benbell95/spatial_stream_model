#!/usr/bin/env Rscript
# Automated unit tests for R modules (s3n_config, s3n_transforms, s3n_covariates,
# s3n_plot_functions, s3n_stats_functions, 06 compute_snapped_reach_point) and
# C++ neighbor / measurement-error prediction kernel (src/s3n_fast_neighbors.so).

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(sprintf("ASSERTION FAILED: %s", msg), call. = FALSE)
}

assert_near <- function(a, b, tol = 1e-4, msg = "") {
  ok <- all(is.finite(a) & is.finite(b) & abs(a - b) <= tol)
  if (!ok) {
    stop(sprintf("ASSERTION FAILED (%s): got %s, expected %s (tol=%g)",
                 msg, paste(a, collapse = ","), paste(b, collapse = ","), tol),
         call. = FALSE)
  }
}

cat("=== 1. Testing R Module Loading & Determinand/Source Configs ===\n")
source("scripts/s3n_plot_functions.R")
source("scripts/s3n_stats_functions.R")

cat_dets <- load_determinand_catalog(".")
assert_true(length(cat_dets) >= 50L, "Expected >= 50 determinands in configs/determinands.json")
m_no3 <- get_determinand_meta("0117")
assert_true(grepl("Nitrate", m_no3[["name"]], ignore.case = TRUE), "0117 should resolve to Nitrate")
m_whpt <- get_determinand_meta("whpt_aspt")
assert_true(m_whpt[["unit"]] == "ASPT", "whpt_aspt unit should be ASPT")

assert_true(resolve_source_tag("rsn") == "rsn", "resolve_source_tag rsn")
assert_true(resolve_source_tag("waterblitz_tobit+rsn") == "rsn+waterblitz_tobit", "combo order sorting")
assert_true(resolve_source_tag("rsn,waterblitz-tobit") == "rsn+waterblitz_tobit", "hyphen/comma normalization")
assert_true(grepl("RSN", format_source_abbr("rsn+waterblitz_tobit")), "format_source_abbr contains RSN")

cat("  [PASS] Configs and source tag helpers\n")

cat("=== 2. Testing Transformations, Variance Delta-Method & Normalisation ===\n")
y_pos <- c(0.5, 1.0, 2.5, 5.0, 10.0)
y_prop <- c(0.05, 0.25, 0.50, 0.75, 0.95)

for (tt in c("log", "log10", "log1p", "sqrt", "cbrt", "fourthrt", "reciprocal", "boxcox", "yeojohnson", "none")) {
  tr <- build_transformer(tt, y_obs = y_pos, lambda = if (tt %in% c("boxcox", "yeojohnson")) 0.5 else NA_real_)
  z <- tr$fwd(y_pos)
  y_back <- tr$inv(z)
  assert_near(y_back, y_pos, tol = 1e-4, msg = paste("roundtrip", tt))
  if (tt != "none") {
    v_sc <- calc_response_meas_variance(y_pos, rep(0.2, length(y_pos)), tr)
    assert_true(all(is.finite(v_sc) & v_sc >= 0), paste("delta variance finite/non-neg for", tt))
  }
}

for (tt in c("asin", "logit")) {
  tr <- build_transformer(tt, y_obs = y_prop)
  z <- tr$fwd(y_prop)
  y_back <- tr$inv(z)
  assert_near(y_back, y_prop, tol = 1e-4, msg = paste("roundtrip", tt))
  v_sc <- calc_response_meas_variance(y_prop, rep(0.02, length(y_prop)), tr)
  assert_true(all(is.finite(v_sc) & v_sc >= 0), paste("delta variance finite/non-neg for", tt))
}

# Aliases + normalisation columns
df_legacy <- data.frame(
  RivID = 1:5,
  has_obs = c(1L, 0L, 1L, 0L, 1L),
  pred_log_combined = c(1.0, 1.5, 2.0, 2.5, 3.0),
  obs_log_nitrate = c(1.1, NA_real_, 2.1, NA_real_, 2.9),
  nitrate_mean = c(3.0, NA_real_, 8.0, NA_real_, 18.0)
)
df_norm <- normalize_pred_df(df_legacy)
assert_true("pred_model_sc" %in% names(df_norm), "normalize_pred_df aliases pred_log_combined -> pred_model_sc")
assert_true("obs_model_sc" %in% names(df_norm), "normalize_pred_df aliases obs_log_nitrate -> obs_model_sc")
assert_true("obs_value" %in% names(df_norm), "normalize_pred_df aliases nitrate_mean -> obs_value")

df_norm <- add_normalised_columns(df_norm, ref_df = df_norm)
assert_near(mean(df_norm$pred_zscore), 0.0, tol = 1e-4, msg = "z-score mean ~ 0")
assert_near(sd(df_norm$pred_zscore), 1.0, tol = 1e-3, msg = "z-score sd ~ 1")
assert_near(range(df_norm$pred_minmax), c(0.0, 1.0), tol = 1e-4, msg = "minmax range [0, 1]")
assert_true(!any(is.na(df_norm$obs_zscore[c(1L, 3L, 5L)])), "obs_zscore computed for observed reaches")
assert_true(all(is.na(df_norm$obs_zscore[c(2L, 4L)])), "obs_zscore NA for unobserved reaches")

cat("  [PASS] All 12 transforms, delta-method variances, aliases, and normalisation\n")

cat("=== 3. Testing Stats Helpers & Snapping Geometry ===\n")
mdw <- new_md_writer(2L)
mdw$add_line("Line %d", 1L)
mdw$add_line("Line %d", 2L)
mdw$add_line("Line %d", 3L)  # triggers capacity doubling
assert_true(identical(mdw$get_lines(), c("Line 1", "Line 2", "Line 3")), "new_md_writer growth")

# Parse only compute_snapped_reach_point without invoking top-level CLI in 06_station_snapping_diagnostics.R
lines06 <- readLines("scripts/06_station_snapping_diagnostics.R", warn = FALSE)
idx_fn <- grep("^compute_snapped_reach_point <- function", lines06)[1]
idx_end <- grep("^# Helper: Load and merge", lines06)[1] - 1L
env06 <- new.env(parent = globalenv())
eval(parse(text = lines06[idx_fn:idx_end]), envir = env06)

snap_out <- env06$compute_snapped_reach_point(
  st_x = c(50.0, 50.0),
  st_y = c(10.0, 0.0),    # second point lies directly on chord (v_norm == 0)
  sx   = c(0.0,  0.0),
  sy   = c(0.0,  0.0),
  mx   = c(50.0, 50.0),
  my   = c(0.0,  0.0),
  ex   = c(100.0, 100.0),
  ey   = c(0.0,  0.0),
  snap_dist_m = c(10.0, 0.0)
)
assert_near(snap_out$x, c(50.0, 50.0), tol = 1e-4, msg = "snapped X")
assert_near(snap_out$y, c(0.0, 0.0), tol = 1e-4, msg = "snapped Y including on-chord case")

cat("  [PASS] Buffered markdown writer and station chord projection\n")

cat("=== 4. Testing C++ Shared Library (s3n_fast_neighbors.so) ===\n")
if (!is.loaded("C_compute_s3n_neighbors")) {
  dyn.load("src/s3n_fast_neighbors.so")
}
assert_true(is.loaded("C_build_tree_tin_tout"), "C_build_tree_tin_tout loaded")
assert_true(is.loaded("C_compute_s3n_neighbors"), "C_compute_s3n_neighbors loaded")
assert_true(is.loaded("C_s3n_predict_meas_error"), "C_s3n_predict_meas_error loaded")

# 6-reach bifurcating tree (first 4 reaches observed, last 2 unobserved):
# Index (0-based): 0..5. Topology:
#   mouth = 2 (CatchID=1, d2m=2.0, up=1.5, parent_dn=-1)
#   trib1 = 0 -> dn=2 (CatchID=1, d2m=4.0, up=3.5, afv=log(0.6))
#   trib2 = 1 -> dn=2 (CatchID=1, d2m=4.0, up=3.5, afv=log(0.4))
#   head1 = 3 -> dn=0 (CatchID=1, d2m=6.0, up=5.5, afv=log(0.6))
#   pred1 = 4 -> dn=0 (CatchID=1, d2m=5.0, up=4.5, afv=log(0.6))
#   pred2 = 5 -> dn=1 (CatchID=1, d2m=5.0, up=4.5, afv=log(0.4))
n_obs <- 4L
q_tot <- 6L
m_nbr <- 3L
parent_0 <- c(2L, 2L, -1L, 0L, 0L, 1L)
cid      <- rep(1L, q_tot)
d2m      <- c(4.0, 4.0, 2.0, 6.0, 5.0, 5.0)
up_dist  <- c(3.5, 3.5, 1.5, 5.5, 4.5, 4.5)
log_afv  <- c(log(0.6), log(0.4), 0.0, log(0.6), log(0.6), log(0.4))
mx       <- c(1000, 1200, 1100, 900, 950, 1250)
my       <- c(2000, 2000, 1000, 3000, 2500, 2500)

tt <- .Call("C_build_tree_tin_tout", as.integer(parent_0))
assert_true(length(tt[[1]]) == q_tot && length(tt[[2]]) == q_tot, "tin/tout vector lengths")
assert_true(all(tt[[2]] >= tt[[1]]), "tout >= tin for every tree node")

nb_tu <- .Call("C_compute_s3n_neighbors",
               as.integer(n_obs), as.integer(m_nbr), as.integer(parent_0),
               as.integer(cid), as.integer(tt[[1]]), as.integer(tt[[2]]),
               as.double(d2m), as.double(up_dist), as.double(log_afv),
               as.double(mx), as.double(my), 0L)
nb_td <- .Call("C_compute_s3n_neighbors",
               as.integer(n_obs), as.integer(m_nbr), as.integer(parent_0),
               as.integer(cid), as.integer(tt[[1]]), as.integer(tt[[2]]),
               as.double(d2m), as.double(up_dist), as.double(log_afv),
               as.double(mx), as.double(my), 1L)

assert_true(!is.null(nb_tu$obs_neighbors) && !is.null(nb_tu$pred_neighbors), "TU neighbor struct complete")
assert_true(!is.null(nb_td$obs_neighbors) && !is.null(nb_td$pred_neighbors), "TD neighbor struct complete")

X_all <- cbind(1.0, (up_dist - mean(up_dist)) / sd(up_dist))
X_obs <- X_all[1:n_obs, , drop = FALSE]
y_obs <- matrix(c(2.1, 1.4, 1.8, 2.5), ncol = 1L)
v_meas <- c(0.05, 0.05, 0.02, 0.08)
beta <- c(1.9, 0.3)
theta <- c(0.4, 0.05, 0.2)

pr_me <- .Call("C_s3n_predict_meas_error",
               X_obs, y_obs, v_meas, X_all, beta, theta, as.integer(m_nbr),
               nb_tu$pred_neighbors$nnIndx, nb_tu$pred_neighbors$nnDist, nb_tu$pred_neighbors$nnWght,
               nb_tu$obsobs_dist, nb_tu$obsobs_wt)
assert_true(length(pr_me$prediction) == q_tot, "ME prediction length == q_tot")
assert_true(all(is.finite(pr_me$prediction)), "ME predictions finite")
assert_true(all(pr_me$prediction.ci[, 2] > pr_me$prediction.ci[, 1]), "ME prediction CI upper > lower")

cat("  [PASS] C++ tree Euler tour, TU/TD neighbor builders, and ME kriging predictor\n")
cat("\nALL R & C++ UNIT TESTS PASSED!\n")
