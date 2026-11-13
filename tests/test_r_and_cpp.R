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

# Canonical prediction schema + normalisation columns
df_canon <- data.frame(
  RivID         = 1:5,
  has_obs       = c(1L, 0L, 1L, 0L, 1L),
  pred_model_sc = c(1.0, 1.5, 2.0, 2.5, 3.0),
  obs_model_sc  = c(1.1, NA_real_, 2.1, NA_real_, 2.9),
  obs_value     = c(3.0, NA_real_, 8.0, NA_real_, 18.0)
)
df_norm <- add_normalised_columns(df_canon, ref_df = df_canon)
assert_near(mean(df_norm$pred_zscore), 0.0, tol = 1e-4, msg = "z-score mean ~ 0")
assert_near(sd(df_norm$pred_zscore), 1.0, tol = 1e-3, msg = "z-score sd ~ 1")
assert_near(range(df_norm$pred_minmax), c(0.0, 1.0), tol = 1e-4, msg = "minmax range [0, 1]")
assert_true(!any(is.na(df_norm$obs_zscore[c(1L, 3L, 5L)])), "obs_zscore computed for observed reaches")
assert_true(all(is.na(df_norm$obs_zscore[c(2L, 4L)])), "obs_zscore NA for unobserved reaches")

# Yeo-Johnson with negative values (tests Jacobian sign), asin percentage scale=100, and unknown transform guard
y_neg <- c(-5.0, -1.0, 0.0, 1.5, 6.0)
tr_yj <- build_transformer("yeojohnson", y_obs = y_neg)
assert_true(is.finite(tr_yj$lambda), "YJ profile ML lambda finite with negative values")
assert_near(tr_yj$inv(tr_yj$fwd(y_neg)), y_neg, tol = 1e-4, msg = "YJ negative-value roundtrip")
v_yj <- calc_response_meas_variance(y_neg, rep(0.25, length(y_neg)), tr_yj)
assert_true(all(is.finite(v_yj) & v_yj >= 0), "YJ delta variance finite/non-neg")

y_pct <- c(5.0, 25.0, 50.0, 75.0, 95.0)
tr_asin_pct <- build_transformer("asin", y_obs = y_pct)
assert_near(tr_asin_pct$scale, 100.0, tol = 1e-6, msg = "asin auto-detects scale = 100 for percentages")
assert_near(tr_asin_pct$inv(tr_asin_pct$fwd(y_pct)), y_pct, tol = 1e-4, msg = "asin percentage roundtrip")
v_asin_pct <- calc_response_meas_variance(y_pct, rep(4.0, length(y_pct)), tr_asin_pct)
assert_true(all(is.finite(v_asin_pct) & v_asin_pct >= 0), "asin percentage delta variance finite/non-neg")

bad_err <- tryCatch({ build_transformer("definitely_invalid_transform"); FALSE }, error = function(e) TRUE)
assert_true(bad_err, "build_transformer rejects unknown transform names")

cat("  [PASS] All 12 transforms, delta-method variances, canonical schemas, and normalisation\n")

cat("=== 2b. Testing Covariate Specs & Design Matrix Imputation ===\n")
cov_aq <- apply_covariate_specs("+aquifer", c("Elevation_z", "Slope_z", "Is_Chalk_Limestone", "Is_Sandstone"))
assert_true("Is_Aquifer" %in% cov_aq && !("Is_Chalk_Limestone" %in% cov_aq) && !("Is_Sandstone" %in% cov_aq),
            "+aquifer substitutes Is_Chalk_Limestone/Is_Sandstone with Is_Aquifer")
cov_mixed <- apply_covariate_specs(c("elevation", "chalk_limestone", "+aquifer"), c("Slope_z"))
assert_true(identical(cov_mixed, c("Elevation_z", "Is_Aquifer")),
            "mixed explicit + modifier spec resolves cleanly")
df_cov <- data.frame(
  Elevation          = c(50, NA_real_, 150, 200),
  Slope              = c(1.0, 2.0, 3.0, 4.0),
  US_Accum_km        = c(5.0, 10.0, 20.0, 40.0),
  Is_Chalk_Limestone = c(1L, 0L, 1L, 0L),
  Is_Sandstone       = c(0L, 1L, 0L, 0L),
  Is_Peat            = c(0L, 0L, 1L, 0L),
  Is_Aquifer         = c(1L, 0L, 1L, 0L),
  mid_X              = c(400000, 410000, 420000, 430000),
  mid_Y              = c(300000, 310000, 320000, 330000),
  has_obs            = c(1L, 0L, 1L, 1L)
)
dm <- build_design_matrix(df_cov, include_peat = TRUE, var_specs = c("Elevation_z", "Slope_z"))
assert_true("Elevation_z" %in% colnames(dm) && all(is.finite(dm)),
            "build_design_matrix imputes isolated NA without dropping covariate")
cat("  [PASS] Covariate spec resolution and design matrix NA imputation\n")

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

cat(sprintf("=== 4. Testing C++ Shared Library (s3n_fast_neighbors%s) ===\n", .Platform$dynlib.ext))
if (!is.loaded("C_compute_s3n_neighbors")) {
  dyn.load(file.path("src", paste0("s3n_fast_neighbors", .Platform$dynlib.ext)))
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

cat("=== 5. Testing New Core S3N Improvements (Multi-Kernel, Physical AFV, Joint TU+TD+Euc NNGP) ===\n")
source("scripts/s3n_model_functions.R")
assert_true(is.loaded("C_compute_log_afv"), "C_compute_log_afv loaded")
assert_true(is.loaded("C_s3n_nngp_loglik"), "C_s3n_nngp_loglik loaded")
assert_true(is.loaded("C_s3n_joint_nngp"), "C_s3n_joint_nngp loaded")

# 5a. Physical AFV Confluence Weighting (C_compute_log_afv & compute_reach_log_afv)
w_accum <- c(6.0, 4.0, 10.0, 3.0, 3.0, 4.0)
afv_cpp <- .Call("C_compute_log_afv", as.integer(parent_0), as.double(w_accum), 1L)
assert_near(afv_cpp[3], 0.0, tol = 1e-8, msg = "Outlet reach (index 2) has log_AFV == 0")
assert_near(exp(afv_cpp[1]) + exp(afv_cpp[2]), 1.0, tol = 1e-8, msg = "Confluence child AFVs sum to 1.0 at mouth")
assert_near(afv_cpp[1], log(0.6), tol = 1e-8, msg = "Child 0 log_AFV == log(0.6)")
assert_near(afv_cpp[2], log(0.4), tol = 1e-8, msg = "Child 1 log_AFV == log(0.4)")
assert_near(exp(afv_cpp[4]) + exp(afv_cpp[5]), exp(afv_cpp[1]), tol = 1e-8, msg = "Upstream split conserves parent AFV")

df_toy <- data.frame(
  RivID = 1:6,
  parent_dn = c(3L, 3L, 0L, 1L, 1L, 2L),
  US_Accum_km = w_accum,
  Shreve = c(2L, 1L, 3L, 1L, 1L, 1L),
  Elevation = c(120, 80, 20, 250, 220, 150),
  Is_Aquifer = c(1L, 0L, 0L, 1L, 1L, 0L)
)
afv_len <- compute_reach_log_afv(df_toy, "length", parent_0 = parent_0)
afv_shr <- compute_reach_log_afv(df_toy, "shreve", parent_0 = parent_0)
afv_hyd <- compute_reach_log_afv(df_toy, "hydro_yield", parent_0 = parent_0)
assert_near(afv_len[1], log(0.6), tol = 1e-6, msg = "length AFV scheme matches US_Accum_km split")
assert_near(afv_shr[1], log(2.0 / 3.0), tol = 1e-6, msg = "shreve AFV scheme matches 2/(2+1) split")
assert_true(afv_hyd[1] > afv_len[1], "hydro_yield AFV scheme boosts high-elevation aquifer tributary")

# 5b. Multi-Kernel NNGP Log-Likelihood & Prediction (Exponential, Spherical, Mariah, Linear, Auto)
for (kid in 0:3) {
  ll_res <- .Call("C_s3n_nngp_loglik", X_obs, as.double(y_obs), nb_tu$obs_neighbors, 0.15, 0.25, as.integer(kid))
  assert_true(is.finite(ll_res$neg2loglik), paste0("C_s3n_nngp_loglik finite for kernel_id=", kid))
  assert_true(length(ll_res$Beta) == ncol(X_obs), paste0("GLS Beta length for kernel_id=", kid))
  assert_true(ll_res$sigma.sq > 0, paste0("Positive sigma.sq for kernel_id=", kid))

  pr_k <- .Call("C_s3n_predict_meas_error",
                X_obs, y_obs, v_meas, X_all, beta, c(0.4, 0.05, 0.25, kid), as.integer(m_nbr),
                nb_tu$pred_neighbors$nnIndx, nb_tu$pred_neighbors$nnDist, nb_tu$pred_neighbors$nnWght,
                nb_tu$obsobs_dist, nb_tu$obsobs_wt)
  assert_true(all(is.finite(pr_k$prediction)) && all(pr_k$se > 0), paste0("Prediction finite for kernel_id=", kid))
}

fit_auto <- fit_nngp_cpp(X_obs, as.double(y_obs), nb_tu$obs_neighbors, var(as.double(y_obs)), kernel_name = "auto")
assert_true(fit_auto$kernel_name %in% c("exponential", "spherical", "mariah", "linear"), "Auto kernel selection valid")
assert_true(length(fit_auto$auto_aic_table) == 4L, "Auto kernel AIC table has 4 kernels")

# 5c. Joint Additive NNGP (Tail-Up + Tail-Down + Euclidean)
params_joint <- c(0.25, 0.15, 0.05, 0.04, 0.20, 0.15, 0.10, 0.0)
jl_res <- .Call("C_s3n_joint_nngp",
                X_obs, as.double(y_obs), v_meas, X_obs, numeric(ncol(X_obs)),
                nb_tu, nb_td, as.double(mx), as.double(my),
                as.integer(m_nbr), as.double(params_joint), 0L)
assert_true(is.finite(jl_res$neg2loglik) && is.finite(jl_res$neg2loglik_exact), "Joint NNGP log-likelihood finite")
assert_true(all(is.finite(jl_res$Beta)), "Joint NNGP GLS Beta finite")

jp_res <- .Call("C_s3n_joint_nngp",
                X_obs, as.double(y_obs), v_meas, X_all, as.double(jl_res$Beta),
                nb_tu, nb_td, as.double(mx), as.double(my),
                as.integer(m_nbr), as.double(params_joint), 1L)
assert_true(length(jp_res$prediction) == q_tot && all(is.finite(jp_res$prediction)), "Joint NNGP predictions finite")
assert_true(all(jp_res$se > 0) && all(jp_res$prediction.ci[, 2] > jp_res$prediction.ci[, 1]), "Joint NNGP SE & CI valid")

fit_td_toy <- fit_nngp_cpp(X_obs, as.double(y_obs), nb_td$obs_neighbors, var(as.double(y_obs)), kernel_name = "mariah")
jf_full <- fit_and_predict_joint_s3n(
  X_obs = X_obs, y_vec_obs = as.double(y_obs), v_meas = v_meas, has_me = TRUE,
  X_all = X_all, nb_tu = nb_tu, nb_td = nb_td, mid_x = mx, mid_y = my,
  m_use = m_nbr, fit_tu = fit_auto, fit_td = fit_td_toy, kernel_id = 2L
)
assert_true(all(is.finite(jf_full$prediction)) && length(jf_full$prediction) == q_tot, "fit_and_predict_joint_s3n complete")

# 5d. End-to-end fit_s3n_domain with auto-kernel, joint TU+TD+Euc kriging, and hydro_yield AFV
set.seed(42)
n_syn <- 24L
df_syn <- data.frame(
  RivID             = seq_len(n_syn),
  CatchID           = rep(c(1L, 2L), each = 12L),
  parent_dn         = c(0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L,
                        0L, 13L, 13L, 14L, 14L, 15L, 15L, 16L, 16L, 17L, 17L, 18L),
  Dist2Mth_km       = rep(c(0, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 8), 2),
  up_dist_mid_km    = rep(c(1, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 9), 2),
  obs_up_dist_km    = rep(c(1, 3, 3, 5, 5, 5, 5, 7, 7, 7, 7, 9), 2),
  US_Accum_km       = rep(c(30, 14, 15, 6, 7, 7, 7, 2, 3, 3, 3, 1), 2),
  log_AFV           = rep(0.0, n_syn),
  Shreve            = rep(c(6L, 3L, 3L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), 2),
  Strahler          = rep(c(3L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), 2),
  start_X           = seq(400000, by = 500, length.out = n_syn),
  start_Y           = seq(300000, by = 500, length.out = n_syn),
  mid_X             = seq(400250, by = 500, length.out = n_syn),
  mid_Y             = seq(300250, by = 500, length.out = n_syn),
  end_X             = seq(400500, by = 500, length.out = n_syn),
  end_Y             = seq(300500, by = 500, length.out = n_syn),
  Elevation         = seq(20, 250, length.out = n_syn),
  Slope             = rep(2.5, n_syn),
  Sinuosity         = rep(1.15, n_syn),
  SrcElevation      = seq(100, 330, length.out = n_syn),
  Is_Chalk_Limestone = rep(c(0L, 1L), length.out = n_syn),
  Is_Sandstone      = rep(c(1L, 0L), length.out = n_syn),
  Is_Aquifer        = rep(c(1L, 0L), length.out = n_syn),
  Is_Peat           = rep(0L, n_syn),
  has_obs           = c(rep(1L, 18L), rep(0L, 6L)),
  obs_val           = NA_real_
)
df_syn$obs_val[1:18] <- exp(1.5 - 0.003 * df_syn$Elevation[1:18] + 0.4 * df_syn$Is_Aquifer[1:18] + rnorm(18, 0, 0.1))

sink(tempfile())
dom_out <- fit_s3n_domain(
  df = df_syn, domain_id = "test", domain_label = "Synthetic Basin",
  trans_obj = build_transformer("log"), var_specs = c("Elevation", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "auto", combine_mode = "joint", afv_scheme = "hydro_yield"
)
sink()
assert_true(!is.null(dom_out), "fit_s3n_domain returned non-null result")
assert_true(nrow(dom_out$cv_summary) == 4L && all(is.finite(dom_out$cv_summary$RMSE)), "5-fold CV summary valid")
assert_true("S3N Joint (TU + TD + Euc)" %in% dom_out$param_df$Model, "Joint model parameters recorded in param_df")

cat("  [PASS] Multi-kernel NNGP, physical AFV confluence schemes, and Joint Additive TU+TD+Euc kriging\n")

cat("=== 6. Testing Non-Linear / Spline Fixed Effects (GAM + S3N) & Tree-Ensemble Regression Kriging (RF / GBM + S3N) ===\n")
assert_true(parse_trend_arg("--trend=gam") == "gam", "parse_trend_arg gam")
assert_true(parse_trend_arg("--gam-s3n") == "gam", "parse_trend_arg --gam-s3n alias")
assert_true(parse_trend_arg("--trend=spline") == "spline", "parse_trend_arg spline")
assert_true(parse_trend_arg("--trend=rf") == "rf", "parse_trend_arg rf")
assert_true(parse_trend_arg("--trend=gbm") == "gbm", "parse_trend_arg gbm")

X_syn_all <- build_design_matrix(df_syn, include_peat = FALSE, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"))
X_syn_obs <- X_syn_all[1:18, , drop = FALSE]
y_syn_obs <- log(df_syn$obs_val[1:18]) + 0.35 * sin(X_syn_obs[, "Elevation_z"] * 1.5)

for (t_mode in c("spline", "gam", "rf", "gbm")) {
  tf_res <- fit_nonlinear_trend(
    X_tr = X_syn_obs, y_tr_vec = y_syn_obs, X_eval = X_syn_all,
    trend_mode = t_mode, spline_df = 3L, n_trees = 30L
  )
  assert_true(tf_res$mode == t_mode, paste0("fit_nonlinear_trend mode=", t_mode))
  assert_true(length(tf_res$g_perp_tr) == 18L && all(is.finite(tf_res$g_perp_tr)), paste0("g_perp_tr finite for ", t_mode))
  assert_true(length(tf_res$g_perp_eval) == n_syn && all(is.finite(tf_res$g_perp_eval)), paste0("g_perp_eval finite for ", t_mode))
  # Verify orthogonality: X_syn_obs^T %*% g_perp_tr == 0
  orth_max <- max(abs(crossprod(X_syn_obs, tf_res$g_perp_tr)))
  assert_true(orth_max < 1e-6, paste0("g_perp_tr orthogonal to X_obs for ", t_mode))
  assert_true(!is.null(tf_res$info_df) && nrow(tf_res$info_df) > 0L, paste0("info_df populated for ", t_mode))
}

for (t_mode in c("gam", "rf", "gbm")) {
  sink(tempfile())
  dom_t <- fit_s3n_domain(
    df = df_syn, domain_id = "test", domain_label = "Synthetic Basin",
    trans_obj = build_transformer("log"), var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
    m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
    kernel = "exponential", combine_mode = "joint", afv_scheme = "length",
    trend = t_mode, spline_df = 3L, n_trees = 25L
  )
  sink()
  assert_true(!is.null(dom_t), paste0("fit_s3n_domain succeeded for trend=", t_mode))
  assert_true("pred_trend_sc" %in% names(dom_t$pred_df) && all(is.finite(dom_t$pred_df$pred_trend_sc)),
              paste0("pred_trend_sc exported for trend=", t_mode))
  assert_true(any(grepl("^S3N Trend", dom_t$param_df$Model)), paste0("S3N Trend recorded in param_df for trend=", t_mode))
}

cat("  [PASS] Spline, GAM + S3N, Random Forest + S3N, and Gradient Boosted Trees + S3N regression kriging\n")

# ==============================================================================
# 7. Validation Modes, Regulatory Exceedance Risk & Optimal Monitoring Design
# ==============================================================================
cat("[Test 7/7] Testing CV modes, exceedance probability, and C++ optimal design...\n")

assert_true(parse_cv_mode_arg("--cv-mode=catchment") == "catchment", "parse_cv_mode_arg catchment")
assert_true(parse_cv_mode_arg("--cv-mode=strahler") == "strahler", "parse_cv_mode_arg strahler")
assert_true(parse_cv_mode_arg("--cv-mode=spatial_block") == "spatial_block", "parse_cv_mode_arg spatial_block")
assert_true(abs(parse_threshold_arg("--threshold=11.3") - 11.3) < 1e-9, "parse_threshold_arg")
assert_true(identical(parse_threshold_arg("--threshold=uk"), "uk"), "parse_threshold_arg uk")
uk_no3 <- resolve_threshold_spec("uk", "above", "0117")
assert_true(abs(uk_no3$val - 11.3) < 1e-9 && uk_no3$dir == "above", "UK threshold preset 0117")
uk_do <- resolve_threshold_spec("uk", "above", "9924")
assert_true(abs(uk_do$val - 6.0) < 1e-9 && uk_do$dir == "below", "UK threshold preset 9924 auto-below")
assert_true(parse_threshold_dir_arg("--threshold-dir=below") == "below", "parse_threshold_dir_arg below")
assert_true(parse_design_sites_arg("--design-sites=5") == 5L, "parse_design_sites_arg")

df_obs_syn <- df_syn[df_syn$has_obs == 1L, ]
df_obs_syn$CatchID <- rep(c(101L, 102L, 103L, 104L), length.out = nrow(df_obs_syn))
df_obs_syn$Strahler <- rep(c(1L, 2L, 3L, 4L), length.out = nrow(df_obs_syn))

cv_rand <- build_cv_folds(df_obs_syn, cv_mode = "random", k_folds = 5L)
assert_true(cv_rand$actual_mode == "random" && length(unique(cv_rand$folds)) == 5L, "random CV 5 folds")

cv_cat <- build_cv_folds(df_obs_syn, cv_mode = "catchment", k_folds = 5L)
assert_true(cv_cat$actual_mode == "catchment", "catchment CV mode active")
# Verify each CatchID maps to exactly one fold
for (cid in unique(df_obs_syn$CatchID)) {
  assert_true(length(unique(cv_cat$folds[df_obs_syn$CatchID == cid])) == 1L,
              paste0("CatchID ", cid, " stays in a single fold"))
}

cv_str <- build_cv_folds(df_obs_syn, cv_mode = "strahler")
assert_true(cv_str$actual_mode == "strahler", "strahler CV mode active")
assert_true(all(cv_str$folds[df_obs_syn$Strahler >= 3L] == 0L), "strahler trains on mainstems >= 3")
assert_true(all(cv_str$folds[df_obs_syn$Strahler <= 2L] == 1L), "strahler tests on headwaters <= 2")

cv_sb <- build_cv_folds(df_obs_syn, cv_mode = "spatial_block", k_folds = 4L)
assert_true(cv_sb$actual_mode == "spatial_block" && length(unique(cv_sb$folds)) == 4L, "spatial_block CV 4 folds")

# Exceedance risk tests (above & below)
tr_log <- build_transformer("log")
ex_above <- compute_exceedance_risk(c(log(5.0), log(10.0), log(20.0)), c(0.1, 0.2, 0.1), 10.0, tr_log, dir = "above")
assert_true(ex_above$prob[1] < 0.05 && ex_above$tier[1] == "Confidently Compliant", "above tier low")
assert_true(abs(ex_above$prob[2] - 0.5) < 1e-6 && ex_above$tier[2] == "Uncertain / Priority Monitoring", "above tier mid")
assert_true(ex_above$prob[3] > 0.95 && ex_above$tier[3] == "Confidently Exceeding", "above tier high")

ex_below <- compute_exceedance_risk(c(log(5.0), log(20.0)), c(0.1, 0.1), 10.0, tr_log, dir = "below")
assert_true(ex_below$prob[1] > 0.95 && ex_below$tier[1] == "Confidently Exceeding", "below tier exceeding")
assert_true(ex_below$prob[2] < 0.05 && ex_below$tier[2] == "Confidently Compliant", "below tier compliant")

# End-to-end fit_s3n_domain with catchment CV + threshold + optimal monitoring design
df_syn_cat <- df_syn
df_syn_cat$CatchID <- rep(c(101L, 102L, 103L, 104L), length.out = nrow(df_syn_cat))
sink(tempfile())
dom_val <- fit_s3n_domain(
  df = df_syn_cat, domain_id = "test", domain_label = "Synthetic Basin",
  trans_obj = tr_log, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "exponential", combine_mode = "joint", afv_scheme = "length",
  trend = "linear", cv_mode = "catchment", threshold = 5.0,
  threshold_dir = "above", design_sites = 3L
)
sink()
assert_true(!is.null(dom_val), "fit_s3n_domain validation+risk+design succeeded")
assert_true(all(c("exceed_prob", "compliance_tier") %in% names(dom_val$pred_df)), "pred_df has exceedance columns")
assert_true(!is.null(dom_val$exceed_summary) && sum(dom_val$exceed_summary$Reach_Count) == nrow(df_syn_cat),
            "exceed_summary counts sum to total reaches")
assert_true(!is.null(dom_val$design_df) && nrow(dom_val$design_df) == 3L, "design_df ranked 3 sites")
assert_true(all(dom_val$design_df$Var_Reduced >= 0) && all(diff(dom_val$design_df$Cum_Pct_Reduced) >= -1e-9),
            "design_df cumulative variance reduction is non-decreasing")
assert_true(all(dom_val$cv_summary$CV_Mode == "catchment"), "cv_summary records CV_Mode=catchment")

# Verify exceedance risk plotting functions render cleanly to PNG
tmp_ex_eng <- tempfile(fileext = ".png")
tmp_ex_bas <- tempfile(fileext = ".png")
plot_fig_exceedance_england(dom_val$pred_df, "0117", "Nitrate as N", "mg/L", tmp_ex_eng,
                            threshold_val = 5.0, threshold_dir = "above")
plot_fig_exceedance_basins(dom_val$pred_df, dom_val$pred_df, dom_val$pred_df,
                           "0117", "Nitrate as N", "mg/L", tmp_ex_bas,
                           threshold_val = 5.0, threshold_dir = "above")
assert_true(file.exists(tmp_ex_eng) && file.info(tmp_ex_eng)$size > 1000, "exceedance england PNG rendered")
assert_true(file.exists(tmp_ex_bas) && file.info(tmp_ex_bas)$size > 1000, "exceedance basins PNG rendered")
tmp_opt_des <- tempfile(fileext = ".png")
plot_fig_optimal_design(dom_val$pred_df, dom_val$design_df, "Synthetic Basin",
                        "0117", "Nitrate as N", "mg/L", tr_log, tmp_opt_des)
assert_true(file.exists(tmp_opt_des) && file.info(tmp_opt_des)$size > 1000, "optimal design PNG rendered")
unlink(c(tmp_ex_eng, tmp_ex_bas, tmp_opt_des))

cat("  [PASS] Hydrological CV modes, regulatory exceedance tiers, plots, and C++ optimal monitoring design\n")

# ==============================================================================
# 8. Spatial Tobit EM Censored Kriging & Spatio-Temporal Borrowing / Trend Mapping
# ==============================================================================
cat("[Test 8/8] Testing Spatial Tobit EM censored kriging & Spatio-Temporal borrowing/trends...\n")

assert_true(parse_st_borrow_arg("--st-borrow") == TRUE, "parse_st_borrow_arg TRUE")
assert_true(parse_st_borrow_arg(character(0)) == FALSE, "parse_st_borrow_arg FALSE")

tn_res <- trunc_norm_below_stats(mu = c(0.0, 2.0), sd_val = c(1.0, 1.0), upper = c(0.5, 0.5))
assert_true(all(tn_res$mean <= 0.5 + 1e-9) && tn_res$mean[1] < 0.0, "trunc_norm_below_stats upper bound")
assert_true(all(tn_res$var > 0 & tn_res$var < 1.0), "trunc_norm_below_stats variance < unconditional")

# Spatial Tobit EM Censored Kriging on synthetic basin with 8 censored reaches
df_syn_em <- df_syn
df_syn_em$is_censored <- 0L
df_syn_em$lod_sc <- NA_real_
cens_rows <- which(df_syn_em$has_obs == 1L)[1:8]
df_syn_em$is_censored[cens_rows] <- 1L
df_syn_em$lod_sc[cens_rows] <- log(df_syn_em$obs_val[cens_rows] + 0.5)

sink(tempfile())
dom_em <- fit_s3n_domain(
  df = df_syn_em, domain_id = "test", domain_label = "Synthetic Tobit Basin",
  trans_obj = tr_log, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "exponential", combine_mode = "joint", afv_scheme = "length",
  censor_mode = "em"
)
sink()
assert_true(!is.null(dom_em), "fit_s3n_domain censor_mode=em succeeded")
assert_true(grepl("\\[Tobit-EM\\]", dom_em$cv_summary$Spatial_Method[1]), "Spatial_Method tagged [Tobit-EM]")
assert_true(any(dom_em$param_df$Model == "Spatial Tobit EM"), "Spatial Tobit EM logged in param_df")

# Spatio-Temporal Dynamic Borrowing & Trend Mapping across 3 synthetic periods
mk_slice <- function(shift, se_mult = 1.0) {
  d <- dom_val$pred_df
  d$pred_model_sc <- d$pred_model_sc + shift
  d$se_model_sc   <- pmax(0.05, d$se_model_sc * se_mult)
  d$obs_model_sc  <- ifelse(d$has_obs == 1L, d$obs_model_sc + shift, NA_real_)
  list(eng_pred = d, hum_pred = d, tha_pred = d, sev_pred = d)
}
sl_list <- list(mk_slice(0.0, 1.0), mk_slice(0.3, 1.5), mk_slice(0.6, 1.0))
st_out  <- apply_spatiotemporal_borrowing(sl_list, tr_log, threshold_val = 5.0, threshold_dir = "above")
assert_true(is.finite(st_out$rho_temporal) && st_out$rho_temporal >= 0.15, "temporal rho_t estimated")
assert_true(mean(st_out$slices[[2]]$eng_pred$se_model_sc) < mean(sl_list[[2]]$eng_pred$se_model_sc),
            "dynamic borrowing reduces SE in sparser middle slice")

tr_df <- compute_spatiotemporal_trends(st_out$slices, tr_log)
assert_true(!is.null(tr_df) && nrow(tr_df) == nrow(df_syn), "compute_spatiotemporal_trends row count")
assert_true(all(tr_df$trend_slope_sc > 0), "positive temporal trend detected across increasing slices")

tmp_fig6 <- tempfile(fileext = ".png")
plot_fig6_spatiotemporal_trends(tr_df, "year", c("Year 2021", "Year 2022", "Year 2023"),
                                "0117", "Nitrate as N", "mg/L", tr_log, tmp_fig6,
                                borrowed = TRUE, rho_t = st_out$rho_temporal)
assert_true(file.exists(tmp_fig6) && file.info(tmp_fig6)$size > 1000, "Figure 6 spatio-temporal trend PNG rendered")
unlink(tmp_fig6)

# Period-to-period transitions (Figure 7) for season (cyclic) & year (non-cyclic)
seas_list <- list(mk_slice(0.4, 1.0), mk_slice(0.1, 1.0), mk_slice(-0.3, 1.0), mk_slice(0.2, 1.0))
seas_lbls <- c("Season: winter", "Season: spring", "Season: summer", "Season: autumn")
seas_tr   <- compute_spatiotemporal_transitions(seas_list, seas_lbls, tr_log, cyclic = TRUE)
assert_true(!is.null(seas_tr) && length(seas_tr$items) == 4L, "4 cyclic seasonal transitions computed")
assert_true(seas_tr$items[[1]]$transition_label == "Winter -> Spring", "First transition is Winter -> Spring")
assert_true(seas_tr$items[[4]]$transition_label == "Autumn -> Winter", "Fourth cyclic transition is Autumn -> Winter")

yr_tr <- compute_spatiotemporal_transitions(sl_list, c("Year 2021", "Year 2022", "Year 2023"), tr_log, cyclic = FALSE)
assert_true(!is.null(yr_tr) && length(yr_tr$items) == 2L, "2 non-cyclic yearly transitions computed")
assert_true(yr_tr$items[[1]]$transition_label == "Year 2021 -> Year 2022", "First yearly transition")

tmp_fig7_s <- tempfile(fileext = ".png")
tmp_fig7_y <- tempfile(fileext = ".png")
plot_fig7_spatiotemporal_transitions(seas_tr$items, "season", "0117", "Nitrate as N", "mg/L",
                                     tr_log, tmp_fig7_s, borrowed = FALSE)
plot_fig7_spatiotemporal_transitions(yr_tr$items, "year", "0117", "Nitrate as N", "mg/L",
                                     tr_log, tmp_fig7_y, borrowed = TRUE, rho_t = st_out$rho_temporal)
assert_true(file.exists(tmp_fig7_s) && file.info(tmp_fig7_s)$size > 1000, "Figure 7 seasonal transition grid PNG rendered")
assert_true(file.exists(tmp_fig7_y) && file.info(tmp_fig7_y)$size > 1000, "Figure 7 yearly transition grid PNG rendered")
unlink(c(tmp_fig7_s, tmp_fig7_y))

cat("  [PASS] Spatial Tobit EM censored kriging, dynamic temporal borrowing & trend mapping\n")

# ==============================================================================
# 9. Multivariate Stream Co-Kriging (Factor-LMC Joint NNGP & Co-Located)
# ==============================================================================
cat("[9/9] Testing Multivariate Stream Co-Kriging (Factor-LMC Joint NNGP & Dual CV)...\n")

source("scripts/s3n_multivariate_functions.R")

parsed_u <- parse_cokrig_targets("0117,0180")
assert_true(length(parsed_u) == 2L && parsed_u[[1]]$dataset == "wq" && parsed_u[[2]]$code == "0180",
            "parse_cokrig_targets unprefixed defaults to wq")
parsed_m <- parse_cokrig_targets("inverts:whpt_aspt,wq:0117,wq:0180,wq:9924")
assert_true(length(parsed_m) == 4L && parsed_m[[1]]$dataset == "inverts" && parsed_m[[4]]$code == "9924",
            "parse_cokrig_targets handles mixed cross-dataset targets")

# Synthetic P=3 heterotopic network (reusing df_toy-style stream tree)
set.seed(2026)
q_mv <- 60L
n_mv <- 40L
df_mv <- data.frame(
  RivID          = seq_len(q_mv),
  CatchID        = rep(901L, q_mv),
  parent_dn      = c(0L, seq_len(q_mv - 1L)),
  Shreve         = pmax(1L, as.integer(round(seq(20, 1, length.out = q_mv)))),
  Strahler       = rep(2L, q_mv),
  Dist2Mth_km    = seq(0.5, 30.0, length.out = q_mv),
  up_dist_mid_km = seq(1.0, 30.5, length.out = q_mv),
  US_Accum_km    = seq(150.0, 5.0, length.out = q_mv),
  log_AFV        = seq(0.0, -2.5, length.out = q_mv),
  Elevation      = seq(10, 220, length.out = q_mv),
  Slope          = rep(2.0, q_mv),
  BFIHOST        = rep(0.55, q_mv),
  PROPWET        = rep(0.40, q_mv),
  SAAR           = rep(900, q_mv),
  URBEXT2000     = seq(0.18, 0.01, length.out = q_mv),
  ARABLE_PCT     = seq(35, 5, length.out = q_mv),
  PASTURE_PCT    = rep(25, q_mv),
  WOODLAND_PCT   = rep(15, q_mv),
  PEAT_PCT       = rep(2, q_mv),
  Is_Aquifer     = rep(0L, q_mv),
  start_X        = rep(400000, q_mv),
  start_Y        = seq(300000, 329500, length.out = q_mv),
  mid_X          = rep(400000, q_mv),
  mid_Y          = seq(300250, 329750, length.out = q_mv),
  end_X          = rep(400000, q_mv),
  end_Y          = seq(300500, 330000, length.out = q_mv),
  stringsAsFactors = FALSE
)
obs_mv_rows <- seq(1L, 55L, length.out = n_mv)
obs_mv_rows <- unique(as.integer(round(obs_mv_rows)))
if (length(obs_mv_rows) < n_mv) obs_mv_rows <- seq_len(n_mv)

# Latent shared spatial driver along downstream/upstream gradient -> strong cross-correlation
lat_f <- sin(seq(0, 2.2 * pi, length.out = length(obs_mv_rows))) + 0.6 * cos(seq(0, pi, length.out = length(obs_mv_rows)))
Y_mv_raw <- cbind(
  log(pmax(0.2, 5.0 + 2.2 * lat_f + rnorm(length(obs_mv_rows), 0, 0.25))),
  log(pmax(0.02, 0.25 + 0.14 * lat_f + rnorm(length(obs_mv_rows), 0, 0.025))),
  5.8 - 0.85 * lat_f + rnorm(length(obs_mv_rows), 0, 0.18)
)
# Introduce heterotopic missingness (target 1 missing at 20% of sites)
miss_t1 <- seq(3L, length(obs_mv_rows), by = 5L)
Y_mv_raw[miss_t1, 1] <- NA_real_

targets_mv_meta <- list(
  list(code = "0117", dataset = "wq", slug = "wq_0117", label = "Nitrate", name = "Nitrate-N", unit = "mg/L"),
  list(code = "0180", dataset = "wq", slug = "wq_0180", label = "Orthophosphate", name = "Ortho-P", unit = "mg/L"),
  list(code = "whpt_aspt", dataset = "inverts", slug = "inverts_whpt_aspt", label = "WHPT ASPT", name = "WHPT ASPT", unit = "score")
)
trans_mv_list <- list(
  build_transformer("log"),
  build_transformer("log"),
  build_transformer("none")
)

# Direct C++ Engine Test (Modes 0, 1, 2)
X_mv_test <- cbind(1, scale(df_mv$Elevation[obs_mv_rows]))
par_mv    <- match(df_mv$parent_dn, df_mv$RivID, nomatch = 0L)
p0_mv     <- ifelse(df_mv$parent_dn == 0 | par_mv == 0L, -1L, par_mv - 1L)
tt_mv     <- .Call("C_build_tree_tin_tout", as.integer(p0_mv))
nb_mv_tu  <- .Call("C_compute_s3n_neighbors", as.integer(length(obs_mv_rows)), 5L, as.integer(p0_mv),
                   as.integer(df_mv$CatchID), as.integer(tt_mv[[1]]), as.integer(tt_mv[[2]]),
                   as.double(df_mv$Dist2Mth_km), as.double(df_mv$up_dist_mid_km), as.double(df_mv$log_AFV),
                   as.double(df_mv$mid_X), as.double(df_mv$mid_Y), 0L)
nb_mv_td  <- .Call("C_compute_s3n_neighbors", as.integer(length(obs_mv_rows)), 5L, as.integer(p0_mv),
                   as.integer(df_mv$CatchID), as.integer(tt_mv[[1]]), as.integer(tt_mv[[2]]),
                   as.double(df_mv$Dist2Mth_km), as.double(df_mv$up_dist_mid_km), as.double(df_mv$log_AFV),
                   as.double(df_mv$mid_X), as.double(df_mv$mid_Y), 1L)
Y_std_test <- scale(Y_mv_raw)
A_tu_t <- matrix(c(0.5, 0.35, -0.35, 0.35, 0.5, -0.35, -0.35, -0.35, 0.5), 3, 3)
A_td_t <- diag(0.3, 3)
A_eu_t <- matrix(0.0, 3, 3)
Psi_t  <- diag(0.15, 3)
phis_t <- c(0.1, 0.1, 0.1, 1.0)

m0_res <- .Call("C_s3n_multivariate_nngp", X_mv_test, Y_std_test, X_mv_test, matrix(0, ncol(X_mv_test), 3),
                nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X[obs_mv_rows]), as.double(df_mv$mid_Y[obs_mv_rows]), 5L,
                A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 0L)
assert_true(is.finite(m0_res$neg2loglik) && identical(dim(m0_res$Beta), c(2L, 3L)),
            "C_s3n_multivariate_nngp Mode 0 produces finite -2*logL and Beta matrix")

m1_res <- .Call("C_s3n_multivariate_nngp", X_mv_test, Y_std_test, cbind(1, scale(df_mv$Elevation)), m0_res$Beta,
                nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X), as.double(df_mv$mid_Y), 5L,
                A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 1L)
assert_true(identical(dim(m1_res$prediction), c(q_mv, 3L)) && all(is.finite(m1_res$prediction)) && all(m1_res$se > 0),
            "C_s3n_multivariate_nngp Mode 1 produces finite predictions and positive SEs across P=3 targets")

m2_res <- .Call("C_s3n_multivariate_nngp", X_mv_test, Y_std_test, X_mv_test, m0_res$Beta,
                nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X[obs_mv_rows]), as.double(df_mv$mid_Y[obs_mv_rows]), 5L,
                A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 2L)
ok_t1 <- which(is.finite(Y_std_test[, 1]))
r2_aux_direct <- 1 - sum((Y_std_test[ok_t1, 1] - m2_res$aux_pred[ok_t1, 1])^2) / sum((Y_std_test[ok_t1, 1] - mean(Y_std_test[ok_t1, 1]))^2)
assert_true(r2_aux_direct > 0.65, sprintf("Mode 2 Aux-Aided LOO CV accurately exploits cross-variable correlation (R2=%.3f)", r2_aux_direct))

# End-to-end Factor-LMC NNGP (K=3 full rank & K=2 reduced latent factors) + Colocated fallback
fit_lmc_full <- fit_multivariate_s3n_domain(
  df_reaches      = df_mv,
  Y_obs_raw_sc    = Y_mv_raw,
  obs_up_dist_obs = df_mv$up_dist_mid_km[obs_mv_rows],
  obs_reach_rows  = obs_mv_rows,
  targets_meta    = targets_mv_meta,
  trans_list      = trans_mv_list,
  domain_id       = "mv_toy",
  domain_label    = "Toy Multivariate Basin",
  m_neighbors     = 5L,
  k_factors       = 3L,
  method          = "lmc",
  thresholds      = c(11.3, 0.2, 6.0),
  threshold_dirs  = c("above", "above", "below"),
  n_omp           = 1L
)
assert_true(length(fit_lmc_full$pred_dfs) == 3L && all(c("pred_combined", "exceed_prob") %in% names(fit_lmc_full$pred_dfs[[1]])),
            "fit_multivariate_s3n_domain returns full-network predictions + exceedance tiers per target")
assert_true(all(c("Univariate S3N", "Co-Kriging (All-Holdout)", "Co-Kriging (Aux-Aided)") %in% fit_lmc_full$cv_compare_df$Model),
            "Dual CV summary contains Univariate, All-Holdout Co-Kriging, and Aux-Aided Co-Kriging")
assert_true(all(abs(fit_lmc_full$cross_cov_df$Corr_TU) <= 1.0 + 1e-6 & abs(fit_lmc_full$cross_cov_df$Corr_TD) <= 1.0 + 1e-6),
            "Extracted cross-correlations are valid within [-1, 1]")

fit_lmc_k2 <- fit_multivariate_s3n_domain(
  df_reaches      = df_mv,
  Y_obs_raw_sc    = Y_mv_raw,
  obs_up_dist_obs = df_mv$up_dist_mid_km[obs_mv_rows],
  obs_reach_rows  = obs_mv_rows,
  targets_meta    = targets_mv_meta,
  trans_list      = trans_mv_list,
  domain_id       = "mv_toy",
  domain_label    = "Toy Multivariate Basin (K=2)",
  m_neighbors     = 5L,
  k_factors       = 2L,
  method          = "lmc",
  n_omp           = 1L
)
assert_true(nrow(fit_lmc_k2$cv_compare_df) == 9L, "Factor-reduced K=2 LMC fits and evaluates all 3 targets cleanly")

fit_coloc <- fit_multivariate_s3n_domain(
  df_reaches      = df_mv,
  Y_obs_raw_sc    = Y_mv_raw,
  obs_up_dist_obs = df_mv$up_dist_mid_km[obs_mv_rows],
  obs_reach_rows  = obs_mv_rows,
  targets_meta    = targets_mv_meta,
  trans_list      = trans_mv_list,
  domain_id       = "mv_toy",
  domain_label    = "Toy Multivariate Basin (Colocated)",
  m_neighbors     = 5L,
  k_factors       = 3L,
  method          = "colocated",
  n_omp           = 1L
)
assert_true(all(is.finite(fit_coloc$pred_dfs[[1]]$pred_combined)), "Co-located residual co-kriging fallback produces finite predictions")

tmp_fig_ck <- tempfile(fileext = ".png")
plot_fig_cokrig_diagnostics(
  cross_cov_df  = fit_lmc_full$cross_cov_df,
  cv_compare_df = fit_lmc_full$cv_compare_df,
  var_names     = vapply(targets_mv_meta, function(t) t$label, character(1)),
  out_file      = tmp_fig_ck,
  network_label = "Test Network",
  source_label  = "RSN"
)
assert_true(file.exists(tmp_fig_ck) && file.info(tmp_fig_ck)$size > 5000,
            "Co-kriging cross-correlation & dual CV diagnostic PNG rendered cleanly")
unlink(tmp_fig_ck)

assert_true(!is.null(fit_lmc_full$coverage_df) && fit_lmc_full$coverage_df$Aux_Gap_Reaches[1] > 0L,
            "Heterotopic gap-filling reach counts tracked accurately in coverage_df")

tmp_rep_ck <- tempfile(fileext = ".md")
rep_lines <- generate_cokriging_report(
  cv_all       = fit_lmc_full$cv_compare_df,
  cc_all       = fit_lmc_full$cross_cov_df,
  coverage_df  = fit_lmc_full$coverage_df,
  exceed_df    = fit_lmc_full$exceed_summary,
  out_file     = tmp_rep_ck,
  method_label = "Factor-LMC Joint NNGP (K=3)",
  cv_mode      = "random"
)
assert_true(file.exists(tmp_rep_ck) && any(grepl("Executive Verdict", rep_lines)) && any(grepl("Where Co-Kriging Helps", rep_lines)),
            "Co-kriging markdown evaluation report written cleanly")
unlink(tmp_rep_ck)

cat("  [PASS] Multivariate Stream Co-Kriging (Factor-LMC Joint NNGP & Co-Located)\n")

# ==============================================================================
# 10. Figure Format (PNG / PDF / SVG), DPI & Size Scaling Options
# ==============================================================================
cat("[10/10] Testing Figure Format (png/pdf/svg), DPI, and Size Scaling...\n")

assert_true(abs(parse_fig_size_scale("150%") - 1.50) < 1e-9, "parse_fig_size_scale handles '150%'")
assert_true(abs(parse_fig_size_scale("+25%") - 1.25) < 1e-9, "parse_fig_size_scale handles '+25%'")
assert_true(abs(parse_fig_size_scale("-20%") - 0.80) < 1e-9, "parse_fig_size_scale handles '-20%'")
assert_true(abs(parse_fig_size_scale("0.75") - 0.75) < 1e-9, "parse_fig_size_scale handles '0.75'")

read_png_dims <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  hdr <- readBin(con, what = "raw", n = 24L)
  stopifnot(identical(hdr[1:4], as.raw(c(0x89, 0x50, 0x4e, 0x47))))
  w <- sum(as.integer(hdr[17:20]) * c(16777216, 65536, 256, 1))
  h <- sum(as.integer(hdr[21:24]) * c(16777216, 65536, 256, 1))
  c(width = w, height = h)
}

var_lbls_t <- vapply(targets_mv_meta, function(t) t$label, character(1))

# 10a. PNG baseline (200 DPI, 100%) vs 300 DPI vs +50% size
tmp_base <- tempfile(fileext = ".png")
configure_plot_options(args = c("--fig-format=png", "--dpi=200", "--fig-size=100%"))
plot_fig_cokrig_diagnostics(fit_lmc_full$cross_cov_df, fit_lmc_full$cv_compare_df, var_lbls_t, tmp_base)
dims_base <- read_png_dims(tmp_base)

tmp_dpi <- tempfile(fileext = ".png")
configure_plot_options(args = c("--fig-format=png", "--dpi=300", "--fig-size=100%"))
plot_fig_cokrig_diagnostics(fit_lmc_full$cross_cov_df, fit_lmc_full$cv_compare_df, var_lbls_t, tmp_dpi)
dims_dpi <- read_png_dims(tmp_dpi)
assert_true(abs(dims_dpi["width"] / dims_base["width"] - 1.5) < 0.02 &&
            abs(dims_dpi["height"] / dims_base["height"] - 1.5) < 0.02,
            sprintf("DPI=300 scales PNG pixel dimensions 1.5x (%dx%d -> %dx%d)",
                    dims_base["width"], dims_base["height"], dims_dpi["width"], dims_dpi["height"]))

tmp_sz <- tempfile(fileext = ".png")
configure_plot_options(args = c("--fig-format=png", "--dpi=200", "--fig-size=+50%"))
plot_fig_cokrig_diagnostics(fit_lmc_full$cross_cov_df, fit_lmc_full$cv_compare_df, var_lbls_t, tmp_sz)
dims_sz <- read_png_dims(tmp_sz)
assert_true(abs(dims_sz["width"] / dims_base["width"] - 1.5) < 0.02 &&
            abs(dims_sz["height"] / dims_base["height"] - 1.5) < 0.02,
            sprintf("--fig-size=+50%% scales PNG dimensions 1.5x (%dx%d -> %dx%d)",
                    dims_base["width"], dims_base["height"], dims_sz["width"], dims_sz["height"]))

# 10b. Vector PDF via --fig-format=pdf (auto-rewrites .png -> .pdf)
tmp_pdf_stem <- tempfile(fileext = ".png")
configure_plot_options(args = c("--fig-format=pdf", "--fig-size=80%"))
plot_fig_cokrig_diagnostics(fit_lmc_full$cross_cov_df, fit_lmc_full$cv_compare_df, var_lbls_t, tmp_pdf_stem)
tmp_pdf <- sub("\\.png$", ".pdf", tmp_pdf_stem)
assert_true(file.exists(tmp_pdf) && file.info(tmp_pdf)$size > 2000, "Vector PDF figure created")
pdf_hdr <- readLines(tmp_pdf, n = 1L, warn = FALSE)
assert_true(grepl("^%PDF-", pdf_hdr), "PDF file has valid %PDF- magic header")

# 10c. Vector SVG via --fig-format=svg (auto-rewrites .png -> .svg)
tmp_svg_stem <- tempfile(fileext = ".png")
configure_plot_options(args = c("--fig-format=svg", "--fig-size=120%"))
plot_fig_cokrig_diagnostics(fit_lmc_full$cross_cov_df, fit_lmc_full$cv_compare_df, var_lbls_t, tmp_svg_stem)
tmp_svg <- sub("\\.png$", ".svg", tmp_svg_stem)
assert_true(file.exists(tmp_svg) && file.info(tmp_svg)$size > 2000, "Vector SVG figure created")
svg_txt <- paste(readLines(tmp_svg, n = 10L, warn = FALSE), collapse = "\n")
assert_true(grepl("<svg", svg_txt, fixed = TRUE), "SVG file contains valid <svg root element")

# Reset global plot options & clean up temp files
options(s3n.fig_format = "png", s3n.fig_dpi = NULL, s3n.fig_scale = 1.0)
unlink(c(tmp_base, tmp_dpi, tmp_sz, tmp_pdf, tmp_svg))

cat("  [PASS] Figure formats (png/pdf/svg), DPI, and percentage size scaling\n")

# ==============================================================================
# 11. Multi-Source Preference (--prefer-source), Weights & Weighted Multivariate NNGP
# ==============================================================================
cat("[11/11] Testing Multi-Source Priority, Weights, SE, and Weighted C++ Multivariate NNGP...\n")

pref_test <- parse_prefer_source_arg(c("--prefer-source=RSN,WaterBlitz"))
assert_true(identical(pref_test, c("rsn", "waterblitz")), "parse_prefer_source_arg parses ordered lowercase sources")

sw_frac <- parse_source_kv_arg(c("--source-weight=rsn:1,waterblitz:0.25"), "--source-weight=", min_val = 1e-4, normalize_gt1 = TRUE)
assert_near(unname(sw_frac["rsn"]), 1.0, tol = 1e-8, msg = "source-weight rsn=1")
assert_near(unname(sw_frac["waterblitz"]), 0.25, tol = 1e-8, msg = "source-weight waterblitz=0.25")

sw_gt1 <- parse_source_kv_arg(c("--source-weight=rsn:5,waterblitz:1"), "--source-weight=", min_val = 1e-4, normalize_gt1 = TRUE)
assert_near(unname(sw_gt1["rsn"]), 1.0, tol = 1e-8, msg = "source-weight normalizes max > 1 to 1.0")
assert_near(unname(sw_gt1["waterblitz"]), 0.2, tol = 1e-8, msg = "source-weight normalizes 1/5 to 0.2")

df_prio <- data.frame(
  RivID            = c(10L, 10L, 10L, 20L, 20L, 30L),
  data_source      = c("RSN", "WATERBLITZ_TOBIT", "WATERBLITZ_TOBIT", "WATERBLITZ_TOBIT", "SSN", "WATERBLITZ_TOBIT"),
  val              = c(4.0, 20.0, 22.0, 8.0, 12.0, 9.0),
  y_sample         = log(c(4.0, 20.0, 22.0, 8.0, 12.0, 9.0)),
  obs_up_dist_km   = c(10, 10, 10, 20, 20, 30),
  meas_uncertainty = NA_real_,
  meas_method      = "",
  stringsAsFactors = FALSE
)

pf_res <- filter_samples_by_source_priority(df_prio, c("rsn", "ssn"))
assert_true(pf_res$n_dropped == 3L && pf_res$n_reaches == 2L,
            sprintf("filter_samples_by_source_priority drops 3 lower-priority rows across 2 reaches (got %d/%d)", pf_res$n_dropped, pf_res$n_reaches))
assert_true(all(pf_res$df$RivID == c(10L, 20L, 30L)) &&
            identical(pf_res$df$data_source, c("RSN", "SSN", "WATERBLITZ_TOBIT")),
            "Priority filter keeps RSN on RivID 10, SSN on RivID 20, and retains WaterBlitz on RivID 30")

# Repeat WaterBlitz visits (10 rows on RivID 10 if pooled) should not outvote 1 RSN visit when waterblitz weight=0.25
tr_log_t <- build_transformer("log")
agg_wt <- aggregate_samples_by_reach_weighted(
  df_prio[df_prio$RivID %in% c(10L, 30L), ],
  tr_log_t,
  src_weights = c(rsn = 1.0, waterblitz = 0.25),
  src_se      = c(rsn = 0.02, waterblitz = 0.8)
)
# RivID 10: mean(RSN)=4.0 (w=1), mean(WB)=21.0 (w=0.25) -> (1*4 + 0.25*21)/1.25 = 7.4 (vs unweighted 15.33)
assert_near(agg_wt$val[agg_wt$RivID == 10L], 7.4, tol = 1e-6, msg = "Source-grouped weighted mean prevents frequency bias")
assert_near(agg_wt$source_weight[agg_wt$RivID == 10L], 1.0, tol = 1e-6, msg = "Mixed reach retains max weight 1.0")
assert_near(agg_wt$source_weight[agg_wt$RivID == 30L], 0.25, tol = 1e-6, msg = "WaterBlitz-only reach gets weight 0.25")
assert_near(agg_wt$meas_uncertainty[agg_wt$RivID == 30L], 0.8, tol = 1e-6, msg = "Source SE applied cleanly")

# Verify C_s3n_multivariate_nngp_weighted matches unweighted when V_obs = 0 and responds cleanly when V_obs > 0
V_zero <- matrix(0.0, nrow = nrow(Y_std_test), ncol = 3L)
m0_wt0 <- .Call("C_s3n_multivariate_nngp_weighted", X_mv_test, Y_std_test, X_mv_test, matrix(0, ncol(X_mv_test), 3),
                nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X[obs_mv_rows]), as.double(df_mv$mid_Y[obs_mv_rows]), 5L,
                A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 0L, V_zero)
assert_near(m0_wt0$neg2loglik, m0_res$neg2loglik, tol = 1e-8, msg = "Zero-V weighted NNGP Mode 0 matches exact unweighted NLL")

m1_wt0 <- .Call("C_s3n_multivariate_nngp_weighted", X_mv_test, Y_std_test, cbind(1, scale(df_mv$Elevation)), m0_res$Beta,
                nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X), as.double(df_mv$mid_Y), 5L,
                A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 1L, V_zero)
assert_near(m1_wt0$prediction, m1_res$prediction, tol = 1e-8, msg = "Zero-V weighted NNGP Mode 1 matches unweighted predictions")

V_hi <- matrix(0.0, nrow = nrow(Y_std_test), ncol = 3L)
V_hi[seq(1L, nrow(Y_std_test), by = 2L), 1L] <- 5.0
m1_wthi <- .Call("C_s3n_multivariate_nngp_weighted", X_mv_test, Y_std_test, cbind(1, scale(df_mv$Elevation)), m0_res$Beta,
                 nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X), as.double(df_mv$mid_Y), 5L,
                 A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 1L, V_hi)
assert_true(all(m1_wthi$se[seq(1L, nrow(Y_std_test), by = 2L), 1L] > m1_res$se[seq(1L, nrow(Y_std_test), by = 2L), 1L] + 1e-3),
            "Down-weighted (V_obs > 0) reaches strictly exhibit larger posterior kriging SE than unweighted reaches")

m2_wthi <- .Call("C_s3n_multivariate_nngp_weighted", X_mv_test, Y_std_test, X_mv_test, m0_res$Beta,
                 nb_mv_tu, nb_mv_td, as.double(df_mv$mid_X[obs_mv_rows]), as.double(df_mv$mid_Y[obs_mv_rows]), 5L,
                 A_tu_t, A_td_t, A_eu_t, Psi_t, phis_t, 2L, V_hi)
assert_true(all(is.finite(m2_wthi$aux_pred[ok_t1, 1])) && all(m2_wthi$aux_se[ok_t1, 1] > 0),
            "Weighted NNGP Mode 2 Aux-Aided LOO CV runs cleanly with V_obs > 0")

# End-to-end fit_multivariate_s3n_domain with W_obs_mat & V_obs_phys_sc (both lmc & colocated)
W_test_mat <- matrix(1.0, nrow = nrow(Y_mv_raw), ncol = 3L)
W_test_mat[seq(1L, nrow(Y_mv_raw), by = 2L), 1L] <- 0.2
V_test_phys <- matrix(0.0, nrow = nrow(Y_mv_raw), ncol = 3L)
V_test_phys[seq(1L, nrow(Y_mv_raw), by = 2L), 2L] <- 0.25

fit_mv_wt <- fit_multivariate_s3n_domain(
  df_reaches      = df_mv,
  Y_obs_raw_sc    = Y_mv_raw,
  obs_up_dist_obs = df_mv$up_dist_mid_km[obs_mv_rows],
  obs_reach_rows  = obs_mv_rows,
  targets_meta    = targets_mv_meta,
  trans_list      = trans_mv_list,
  domain_id       = "mv_toy",
  domain_label    = "Toy Weighted Multivariate Basin",
  m_neighbors     = 5L,
  k_factors       = 3L,
  method          = "lmc",
  n_omp           = 1L,
  V_obs_phys_sc   = V_test_phys,
  W_obs_mat       = W_test_mat
)
assert_true(all(is.finite(fit_mv_wt$pred_dfs[[1]]$pred_combined)) && all(fit_mv_wt$pred_dfs[[1]]$se_model_sc > 0),
            "fit_multivariate_s3n_domain (lmc) runs end-to-end with source weights and physical SEs")

fit_coloc_wt <- fit_multivariate_s3n_domain(
  df_reaches      = df_mv,
  Y_obs_raw_sc    = Y_mv_raw,
  obs_up_dist_obs = df_mv$up_dist_mid_km[obs_mv_rows],
  obs_reach_rows  = obs_mv_rows,
  targets_meta    = targets_mv_meta,
  trans_list      = trans_mv_list,
  domain_id       = "mv_toy",
  domain_label    = "Toy Weighted Colocated Basin",
  m_neighbors     = 5L,
  k_factors       = 3L,
  method          = "colocated",
  n_omp           = 1L,
  V_obs_phys_sc   = V_test_phys,
  W_obs_mat       = W_test_mat
)
assert_true(all(is.finite(fit_coloc_wt$pred_dfs[[1]]$pred_combined)),
            "fit_multivariate_s3n_domain (colocated) runs end-to-end with source weights and physical SEs")

# Single-determinand fit_s3n_domain with obs_source_weight (paired comparison on identical reaches)
df_syn_unw <- df_syn
df_syn_unw$obs_source_weight <- 1.0
df_syn_sw <- df_syn
obs_sw_idx <- which(df_syn_sw$has_obs == 1L)
df_syn_sw$obs_source_weight <- 1.0
low_w_idx <- obs_sw_idx[seq(1L, length(obs_sw_idx), by = 2L)]
df_syn_sw$obs_source_weight[low_w_idx] <- 0.2

sink(tempfile())
dom_unw_v <- fit_s3n_domain(
  df = df_syn_unw, domain_id = "test", domain_label = "Synthetic Unweighted Basin",
  trans_obj = tr_log, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "exponential", combine_mode = "precision", afv_scheme = "length"
)
dom_sw_v <- fit_s3n_domain(
  df = df_syn_sw, domain_id = "test", domain_label = "Synthetic Source-Weighted Basin",
  trans_obj = tr_log, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "exponential", combine_mode = "precision", afv_scheme = "length"
)
dom_sw_j <- fit_s3n_domain(
  df = df_syn_sw, domain_id = "test", domain_label = "Synthetic Source-Weighted Joint",
  trans_obj = tr_log, var_specs = c("Elevation", "US_Accum_km", "Is_Aquifer"),
  m_neighbors = 5L, n_boot = 10L, n_omp = 1L,
  kernel = "exponential", combine_mode = "joint", afv_scheme = "length"
)
sink()
assert_true(!is.null(dom_sw_j) && all(is.finite(dom_sw_j$pred_df$pred_combined)) && all(dom_sw_j$pred_df$se_model_sc > 0),
            "fit_s3n_domain (joint) runs cleanly with obs_source_weight")
rivs_low <- df_syn_sw$RivID[low_w_idx]
se_low_unw <- mean(dom_unw_v$pred_df$se_model_sc[match(rivs_low, dom_unw_v$pred_df$RivID)])
se_low_wt  <- mean(dom_sw_v$pred_df$se_model_sc[match(rivs_low,  dom_sw_v$pred_df$RivID)])
assert_true(se_low_wt > se_low_unw,
            sprintf("Down-weighted (w=0.2) reaches have larger posterior SE than when unweighted (%.5f > %.5f)",
                    se_low_wt, se_low_unw))
assert_true(any(dom_sw_j$param_df$Parameter == "mean_source_penalty_var" & dom_sw_j$param_df$Estimate > 0),
            "param_df logs positive mean_source_penalty_var when source weights < 1")

cat("  [PASS] Multi-source priority filtering, reliability weights, source SEs, and weighted C++ NNGP\n")
cat("\nALL R & C++ UNIT TESTS PASSED!\n")
