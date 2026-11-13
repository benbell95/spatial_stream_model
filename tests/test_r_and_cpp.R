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
unlink(c(tmp_ex_eng, tmp_ex_bas))

cat("  [PASS] Hydrological CV modes, regulatory exceedance tiers, plots, and C++ optimal monitoring design\n")
cat("\nALL R & C++ UNIT TESTS PASSED!\n")
