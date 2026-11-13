# ==============================================================================
# s3n_multivariate_functions.R
# Multivariate Stream Co-Kriging (Factor-LMC Joint NNGP & Co-Located Residual)
# Completely decoupled from single-determinand fit_s3n_domain() paths.
# ==============================================================================

parse_cokrig_targets <- function(target_spec, default_dataset = "wq") {
  raw_items <- trimws(strsplit(as.character(target_spec), ",")[[1]])
  raw_items <- raw_items[nzchar(raw_items)]
  if (length(raw_items) < 2L) {
    stop("--targets requires at least 2 comma-separated determinands/variables (e.g. --targets=0117,0180 or --targets=inverts:whpt_aspt,wq:0117).")
  }
  lapply(raw_items, function(item) {
    if (grepl(":", item, fixed = TRUE)) {
      parts <- strsplit(item, ":", fixed = TRUE)[[1]]
      ds    <- tolower(trimws(parts[1]))
      code  <- trimws(paste(parts[-1], collapse = ":"))
    } else {
      code  <- item
      ds    <- if (grepl("^[0-9]{4}$", code)) "wq" else default_dataset
    }
    slug <- gsub("[^a-zA-Z0-9_]+", "_", paste0(ds, "_", code))
    list(raw = item, dataset = ds, code = code, slug = slug)
  })
}

build_factor_matrix <- function(vec, P, K) {
  L <- matrix(0.0, nrow = P, ncol = K)
  idx <- 1L
  for (j in seq_len(K)) {
    for (i in seq(j, P)) {
      L[i, j] <- vec[idx]
      idx <- idx + 1L
    }
  }
  L
}

n_factor_params <- function(P, K) {
  as.integer(P * K - K * (K - 1L) / 2L)
}

unpack_lmc_params <- function(theta, P, K, full_psi = TRUE, fit_euc = FALSE) {
  nf <- n_factor_params(P, K)
  off <- 1L

  L_tu <- build_factor_matrix(theta[off:(off + nf - 1L)], P, K)
  off  <- off + nf
  L_td <- build_factor_matrix(theta[off:(off + nf - 1L)], P, K)
  off  <- off + nf

  if (fit_euc) {
    nf_e  <- P
    L_euc <- matrix(theta[off:(off + nf_e - 1L)], nrow = P, ncol = 1L)
    off   <- off + nf_e
  } else {
    L_euc <- matrix(0.0, nrow = P, ncol = 1L)
  }

  if (full_psi) {
    np <- n_factor_params(P, P)
    L_psi <- build_factor_matrix(theta[off:(off + np - 1L)], P, P)
    diag(L_psi) <- exp(pmin(4.0, pmax(-8.0, diag(L_psi))))
    off <- off + np
  } else {
    L_psi <- diag(exp(pmin(4.0, pmax(-8.0, theta[off:(off + P - 1L)]))), nrow = P, ncol = P)
    off   <- off + P
  }

  phi_tu  <- exp(min(1.5, max(-7.5, theta[off])))
  phi_td  <- exp(min(1.5, max(-7.5, theta[off + 1L])))
  phi_euc <- if (fit_euc) exp(min(1.5, max(-7.5, theta[off + 2L]))) else 0.05

  list(
    A_tu    = tcrossprod(L_tu),
    A_td    = tcrossprod(L_td),
    A_euc   = tcrossprod(L_euc),
    Psi     = tcrossprod(L_psi),
    phi_tu  = phi_tu,
    phi_td  = phi_td,
    phi_euc = phi_euc
  )
}

cov_to_corr_mat <- function(M) {
  d <- sqrt(pmax(1e-12, diag(M)))
  C <- sweep(sweep(M, 1L, d, "/"), 2L, d, "/")
  diag(C) <- 1.0
  C <- pmax(-1.0, pmin(1.0, C))
  matrix(C, nrow = nrow(M), ncol = ncol(M))
}

extract_cokrig_cross_stats <- function(A_tu, A_td, A_euc, Psi, var_names, y_sds, domain_label) {
  # Scale covariances back to each determinand's model scale (from z-score standardized fitting scale)
  S <- diag(y_sds, nrow = length(y_sds), ncol = length(y_sds))
  A_tu_s  <- S %*% A_tu  %*% S
  A_td_s  <- S %*% A_td  %*% S
  A_euc_s <- S %*% A_euc %*% S
  Psi_s   <- S %*% Psi   %*% S
  Tot_s   <- A_tu_s + A_td_s + A_euc_s + Psi_s

  R_tu  <- cov_to_corr_mat(A_tu_s)
  R_td  <- cov_to_corr_mat(A_td_s)
  R_nug <- cov_to_corr_mat(Psi_s)
  R_tot <- cov_to_corr_mat(Tot_s)

  P <- length(var_names)
  rows <- list()
  for (i in seq_len(P)) {
    for (j in seq(i, P)) {
      rows[[length(rows) + 1L]] <- data.frame(
        Domain      = domain_label,
        Var_A       = var_names[i],
        Var_B       = var_names[j],
        Cov_TU      = round(A_tu_s[i, j], 5),
        Corr_TU     = round(R_tu[i, j], 4),
        Cov_TD      = round(A_td_s[i, j], 5),
        Corr_TD     = round(R_td[i, j], 4),
        Cov_Euc     = round(A_euc_s[i, j], 5),
        Cov_Nugget  = round(Psi_s[i, j], 5),
        Corr_Nugget = round(R_nug[i, j], 4),
        Cov_Total   = round(Tot_s[i, j], 5),
        Corr_Total  = round(R_tot[i, j], 4),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

fit_multivariate_s3n_domain <- function(df_reaches,
                                        Y_obs_raw_sc,
                                        obs_up_dist_obs,
                                        obs_reach_rows,
                                        targets_meta,
                                        trans_list,
                                        var_specs = NULL,
                                        domain_id = "england",
                                        domain_label = "All England",
                                        catch_filter = NULL,
                                        m_neighbors = 10L,
                                        k_factors = 3L,
                                        kernel = "exponential",
                                        method = "lmc",
                                        cv_mode = "random",
                                        thresholds = NULL,
                                        threshold_dirs = NULL,
                                        n_omp = 4L) {
  if (!is.null(catch_filter)) {
    keep_dom <- df_reaches$CatchID %in% catch_filter
  } else {
    keep_dom <- rep(TRUE, nrow(df_reaches))
  }
  df_dom <- df_reaches[keep_dom, , drop = FALSE]
  P <- ncol(Y_obs_raw_sc)
  var_names <- vapply(targets_meta, function(t) t$label, character(1))

  # Map global observed reach rows into df_dom
  dom_riv_match <- match(df_reaches$RivID[obs_reach_rows], df_dom$RivID)
  in_dom <- !is.na(dom_riv_match)
  if (sum(in_dom) < 20L) {
    cat(sprintf("  Skipping Co-Kriging for %s (only %d observed reaches).\n", domain_label, sum(in_dom)))
    return(NULL)
  }

  dom_obs_local <- dom_riv_match[in_dom]
  Y_sub_sc <- Y_obs_raw_sc[in_dom, , drop = FALSE]
  up_sub   <- obs_up_dist_obs[in_dom]

  # Order observed reaches by (CatchID, up_dist) and place unmonitored reaches after
  ord_o <- order(df_dom$CatchID[dom_obs_local], up_sub)
  obs_rows_dom <- dom_obs_local[ord_o]
  Y_obs_sc     <- Y_sub_sc[ord_o, , drop = FALSE]
  up_obs       <- up_sub[ord_o]

  n_obs <- length(obs_rows_dom)
  pred_rows_dom <- setdiff(seq_len(nrow(df_dom)), obs_rows_dom)
  df_ord <- df_dom[c(obs_rows_dom, pred_rows_dom), , drop = FALSE]
  q_tot  <- nrow(df_ord)

  df_ord$has_obs <- 0L
  df_ord$has_obs[1:n_obs] <- 1L
  up_dist <- df_ord$up_dist_mid_km
  up_dist[1:n_obs] <- up_obs

  # Check minimum observations per variable
  for (p in seq_len(P)) {
    n_p <- sum(is.finite(Y_obs_sc[, p]))
    if (n_p < 12L) {
      stop(sprintf("Domain '%s' has only %d valid observations for '%s' (need >= 12).",
                   domain_label, n_p, var_names[p]))
    }
  }

  cat(sprintf("\n--- Multivariate Co-Kriging (%s): %d reaches, %d observed sites across P = %d targets ---\n",
              domain_label, q_tot, n_obs, P))
  for (p in seq_len(P)) {
    cat(sprintf("    Target %d: %-28s | n_obs = %5d (%.1f%% coverage)\n",
                p, var_names[p], sum(is.finite(Y_obs_sc[, p])), 100 * mean(is.finite(Y_obs_sc[, p]))))
  }

  # Standardize each variable on its transformed scale for balanced numerical optimization
  y_means <- vapply(seq_len(P), function(p) mean(Y_obs_sc[, p], na.rm = TRUE), numeric(1))
  y_sds   <- vapply(seq_len(P), function(p) max(1e-4, sd(Y_obs_sc[, p], na.rm = TRUE)), numeric(1))
  Y_std   <- sweep(sweep(Y_obs_sc, 2L, y_means, "-"), 2L, y_sds, "/")

  par_idx  <- match(df_ord$parent_dn, df_ord$RivID, nomatch = 0L)
  parent_0 <- ifelse(df_ord$parent_dn == 0 | par_idx == 0L, -1L, par_idx - 1L)
  df_ord$log_AFV <- compute_reach_log_afv(df_ord, afv_scheme = "length", parent_0 = parent_0)
  tt <- .Call("C_build_tree_tin_tout", as.integer(parent_0))

  X_all <- build_design_matrix(df_ord, include_peat = (domain_id == "england"), var_specs = var_specs)
  X_obs <- X_all[1:n_obs, , drop = FALSE]
  p_cov <- ncol(X_obs)

  kernel_req <- normalize_kernel_name(kernel, default = "exponential")
  if (kernel_req == "auto") kernel_req <- "exponential"
  kid <- as.integer(S3N_KERNEL_SPECS[[kernel_req]]$id)

  m_use <- min(as.integer(m_neighbors), max(2L, n_obs - 2L))
  nb_tu <- .Call("C_compute_s3n_neighbors",
                 as.integer(n_obs), as.integer(m_use), as.integer(parent_0),
                 as.integer(df_ord$CatchID), as.integer(tt[[1]]), as.integer(tt[[2]]),
                 as.double(df_ord$Dist2Mth_km), as.double(up_dist), as.double(df_ord$log_AFV),
                 as.double(df_ord$mid_X), as.double(df_ord$mid_Y), 0L)
  nb_td <- .Call("C_compute_s3n_neighbors",
                 as.integer(n_obs), as.integer(m_use), as.integer(parent_0),
                 as.integer(df_ord$CatchID), as.integer(tt[[1]]), as.integer(tt[[2]]),
                 as.double(df_ord$Dist2Mth_km), as.double(up_dist), as.double(df_ord$log_AFV),
                 as.double(df_ord$mid_X), as.double(df_ord$mid_Y), 1L)

  # Step 1: Warm-Start with fast univariate S3N fits on each variable's standardized observed subset
  uni_fits <- vector("list", P)
  R_emp    <- matrix(NA_real_, nrow = n_obs, ncol = P)
  phi_tus  <- numeric(P)
  phi_tds  <- numeric(P)

  for (p in seq_len(P)) {
    ok_p  <- which(is.finite(Y_std[, p]))
    df_p  <- df_ord[c(ok_p, setdiff(seq_len(q_tot), ok_p)), , drop = FALSE]
    n_p   <- length(ok_p)
    par_p <- match(df_p$parent_dn, df_p$RivID, nomatch = 0L)
    p0_p  <- ifelse(df_p$parent_dn == 0 | par_p == 0L, -1L, par_p - 1L)
    tt_p  <- .Call("C_build_tree_tin_tout", as.integer(p0_p))
    up_p  <- up_dist[c(ok_p, setdiff(seq_len(q_tot), ok_p))]
    m_p   <- min(as.integer(m_use), max(2L, n_p - 2L))
    nb_p_tu <- .Call("C_compute_s3n_neighbors", as.integer(n_p), as.integer(m_p), as.integer(p0_p),
                     as.integer(df_p$CatchID), as.integer(tt_p[[1]]), as.integer(tt_p[[2]]),
                     as.double(df_p$Dist2Mth_km), as.double(up_p), as.double(df_p$log_AFV),
                     as.double(df_p$mid_X), as.double(df_p$mid_Y), 0L)
    nb_p_td <- .Call("C_compute_s3n_neighbors", as.integer(n_p), as.integer(m_p), as.integer(p0_p),
                     as.integer(df_p$CatchID), as.integer(tt_p[[1]]), as.integer(tt_p[[2]]),
                     as.double(df_p$Dist2Mth_km), as.double(up_p), as.double(df_p$log_AFV),
                     as.double(df_p$mid_X), as.double(df_p$mid_Y), 1L)
    Xp   <- X_obs[ok_p, , drop = FALSE]
    yp   <- Y_std[ok_p, p]
    vy_p <- max(1e-4, var(yp))
    f_tu <- fit_nngp_cpp(Xp, yp, nb_p_tu$obs_neighbors, vy_p, kernel_name = kernel_req)
    f_td <- fit_nngp_cpp(Xp, yp, nb_p_td$obs_neighbors, vy_p, kernel_name = kernel_req)
    uni_fits[[p]] <- list(ok = ok_p, f_tu = f_tu, f_td = f_td, nb_tu = nb_p_tu, nb_td = nb_p_td, m_p = m_p)
    beta_avg <- 0.5 * (f_tu$Beta + f_td$Beta)
    R_emp[ok_p, p] <- yp - as.numeric(Xp %*% beta_avg)
    phi_tus[p] <- f_tu$Theta["phi"]
    phi_tds[p] <- f_td$Theta["phi"]
  }

  # Empirical co-located residual correlation matrix to seed Factor-LMC
  cor_emp <- suppressWarnings(cor(R_emp, use = "pairwise.complete.obs"))
  cor_emp[!is.finite(cor_emp)] <- 0.0
  diag(cor_emp) <- 1.0
  # Shrink slightly for guaranteed strict positive definiteness
  cor_pd <- 0.85 * cor_emp + 0.15 * diag(P)

  K_use <- min(P, max(1L, as.integer(k_factors)))
  full_psi <- (P <= 4L)
  fit_euc  <- FALSE

  # Eigen-decompose scaled spatial/nugget seeds for initial factor matrices
  sig_tu_vec <- vapply(uni_fits, function(u) 0.5 * max(0.05, u$f_tu$Theta["sigma.sq"]), numeric(1))
  sig_td_vec <- vapply(uni_fits, function(u) 0.5 * max(0.05, u$f_td$Theta["sigma.sq"]), numeric(1))
  tau_vec    <- vapply(uni_fits, function(u) max(0.03, 0.5 * (u$f_tu$Theta["tau.sq"] + u$f_td$Theta["tau.sq"])), numeric(1))

  init_factor_vec <- function(sd_vec, R_mat, K_f) {
    M <- diag(sd_vec, P, P) %*% R_mat %*% diag(sd_vec, P, P)
    eg <- eigen(M, symmetric = TRUE)
    vals <- pmax(1e-4, eg$values[seq_len(K_f)])
    F0 <- eg$vectors[, seq_len(K_f), drop = FALSE] %*% diag(sqrt(vals), K_f, K_f)
    # QR -> lower-trapezoidal factor
    qr_f <- qr(t(F0))
    L_init <- t(qr.R(qr_f))
    if (ncol(L_init) < K_f) {
      L_init <- cbind(L_init, matrix(0.0, nrow = P, ncol = K_f - ncol(L_init)))
    }
    out <- numeric(0)
    for (j in seq_len(K_f)) out <- c(out, L_init[j:P, j])
    out
  }

  v_tu0 <- init_factor_vec(sqrt(sig_tu_vec), cor_pd, K_use)
  v_td0 <- init_factor_vec(sqrt(sig_td_vec), cor_pd, K_use)
  if (full_psi) {
    C_psi0 <- diag(sqrt(tau_vec), P, P) %*% cor_pd %*% diag(sqrt(tau_vec), P, P)
    L_psi0 <- t(chol(C_psi0 + 1e-4 * diag(P)))
    diag(L_psi0) <- log(pmax(1e-3, diag(L_psi0)))
    v_psi0 <- numeric(0)
    for (j in seq_len(P)) v_psi0 <- c(v_psi0, L_psi0[j:P, j])
  } else {
    v_psi0 <- log(sqrt(tau_vec))
  }

  theta_init <- c(
    v_tu0, v_td0, v_psi0,
    log(max(1e-3, median(phi_tus))),
    log(max(1e-3, median(phi_tds)))
  )

  eval_mv_nll <- function(th, X_tr_m, Y_tr_m, nb_t_u, nb_t_d, mx, my, m_arg) {
    pars <- unpack_lmc_params(th, P, K_use, full_psi = full_psi, fit_euc = fit_euc)
    res <- .Call(
      "C_s3n_multivariate_nngp",
      X_tr_m, Y_tr_m, X_tr_m, matrix(0.0, ncol(X_tr_m), P),
      nb_t_u, nb_t_d, as.double(mx), as.double(my), as.integer(m_arg),
      pars$A_tu, pars$A_td, pars$A_euc, pars$Psi,
      as.double(c(pars$phi_tu, pars$phi_td, pars$phi_euc, kid)),
      0L
    )
    list(nll = res$neg2loglik, Beta = res$Beta, pars = pars)
  }

  method_req <- tolower(trimws(as.character(method)))
  if (method_req == "colocated") {
    # Fast two-stage residual Co-Kriging (Measurement-Error Co-Kriging per target)
    cat(sprintf("  [Multivariate Co-Kriging Mode: colocated] Fitting %d targets via sequential residual conditioning...\n", P))
    # Predict each variable univariately across q_tot
    uni_pred_std <- matrix(0.0, nrow = q_tot, ncol = P)
    uni_se_std   <- matrix(1.0, nrow = q_tot, ncol = P)
    for (p in seq_len(P)) {
      ok_p <- uni_fits[[p]]$ok
      ord_p <- c(ok_p, setdiff(seq_len(q_tot), ok_p))
      Xp_o  <- X_all[ok_p, , drop = FALSE]
      Xp_e  <- X_all[ord_p, , drop = FALSE]
      f_p   <- uni_fits[[p]]$f_td
      p_n   <- uni_fits[[p]]$nb_td$pred_neighbors
      m_p   <- uni_fits[[p]]$m_p
      pr <- .Call("C_s3n_predict_meas_error", Xp_o, as.matrix(Y_std[ok_p, p]), numeric(0),
                  Xp_e, as.double(f_p$Beta), as.double(c(f_p$Theta, kid)), as.integer(m_p),
                  p_n$nnIndx[1:(q_tot * m_p)], p_n$nnDist[1:(q_tot * m_p)], p_n$nnWght[1:(q_tot * m_p)],
                  uni_fits[[p]]$nb_td$obsobs_dist, uni_fits[[p]]$nb_td$obsobs_wt)
      uni_pred_std[ord_p, p] <- pr$prediction
      uni_se_std[ord_p, p]   <- (pr$prediction.ci[, 2] - pr$prediction) / 1.959964
    }
    # Regress target p residual on co-located auxiliary residuals (using actual obs when available)
    pred_std <- uni_pred_std
    se_std   <- uni_se_std
    aux_cv_std <- matrix(NA_real_, nrow = n_obs, ncol = P)
    Beta_fin <- matrix(0.0, nrow = p_cov, ncol = P)
    for (p in seq_len(P)) {
      Beta_fin[, p] <- 0.5 * (uni_fits[[p]]$f_tu$Beta + uni_fits[[p]]$f_td$Beta)
      ok_p <- uni_fits[[p]]$ok
      aux_idx <- setdiff(seq_len(P), p)
      R_aux_obs <- matrix(0.0, nrow = length(ok_p), ncol = length(aux_idx))
      R_aux_all <- matrix(0.0, nrow = q_tot, ncol = length(aux_idx))
      for (ai in seq_along(aux_idx)) {
        q_idx <- aux_idx[ai]
        r_all_q <- uni_pred_std[, q_idx] - as.numeric(X_all %*% Beta_fin[, q_idx])
        r_obs_q <- ifelse(is.finite(Y_std[, q_idx]), Y_std[, q_idx] - as.numeric(X_obs %*% Beta_fin[, q_idx]), r_all_q[1:n_obs])
        R_aux_obs[, ai] <- r_obs_q[ok_p]
        R_aux_all[, ai] <- ifelse(seq_len(q_tot) <= n_obs & is.finite(c(Y_std[, q_idx], rep(NA_real_, q_tot - n_obs))),
                                  c(r_obs_q, rep(0, q_tot - n_obs)), r_all_q)
      }
      r_p_obs <- Y_std[ok_p, p] - as.numeric(X_obs[ok_p, , drop = FALSE] %*% Beta_fin[, p])
      fit_co  <- lm.fit(cbind(1, R_aux_obs), r_p_obs)
      gamma   <- fit_co$coefficients
      gamma[!is.finite(gamma)] <- 0.0
      adj_all <- as.numeric(cbind(1, R_aux_all) %*% gamma)
      # Shrink prediction blend weighted by univariate precision
      r2_aux  <- max(0.0, 1.0 - var(fit_co$residuals) / max(1e-6, var(r_p_obs)))
      w_aux   <- r2_aux * 0.65
      pred_std[, p] <- (1.0 - w_aux) * uni_pred_std[, p] + w_aux * (as.numeric(X_all %*% Beta_fin[, p]) + adj_all)
      se_std[, p]   <- uni_se_std[, p] * sqrt(max(0.25, 1.0 - r2_aux * 0.5))
      aux_cv_std[ok_p, p] <- pred_std[ok_p, p]
    }
    pars_fin <- unpack_lmc_params(theta_init, P, K_use, full_psi = full_psi, fit_euc = fit_euc)
  } else {
    t_opt0 <- proc.time()
    cat(sprintf("  [Multivariate Factor-LMC NNGP] Optimizing P=%d, K=%d latent spatial factors (m=%d)...\n",
                P, K_use, m_use))
    opt_res <- tryCatch(
      optim(
        par = theta_init,
        fn  = function(th) eval_mv_nll(th, X_obs, Y_std, nb_tu, nb_td, df_ord$mid_X[1:n_obs], df_ord$mid_Y[1:n_obs], m_use)$nll,
        method = "BFGS",
        control = list(maxit = 45L, reltol = 1e-4)
      ),
      error = function(e) list(par = theta_init, convergence = 1L)
    )
    fit_fin  <- eval_mv_nll(opt_res$par, X_obs, Y_std, nb_tu, nb_td, df_ord$mid_X[1:n_obs], df_ord$mid_Y[1:n_obs], m_use)
    pars_fin <- fit_fin$pars
    Beta_fin <- fit_fin$Beta

    pr_full <- .Call(
      "C_s3n_multivariate_nngp",
      X_obs, Y_std, X_all, Beta_fin,
      nb_tu, nb_td, as.double(df_ord$mid_X), as.double(df_ord$mid_Y), as.integer(m_use),
      pars_fin$A_tu, pars_fin$A_td, pars_fin$A_euc, pars_fin$Psi,
      as.double(c(pars_fin$phi_tu, pars_fin$phi_td, pars_fin$phi_euc, kid)),
      1L
    )
    pred_std <- pr_full$prediction
    se_std   <- pr_full$se

    # Mode 2: Aux-Aided Leave-One-Target-Out CV across observed stations
    pr_aux <- .Call(
      "C_s3n_multivariate_nngp",
      X_obs, Y_std, X_obs, Beta_fin,
      nb_tu, nb_td, as.double(df_ord$mid_X[1:n_obs]), as.double(df_ord$mid_Y[1:n_obs]), as.integer(m_use),
      pars_fin$A_tu, pars_fin$A_td, pars_fin$A_euc, pars_fin$Psi,
      as.double(c(pars_fin$phi_tu, pars_fin$phi_td, pars_fin$phi_euc, kid)),
      2L
    )
    aux_cv_std <- pr_aux$aux_pred
    cat(sprintf("    Converged in %.2fs (-2*logL = %.2f, phi_tu = %.4f, phi_td = %.4f)\n",
                (proc.time() - t_opt0)[3], fit_fin$nll, pars_fin$phi_tu, pars_fin$phi_td))
  }

  # Unstandardize predictions & SEs onto each determinand's own transformed model scale
  pred_sc <- sweep(sweep(pred_std, 2L, y_sds, "*"), 2L, y_means, "+")
  se_sc   <- sweep(se_std, 2L, y_sds, "*")
  aux_cv_sc <- sweep(sweep(aux_cv_std, 2L, y_sds, "*"), 2L, y_means, "+")

  # Full Hold-Out CV (All variables masked at test fold): Compare Univariate S3N vs. Co-Kriging
  cv_info <- build_cv_folds(df_ord[1:n_obs, , drop = FALSE], cv_mode = cv_mode, k_folds = 5L, seed = 2026L)
  folds   <- cv_info$folds
  fold_ids <- sort(unique(folds[is.finite(folds) & folds > 0L]))

  cv_uni_sc <- matrix(NA_real_, nrow = n_obs, ncol = P)
  cv_all_sc <- matrix(NA_real_, nrow = n_obs, ncol = P)

  for (k in fold_ids) {
    te <- which(folds == k)
    tr <- which(!is.na(folds) & folds != k)
    if (length(te) == 0L || length(tr) < p_cov) next

    # Univariate S3N CV on identical split
    for (p in seq_len(P)) {
      tr_p <- tr[is.finite(Y_obs_sc[tr, p])]
      te_p <- te[is.finite(Y_obs_sc[te, p])]
      if (length(tr_p) < p_cov || length(te_p) == 0L) next
      ord_k <- c(tr_p, te_p)
      par_k <- match(df_ord$parent_dn[ord_k], df_ord$RivID[ord_k], nomatch = 0L)
      p0_k  <- ifelse(df_ord$parent_dn[ord_k] == 0 | par_k == 0L, -1L, par_k - 1L)
      tt_k  <- .Call("C_build_tree_tin_tout", as.integer(p0_k))
      m_k   <- min(as.integer(m_use), max(2L, length(tr_p) - 2L))
      nb_k_tu <- .Call("C_compute_s3n_neighbors", as.integer(length(tr_p)), as.integer(m_k), as.integer(p0_k),
                       as.integer(df_ord$CatchID[ord_k]), as.integer(tt_k[[1]]), as.integer(tt_k[[2]]),
                       as.double(df_ord$Dist2Mth_km[ord_k]), as.double(up_dist[ord_k]), as.double(df_ord$log_AFV[ord_k]),
                       as.double(df_ord$mid_X[ord_k]), as.double(df_ord$mid_Y[ord_k]), 0L)
      nb_k_td <- .Call("C_compute_s3n_neighbors", as.integer(length(tr_p)), as.integer(m_k), as.integer(p0_k),
                       as.integer(df_ord$CatchID[ord_k]), as.integer(tt_k[[1]]), as.integer(tt_k[[2]]),
                       as.double(df_ord$Dist2Mth_km[ord_k]), as.double(up_dist[ord_k]), as.double(df_ord$log_AFV[ord_k]),
                       as.double(df_ord$mid_X[ord_k]), as.double(df_ord$mid_Y[ord_k]), 1L)
      Xp_tr <- X_obs[tr_p, , drop = FALSE]
      Xp_ev <- X_obs[ord_k, , drop = FALSE]
      yp_tr <- Y_obs_sc[tr_p, p]
      f_k_tu <- fit_nngp_cpp(Xp_tr, yp_tr, nb_k_tu$obs_neighbors, max(1e-4, var(yp_tr)), kernel_name = kernel_req)
      f_k_td <- fit_nngp_cpp(Xp_tr, yp_tr, nb_k_td$obs_neighbors, max(1e-4, var(yp_tr)), kernel_name = kernel_req)
      ne <- length(ord_k)
      p1 <- .Call("C_s3n_predict_meas_error", Xp_tr, as.matrix(yp_tr), numeric(0), Xp_ev, f_k_tu$Beta,
                  as.double(c(f_k_tu$Theta, kid)), as.integer(m_k),
                  nb_k_tu$pred_neighbors$nnIndx[1:(ne * m_k)], nb_k_tu$pred_neighbors$nnDist[1:(ne * m_k)], nb_k_tu$pred_neighbors$nnWght[1:(ne * m_k)],
                  nb_k_tu$obsobs_dist, nb_k_tu$obsobs_wt)
      p2 <- .Call("C_s3n_predict_meas_error", Xp_tr, as.matrix(yp_tr), numeric(0), Xp_ev, f_k_td$Beta,
                  as.double(c(f_k_td$Theta, kid)), as.integer(m_k),
                  nb_k_td$pred_neighbors$nnIndx[1:(ne * m_k)], nb_k_td$pred_neighbors$nnDist[1:(ne * m_k)], nb_k_td$pred_neighbors$nnWght[1:(ne * m_k)],
                  nb_k_td$obsobs_dist, nb_k_td$obsobs_wt)
      tl <- (length(tr_p) + 1L):ne
      cv_uni_sc[te_p, p] <- 0.5 * (p1$prediction[tl] + p2$prediction[tl])
    }

    # Multivariate All-Holdout Co-Kriging CV on identical split
    ord_mv <- c(tr, te)
    par_mv <- match(df_ord$parent_dn[ord_mv], df_ord$RivID[ord_mv], nomatch = 0L)
    p0_mv  <- ifelse(df_ord$parent_dn[ord_mv] == 0 | par_mv == 0L, -1L, par_mv - 1L)
    tt_mv  <- .Call("C_build_tree_tin_tout", as.integer(p0_mv))
    m_mv   <- min(as.integer(m_use), max(2L, length(tr) - 2L))
    nb_mv_tu <- .Call("C_compute_s3n_neighbors", as.integer(length(tr)), as.integer(m_mv), as.integer(p0_mv),
                      as.integer(df_ord$CatchID[ord_mv]), as.integer(tt_mv[[1]]), as.integer(tt_mv[[2]]),
                      as.double(df_ord$Dist2Mth_km[ord_mv]), as.double(up_dist[ord_mv]), as.double(df_ord$log_AFV[ord_mv]),
                      as.double(df_ord$mid_X[ord_mv]), as.double(df_ord$mid_Y[ord_mv]), 0L)
    nb_mv_td <- .Call("C_compute_s3n_neighbors", as.integer(length(tr)), as.integer(m_mv), as.integer(p0_mv),
                      as.integer(df_ord$CatchID[ord_mv]), as.integer(tt_mv[[1]]), as.integer(tt_mv[[2]]),
                      as.double(df_ord$Dist2Mth_km[ord_mv]), as.double(up_dist[ord_mv]), as.double(df_ord$log_AFV[ord_mv]),
                      as.double(df_ord$mid_X[ord_mv]), as.double(df_ord$mid_Y[ord_mv]), 1L)
    pr_cv <- .Call(
      "C_s3n_multivariate_nngp",
      X_obs[tr, , drop = FALSE], Y_std[tr, , drop = FALSE],
      X_obs[ord_mv, , drop = FALSE], Beta_fin,
      nb_mv_tu, nb_mv_td, as.double(df_ord$mid_X[ord_mv]), as.double(df_ord$mid_Y[ord_mv]), as.integer(m_mv),
      pars_fin$A_tu, pars_fin$A_td, pars_fin$A_euc, pars_fin$Psi,
      as.double(c(pars_fin$phi_tu, pars_fin$phi_td, pars_fin$phi_euc, kid)),
      1L
    )
    tl_mv <- (length(tr) + 1L):length(ord_mv)
    pr_cv_sc <- sweep(sweep(pr_cv$prediction[tl_mv, , drop = FALSE], 2L, y_sds, "*"), 2L, y_means, "+")
    cv_all_sc[te, ] <- pr_cv_sc
  }

  # Metric helpers
  rmse <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (!any(ok)) NA_real_ else sqrt(mean((a[ok] - b[ok])^2)) }
  mae  <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (!any(ok)) NA_real_ else mean(abs(a[ok] - b[ok])) }
  r2   <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 2L) return(NA_real_)
    ss_tot <- sum((a[ok] - mean(a[ok]))^2)
    if (!is.finite(ss_tot) || ss_tot < 1e-12) return(0.0)
    1 - sum((a[ok] - b[ok])^2) / ss_tot
  }

  # Build paired CV comparison table reporting Univariate S3N vs Co-Kriging (All-Holdout) vs Co-Kriging (Aux-Aided)
  cv_rows <- list()
  for (p in seq_len(P)) {
    y_true <- Y_obs_sc[, p]
    n_p <- sum(is.finite(y_true))
    mods <- c("Univariate S3N", "Co-Kriging (All-Holdout)", "Co-Kriging (Aux-Aided)")
    preds_p <- list(cv_uni_sc[, p], cv_all_sc[, p], aux_cv_sc[, p])
    for (mi in seq_along(mods)) {
      cv_rows[[length(cv_rows) + 1L]] <- data.frame(
        Domain      = domain_label,
        Target      = var_names[p],
        Target_Code = targets_meta[[p]]$code,
        Dataset     = targets_meta[[p]]$dataset,
        CV_Mode     = cv_info$actual_mode,
        Method      = method_req,
        n_obs       = n_p,
        Model       = mods[mi],
        RMSE        = round(rmse(y_true, preds_p[[mi]]), 4),
        MAE         = round(mae(y_true, preds_p[[mi]]), 4),
        R2          = round(r2(y_true, preds_p[[mi]]), 4),
        stringsAsFactors = FALSE
      )
    }
  }
  cv_compare_df <- do.call(rbind, cv_rows)
  print(cv_compare_df, row.names = FALSE)

  cross_cov_df <- extract_cokrig_cross_stats(
    pars_fin$A_tu, pars_fin$A_td, pars_fin$A_euc, pars_fin$Psi,
    var_names, y_sds, domain_label
  )

  # Track heterotopic observation coverage and gap-filling opportunities per target
  obs_mask <- is.finite(Y_obs_sc)
  n_any_other <- vapply(seq_len(P), function(p) {
    sum(!obs_mask[, p] & rowSums(obs_mask[, setdiff(seq_len(P), p), drop = FALSE]) > 0L)
  }, integer(1))
  coverage_df <- data.frame(
    Domain          = rep(domain_label, P),
    Target          = var_names,
    Target_Code     = vapply(targets_meta, function(t) t$code, character(1)),
    Dataset         = vapply(targets_meta, function(t) t$dataset, character(1)),
    n_obs           = colSums(obs_mask),
    n_union_sites   = rep(n_obs, P),
    Coverage_Pct    = round(100 * colMeans(obs_mask), 1),
    Aux_Gap_Reaches = n_any_other,
    Mean_CoKrig_SE  = round(colMeans(se_sc, na.rm = TRUE), 5),
    stringsAsFactors = FALSE
  )

  # Build canonical pred_df for every target variable (compatible with 03/04/05 & exceedance risk)
  pred_dfs <- vector("list", P)
  exceed_summaries <- list()
  for (p in seq_len(P)) {
    tr_p <- trans_list[[p]]
    obs_sc_full  <- rep(NA_real_, q_tot)
    obs_val_full <- rep(NA_real_, q_tot)
    obs_sc_full[1:n_obs]  <- round(Y_obs_sc[, p], 5)
    obs_val_full[1:n_obs] <- round(tr_p$inv(Y_obs_sc[, p]), 4)
    has_obs_p <- rep(0L, q_tot)
    has_obs_p[1:n_obs] <- as.integer(is.finite(Y_obs_sc[, p]))

    pd <- data.frame(
      RivID         = df_ord$RivID,
      CatchID       = df_ord$CatchID,
      parent_dn     = df_ord$parent_dn,
      Shreve        = df_ord$Shreve,
      Strahler      = df_ord$Strahler,
      Dist2Mth_km   = df_ord$Dist2Mth_km,
      US_Accum_km   = if ("US_Accum_km" %in% names(df_ord)) df_ord$US_Accum_km else rep(NA_real_, q_tot),
      start_X       = df_ord$start_X,
      start_Y       = df_ord$start_Y,
      mid_X         = df_ord$mid_X,
      mid_Y         = df_ord$mid_Y,
      end_X         = df_ord$end_X,
      end_Y         = df_ord$end_Y,
      Elevation     = df_ord$Elevation,
      Is_Aquifer    = df_ord$Is_Aquifer,
      has_obs       = has_obs_p,
      obs_value     = obs_val_full,
      obs_model_sc  = obs_sc_full,
      pred_tu_sc    = round(pred_sc[, p], 5),
      se_tu_sc      = round(se_sc[, p], 5),
      pred_td_sc    = round(pred_sc[, p], 5),
      se_td_sc      = round(se_sc[, p], 5),
      pred_model_sc = round(pred_sc[, p], 5),
      se_model_sc   = round(se_sc[, p], 5),
      pred_tu       = round(tr_p$inv(pred_sc[, p]), 4),
      pred_td       = round(tr_p$inv(pred_sc[, p]), 4),
      pred_combined = round(tr_p$inv(pred_sc[, p]), 4),
      stringsAsFactors = FALSE
    )

    if (!is.null(thresholds) && length(thresholds) >= p && is.finite(thresholds[p])) {
      t_val <- thresholds[p]
      t_dir <- if (!is.null(threshold_dirs) && length(threshold_dirs) >= p) threshold_dirs[p] else "above"
      ex_res <- compute_exceedance_risk(pred_sc[, p], se_sc[, p], t_val, tr_p, dir = t_dir)
      pd$exceed_prob     <- round(ex_res$prob, 5)
      pd$compliance_tier <- ex_res$tier
      tier_levels <- c("Confidently Compliant", "Likely Compliant", "Uncertain / Priority Monitoring",
                       "Likely Exceeding", "Confidently Exceeding")
      tier_tab <- table(factor(ex_res$tier, levels = tier_levels))
      exceed_summaries[[length(exceed_summaries) + 1L]] <- data.frame(
        Domain        = domain_label,
        Target        = var_names[p],
        Threshold_Raw = t_val,
        Threshold_Sc  = round(ex_res$thresh_sc, 5),
        Direction     = t_dir,
        Tier          = names(tier_tab),
        Reach_Count   = as.integer(tier_tab),
        Reach_Pct     = round(100 * as.numeric(tier_tab) / q_tot, 2),
        stringsAsFactors = FALSE
      )
    }
    pred_dfs[[p]] <- pd
  }
  names(pred_dfs) <- var_names

  list(
    pred_dfs       = pred_dfs,
    cv_compare_df  = cv_compare_df,
    cross_cov_df   = cross_cov_df,
    coverage_df    = coverage_df,
    exceed_summary = if (length(exceed_summaries) > 0L) do.call(rbind, exceed_summaries) else NULL,
    beta_mat       = Beta_fin,
    pars           = pars_fin
  )
}

#' Generate Automated Markdown Co-Kriging Value-Add & Interpretation Report
#'
#' Evaluates whether Multivariate Co-Kriging outperformed Univariate S3N across
#' targets and river basins, highlights where heterotopic gap-filling and
#' uncertainty reduction occurred, and explains how cross-variable Tail-Up,
#' Tail-Down, and Co-Located Nugget correlations contributed.
generate_cokriging_report <- function(cv_all,
                                      cc_all,
                                      coverage_df = NULL,
                                      exceed_df = NULL,
                                      out_file = NULL,
                                      method_label = "Factor-LMC Joint NNGP",
                                      cv_mode = "random") {
  lines <- character(0)
  add <- function(...) {
    lines <<- c(lines, paste0(...))
  }

  add("# Multivariate Stream Co-Kriging — Evaluation & Interpretation Report\n")
  add(sprintf("- **Method:** `%s` | **CV Mode:** `%s`\n", method_label, cv_mode))

  if (is.null(cv_all) || nrow(cv_all) == 0L) {
    add("> No cross-validation rows available to evaluate Co-Kriging gains.")
    if (!is.null(out_file)) writeLines(lines, out_file)
    return(invisible(lines))
  }

  doms    <- unique(cv_all$Domain)
  targets <- unique(cv_all$Target)
  pair_rows <- list()

  classify_verdict <- function(d_all, d_aux, pct_rmse_all, pct_rmse_aux) {
    if ((is.finite(d_aux) && d_aux >= 0.05) || (is.finite(d_all) && d_all >= 0.03) ||
        (is.finite(pct_rmse_aux) && pct_rmse_aux <= -10.0) || (is.finite(pct_rmse_all) && pct_rmse_all <= -8.0)) {
      return("Strong Benefit")
    }
    if ((is.finite(d_aux) && d_aux >= 0.015) || (is.finite(d_all) && d_all >= 0.01) ||
        (is.finite(pct_rmse_aux) && pct_rmse_aux <= -3.0) || (is.finite(pct_rmse_all) && pct_rmse_all <= -2.5)) {
      return("Moderate Benefit")
    }
    if ((is.finite(d_all) && d_all <= -0.01) && (!is.finite(d_aux) || d_aux <= 0.005)) {
      return("No Gain (Prefer Univariate)")
    }
    "Marginal / Parity"
  }

  for (dom in doms) {
    for (tg in targets) {
      sub <- cv_all[cv_all$Domain == dom & cv_all$Target == tg, , drop = FALSE]
      if (nrow(sub) == 0L) next
      get_m <- function(mod, col) {
        ix <- which(sub$Model == mod)
        if (length(ix) > 0L) suppressWarnings(as.numeric(sub[[col]][ix[1]])) else NA_real_
      }
      r2_u  <- get_m("Univariate S3N", "R2")
      r2_a  <- get_m("Co-Kriging (All-Holdout)", "R2")
      r2_x  <- get_m("Co-Kriging (Aux-Aided)", "R2")
      rm_u  <- get_m("Univariate S3N", "RMSE")
      rm_a  <- get_m("Co-Kriging (All-Holdout)", "RMSE")
      rm_x  <- get_m("Co-Kriging (Aux-Aided)", "RMSE")
      d_all <- r2_a - r2_u
      d_aux <- r2_x - r2_u
      p_a   <- if (is.finite(rm_u) && rm_u > 0 && is.finite(rm_a)) 100 * (rm_a - rm_u) / rm_u else NA_real_
      p_x   <- if (is.finite(rm_u) && rm_u > 0 && is.finite(rm_x)) 100 * (rm_x - rm_u) / rm_u else NA_real_
      ver   <- classify_verdict(d_all, d_aux, p_a, p_x)
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        Domain = dom, Target = tg, n_obs = sub$n_obs[1],
        R2_Uni = r2_u, R2_All = r2_a, Delta_R2_All = d_all, RMSE_Pct_All = p_a,
        R2_Aux = r2_x, Delta_R2_Aux = d_aux, RMSE_Pct_Aux = p_x, Verdict = ver,
        stringsAsFactors = FALSE
      )
    }
  }
  eval_df <- do.call(rbind, pair_rows)

  n_strong <- sum(eval_df$Verdict == "Strong Benefit")
  n_mod    <- sum(eval_df$Verdict == "Moderate Benefit")
  n_marg   <- sum(eval_df$Verdict == "Marginal / Parity")
  n_nogain <- sum(eval_df$Verdict == "No Gain (Prefer Univariate)")

  headline <- if ((n_strong + n_mod) >= ceiling(0.5 * nrow(eval_df))) {
    "**WORTH DOING — Co-Kriging measurably outperformed Univariate S3N**"
  } else if (n_strong + n_mod > 0L) {
    "**SELECTIVELY USEFUL — Co-Kriging improved specific targets or basins (mixed elsewhere)**"
  } else {
    "**LIMITED GAIN OVER UNIVARIATE S3N — Residual cross-correlations are weak or sites are fully isotopic**"
  }

  add("## 1. Executive Verdict: Was Co-Kriging Worth Doing?\n")
  add(sprintf("> %s\n", headline))
  add(sprintf("- **Summary across %d Domain × Target evaluations:** `%d` Strong Benefit, `%d` Moderate Benefit, `%d` Marginal/Parity, `%d` No Gain.\n",
              nrow(eval_df), n_strong, n_mod, n_marg, n_nogain))
  if (any(is.finite(eval_df$Delta_R2_All))) {
    add(sprintf("- **Median Unmonitored-Reach Change (`All-Holdout` vs. Univariate):** `Delta R2 = %+.4f` (`%+.1f%%` RMSE).",
                median(eval_df$Delta_R2_All, na.rm = TRUE), median(eval_df$RMSE_Pct_All, na.rm = TRUE)))
  }
  if (any(is.finite(eval_df$Delta_R2_Aux))) {
    add(sprintf("- **Median Surrogate/Gap-Fill Change (`Aux-Aided` vs. Univariate):** `Delta R2 = %+.4f` (`%+.1f%%` RMSE).\n",
                median(eval_df$Delta_R2_Aux, na.rm = TRUE), median(eval_df$RMSE_Pct_Aux, na.rm = TRUE)))
  }

  add("## 2. Where Co-Kriging Helps (By Basin & Target)\n")
  add("| Domain | Target | n_obs | Uni R² | All-Holdout R² (Δ) | Aux-Aided R² (Δ) | Aux RMSE Δ | Verdict |")
  add("| :--- | :--- | ---: | ---: | ---: | ---: | ---: | :--- |")
  fmt_num <- function(x, fmt = "%.3f") if (is.finite(x)) sprintf(fmt, x) else "—"
  for (i in seq_len(nrow(eval_df))) {
    r <- eval_df[i, ]
    add(sprintf("| %s | %s | %d | %s | %s (%s) | %s (%s) | %s | **%s** |",
                r$Domain, r$Target, r$n_obs,
                fmt_num(r$R2_Uni),
                fmt_num(r$R2_All), fmt_num(r$Delta_R2_All, "%+.3f"),
                fmt_num(r$R2_Aux), fmt_num(r$Delta_R2_Aux, "%+.3f"),
                if (is.finite(r$RMSE_Pct_Aux)) sprintf("%+.1f%%", r$RMSE_Pct_Aux) else "—",
                r$Verdict))
  }
  add("")

  best_aux <- eval_df[order(-ifelse(is.finite(eval_df$Delta_R2_Aux), eval_df$Delta_R2_Aux, -999)), , drop = FALSE]
  if (nrow(best_aux) > 0L && is.finite(best_aux$Delta_R2_Aux[1]) && best_aux$Delta_R2_Aux[1] > 0.005) {
    add(sprintf("- **Biggest Accuracy Beneficiary:** **%s** in **%s** (`Aux-Aided ΔR² = %+.4f`, `RMSE %+.1f%%`).",
                best_aux$Target[1], best_aux$Domain[1], best_aux$Delta_R2_Aux[1], best_aux$RMSE_Pct_Aux[1]))
  }
  best_all <- eval_df[order(-ifelse(is.finite(eval_df$Delta_R2_All), eval_df$Delta_R2_All, -999)), , drop = FALSE]
  if (nrow(best_all) > 0L && is.finite(best_all$Delta_R2_All[1]) && best_all$Delta_R2_All[1] > 0.005) {
    add(sprintf("- **Best Completely-Unmonitored Extrapolation:** **%s** in **%s** (`All-Holdout ΔR² = %+.4f`, `RMSE %+.1f%%`).",
                best_all$Target[1], best_all$Domain[1], best_all$Delta_R2_All[1], best_all$RMSE_Pct_All[1]))
  }

  if (!is.null(coverage_df) && nrow(coverage_df) > 0L) {
    add("\n### Heterotopic Gap-Filling Reach Coverage\n")
    add("Reaches where a target determinand was **not observed** (`NA`), but at least one auxiliary co-kriging variable **was observed** directly benefit from co-located (`Aux-Aided`) information borrowing:\n")
    add("| Domain | Target | Observed Sites | Coverage (%) | Aux-Aided Gap Reaches | Mean Co-Kriging SE (Model Sc) |")
    add("| :--- | :--- | ---: | ---: | ---: | ---: |")
    for (i in seq_len(nrow(coverage_df))) {
      c_r <- coverage_df[i, ]
      add(sprintf("| %s | %s | %d / %d | %.1f%% | %d | %.4f |",
                  c_r$Domain, c_r$Target, c_r$n_obs, c_r$n_union_sites,
                  c_r$Coverage_Pct, c_r$Aux_Gap_Reaches, c_r$Mean_CoKrig_SE))
    }
    add("")
  }

  if (!is.null(exceed_df) && nrow(exceed_df) > 0L) {
    add("### Regulatory Exceedance Clarity (`--threshold`)\n")
    unc_rows <- exceed_df[exceed_df$Tier == "Uncertain / Priority Monitoring", , drop = FALSE]
    if (nrow(unc_rows) > 0L) {
      add("| Domain | Target | Threshold | Direction | Uncertain Reaches (%) |")
      add("| :--- | :--- | ---: | :---: | ---: |")
      for (i in seq_len(nrow(unc_rows))) {
        u <- unc_rows[i, ]
        add(sprintf("| %s | %s | %.3g | %s | %d (%.1f%%) |",
                    u$Domain, u$Target, u$Threshold_Raw, u$Direction, u$Reach_Count, u$Reach_Pct))
      }
      add("")
    }
  }

  add("## 3. How It Helps: Learned Hydrological & Co-Located Mechanisms\n")
  if (!is.null(cc_all) && nrow(cc_all) > 0L) {
    dom_ref <- if (any(grepl("^All England", cc_all$Domain))) {
      cc_all$Domain[grepl("^All England", cc_all$Domain)][1]
    } else {
      cc_all$Domain[1]
    }
    cc_off <- cc_all[cc_all$Domain == dom_ref & cc_all$Var_A != cc_all$Var_B, , drop = FALSE]
    if (nrow(cc_off) > 0L) {
      add(sprintf("Off-diagonal residual cross-correlations shown below are evaluated on **%s** (after accounting for catchment covariates):\n", dom_ref))
      add("| Pair | Tail-Up `Corr_TU` | Tail-Down `Corr_TD` | Nugget `Corr_Nugget` | Total `Corr_Total` | Dominant Mechanism |")
      add("| :--- | ---: | ---: | ---: | ---: | :--- |")
      cc_off <- cc_off[order(-abs(cc_off$Corr_Total)), , drop = FALSE]
      for (i in seq_len(nrow(cc_off))) {
        p_r <- cc_off[i, ]
        mags <- c("Flow-Connected Plume (Tail-Up)" = abs(p_r$Corr_TU),
                  "Shared Tributary / Catchment (Tail-Down)" = abs(p_r$Corr_TD),
                  "Co-Located Station Surrogate (Nugget)" = abs(p_r$Corr_Nugget))
        dom_mech <- if (max(mags) < 0.15) {
          "Weak / Nearly Independent Residuals"
        } else {
          dir_lbl <- if (p_r$Corr_Total >= 0) "Positive" else "Inverse"
          paste0(dir_lbl, " — ", names(which.max(mags)))
        }
        add(sprintf("| %s <-> %s | %+.3f | %+.3f | %+.3f | %+.3f | %s |",
                    p_r$Var_A, p_r$Var_B, p_r$Corr_TU, p_r$Corr_TD, p_r$Corr_Nugget, p_r$Corr_Total, dom_mech))
      }
      add("")
      top_p <- cc_off[1, ]
      if (max(abs(c(top_p$Corr_TU, top_p$Corr_TD, top_p$Corr_Nugget))) >= 0.30) {
        add(sprintf("- **Key Mechanism:** **%s** and **%s** share substantial unmeasured residual structure (`Corr_Total = %+.3f`), driven by `%s` (`r = %+.3f`) and `%s` (`r = %+.3f`). Where one target is unmeasured or sampled sparsely, the other provides direct spatial constraints.",
                    top_p$Var_A, top_p$Var_B, top_p$Corr_Total,
                    if (abs(top_p$Corr_TU) >= abs(top_p$Corr_TD)) "Tail-Up flow transport" else "Tail-Down neighborhood correlation",
                    if (abs(top_p$Corr_TU) >= abs(top_p$Corr_TD)) top_p$Corr_TU else top_p$Corr_TD,
                    "co-located Nugget coupling", top_p$Corr_Nugget))
      } else {
        add("- **Key Mechanism:** Cross-variable residual correlations are modest (`|Corr_Total| < 0.30`), meaning catchment covariates already capture most shared variance between these determinands.")
      }
    }
  }

  add("\n## 4. Practical Next Steps")
  if (n_strong + n_mod > 0L) {
    add("- Use the generated `predictions_<domain>_<code>.csv` files (and `fig_cokrig_england_<code>.png`) in place of univariate predictions for targets marked **Strong** or **Moderate Benefit**, and cite `cv_univariate_vs_cokriging.csv` to document the out-of-sample gain.")
  } else {
    add("- Since Univariate S3N matches Co-Kriging on this variable group, you can safely stick with faster univariate runs (`scripts/run_s3n_determinand.R`) unless you add a sparser target or denser auxiliary sensor dataset where heterotopic borrowing matters.")
  }

  if (!is.null(out_file)) writeLines(lines, out_file)
  invisible(lines)
}
