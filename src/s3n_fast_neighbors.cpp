#include <R.h>
#include <Rinternals.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

struct Cand {
    bool fc;
    double dist;
    double wt;
    int idx; // 0-based obs index
    bool operator<(const Cand& o) const {
        if (fc != o.fc) return fc > o.fc; // FC = true comes first
        if (dist != o.dist) return dist < o.dist;
        return idx < o.idx;
    }
};

extern "C" {

// Computes DFS entry/exit times (tin and tout) for a directed river forest
// given 0-based downstream parent index parent_dn (-1 for river mouth roots).
// Reach b is upstream of (or equal to) reach a iff tin[a] <= tin[b] && tout[b] <= tout[a].
SEXP C_build_tree_tin_tout(SEXP parent_dn_r) {
    int q = Rf_length(parent_dn_r);
    const int* parent_dn = INTEGER(parent_dn_r);
    std::vector<std::vector<int>> children(q);
    std::vector<int> roots;
    for (int i = 0; i < q; ++i) {
        int p = parent_dn[i];
        if (p < 0 || p >= q) {
            roots.push_back(i);
        } else {
            children[p].push_back(i);
        }
    }
    SEXP tin_r = PROTECT(Rf_allocVector(INTSXP, q));
    SEXP tout_r = PROTECT(Rf_allocVector(INTSXP, q));
    int* tin = INTEGER(tin_r);
    int* tout = INTEGER(tout_r);

    int timer = 0;
    struct Frame { int u; int child_idx; };
    std::vector<Frame> stack;
    for (int r : roots) {
        stack.push_back({r, 0});
        tin[r] = ++timer;
        while (!stack.empty()) {
            Frame& top = stack.back();
            int u = top.u;
            if (top.child_idx < (int)children[u].size()) {
                int v = children[u][top.child_idx++];
                tin[v] = ++timer;
                stack.push_back({v, 0});
            } else {
                tout[u] = timer;
                stack.pop_back();
            }
        }
    }

    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, tin_r);
    SET_VECTOR_ELT(res, 1, tout_r);
    UNPROTECT(3);
    return res;
}

// Computes obsobs_dist (n x n), obsobs_wt (n x n), obs_neighbors, and pred_neighbors
// model_type:
//   0 = Tail-Up (flow-connected covariance with additive function weights sqrt(AFV_up/AFV_dn))
//   1 = Tail-Down (hydrologic stream distance covariance with unit weights within same river network)
SEXP C_compute_s3n_neighbors(
    SEXP n_r, SEXP m_r,
    SEXP parent_dn_r, SEXP net_id_r,
    SEXP tin_r, SEXP tout_r,
    SEXP d2m_km_r, SEXP up_dist_km_r, SEXP log_afv_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP model_type_r
) {
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int model_type = INTEGER(model_type_r)[0];
    int q = Rf_length(parent_dn_r);

    const int* parent_dn = INTEGER(parent_dn_r);
    const int* net_id = INTEGER(net_id_r);
    const int* tin = INTEGER(tin_r);
    const int* tout = INTEGER(tout_r);
    const double* d2m_km = REAL(d2m_km_r);
    const double* up_dist_km = REAL(up_dist_km_r);
    const double* log_afv = REAL(log_afv_r);
    const double* mid_x = REAL(mid_x_r);
    const double* mid_y = REAL(mid_y_r);

    auto calc_pair = [&](int a, int b, bool& fc, double& dist, double& wt) {
        if (a == b) {
            fc = true;
            dist = 0.0;
            wt = 1.0;
            return;
        }
        if (net_id[a] == net_id[b]) {
            int ta = tin[a], tb = tin[b];
            if ((ta <= tb && tout[b] <= tout[a]) || (tb <= ta && tout[a] <= tout[b])) {
                fc = true;
                dist = std::max(1e-4, std::fabs(up_dist_km[a] - up_dist_km[b]));
                if (model_type == 1) {
                    wt = std::min(1.0 - 1e-6, std::exp(-0.001 * dist));
                } else {
                    double afv_up = std::min(log_afv[a], log_afv[b]);
                    double afv_dn = std::max(log_afv[a], log_afv[b]);
                    wt = std::min(1.0 - 1e-6, std::exp(0.5 * (afv_up - afv_dn) - 0.001 * dist));
                }
            } else {
                int anc = parent_dn[a];
                int hops = 0;
                while (anc >= 0 && anc < q && hops < q && !(tin[anc] <= tb && tout[b] <= tout[anc])) {
                    anc = parent_dn[anc];
                    ++hops;
                }
                double ncj_up = (anc >= 0 && anc < q && hops < q) ? d2m_km[anc] : 0.0;
                dist = std::max(1e-4, up_dist_km[a] + up_dist_km[b] - 2.0 * ncj_up);
                if (model_type == 1) {
                    fc = true;
                    wt = std::min(1.0 - 1e-6, std::exp(-0.001 * dist));
                } else {
                    fc = false;
                    wt = 0.0;
                }
            }
        } else {
            fc = false;
            wt = 0.0;
            double dx = (mid_x[a] - mid_x[b]) * 0.001;
            double dy = (mid_y[a] - mid_y[b]) * 0.001;
            dist = 1e6 + std::sqrt(dx * dx + dy * dy);
        }
    };

    // 1. Compute symmetric n x n obsobs_dist and obsobs_wt
    SEXP obsobs_dist_r = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    SEXP obsobs_wt_r = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    double* obsobs_dist = REAL(obsobs_dist_r);
    double* obsobs_wt = REAL(obsobs_wt_r);

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < n; ++i) {
        obsobs_dist[i + i * n] = 0.0;
        obsobs_wt[i + i * n] = 1.0;
        for (int j = 0; j < i; ++j) {
            bool fc; double dist, wt;
            calc_pair(i, j, fc, dist, wt);
            obsobs_dist[i + j * n] = dist;
            obsobs_dist[j + i * n] = dist;
            obsobs_wt[i + j * n] = wt;
            obsobs_wt[j + i * n] = wt;
        }
    }

    // 2. Build obs_neighbors for estimation
    std::vector<int> numnn(n), ind1stnn(n), numDelem(n), ind1stD(n);
    int total_d = 0, total_D = 0;
    for (int i = 0; i < n; ++i) {
        numnn[i] = (i < m) ? i : m;
        ind1stnn[i] = total_d;
        total_d += numnn[i];
        numDelem[i] = numnn[i] * numnn[i];
        ind1stD[i] = total_D;
        total_D += numDelem[i];
    }

    SEXP nnIndxLU_r = PROTECT(Rf_allocVector(INTSXP, 2 * n));
    SEXP CIndx_r = PROTECT(Rf_allocVector(INTSXP, 2 * n));
    int* nnIndxLU = INTEGER(nnIndxLU_r);
    int* CIndx = INTEGER(CIndx_r);
    for (int i = 0; i < n; ++i) {
        nnIndxLU[i] = ind1stnn[i];
        nnIndxLU[n + i] = numnn[i];
        CIndx[i] = ind1stD[i];
        CIndx[n + i] = numDelem[i];
    }

    SEXP d_r = PROTECT(Rf_allocVector(REALSXP, total_d));
    SEXP nnWt_r = PROTECT(Rf_allocVector(REALSXP, total_d));
    SEXP nnIndx_r = PROTECT(Rf_allocVector(INTSXP, total_d));
    SEXP D_r = PROTECT(Rf_allocVector(REALSXP, total_D));
    SEXP Dwt_r = PROTECT(Rf_allocVector(REALSXP, total_D));

    double* d_vec = REAL(d_r);
    double* nnWt_vec = REAL(nnWt_r);
    int* nnIndx_vec = INTEGER(nnIndx_r);
    double* D_vec = REAL(D_r);
    double* Dwt_vec = REAL(Dwt_r);

    std::fill(D_vec, D_vec + total_D, 0.0);
    std::fill(Dwt_vec, Dwt_vec + total_D, 0.0);

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 1; i < n; ++i) {
        std::vector<Cand> cands(i);
        for (int j = 0; j < i; ++j) {
            double wt = obsobs_wt[i + j * n];
            double dist = obsobs_dist[i + j * n];
            cands[j] = {wt > 0.0, dist, wt, j};
        }
        int mi = numnn[i];
        if (i > mi) {
            std::partial_sort(cands.begin(), cands.begin() + mi, cands.end());
        } else {
            std::sort(cands.begin(), cands.end());
        }
        int base_d = ind1stnn[i];
        for (int k = 0; k < mi; ++k) {
            nnIndx_vec[base_d + k] = cands[k].idx;
            d_vec[base_d + k] = cands[k].dist;
            nnWt_vec[base_d + k] = cands[k].wt;
        }
        int base_D = ind1stD[i];
        for (int k = 0; k < mi; ++k) {
            int nk = cands[k].idx;
            for (int l = 0; l <= k; ++l) {
                int nl = cands[l].idx;
                int pos = base_D + l * mi + k;
                if (l == k) {
                    D_vec[pos] = 0.0;
                    Dwt_vec[pos] = 1.0;
                } else {
                    D_vec[pos] = obsobs_dist[nk + nl * n];
                    Dwt_vec[pos] = obsobs_wt[nk + nl * n];
                }
            }
        }
    }

    // 3. Build pred_neighbors for all q locations (0..n-1 are obs, n..q-1 are unobserved preds)
    // Pre-index observations by net_id: since same-network pairs have dist < 1e6 (vs >= 1e6 cross-network),
    // any net_id with >= m observations finds its exact top-m neighbors without scanning other basins.
    std::unordered_map<int, std::vector<int>> obs_by_net;
    obs_by_net.reserve(256);
    for (int j = 0; j < n; ++j) {
        obs_by_net[net_id[j]].push_back(j);
    }

    SEXP pred_nnIndx_r = PROTECT(Rf_allocVector(INTSXP, q * m));
    SEXP pred_nnDist_r = PROTECT(Rf_allocVector(REALSXP, q * m));
    SEXP pred_nnWght_r = PROTECT(Rf_allocVector(REALSXP, q * m));
    int* pred_nnIndx = INTEGER(pred_nnIndx_r);
    double* pred_nnDist = REAL(pred_nnDist_r);
    double* pred_nnWght = REAL(pred_nnWght_r);

    int m_use_pred = (n < m) ? n : m;

    #pragma omp parallel
    {
        std::vector<Cand> cands;
        cands.reserve(n);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < q; ++i) {
            auto it_net = obs_by_net.find(net_id[i]);
            bool use_same_net_only = (it_net != obs_by_net.end() && (int)it_net->second.size() >= m_use_pred && n > 0);
            cands.clear();
            if (use_same_net_only) {
                const std::vector<int>& same_obs = it_net->second;
                for (int j : same_obs) {
                    if (i < n) {
                        double wt = obsobs_wt[i + j * n];
                        double dist = obsobs_dist[i + j * n];
                        cands.push_back({wt > 0.0, dist, wt, j});
                    } else {
                        bool fc; double dist, wt;
                        calc_pair(i, j, fc, dist, wt);
                        cands.push_back({fc, dist, wt, j});
                    }
                }
            } else {
                for (int j = 0; j < n; ++j) {
                    if (i < n) {
                        double wt = obsobs_wt[i + j * n];
                        double dist = obsobs_dist[i + j * n];
                        cands.push_back({wt > 0.0, dist, wt, j});
                    } else {
                        bool fc; double dist, wt;
                        calc_pair(i, j, fc, dist, wt);
                        cands.push_back({fc, dist, wt, j});
                    }
                }
            }
            if ((int)cands.size() > m_use_pred) {
                std::partial_sort(cands.begin(), cands.begin() + m_use_pred, cands.end());
            } else {
                std::sort(cands.begin(), cands.end());
            }
            int base_p = i * m;
            for (int k = 0; k < m; ++k) {
                if (k < m_use_pred) {
                    pred_nnIndx[base_p + k] = cands[k].idx;
                    pred_nnDist[base_p + k] = cands[k].dist;
                    pred_nnWght[base_p + k] = cands[k].wt;
                } else {
                    pred_nnIndx[base_p + k] = (m_use_pred > 0) ? cands[m_use_pred - 1].idx : 0;
                    pred_nnDist[base_p + k] = 1e6;
                    pred_nnWght[base_p + k] = 0.0;
                }
            }
        }
    }

    // Package obs_neighbors list
    SEXP obs_nb = PROTECT(Rf_allocVector(VECSXP, 8));
    SEXP obs_nb_names = PROTECT(Rf_allocVector(STRSXP, 8));
    SET_VECTOR_ELT(obs_nb, 0, d_r);          SET_STRING_ELT(obs_nb_names, 0, Rf_mkChar("d"));
    SET_VECTOR_ELT(obs_nb, 1, nnWt_r);       SET_STRING_ELT(obs_nb_names, 1, Rf_mkChar("nnWt"));
    SET_VECTOR_ELT(obs_nb, 2, nnIndx_r);     SET_STRING_ELT(obs_nb_names, 2, Rf_mkChar("nnIndx"));
    SET_VECTOR_ELT(obs_nb, 3, nnIndxLU_r);   SET_STRING_ELT(obs_nb_names, 3, Rf_mkChar("nnIndxLU"));
    SET_VECTOR_ELT(obs_nb, 4, D_r);          SET_STRING_ELT(obs_nb_names, 4, Rf_mkChar("D"));
    SET_VECTOR_ELT(obs_nb, 5, Dwt_r);        SET_STRING_ELT(obs_nb_names, 5, Rf_mkChar("Dwt"));
    SET_VECTOR_ELT(obs_nb, 6, CIndx_r);      SET_STRING_ELT(obs_nb_names, 6, Rf_mkChar("CIndx"));
    SET_VECTOR_ELT(obs_nb, 7, Rf_ScalarInteger(total_D)); SET_STRING_ELT(obs_nb_names, 7, Rf_mkChar("Length.D"));
    Rf_setAttrib(obs_nb, R_NamesSymbol, obs_nb_names);

    // Package pred_neighbors list
    SEXP pred_nb = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP pred_nb_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_VECTOR_ELT(pred_nb, 0, pred_nnIndx_r); SET_STRING_ELT(pred_nb_names, 0, Rf_mkChar("nnIndx"));
    SET_VECTOR_ELT(pred_nb, 1, pred_nnDist_r); SET_STRING_ELT(pred_nb_names, 1, Rf_mkChar("nnDist"));
    SET_VECTOR_ELT(pred_nb, 2, pred_nnWght_r); SET_STRING_ELT(pred_nb_names, 2, Rf_mkChar("nnWght"));
    Rf_setAttrib(pred_nb, R_NamesSymbol, pred_nb_names);

    // Top-level return list
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_VECTOR_ELT(out, 0, obsobs_dist_r); SET_STRING_ELT(out_names, 0, Rf_mkChar("obsobs_dist"));
    SET_VECTOR_ELT(out, 1, obsobs_wt_r);   SET_STRING_ELT(out_names, 1, Rf_mkChar("obsobs_wt"));
    SET_VECTOR_ELT(out, 2, obs_nb);        SET_STRING_ELT(out_names, 2, Rf_mkChar("obs_neighbors"));
    SET_VECTOR_ELT(out, 3, pred_nb);       SET_STRING_ELT(out_names, 3, Rf_mkChar("pred_neighbors"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);

    UNPROTECT(18);
    return out;
}

// Solves C * x = b for an m x m symmetric positive definite matrix C using Cholesky.
// If C has zero or negative pivots, adds progressive jitter to diagonal.
static bool solve_chol_stack(int m, const double* C, const double* b, double* x) {
    if (m <= 0 || m > 64) return false;
    double L[64 * 64];
    double jitter = 0.0;
    for (int attempt = 0; attempt < 5; ++attempt) {
        bool ok = true;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j <= i; ++j) {
                double s = C[i * m + j];
                if (i == j) s += jitter;
                for (int k = 0; k < j; ++k) {
                    s -= L[i * m + k] * L[j * m + k];
                }
                if (i == j) {
                    if (s <= 1e-12) {
                        ok = false;
                        break;
                    }
                    L[i * m + j] = std::sqrt(s);
                } else {
                    L[i * m + j] = s / L[j * m + j];
                }
            }
            if (!ok) break;
        }
        if (ok) {
            double y[64];
            for (int i = 0; i < m; ++i) {
                double s = b[i];
                for (int k = 0; k < i; ++k) {
                    s -= L[i * m + k] * y[k];
                }
                y[i] = s / L[i * m + i];
            }
            for (int i = m - 1; i >= 0; --i) {
                double s = y[i];
                for (int k = i + 1; k < m; ++k) {
                    s -= L[k * m + i] * x[k];
                }
                x[i] = s / L[i * m + i];
            }
            return true;
        }
        jitter = (jitter == 0.0) ? 1e-6 : jitter * 10.0;
    }
    for (int i = 0; i < m; ++i) {
        double d = C[i * m + i] + 1e-4;
        x[i] = b[i] / d;
    }
    return false;
}

// =============================================================================
// Stream & Euclidean Covariance Kernel Evaluator
// kernel_id:
//   0 = exponential: exp(-phi * d)
//   1 = spherical:   (1 - 1.5*(phi*d) + 0.5*(phi*d)^3) for phi*d < 1, else 0
//   2 = mariah:      log(phi*d + 1) / (phi*d) for d > 0, 1 for d == 0
//   3 = linear:      max(0, 1 - phi * d)
// =============================================================================
inline double eval_stream_kernel(double dist, double phi, int kernel_id) {
    if (dist <= 0.0) return 1.0;
    if (dist >= 1e5 || phi <= 0.0) return 0.0;
    double u = phi * dist;
    switch (kernel_id) {
        case 1: // spherical
            if (u >= 1.0) return 0.0;
            return 1.0 - 1.5 * u + 0.5 * u * u * u;
        case 2: // mariah
            if (u < 1e-6) return 1.0 - 0.5 * u;
            return std::log(u + 1.0) / u;
        case 3: // linear-with-sill
            if (u >= 1.0) return 0.0;
            return 1.0 - u;
        case 0: // exponential (default)
        default:
            return std::exp(-u);
    }
}

// Fast NNGP Kriging Prediction with Heteroskedastic Measurement Error (Approach A)
// Solves: [sigma^2 * R + tau_micro^2 * I + diag(v_meas)] * lambda = c
// for each reach i in 0..q-1.
// Theta_r: 3-vector (sigma.sq, tau.sq, phi) or 4-vector (..., kernel_id)
SEXP C_s3n_predict_meas_error(
    SEXP X_r,            // n x p matrix of observed covariates
    SEXP y_r,            // n x 1 matrix of observed response
    SEXP v_meas_r,       // n vector of measurement error variances (>= 0)
    SEXP X0_r,           // q x p matrix of target reach covariates
    SEXP Beta_r,         // p vector of fixed effect coefficients
    SEXP Theta_r,        // 3 or 4 vector: sigma.sq, tau.sq (microscale), phi, [kernel_id]
    SEXP m_r,            // integer scalar: number of neighbors m
    SEXP pred_nnIndx_r,  // (q * m) 0-based neighbor indices
    SEXP pred_nnDist_r,  // (q * m) target-to-neighbor distances
    SEXP pred_nnWght_r,  // (q * m) target-to-neighbor weights
    SEXP obsobs_dist_r,  // n x n neighbor pairwise distances
    SEXP obsobs_wt_r     // n x n neighbor pairwise weights
) {
    int n = Rf_nrows(X_r);
    int p = Rf_ncols(X_r);
    int q = Rf_nrows(X0_r);
    int m_stride = INTEGER(m_r)[0];
    int m = m_stride;
    if (m > 64) m = 64;
    if (m > n) m = n;

    const double* X = REAL(X_r);
    const double* y = REAL(y_r);
    const double* v_meas = (Rf_length(v_meas_r) > 0) ? REAL(v_meas_r) : nullptr;
    const double* X0 = REAL(X0_r);
    const double* Beta = REAL(Beta_r);
    const double* Theta = REAL(Theta_r);
    double sigma_sq = Theta[0];
    double tau_sq = Theta[1];
    double phi = Theta[2];
    int kernel_id = (Rf_length(Theta_r) >= 4) ? (int)std::round(Theta[3]) : 0;

    const int* pred_nnIndx = INTEGER(pred_nnIndx_r);
    const double* pred_nnDist = REAL(pred_nnDist_r);
    const double* pred_nnWght = REAL(pred_nnWght_r);
    const double* obsobs_dist = REAL(obsobs_dist_r);
    const double* obsobs_wt = REAL(obsobs_wt_r);

    SEXP pred_r = PROTECT(Rf_allocVector(REALSXP, q));
    SEXP se_r = PROTECT(Rf_allocVector(REALSXP, q));
    SEXP ci_r = PROTECT(Rf_allocMatrix(REALSXP, q, 2));
    double* pred = REAL(pred_r);
    double* se = REAL(se_r);
    double* ci = REAL(ci_r);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < q; ++i) {
        // 1. Mean prediction at reach i
        double mu_0 = 0.0;
        for (int b = 0; b < p; ++b) {
            mu_0 += X0[i + b * q] * Beta[b];
        }

        // 2. Neighbor indices, distances, weights
        int base_p = i * m_stride;
        if (m <= 0 || n <= 0) {
            double s_pred = std::sqrt(std::max(1e-8, sigma_sq + tau_sq));
            pred[i] = mu_0;
            se[i] = s_pred;
            ci[i] = mu_0 - 1.959964 * s_pred;
            ci[i + q] = mu_0 + 1.959964 * s_pred;
            continue;
        }
        double C[64 * 64];
        double c[64];
        double lambda[64];
        double resid[64];

        for (int k = 0; k < m; ++k) {
            int idx_k = pred_nnIndx[base_p + k];
            if (idx_k < 0) idx_k = 0;
            if (idx_k >= n) idx_k = n - 1;

            double d_ik = pred_nnDist[base_p + k];
            double w_ik = pred_nnWght[base_p + k];
            c[k] = (w_ik > 0.0 && d_ik < 1e5)
                   ? (sigma_sq * w_ik * eval_stream_kernel(d_ik, phi, kernel_id))
                   : 0.0;

            double mu_k = 0.0;
            for (int b = 0; b < p; ++b) {
                mu_k += X[idx_k + b * n] * Beta[b];
            }
            resid[k] = y[idx_k] - mu_k;
        }

        // 3. Build m x m covariance matrix C
        for (int k = 0; k < m; ++k) {
            int idx_k = pred_nnIndx[base_p + k];
            if (idx_k < 0) idx_k = 0;
            if (idx_k >= n) idx_k = n - 1;

            for (int l = 0; l <= k; ++l) {
                int idx_l = pred_nnIndx[base_p + l];
                if (idx_l < 0) idx_l = 0;
                if (idx_l >= n) idx_l = n - 1;

                double D_kl = obsobs_dist[idx_k + idx_l * n];
                double W_kl = obsobs_wt[idx_k + idx_l * n];
                double cov_val = (W_kl > 0.0 && D_kl < 1e5)
                                 ? (sigma_sq * W_kl * eval_stream_kernel(D_kl, phi, kernel_id))
                                 : 0.0;

                if (k == l) {
                    cov_val += tau_sq + (v_meas ? v_meas[idx_k] : 0.0);
                }
                C[k * m + l] = cov_val;
                C[l * m + k] = cov_val;
            }
        }

        // 4. Solve C * lambda = c
        solve_chol_stack(m, C, c, lambda);

        // 5. Compute spatial prediction
        double delta_y = 0.0;
        double lambda_dot_c = 0.0;
        for (int k = 0; k < m; ++k) {
            delta_y += lambda[k] * resid[k];
            lambda_dot_c += lambda[k] * c[k];
        }

        double y_pred = mu_0 + delta_y;
        double v_pred = sigma_sq + tau_sq - lambda_dot_c;
        if (v_pred < 1e-8) v_pred = 1e-8;
        double s_pred = std::sqrt(v_pred);

        pred[i] = y_pred;
        se[i] = s_pred;
        ci[i] = y_pred - 1.959964 * s_pred;
        ci[i + q] = y_pred + 1.959964 * s_pred;
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_VECTOR_ELT(out, 0, pred_r); SET_STRING_ELT(out_names, 0, Rf_mkChar("prediction"));
    SET_VECTOR_ELT(out, 1, ci_r);   SET_STRING_ELT(out_names, 1, Rf_mkChar("prediction.ci"));
    SET_VECTOR_ELT(out, 2, se_r);   SET_STRING_ELT(out_names, 2, Rf_mkChar("se"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);

    UNPROTECT(5);
    return out;
}

// =============================================================================
// C_compute_log_afv:
// Fast O(N) C++ computation of log Additive Function Values (log_afv) for a
// directed stream network using arbitrary physical reach weights (e.g. reach
// length, Shreve magnitude, or hydrological yield proxy).
//
// Args:
//   parent_dn_r : integer vector (length N), 0-based downstream parent index (-1 for outlet)
//   local_w_r   : numeric vector (length N), positive local reach weight (or accumulated weight)
//   is_accum_r  : integer scalar (1 if local_w_r is already accumulated upstream mass, 0 if local)
// =============================================================================
SEXP C_compute_log_afv(SEXP parent_dn_r, SEXP local_w_r, SEXP is_accum_r) {
    int N = Rf_length(parent_dn_r);
    const int* parent_dn = INTEGER(parent_dn_r);
    const double* local_w = REAL(local_w_r);
    int is_accum = Rf_asInteger(is_accum_r);

    std::vector<int> head(N, -1);
    std::vector<int> next_child(N, -1);
    std::vector<int> roots;
    roots.reserve(1024);

    for (int i = 0; i < N; ++i) {
        int p = parent_dn[i];
        if (p >= 0 && p < N && p != i) {
            next_child[i] = head[p];
            head[p] = i;
        } else {
            roots.push_back(i);
        }
    }

    std::vector<int> order;
    order.reserve(N);
    std::vector<int> stack = roots;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        order.push_back(u);
        for (int c = head[u]; c != -1; c = next_child[c]) {
            stack.push_back(c);
        }
    }
    if ((int)order.size() < N) {
        std::vector<char> seen(N, 0);
        for (int u : order) seen[u] = 1;
        for (int i = 0; i < N; ++i) {
            if (!seen[i]) order.push_back(i);
        }
    }

    std::vector<double> mass(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double w = local_w[i];
        if (!std::isfinite(w) || w < 1e-4) w = 1e-4;
        mass[i] = w;
    }
    if (!is_accum) {
        for (int idx = (int)order.size() - 1; idx >= 0; --idx) {
            int u = order[idx];
            double sum_c = 0.0;
            for (int c = head[u]; c != -1; c = next_child[c]) {
                sum_c += mass[c];
            }
            mass[u] += sum_c;
        }
    }

    SEXP log_afv_r = PROTECT(Rf_allocVector(REALSXP, N));
    double* log_afv = REAL(log_afv_r);
    for (int i = 0; i < N; ++i) log_afv[i] = 0.0;

    for (int u : order) {
        double sum_c = 0.0;
        for (int c = head[u]; c != -1; c = next_child[c]) {
            sum_c += mass[c];
        }
        if (sum_c > 0.0) {
            for (int c = head[u]; c != -1; c = next_child[c]) {
                double prop = mass[c] / sum_c;
                if (prop < 1e-12) prop = 1e-12;
                log_afv[c] = log_afv[u] + std::log(prop);
            }
        }
    }

    UNPROTECT(1);
    return log_afv_r;
}

// =============================================================================
// C_s3n_nngp_loglik:
// Fast OpenMP NNGP profile log-likelihood and GLS solver for any stream kernel
// (0 = exponential, 1 = spherical, 2 = mariah, 3 = linear).
// =============================================================================
SEXP C_s3n_nngp_loglik(
    SEXP X_r, SEXP y_r, SEXP obs_nb_r,
    SEXP eta_r, SEXP phi_r, SEXP kernel_id_r
) {
    int n = Rf_nrows(X_r);
    int p = Rf_ncols(X_r);
    const double* X = REAL(X_r);
    const double* y = REAL(y_r);
    double eta = std::max(0.0, Rf_asReal(eta_r));
    double phi = std::max(1e-6, Rf_asReal(phi_r));
    int kernel_id = Rf_asInteger(kernel_id_r);

    const double* d        = REAL(VECTOR_ELT(obs_nb_r, 0));
    const double* nnWt     = REAL(VECTOR_ELT(obs_nb_r, 1));
    const int*    nnIndx   = INTEGER(VECTOR_ELT(obs_nb_r, 2));
    const int*    nnIndxLU = INTEGER(VECTOR_ELT(obs_nb_r, 3));
    const double* D        = REAL(VECTOR_ELT(obs_nb_r, 4));
    const double* Dwt      = REAL(VECTOR_ELT(obs_nb_r, 5));
    const int*    CIndx    = INTEGER(VECTOR_ELT(obs_nb_r, 6));

    std::vector<double> y_tilde(n, 0.0);
    std::vector<double> X_tilde((size_t)n * p, 0.0);
    std::vector<double> log_f(n, 0.0);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        int start = nnIndxLU[i];
        int k     = nnIndxLU[n + i];
        if (k > 64) k = 64;
        double diag_i = 1.0 + eta;

        if (k <= 0) {
            double f_i = std::max(1e-8, diag_i);
            log_f[i] = std::log(f_i);
            double inv_sqrt_f = 1.0 / std::sqrt(f_i);
            y_tilde[i] = y[i] * inv_sqrt_f;
            for (int j = 0; j < p; ++j) {
                X_tilde[i + j * n] = X[i + j * n] * inv_sqrt_f;
            }
            continue;
        }

        double c_vec[64];
        double C_mat[64 * 64];
        double b_vec[64];

        int valid_k = 0;
        for (int l = 0; l < k; ++l) {
            double w_il = nnWt[start + l];
            double d_il = d[start + l];
            c_vec[l] = (w_il > 0.0 && d_il < 1e5)
                       ? (w_il * eval_stream_kernel(d_il, phi, kernel_id))
                       : 0.0;
            if (c_vec[l] > 0.0) valid_k++;
        }

        if (valid_k == 0) {
            double f_i = std::max(1e-8, diag_i);
            log_f[i] = std::log(f_i);
            double inv_sqrt_f = 1.0 / std::sqrt(f_i);
            y_tilde[i] = y[i] * inv_sqrt_f;
            for (int j = 0; j < p; ++j) {
                X_tilde[i + j * n] = X[i + j * n] * inv_sqrt_f;
            }
            continue;
        }

        int c_start = CIndx[i];
        int k_full  = nnIndxLU[n + i];
        for (int l1 = 0; l1 < k; ++l1) {
            C_mat[l1 * k + l1] = 1.0 + eta;
            for (int l2 = l1 + 1; l2 < k; ++l2) {
                int idx = c_start + l1 * k_full + l2;
                double w_12 = Dwt[idx];
                double d_12 = D[idx];
                double cov = (w_12 > 0.0 && d_12 < 1e5)
                             ? (w_12 * eval_stream_kernel(d_12, phi, kernel_id))
                             : 0.0;
                C_mat[l1 * k + l2] = cov;
                C_mat[l2 * k + l1] = cov;
            }
        }

        solve_chol_stack(k, C_mat, c_vec, b_vec);

        double c_dot_b = 0.0;
        double y_hat_nb = 0.0;
        for (int l = 0; l < k; ++l) {
            int nb = nnIndx[start + l];
            c_dot_b  += c_vec[l] * b_vec[l];
            y_hat_nb += b_vec[l] * y[nb];
        }
        double f_i = std::max(1e-8, diag_i - c_dot_b);
        log_f[i] = std::log(f_i);
        double inv_sqrt_f = 1.0 / std::sqrt(f_i);

        y_tilde[i] = (y[i] - y_hat_nb) * inv_sqrt_f;
        for (int j = 0; j < p; ++j) {
            double x_hat_nb = 0.0;
            for (int l = 0; l < k; ++l) {
                int nb = nnIndx[start + l];
                x_hat_nb += b_vec[l] * X[nb + j * n];
            }
            X_tilde[i + j * n] = (X[i + j * n] - x_hat_nb) * inv_sqrt_f;
        }
    }

    std::vector<double> XtX(p * p, 0.0);
    std::vector<double> Xty(p, 0.0);
    for (int j1 = 0; j1 < p; ++j1) {
        const double* col1 = &X_tilde[j1 * n];
        double s_y = 0.0;
        for (int i = 0; i < n; ++i) s_y += col1[i] * y_tilde[i];
        Xty[j1] = s_y;
        for (int j2 = j1; j2 < p; ++j2) {
            const double* col2 = &X_tilde[j2 * n];
            double s_xx = 0.0;
            for (int i = 0; i < n; ++i) s_xx += col1[i] * col2[i];
            XtX[j1 * p + j2] = s_xx;
            XtX[j2 * p + j1] = s_xx;
        }
    }

    SEXP beta_r = PROTECT(Rf_allocVector(REALSXP, p));
    double* beta = REAL(beta_r);
    solve_chol_stack(p, XtX.data(), Xty.data(), beta);

    double rss = 0.0;
    double sum_log_f = 0.0;
    for (int i = 0; i < n; ++i) {
        double mu_i = 0.0;
        for (int j = 0; j < p; ++j) mu_i += X_tilde[i + j * n] * beta[j];
        double res = y_tilde[i] - mu_i;
        rss += res * res;
        sum_log_f += log_f[i];
    }

    double sigma_sq = std::max(1e-10, rss / (double)n);
    double tau_sq   = eta * sigma_sq;
    double neg2ll   = (double)n * std::log(sigma_sq) + sum_log_f + (double)n;

    SEXP cov_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, p));
    double* cov_beta = REAL(cov_beta_r);
    for (int j = 0; j < p; ++j) {
        std::vector<double> e_j(p, 0.0);
        std::vector<double> col_j(p, 0.0);
        e_j[j] = 1.0;
        solve_chol_stack(p, XtX.data(), e_j.data(), col_j.data());
        for (int i = 0; i < p; ++i) {
            cov_beta[i + j * p] = col_j[i] * sigma_sq;
        }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(neg2ll));
    SET_VECTOR_ELT(out, 1, beta_r);
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(sigma_sq));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(tau_sq));
    SET_VECTOR_ELT(out, 4, cov_beta_r);

    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, Rf_mkChar("neg2loglik"));
    SET_STRING_ELT(names, 1, Rf_mkChar("Beta"));
    SET_STRING_ELT(names, 2, Rf_mkChar("sigma.sq"));
    SET_STRING_ELT(names, 3, Rf_mkChar("tau.sq"));
    SET_STRING_ELT(names, 4, Rf_mkChar("cov.beta"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(4);
    return out;
}

// =============================================================================
// C_s3n_joint_nngp:
// Joint Additive Spatial Stream Network Model (Ver Hoef & Peterson 2010):
//   C_joint(i, j) = sigma^2_TU * w_TU(i,j) * K_k(d_TU; phi_TU)
//                 + sigma^2_TD * w_TD(i,j) * K_k(d_TD; phi_TD)
//                 + sigma^2_Euc * K_k(d_Euc; phi_Euc)
//                 + (tau^2 + v_meas[i]) * I(i == j)
//
// Supports two modes via mode_r:
//   mode = 0 : Joint NNGP profile log-likelihood + GLS Beta & Cov(Beta) on n obs
//   mode = 1 : Joint NNGP universal kriging prediction across all q reaches
//
// Args:
//   X_r          : n x p matrix of observed covariates
//   y_r          : n vector of observed response
//   v_meas_r     : n vector of measurement error variances (or empty)
//   X0_r         : q x p matrix of target reach covariates (used in mode 1)
//   Beta_r       : p vector of GLS fixed-effect coefficients (used in mode 1)
//   nb_tu_r      : list from C_compute_s3n_neighbors(..., model_type = 0)
//   nb_td_r      : list from C_compute_s3n_neighbors(..., model_type = 1)
//   mid_x_r      : length q vector of reach midpoint X (meters)
//   mid_y_r      : length q vector of reach midpoint Y (meters)
//   m_r          : integer scalar number of NNGP neighbors m
//   params_r     : 8-vector: (sigma_sq_tu, sigma_sq_td, sigma_sq_euc, tau_sq,
//                             phi_tu, phi_td, phi_euc, kernel_id)
//   mode_r       : integer scalar (0 = loglik/GLS fit, 1 = kriging prediction)
// =============================================================================
SEXP C_s3n_joint_nngp(
    SEXP X_r, SEXP y_r, SEXP v_meas_r,
    SEXP X0_r, SEXP Beta_r,
    SEXP nb_tu_r, SEXP nb_td_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP m_r, SEXP params_r, SEXP mode_r
) {
    int n = Rf_nrows(X_r);
    int p = Rf_ncols(X_r);
    int m_stride = Rf_asInteger(m_r);
    int m = m_stride;
    if (m > 64) m = 64;
    if (m > n) m = n;
    int mode = Rf_asInteger(mode_r);

    const double* X = REAL(X_r);
    const double* y = REAL(y_r);
    const double* v_meas = (Rf_length(v_meas_r) == n) ? REAL(v_meas_r) : nullptr;
    const double* mid_x = REAL(mid_x_r);
    const double* mid_y = REAL(mid_y_r);

    const double* params = REAL(params_r);
    double sig_tu  = std::max(0.0, params[0]);
    double sig_td  = std::max(0.0, params[1]);
    double sig_euc = std::max(0.0, params[2]);
    double tau_sq  = std::max(0.0, params[3]);
    double phi_tu  = std::max(1e-6, params[4]);
    double phi_td  = std::max(1e-6, params[5]);
    double phi_euc = std::max(1e-6, params[6]);
    int kernel_id  = (Rf_length(params_r) >= 8) ? (int)std::round(params[7]) : 0;

    double total_sill = sig_tu + sig_td + sig_euc + tau_sq;
    if (total_sill < 1e-10) total_sill = 1e-10;

    // Extract obs x obs pairwise distances & weights from nb_tu and nb_td
    const double* tu_oo_dist = REAL(VECTOR_ELT(nb_tu_r, 0));
    const double* tu_oo_wt   = REAL(VECTOR_ELT(nb_tu_r, 1));
    const double* td_oo_dist = REAL(VECTOR_ELT(nb_td_r, 0));
    const double* td_oo_wt   = REAL(VECTOR_ELT(nb_td_r, 1));

    auto eval_obs_pair_cov = [&](int idx1, int idx2) -> double {
        double cov = 0.0;
        if (sig_tu > 0.0) {
            double w = tu_oo_wt[idx1 + idx2 * n];
            double d = tu_oo_dist[idx1 + idx2 * n];
            if (w > 0.0 && d < 1e5) {
                cov += sig_tu * w * eval_stream_kernel(d, phi_tu, kernel_id);
            }
        }
        if (sig_td > 0.0) {
            double w = td_oo_wt[idx1 + idx2 * n];
            double d = td_oo_dist[idx1 + idx2 * n];
            if (w > 0.0 && d < 1e5) {
                cov += sig_td * w * eval_stream_kernel(d, phi_td, kernel_id);
            }
        }
        if (sig_euc > 0.0) {
            double dx = (mid_x[idx1] - mid_x[idx2]) * 0.001;
            double dy = (mid_y[idx1] - mid_y[idx2]) * 0.001;
            double d_euc = std::hypot(dx, dy);
            cov += sig_euc * eval_stream_kernel(d_euc, phi_euc, kernel_id);
        }
        return cov;
    };

    if (mode == 0) {
        // Mode 0: Joint NNGP log-likelihood and GLS Beta estimation on 0..n-1
        std::vector<double> y_tilde(n, 0.0);
        std::vector<double> X_tilde((size_t)n * p, 0.0);
        std::vector<double> log_f(n, 0.0);

        #pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < n; ++i) {
            double diag_i = total_sill + (v_meas ? v_meas[i] : 0.0);
            if (i == 0 || m <= 0) {
                double f_i = std::max(1e-8, diag_i);
                log_f[i] = std::log(f_i);
                double inv_sqrt_f = 1.0 / std::sqrt(f_i);
                y_tilde[i] = y[i] * inv_sqrt_f;
                for (int j = 0; j < p; ++j) {
                    X_tilde[i + j * n] = X[i + j * n] * inv_sqrt_f;
                }
                continue;
            }

            // Rank predecessor observations j < i by joint covariance with i
            int k = (i < m) ? i : m;
            std::vector<std::pair<double, int>> cands(i);
            for (int j = 0; j < i; ++j) {
                double c_ij = eval_obs_pair_cov(i, j);
                cands[j] = {-c_ij, j}; // sort descending by covariance
            }
            if (i > k) {
                std::partial_sort(cands.begin(), cands.begin() + k, cands.end());
            } else {
                std::sort(cands.begin(), cands.end());
            }

            double c_vec[64];
            double C_mat[64 * 64];
            double b_vec[64];
            int nb_idx[64];
            int valid_k = 0;

            for (int l = 0; l < k; ++l) {
                nb_idx[l] = cands[l].second;
                c_vec[l]  = -cands[l].first;
                if (c_vec[l] > 1e-12) valid_k++;
            }

            if (valid_k == 0) {
                double f_i = std::max(1e-8, diag_i);
                log_f[i] = std::log(f_i);
                double inv_sqrt_f = 1.0 / std::sqrt(f_i);
                y_tilde[i] = y[i] * inv_sqrt_f;
                for (int j = 0; j < p; ++j) {
                    X_tilde[i + j * n] = X[i + j * n] * inv_sqrt_f;
                }
                continue;
            }

            for (int l1 = 0; l1 < k; ++l1) {
                int j1 = nb_idx[l1];
                C_mat[l1 * k + l1] = total_sill + (v_meas ? v_meas[j1] : 0.0);
                for (int l2 = l1 + 1; l2 < k; ++l2) {
                    int j2 = nb_idx[l2];
                    double cov = eval_obs_pair_cov(j1, j2);
                    C_mat[l1 * k + l2] = cov;
                    C_mat[l2 * k + l1] = cov;
                }
            }

            solve_chol_stack(k, C_mat, c_vec, b_vec);

            double c_dot_b = 0.0;
            double y_hat_nb = 0.0;
            for (int l = 0; l < k; ++l) {
                int nb = nb_idx[l];
                c_dot_b  += c_vec[l] * b_vec[l];
                y_hat_nb += b_vec[l] * y[nb];
            }
            double f_i = std::max(1e-8, diag_i - c_dot_b);
            log_f[i] = std::log(f_i);
            double inv_sqrt_f = 1.0 / std::sqrt(f_i);

            y_tilde[i] = (y[i] - y_hat_nb) * inv_sqrt_f;
            for (int j = 0; j < p; ++j) {
                double x_hat_nb = 0.0;
                for (int l = 0; l < k; ++l) {
                    int nb = nb_idx[l];
                    x_hat_nb += b_vec[l] * X[nb + j * n];
                }
                X_tilde[i + j * n] = (X[i + j * n] - x_hat_nb) * inv_sqrt_f;
            }
        }

        std::vector<double> XtX(p * p, 0.0);
        std::vector<double> Xty(p, 0.0);
        for (int j1 = 0; j1 < p; ++j1) {
            const double* col1 = &X_tilde[j1 * n];
            double s_y = 0.0;
            for (int i = 0; i < n; ++i) s_y += col1[i] * y_tilde[i];
            Xty[j1] = s_y;
            for (int j2 = j1; j2 < p; ++j2) {
                const double* col2 = &X_tilde[j2 * n];
                double s_xx = 0.0;
                for (int i = 0; i < n; ++i) s_xx += col1[i] * col2[i];
                XtX[j1 * p + j2] = s_xx;
                XtX[j2 * p + j1] = s_xx;
            }
        }

        SEXP beta_r = PROTECT(Rf_allocVector(REALSXP, p));
        double* beta = REAL(beta_r);
        solve_chol_stack(p, XtX.data(), Xty.data(), beta);

        double rss = 0.0;
        double sum_log_f = 0.0;
        for (int i = 0; i < n; ++i) {
            double mu_i = 0.0;
            for (int j = 0; j < p; ++j) mu_i += X_tilde[i + j * n] * beta[j];
            double res = y_tilde[i] - mu_i;
            rss += res * res;
            sum_log_f += log_f[i];
        }

        // Exact NNGP -2*loglik for fixed variance components (or profiled scale factor)
        double scale = std::max(1e-10, rss / (double)n);
        double neg2ll_profiled = (double)n * std::log(scale) + sum_log_f + (double)n;
        double neg2ll_exact    = rss + sum_log_f;

        SEXP cov_beta_r = PROTECT(Rf_allocMatrix(REALSXP, p, p));
        double* cov_beta = REAL(cov_beta_r);
        for (int j = 0; j < p; ++j) {
            std::vector<double> e_j(p, 0.0);
            std::vector<double> col_j(p, 0.0);
            e_j[j] = 1.0;
            solve_chol_stack(p, XtX.data(), e_j.data(), col_j.data());
            for (int i = 0; i < p; ++i) {
                cov_beta[i + j * p] = col_j[i] * scale;
            }
        }

        SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
        SET_VECTOR_ELT(out, 0, Rf_ScalarReal(neg2ll_profiled));
        SET_VECTOR_ELT(out, 1, Rf_ScalarReal(neg2ll_exact));
        SET_VECTOR_ELT(out, 2, Rf_ScalarReal(scale));
        SET_VECTOR_ELT(out, 3, beta_r);
        SET_VECTOR_ELT(out, 4, cov_beta_r);

        SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
        SET_STRING_ELT(names, 0, Rf_mkChar("neg2loglik"));
        SET_STRING_ELT(names, 1, Rf_mkChar("neg2loglik_exact"));
        SET_STRING_ELT(names, 2, Rf_mkChar("scale"));
        SET_STRING_ELT(names, 3, Rf_mkChar("Beta"));
        SET_STRING_ELT(names, 4, Rf_mkChar("cov.beta"));
        Rf_setAttrib(out, R_NamesSymbol, names);

        UNPROTECT(4);
        return out;
    }

    // Mode 1: Joint NNGP Kriging Prediction across all q reaches
    int q = Rf_nrows(X0_r);
    const double* X0 = REAL(X0_r);
    const double* Beta = REAL(Beta_r);

    SEXP pred_nb_tu = VECTOR_ELT(nb_tu_r, 3);
    SEXP pred_nb_td = VECTOR_ELT(nb_td_r, 3);
    const int*    tu_p_idx  = INTEGER(VECTOR_ELT(pred_nb_tu, 0));
    const double* tu_p_dist = REAL(VECTOR_ELT(pred_nb_tu, 1));
    const double* tu_p_wt   = REAL(VECTOR_ELT(pred_nb_tu, 2));
    const int*    td_p_idx  = INTEGER(VECTOR_ELT(pred_nb_td, 0));
    const double* td_p_dist = REAL(VECTOR_ELT(pred_nb_td, 1));
    const double* td_p_wt   = REAL(VECTOR_ELT(pred_nb_td, 2));

    // Precompute observed residuals r_j = y_j - X_j * Beta
    std::vector<double> obs_resid(n, 0.0);
    for (int j = 0; j < n; ++j) {
        double mu_j = 0.0;
        for (int b = 0; b < p; ++b) mu_j += X[j + b * n] * Beta[b];
        obs_resid[j] = y[j] - mu_j;
    }

    SEXP pred_r = PROTECT(Rf_allocVector(REALSXP, q));
    SEXP se_r   = PROTECT(Rf_allocVector(REALSXP, q));
    SEXP ci_r   = PROTECT(Rf_allocMatrix(REALSXP, q, 2));
    double* pred = REAL(pred_r);
    double* se   = REAL(se_r);
    double* ci   = REAL(ci_r);

    struct JointCand {
        int idx;
        double c_joint;
    };

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < q; ++i) {
        double mu_0 = 0.0;
        for (int b = 0; b < p; ++b) {
            mu_0 += X0[i + b * q] * Beta[b];
        }

        if (m <= 0 || n <= 0) {
            double s_pred = std::sqrt(std::max(1e-8, total_sill));
            pred[i] = mu_0;
            se[i]   = s_pred;
            ci[i]   = mu_0 - 1.959964 * s_pred;
            ci[i + q] = mu_0 + 1.959964 * s_pred;
            continue;
        }

        int base_p = i * m_stride;
        // Collect unique candidate observations from TU neighbors, TD neighbors,
        // and (when sig_euc > 0) top Euclidean neighbors
        int cand_ids[192];
        double cand_tu_d[192], cand_tu_w[192];
        double cand_td_d[192], cand_td_w[192];
        int n_cand = 0;

        auto find_or_add = [&](int obs_j) -> int {
            if (obs_j < 0 || obs_j >= n) return -1;
            for (int c = 0; c < n_cand; ++c) {
                if (cand_ids[c] == obs_j) return c;
            }
            if (n_cand >= 192) return -1;
            int c = n_cand++;
            cand_ids[c]  = obs_j;
            cand_tu_d[c] = 1e6; cand_tu_w[c] = 0.0;
            cand_td_d[c] = 1e6; cand_td_w[c] = 0.0;
            return c;
        };

        for (int k = 0; k < m; ++k) {
            int j_tu = tu_p_idx[base_p + k];
            int c = find_or_add(j_tu);
            if (c >= 0) {
                cand_tu_d[c] = tu_p_dist[base_p + k];
                cand_tu_w[c] = tu_p_wt[base_p + k];
                if (cand_tu_w[c] > 0.0 && cand_tu_d[c] < 1e5) {
                    cand_td_d[c] = cand_tu_d[c];
                    cand_td_w[c] = 1.0;
                }
            }
        }
        for (int k = 0; k < m; ++k) {
            int j_td = td_p_idx[base_p + k];
            int c = find_or_add(j_td);
            if (c >= 0) {
                double d_td = td_p_dist[base_p + k];
                double w_td = td_p_wt[base_p + k];
                if (w_td > 0.0 && d_td < 1e5) {
                    cand_td_d[c] = d_td;
                    cand_td_w[c] = w_td;
                }
            }
        }

        if (sig_euc > 0.0) {
            // Find up to m closest Euclidean observations to reach i
            double best_d2[64];
            int best_j[64];
            int n_best = 0;
            double xi = mid_x[i], yi = mid_y[i];
            for (int j = 0; j < n; ++j) {
                double dx = (xi - mid_x[j]) * 0.001;
                double dy = (yi - mid_y[j]) * 0.001;
                double d2 = dx * dx + dy * dy;
                if (n_best < m) {
                    int pos = n_best++;
                    while (pos > 0 && best_d2[pos - 1] > d2) {
                        best_d2[pos] = best_d2[pos - 1];
                        best_j[pos]  = best_j[pos - 1];
                        --pos;
                    }
                    best_d2[pos] = d2;
                    best_j[pos]  = j;
                } else if (d2 < best_d2[m - 1]) {
                    int pos = m - 1;
                    while (pos > 0 && best_d2[pos - 1] > d2) {
                        best_d2[pos] = best_d2[pos - 1];
                        best_j[pos]  = best_j[pos - 1];
                        --pos;
                    }
                    best_d2[pos] = d2;
                    best_j[pos]  = j;
                }
            }
            for (int k = 0; k < n_best; ++k) {
                find_or_add(best_j[k]);
            }
        }

        // Compute joint target-to-candidate covariance for each unique candidate
        JointCand ranked[192];
        double xi = mid_x[i], yi = mid_y[i];
        for (int c = 0; c < n_cand; ++c) {
            int j = cand_ids[c];
            double cov = 0.0;
            if (sig_tu > 0.0 && cand_tu_w[c] > 0.0 && cand_tu_d[c] < 1e5) {
                cov += sig_tu * cand_tu_w[c] * eval_stream_kernel(cand_tu_d[c], phi_tu, kernel_id);
            }
            if (sig_td > 0.0 && cand_td_w[c] > 0.0 && cand_td_d[c] < 1e5) {
                cov += sig_td * cand_td_w[c] * eval_stream_kernel(cand_td_d[c], phi_td, kernel_id);
            }
            if (sig_euc > 0.0) {
                double dx = (xi - mid_x[j]) * 0.001;
                double dy = (yi - mid_y[j]) * 0.001;
                double d_euc = std::hypot(dx, dy);
                cov += sig_euc * eval_stream_kernel(d_euc, phi_euc, kernel_id);
            }
            ranked[c] = {j, cov};
        }

        std::sort(ranked, ranked + n_cand, [](const JointCand& a, const JointCand& b) {
            if (a.c_joint != b.c_joint) return a.c_joint > b.c_joint;
            return a.idx < b.idx;
        });

        int k_use = (n_cand < m) ? n_cand : m;
        double C[64 * 64];
        double c_vec[64];
        double lambda[64];

        for (int k = 0; k < k_use; ++k) {
            c_vec[k] = ranked[k].c_joint;
            int idx_k = ranked[k].idx;
            for (int l = 0; l <= k; ++l) {
                int idx_l = ranked[l].idx;
                double cov_val = eval_obs_pair_cov(idx_k, idx_l);
                if (k == l) {
                    cov_val = total_sill + (v_meas ? v_meas[idx_k] : 0.0);
                }
                C[k * k_use + l] = cov_val;
                C[l * k_use + k] = cov_val;
            }
        }

        solve_chol_stack(k_use, C, c_vec, lambda);

        double delta_y = 0.0;
        double lambda_dot_c = 0.0;
        for (int k = 0; k < k_use; ++k) {
            int idx_k = ranked[k].idx;
            delta_y      += lambda[k] * obs_resid[idx_k];
            lambda_dot_c += lambda[k] * c_vec[k];
        }

        double y_pred = mu_0 + delta_y;
        double v_pred = total_sill - lambda_dot_c;
        if (v_pred < 1e-8) v_pred = 1e-8;
        double s_pred = std::sqrt(v_pred);

        pred[i]   = y_pred;
        se[i]     = s_pred;
        ci[i]     = y_pred - 1.959964 * s_pred;
        ci[i + q] = y_pred + 1.959964 * s_pred;
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_VECTOR_ELT(out, 0, pred_r); SET_STRING_ELT(out_names, 0, Rf_mkChar("prediction"));
    SET_VECTOR_ELT(out, 1, ci_r);   SET_STRING_ELT(out_names, 1, Rf_mkChar("prediction.ci"));
    SET_VECTOR_ELT(out, 2, se_r);   SET_STRING_ELT(out_names, 2, Rf_mkChar("se"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);

    UNPROTECT(5);
    return out;
}

}



