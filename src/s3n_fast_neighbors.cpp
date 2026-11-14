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
    std::fill(tin, tin + q, 0);
    std::fill(tout, tout + q, 0);

    int timer = 0;
    struct Frame { int u; int child_idx; };
    std::vector<Frame> stack;
    auto traverse_root = [&](int r) {
        stack.clear();
        stack.push_back({r, 0});
        tin[r] = ++timer;
        while (!stack.empty()) {
            Frame& top = stack.back();
            int u = top.u;
            if (top.child_idx < (int)children[u].size()) {
                int v = children[u][top.child_idx++];
                if (tin[v] == 0) {
                    tin[v] = ++timer;
                    stack.push_back({v, 0});
                }
            } else {
                tout[u] = timer;
                stack.pop_back();
            }
        }
    };
    for (int r : roots) {
        if (tin[r] == 0) traverse_root(r);
    }
    for (int i = 0; i < q; ++i) {
        if (tin[i] == 0) traverse_root(i);
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
    int q_pred = (Rf_length(n_r) >= 2) ? std::min(q, std::max(0, INTEGER(n_r)[1])) : q;

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

    SEXP pred_nnIndx_r = PROTECT(Rf_allocVector(INTSXP, q_pred * m));
    SEXP pred_nnDist_r = PROTECT(Rf_allocVector(REALSXP, q_pred * m));
    SEXP pred_nnWght_r = PROTECT(Rf_allocVector(REALSXP, q_pred * m));
    int* pred_nnIndx = INTEGER(pred_nnIndx_r);
    double* pred_nnDist = REAL(pred_nnDist_r);
    double* pred_nnWght = REAL(pred_nnWght_r);

    int m_use_pred = (n < m) ? n : m;

    if (q_pred > 0) {
        #pragma omp parallel
        {
            std::vector<Cand> cands;
            cands.reserve(n);

            #pragma omp for schedule(dynamic, 64)
            for (int i = 0; i < q_pred; ++i) {
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
    if (m <= 0) return false;
    std::fill(x, x + m, 0.0);
    double L_stack[64 * 64];
    std::vector<double> L_heap;
    double* L = L_stack;
    if (m > 64) {
        L_heap.resize((size_t)m * m);
        L = L_heap.data();
    }
    double jitter = 0.0;
    for (int attempt = 0; attempt < 5; ++attempt) {
        bool ok = true;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j <= i; ++j) {
                double s = C[(size_t)i * m + j];
                if (i == j) s += jitter;
                for (int k = 0; k < j; ++k) {
                    s -= L[(size_t)i * m + k] * L[(size_t)j * m + k];
                }
                if (i == j) {
                    if (s <= 1e-12) {
                        ok = false;
                        break;
                    }
                    L[(size_t)i * m + j] = std::sqrt(s);
                } else {
                    L[(size_t)i * m + j] = s / L[(size_t)j * m + j];
                }
            }
            if (!ok) break;
        }
        if (ok) {
            std::vector<double> y_buf;
            double y_stack[64];
            double* y = (m <= 64) ? y_stack : (y_buf.resize(m), y_buf.data());
            for (int i = 0; i < m; ++i) {
                double s = b[i];
                for (int k = 0; k < i; ++k) s -= L[(size_t)i * m + k] * y[k];
                y[i] = s / L[(size_t)i * m + i];
            }
            for (int i = m - 1; i >= 0; --i) {
                double s = y[i];
                for (int k = i + 1; k < m; ++k) s -= L[(size_t)k * m + i] * x[k];
                x[i] = s / L[(size_t)i * m + i];
            }
            return true;
        }
        jitter = (attempt == 0) ? 1e-6 : jitter * 100.0;
    }
    for (int i = 0; i < m; ++i) {
        double d = std::max(1e-8, C[(size_t)i * m + i] + 1e-4);
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
    int kid_tu     = (Rf_length(params_r) >= 8)  ? (int)std::round(params[7]) : 0;
    int kid_td     = (Rf_length(params_r) >= 9)  ? (int)std::round(params[8]) : kid_tu;
    int kid_euc    = (Rf_length(params_r) >= 10) ? (int)std::round(params[9]) : kid_tu;

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
                cov += sig_tu * w * eval_stream_kernel(d, phi_tu, kid_tu);
            }
        }
        if (sig_td > 0.0) {
            double w = td_oo_wt[idx1 + idx2 * n];
            double d = td_oo_dist[idx1 + idx2 * n];
            if (w > 0.0 && d < 1e5) {
                cov += sig_td * w * eval_stream_kernel(d, phi_td, kid_td);
            }
        }
        if (sig_euc > 0.0) {
            double dx = (mid_x[idx1] - mid_x[idx2]) * 0.001;
            double dy = (mid_y[idx1] - mid_y[idx2]) * 0.001;
            double d_euc = std::hypot(dx, dy);
            cov += sig_euc * eval_stream_kernel(d_euc, phi_euc, kid_euc);
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
                    cand_td_w[c] = std::min(1.0 - 1e-6, std::exp(-0.001 * cand_tu_d[c]));
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
                cov += sig_tu * cand_tu_w[c] * eval_stream_kernel(cand_tu_d[c], phi_tu, kid_tu);
            }
            if (sig_td > 0.0 && cand_td_w[c] > 0.0 && cand_td_d[c] < 1e5) {
                cov += sig_td * cand_td_w[c] * eval_stream_kernel(cand_td_d[c], phi_td, kid_td);
            }
            if (sig_euc > 0.0) {
                double dx = (xi - mid_x[j]) * 0.001;
                double dy = (yi - mid_y[j]) * 0.001;
                double d_euc = std::hypot(dx, dy);
                cov += sig_euc * eval_stream_kernel(d_euc, phi_euc, kid_euc);
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

// =============================================================================
// C_s3n_optimal_design:
// Fast OpenMP Greedy Variance-Reduction (Active-Learning) Monitoring Network
// Design. Given current kriging SEs across the river network and fitted
// stream covariance parameters, sequentially selects K unmonitored reaches
// that maximize expected reduction in domain-wide kriging variance.
//
// Args:
//   parent_dn_r, net_id_r, tin_r, tout_r, d2m_km_r, up_dist_km_r, mid_x_r, mid_y_r
//   se_curr_r   : length q vector of current kriging SEs
//   params_r    : 5-vector: (sigma_sq_stream, sigma_sq_euc, tau_sq, phi, kernel_id)
//   cand_idx_r  : 0-based integer vector of candidate unmonitored reach indices
//   eval_idx_r  : 0-based integer vector of target evaluation reach indices
//   k_design_r  : integer scalar K (number of new monitoring stations to select)
// =============================================================================
SEXP C_s3n_optimal_design(
    SEXP parent_dn_r, SEXP net_id_r,
    SEXP tin_r, SEXP tout_r,
    SEXP d2m_km_r, SEXP up_dist_km_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP se_curr_r, SEXP params_r,
    SEXP cand_idx_r, SEXP eval_idx_r,
    SEXP k_design_r
) {
    int q      = Rf_length(parent_dn_r);
    int n_cand = Rf_length(cand_idx_r);
    int n_eval = Rf_length(eval_idx_r);
    int K      = std::min(Rf_asInteger(k_design_r), n_cand);

    const int*    parent_dn  = INTEGER(parent_dn_r);
    const int*    net_id     = INTEGER(net_id_r);
    const int*    tin        = INTEGER(tin_r);
    const int*    tout       = INTEGER(tout_r);
    const double* d2m_km     = REAL(d2m_km_r);
    const double* up_dist_km = REAL(up_dist_km_r);
    const double* mid_x      = REAL(mid_x_r);
    const double* mid_y      = REAL(mid_y_r);
    const double* se_curr    = REAL(se_curr_r);
    const int*    cand_idx   = INTEGER(cand_idx_r);
    const int*    eval_idx   = INTEGER(eval_idx_r);

    const double* params = REAL(params_r);
    double sig_str = std::max(0.0, params[0]);
    double sig_euc = std::max(0.0, params[1]);
    double tau_sq  = std::max(1e-6, params[2]);
    double phi     = std::max(1e-6, params[3]);
    int kernel_id  = (Rf_length(params_r) >= 5) ? (int)std::round(params[4]) : 0;
    double sig_tot = std::max(1e-8, sig_str + sig_euc);
    double w_str   = sig_str / sig_tot;
    double w_euc   = sig_euc / sig_tot;

    auto pair_corr = [&](int a, int b) -> double {
        if (a == b) return 1.0;
        double r_str = 0.0;
        if (w_str > 0.0 && net_id[a] == net_id[b]) {
            double dist = 1e6;
            if ((tin[a] <= tin[b] && tout[b] <= tout[a]) ||
                (tin[b] <= tin[a] && tout[a] <= tout[b])) {
                dist = std::max(1e-4, std::abs(up_dist_km[a] - up_dist_km[b]));
            } else {
                int anc = a;
                int guard = 0;
                while (anc >= 0 && anc < q && !(tin[anc] <= tin[b] && tout[b] <= tout[anc]) && guard < q) {
                    anc = parent_dn[anc];
                    guard++;
                }
                if (anc >= 0 && anc < q && tin[anc] <= tin[b] && tout[b] <= tout[anc]) {
                    dist = std::max(1e-4, up_dist_km[a] + up_dist_km[b] - 2.0 * d2m_km[anc]);
                }
            }
            if (dist < 1e5) {
                r_str = eval_stream_kernel(dist, phi, kernel_id);
            }
        }
        double r_euc = 0.0;
        if (w_euc > 0.0) {
            double dx = (mid_x[a] - mid_x[b]) * 0.001;
            double dy = (mid_y[a] - mid_y[b]) * 0.001;
            r_euc = eval_stream_kernel(std::hypot(dx, dy), phi, kernel_id);
        }
        return w_str * r_str + w_euc * r_euc;
    };

    std::vector<double> v_eval(n_eval, 0.0);
    double base_var_sum = 0.0;
    for (int e = 0; e < n_eval; ++e) {
        double s = se_curr[eval_idx[e]];
        v_eval[e] = std::max(0.0, s * s);
        base_var_sum += v_eval[e];
    }
    if (base_var_sum < 1e-12) base_var_sum = 1e-12;

    std::vector<double> v_cand(n_cand, 0.0);
    std::vector<char>   chosen(n_cand, 0);
    for (int c = 0; c < n_cand; ++c) {
        double s = se_curr[cand_idx[c]];
        v_cand[c] = std::max(0.0, s * s);
    }

    SEXP sel_pos_r    = PROTECT(Rf_allocVector(INTSXP, K));
    SEXP dv_r         = PROTECT(Rf_allocVector(REALSXP, K));
    SEXP pct_r        = PROTECT(Rf_allocVector(REALSXP, K));
    SEXP cum_pct_r    = PROTECT(Rf_allocVector(REALSXP, K));
    SEXP se_before_r  = PROTECT(Rf_allocVector(REALSXP, K));
    SEXP se_after_r   = PROTECT(Rf_allocVector(REALSXP, K));

    int*    sel_pos   = INTEGER(sel_pos_r);
    double* dv_out    = REAL(dv_r);
    double* pct_out   = REAL(pct_r);
    double* cum_out   = REAL(cum_pct_r);
    double* se_bef    = REAL(se_before_r);
    double* se_aft    = REAL(se_after_r);
    std::fill(sel_pos, sel_pos + K, 0);
    std::fill(dv_out,  dv_out + K,  0.0);
    std::fill(pct_out, pct_out + K, 0.0);
    std::fill(cum_out, cum_out + K, 0.0);
    std::fill(se_bef,  se_bef + K,  0.0);
    std::fill(se_aft,  se_aft + K,  0.0);

    // Precompute candidate x eval correlation matrix once when memory budget permits (<= 80 MB)
    bool use_cache = ((size_t)n_cand * (size_t)n_eval <= 20000000ULL);
    std::vector<float> R_ce;
    if (use_cache && n_cand > 0 && n_eval > 0) {
        R_ce.resize((size_t)n_cand * n_eval, 0.0f);
        #pragma omp parallel for schedule(dynamic, 4)
        for (int c = 0; c < n_cand; ++c) {
            int u = cand_idx[c];
            int net_u = net_id[u];
            float* row_c = &R_ce[(size_t)c * n_eval];
            for (int e = 0; e < n_eval; ++e) {
                int ev = eval_idx[e];
                if (w_euc == 0.0 && net_id[ev] != net_u) continue;
                row_c[e] = (float)pair_corr(ev, u);
            }
        }
    }

    double cum_dv = 0.0;

    for (int step = 0; step < K; ++step) {
        double mean_se_b = 0.0;
        for (int e = 0; e < n_eval; ++e) mean_se_b += std::sqrt(v_eval[e]);
        mean_se_b /= std::max(1, n_eval);

        std::vector<double> cand_gains(n_cand, -1.0);
        #pragma omp parallel for schedule(dynamic, 4)
        for (int c = 0; c < n_cand; ++c) {
            if (chosen[c] || v_cand[c] <= 1e-10) continue;
            int u = cand_idx[c];
            int net_u = net_id[u];
            double vc = v_cand[c];
            double denom = vc + tau_sq;
            double total_dv = 0.0;
            const float* row_c = use_cache ? &R_ce[(size_t)c * n_eval] : nullptr;
            for (int e = 0; e < n_eval; ++e) {
                if (v_eval[e] <= 1e-10) continue;
                double r = row_c ? (double)row_c[e] : ((w_euc == 0.0 && net_id[eval_idx[e]] != net_u) ? 0.0 : pair_corr(eval_idx[e], u));
                if (r <= 1e-6) continue;
                double cov_ec = std::sqrt(v_eval[e] * vc) * r;
                double dv = (cov_ec * cov_ec) / denom;
                if (dv > v_eval[e]) dv = v_eval[e];
                total_dv += dv;
            }
            cand_gains[c] = total_dv;
        }

        int best_c = -1;
        double best_dv = -1.0;
        for (int c = 0; c < n_cand; ++c) {
            if (!chosen[c] && cand_gains[c] > best_dv) {
                best_dv = cand_gains[c];
                best_c = c;
            }
        }
        if (best_c < 0) break;

        chosen[best_c] = 1;
        int u_best = cand_idx[best_c];
        int net_best = net_id[u_best];
        double vc_best = v_cand[best_c];
        double denom_best = vc_best + tau_sq;

        // Update remaining evaluation and candidate variances in parallel
        const float* best_row = use_cache ? &R_ce[(size_t)best_c * n_eval] : nullptr;
        #pragma omp parallel for schedule(static)
        for (int e = 0; e < n_eval; ++e) {
            if (v_eval[e] <= 1e-10) continue;
            double r = best_row ? (double)best_row[e] : ((w_euc == 0.0 && net_id[eval_idx[e]] != net_best) ? 0.0 : pair_corr(eval_idx[e], u_best));
            if (r <= 1e-6) continue;
            double cov_ec = std::sqrt(v_eval[e] * vc_best) * r;
            double dv = (cov_ec * cov_ec) / denom_best;
            v_eval[e] = std::max(0.0, v_eval[e] - dv);
        }
        #pragma omp parallel for schedule(static)
        for (int c = 0; c < n_cand; ++c) {
            if (chosen[c] || v_cand[c] <= 1e-10) continue;
            if (w_euc == 0.0 && net_id[cand_idx[c]] != net_best) continue;
            double r = pair_corr(cand_idx[c], u_best);
            if (r <= 1e-6) continue;
            double cov_cc = std::sqrt(v_cand[c] * vc_best) * r;
            double dv = (cov_cc * cov_cc) / denom_best;
            v_cand[c] = std::max(0.0, v_cand[c] - dv);
        }

        double mean_se_a = 0.0;
        for (int e = 0; e < n_eval; ++e) mean_se_a += std::sqrt(v_eval[e]);
        mean_se_a /= std::max(1, n_eval);

        cum_dv += std::max(0.0, best_dv);
        sel_pos[step] = best_c + 1; // 1-based position into cand_idx_r
        dv_out[step]  = std::max(0.0, best_dv);
        pct_out[step] = 100.0 * std::max(0.0, best_dv) / base_var_sum;
        cum_out[step] = 100.0 * cum_dv / base_var_sum;
        se_bef[step]  = mean_se_b;
        se_aft[step]  = mean_se_a;
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
    SET_VECTOR_ELT(out, 0, sel_pos_r);   SET_STRING_ELT(nms, 0, Rf_mkChar("cand_pos"));
    SET_VECTOR_ELT(out, 1, dv_r);        SET_STRING_ELT(nms, 1, Rf_mkChar("var_reduced"));
    SET_VECTOR_ELT(out, 2, pct_r);       SET_STRING_ELT(nms, 2, Rf_mkChar("pct_var_reduced"));
    SET_VECTOR_ELT(out, 3, cum_pct_r);   SET_STRING_ELT(nms, 3, Rf_mkChar("cum_pct_reduced"));
    SET_VECTOR_ELT(out, 4, se_before_r); SET_STRING_ELT(nms, 4, Rf_mkChar("mean_se_before"));
    SET_VECTOR_ELT(out, 5, se_after_r);  SET_STRING_ELT(nms, 5, Rf_mkChar("mean_se_after"));
    Rf_setAttrib(out, R_NamesSymbol, nms);

    UNPROTECT(8);
    return out;
}

// =============================================================================
// Heap-backed Cholesky factorization and multi-RHS solver (supports d up to ~512)
// C is row-major n x n symmetric positive semi-definite.
// On success returns L (lower triangular, row-major n x n).
// =============================================================================
static bool chol_factor_dyn(int n, const double* C, std::vector<double>& L) {
    L.assign((size_t)n * n, 0.0);
    if (n <= 0) return true;
    double trace_avg = 0.0;
    for (int i = 0; i < n; ++i) trace_avg += std::max(0.0, C[(size_t)i * n + i]);
    trace_avg = std::max(1e-6, trace_avg / n);
    double jitter = 1e-8 * trace_avg;

    for (int attempt = 0; attempt < 6; ++attempt) {
        bool ok = true;
        std::fill(L.begin(), L.end(), 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double s = C[(size_t)i * n + j];
                if (i == j) s += jitter;
                const double* Li = &L[(size_t)i * n];
                const double* Lj = &L[(size_t)j * n];
                for (int k = 0; k < j; ++k) s -= Li[k] * Lj[k];
                if (i == j) {
                    if (s <= 1e-12 * trace_avg) {
                        ok = false;
                        break;
                    }
                    L[(size_t)i * n + j] = std::sqrt(s);
                } else {
                    L[(size_t)i * n + j] = s / L[(size_t)j * n + j];
                }
            }
            if (!ok) break;
        }
        if (ok) return true;
        jitter *= 10.0;
    }
    for (int i = 0; i < n; ++i) {
        L[(size_t)i * n + i] = std::sqrt(std::max(1e-6, C[(size_t)i * n + i] + 1e-4));
    }
    return false;
}

static void chol_forward_sub(int n, const double* L, const double* b, double* y) {
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        const double* Li = &L[(size_t)i * n];
        for (int k = 0; k < i; ++k) s -= Li[k] * y[k];
        y[i] = s / Li[i];
    }
}

static void chol_back_sub(int n, const double* L, const double* y, double* x) {
    for (int i = n - 1; i >= 0; --i) {
        double s = y[i];
        for (int k = i + 1; k < n; ++k) s -= L[(size_t)k * n + i] * x[k];
        x[i] = s / L[(size_t)i * n + i];
    }
}

// =============================================================================
// Multivariate Factor-LMC NNGP Engine (Heterotopic TU + TD + Euclidean)
//
// Modes:
//   0 = Exact Multivariate NNGP -2*loglik & Spatial GLS Beta_mat (p_cov x P)
//   1 = Multivariate Co-Kriging Prediction across all q reaches (pred, se)
//   2 = Aux-Aided Leave-One-Target-Out CV across n_obs reaches (aux_pred, aux_se)
// =============================================================================
static SEXP s3n_multivariate_nngp_impl(
    SEXP X_obs_r, SEXP Y_obs_r,
    SEXP X0_r, SEXP Beta_r,
    SEXP nb_tu_r, SEXP nb_td_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP m_r,
    SEXP A_tu_r, SEXP A_td_r, SEXP A_euc_r, SEXP Psi_r,
    SEXP phis_r, SEXP mode_r,
    const double* V_obs, bool has_v
) {
    int n     = Rf_nrows(X_obs_r);
    int p_cov = Rf_ncols(X_obs_r);
    int P     = Rf_ncols(Y_obs_r);
    int m     = std::min(32, std::max(1, Rf_asInteger(m_r)));
    int mode  = Rf_asInteger(mode_r);

    const double* X_obs = REAL(X_obs_r);
    const double* Y_obs = REAL(Y_obs_r);
    const double* mid_x = REAL(mid_x_r);
    const double* mid_y = REAL(mid_y_r);

    const double* A_tu  = REAL(A_tu_r);
    const double* A_td  = REAL(A_td_r);
    const double* A_euc = REAL(A_euc_r);
    const double* Psi   = REAL(Psi_r);

    const double* phis = REAL(phis_r);
    double phi_tu  = std::max(1e-6, phis[0]);
    double phi_td  = std::max(1e-6, phis[1]);
    double phi_euc = std::max(1e-6, phis[2]);
    int kernel_id  = (Rf_length(phis_r) >= 4) ? (int)std::round(phis[3]) : 0;

    bool has_tu  = false, has_td = false, has_euc = false;
    for (int a = 0; a < P; ++a) {
        if (A_tu[a + a * P]  > 1e-12) has_tu  = true;
        if (A_td[a + a * P]  > 1e-12) has_td  = true;
        if (A_euc[a + a * P] > 1e-12) has_euc = true;
    }

    const double* tu_oo_dist = REAL(VECTOR_ELT(nb_tu_r, 0));
    const double* tu_oo_wt   = REAL(VECTOR_ELT(nb_tu_r, 1));
    const double* td_oo_dist = REAL(VECTOR_ELT(nb_td_r, 0));
    const double* td_oo_wt   = REAL(VECTOR_ELT(nb_td_r, 1));

    // Precompute spatial scalar kernel weights between observed reaches u1, u2
    auto obs_spatial_weights = [&](int u1, int u2, double& wtu_k, double& wtd_k, double& weuc_k) {
        if (u1 == u2) {
            wtu_k = has_tu ? 1.0 : 0.0;
            wtd_k = has_td ? 1.0 : 0.0;
            weuc_k = has_euc ? 1.0 : 0.0;
            return;
        }
        wtu_k = 0.0; wtd_k = 0.0; weuc_k = 0.0;
        if (has_tu) {
            double w = tu_oo_wt[u1 + u2 * n];
            double d = tu_oo_dist[u1 + u2 * n];
            if (w > 0.0 && d < 1e5) wtu_k = w * eval_stream_kernel(d, phi_tu, kernel_id);
        }
        if (has_td) {
            double w = td_oo_wt[u1 + u2 * n];
            double d = td_oo_dist[u1 + u2 * n];
            if (w > 0.0 && d < 1e5) wtd_k = w * eval_stream_kernel(d, phi_td, kernel_id);
        }
        if (has_euc) {
            double dx = (mid_x[u1] - mid_x[u2]) * 0.001;
            double dy = (mid_y[u1] - mid_y[u2]) * 0.001;
            weuc_k = eval_stream_kernel(std::hypot(dx, dy), phi_euc, kernel_id);
        }
    };

    auto obs_scalar_cov = [&](int u1, int a1, int u2, int a2,
                              double wtu_k, double wtd_k, double weuc_k) -> double {
        double cov = 0.0;
        int idx = a1 + a2 * P;
        if (has_tu && wtu_k > 0.0)   cov += wtu_k  * A_tu[idx];
        if (has_td && wtd_k > 0.0)   cov += wtd_k  * A_td[idx];
        if (has_euc && weuc_k > 0.0) cov += weuc_k * A_euc[idx];
        if (u1 == u2)                cov += Psi[idx];
        return cov;
    };

    if (mode == 0) {
        int P_beta = P * p_cov;
        int n_thr = 1;
        #ifdef _OPENMP
        n_thr = std::max(1, omp_get_max_threads());
        #endif
        std::vector<std::vector<double>> XtX_thr(n_thr, std::vector<double>((size_t)P_beta * P_beta, 0.0));
        std::vector<std::vector<double>> Xty_thr(n_thr, std::vector<double>(P_beta, 0.0));
        std::vector<double> sum_logd_thr(n_thr, 0.0);
        std::vector<double> yt_sum2_thr(n_thr, 0.0);
        std::vector<int> n_eff_thr(n_thr, 0);

        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            std::vector<double> L_N, L_D;

            #pragma omp for schedule(dynamic, 8)
            for (int i = 0; i < n; ++i) {
                std::vector<int> obs_b;
                for (int b = 0; b < P; ++b) {
                    if (R_FINITE(Y_obs[i + b * n])) obs_b.push_back(b);
                }
                int di = (int)obs_b.size();
                if (di == 0) continue;

                // Gather sorted predecessor reaches j < i so sparse heterotopic targets get up to m neighbors
                std::vector<std::pair<double, int>> cands;
                if (i > 0) {
                    cands.reserve(i);
                    for (int j = 0; j < i; ++j) {
                        double wtu_k, wtd_k, weuc_k;
                        obs_spatial_weights(i, j, wtu_k, wtd_k, weuc_k);
                        double score = wtu_k + wtd_k + weuc_k;
                        if (score > 1e-10) cands.push_back({-score, j});
                    }
                    std::sort(cands.begin(), cands.end());
                }

                // Flatten observed neighbor entries (u, a) capping each target a at m observations
                std::vector<int> nb_u, nb_a, var_cnt(P, 0);
                std::vector<double> nb_wtu, nb_wtd, nb_weuc;
                for (size_t c = 0; c < cands.size(); ++c) {
                    int u = cands[c].second;
                    double wtu_k, wtd_k, weuc_k;
                    obs_spatial_weights(i, u, wtu_k, wtd_k, weuc_k);
                    bool all_full = true;
                    for (int a = 0; a < P; ++a) {
                        if (var_cnt[a] < m) {
                            all_full = false;
                            if (R_FINITE(Y_obs[u + a * n])) {
                                nb_u.push_back(u);
                                nb_a.push_back(a);
                                nb_wtu.push_back(wtu_k);
                                nb_wtd.push_back(wtd_k);
                                nb_weuc.push_back(weuc_k);
                                var_cnt[a]++;
                            }
                        }
                    }
                    if (all_full) break;
                }
                int K = (int)nb_u.size();

                // Block covariance C_ii (di x di)
                std::vector<double> C_ii((size_t)di * di, 0.0);
                for (int r = 0; r < di; ++r) {
                    for (int s = 0; s <= r; ++s) {
                        int b1 = obs_b[r], b2 = obs_b[s];
                        double v = A_tu[b1 + b2 * P] + A_td[b1 + b2 * P] +
                                   A_euc[b1 + b2 * P] + Psi[b1 + b2 * P];
                        if (has_v && r == s) {
                            double vi = V_obs[i + b1 * n];
                            if (R_FINITE(vi) && vi > 0.0) v += vi;
                        }
                        C_ii[r * di + s] = v;
                        C_ii[s * di + r] = v;
                    }
                }

                // Residual block dy_i (di) and dX_i (di x P_beta)
                std::vector<double> dy(di, 0.0);
                std::vector<double> dX((size_t)di * P_beta, 0.0);
                for (int r = 0; r < di; ++r) {
                    int br = obs_b[r];
                    dy[r] = Y_obs[i + br * n];
                    for (int pc = 0; pc < p_cov; ++pc) {
                        dX[r + (br * p_cov + pc) * di] = X_obs[i + pc * n];
                    }
                }

                std::vector<double> D_i = C_ii;
                if (K > 0) {
                    std::vector<double> C_NN((size_t)K * K, 0.0);
                    for (int s1 = 0; s1 < K; ++s1) {
                        int u1 = nb_u[s1], a1 = nb_a[s1];
                        double diag_v = A_tu[a1 + a1 * P] + A_td[a1 + a1 * P] +
                                        A_euc[a1 + a1 * P] + Psi[a1 + a1 * P];
                        if (has_v) {
                            double vi = V_obs[u1 + a1 * n];
                            if (R_FINITE(vi) && vi > 0.0) diag_v += vi;
                        }
                        C_NN[(size_t)s1 * K + s1] = diag_v;
                        for (int s2 = s1 + 1; s2 < K; ++s2) {
                            int u2 = nb_u[s2], a2 = nb_a[s2];
                            double w1, w2, w3;
                            obs_spatial_weights(u1, u2, w1, w2, w3);
                            double cov = obs_scalar_cov(u1, a1, u2, a2, w1, w2, w3);
                            C_NN[(size_t)s1 * K + s2] = cov;
                            C_NN[(size_t)s2 * K + s1] = cov;
                        }
                    }
                    chol_factor_dyn(K, C_NN.data(), L_N);

                    std::vector<double> tmp(K), U_i((size_t)K * di, 0.0);
                    std::vector<double> c_col(K);
                    for (int r = 0; r < di; ++r) {
                        int br = obs_b[r];
                        for (int s = 0; s < K; ++s) {
                            c_col[s] = obs_scalar_cov(nb_u[s], nb_a[s], i, br, nb_wtu[s], nb_wtd[s], nb_weuc[s]);
                        }
                        double* col_Ur = &U_i[(size_t)r * K];
                        chol_forward_sub(K, L_N.data(), c_col.data(), col_Ur);
                        chol_back_sub(K, L_N.data(), col_Ur, tmp.data());
                        for (int s = 0; s < K; ++s) {
                            dy[r] -= tmp[s] * Y_obs[nb_u[s] + nb_a[s] * n];
                            for (int pc = 0; pc < p_cov; ++pc) {
                                dX[r + (nb_a[s] * p_cov + pc) * di] -= tmp[s] * X_obs[nb_u[s] + pc * n];
                            }
                        }
                    }
                    // D_i = C_ii - B_i^T * C_Ni = C_ii - U_i^T * U_i
                    for (int r1 = 0; r1 < di; ++r1) {
                        const double* u1_ptr = &U_i[(size_t)r1 * K];
                        for (int r2 = r1; r2 < di; ++r2) {
                            const double* u2_ptr = &U_i[(size_t)r2 * K];
                            double dot = 0.0;
                            for (int s = 0; s < K; ++s) dot += u1_ptr[s] * u2_ptr[s];
                            D_i[r1 * di + r2] -= dot;
                            D_i[r2 * di + r1] = D_i[r1 * di + r2];
                        }
                    }
                }

                chol_factor_dyn(di, D_i.data(), L_D);
                double logd = 0.0;
                for (int r = 0; r < di; ++r) logd += 2.0 * std::log(L_D[r * di + r]);
                sum_logd_thr[tid] += logd;
                n_eff_thr[tid]    += di;

                std::vector<double> yt(di), Xt((size_t)di * P_beta);
                chol_forward_sub(di, L_D.data(), dy.data(), yt.data());
                for (int r = 0; r < di; ++r) yt_sum2_thr[tid] += yt[r] * yt[r];

                std::vector<double> col_x(di), col_xt(di);
                for (int b = 0; b < P_beta; ++b) {
                    for (int r = 0; r < di; ++r) col_x[r] = dX[r + b * di];
                    chol_forward_sub(di, L_D.data(), col_x.data(), col_xt.data());
                    for (int r = 0; r < di; ++r) Xt[r + b * di] = col_xt[r];
                }

                for (int b1 = 0; b1 < P_beta; ++b1) {
                    double sy = 0.0;
                    for (int r = 0; r < di; ++r) sy += Xt[r + b1 * di] * yt[r];
                    Xty_thr[tid][b1] += sy;
                    for (int b2 = b1; b2 < P_beta; ++b2) {
                        double sx = 0.0;
                        for (int r = 0; r < di; ++r) sx += Xt[r + b1 * di] * Xt[r + b2 * di];
                        XtX_thr[tid][(size_t)b1 * P_beta + b2] += sx;
                        if (b1 != b2) XtX_thr[tid][(size_t)b2 * P_beta + b1] += sx;
                    }
                }
            }
        }

        std::vector<double> XtX((size_t)P_beta * P_beta, 0.0), Xty(P_beta, 0.0);
        double sum_logd = 0.0, yt_sum2 = 0.0;
        int n_eff = 0;
        for (int t = 0; t < n_thr; ++t) {
            sum_logd += sum_logd_thr[t];
            yt_sum2  += yt_sum2_thr[t];
            n_eff    += n_eff_thr[t];
            for (int b = 0; b < P_beta; ++b) Xty[b] += Xty_thr[t][b];
            for (size_t idx = 0; idx < XtX.size(); ++idx) XtX[idx] += XtX_thr[t][idx];
        }

        std::vector<double> L_B, y_b(P_beta), beta_vec(P_beta, 0.0);
        chol_factor_dyn(P_beta, XtX.data(), L_B);
        chol_forward_sub(P_beta, L_B.data(), Xty.data(), y_b.data());
        chol_back_sub(P_beta, L_B.data(), y_b.data(), beta_vec.data());

        // Exact GLS RSS = ||y_tilde - X_tilde * beta||^2 = yt_sum2 - 2 * beta^T Xty + beta^T XtX beta
        double beta_Xty = 0.0, beta_XtX_beta = 0.0;
        for (int b1 = 0; b1 < P_beta; ++b1) {
            beta_Xty += beta_vec[b1] * Xty[b1];
            for (int b2 = 0; b2 < P_beta; ++b2) {
                beta_XtX_beta += beta_vec[b1] * XtX[(size_t)b1 * P_beta + b2] * beta_vec[b2];
            }
        }
        double rss = std::max(1e-10, yt_sum2 - 2.0 * beta_Xty + beta_XtX_beta);
        double neg2ll = rss + sum_logd;

        SEXP beta_mat_r = PROTECT(Rf_allocMatrix(REALSXP, p_cov, P));
        double* bm = REAL(beta_mat_r);
        for (int v = 0; v < P; ++v) {
            for (int pc = 0; pc < p_cov; ++pc) bm[pc + v * p_cov] = beta_vec[v * p_cov + pc];
        }

        SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
        SET_VECTOR_ELT(out, 0, Rf_ScalarReal(neg2ll));
        SET_VECTOR_ELT(out, 1, Rf_ScalarReal(rss / std::max(1, n_eff)));
        SET_VECTOR_ELT(out, 2, beta_mat_r);
        SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(n_eff));
        SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(nms, 0, Rf_mkChar("neg2loglik"));
        SET_STRING_ELT(nms, 1, Rf_mkChar("whitened_var"));
        SET_STRING_ELT(nms, 2, Rf_mkChar("Beta"));
        SET_STRING_ELT(nms, 3, Rf_mkChar("n_eff"));
        Rf_setAttrib(out, R_NamesSymbol, nms);
        UNPROTECT(3);
        return out;
    }

    // Precompute observed residual matrix r_obs (n x P) for Modes 1 & 2
    const double* Beta = REAL(Beta_r);
    std::vector<double> r_obs((size_t)n * P, NA_REAL);
    for (int i = 0; i < n; ++i) {
        for (int b = 0; b < P; ++b) {
            if (!R_FINITE(Y_obs[i + b * n])) continue;
            double mu = 0.0;
            for (int pc = 0; pc < p_cov; ++pc) mu += X_obs[i + pc * n] * Beta[pc + b * p_cov];
            r_obs[i + b * n] = Y_obs[i + b * n] - mu;
        }
    }

    if (mode == 1) {
        // Mode 1: Full network multivariate Co-Kriging prediction across all q reaches
        int q = Rf_nrows(X0_r);
        const double* X0 = REAL(X0_r);

        SEXP pred_nb_tu = VECTOR_ELT(nb_tu_r, 3);
        SEXP pred_nb_td = VECTOR_ELT(nb_td_r, 3);
        const int*    tu_p_idx  = INTEGER(VECTOR_ELT(pred_nb_tu, 0));
        const double* tu_p_dist = REAL(VECTOR_ELT(pred_nb_tu, 1));
        const double* tu_p_wt   = REAL(VECTOR_ELT(pred_nb_tu, 2));
        const int*    td_p_idx  = INTEGER(VECTOR_ELT(pred_nb_td, 0));
        const double* td_p_dist = REAL(VECTOR_ELT(pred_nb_td, 1));
        const double* td_p_wt   = REAL(VECTOR_ELT(pred_nb_td, 2));
        int m_stride = Rf_length(VECTOR_ELT(pred_nb_tu, 0)) / std::max(1, q);

        SEXP pred_r = PROTECT(Rf_allocMatrix(REALSXP, q, P));
        SEXP se_r   = PROTECT(Rf_allocMatrix(REALSXP, q, P));
        double* pred = REAL(pred_r);
        double* se   = REAL(se_r);

        #pragma omp parallel
        {
            std::vector<double> L_N;
            #pragma omp for schedule(dynamic, 32)
            for (int i = 0; i < q; ++i) {
                int base_p = i * m_stride;
                int cand_ids[128];
                double cand_wtu[128], cand_wtd[128], cand_weuc[128], cand_sc[128];
                int n_cand = 0;

                auto add_cand = [&](int u) -> int {
                    if (u < 0 || u >= n) return -1;
                    for (int c = 0; c < n_cand; ++c) if (cand_ids[c] == u) return c;
                    if (n_cand >= 128) return -1;
                    int c = n_cand++;
                    cand_ids[c] = u;
                    cand_wtu[c] = 0.0; cand_wtd[c] = 0.0; cand_weuc[c] = 0.0; cand_sc[c] = 0.0;
                    return c;
                };

                for (int k = 0; k < m_stride; ++k) {
                    int c = add_cand(tu_p_idx[base_p + k]);
                    if (c >= 0) {
                        double d = tu_p_dist[base_p + k], w = tu_p_wt[base_p + k];
                        if (w > 0.0 && d < 1e5 && has_tu) cand_wtu[c] = w * eval_stream_kernel(d, phi_tu, kernel_id);
                        if (w > 0.0 && d < 1e5 && has_td) {
                            double w_td_raw = std::min(1.0 - 1e-6, std::exp(-0.001 * d));
                            cand_wtd[c] = w_td_raw * eval_stream_kernel(d, phi_td, kernel_id);
                        }
                    }
                }
                for (int k = 0; k < m_stride; ++k) {
                    int c = add_cand(td_p_idx[base_p + k]);
                    if (c >= 0) {
                        double d = td_p_dist[base_p + k], w = td_p_wt[base_p + k];
                        if (w > 0.0 && d < 1e5 && has_td) {
                            cand_wtd[c] = std::max(cand_wtd[c], w * eval_stream_kernel(d, phi_td, kernel_id));
                        }
                    }
                }
                if (has_euc) {
                    for (int c = 0; c < n_cand; ++c) {
                        int u = cand_ids[c];
                        cand_weuc[c] = eval_stream_kernel(
                            std::hypot((mid_x[i] - mid_x[u]) * 0.001, (mid_y[i] - mid_y[u]) * 0.001),
                            phi_euc, kernel_id);
                    }
                    double best_d2[32]; int best_j[32]; int nb_e = 0;
                    int e_cap = std::min(m, 32);
                    for (int j = 0; j < n; ++j) {
                        double d2 = std::pow((mid_x[i] - mid_x[j]) * 0.001, 2) + std::pow((mid_y[i] - mid_y[j]) * 0.001, 2);
                        if (nb_e < e_cap) {
                            int pos = nb_e++;
                            while (pos > 0 && best_d2[pos - 1] > d2) { best_d2[pos] = best_d2[pos - 1]; best_j[pos] = best_j[pos - 1]; --pos; }
                            best_d2[pos] = d2; best_j[pos] = j;
                        } else if (d2 < best_d2[nb_e - 1]) {
                            int pos = nb_e - 1;
                            while (pos > 0 && best_d2[pos - 1] > d2) { best_d2[pos] = best_d2[pos - 1]; best_j[pos] = best_j[pos - 1]; --pos; }
                            best_d2[pos] = d2; best_j[pos] = j;
                        }
                    }
                    for (int e = 0; e < nb_e; ++e) {
                        int c = add_cand(best_j[e]);
                        if (c >= 0) cand_weuc[c] = eval_stream_kernel(std::sqrt(best_d2[e]), phi_euc, kernel_id);
                    }
                }

                std::vector<std::pair<double, int>> ranked;
                ranked.reserve(n_cand);
                for (int c = 0; c < n_cand; ++c) {
                    cand_sc[c] = cand_wtu[c] + cand_wtd[c] + cand_weuc[c];
                    if (cand_sc[c] > 1e-10) ranked.push_back({-cand_sc[c], c});
                }
                std::sort(ranked.begin(), ranked.end());

                std::vector<int> nb_u, nb_a, var_cnt(P, 0);
                std::vector<double> nb_wtu, nb_wtd, nb_weuc, r_vec;
                for (size_t r = 0; r < ranked.size(); ++r) {
                    int c = ranked[r].second;
                    int u = cand_ids[c];
                    bool all_full = true;
                    for (int a = 0; a < P; ++a) {
                        if (var_cnt[a] < m) {
                            all_full = false;
                            double rv = r_obs[u + a * n];
                            if (R_FINITE(rv)) {
                                nb_u.push_back(u); nb_a.push_back(a);
                                nb_wtu.push_back(cand_wtu[c]); nb_wtd.push_back(cand_wtd[c]); nb_weuc.push_back(cand_weuc[c]);
                                r_vec.push_back(rv);
                                var_cnt[a]++;
                            }
                        }
                    }
                    if (all_full) break;
                }

                int K = (int)nb_u.size();
                std::vector<double> alpha(K, 0.0);
                if (K > 0) {
                    std::vector<double> C_NN((size_t)K * K, 0.0);
                    for (int s1 = 0; s1 < K; ++s1) {
                        int u1 = nb_u[s1], a1 = nb_a[s1];
                        double diag_v = A_tu[a1 + a1 * P] + A_td[a1 + a1 * P] + A_euc[a1 + a1 * P] + Psi[a1 + a1 * P];
                        if (has_v) {
                            double vi = V_obs[u1 + a1 * n];
                            if (R_FINITE(vi) && vi > 0.0) diag_v += vi;
                        }
                        C_NN[(size_t)s1 * K + s1] = diag_v;
                        for (int s2 = s1 + 1; s2 < K; ++s2) {
                            int u2 = nb_u[s2], a2 = nb_a[s2];
                            double w1, w2, w3;
                            obs_spatial_weights(u1, u2, w1, w2, w3);
                            double cov = obs_scalar_cov(u1, a1, u2, a2, w1, w2, w3);
                            C_NN[(size_t)s1 * K + s2] = cov;
                            C_NN[(size_t)s2 * K + s1] = cov;
                        }
                    }
                    chol_factor_dyn(K, C_NN.data(), L_N);
                    std::vector<double> yt_n(K);
                    chol_forward_sub(K, L_N.data(), r_vec.data(), yt_n.data());
                    chol_back_sub(K, L_N.data(), yt_n.data(), alpha.data());
                }

                std::vector<double> c_b(K), v_b(K);
                for (int b = 0; b < P; ++b) {
                    double mu_b = 0.0;
                    for (int pc = 0; pc < p_cov; ++pc) mu_b += X0[i + pc * q] * Beta[pc + b * p_cov];
                    double var_b = A_tu[b + b * P] + A_td[b + b * P] + A_euc[b + b * P] + Psi[b + b * P];
                    double shift = 0.0, red = 0.0;
                    if (K > 0) {
                        for (int s = 0; s < K; ++s) {
                            double w_tu_s  = nb_wtu[s];
                            double w_td_s  = nb_wtd[s];
                            double w_euc_s = nb_weuc[s];
                            double cov = 0.0;
                            int idx = nb_a[s] + b * P;
                            if (has_tu  && w_tu_s  > 0.0) cov += w_tu_s  * A_tu[idx];
                            if (has_td  && w_td_s  > 0.0) cov += w_td_s  * A_td[idx];
                            if (has_euc && w_euc_s > 0.0) cov += w_euc_s * A_euc[idx];
                            if (nb_u[s] == i && nb_a[s] != b) cov += Psi[idx];
                            c_b[s] = cov;
                            shift += cov * alpha[s];
                        }
                        chol_forward_sub(K, L_N.data(), c_b.data(), v_b.data());
                        for (int s = 0; s < K; ++s) red += v_b[s] * v_b[s];
                    }
                    pred[i + b * q] = mu_b + shift;
                    se[i + b * q]   = std::sqrt(std::max(1e-8, var_b - red));
                }
            }
        }

        SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(out, 0, pred_r);
        SET_VECTOR_ELT(out, 1, se_r);
        SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(nms, 0, Rf_mkChar("prediction"));
        SET_STRING_ELT(nms, 1, Rf_mkChar("se"));
        Rf_setAttrib(out, R_NamesSymbol, nms);
        UNPROTECT(4);
        return out;
    }

    // Mode 2: Aux-Aided Leave-One-Target-Out CV across n observed stations
    SEXP aux_pred_r = PROTECT(Rf_allocMatrix(REALSXP, n, P));
    SEXP aux_se_r   = PROTECT(Rf_allocMatrix(REALSXP, n, P));
    double* aux_pred = REAL(aux_pred_r);
    double* aux_se   = REAL(aux_se_r);
    for (size_t k = 0; k < (size_t)n * P; ++k) {
        aux_pred[k] = NA_REAL;
        aux_se[k]   = NA_REAL;
    }

    #pragma omp parallel
    {
        std::vector<double> L_N;
        #pragma omp for schedule(dynamic, 8)
        for (int i = 0; i < n; ++i) {
            std::vector<std::pair<double, int>> cands;
            cands.reserve(n);
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                double w1, w2, w3;
                obs_spatial_weights(i, j, w1, w2, w3);
                double sc = w1 + w2 + w3;
                if (sc > 1e-10) cands.push_back({-sc, j});
            }
            std::sort(cands.begin(), cands.end());

            for (int b = 0; b < P; ++b) {
                if (!R_FINITE(Y_obs[i + b * n])) continue;
                std::vector<int> nb_u, nb_a, var_cnt(P, 0);
                std::vector<double> nb_w1, nb_w2, nb_w3, r_vec;

                // Co-located auxiliary observations at reach i itself (a != b)
                for (int a = 0; a < P; ++a) {
                    if (a == b || !R_FINITE(r_obs[i + a * n])) continue;
                    nb_u.push_back(i); nb_a.push_back(a);
                    nb_w1.push_back(has_tu ? 1.0 : 0.0);
                    nb_w2.push_back(has_td ? 1.0 : 0.0);
                    nb_w3.push_back(has_euc ? 1.0 : 0.0);
                    r_vec.push_back(r_obs[i + a * n]);
                }
                // Neighboring station observations (u != i), capped at m per target a
                for (size_t c = 0; c < cands.size(); ++c) {
                    int u = cands[c].second;
                    double w1, w2, w3;
                    obs_spatial_weights(i, u, w1, w2, w3);
                    bool all_full = true;
                    for (int a = 0; a < P; ++a) {
                        if (var_cnt[a] < m) {
                            all_full = false;
                            double rv = r_obs[u + a * n];
                            if (R_FINITE(rv)) {
                                nb_u.push_back(u); nb_a.push_back(a);
                                nb_w1.push_back(w1); nb_w2.push_back(w2); nb_w3.push_back(w3);
                                r_vec.push_back(rv);
                                var_cnt[a]++;
                            }
                        }
                    }
                    if (all_full) break;
                }

                double mu_b = 0.0;
                for (int pc = 0; pc < p_cov; ++pc) mu_b += X_obs[i + pc * n] * Beta[pc + b * p_cov];
                double var_b = A_tu[b + b * P] + A_td[b + b * P] + A_euc[b + b * P] + Psi[b + b * P];

                int K = (int)nb_u.size();
                if (K == 0) {
                    aux_pred[i + b * n] = mu_b;
                    aux_se[i + b * n]   = std::sqrt(std::max(1e-8, var_b));
                    continue;
                }

                std::vector<double> C_NN((size_t)K * K, 0.0), c_b(K), alpha(K), tmp(K);
                for (int s1 = 0; s1 < K; ++s1) {
                    int u1 = nb_u[s1], a1 = nb_a[s1];
                    double diag_v = A_tu[a1 + a1 * P] + A_td[a1 + a1 * P] + A_euc[a1 + a1 * P] + Psi[a1 + a1 * P];
                    if (has_v) {
                        double vi = V_obs[u1 + a1 * n];
                        if (R_FINITE(vi) && vi > 0.0) diag_v += vi;
                    }
                    C_NN[(size_t)s1 * K + s1] = diag_v;
                    c_b[s1] = obs_scalar_cov(u1, a1, i, b, nb_w1[s1], nb_w2[s1], nb_w3[s1]);
                    for (int s2 = s1 + 1; s2 < K; ++s2) {
                        int u2 = nb_u[s2], a2 = nb_a[s2];
                        double w1, w2, w3;
                        obs_spatial_weights(u1, u2, w1, w2, w3);
                        double cov = obs_scalar_cov(u1, a1, u2, a2, w1, w2, w3);
                        C_NN[(size_t)s1 * K + s2] = cov;
                        C_NN[(size_t)s2 * K + s1] = cov;
                    }
                }
                chol_factor_dyn(K, C_NN.data(), L_N);
                chol_forward_sub(K, L_N.data(), r_vec.data(), tmp.data());
                chol_back_sub(K, L_N.data(), tmp.data(), alpha.data());

                double shift = 0.0;
                for (int s = 0; s < K; ++s) shift += c_b[s] * alpha[s];
                chol_forward_sub(K, L_N.data(), c_b.data(), tmp.data());
                double red = 0.0;
                for (int s = 0; s < K; ++s) red += tmp[s] * tmp[s];

                aux_pred[i + b * n] = mu_b + shift;
                aux_se[i + b * n]   = std::sqrt(std::max(1e-8, var_b - red));
            }
        }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, aux_pred_r);
    SET_VECTOR_ELT(out, 1, aux_se_r);
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("aux_pred"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("aux_se"));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(4);
    return out;
}

SEXP C_s3n_multivariate_nngp(
    SEXP X_obs_r, SEXP Y_obs_r,
    SEXP X0_r, SEXP Beta_r,
    SEXP nb_tu_r, SEXP nb_td_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP m_r,
    SEXP A_tu_r, SEXP A_td_r, SEXP A_euc_r, SEXP Psi_r,
    SEXP phis_r, SEXP mode_r
) {
    return s3n_multivariate_nngp_impl(
        X_obs_r, Y_obs_r, X0_r, Beta_r, nb_tu_r, nb_td_r,
        mid_x_r, mid_y_r, m_r, A_tu_r, A_td_r, A_euc_r, Psi_r, phis_r, mode_r,
        nullptr, false
    );
}

SEXP C_s3n_multivariate_nngp_weighted(
    SEXP X_obs_r, SEXP Y_obs_r,
    SEXP X0_r, SEXP Beta_r,
    SEXP nb_tu_r, SEXP nb_td_r,
    SEXP mid_x_r, SEXP mid_y_r,
    SEXP m_r,
    SEXP A_tu_r, SEXP A_td_r, SEXP A_euc_r, SEXP Psi_r,
    SEXP phis_r, SEXP mode_r, SEXP V_obs_r
) {
    int n = Rf_nrows(X_obs_r);
    int P = Rf_ncols(Y_obs_r);
    bool has_v = (!Rf_isNull(V_obs_r) && Rf_isReal(V_obs_r) && Rf_length(V_obs_r) == n * P);
    const double* V_obs = has_v ? REAL(V_obs_r) : nullptr;
    return s3n_multivariate_nngp_impl(
        X_obs_r, Y_obs_r, X0_r, Beta_r, nb_tu_r, nb_td_r,
        mid_x_r, mid_y_r, m_r, A_tu_r, A_td_r, A_euc_r, Psi_r, phis_r, mode_r,
        V_obs, has_v
    );
}

}



