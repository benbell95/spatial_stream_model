#include <R.h>
#include <Rinternals.h>
#include <cmath>
#include <vector>
#include <algorithm>
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
                    wt = 1.0;
                } else {
                    double afv_up = std::min(log_afv[a], log_afv[b]);
                    double afv_dn = std::max(log_afv[a], log_afv[b]);
                    wt = std::exp(0.5 * (afv_up - afv_dn));
                }
            } else {
                int anc = parent_dn[a];
                while (anc >= 0 && !(tin[anc] <= tb && tout[b] <= tout[anc])) {
                    anc = parent_dn[anc];
                }
                double ncj_up = (anc >= 0) ? d2m_km[anc] : 0.0;
                dist = std::max(1e-4, up_dist_km[a] + up_dist_km[b] - 2.0 * ncj_up);
                if (model_type == 1) {
                    fc = true;
                    wt = 1.0;
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
    SEXP pred_nnIndx_r = PROTECT(Rf_allocVector(INTSXP, q * m));
    SEXP pred_nnDist_r = PROTECT(Rf_allocVector(REALSXP, q * m));
    SEXP pred_nnWght_r = PROTECT(Rf_allocVector(REALSXP, q * m));
    int* pred_nnIndx = INTEGER(pred_nnIndx_r);
    double* pred_nnDist = REAL(pred_nnDist_r);
    double* pred_nnWght = REAL(pred_nnWght_r);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < q; ++i) {
        std::vector<Cand> cands(n);
        if (i < n) {
            for (int j = 0; j < n; ++j) {
                double wt = obsobs_wt[i + j * n];
                double dist = obsobs_dist[i + j * n];
                cands[j] = {wt > 0.0, dist, wt, j};
            }
        } else {
            for (int j = 0; j < n; ++j) {
                bool fc; double dist, wt;
                calc_pair(i, j, fc, dist, wt);
                cands[j] = {fc, dist, wt, j};
            }
        }
        std::partial_sort(cands.begin(), cands.begin() + m, cands.end());
        int base_p = i * m;
        for (int k = 0; k < m; ++k) {
            pred_nnIndx[base_p + k] = cands[k].idx;
            pred_nnDist[base_p + k] = cands[k].dist;
            pred_nnWght[base_p + k] = cands[k].wt;
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

}

