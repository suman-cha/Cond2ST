// mmd_kernels.cpp — Rcpp acceleration for R/mmd.R.
//
// Mirrors the pure-R implementations of MMDl, MMDb, and quadratic_MMD in
// experiments/R/mmd.R. Numerics are not byte-for-byte with the R path;
// summation order of the quadratic forms differs. Agreement is held to
// double-precision round-off (all.equal tolerance ~1e-12) and verified by
// experiments/tests/test_mmd_cpp_equiv.R.
//
// Loaded automatically by R/mmd.R when this file is found via
// Rcpp::sourceCpp. Must keep its exported function names in sync with the
// runtime checks in R/mmd.R:
//   - mmd_linear_cpp       (dispatch guard in MMDl)
//   - mmd_block_cpp        (dispatch guard in MMDb)
//   - bootstrap_mmd_cpp    (dispatch guard in quadratic_MMD)

#include <Rcpp.h>
using namespace Rcpp;

static inline double gkernel_rows(const NumericMatrix& X, int i,
                                  const NumericMatrix& Y, int j, double h) {
    const int d = X.ncol();
    double dist = 0.0;
    for (int k = 0; k < d; ++k) {
        const double diff = X(i, k) - Y(j, k);
        dist += diff * diff;
    }
    return std::exp(-dist / (2.0 * h * h));
}

static inline double gkernel_scalar(double a, double b, double h) {
    const double diff = a - b;
    return std::exp(-(diff * diff) / (2.0 * h * h));
}

// [[Rcpp::export]]
double mmd_linear_cpp(NumericMatrix x12, NumericMatrix x22,
                      NumericVector y12, NumericVector y22,
                      NumericVector r_X, double h_x, double h_y) {
    const int n = y12.size();
    const int m = n / 2;
    if (m < 2) stop("mmd_linear_cpp: n/2 must be >= 2");

    NumericVector S(m);
    for (int i = 0; i < m; ++i) {
        const int ipm = i + m;
        const double k_zz = gkernel_rows(x12, i, x12, ipm, h_x) * gkernel_scalar(y12[i],  y12[ipm], h_y);
        const double k_ww = gkernel_rows(x22, i, x22, ipm, h_x) * gkernel_scalar(y22[i],  y22[ipm], h_y);
        const double k_wz = gkernel_rows(x22, i, x12, ipm, h_x) * gkernel_scalar(y22[i],  y12[ipm], h_y);
        const double k_zw = gkernel_rows(x12, i, x22, ipm, h_x) * gkernel_scalar(y12[i],  y22[ipm], h_y);
        S[i] = k_zz + r_X[i] * r_X[ipm] * k_ww - r_X[i] * k_wz - r_X[ipm] * k_zw;
    }

    double sum_S = 0.0;
    for (int i = 0; i < m; ++i) sum_S += S[i];
    const double S_bar = sum_S / static_cast<double>(m);

    double sum_sq = 0.0;
    for (int i = 0; i < m; ++i) {
        const double d = S[i] - S_bar;
        sum_sq += d * d;
    }
    const double sigma_hat = sum_sq / static_cast<double>(m - 1);

    return std::sqrt(static_cast<double>(m)) * S_bar / std::sqrt(sigma_hat);
}

// [[Rcpp::export]]
double mmd_block_cpp(NumericMatrix x12, NumericMatrix x22,
                     NumericVector y12, NumericVector y22,
                     NumericVector r_X, int B_size, double h_x, double h_y) {
    const int n = y12.size();
    const int S = n / B_size;
    if (S < 2) stop("mmd_block_cpp: block size too large for given sample size");

    NumericVector MMD_values(S);
    const double n_pairs = static_cast<double>(B_size) * (B_size - 1) / 2.0;

    for (int b = 0; b < S; ++b) {
        const int off = b * B_size;
        double sum_k = 0.0;
        for (int i = 0; i < B_size - 1; ++i) {
            const int ii = off + i;
            for (int j = i + 1; j < B_size; ++j) {
                const int jj = off + j;
                const double k_zz = gkernel_rows(x12, ii, x12, jj, h_x) * gkernel_scalar(y12[ii], y12[jj], h_y);
                const double k_ww = gkernel_rows(x22, ii, x22, jj, h_x) * gkernel_scalar(y22[ii], y22[jj], h_y);
                const double k_wz = gkernel_rows(x22, ii, x12, jj, h_x) * gkernel_scalar(y22[ii], y12[jj], h_y);
                const double k_zw = gkernel_rows(x12, ii, x22, jj, h_x) * gkernel_scalar(y12[ii], y22[jj], h_y);
                sum_k += k_zz + r_X[ii] * r_X[jj] * k_ww - r_X[ii] * k_wz - r_X[jj] * k_zw;
            }
        }
        MMD_values[b] = sum_k / n_pairs;
    }

    double sum_M = 0.0;
    for (int b = 0; b < S; ++b) sum_M += MMD_values[b];
    const double S_bar = sum_M / static_cast<double>(S);

    double sum_sq = 0.0;
    for (int b = 0; b < S; ++b) {
        const double d = MMD_values[b] - S_bar;
        sum_sq += d * d;
    }
    const double sigma_hat = std::sqrt(sum_sq / static_cast<double>(S - 1));

    if (sigma_hat > 0.0) return std::sqrt(static_cast<double>(S)) * S_bar / sigma_hat;
    return 0.0;
}

// [[Rcpp::export]]
List bootstrap_mmd_cpp(NumericMatrix x12, NumericMatrix x22,
                       NumericVector y12, NumericVector y22,
                       NumericVector r_X, double h_x, double h_y, int B) {
    const int n = y12.size();

    // Build H_hat (symmetric, zero diagonal).
    NumericMatrix H_hat(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double k_zz = gkernel_rows(x12, i, x12, j, h_x) * gkernel_scalar(y12[i], y12[j], h_y);
            const double k_ww = gkernel_rows(x22, i, x22, j, h_x) * gkernel_scalar(y22[i], y22[j], h_y);
            const double k_wz = gkernel_rows(x22, i, x12, j, h_x) * gkernel_scalar(y22[i], y12[j], h_y);
            const double k_zw = gkernel_rows(x12, i, x22, j, h_x) * gkernel_scalar(y12[i], y22[j], h_y);
            H_hat(i, j) = k_zz + r_X[i] * r_X[j] * k_ww - r_X[i] * k_wz - r_X[j] * k_zw;
        }
    }

    double sum_H = 0.0;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            sum_H += H_hat(i, j);                    // column-major sum order
    const double denom = static_cast<double>(n) * (n - 1);
    const double obs_stat = sum_H / denom;

    NumericVector bootstrap_stats(B);
    GetRNGstate();
    NumericVector W(n);
    for (int b = 0; b < B; ++b) {
        for (int i = 0; i < n; ++i) W[i] = ::norm_rand();
        double s = 0.0;
        for (int j = 0; j < n; ++j) {
            const double Wj = W[j];
            for (int i = 0; i < n; ++i) {
                s += W[i] * Wj * H_hat(i, j);
            }
        }
        bootstrap_stats[b] = s / denom;
    }
    PutRNGstate();

    return List::create(Named("obs_stat") = obs_stat,
                        Named("bootstrap_stats") = bootstrap_stats);
}
