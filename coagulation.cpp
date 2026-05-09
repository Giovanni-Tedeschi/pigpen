#include "coagulation.h"
#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm>

// ---------------------------------------------------------------
// helpers
// ---------------------------------------------------------------
static inline double theta_fn(double x) { return x >= 0.0 ? 1.0 : 0.0; }
static inline int    delta_fn(int i, int j) { return i == j ? 1 : 0; }

int CoagKernel::find_bin(const std::vector<double>& m, double val)
{
    auto it = std::lower_bound(m.begin(), m.end(), val);
    if (it != m.end() && *it == val)
        return (int)std::distance(m.begin(), it);
    return (int)std::distance(m.begin(), it) - 1;
}

// ---------------------------------------------------------------
// MassGrid::build
// ---------------------------------------------------------------
void MassGrid::build_from_grid(const std::vector<double>& m_in,
                                const std::vector<double>& dm_in)
{
    m  = m_in;
    dm = dm_in;
    N_m   = (int)m.size();
    m_min = m.front();
    m_max = m.back();

    a  = std::log10(m[1] / m[0]);
    ce = (int)(-1.0 / a * std::log10(1.0 - std::pow(10.0, -a))) + 1;

    mi_m_mim1  .resize(N_m);
    mip1_m_mi  .resize(N_m);
    mip1_m_mim1.resize(N_m);
    for (int i = 0; i < N_m; ++i) {
        mi_m_mim1  [i] = m[i] * (1.0 - std::pow(10.0, -a));
        mip1_m_mi  [i] = m[i] * (std::pow(10.0,  a) - 1.0);
        mip1_m_mim1[i] = m[i] * (std::pow(10.0,  a) - std::pow(10.0, -a));
    }

    printf("[MassGrid] N_m=%d  a=%.4f  ce=%d  m_min=%.3e  m_max=%.3e\n",
           N_m, a, ce, m_min, m_max);
}

// ---------------------------------------------------------------
// CoagKernel::build  (top-level)
// ---------------------------------------------------------------
void CoagKernel::build(const MassGrid& mg)
{
    N_m = mg.N_m;
    int N2 = N_m * N_m;
    int N3 = N2 * N_m;

    // --- K (mass-split weights) ---
    std::vector<double> K(N3, 0.0);
    build_K(mg, K);

    // --- D, E (coagulation geometry) ---
    std::vector<double> D(N2, 0.0), E(N2, 0.0);
    build_DE(mg, D, E);

    // --- C_tilde then C ---
    std::vector<double> C_tilde(N3, 0.0);
    build_C_tilde(mg, K, D, E, C_tilde);

    std::vector<double> C(N3, 0.0);
    build_C_full(N_m, C_tilde, C);

    // --- sparse structures ---
    build_C_sparse_internal(mg, C);
    build_correction_sparse_internal(mg, K);
}

// ---------------------------------------------------------------
// Build K: distribute m[i]+m[j] onto the two nearest bins
// ---------------------------------------------------------------
void CoagKernel::build_K(const MassGrid& mg, std::vector<double>& K) const
{
    for (int i = 0; i < N_m; ++i) {
        for (int j = 0; j < N_m; ++j) {
            double m_coll = mg.m[i] + mg.m[j];
            int k = find_bin(mg.m, m_coll);
            if (k >= 0 && k < N_m - 1) {
                double epsilon     = (mg.m[k+1] - m_coll) / (mg.m[k+1] - mg.m[k]);
                K[idx3(i,j,k)]   = epsilon;
                K[idx3(i,j,k+1)] = 1.0 - epsilon;
            }
        }
    }
}

// ---------------------------------------------------------------
// Build D, E
// ---------------------------------------------------------------
void CoagKernel::build_DE(const MassGrid& mg,
                           std::vector<double>& D,
                           std::vector<double>& E) const
{
    int ce = mg.ce;
    for (int i = 0; i < N_m; ++i) {
        for (int j = 0; j < N_m; ++j) {
            if (i <= j + 1 - ce) {
                D[idx2(i,j)] = -mg.m[i] / mg.mip1_m_mi[j];
            } else {
                D[idx2(i,j)] = -1.0;
            }

            if (i <= j - ce) {
                E[idx2(i,j)] = mg.m[i] / mg.mi_m_mim1[j];
            } else {
                E[idx2(i,j)] = (1.0 - (mg.m[i] - mg.mi_m_mim1[j]) / mg.mip1_m_mi[j])
                              * theta_fn(mg.mip1_m_mim1[j] - mg.m[i]);
            }
        }
    }
}

// ---------------------------------------------------------------
// Build C_tilde
// ---------------------------------------------------------------
void CoagKernel::build_C_tilde(const MassGrid& mg,
                                const std::vector<double>& K,
                                const std::vector<double>& D,
                                const std::vector<double>& E,
                                std::vector<double>& C_tilde) const
{
    for (int i = 0; i < N_m; ++i) {
        for (int j = 0; j < N_m; ++j) {
            if (mg.m[i] + mg.m[j] > mg.m[N_m-1]) continue;
            for (int k = 0; k < N_m; ++k) {
                double val = 0.0;
                val += 0.5 * K[idx3(i,j,k)] * delta_fn(i,j);
                val += D[idx2(i,j)] * delta_fn(j,k);
                val += K[idx3(i,j,k)] * theta_fn(k - i - 1.5) * theta_fn(i - j - 0.5);
                if (i + 1 < N_m)
                    val += E[idx2(j,i+1)] * delta_fn(i, k-1) * theta_fn(k - j - 1.5);
                C_tilde[idx3(i,j,k)] = val;
            }
        }
    }
}

// ---------------------------------------------------------------
// Symmetrise: C[i,j,k] = C_tilde[i,j,k] + C_tilde[j,i,k]  (i>=j)
// ---------------------------------------------------------------
void CoagKernel::build_C_full(int N_m,
                               const std::vector<double>& C_tilde,
                               std::vector<double>& C) const
{
    for (int i = 0; i < N_m; ++i)
        for (int j = 0; j <= i; ++j)
            for (int k = 0; k < N_m; ++k)
                C[idx3(i,j,k)] = C_tilde[idx3(i,j,k)] + C_tilde[idx3(j,i,k)];
}

// ---------------------------------------------------------------
// Build sparse C (per output bin k, list all (i,j) pairs)
// ---------------------------------------------------------------
void CoagKernel::build_C_sparse_internal(const MassGrid& mg,
                                          const std::vector<double>& C)
{
    C_sparse_offsets.assign(N_m + 1, 0);

    for (int k = 0; k < N_m; ++k)
        for (int i = 0; i < N_m; ++i)
            for (int j = 0; j <= i; ++j)
                if (C[idx3(i,j,k)] != 0.0)
                    C_sparse_offsets[k+1]++;

    for (int k = 0; k < N_m; ++k)
        C_sparse_offsets[k+1] += C_sparse_offsets[k];

    C_sparse_total = C_sparse_offsets[N_m];
    C_sparse_i  .resize(C_sparse_total);
    C_sparse_j  .resize(C_sparse_total);
    C_sparse_val.resize(C_sparse_total);

    std::vector<int> cnt(N_m, 0);
    for (int k = 0; k < N_m; ++k) {
        for (int i = 0; i < N_m; ++i) {
            for (int j = 0; j <= i; ++j) {
                double cval = C[idx3(i,j,k)];
                if (cval == 0.0) continue;
                // Match dustcpp_mass.cpp: store C * pStick / (1+delta)
                // pStick=1 for pure coagulation
                double sym = (i == j) ? 0.5 : 1.0;
                int pos = C_sparse_offsets[k] + cnt[k]++;
                C_sparse_i  [pos] = i;
                C_sparse_j  [pos] = j;
                C_sparse_val[pos] = cval * sym;  // = C/(1+delta), same as dustcpp
            }
        }
    }
    printf("[CoagKernel] C_sparse: %d nonzeros\n", C_sparse_total);
}

// ---------------------------------------------------------------
// Build sparse correction (mass-flux limiter)
// ---------------------------------------------------------------
void CoagKernel::build_correction_sparse_internal(const MassGrid& mg,
                                                   const std::vector<double>& K)
{
    corr_offsets.assign(N_m + 1, 0);

    for (int i = 0; i < N_m; ++i) {
        for (int j = 0; j <= i; ++j) {
            double m_coll = mg.m[i] + mg.m[j];
            int k_low = find_bin(mg.m, m_coll);
            if (k_low < 0 || k_low >= N_m - 1) continue;
            corr_offsets[k_low + 1]++;
        }
    }
    for (int k = 0; k < N_m; ++k)
        corr_offsets[k+1] += corr_offsets[k];

    corr_total = corr_offsets[N_m];
    corr_pair_i.resize(corr_total);
    corr_pair_j.resize(corr_total);
    corr_eps   .resize(corr_total);
    corr_f_low .resize(corr_total);
    corr_f_high.resize(corr_total);

    std::vector<int> cnt(N_m, 0);
    for (int i = 0; i < N_m; ++i) {
        for (int j = 0; j <= i; ++j) {
            double m_coll = mg.m[i] + mg.m[j];
            int k_low = find_bin(mg.m, m_coll);
            if (k_low < 0 || k_low >= N_m - 1) continue;

            int pos = corr_offsets[k_low] + cnt[k_low]++;
            corr_pair_i[pos] = i;
            corr_pair_j[pos] = j;
            corr_eps   [pos] = K[idx3(i,j,k_low)];

            double f_low  = m_coll - mg.m[k_low];
            double f_high = m_coll - mg.m[k_low+1];
            corr_f_low [pos] = std::max(0.0, std::min(1.0,  f_low));
            corr_f_high[pos] = std::min(0.0, std::max(-1.0, f_high));
        }
    }
    printf("[CoagKernel] correction_sparse: %d nonzeros\n", corr_total);
}

// ---------------------------------------------------------------
// compute_drhodt  (pure coagulation, no fragmentation)
// collision_kernel = 1 everywhere for now
// ---------------------------------------------------------------
void CoagKernel::compute_drhodt(const MassGrid& mg,
                                 const std::vector<double>& rho,
                                 std::vector<double>& drhodt) const
{
    drhodt.assign(N_m, 0.0);

    for (int k = 0; k < N_m; ++k) {
        double val = 0.0;
        int start = C_sparse_offsets[k];
        int end   = C_sparse_offsets[k+1];
        for (int ci = start; ci < end; ++ci) {
            int i    = C_sparse_i[ci];
            int j    = C_sparse_j[ci];
            double w = C_sparse_val[ci];
            // Rstick = kernel(=1) / (1+delta), matching dustcpp_mass.cpp
            double Rstick = (i == j) ? 0.5 : 1.0;
            val += w * Rstick * rho[i] * rho[j] * mg.dm[i] * mg.dm[j] / (mg.m[i] * mg.m[j]);
        }
        drhodt[k] = val * mg.m[k] / mg.dm[k];
    }

    apply_correction(mg, rho, drhodt);
}

// ---------------------------------------------------------------
// apply_correction  (same flux-limiter as dustcpp_mass.cpp)
// ---------------------------------------------------------------
void CoagKernel::apply_correction(const MassGrid& mg,
                                   const std::vector<double>& rho,
                                   std::vector<double>& drhodt) const
{
    std::vector<double> slope_S(N_m, 0.0);

    for (int k_low = 0; k_low < N_m - 1; ++k_low) {
        int start = corr_offsets[k_low];
        int end   = corr_offsets[k_low+1];
        for (int ci = start; ci < end; ++ci) {
            int i       = corr_pair_i[ci];
            int j       = corr_pair_j[ci];
            double eps  = corr_eps[ci];
            double f_lo = corr_f_low[ci];
            double f_hi = corr_f_high[ci];

            // pStick=1, kernel=1, sym factor 1/(1+delta) already in eps via K
            double rate = rho[i] * rho[j] * mg.dm[i] * mg.dm[j] / (mg.m[i] * mg.m[j]);

            slope_S[k_low]   += eps         * rate * f_lo  / mg.dm[k_low];
            slope_S[k_low+1] += (1.0 - eps) * rate * f_hi  / mg.dm[k_low+1];
        }
    }

    std::vector<double> T_mass(N_m + 1, 0.0);
    for (int k = 1; k < N_m; ++k) {
        double s0  = slope_S[k-1];
        double s1  = slope_S[k];
        double raw = s0 + s1;
        double T_number;
        if (raw > 0.0)
            T_number = std::min(raw, std::max(0.0, s0));
        else
            T_number = std::max(raw, std::min(0.0, s1));

        T_mass[k] = (T_number >= 0.0) ? T_number * mg.m[k-1]
                                       : T_number * mg.m[k];
    }

    for (int k = 0; k < N_m; ++k)
        drhodt[k] += (T_mass[k] - T_mass[k+1]) / mg.dm[k];
}