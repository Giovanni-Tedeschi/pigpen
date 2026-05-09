#ifndef COAGULATION_H
#define COAGULATION_H

#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>

// ---------------------------------------------------------------
// MassGrid: logarithmic mass grid and associated geometry
// ---------------------------------------------------------------
struct MassGrid {
    double m_min;
    double m_max;
    int    bins_per_decade;
    int    N_m;

    double a;   // log10 spacing = 1/bins_per_decade
    int    ce;  // coagulation boundary index

    std::vector<double> m;
    std::vector<double> dm;
    std::vector<double> mi_m_mim1;    // m[i] - m[i-1]
    std::vector<double> mip1_m_mi;    // m[i+1] - m[i]
    std::vector<double> mip1_m_mim1;  // m[i+1] - m[i-1]

    MassGrid() : m_min(0), m_max(0), bins_per_decade(0), N_m(0), a(0), ce(0) {}

    void build_from_grid(const std::vector<double>& m_in,
                         const std::vector<double>& dm_in);  // build from existing m, dm
};

// ---------------------------------------------------------------
// CoagKernel: sparse coagulation matrices for drhodt
// ---------------------------------------------------------------
struct CoagKernel {
    int N_m;

    // Sparse C matrix (coagulation gain)
    std::vector<int>    C_sparse_offsets;
    std::vector<int>    C_sparse_i;
    std::vector<int>    C_sparse_j;
    std::vector<double> C_sparse_val;
    int C_sparse_total = 0;

    // Sparse correction (mass-conserving flux limiter)
    std::vector<int>    corr_offsets;
    std::vector<int>    corr_pair_i;
    std::vector<int>    corr_pair_j;
    std::vector<double> corr_eps;
    std::vector<double> corr_f_low;
    std::vector<double> corr_f_high;
    int corr_total = 0;

    CoagKernel() : N_m(0) {}

    // Build everything from a MassGrid
    void build(const MassGrid& mg);

    // Compute drhodt for a single cell's rho array
    void compute_drhodt(const MassGrid& mg,
                        const std::vector<double>& rho,
                        std::vector<double>& drhodt) const;

private:
    // Internals used during build
    static int find_bin(const std::vector<double>& m, double val);

    void build_K(const MassGrid& mg,
                 std::vector<double>& K) const;

    void build_DE(const MassGrid& mg,
                  std::vector<double>& D,
                  std::vector<double>& E) const;

    void build_C_tilde(const MassGrid& mg,
                       const std::vector<double>& K,
                       const std::vector<double>& D,
                       const std::vector<double>& E,
                       std::vector<double>& C_tilde) const;

    void build_C_full(int N_m,
                      const std::vector<double>& C_tilde,
                      std::vector<double>& C) const;

    void build_C_sparse_internal(const MassGrid& mg,
                                 const std::vector<double>& C);

    void build_correction_sparse_internal(const MassGrid& mg,
                                          const std::vector<double>& K);

    void apply_correction(const MassGrid& mg,
                          const std::vector<double>& rho,
                          std::vector<double>& drhodt) const;

    // index helpers
    inline int idx2(int i, int j) const { return i * N_m + j; }
    inline int idx3(int i, int j, int k) const { return (i * N_m + j) * N_m + k; }
};

#endif // COAGULATION_H