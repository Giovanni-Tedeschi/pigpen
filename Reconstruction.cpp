#include "classes.h"

// Helper lambdas at the top of compute_slopes
auto minmod = [](double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return std::abs(a) < std::abs(b) ? a : b;
};

// van Leer limiter
auto van_leer = [](double a, double b) -> double {
    if (a * b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
};

// MC (monotonized central) limiter
auto MC = [](double a, double b) -> double {
    if (a * b <= 0.0) return 0.0;
    double c = 0.5 * (a + b);
    return std::abs(c) < 2.0*std::abs(a) && std::abs(c) < 2.0*std::abs(b)
           ? c : (std::abs(a) < std::abs(b) ? 2.0*a : 2.0*b);
};

void compute_slopes(std::vector<Cell> &c, Params p, const Grid& g)
{
    // X-slopes (always)
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = p.N_ghost - 1; i < g.Nx + p.N_ghost + 1; i++)
    {
        int flat  = g.flat_idx(i,     j + p.N_ghost, k + p.N_ghost);
        int flatL = g.flat_idx(i - 1, j + p.N_ghost, k + p.N_ghost);
        int flatR = g.flat_idx(i + 1, j + p.N_ghost, k + p.N_ghost);

        for (size_t s = 0; s < c[flat].U.size(); ++s)
        for (size_t var = 0; var < c[flat].U[s].size(); ++var)
        {
            double dR = c[flatR].U[s][var] - c[flat].U[s][var];
            double dL = c[flat].U[s][var]  - c[flatL].U[s][var];
            c[flat].dU[s][var] = 0.5 * van_leer(dL, dR);   // swap to van_leer or minmod as needed
        }
    }
    // Y-slopes (only if N_dims >= 2)
    if (p.N_dims >= 2)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < g.Nz; k++)
        for (int i = 0; i < g.Nx; i++)
        for (int j = p.N_ghost - 1; j < g.Ny + p.N_ghost + 1; j++)
        {
            int flat  = g.flat_idx(i + p.N_ghost, j,     k + p.N_ghost);
            int flatL = g.flat_idx(i + p.N_ghost, j - 1, k + p.N_ghost);
            int flatR = g.flat_idx(i + p.N_ghost, j + 1, k + p.N_ghost);

            for (size_t s = 0; s < c[flat].U.size(); ++s)
            for (size_t var = 0; var < c[flat].U[s].size(); ++var)
            {
                double dR = c[flatR].U[s][var] - c[flat].U[s][var];
                double dL = c[flat].U[s][var]  - c[flatL].U[s][var];
                c[flat].dUy[s][var] = 0.5 * van_leer(dL, dR);   // swap to van_leer or minmod as needed
            }
        }
    }

    // Z-slopes (only if N_dims == 3)
    if (p.N_dims == 3)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        for (int k = p.N_ghost - 1; k < g.Nz + p.N_ghost + 1; k++)
        {
            int flat  = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k    );
            int flatL = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k - 1);
            int flatR = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + 1);

            for (size_t s = 0; s < c[flat].U.size(); ++s)
            for (size_t var = 0; var < c[flat].U[s].size(); ++var)
            {
                double dR = c[flatR].U[s][var] - c[flat].U[s][var];
                double dL = c[flat].U[s][var]  - c[flatL].U[s][var];
                c[flat].dUz[s][var] = 0.5 * van_leer(dL, dR);   // swap to van_leer or minmod as needed
            }
        }
    }
}


void reconstruct_cell_pair(Cell &Left, Cell &Right, int sign)
{
    for (size_t s = 0; s < Left.U.size(); ++s)
    for (size_t var = 0; var < Left.U[s].size(); ++var)
    {
        Left.U[s][var]  += Left.dU[s][var]  * sign;
        Right.U[s][var] -= Right.dU[s][var] * sign;
    }

    Left.get_W_from_U();
    Right.get_W_from_U();
}