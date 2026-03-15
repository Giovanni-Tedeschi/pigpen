#include "classes.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p, const Grid& g)
{
    // ---- Y boundaries first (only if N_dims >= 2) ----
    // Must be done before X so that X ghost filling also gets correct Y ghosts
    if (p.N_dims >= 2)
    {
        for (int k = 0; k < g.Nz_tot; k++)
        for (int i = p.N_ghost; i < g.Nx + p.N_ghost; i++)  // active X only
        {
            if (p.BC == 0) // Transmissive
            {
                for (int j = 0; j < p.N_ghost; j++)
                {
                    c[g.flat_idx(i, j, k)]                = c[g.flat_idx(i, p.N_ghost, k)];
                    c[g.flat_idx(i, g.Ny_tot-1-j, k)]     = c[g.flat_idx(i, g.Ny+p.N_ghost-1, k)];
                }
            }
            else // Periodic
            {
                for (int j = 0; j < p.N_ghost; j++)
                {
                    c[g.flat_idx(i, p.N_ghost-1-j, k)]    = c[g.flat_idx(i, g.Ny+p.N_ghost-1-j, k)];
                    c[g.flat_idx(i, g.Ny+p.N_ghost+j, k)] = c[g.flat_idx(i, p.N_ghost+j, k)];
                }
            }
        }
    }

    // ---- Z boundaries (only if N_dims == 3) ----
    // Must be done before X so that X ghost filling also gets correct Z ghosts
    if (p.N_dims == 3)
    {
        for (int j = p.N_ghost; j < g.Ny + p.N_ghost; j++)  // active Y only
        for (int i = p.N_ghost; i < g.Nx + p.N_ghost; i++)  // active X only
        {
            if (p.BC == 0) // Transmissive
            {
                for (int k = 0; k < p.N_ghost; k++)
                {
                    c[g.flat_idx(i, j, k)]                = c[g.flat_idx(i, j, p.N_ghost)];
                    c[g.flat_idx(i, j, g.Nz_tot-1-k)]     = c[g.flat_idx(i, j, g.Nz+p.N_ghost-1)];
                }
            }
            else // Periodic
            {
                for (int k = 0; k < p.N_ghost; k++)
                {
                    c[g.flat_idx(i, j, p.N_ghost-1-k)]    = c[g.flat_idx(i, j, g.Nz+p.N_ghost-1-k)];
                    c[g.flat_idx(i, j, g.Nz+p.N_ghost+k)] = c[g.flat_idx(i, j, p.N_ghost+k)];
                }
            }
        }
    }

    // ---- X boundaries last (fills corners correctly) ----
    for (int k = 0; k < g.Nz_tot; k++)
    for (int j = 0; j < g.Ny_tot; j++)
    {
        if (p.BC == 0) // Transmissive
        {
            for (int i = 0; i < p.N_ghost; i++)
            {
                c[g.flat_idx(i, j, k)]              = c[g.flat_idx(p.N_ghost, j, k)];
                c[g.flat_idx(g.Nx_tot-1-i, j, k)]   = c[g.flat_idx(g.Nx+p.N_ghost-1, j, k)];
            }
        }
        else // Periodic
        {
            for (int i = 0; i < p.N_ghost; i++)
            {
                c[g.flat_idx(p.N_ghost-1-i, j, k)]  = c[g.flat_idx(g.Nx+p.N_ghost-1-i, j, k)];
                c[g.flat_idx(g.Nx+p.N_ghost+i, j, k)] = c[g.flat_idx(p.N_ghost+i, j, k)];
            }
        }

        if (p.Omega0 != 0)
        {
            for (int i = 0; i < p.N_ghost; i++)
            {
                for (size_t s = 0; s < c[0].U.size(); ++s)
                {
                    c[g.flat_idx(p.N_ghost-1-i, j, k)].U[s][idx.vy]      += p.q * p.L * p.Omega0 * c[g.flat_idx(g.Nx+p.N_ghost-1-i, j, k)].U[s][idx.rho];
                    c[g.flat_idx(g.Nx+p.N_ghost+i, j, k)].U[s][idx.vy]   -= p.q * p.L * p.Omega0 * c[g.flat_idx(p.N_ghost+i, j, k)].U[s][idx.rho];
                }
                c[g.flat_idx(p.N_ghost-1-i, j, k)].get_W_from_U();
                c[g.flat_idx(g.Nx+p.N_ghost+i, j, k)].get_W_from_U();
            }
        }
    }
}