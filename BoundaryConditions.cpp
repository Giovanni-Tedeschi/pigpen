#include "classes.h"

// BC codes
static constexpr int BC_TRANSMISSIVE = 0;
static constexpr int BC_PERIODIC     = 1;
static constexpr int BC_REFLECTIVE   = 2;

static inline int resolve_bc(int face_bc, int fallback_bc)
{
    return (face_bc >= BC_TRANSMISSIVE && face_bc <= BC_REFLECTIVE) ? face_bc : fallback_bc;
}

static inline void reflect_normal_momentum(Cell& cell, int normal_dim, const Params& p)
{
    // Flip normal momentum for gas + dust in conserved variables U
    for (int s = 0; s <= p.N_dust; ++s)
    {
        if (normal_dim == 0) cell.U[s][idx.vx] *= -1.0;
        if (normal_dim == 1 && p.N_dims >= 2) cell.U[s][idx.vy] *= -1.0;
        if (normal_dim == 2 && p.N_dims == 3) cell.U[s][idx.vz] *= -1.0;
    }

    // If PTC is active, flip mixed stress terms with one normal index
    if (p.PTC == 1)
    {
        for (int s = 1; s <= p.N_dust; ++s)
        {
            if (normal_dim == 0)
            {
                if (p.N_dims >= 2) cell.U[s][idx.s12] *= -1.0;
                if (p.N_dims == 3) cell.U[s][idx.s13] *= -1.0;
            }
            else if (normal_dim == 1 && p.N_dims >= 2)
            {
                cell.U[s][idx.s12] *= -1.0;
                if (p.N_dims == 3) cell.U[s][idx.s23] *= -1.0;
            }
            else if (normal_dim == 2 && p.N_dims == 3)
            {
                cell.U[s][idx.s13] *= -1.0;
                cell.U[s][idx.s23] *= -1.0;
            }
        }
    }

    cell.get_W_from_U();
}

void apply_boundary_conditions(std::vector<Cell> &c, Params p, const Grid& g)
{
    const int bc_fallback = static_cast<int>(p.BC);

    const int bc_xmin = resolve_bc(p.BC_xmin, bc_fallback);
    const int bc_xmax = resolve_bc(p.BC_xmax, bc_fallback);
    const int bc_ymin = resolve_bc(p.BC_ymin, bc_fallback);
    const int bc_ymax = resolve_bc(p.BC_ymax, bc_fallback);
    const int bc_zmin = resolve_bc(p.BC_zmin, bc_fallback);
    const int bc_zmax = resolve_bc(p.BC_zmax, bc_fallback);

    // ---- Y boundaries first (if 2D/3D) ----
    if (p.N_dims >= 2)
    {
        for (int k = 0; k < g.Nz_tot; ++k)
        for (int i = p.N_ghost; i < g.Nx + p.N_ghost; ++i) // active X only
        {
            for (int gg = 0; gg < p.N_ghost; ++gg)
            {
                const int jL_g = p.N_ghost - 1 - gg;
                const int jR_g = g.Ny + p.N_ghost + gg;

                // ymin
                if (bc_ymin == BC_TRANSMISSIVE)
                {
                    c[g.flat_idx(i, jL_g, k)] = c[g.flat_idx(i, p.N_ghost, k)];
                }
                else if (bc_ymin == BC_PERIODIC)
                {
                    c[g.flat_idx(i, jL_g, k)] = c[g.flat_idx(i, g.Ny + p.N_ghost - 1 - gg, k)];
                }
                else // reflective
                {
                    c[g.flat_idx(i, jL_g, k)] = c[g.flat_idx(i, p.N_ghost + gg, k)];
                    reflect_normal_momentum(c[g.flat_idx(i, jL_g, k)], 1, p);
                }

                // ymax
                if (bc_ymax == BC_TRANSMISSIVE)
                {
                    c[g.flat_idx(i, jR_g, k)] = c[g.flat_idx(i, g.Ny + p.N_ghost - 1, k)];
                }
                else if (bc_ymax == BC_PERIODIC)
                {
                    c[g.flat_idx(i, jR_g, k)] = c[g.flat_idx(i, p.N_ghost + gg, k)];
                }
                else // reflective
                {
                    c[g.flat_idx(i, jR_g, k)] = c[g.flat_idx(i, g.Ny + p.N_ghost - 1 - gg, k)];
                    reflect_normal_momentum(c[g.flat_idx(i, jR_g, k)], 1, p);
                }
            }
        }
    }

    // ---- Z boundaries second (if 3D) ----
    if (p.N_dims == 3)
    {
        for (int j = p.N_ghost; j < g.Ny + p.N_ghost; ++j) // active Y only
        for (int i = p.N_ghost; i < g.Nx + p.N_ghost; ++i) // active X only
        {
            for (int gg = 0; gg < p.N_ghost; ++gg)
            {
                const int kL_g = p.N_ghost - 1 - gg;
                const int kR_g = g.Nz + p.N_ghost + gg;

                // zmin
                if (bc_zmin == BC_TRANSMISSIVE)
                {
                    c[g.flat_idx(i, j, kL_g)] = c[g.flat_idx(i, j, p.N_ghost)];
                }
                else if (bc_zmin == BC_PERIODIC)
                {
                    c[g.flat_idx(i, j, kL_g)] = c[g.flat_idx(i, j, g.Nz + p.N_ghost - 1 - gg)];
                }
                else // reflective
                {
                    c[g.flat_idx(i, j, kL_g)] = c[g.flat_idx(i, j, p.N_ghost + gg)];
                    reflect_normal_momentum(c[g.flat_idx(i, j, kL_g)], 2, p);
                }

                // zmax
                if (bc_zmax == BC_TRANSMISSIVE)
                {
                    c[g.flat_idx(i, j, kR_g)] = c[g.flat_idx(i, j, g.Nz + p.N_ghost - 1)];
                }
                else if (bc_zmax == BC_PERIODIC)
                {
                    c[g.flat_idx(i, j, kR_g)] = c[g.flat_idx(i, j, p.N_ghost + gg)];
                }
                else // reflective
                {
                    c[g.flat_idx(i, j, kR_g)] = c[g.flat_idx(i, j, g.Nz + p.N_ghost - 1 - gg)];
                    reflect_normal_momentum(c[g.flat_idx(i, j, kR_g)], 2, p);
                }
            }
        }
    }

    // ---- X boundaries last (fills corners) ----
    for (int k = 0; k < g.Nz_tot; ++k)
    for (int j = 0; j < g.Ny_tot; ++j)
    {
        for (int gg = 0; gg < p.N_ghost; ++gg)
        {
            const int iL_g = p.N_ghost - 1 - gg;
            const int iR_g = g.Nx + p.N_ghost + gg;

            // xmin
            if (bc_xmin == BC_TRANSMISSIVE)
            {
                c[g.flat_idx(iL_g, j, k)] = c[g.flat_idx(p.N_ghost, j, k)];
            }
            else if (bc_xmin == BC_PERIODIC)
            {
                c[g.flat_idx(iL_g, j, k)] = c[g.flat_idx(g.Nx + p.N_ghost - 1 - gg, j, k)];
            }
            else // reflective
            {
                c[g.flat_idx(iL_g, j, k)] = c[g.flat_idx(p.N_ghost + gg, j, k)];
                reflect_normal_momentum(c[g.flat_idx(iL_g, j, k)], 0, p);
            }

            // xmax
            if (bc_xmax == BC_TRANSMISSIVE)
            {
                c[g.flat_idx(iR_g, j, k)] = c[g.flat_idx(g.Nx + p.N_ghost - 1, j, k)];
            }
            else if (bc_xmax == BC_PERIODIC)
            {
                c[g.flat_idx(iR_g, j, k)] = c[g.flat_idx(p.N_ghost + gg, j, k)];
            }
            else // reflective
            {
                c[g.flat_idx(iR_g, j, k)] = c[g.flat_idx(g.Nx + p.N_ghost - 1 - gg, j, k)];
                reflect_normal_momentum(c[g.flat_idx(iR_g, j, k)], 0, p);
            }
        }

        // Shearing-box correction only meaningful with periodic X/X
        if (p.Omega0 != 0 && bc_xmin == BC_PERIODIC && bc_xmax == BC_PERIODIC)
        {
            for (int i = 0; i < p.N_ghost; i++)
            {
                for (size_t s = 0; s < c[0].U.size(); ++s)
                {
                    c[g.flat_idx(p.N_ghost-1-i, j, k)].U[s][idx.vy]    += p.q * p.L * p.Omega0 * c[g.flat_idx(g.Nx+p.N_ghost-1-i, j, k)].U[s][idx.rho];
                    c[g.flat_idx(g.Nx+p.N_ghost+i, j, k)].U[s][idx.vy] -= p.q * p.L * p.Omega0 * c[g.flat_idx(p.N_ghost+i, j, k)].U[s][idx.rho];
                }
                c[g.flat_idx(p.N_ghost-1-i, j, k)].get_W_from_U();
                c[g.flat_idx(g.Nx+p.N_ghost+i, j, k)].get_W_from_U();
            }
        }
    }

    // INJECT DUST (unchanged)
    for (int k = p.N_ghost; k < g.Nz + p.N_ghost; ++k)
    for (int j = p.N_ghost; j < g.Ny + p.N_ghost; ++j)
    {
        double inject_y1 = 0.5;
        double inject_y2 = 1.5;
        double inject_half_width = 0.11;
        double inject_rho = 1.0;
        double inject_vx = 0.2;
        double inject_vy = 0.0;

        const int ref = g.flat_idx(p.N_ghost, j, k);
        const double y = c[ref].y_center;

        const bool in_jet =
            (std::abs(y - inject_y1) < inject_half_width) ||
            (std::abs(y - inject_y2) < inject_half_width);

        for (int i = 0; i < p.N_ghost; ++i)
            {
                const int flat = g.flat_idx(i, j, k);
                for (int s = 1; s <= p.N_dust; ++s)
                {
                    if (in_jet){
                        c[flat].U[s][idx.rho] = inject_rho;
                        c[flat].U[s][idx.vx]  = inject_rho * inject_vx;
                        if (p.N_dims >= 2) c[flat].U[s][idx.vy] = inject_rho * inject_vy;

                        if (p.PTC == 1)
                        {
                            c[flat].U[s][idx.s11] = inject_rho * (1e-8 + pow(inject_vx, 2));
                            if (p.N_dims >= 2)
                            {
                                c[flat].U[s][idx.s12] = inject_rho * (inject_vx * inject_vy);
                                c[flat].U[s][idx.s22] = inject_rho * (1e-8 + pow(inject_vy, 2));
                            }
                        }
                    }else{
                        c[flat].U[s][idx.rho] = 1e-8;
                        c[flat].U[s][idx.vx]  = 0.0;
                        c[flat].U[s][idx.vy]  = 0.0;
                        if (p.PTC == 1)
                        {
                            c[flat].U[s][idx.s11] = 1e-8 * 1e-8;
                            c[flat].U[s][idx.s12] = 0.0;
                            c[flat].U[s][idx.s22] = 1e-8 * 1e-8;
                        }
                    }
                }
                
                c[flat].get_W_from_U();
            }
    }
}