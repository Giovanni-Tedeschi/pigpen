#include "classes.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p)
{
    if (p.BC == 0){ // Transmissive boundaries
        for(int i = 0; i < p.N_ghost; i++){
            c[i] = c[p.N_ghost];
            c[p.N_cells+p.N_ghost+i] = c[p.N_cells+p.N_ghost-1];
        }
    }
    else // Periodic boundaries
    {
        for(int i = 0; i < p.N_ghost; i++){
            c[p.N_ghost-1-i] = c[p.N_cells+p.N_ghost-1-i];
            c[p.N_cells + p.N_ghost+i] = c[p.N_ghost+i];
        }
    }

    if(p.Omega0 != 0){
        for(int i = 0; i < p.N_ghost; i++){
            for (size_t k = 0; k < c[0].U.size(); ++k) {
                c[p.N_ghost-1-i].U[k][idx.vy]           += p.q * p.L * p.Omega0 * c[p.N_cells+p.N_ghost-1-i].U[k][idx.rho];
                c[p.N_cells + p.N_ghost+i].U[k][idx.vy] -= p.q * p.L * p.Omega0 * c[p.N_ghost+i].U[k][idx.rho];
            }
            c[p.N_ghost-1-i].get_W_from_U();
            c[p.N_cells + p.N_ghost+i].get_W_from_U();
        }
    }
}