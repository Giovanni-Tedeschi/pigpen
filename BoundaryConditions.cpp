#include "classes.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p)
{
    if (p.BC == 0) // Transmissive boundaries
    {
        c[0] = c[1];
        c[p.N_cells + 1] = c[p.N_cells];
    }
    else // Periodic boundaries
    {
        c[0] = c[p.N_cells];
        c[p.N_cells + 1] = c[1];
    }
}