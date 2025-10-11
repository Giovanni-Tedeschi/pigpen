#include "classes.h"

void compute_slopes(std::vector<Cell> &c, Params p)
{
    double denom;
    for (int i = p.N_ghost-1; i < p.N_cells + p.N_ghost+1; i++)
    {
        for (size_t k = 0; k < c[i].U.size(); ++k) {
            for (size_t l = 0; l < c[i].U[k].size(); ++l) {
                c[i].dU[k][l] = std::max((c[i+1].U[k][l] - c[i].U[k][l]) * (c[i].U[k][l] - c[i-1].U[k][l]), 0.0) / (c[i+1].U[k][l] - c[i-1].U[k][l] + 1e-20);
            }
        }
    }
}

void reconstruct_cell_pair(Cell &Left, Cell &Right, int sign)
{
    for (size_t k = 0; k < Left.U.size(); ++k) {
        for (size_t l = 0; l < Left.U[k].size(); ++l) {
            Left.U[k][l] += Left.dU[k][l] * sign;
            Right.U[k][l] -= Right.dU[k][l] * sign;
        }
    }

    Left.get_W_from_U();
    Right.get_W_from_U();
}