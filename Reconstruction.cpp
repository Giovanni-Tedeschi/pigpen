#include "classes.h"

void compute_slopes(std::vector<Cell> &c, Params p)
{
    for (int i = 1; i <= p.N_cells; i++)
    {
        c[i].dU[0][0] = std::max((c[i+1].U[0][0] - c[i].U[0][0]) * (c[i].U[0][0] - c[i-1].U[0][0]), 0.0) / (c[i+1].U[0][0] - c[i-1].U[0][0] + 1e-20); 
        c[i].dU[0][1] = std::max((c[i+1].U[0][1] - c[i].U[0][1]) * (c[i].U[0][1] - c[i-1].U[0][1]), 0.0) / (c[i+1].U[0][1] - c[i-1].U[0][1] + 1e-20);
        c[i].dU[0][2] = std::max((c[i+1].U[0][2] - c[i].U[0][2]) * (c[i].U[0][2] - c[i-1].U[0][2]), 0.0) / (c[i+1].U[0][2] - c[i-1].U[0][2] + 1e-20);
        for(int j = 1; j <= p.N_dust; j++){
            c[i].dU[j][0] = std::max((c[i+1].U[j][0] - c[i].U[j][0]) * (c[i].U[j][0] - c[i-1].U[j][0]), 0.0) / (c[i+1].U[j][0] - c[i-1].U[j][0] + 1e-20);
            c[i].dU[j][1] = std::max((c[i+1].U[j][1] - c[i].U[j][1]) * (c[i].U[j][1] - c[i-1].U[j][1]), 0.0) / (c[i+1].U[j][1] - c[i-1].U[j][1] + 1e-20);
        }
    }
}

void reconstruct_cell_pair(Cell &Left, Cell &Right, int N_dust, int sign)
{
    Left.U[0][0] += Left.dU[0][0] * sign;
    Left.U[0][1] += Left.dU[0][1] * sign;
    Left.U[0][2] += Left.dU[0][2] * sign;
    for(int j = 1; j <= N_dust; j++){
        Left.U[j][0] += Left.dU[j][0] * sign;
        Left.U[j][1] += Left.dU[j][1] * sign;
    }

    Right.U[0][0] -= Right.dU[0][0] * sign;
    Right.U[0][1] -= Right.dU[0][1] * sign;
    Right.U[0][2] -= Right.dU[0][2] * sign;
    for(int j = 1; j <= N_dust; j++){
        Right.U[j][0] -= Right.dU[j][0] * sign;
        Right.U[j][1] -= Right.dU[j][1] * sign;
    }

    Left.get_W_from_U();
    Right.get_W_from_U();
}