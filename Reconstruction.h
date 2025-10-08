#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "classes.h"

void compute_slopes(std::vector<Cell> &c, Params p);
void reconstruct_cell_pair(Cell &Left, Cell &Right, int N_dust, int sign);

#endif // RECONSTRUCTION_H