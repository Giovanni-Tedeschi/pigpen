#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "classes.h"

void compute_slopes(std::vector<Cell> &c, Params p, const Grid& g);
void reconstruct_cell_pair(Cell &Left, Cell &Right, int sign);

#endif // RECONSTRUCTION_H