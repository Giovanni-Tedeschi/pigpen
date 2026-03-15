#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "classes.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p, const Grid& g);

#endif // BOUNDARY_CONDITIONS_H