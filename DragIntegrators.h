#ifndef DRAG_INTEGRATORS_H
#define DRAG_INTEGRATORS_H

#include "classes.h"

void RK_K(Cell &c, Params p, double dt, double gamma1, double gamma2, double beta1, double beta2);
void integrate_drag_RK(std::vector<Cell> &c, Params p, double dt);

void MDIRK_K(Cell &c, Params p, double dt, double gamma);
void integrate_drag_MDIRK(std::vector<Cell> &c, Params p, double dt);

#endif // DRAG_INTEGRATORS_H