#ifndef RIEMANN_SOLVERS_H
#define RIEMANN_SOLVERS_H

#include "classes.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p);
void compute_fluxes(std::vector<Cell> &c, Params p, const Grid& g);
void get_dust_flux(Cell& Left, Cell& Right);
void get_hll_flux(Cell& Left, Cell& Right);
void get_hllc_flux(Cell &Left, Cell &Right, double GAMMA);
void get_exact_flux(Cell& Left, Cell& Right, double GAMMA);

double fK(double p,  const VarArray& WL, const VarArray& WR, int K, double GAMMA);
double fprimeK(double p,  const VarArray& WL, const VarArray& WR, int K, double GAMMA);
double get_pstar( const VarArray& WL, const VarArray& WR, double GAMMA);

#endif // RIEMANN_SOLVERS_H