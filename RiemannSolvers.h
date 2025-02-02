#ifndef RIEMANN_SOLVERS_H
#define RIEMANN_SOLVERS_H

#include "classes.h"

void get_dust_flux(Cell& Left, Cell& Right);
void get_hll_flux(Cell& Left, Cell& Right);
void get_exact_flux(Cell& Left, Cell& Right, double GAMMA);
double fK(double p,  std::vector<double> WL,  std::vector<double> WR, int K, double GAMMA);
double fprimeK(double p,  std::vector<double> WL,  std::vector<double> WR, int K, double GAMMA);
double get_pstar( std::vector<double> WL,  std::vector<double> WR, double GAMMA);

#endif // RIEMANN_SOLVERS_H