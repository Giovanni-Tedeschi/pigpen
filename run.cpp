#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <sstream>
#include <cfloat>
#include <omp.h>
#include <iomanip>
#include <unordered_map>
#include "classes.h"
#include "RiemannSolvers.h"
#include "IO.h"

void apply_boundary_conditions(std::vector<Cell> &c, Params p)
{
    if (p.BC == 0)
    {
        for (int i = 0; i < p.N_vars; i++)
        {
            c[0].W[i] = c[1].W[i];
            c[0].U[i] = c[1].U[i];
            c[0].F[i] = c[1].F[i];
            c[p.N_cells + 1].W[i] = c[p.N_cells].W[i];
            c[p.N_cells + 1].U[i] = c[p.N_cells].U[i];
            c[p.N_cells + 1].F[i] = c[p.N_cells].F[i];
        }
    }
    else
    {
        for (int i = 0; i < p.N_vars; i++)
        {
            c[0].W[i] = c[p.N_cells].W[i];
            c[0].U[i] = c[p.N_cells].U[i];
            c[0].F[i] = c[p.N_cells].F[i];
            c[p.N_cells + 1].W[i] = c[1].W[i];
            c[p.N_cells + 1].U[i] = c[1].U[i];
            c[p.N_cells + 1].F[i] = c[1].F[i];
        }
    }
}

void compute_fluxes(std::vector<Cell> &c, Params p)
{
    for (int i = 1; i <= p.N_cells; i++)
    {
        c[i].get_F();
    }

    apply_boundary_conditions(c, p);

    for (int i = 0; i <= p.N_cells; i++)
    {
        if(p.RiemannSolver == 0){
            get_exact_flux(c[i], c[i + 1], p.GAMMA);
        }else{
            get_hll_flux(c[i], c[i + 1]);
        }
        get_dust_flux(c[i], c[i + 1]);
    }
}

void update_variables(std::vector<Cell> &c, Params p, Vars v)
{
    for (int i = 1; i <= p.N_cells; i++)
    {
        for (int j = 0; j < p.N_vars; j++)
            c[i].U[j] += v.dt / p.dx * (c[i].FL[j] - c[i].FR[j]);
        c[i].get_W_from_U();
    }
}

void integrate_drag(std::vector<Cell> &c, Params p, Vars v)
{
    // Find minimum stopping time
    double ts_i;
    double ts = 1e40;

    for (int i = 0; i < p.N_cells; i++)
    {
        ts_i = (c[i].W[0] * c[i].W[3]) / (p.K * (c[i].W[0] + c[i].W[3]));
        ts = std::min(ts, ts_i);
    }

    double drag;
    if(v.dt < ts){
        for (int i = 0; i < p.N_cells; i++)
        {
            drag = -p.K * v.dt * (c[i].W[4] - c[i].W[1]);
            c[i].U[1] -= drag;
            c[i].U[4] += drag;
            c[i].get_W_from_U();
        }
    }else{
        int M = 1;
        do{
          M += 1;
        }while(v.dt/M > ts);

        for(int iM=0; iM<M; iM++){
          for(int i=0; i<p.N_cells; i++){
            drag = -p.K * v.dt * (c[i].W[4] - c[i].W[1]);
            c[i].U[1] -= drag;
            c[i].U[4] += drag;
            c[i].get_W_from_U();
          }
        }
    }
}

void find_dt(std::vector<Cell> c, Params p, Vars &v)
{
    double max_sig_vel = 1e-40;
    for (int i = 1; i <= p.N_cells; i++)
    {
        max_sig_vel = std::max(max_sig_vel, c[i].get_vsig());
    }

    v.dt = p.CFL * p.dx / max_sig_vel;
    //v.dt = 0.00001;
}

int main(int argc, char *argv[])
{
    // Read the param file name as terminal input
    std::string param_file = argv[1];

    Params p = read_param(param_file);

    Vars v;

    std::vector<Cell> c = read_ic(p);

    while (v.t < p.t_max)
    {
        compute_fluxes(c, p);

        find_dt(c, p, v);

        update_variables(c, p, v);

        integrate_drag(c, p, v);

        v.t += v.dt;

        write_output(c, p, v);
    }
}
