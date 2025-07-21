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
#include "DragIntegrators.h"
#include "BoundaryConditions.h"
#include "IO.h"


void integrate_external_force(std::vector<Cell> &c, Params p, double dt){
    for (int i = 1; i <= p.N_cells; i++)
    {
        c[i].U[0][1] += p.g0 * dt;
        c[i].get_W_from_U();
    }
}

void update_variables(std::vector<Cell> &c, Params p, double dt)
{
    for (int i = 1; i <= p.N_cells; i++)
    {
        c[i].U[0][0] += dt / p.dx * (c[i].FL[0][0] - c[i].FR[0][0]);
        c[i].U[0][1] += dt / p.dx * (c[i].FL[0][1] - c[i].FR[0][1]);
        c[i].U[0][2] += dt / p.dx * (c[i].FL[0][2] - c[i].FR[0][2]);
        for(int j = 1; j <= p.N_dust; j++){
            c[i].U[j][0] += dt / p.dx * (c[i].FL[j][0] - c[i].FR[j][0]);
            c[i].U[j][1] += dt / p.dx * (c[i].FL[j][1] - c[i].FR[j][1]);
        }
        c[i].get_W_from_U();
    }
}


void find_dt(std::vector<Cell> c, Params p, Vars &v)
{
    if(p.const_dt < 0.){
        double max_sig_vel = 1e-40;
        for (int i = 1; i <= p.N_cells; i++)
        {
            max_sig_vel = std::max(max_sig_vel, c[i].get_vsig());
        }

        v.dt = p.CFL * p.dx / max_sig_vel;
    }else{
        v.dt = p.const_dt;
    }
}

void do_integration_step(std::vector<Cell> &c, Params p, Vars v){

    if(p.DragIntegrator == 1){
        // DHD
        integrate_drag_RK(c, p, v.dt/2);

        compute_fluxes(c, p);
        update_variables(c, p, v.dt);
        integrate_external_force(c, p, v.dt);

        integrate_drag_RK(c, p, v.dt/2);

    }else if(p.DragIntegrator == 2){
        // DHDHD
        integrate_drag_RK(c, p, v.dt/4);

        compute_fluxes(c, p);
        update_variables(c, p, v.dt/2);
        integrate_external_force(c, p, v.dt/2);

        integrate_drag_RK(c, p, v.dt/2);

        compute_fluxes(c, p);
        update_variables(c, p, v.dt/2);
        integrate_external_force(c, p, v.dt/2);

        integrate_drag_RK(c, p, v.dt/4);
    }else{
        // MDIRK
        integrate_drag_MDIRK(c, p, v.dt);
    }
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
        find_dt(c, p, v);

        do_integration_step(c, p, v);

        v.t += v.dt;

        write_output(c, p, v);
    }

}
