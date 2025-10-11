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
#include "indices.h"
#include "IO.h"


void integrate_external_force(std::vector<Cell> &c, Params p, double dt){
    for (int i = p.N_ghost; i < p.N_cells + p.N_ghost; i++)
    {
        if(p.g0 != 0.0){
            // Constant forcing term
            c[i].U[0][idx.vx] += p.g0 * dt;
        }

        if(p.Omega0 != 0.0){
            double chi0 = 0.005;
            for (size_t k = 0; k < c[i].U.size(); ++k) {
                double vx_old = c[i].Uold[k][idx.vx];

            // Centrifugal force
                c[i].U[k][idx.vx] +=  2. * p.q * pow(p.Omega0,2) * c[i].x_center * c[i].Uold[k][idx.rho] * dt;
            // Coriolis force
                c[i].U[k][idx.vx] += 2. * p.Omega0 * c[i].Uold[k][idx.vy] * dt;
                c[i].U[k][idx.vy] -= 2. * p.Omega0 * vx_old * dt;
            }
            // Pressure gradient force
            c[i].U[0][idx.vx] += chi0 * p.Omega0 * c[i].Uold[0][idx.rho] * dt;
        }
        c[i].get_W_from_U();
    }
}


void update_variables(std::vector<Cell> &c, Params p, double dt)
{
    for (int i = p.N_ghost; i < p.N_cells + p.N_ghost; i++)
    {
        for (size_t k = 0; k < c[i].U.size(); ++k) {
            for (size_t l = 0; l < c[i].U[k].size(); ++l) {
                c[i].U[k][l] += dt / p.dx * (c[i].FL[k][l] - c[i].FR[k][l]);
            }
        }
        c[i].get_W_from_U();
    }
}


void find_dt(std::vector<Cell> c, Params p, Vars &v)
{
    if(p.const_dt < 0.){
        double max_sig_vel = 1e-40;
        for (int i = p.N_ghost; i < p.N_cells + p.N_ghost; i++)
        {
            max_sig_vel = std::max(max_sig_vel, c[i].get_vsig());
        }

        v.dt = p.CFL * p.dx / max_sig_vel;
    }else{
        v.dt = p.const_dt;
    }
}

void save_old_state(std::vector<Cell> &c, Params p){
    for (int i = 0; i < p.N_cells + 2*p.N_ghost; i++)
    {
        for (size_t k = 0; k < c[i].U.size(); ++k) {
            for (size_t l = 0; l < c[i].U[k].size(); ++l) {
                c[i].Uold[k][l] = c[i].U[k][l];
            }
        }   
    }
}

void do_integration_step(std::vector<Cell> &c, Params p, Vars v){

    if(p.DragIntegrator == 1){
        // DHD
        integrate_drag_RK(c, p, v.dt/2);

        save_old_state(c, p);

        compute_fluxes(c, p);
        update_variables(c, p, v.dt);
        integrate_external_force(c, p, v.dt);
        apply_boundary_conditions(c, p);

        integrate_drag_RK(c, p, v.dt/2);

    }else if(p.DragIntegrator == 2){
        // DHDHD
        integrate_drag_RK(c, p, v.dt/4);

        save_old_state(c, p);
        compute_fluxes(c, p);
        update_variables(c, p, v.dt/2);
        integrate_external_force(c, p, v.dt/2);
        apply_boundary_conditions(c, p);

        integrate_drag_RK(c, p, v.dt/2);

        save_old_state(c, p);
        compute_fluxes(c, p);
        update_variables(c, p, v.dt/2);
        integrate_external_force(c, p, v.dt/2);
        apply_boundary_conditions(c, p);

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

    write_output(c, p, v);

    while (v.t < p.t_max)
    {
        find_dt(c, p, v);

        do_integration_step(c, p, v);

        apply_boundary_conditions(c, p);
        
        v.t += v.dt;

        write_output(c, p, v);
    }

}
