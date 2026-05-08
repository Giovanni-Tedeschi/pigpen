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


void integrate_external_force(std::vector<Cell> &c, Params p, const Grid& g, double dt)
{
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        if (p.g0 != 0.0)
        {
            // Constant forcing term
            c[flat].U[0][idx.vx] += p.g0 * dt;
        }

        if (p.Omega0 != 0.0)
        {
            double chi0 = 0.005;
            for (size_t s = 0; s < c[flat].U.size(); ++s)
            {
                double vx_old = c[flat].Uold[s][idx.vx];

                // Centrifugal force
                c[flat].U[s][idx.vx] +=  2. * p.q * pow(p.Omega0, 2) * c[flat].x_center * c[flat].Uold[s][idx.rho] * dt;
                // Coriolis force
                c[flat].U[s][idx.vx] += 2. * p.Omega0 * c[flat].Uold[s][idx.vy] * dt;
                c[flat].U[s][idx.vy] -= 2. * p.Omega0 * vx_old * dt;
            }
            // Pressure gradient force
            c[flat].U[0][idx.vx] += chi0 * p.Omega0 * c[flat].Uold[0][idx.rho] * dt;
        }
        c[flat].get_W_from_U();
    }
}



void update_variables(std::vector<Cell> &c, Params p, const Grid& g, double dt)
{
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);
        int smin = (p.freeze_gas == 1) ? 1 : 0; // if freeze_gas, skip gas (s=0) when updating U
        for (size_t s = smin; s < c[flat].U.size(); ++s){
            for (size_t var = 0; var < c[flat].U[s].size(); ++var)
            {
                c[flat].U[s][var] += dt / p.dx * (c[flat].FL[s][var] - c[flat].FR[s][var]);
                if(p.N_dims >= 2) c[flat].U[s][var] += dt / p.dy * (c[flat].FLy[s][var] - c[flat].FRy[s][var]);
                if(p.N_dims == 3) c[flat].U[s][var] += dt / p.dz * (c[flat].FLz[s][var] - c[flat].FRz[s][var]);
            }
        }
        c[flat].get_W_from_U();

        for (size_t s = 0; s < c[flat].U.size(); ++s){
            if(c[flat].W[s][idx.rho] < 1e-8) c[flat].W[s][idx.rho] = 1e-8;
            if(p.PTC == 1 && s > 0){
                if(c[flat].W[s][idx.s11] < 1e-8) c[flat].W[s][idx.s11] = 1e-8;
                if(p.N_dims >= 2){
                    if(c[flat].W[s][idx.s22] < 1e-8) c[flat].W[s][idx.s22] = 1e-8;
                }
            }
        }
        c[flat].get_U_from_W();
    }
}

void find_dt(const std::vector<Cell> &c, Params p, const Grid& g, Vars &v)
{
    if (p.const_dt < 0.)
    {
        double max_sig_vel = 1e-40;

        #pragma omp parallel for collapse(3) schedule(static) reduction(max:max_sig_vel)
        for (int k = 0; k < g.Nz; k++)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);
            max_sig_vel = std::max(max_sig_vel, c[flat].get_vsig());
        }

        // CFL across all active dimensions
        double min_dx = p.dx;
        if (p.N_dims >= 2) min_dx = std::min(min_dx, p.dy);
        if (p.N_dims == 3) min_dx = std::min(min_dx, p.dz);

        v.dt = p.CFL * min_dx / max_sig_vel;
    }
    else
    {
        v.dt = p.const_dt;
    }
}

void save_old_state(std::vector<Cell> &c, const Grid& g)
{
    #pragma omp parallel for schedule(static)
    for (int flat = 0; flat < g.size(); flat++)
    for (size_t s = 0; s < c[flat].U.size(); ++s)
    for (size_t var = 0; var < c[flat].U[s].size(); ++var)
        c[flat].Uold[s][var] = c[flat].U[s][var];
}


void do_integration_step(std::vector<Cell> &c, Params p, const Grid& g, Vars v)
{
    if (p.DragIntegrator == 1)
    {
        // DHD
        integrate_drag_RK(c, p, g, v.dt/2);

        save_old_state(c, g);

        compute_fluxes(c, p, g);
        update_variables(c, p, g, v.dt);
        integrate_external_force(c, p, g, v.dt);

        integrate_drag_RK(c, p, g, v.dt/2);
    }
    else if (p.DragIntegrator == 2)
    {
        // DHDHD
        integrate_drag_RK(c, p, g, v.dt/4);

        save_old_state(c, g);
        compute_fluxes(c, p, g);
        update_variables(c, p, g, v.dt/2);
        integrate_external_force(c, p, g, v.dt/2);
        apply_boundary_conditions(c, p, g);

        integrate_drag_RK(c, p, g, v.dt/2);

        save_old_state(c, g);
        compute_fluxes(c, p, g);
        update_variables(c, p, g, v.dt/2);
        integrate_external_force(c, p, g, v.dt/2);
        apply_boundary_conditions(c, p, g);

        integrate_drag_RK(c, p, g, v.dt/4);
    }
    else
    {
        // MDIRK
        integrate_drag_MDIRK(c, p, g, v.dt);
    }
}



int main(int argc, char *argv[])
{
    // Read the param file name as terminal input
    std::string param_file = argv[1];

    Params p = read_param(param_file);

    Vars v;

    Grid g;  // will be fully built inside read_ic once N_cells is known

    std::vector<Cell> c = read_ic_hdf5(p, g);

    write_output_hdf5(c, p, v, g);

    while (v.t < p.t_max)
    {
        find_dt(c, p, g, v);

        do_integration_step(c, p, g, v);

        apply_boundary_conditions(c, p, g);

        v.t += v.dt;

        write_output_hdf5(c, p, v, g);
    }
}