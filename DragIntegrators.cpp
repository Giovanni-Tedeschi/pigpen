#include "RiemannSolvers.h"
#include "indices.h"

void RK_K(Cell &c, Params p, double dt, double gamma1, double gamma2, double beta1, double beta2){
    double A1x=0., A1y=0., A1z=0., B1=0, C1=0, D1x=0, D1y=0, D1z=0, E1=0, F1=0, G1x=0, G1y=0, G1z=0, H1=0, L1=0, N1=0;
    double A2x=0., A2y=0., A2z=0., B2=0, C2=0, D2x=0, D2y=0, D2z=0, E2=0, F2=0, G2x=0, G2y=0, G2z=0, H2=0, L2=0, N2=0;
    double lambda=0., delta1 = 0., delta2=0.;
    double eps = 0.;
    int idust = 0.;
    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        c.alpha[j] = p.K[j] / c.U[idust][idx.rho];
        eps = c.U[idust][idx.rho] / c.U[0][idx.rho];

        lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));
        delta1 = 1. / (1. + gamma1 * dt * c.alpha[j]);
        delta2 = 1. / (1. + gamma2 * dt * c.alpha[j]);

        A1x += c.alpha[j] * c.U[idust][idx.vx] * delta1;
        A2x += c.alpha[j] * c.U[idust][idx.vx] * delta2;
        
        B1 += c.alpha[j] * eps * delta1;
        B2 += c.alpha[j] * eps * delta2;

        C1 += pow(c.alpha[j],2)*eps*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
        C2 += pow(c.alpha[j],2)*eps*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;

        D1x += pow(c.alpha[j],2)*c.U[idust][idx.vx]*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
        D2x += pow(c.alpha[j],2)*c.U[idust][idx.vx]*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;
        
        E1 += pow(c.alpha[j],2)*eps*delta1*lambda;
        E2 += pow(c.alpha[j],2)*eps*delta2*lambda;

        F1 += pow(c.alpha[j],2)*eps * (gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2))*delta1*lambda;
        F2 += pow(c.alpha[j],2)*eps * (gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2))*delta2*lambda;

        if(p.N_dims >= 2){
            A1y += c.alpha[j] * c.U[idust][idx.vy] * delta1;
            A2y += c.alpha[j] * c.U[idust][idx.vy] * delta2;
            D1y += pow(c.alpha[j],2)*c.U[idust][idx.vy]*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
            D2y += pow(c.alpha[j],2)*c.U[idust][idx.vy]*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;
        }

        if(p.N_dims == 3){
            A1z += c.alpha[j] * c.U[idust][idx.vz] * delta1;
            A2z += c.alpha[j] * c.U[idust][idx.vz] * delta2;
            D1z += pow(c.alpha[j],2)*c.U[idust][idx.vz]*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
            D2z += pow(c.alpha[j],2)*c.U[idust][idx.vz]*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;
        }
    }

    G1x = A1x - beta1*dt*D1x;
    G2x = A2x - beta2*dt*D2x;

    G1y = A1y - beta1*dt*D1y;
    G2y = A2y - beta2*dt*D2y;

    G1z = A1z - beta1*dt*D1z;
    G2z = A2z - beta2*dt*D2z;

    H1 = B1 - beta1*dt*C1;
    H2 = B2 - beta2*dt*C2;

    L1 = B1 - dt*F1;
    L2 = B2 - dt*F2;

    N1 = 1. + gamma1*dt*B1 - beta1*beta2*dt*dt*E1;
    N2 = 1. + gamma2*dt*B2 - beta1*beta2*dt*dt*E2;

    c.K1[0][idx.rho] = 0.;
    c.K1[0][idx.vx] = (-G1x*N2 + H1*N2*c.U[0][idx.vx] + beta1*dt*L1*(G2x-H2*c.U[0][idx.vx])) / (beta1*beta2*dt*dt*L1*L2 - N1*N2);
    c.K1[0][idx.P] = 0.;

    c.K2[0][idx.rho] = 0.;
    c.K2[0][idx.vx] = (G2x - c.U[0][idx.vx]*H2 - c.K1[0][idx.vx]*beta2*dt*L2)/N2;
    c.K2[0][idx.P] = 0.;

    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        eps = c.U[idust][idx.rho]/c.U[0][idx.rho];
        lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));

        c.K1[idust][idx.rho] = 0.;
        c.K1[idust][idx.vx] = c.alpha[j]*lambda*((c.U[0][idx.vx]*eps - c.U[idust][idx.vx])*(1.+c.alpha[j]*dt*(gamma2-beta1)) + c.K1[0][idx.vx]*eps*dt*(gamma1+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K2[0][idx.vx]*eps*dt*beta1);
        
        c.K2[idust][idx.rho] = 0.;
        c.K2[idust][idx.vx] = c.alpha[j]*lambda*((c.U[0][idx.vx]*eps - c.U[idust][idx.vx])*(1.+c.alpha[j]*dt*(gamma1-beta2)) + c.K2[0][idx.vx]*eps*dt*(gamma2+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K1[0][idx.vx]*eps*dt*beta2);
    }
    
    if(p.N_dims >= 2){
        c.K1[0][idx.vy] = (-G1y*N2 + H1*N2*c.U[0][idx.vy] + beta1*dt*L1*(G2y-H2*c.U[0][idx.vy])) / (beta1*beta2*dt*dt*L1*L2 - N1*N2);
        c.K2[0][idx.vy] = (G2y - c.U[0][idx.vy]*H2 - c.K1[0][idx.vy]*beta2*dt*L2)/N2;
        for(int j=0; j<p.N_dust;j++){
            idust = j+1;
            eps = c.U[idust][idx.rho]/c.U[0][idx.rho];
            lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));
            c.K1[idust][idx.vy] = c.alpha[j]*lambda*((c.U[0][idx.vy]*eps - c.U[idust][idx.vy])*(1.+c.alpha[j]*dt*(gamma2-beta1)) + c.K1[0][idx.vy]*eps*dt*(gamma1+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K2[0][idx.vy]*eps*dt*beta1);
            c.K2[idust][idx.vy] = c.alpha[j]*lambda*((c.U[0][idx.vy]*eps - c.U[idust][idx.vy])*(1.+c.alpha[j]*dt*(gamma1-beta2)) + c.K2[0][idx.vy]*eps*dt*(gamma2+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K1[0][idx.vy]*eps*dt*beta2);
        }
    }

    if(p.N_dims == 3){
        c.K1[0][idx.vz] = (-G1z*N2 + H1*N2*c.U[0][idx.vz] + beta1*dt*L1*(G2z-H2*c.U[0][idx.vz])) / (beta1*beta2*dt*dt*L1*L2 - N1*N2);
        c.K2[0][idx.vz] = (G2z - c.U[0][idx.vz]*H2 - c.K1[0][idx.vz]*beta2*dt*L2)/N2;
        for(int j=0; j<p.N_dust;j++){
            idust = j+1;
            eps = c.U[idust][idx.rho]/c.U[0][idx.rho];
            lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));
            c.K1[idust][idx.vz] = c.alpha[j]*lambda*((c.U[0][idx.vz]*eps - c.U[idust][idx.vz])*(1.+c.alpha[j]*dt*(gamma2-beta1)) + c.K1[0][idx.vz]*eps*dt*(gamma1+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K2[0][idx.vz]*eps*dt*beta1);
            c.K2[idust][idx.vz] = c.alpha[j]*lambda*((c.U[0][idx.vz]*eps - c.U[idust][idx.vz])*(1.+c.alpha[j]*dt*(gamma1-beta2)) + c.K2[0][idx.vz]*eps*dt*(gamma2+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K1[0][idx.vz]*eps*dt*beta2);
        }
    }
}

void integrate_drag_RK(std::vector<Cell> &c, Params p, const Grid& g, double dt)
{
    // Find max stopping time across all active cells
    double ts_max = 0.;

    #pragma omp parallel for collapse(3) schedule(static) reduction(max:ts_max)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);
        for (int s = 1; s <= p.N_dust; s++)
        {
            double ts_i = c[flat].U[s][idx.rho] / p.K[s-1];
            ts_max = std::max(ts_max, ts_i);
        }
    }

    double gamma1, gamma2, beta1, beta2, b;

    if (p.DragIntegrator == 1)
    {
        // PARAMETERS FOR DHD
        if (dt <= ts_max)
        {
            gamma1 = 1.0;
            gamma2 = 0.0;
            b = 1.;
            beta1 = 0.5 - gamma1;
            beta2 = (1. - 3.*gamma1 - 3.*gamma2 + 6.*gamma1*gamma2) / (3. - 6.*gamma1);
        }
        else
        {
            gamma1 = 1.0;
            b = 0.;
            beta2 = gamma1 - 2.;
            gamma2 = 2. - gamma1;
            beta1 = (2. - 2.*gamma1 + gamma1*gamma1) / (2. - gamma1);
        }
    }
    else
    {
        // PARAMETERS FOR DHDHD
        if (dt <= ts_max)
        {
            gamma1 = 1.0;
            gamma2 = 0.0;
            b = 1.;
            beta1 = 0.5 - gamma1;
            beta2 = (1. - 3.*gamma1 - 3.*gamma2 + 6.*gamma1*gamma2) / (3. - 6.*gamma1);
        }
        else
        {
            gamma1 = 1.0;
            b = 1.;
            beta1 = -1. - gamma1;
            gamma2 = 3. - gamma1;
            beta2 = (4. - 3.*gamma1 + pow(gamma1, 2)) / (1. + gamma1);
        }
    }

    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        RK_K(c[flat], p, dt, gamma1, gamma2, beta1, beta2);

        c[flat].U[0][idx.vx] += b * dt * c[flat].K1[0][idx.vx] + (1.-b) * dt * c[flat].K2[0][idx.vx];
        for (int s = 1; s <= p.N_dust; s++)
            c[flat].U[s][idx.vx] += b * dt * c[flat].K1[s][idx.vx] + (1.-b) * dt * c[flat].K2[s][idx.vx];

        if (p.N_dims >= 2)
        {
            c[flat].U[0][idx.vy] += b * dt * c[flat].K1[0][idx.vy] + (1.-b) * dt * c[flat].K2[0][idx.vy];
            for (int s = 1; s <= p.N_dust; s++)
                c[flat].U[s][idx.vy] += b * dt * c[flat].K1[s][idx.vy] + (1.-b) * dt * c[flat].K2[s][idx.vy];
        }

        if (p.N_dims == 3)
        {
            c[flat].U[0][idx.vz] += b * dt * c[flat].K1[0][idx.vz] + (1.-b) * dt * c[flat].K2[0][idx.vz];
            for (int s = 1; s <= p.N_dust; s++)
                c[flat].U[s][idx.vz] += b * dt * c[flat].K1[s][idx.vz] + (1.-b) * dt * c[flat].K2[s][idx.vz];
        }

        c[flat].get_W_from_U();
    }
}

void MDIRK_K(Cell &c, Params p, double dt, double gamma){
    double A = 0.;
    double B = 0.;
    int idust = 0.;
    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        c.alpha[j] = p.K[j] / c.U[idust][0];
        A += (c.U[idust][1] * c.alpha[j]) / (1. + gamma * dt * c.alpha[j]);
        B += c.U[idust][0] / c.U[0][0] * c.alpha[j] / (1. + gamma * dt * c.alpha[j]);
    }
    
    c.K[0][0] = 0.;
    c.K[0][1] = (A - c.U[0][1]*B) / (1. + gamma * dt * B);
    c.K[0][2] = 0.;
    
    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        c.K[idust][0] = 0.;
        c.K[idust][1] = c.alpha[j] / (1.+gamma*dt * c.alpha[j]) * ((c.U[idust][0]/c.U[0][0] * c.U[0][1] - c.U[idust][1])  +   gamma*dt*c.U[idust][0]/c.U[0][0] * c.K[0][1]);
    }
}


void integrate_drag_MDIRK(std::vector<Cell> &c, Params p, const Grid& g, double dt)
{
    compute_fluxes(c, p, g);

    // Save state and compute L operator for all active cells
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        for (int var = 0; var < p.N_var_gas; var++)
        {
            c[flat].Un[0][var] = c[flat].U[0][var];
            c[flat].Ln[0][var] = (c[flat].FL[0][var] - c[flat].FR[0][var]) / p.dx;
        }
        for (int s = 1; s <= p.N_dust; s++)
        for (int var = 0; var < p.N_var_dust; var++)
        {
            c[flat].Un[s][var] = c[flat].U[s][var];
            c[flat].Ln[s][var] = (c[flat].FL[s][var] - c[flat].FR[s][var]) / p.dx;
        }
    }

    // Find max stopping time
    double ts_max = 0.;

    #pragma omp parallel for collapse(3) schedule(static) reduction(max:ts_max)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);
        for (int s = 1; s <= p.N_dust; s++)
        {
            double ts_i = (c[flat].U[0][idx.rho] * c[flat].U[s][idx.rho])
                 / (p.K[s-1] * (c[flat].U[0][idx.rho] + c[flat].U[s][idx.rho]));
            ts_max = std::max(ts_max, ts_i);
        }
    }

    double gamma = 0.5;
    if (dt < ts_max)
        gamma = 1. - 1./sqrt(2.);

    double beta  = 1. - gamma;
    double b2    = gamma;
    double delta = 1. - 1./(2.*gamma);

    // First MDIRK stage
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        for (int var = 0; var < p.N_var_gas; var++)
            c[flat].U[0][var] += gamma * dt * (c[flat].FL[0][var] - c[flat].FR[0][var]) / p.dx;
        for (int s = 1; s <= p.N_dust; s++)
        for (int var = 0; var < p.N_var_dust; var++)
            c[flat].U[s][var] += gamma * dt * (c[flat].FL[s][var] - c[flat].FR[s][var]) / p.dx;

        c[flat].U[0][idx.vx] += gamma * dt * p.g0;

        MDIRK_K(c[flat], p, dt, gamma);

        for (int var = 0; var < p.N_var_gas; var++)
            c[flat].U[0][var] += gamma * dt * c[flat].K[0][var];
        for (int s = 1; s <= p.N_dust; s++)
        for (int var = 0; var < p.N_var_dust; var++)
            c[flat].U[s][var] += gamma * dt * c[flat].K[s][var];

        c[flat].get_W_from_U();
    }

    compute_fluxes(c, p, g);

    // Second MDIRK stage
    #pragma omp parallel for collapse(3) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        for (int var = 0; var < p.N_var_gas; var++)
            c[flat].U[0][var] = c[flat].Un[0][var]
                              + (1.-delta) * dt * (c[flat].FL[0][var] - c[flat].FR[0][var]) / p.dx
                              + delta * dt * c[flat].Ln[0][var]
                              + beta * dt * c[flat].K[0][var];
        for (int s = 1; s <= p.N_dust; s++)
        for (int var = 0; var < p.N_var_dust; var++)
            c[flat].U[s][var] = c[flat].Un[s][var]
                              + (1.-delta) * dt * (c[flat].FL[s][var] - c[flat].FR[s][var]) / p.dx
                              + delta * dt * c[flat].Ln[s][var]
                              + beta * dt * c[flat].K[s][var];

        c[flat].U[0][idx.vx] += dt * p.g0;

        MDIRK_K(c[flat], p, dt, gamma);

        for (int var = 0; var < p.N_var_gas; var++)
            c[flat].U[0][var] += b2 * dt * c[flat].K[0][var];
        for (int s = 1; s <= p.N_dust; s++)
        for (int var = 0; var < p.N_var_dust; var++)
            c[flat].U[s][var] += b2 * dt * c[flat].K[s][var];

        c[flat].get_W_from_U();
    }
}