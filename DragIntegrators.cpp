#include "RiemannSolvers.h"

void RK_K(Cell &c, Params p, double dt, double gamma1, double gamma2, double beta1, double beta2){
    double A1x=0., A1y=0., B1=0, C1=0, D1x=0, D1y=0, E1=0, F1=0, G1x=0, G1y=0, H1=0, L1=0, N1=0;
    double A2x=0., A2y=0., B2=0, C2=0, D2x=0, D2y=0, E2=0, F2=0, G2x=0, G2y=0, H2=0, L2=0, N2=0;
    double lambda=0., delta1 = 0., delta2=0.;
    double eps = 0.;
    int idust = 0.;
    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        c.alpha[j] = p.K[j] / c.U[idust][0];
        eps = c.U[idust][0] / c.U[0][0];

        lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));
        delta1 = 1. / (1. + gamma1 * dt * c.alpha[j]);
        delta2 = 1. / (1. + gamma2 * dt * c.alpha[j]);

        A1x += c.alpha[j] * c.U[idust][1] * delta1;
        A2x += c.alpha[j] * c.U[idust][1] * delta2;

        B1 += c.alpha[j] * eps * delta1;
        B2 += c.alpha[j] * eps * delta2;

        C1 += pow(c.alpha[j],2)*eps*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
        C2 += pow(c.alpha[j],2)*eps*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;

        D1x += pow(c.alpha[j],2)*c.U[idust][1]*(1.+c.alpha[j]*dt*(gamma1-beta2))*delta1*lambda;
        D2x += pow(c.alpha[j],2)*c.U[idust][1]*(1.+c.alpha[j]*dt*(gamma2-beta1))*delta2*lambda;

        E1 += pow(c.alpha[j],2)*eps*delta1*lambda;
        E2 += pow(c.alpha[j],2)*eps*delta2*lambda;

        F1 += pow(c.alpha[j],2)*eps * (gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2))*delta1*lambda;
        F2 += pow(c.alpha[j],2)*eps * (gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2))*delta2*lambda;
    }

    G1x = A1x - beta1*dt*D1x;
    G2x = A2x - beta2*dt*D2x;

    G1y = A1y - beta1*dt*D1y;
    G2y = A2y - beta2*dt*D2y;

    H1 = B1 - beta1*dt*C1;
    H2 = B2 - beta2*dt*C2;

    L1 = B1 - dt*F1;
    L2 = B2 - dt*F2;

    N1 = 1. + gamma1*dt*B1 - beta1*beta2*dt*dt*E1;
    N2 = 1. + gamma2*dt*B2 - beta1*beta2*dt*dt*E2;

    c.K1[0][0] = 0.;
    c.K1[0][1] = (-G1x*N2 + H1*N2*c.U[0][1] + beta1*dt*L1*(G2x-H2*c.U[0][1])) / (beta1*beta2*dt*dt*L1*L2 - N1*N2);
    c.K1[0][2] = 0.;

    c.K2[0][0] = 0.;
    c.K2[0][1] = (G2x - c.U[0][1]*H2 - c.K1[0][1]*beta2*dt*L2)/N2;
    c.K2[0][2] = 0.;

    for(int j=0; j<p.N_dust;j++){
        idust = j+1;
        eps = c.U[idust][0]/c.U[0][0];
        lambda = 1. / (1. + c.alpha[j]*dt*(gamma1 + gamma2 + c.alpha[j]*dt*(gamma1*gamma2 - beta1*beta2)));

        c.K1[idust][0] = 0.;
        c.K1[idust][1] = c.alpha[j]*lambda*((c.U[0][1]*eps - c.U[idust][1])*(1.+c.alpha[j]*dt*(gamma2-beta1)) + c.K1[0][1]*eps*dt*(gamma1+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K2[0][1]*eps*dt*beta1);
        
        c.K2[idust][0] = 0.;
        c.K2[idust][1] = c.alpha[j]*lambda*((c.U[0][1]*eps - c.U[idust][1])*(1.+c.alpha[j]*dt*(gamma1-beta2)) + c.K2[0][1]*eps*dt*(gamma2+c.alpha[j]*dt*(gamma1*gamma2-beta1*beta2)) + c.K1[0][1]*eps*dt*beta2);
    }
}

void integrate_drag_RK(std::vector<Cell> &c, Params p, double dt)
{
    double ts_i;
    double ts_max = 0.;
    for (int i = 1; i <= p.N_cells; i++)
    {
        for(int j = 1; j <= p.N_dust; j++){
            ts_i = c[i].U[j][0] / p.K[j-1];
            ts_max = std::max(ts_max, ts_i);
        }
    }

    double gamma1;
    double gamma2;
    double beta1;
    double beta2;
    double b;
    
    if(p.DragIntegrator == 1){
        // PARAMETERS FOR DHD
        if(dt <= ts_max){
            gamma1 = 1.0;
            gamma2 = 0.0;
            b = 1.;
            beta1 = 0.5 - gamma1;
            beta2 = (1. - 3.*gamma1 - 3.*gamma2 + 6.*gamma1*gamma2)/(3. - 6.*gamma1); 
            
        }else{
            gamma1 = 1.0;
            b = 0.;
            beta2 = gamma1 - 2.;
            gamma2 = 2. - gamma1;
            beta1 = (2. - 2.*gamma1 + gamma1*gamma1)/(2.-gamma1);
        }

    }else{
        // PARAMETERS FOR DHDHD      
        if(dt <= ts_max){
            gamma1 = 1.0;
            gamma2 = 0.0;
            b = 1.;
            beta1 = 0.5 - gamma1;
            beta2 = (1. - 3.*gamma1 - 3.*gamma2 + 6.*gamma1*gamma2)/(3. - 6.*gamma1); 
        }else{
            gamma1 = 1.0;
            b = 1.;
            beta1 = -1. - gamma1;
            gamma2 = 3. - gamma1;
            beta2 = (4. - 3.* gamma1 + pow(gamma1,2))/(1. + gamma1);
        }
    }

    for(int i = 0; i <= p.N_cells; i++){

        RK_K(c[i], p, dt, gamma1, gamma2, beta1, beta2);

        c[i].U[0][1] += b * dt * c[i].K1[0][1] + (1.-b) * dt * c[i].K2[0][1];
        for(int j = 1; j <= p.N_dust; j++){
            c[i].U[j][1] += b * dt * c[i].K1[j][1] + (1.-b) * dt * c[i].K2[j][1];
        }

        c[i].get_W_from_U();
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

void integrate_drag_MDIRK(std::vector<Cell> &c, Params p, double dt)
{
    compute_fluxes(c, p);

    for(int i = 1; i <= p.N_cells; i++){      
        for(int l=0; l<3; l++){
            c[i].Un[0][l] = c[i].U[0][l];    
            c[i].Ln[0][l] = (c[i].FL[0][l] - c[i].FR[0][l]) / p.dx;     
        } 
        for(int j = 1; j <= p.N_dust; j++){
            for(int l=0; l<2; l++){
                c[i].Un[j][l] = c[i].U[j][l];
                c[i].Ln[j][l] = (c[i].FL[j][l] - c[i].FR[j][l]) / p.dx;    
            }
        }
    }
    
    double ts_i;
    double ts_max = 0.;
    for (int i = 1; i <= p.N_cells; i++)
    {
        for(int j = 1; j <= p.N_dust; j++){
            ts_i = (c[i].U[0][0] * c[i].U[j][0]) / (p.K[j-1] * (c[i].U[0][0] + c[i].U[j][0]));
            ts_max = std::max(ts_max, ts_i);
        }
    }
    
    double gamma = 0.5;
    if(dt < ts_max){
        gamma = 1 - 1./sqrt(2);
    }
    
    double beta = 1. - gamma;
    double b2 = gamma;
    double delta = 1. - 1./(2.*gamma);
    
    for(int i = 1; i <= p.N_cells; i++){
        for(int l=0; l<3; l++){
            c[i].U[0][l] +=  gamma * dt * (c[i].FL[0][l] - c[i].FR[0][l]) / p.dx;    
        } 
        for(int j = 1; j <= p.N_dust; j++){
            for(int l=0; l<2; l++){
                c[i].U[j][l] +=  gamma * dt * (c[i].FL[j][l] - c[i].FR[j][l]) / p.dx;
            }
        }
            
        c[i].U[0][1] += gamma * dt * p.g0;
    
        MDIRK_K(c[i], p, dt, gamma);
    
        for(int l=0; l<3; l++){
            c[i].U[0][l] +=  gamma * dt * c[i].K[0][l];
        } 
        for(int j = 1; j <= p.N_dust; j++){
            for(int l=0; l<2; l++){
                c[i].U[j][l] += gamma * dt * c[i].K[j][l];
            }
        }
    
        c[i].get_W_from_U();
    }
    
    compute_fluxes(c, p);
    
    for(int i = 1; i <= p.N_cells; i++){
    
        for(int l=0; l<3; l++){
            c[i].U[0][l] = c[i].Un[0][l] + (1-delta)*dt*(c[i].FL[0][l] - c[i].FR[0][l])/p.dx + delta*dt*c[i].Ln[0][l] + beta * dt * c[i].K[0][l];
        } 
        for(int j = 1; j <= p.N_dust; j++){
            for(int l=0; l<2; l++){
                c[i].U[j][l] = c[i].Un[j][l] + (1-delta)*dt*(c[i].FL[j][l] - c[i].FR[j][l])/p.dx + delta*dt*c[i].Ln[j][l] + beta * dt * c[i].K[j][l];
            }
        }
    
        c[i].U[0][1] += dt * p.g0;
    
        MDIRK_K(c[i], p, dt, gamma);
    
        for(int l=0; l<3; l++){
            c[i].U[0][l] +=  b2 * dt * c[i].K[0][l];
        } 
        for(int j = 1; j <= p.N_dust; j++){
            for(int l=0; l<2; l++){
                c[i].U[j][l] += b2 * dt * c[i].K[j][l];
            }
        }
    
        c[i].get_W_from_U();
    }
}
