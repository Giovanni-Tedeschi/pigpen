#include "indices.h"
#include "RiemannSolvers.h"
#include "BoundaryConditions.h"
#include "Reconstruction.h"
#include <algorithm>
#include <cmath>


void swap_xy(Cell &c)
{
    for (size_t s = 0; s < c.W.size(); s++)
        std::swap(c.W[s][idx.vx], c.W[s][idx.vy]);
    for (size_t s = 0; s < c.U.size(); s++)
        std::swap(c.U[s][idx.vx], c.U[s][idx.vy]);
}

void swap_xz(Cell &c)
{
    for (size_t s = 0; s < c.W.size(); s++)
        std::swap(c.W[s][idx.vx], c.W[s][idx.vz]);
    for (size_t s = 0; s < c.U.size(); s++)
        std::swap(c.U[s][idx.vx], c.U[s][idx.vz]);
}

void swap_xy_flux(Cell &c)
{
    for (size_t s = 0; s < c.FL.size(); s++)
    {
        std::swap(c.FL[s][idx.vx], c.FL[s][idx.vy]);
        std::swap(c.FR[s][idx.vx], c.FR[s][idx.vy]);
    }
}

void swap_xz_flux(Cell &c)
{
    for (size_t s = 0; s < c.FL.size(); s++)
    {
        std::swap(c.FL[s][idx.vx], c.FL[s][idx.vz]);
        std::swap(c.FR[s][idx.vx], c.FR[s][idx.vz]);
    }
}

void compute_fluxes(std::vector<Cell> &c, Params p, const Grid& g)
{
    if (p.apply_reconstruction == 1) compute_slopes(c, p, g);

    // ---- X-sweep ----
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = p.N_ghost - 1; i < g.Nx + p.N_ghost; i++)
    {
        int flatL = g.flat_idx(i,     (p.N_dims >= 2) ? j + p.N_ghost : 0, (p.N_dims == 3) ? k + p.N_ghost : 0);
        int flatR = g.flat_idx(i + 1, (p.N_dims >= 2) ? j + p.N_ghost : 0, (p.N_dims == 3) ? k + p.N_ghost : 0);

        Cell cL = c[flatL];   // work on copies
        Cell cR = c[flatR];

        if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL, cR, 1);

        if (p.RiemannSolver == 0)
            get_exact_flux(cL, cR, p.GAMMA);
        else
            get_hll_flux(cL, cR);

        get_dust_flux(cL, cR);

        // Write fluxes back to original cells
        c[flatL].FR = cL.FR;
        c[flatR].FL = cR.FL;
    }

    // ---- Y-sweep (only if N_dims >= 2) ----
    if (p.N_dims >= 2)
    {
        for (int k = 0; k < g.Nz; k++)
        for (int j = p.N_ghost - 1; j < g.Ny + p.N_ghost; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flatL = g.flat_idx(i + p.N_ghost, j,     (p.N_dims == 3) ? k + p.N_ghost : 0);
            int flatR = g.flat_idx(i + p.N_ghost, j + 1, (p.N_dims == 3) ? k + p.N_ghost : 0);

            Cell cL = c[flatL];   // work on copies
            Cell cR = c[flatR];

            // Swap vx <-> vy on copies only
            swap_xy(cL);
            swap_xy(cR);

            // Slopes were computed in original orientation -- swap dU too
            for (size_t s = 0; s < cL.dU.size(); ++s)
                std::swap(cL.dU[s][idx.vx], cL.dU[s][idx.vy]);
            for (size_t s = 0; s < cR.dU.size(); ++s)
                std::swap(cR.dU[s][idx.vx], cR.dU[s][idx.vy]);

            if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL, cR, 1);

            if (p.RiemannSolver == 0)
                get_exact_flux(cL, cR, p.GAMMA);
            else
                get_hll_flux(cL, cR);

            get_dust_flux(cL, cR);

            // Swap fluxes back before storing
            for (size_t s = 0; s < cL.FR.size(); ++s)
                std::swap(cL.FR[s][idx.vx], cL.FR[s][idx.vy]);
            for (size_t s = 0; s < cR.FL.size(); ++s)
                std::swap(cR.FL[s][idx.vx], cR.FL[s][idx.vy]);

            c[flatL].FR = cL.FR;
            c[flatR].FL = cR.FL;
        }
    }

    // ---- Z-sweep (only if N_dims == 3) ----
    if (p.N_dims == 3)
    {
        for (int k = p.N_ghost - 1; k < g.Nz + p.N_ghost; k++)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flatL = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k    );
            int flatR = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + 1);

            Cell cL = c[flatL];
            Cell cR = c[flatR];

            swap_xz(cL);
            swap_xz(cR);

            for (size_t s = 0; s < cL.dU.size(); ++s)
                std::swap(cL.dU[s][idx.vx], cL.dU[s][idx.vz]);
            for (size_t s = 0; s < cR.dU.size(); ++s)
                std::swap(cR.dU[s][idx.vx], cR.dU[s][idx.vz]);

            if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL, cR, 1);

            if (p.RiemannSolver == 0)
                get_exact_flux(cL, cR, p.GAMMA);
            else if (p.RiemannSolver == 1)
                get_hllc_flux(cL, cR, p.GAMMA);
            else
                get_hll_flux(cL, cR);

            get_dust_flux(cL, cR);

            for (size_t s = 0; s < cL.FR.size(); ++s)
                std::swap(cL.FR[s][idx.vx], cL.FR[s][idx.vz]);
            for (size_t s = 0; s < cR.FL.size(); ++s)
                std::swap(cR.FL[s][idx.vx], cR.FL[s][idx.vz]);

            c[flatL].FR = cL.FR;
            c[flatR].FL = cR.FL;
        }
    }
}


void get_dust_flux(Cell &Left, Cell &Right)
{
    for(int j=1; j<Left.W.size(); j++){

        if ((Left.W[j][idx.vx] > 0.) && (Right.W[j][idx.vx] > 0.))
        {
            Left.FR[j][idx.rho] = Left.W[j][idx.rho] * Left.W[j][idx.vx];
            Left.FR[j][idx.vx] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vx];
            if (Left.N_dims >= 2) Left.FR[j][idx.vy] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vy];
            if (Left.N_dims == 3) Left.FR[j][idx.vz] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vz];
        }
        else if ((Left.W[j][idx.vx] < 0.) && (Right.W[j][idx.vx] < 0.))
        {
            Left.FR[j][idx.rho] = Right.W[j][idx.rho] * Right.W[j][idx.vx];
            Left.FR[j][idx.vx] = Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vx];
            if (Left.N_dims >= 2) Left.FR[j][idx.vy] = Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vy];
            if (Left.N_dims == 3) Left.FR[j][idx.vz] = Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vz];
        }
        else if ((Left.W[j][idx.vx] <= 0.) && (Right.W[j][idx.vx] >= 0.))
        {
            Left.FR[j][idx.rho] = 0.;
            Left.FR[j][idx.vx] = 0.;
            if (Left.N_dims >= 2) Left.FR[j][idx.vy] = 0.;
            if (Left.N_dims == 3) Left.FR[j][idx.vz] = 0.;
        }
        else if ((Left.W[j][idx.vx] > 0.) && (Right.W[j][idx.vx] < 0.))
        {
            Left.FR[j][idx.rho] = Left.W[j][idx.rho] * Left.W[j][idx.vx] + Right.W[j][idx.rho] * Right.W[j][idx.vx];
            Left.FR[j][idx.vx] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vx] + Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vx];
            if (Left.N_dims >= 2) Left.FR[j][idx.vy] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vy] + Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vy];
            if (Left.N_dims == 3) Left.FR[j][idx.vz] = Left.W[j][idx.rho] * Left.W[j][idx.vx] * Left.W[j][idx.vz] + Right.W[j][idx.rho] * Right.W[j][idx.vx] * Right.W[j][idx.vz];
        }

        for (int l = 0; l < Left.N_var_dust; l++)
            Right.FL[j][l] = Left.FR[j][l];
    }
}


void get_hll_flux(Cell &Left, Cell &Right)
{
    double sPlus, sL, sR;
    sPlus = std::max(fabs(Left.W[0][idx.vx]) + sqrt(Left.get_SoundSpeed2()), fabs(Right.W[0][idx.vx]) + sqrt(Right.get_SoundSpeed2()));
    sL = -sPlus;
    sR = sPlus;

    Left.get_F();
    Right.get_F();
    
    if (sL > 0)
    {
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = Left.F[0][j];
    }
    else if (sR < 0)
    {
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = Right.F[0][j];
    }
    else
    {
        for (int j = 0; j < Left.N_var_gas; j++){
            Left.FR[0][j] = (sR * Left.F[0][j] - sL * Right.F[0][j] + sL * sR * (Right.U[0][j] - Left.U[0][j])) / (sR - sL);
        }
    }

    for (int j = 0; j < Left.N_var_gas; j++)
        Right.FL[0][j] = Left.FR[0][j];
}

void get_hllc_flux(Cell &Left, Cell &Right, double GAMMA){
    double rhoL = Left.W[0][idx.rho],  rhoR = Right.W[0][idx.rho];
    double uL   = Left.W[0][idx.vx],   uR   = Right.W[0][idx.vx];
    double pL   = Left.W[0][idx.P],     pR   = Right.W[0][idx.P];
    double aL   = sqrt(Left.get_SoundSpeed2());
    double aR   = sqrt(Right.get_SoundSpeed2());

    // Transverse velocities
    double vyL = (Left.N_dims  >= 2) ? Left.W[0][idx.vy]  : 0.;
    double vyR = (Right.N_dims >= 2) ? Right.W[0][idx.vy] : 0.;
    double vzL = (Left.N_dims  == 3) ? Left.W[0][idx.vz]  : 0.;
    double vzR = (Right.N_dims == 3) ? Right.W[0][idx.vz] : 0.;

    // Energies
    double EL = pL/(GAMMA-1.) + 0.5*rhoL*(uL*uL + vyL*vyL + vzL*vzL);
    double ER = pR/(GAMMA-1.) + 0.5*rhoR*(uR*uR + vyR*vyR + vzR*vzR);

    // Pressure-based wave speed estimates (Toro et al. 1994)
    double pstar = std::max(0.5*(pL+pR) - 0.5*(uR-uL)*0.5*(rhoL+rhoR)*0.5*(aL+aR), 1e-8);

    double qL = (pstar <= pL) ? 1.0 : sqrt(1. + (GAMMA+1.)/(2.*GAMMA)*(pstar/pL - 1.));
    double qR = (pstar <= pR) ? 1.0 : sqrt(1. + (GAMMA+1.)/(2.*GAMMA)*(pstar/pR - 1.));

    double sL = uL - aL*qL;
    double sR = uR + aR*qR;

    // Contact wave speed
    double sStar = (pR - pL + rhoL*uL*(sL - uL) - rhoR*uR*(sR - uR))
                 / (rhoL*(sL - uL) - rhoR*(sR - uR));

    // Helper: builds the HLLC star state U* for one side
    // U*_K = rho_K * (s_K - u_K)/(s_K - s*) * [1, s*, vy_K, vz_K, E_K/rho_K + (s*-u_K)*(s* + p_K/(rho_K*(s_K-u_K)))]
    auto hllc_state = [&](double rhoK, double uK, double vyK, double vzK,
                          double pK, double EK, double sK)
        -> std::vector<double>
    {
        double coeff = rhoK * (sK - uK) / (sK - sStar);
        std::vector<double> Ustar(Left.N_var_gas, 0.);
        Ustar[idx.rho] = coeff;
        Ustar[idx.vx]  = coeff * sStar;
        if (Left.N_dims >= 2) Ustar[idx.vy] = coeff * vyK;
        if (Left.N_dims == 3) Ustar[idx.vz] = coeff * vzK;
        Ustar[idx.P]   = coeff * (EK/rhoK + (sStar - uK)*(sStar + pK/(rhoK*(sK - uK))));
        return Ustar;
    };

    // Compute physical fluxes F = [rho*u, rho*u^2+p, rho*u*vy, rho*u*vz, u*(E+p)]
    auto phys_flux = [&](double rhoK, double uK, double vyK, double vzK,
                         double pK, double EK)
        -> std::vector<double>
    {
        std::vector<double> F(Left.N_var_gas, 0.);
        F[idx.rho] = rhoK * uK;
        F[idx.vx]  = rhoK * uK*uK + pK;
        if (Left.N_dims >= 2) F[idx.vy] = rhoK * uK * vyK;
        if (Left.N_dims == 3) F[idx.vz] = rhoK * uK * vzK;
        F[idx.P]   = uK * (EK + pK);
        return F;
    };

    std::vector<double> FL = phys_flux(rhoL, uL, vyL, vzL, pL, EL);
    std::vector<double> FR = phys_flux(rhoR, uR, vyR, vzR, pR, ER);

    if (sL >= 0.)
    {
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = FL[j];
    }
    else if (sStar >= 0.)
    {
        std::vector<double> UstarL = hllc_state(rhoL, uL, vyL, vzL, pL, EL, sL);
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = FL[j] + sL * (UstarL[j] - Left.U[0][j]);
    }
    else if (sR >= 0.)
    {
        std::vector<double> UstarR = hllc_state(rhoR, uR, vyR, vzR, pR, ER, sR);
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = FR[j] + sR * (UstarR[j] - Right.U[0][j]);
    }
    else
    {
        for (int j = 0; j < Left.N_var_gas; j++)
            Left.FR[0][j] = FR[j];
    }

    for (int j = 0; j < Left.N_var_gas; j++)
        Right.FL[0][j] = Left.FR[0][j];
}


void get_exact_flux(Cell &Left, Cell &Right, double GAMMA){
    // CALCULATE THE GAS-PART OF THE WAVE DIAGRAM
    double Wstar[3];
    double pStar = get_pstar(Left.W[0], Right.W[0], GAMMA);
    double uStar = 0.5*(Left.W[0][1] + Right.W[0][1]) + 0.5*(fK(pStar, Left.W[0], Right.W[0], 1, GAMMA) - fK(pStar, Left.W[0], Right.W[0], -1, GAMMA));
    
    if(uStar > 0.){
        double aL = sqrt(Left.get_SoundSpeed2());
        if(pStar > Left.W[0][2]){
            // LEFT SHOCK WAVE
            double sL = Left.W[0][1] - aL*sqrt((GAMMA+1.)/2./GAMMA*pStar/Left.W[0][2] + (GAMMA-1.)/2./GAMMA);
            if(sL > 0.){
                Wstar[0] = Left.W[0][0];
                Wstar[1] = Left.W[0][1];
                Wstar[2] = Left.W[0][2];
            }else{
                Wstar[0] = Left.W[0][0]* ((pStar / Left.W[0][2] + (GAMMA-1.)/(GAMMA+1.)) / (1. + ((GAMMA-1.)/(GAMMA+1.))*pStar/Left.W[0][2]));
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }
        }else{
            // LEFT RAREFACTION WAVE
            double aStarL = aL * pow(pStar/Left.W[0][2],(GAMMA-1.)/2./GAMMA);
            double sHL = Left.W[0][1] - aL;
            double sTL = uStar - aStarL;
            if(sTL < 0.){
                Wstar[0] = Left.W[0][0] * pow(pStar/Left.W[0][2],1./GAMMA);
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }else if(sHL > 0.){
                Wstar[0] = Left.W[0][0];
                Wstar[1] = Left.W[0][1];
                Wstar[2] = Left.W[0][2];
            }else{
                Wstar[0] = Left.W[0][0] * pow(2./(GAMMA+1.) + (GAMMA-1.)/(GAMMA+1.)/aL * Left.W[0][1], 2./(GAMMA-1.));
                Wstar[1] = 2./(GAMMA+1.)* (aL + (GAMMA-1.)/2. *Left.W[0][1]);
                Wstar[2] = Left.W[0][2] * pow(2./(GAMMA+1.) + (GAMMA-1.)/(GAMMA+1.)/aL * Left.W[0][1], 2.*GAMMA/(GAMMA-1.));
            }
        }
    }else{
        double aR = sqrt(Right.get_SoundSpeed2());
        if(pStar > Right.W[0][2]){
            // RIGHT SHOCK WAVE
            double sR = Right.W[0][1] + aR*sqrt((GAMMA+1.)/2./GAMMA * pStar/Right.W[0][2] + (GAMMA-1.)/2./GAMMA);
            if(sR < 0.){
                Wstar[0] = Right.W[0][0];
                Wstar[1] = Right.W[0][1];
                Wstar[2] = Right.W[0][2];
            }else{
                Wstar[0] = Right.W[0][0] * ((pStar / Right.W[0][2] + (GAMMA-1.)/(GAMMA+1.)) / (1. + ((GAMMA-1.)/(GAMMA+1.))*pStar/Right.W[0][2]));
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }
        }else{
            // RIGHT RAREFACTION WAVE
            double aStarR = aR * pow(pStar/Right.W[0][2],(GAMMA-1.)/2./GAMMA);
            double sHR = Right.W[0][1] + aR;
            double sTR = uStar + aStarR;
            if(sHR < 0.){
                Wstar[0] = Right.W[0][0];
                Wstar[1] = Right.W[0][1];
                Wstar[2] = Right.W[0][2];
            }else if(sTR > 0.){
                Wstar[0] = Right.W[0][0] * pow(pStar/Right.W[0][2],1./GAMMA);
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }else{
                Wstar[0] = Right.W[0][0] * pow(2./(GAMMA+1.) - (GAMMA-1.)/(GAMMA+1.)/aR * Right.W[0][1], 2./(GAMMA-1.));
                Wstar[1] = 2./(GAMMA+1.)* (-aR + (GAMMA-1.)/2. * Right.W[0][1]);
                Wstar[2] = Right.W[0][2] * pow(2./(GAMMA+1.) - (GAMMA-1.)/(GAMMA+1.)/aR * Right.W[0][1], 2.*GAMMA/(GAMMA-1.));
            }
        }
    }

    // COMPUTE THE INTERMEDIATE INTERCELL FLUX FOR THE GAS PART
    double vy_star = (uStar > 0.) ? Left.W[0][idx.vy] : Right.W[0][idx.vy];
    double vz_star = (uStar > 0.) ? Left.W[0][idx.vz] : Right.W[0][idx.vz];

    Left.FR[0][idx.rho] = Wstar[0] * Wstar[1];
    Left.FR[0][idx.vx]  = Wstar[0] * pow(Wstar[1],2) + Wstar[2];
    if (Left.N_dims >= 2) Left.FR[0][idx.vy] = Wstar[0] * Wstar[1] * vy_star;
    if (Left.N_dims == 3) Left.FR[0][idx.vz] = Wstar[0] * Wstar[1] * vz_star;
    Left.FR[0][idx.P]   = Wstar[1] * (0.5*Wstar[0]*(pow(Wstar[1],2) + pow(vy_star,2) + pow(vz_star,2))
                          + Wstar[2]*GAMMA/(GAMMA-1.));

    for (int j = 0; j < Left.N_var_gas; j++)
        Right.FL[0][j] = Left.FR[0][j];
}



double get_pstar(std::vector<double> WL, std::vector<double> WR, double GAMMA){
    double pPV = 0.5*(WL[2] + WR[2]) - 1./8.*(WR[1] - WL[1])*(WL[0] + WR[0])*(sqrt(GAMMA*WL[2]/WL[0]) + sqrt(GAMMA*WR[2]/WR[0]));
    double pguess = std::max(pPV, 1e-8);
    double p = pguess;       // initialise to pguess, not 0
    double pprev = 0.;       // initialise to 0 so the loop is entered
    int iter = 0;
    while(fabs(p - pprev) / (0.5*(p + pprev)) > 1e-8){
        pprev = p;
        p = pprev - (fK(pprev, WL, WR, -1, GAMMA) + fK(pprev, WL, WR, 1, GAMMA) + WR[1] - WL[1])
                  / (fprimeK(pprev, WL, WR, -1, GAMMA) + fprimeK(pprev, WL, WR, 1, GAMMA));
        p = std::max(p, 1e-8);   // clamp to avoid negative pressure
        if (++iter > 100) break;  // safety valve
    }
    return p;
}

double fK(double p,  std::vector<double> WL,  std::vector<double> WR, int K, double GAMMA){
    double rhoK, pK, aK, AK, BK;
    if(K == 1){
        rhoK = WR[0];
        pK = WR[2];
    }else if(K == -1){
        rhoK = WL[0];
        pK = WL[2];
    }
    aK = pow(GAMMA * pK / rhoK, 0.5);
    AK = 2./(GAMMA+1.)/rhoK;
    BK = (GAMMA-1.)/(GAMMA+1.) * pK;
    if(p > pK){
        return (p-pK)*pow(AK/(p+BK),0.5);
    }else if(p <= pK){
        return 2.*aK/(GAMMA-1.)*(pow(p/pK,(GAMMA-1.)/2/GAMMA) - 1.);
    }else{
        return 0.;
    }
}

double fprimeK(double p,  std::vector<double> WL,  std::vector<double> WR, int K, double GAMMA){
    double rhoK, pK, aK, AK, BK;
    if(K == 1){
        rhoK = WR[0];
        pK = WR[2];
    }else if(K == -1){
        rhoK = WL[0];
        pK = WL[2];
    }
    aK = pow(GAMMA * pK / rhoK, 0.5);
    AK = 2./(GAMMA+1.)/rhoK;
    BK = (GAMMA-1.)/(GAMMA+1.) * pK;
    if(p > pK){
        return pow(AK/(p+BK),0.5)*(1 - (p-pK)/(2*(p+BK)));
    }else if(p <= pK){
        return pow(p/pK,-(GAMMA+1.)/2./GAMMA)/(rhoK*aK);
    }else{
        return 0.;
    }
}