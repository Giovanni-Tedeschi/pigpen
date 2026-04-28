#include "indices.h"
#include "RiemannSolvers.h"
#include "BoundaryConditions.h"
#include "Reconstruction.h"
#include <algorithm>
#include <cmath>
#include <omp.h>


void swap_xy(Cell &c)
{
    for (size_t s = 0; s < c.W.size(); s++){
        std::swap(c.W[s][idx.vx], c.W[s][idx.vy]);
        if(c.PTC==1 && s > 0) std::swap(c.W[s][idx.s11], c.W[s][idx.s22]);
    }
    for (size_t s = 0; s < c.U.size(); s++){
        std::swap(c.U[s][idx.vx], c.U[s][idx.vy]);
        if(c.PTC==1 && s > 0) std::swap(c.U[s][idx.s11], c.U[s][idx.s22]);
    }
}

void swap_xz(Cell &c)
{
    for (size_t s = 0; s < c.W.size(); s++){
        std::swap(c.W[s][idx.vx], c.W[s][idx.vz]);
        if(c.PTC==1 && s > 0){
            std::swap(c.W[s][idx.s11], c.W[s][idx.s33]);
            std::swap(c.W[s][idx.s12], c.W[s][idx.s23]);
        }
    }
    for (size_t s = 0; s < c.U.size(); s++){
        std::swap(c.U[s][idx.vx], c.U[s][idx.vz]);
        if(c.PTC==1 && s > 0){
            std::swap(c.U[s][idx.s11], c.U[s][idx.s33]);
            std::swap(c.U[s][idx.s12], c.U[s][idx.s23]);
        }
    }
}

void swap_xy_flux(Cell &c)
{
    for (size_t s = 0; s < c.FL.size(); s++)
    {
        std::swap(c.FL[s][idx.vx], c.FL[s][idx.vy]);
        if(c.PTC==1 && s > 0) std::swap(c.FL[s][idx.s11], c.FL[s][idx.s22]);
        std::swap(c.FR[s][idx.vx], c.FR[s][idx.vy]);
        if(c.PTC==1 && s > 0) std::swap(c.FR[s][idx.s11], c.FR[s][idx.s22]);
    }
}

void swap_xz_flux(Cell &c)
{
    for (size_t s = 0; s < c.FL.size(); s++)
    {
        std::swap(c.FL[s][idx.vx], c.FL[s][idx.vz]);
        if(c.PTC==1 && s > 0){
            std::swap(c.FL[s][idx.s11], c.FL[s][idx.s33]);
            std::swap(c.FL[s][idx.s12], c.FL[s][idx.s23]);
        }
        std::swap(c.FR[s][idx.vx], c.FR[s][idx.vz]);
        if(c.PTC==1 && s > 0){
            std::swap(c.FR[s][idx.s11], c.FR[s][idx.s33]);
            std::swap(c.FR[s][idx.s12], c.FR[s][idx.s23]);
        }
    }
}

void ensure_positivity(Cell &Left, Cell &Right, Params p){
    double rho_min = 1e-8;
    for (size_t s = 0; s < Left.W.size(); ++s)
    {
        // Clamp density in W and U consistently
        if (Left.W[s][idx.rho] < rho_min){
            Left.W[s][idx.rho] = rho_min;
            Left.U[s][idx.rho] = rho_min;   // U[rho] = rho directly
        }
        if (Right.W[s][idx.rho] < rho_min){
            Right.W[s][idx.rho] = rho_min;
            Right.U[s][idx.rho] = rho_min;
        }

        if (s > 0 && p.PTC == 1)
        {
            // Clamp diagonal stress in W; recompute corresponding U component
            // U[s11] = rho * s11, so clamp W then fix U
            if (Left.W[s][idx.s11] < 1e-8){
                Left.W[s][idx.s11] = 1e-8;
                Left.U[s][idx.s11] = Left.W[s][idx.rho] * 1e-8;
            }
            if (Right.W[s][idx.s11] < 1e-8){
                Right.W[s][idx.s11] = 1e-8;
                Right.U[s][idx.s11] = Right.W[s][idx.rho] * 1e-8;
            }
            if (p.N_dims >= 2){
                if (Left.W[s][idx.s22] < 1e-8){
                    Left.W[s][idx.s22] = 1e-8;
                    Left.U[s][idx.s22] = Left.W[s][idx.rho] * 1e-8;
                }
                if (Right.W[s][idx.s22] < 1e-8){
                    Right.W[s][idx.s22] = 1e-8;
                    Right.U[s][idx.s22] = Right.W[s][idx.rho] * 1e-8;
                }
            }
        }
    }
    // Do NOT call get_U_from_W() here — it would overwrite the reconstructed U
}

void compute_fluxes(std::vector<Cell> &c, Params p, const Grid& g)
{
    if (p.apply_reconstruction == 1) compute_slopes(c, p, g);

    // Pre-allocate one cL/cR pair per thread — reused every iteration, no per-call heap alloc
    // Use a representative active cell to get the right sizes
    int ref = g.flat_idx(p.N_ghost, p.N_ghost, p.N_ghost);
    int nthreads;
    #pragma omp parallel
    {
        #pragma omp single
        nthreads = omp_get_num_threads();
    }
    std::vector<Cell> cL(nthreads, c[ref]);
    std::vector<Cell> cR(nthreads, c[ref]);

    // ---- X-sweep ----
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = p.N_ghost - 1; i < g.Nx + p.N_ghost; i++)
    {
        int tid  = omp_get_thread_num();
        int flatL = g.flat_idx(i,     (p.N_dims >= 2) ? j + p.N_ghost : 0, (p.N_dims == 3) ? k + p.N_ghost : 0);
        int flatR = g.flat_idx(i + 1, (p.N_dims >= 2) ? j + p.N_ghost : 0, (p.N_dims == 3) ? k + p.N_ghost : 0);
        
        cL[tid] = c[flatL];
        cR[tid] = c[flatR];

        if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL[tid], cR[tid], 1);
        ensure_positivity(cL[tid], cR[tid], p);

        if (p.RiemannSolver == 0)
            get_exact_flux(cL[tid], cR[tid], p.GAMMA);
        else if (p.RiemannSolver == 1)
            get_hllc_flux(cL[tid], cR[tid], p.GAMMA);
        else
            get_hll_flux(cL[tid], cR[tid]);

        if(p.PTC == 0){
            get_dust_flux(cL[tid], cR[tid]);
        }else{
            get_dust_flux_PTC(cL[tid], cR[tid]);
        }

        // Write fluxes back to original cells
        c[flatL].FR = cL[tid].FR;
        c[flatR].FL = cR[tid].FL;

    }

    // ---- Y-sweep (only if N_dims >= 2) ----
    if (p.N_dims >= 2)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < g.Nz; k++)
        for (int i = 0; i < g.Nx; i++)
        for (int j = p.N_ghost - 1; j < g.Ny + p.N_ghost; j++)
        {
            int tid  = omp_get_thread_num();
            int flatL = g.flat_idx(i + p.N_ghost, j,     (p.N_dims == 3) ? k + p.N_ghost : 0);
            int flatR = g.flat_idx(i + p.N_ghost, j + 1, (p.N_dims == 3) ? k + p.N_ghost : 0);
            
            cL[tid] = c[flatL];   // work on copies
            cR[tid] = c[flatR];

            // Swap vx <-> vy on copies only
            swap_xy(cL[tid]);
            swap_xy(cR[tid]);

            // Slopes were computed in original orientation -- swap dU too
            for (size_t s = 0; s < cL[tid].dU.size(); ++s)
            for (size_t var = 0; var < cL[tid].dU[s].size(); ++var)
            {
                cL[tid].dU[s][var] = c[flatL].dUy[s][var];
                cR[tid].dU[s][var] = c[flatR].dUy[s][var];
            }

            for (size_t s = 0; s < cL[tid].dU.size(); ++s)
            {
                std::swap(cL[tid].dU[s][idx.vx], cL[tid].dU[s][idx.vy]);
                std::swap(cR[tid].dU[s][idx.vx], cR[tid].dU[s][idx.vy]);
                if(cL[tid].PTC==1 && s > 0){
                    std::swap(cL[tid].dU[s][idx.s11], cL[tid].dU[s][idx.s22]);
                    std::swap(cR[tid].dU[s][idx.s11], cR[tid].dU[s][idx.s22]);
                }
            }

            if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL[tid], cR[tid], 1);
            ensure_positivity(cL[tid], cR[tid], p);

            if (p.RiemannSolver == 0)
                get_exact_flux(cL[tid], cR[tid], p.GAMMA);
            else if (p.RiemannSolver == 1)
                get_hllc_flux(cL[tid], cR[tid], p.GAMMA);
            else
                get_hll_flux(cL[tid], cR[tid]);

            if(p.PTC == 0){
                get_dust_flux(cL[tid], cR[tid]);
            }else{
                get_dust_flux_PTC(cL[tid], cR[tid]);
            }

            // Swap fluxes back before storing
            swap_xy_flux(cL[tid]);
            swap_xy_flux(cR[tid]);

            c[flatL].FRy = cL[tid].FR;
            c[flatR].FLy = cR[tid].FL;
        }
    }

    // ---- Z-sweep (only if N_dims == 3) ----
    if (p.N_dims == 3)
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        for (int k = p.N_ghost - 1; k < g.Nz + p.N_ghost; k++)
        {
            int tid  = omp_get_thread_num();
            int flatL = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k    );
            int flatR = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + 1);

            cL[tid] = c[flatL];
            cR[tid] = c[flatR];

            swap_xz(cL[tid]);
            swap_xz(cR[tid]);

            for (size_t s = 0; s < cL[tid].dU.size(); ++s)
            for (size_t var = 0; var < cL[tid].dU[s].size(); ++var)
            {
                cL[tid].dU[s][var] = c[flatL].dUz[s][var];
                cR[tid].dU[s][var] = c[flatR].dUz[s][var];
            }
            for (size_t s = 0; s < cL[tid].dU.size(); ++s)
            {
                std::swap(cL[tid].dU[s][idx.vx], cL[tid].dU[s][idx.vz]);
                std::swap(cR[tid].dU[s][idx.vx], cR[tid].dU[s][idx.vz]);
                if(cL[tid].PTC==1 && s > 0){
                    std::swap(cL[tid].dU[s][idx.s11], cL[tid].dU[s][idx.s33]);
                    std::swap(cL[tid].dU[s][idx.s12], cL[tid].dU[s][idx.s23]);
                    std::swap(cR[tid].dU[s][idx.s11], cR[tid].dU[s][idx.s33]);
                    std::swap(cR[tid].dU[s][idx.s12], cR[tid].dU[s][idx.s23]);
                }
            }

            if (p.apply_reconstruction == 1) reconstruct_cell_pair(cL[tid], cR[tid], 1);
            ensure_positivity(cL[tid], cR[tid], p);

            if (p.RiemannSolver == 0)
                get_exact_flux(cL[tid], cR[tid], p.GAMMA);
            else if (p.RiemannSolver == 1)
                get_hllc_flux(cL[tid], cR[tid], p.GAMMA);
            else
                get_hll_flux(cL[tid], cR[tid]);

            if(p.PTC == 0){
                get_dust_flux(cL[tid], cR[tid]);
            }else{
                get_dust_flux_PTC(cL[tid], cR[tid]);
            }

            swap_xz_flux(cL[tid]);
            swap_xz_flux(cR[tid]);

            c[flatL].FRz = cL[tid].FR;
            c[flatR].FLz = cR[tid].FL;
        }
    }
}


void get_dust_flux(Cell &Left, Cell &Right)
{
    Left.get_dust_F();
    Right.get_dust_F();

    for(int j=1; j<Left.W.size(); j++){
        if ((Left.W[j][idx.vx] > 0.) && (Right.W[j][idx.vx] > 0.))
        {
            for(int l=0; l<Left.N_var_dust; l++){
                Left.FR[j][l] = Left.F[j][l];
            }
        }
        else if ((Left.W[j][idx.vx] < 0.) && (Right.W[j][idx.vx] < 0.))
        {
            for(int l=0; l<Left.N_var_dust; l++){
                Left.FR[j][l] = Right.F[j][l];
            }
        }
        else if ((Left.W[j][idx.vx] <= 0.) && (Right.W[j][idx.vx] >= 0.))
        {
            for(int l=0; l<Left.N_var_dust; l++){
                Left.FR[j][l] = 0.0;
            }
        }
        else if ((Left.W[j][idx.vx] > 0.) && (Right.W[j][idx.vx] < 0.))
        {
            for(int l=0; l<Left.N_var_dust; l++){
                Left.FR[j][l] = Left.F[j][l] + Right.F[j][l];
            }
        }

        for (int l = 0; l < Left.N_var_dust; l++)
            Right.FL[j][l] = Left.FR[j][l];
    }
}


void get_dust_flux_PTC(Cell &Left, Cell &Right)
{
    Left.get_dust_F();
    Right.get_dust_F();

    for(int j=1; j<Left.W.size(); j++){
        if(Left.W[j][idx.s11] < 1e-8) Left.W[j][idx.s11] = 1e-8;
        if(Right.W[j][idx.s11] < 1e-8) Right.W[j][idx.s11] = 1e-8;
        double lmin_left  = Left.W[j][idx.vx]  - sqrt(3.*Left.W[j][idx.s11]);
        double lmax_right = Right.W[j][idx.vx] + sqrt(3.*Right.W[j][idx.s11]);
        double inv_denom = 1./(lmin_left - lmax_right);
        double int_state = 0.0;
        for(int l=0; l<Left.N_var_dust; l++){
            if(fabs(lmin_left - lmax_right) > 1e-8){
                int_state = (lmin_left * Left.U[j][l] - lmax_right * Right.U[j][l]) * inv_denom;
                int_state -= (Left.F[j][l] - Right.F[j][l]) * inv_denom;
            }else{
                int_state = 0.0;
            }
            Left.FR[j][l] = 0.5 * (Left.F[j][l] + Right.F[j][l]) - 0.5 * fabs(lmin_left) * (int_state - Left.U[j][l]) - 0.5*fabs(lmax_right) * (Right.U[j][l] - int_state);
            Right.FL[j][l] = Left.FR[j][l];
        }
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



double get_pstar(const VarArray& WL, const VarArray& WR, double GAMMA){
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

double fK(double p,  const VarArray& WL,  const VarArray& WR, int K, double GAMMA){
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

double fprimeK(double p,  const VarArray& WL,  const VarArray& WR, int K, double GAMMA){
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