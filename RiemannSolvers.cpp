#include "RiemannSolvers.h"
#include <algorithm>
#include <cmath>

void get_dust_flux(Cell &Left, Cell &Right)
{
    double rho_dL = Left.W[3];
    double vel_dL = Left.W[4];
    double rho_dR = Right.W[3];
    double vel_dR = Right.W[4];
    if ((vel_dL >= 0.) && (vel_dR >= 0.))
    {
        Left.FR[3] = rho_dL * vel_dL;
        Left.FR[4] = rho_dL * pow(vel_dL, 2);
    }
    else if ((vel_dL <= 0.) && (vel_dR <= 0.))
    {
        Left.FR[3] = rho_dR * vel_dR;
        Left.FR[4] = rho_dR * pow(vel_dR, 2);
    }
    else if ((vel_dL <= 0.) && (vel_dR >= 0.))
    {
        Left.FR[3] = 0.;
        Left.FR[4] = 0.;
    }
    else if ((vel_dL >= 0.) && (vel_dR <= 0.))
    {
        Left.FR[3] = rho_dL * vel_dL + rho_dR * vel_dR;
        Left.FR[4] = rho_dL * pow(vel_dL, 2) + rho_dR * pow(vel_dR, 2);
    }
    for (int j = 3; j < 5; j++)
        Right.FL[j] = Left.FR[j];
}


void get_hll_flux(Cell &Left, Cell &Right)
{
    double sPlus, sL, sR;
    sPlus = std::max(fabs(Left.W[1]) + sqrt(Left.get_SoundSpeed2()), fabs(Right.W[1]) + sqrt(Right.get_SoundSpeed2()));
    sL = -sPlus;
    sR = sPlus;

    if (sL >= 0)
    {
        for (int j = 0; j < 3; j++)
            Left.FR[j] = Left.F[j];
    }
    else if (sR <= 0)
    {
        for (int j = 0; j < 3; j++)
            Left.FR[j] = Right.F[j];
    }
    else
    {
        for (int j = 0; j < 3; j++)
            Left.FR[j] = (sR * Left.F[j] - sL * Right.F[j] + sL * sR * (Right.U[j] - Left.U[j])) / (sR - sL);
    }

    for (int j = 0; j < 3; j++)
        Right.FL[j] = Left.FR[j];
}


void get_exact_flux(Cell &Left, Cell &Right, double GAMMA){
    // CALCULATE THE GAS-PART OF THE WAVE DIAGRAM
    double Wstar[3];
    double pStar = get_pstar(Left.W, Right.W, GAMMA);
    double uStar = 0.5*(Left.W[1] + Right.W[1]) + 0.5*(fK(pStar, Left.W, Right.W, 1, GAMMA) - fK(pStar, Left.W, Right.W, -1, GAMMA));
    
    if(uStar > 0.){
        double aL = sqrt(Left.get_SoundSpeed2());
        if(pStar > Left.W[2]){
            // LEFT SHOCK WAVE
            double sL = Left.W[1] - aL*sqrt((GAMMA+1.)/2./GAMMA*pStar/Left.W[2] + (GAMMA-1.)/2./GAMMA);
            if(sL > 0.){
                Wstar[0] = Left.W[0];
                Wstar[1] = Left.W[1];
                Wstar[2] = Left.W[2];
            }else{
                Wstar[0] = Left.W[0]* ((pStar / Left.W[2] + (GAMMA-1.)/(GAMMA+1.)) / (1. + ((GAMMA-1.)/(GAMMA+1.))*pStar/Left.W[2]));
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }
        }else{
            // LEFT RAREFACTION WAVE
            double aStarL = aL * pow(pStar/Left.W[2],(GAMMA-1.)/2./GAMMA);
            double sHL = Left.W[1] - aL;
            double sTL = uStar - aStarL;
            if(sTL < 0.){
                Wstar[0] = Left.W[0] * pow(pStar/Left.W[2],1./GAMMA);
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }else if(sHL > 0.){
                Wstar[0] = Left.W[0];
                Wstar[1] = Left.W[1];
                Wstar[2] = Left.W[2];
            }else{
                Wstar[0] = Left.W[0] * pow(2./(GAMMA+1.) + (GAMMA-1.)/(GAMMA+1.)/aL * Left.W[1], 2./(GAMMA-1.));
                Wstar[1] = 2./(GAMMA+1.)* (aL + (GAMMA-1.)/2. *Left.W[1]);
                Wstar[2] = Left.W[2] * pow(2./(GAMMA+1.) + (GAMMA-1.)/(GAMMA+1.)/aL * Left.W[1], 2.*GAMMA/(GAMMA-1.));
            }
        }
    }else{
        double aR = sqrt(Right.get_SoundSpeed2());
        if(pStar > Right.W[2]){
            // RIGHT SHOCK WAVE
            double sR = Right.W[1] + aR*sqrt((GAMMA+1.)/2./GAMMA * pStar/Right.W[2] + (GAMMA-1.)/2./GAMMA);
            if(sR < 0.){
                Wstar[0] = Right.W[0];
                Wstar[1] = Right.W[1];
                Wstar[2] = Right.W[2];
            }else{
                Wstar[0] = Right.W[0] * ((pStar / Right.W[2] + (GAMMA-1.)/(GAMMA+1.)) / (1. + ((GAMMA-1.)/(GAMMA+1.))*pStar/Right.W[2]));
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }
        }else{
            // RIGHT RAREFACTION WAVE
            double aStarR = aR * pow(pStar/Right.W[2],(GAMMA-1.)/2./GAMMA);
            double sHR = Right.W[1] + aR;
            double sTR = uStar + aStarR;
            if(sHR < 0.){
                Wstar[0] = Right.W[0];
                Wstar[1] = Right.W[1];
                Wstar[2] = Right.W[2];
            }else if(sTR > 0.){
                Wstar[0] = Right.W[0] * pow(pStar/Right.W[2],1./GAMMA);
                Wstar[1] = uStar;
                Wstar[2] = pStar;
            }else{
                Wstar[0] = Right.W[0] * pow(2./(GAMMA+1.) - (GAMMA-1.)/(GAMMA+1.)/aR * Right.W[1], 2./(GAMMA-1.));
                Wstar[1] = 2./(GAMMA+1.)* (-aR + (GAMMA-1.)/2. * Right.W[1]);
                Wstar[2] = Right.W[2] * pow(2./(GAMMA+1.) - (GAMMA-1.)/(GAMMA+1.)/aR * Right.W[1], 2.*GAMMA/(GAMMA-1.));
            }
        }
    }

    // COMPUTE THE INTERMEDIATE INTERCELL FLUX FOR THE GAS PART
    Left.FR[0] = Wstar[0] * Wstar[1];
    Left.FR[1] = Wstar[0] * pow(Wstar[1],2) + Wstar[2];
    Left.FR[2] = 0.5 * Wstar[0] * pow(Wstar[1],3) + Wstar[1]*Wstar[2] + Wstar[1]*Wstar[2]/(GAMMA-1.);
    for (int j = 0; j < 3; j++)
        Right.FL[j] = Left.FR[j];
}



double get_pstar( std::vector<double> WL,  std::vector<double> WR, double GAMMA){
    double pPV = 0.5*(WL[2] + WR[2]) - 1./8.*(WR[1] - WL[1])*(WL[0] + WR[0])*(sqrt(GAMMA*WL[2]/WL[0]) + sqrt(GAMMA*WR[2]/WR[0]));
    double pguess = std::max(pPV,1e-8);
    double p = 0.;
    double pprev = pguess;
    while(fabs(p-pprev) / (0.5*(p+pprev)) > 1e-8){
        pprev = pguess;
        p = pguess - (fK(pguess, WL, WR,-1, GAMMA) + fK(pguess, WL, WR,1, GAMMA) + WR[1] - WL[1]) / (fprimeK(pguess, WL, WR, -1, GAMMA) + fprimeK(pguess, WL, WR, 1, GAMMA));
        pguess = p;
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