#ifndef CELL_H
#define CELL_H

#include <vector>
#include <cmath>
#include <sstream>

class Cell{
    public:
        std::vector<double> W;
        std::vector<double> U;
        std::vector<double> F;
        std::vector<double> FL;
        std::vector<double> FR;
        double GAMMA;

        void get_U_from_W(){
            // GAS
            U[0] = W[0];                                               
            U[1] = W[0] * W[1];        
            U[2] = 0.5 * W[0] * W[1] * W[1] + W[2] / (GAMMA-1.);    
            // DUST
            U[2] += 0.5 * W[3] * W[4] * W[4];
            U[3] = W[3];
            U[4] = W[3] * W[4];                                
        }

        void get_W_from_U(){
            // GAS
            W[0] = U[0];                    
            W[1] = U[1] / U[0];    
            W[2] = (U[2] - 0.5 * W[0]*W[1]*W[1]) * (GAMMA-1.);   
            // DUST
            W[3] = U[3];
            W[4] = U[4] / U[3];
            W[2] -= (0.5 * W[3] * W[4] * W[4])* (GAMMA-1.);
        }

        void get_F(){
            // GAS
            F[0] = W[0]*W[1];
            F[1] = W[0]*W[1]*W[1] + W[2];
            F[2] = W[1] * (0.5 * W[0] * W[1] * W[1] + W[2] * GAMMA / (GAMMA-1.)); 
            // DUST
            F[3] = W[3]*W[4];
            F[4] = W[3]*W[4]*W[4];
            F[2] += W[4] * (0.5 * W[3] * W[4] * W[4]);
        }

        double get_SoundSpeed2(){
            return GAMMA*W[2]/W[0];
        }

        double get_vsig(){
            return fabs(W[1]) + sqrt(get_SoundSpeed2()) + fabs(W[4]);
        }
};

class Params{
    public:
        double CFL;
        double GAMMA;
        double t_max;
        double dt_snap;
        double L;
        double dx;
        double BC;
        int RiemannSolver;
        int N_cells;
        int N_dims;
        int N_dust;
        int N_vars;
        double K;
        std::string input_file;
        std::string output_dir;
};

class Vars{
    public:
        double t;
        double dt;
        int k_snap;

        Vars(){
            t = 0.;
            dt = 1e-3;
            k_snap = 0;
        }
};


#endif // CELL_H