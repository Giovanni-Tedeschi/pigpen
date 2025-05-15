#ifndef CELL_H
#define CELL_H

#include <vector>
#include <cmath>
#include <sstream>

class Cell{
    public:
        std::vector<std::vector<double>> W;
        std::vector<std::vector<double>> U;
        std::vector<std::vector<double>> F;
        std::vector<std::vector<double>> FL;
        std::vector<std::vector<double>> FR;
        std::vector<std::vector<double>> K1;
        std::vector<std::vector<double>> K2;
        std::vector<std::vector<double>> K;
        std::vector<std::vector<double>> Ln;
        std::vector<std::vector<double>> Un;
        std::vector<double> alpha;
        double GAMMA;
        int N_dust;

        Cell& operator=(const Cell& other) {
            if (this != &other) {  // Prevent self-assignment
                // Copy all member variables
                W = other.W;
                U = other.U;
                F = other.F;
            }
            return *this;
        }

        void initialize() {
            alpha.resize(N_dust);

            W.resize(N_dust+1);
            U.resize(N_dust+1);
            F.resize(N_dust+1);
            FL.resize(N_dust+1);
            FR.resize(N_dust+1);
            K1.resize(N_dust+1);
            K2.resize(N_dust+1);
            K.resize(N_dust+1);
            Un.resize(N_dust+1);
            Ln.resize(N_dust+1);

            W[0].resize(3);
            U[0].resize(3);
            F[0].resize(3);
            FL[0].resize(3);
            FR[0].resize(3);
            K1[0].resize(3);
            K2[0].resize(3);
            K[0].resize(3);
            Un[0].resize(3);
            Ln[0].resize(3);

            for(int j=1; j<=N_dust; j++){
                W[j].resize(2);
                U[j].resize(2);   
                F[j].resize(2);
                FL[j].resize(2);
                FR[j].resize(2);
                K1[j].resize(2);
                K2[j].resize(2);
                K[j].resize(2);
                Un[j].resize(2);
                Ln[j].resize(2);
            }
        }


        void get_U_from_W(){
            // GAS
            U[0][0] = W[0][0];                                               
            U[0][1] = W[0][0] * W[0][1];        
            U[0][2] = 0.5 * W[0][0] * W[0][1] * W[0][1] + W[0][2] / (GAMMA-1.);    
            // DUST
            for(int j=1; j<=N_dust; j++){
                U[j][0] = W[j][0];
                U[j][1] = W[j][0] * W[j][1];           
                U[0][2] += 0.5 * W[j][0] * W[j][1] * W[j][1];      
            }               
        }

        void get_W_from_U(){
            // GAS
            W[0][0] = U[0][0];                    
            W[0][1] = U[0][1] / U[0][0];    
            W[0][2] = (U[0][2] - 0.5 * W[0][0]*W[0][1]*W[0][1]) * (GAMMA-1.);   
            // DUST
            for(int j=1; j<=N_dust; j++){
                W[j][0] = U[j][0];
                W[j][1] = U[j][1] / U[j][0];
                W[0][2] -= (0.5 * W[j][0] * W[j][1] * W[j][1])* (GAMMA-1.);
            }
        }

        void get_F(){
            // GAS
            F[0][0] = W[0][0]*W[0][1];
            F[0][1] = W[0][0]*W[0][1]*W[0][1] + W[0][2];
            F[0][2] = W[0][1] * (0.5 * W[0][0] * W[0][1] * W[0][1] + W[0][2] * GAMMA / (GAMMA-1.)); 
            // DUST
            for(int j=1; j<=N_dust; j++){
                F[j][0] = W[j][0]*W[j][1];
                F[j][1] = W[j][0]*W[j][1]*W[j][1];
            }
        }

        double get_SoundSpeed2(){
            return GAMMA*W[0][2]/W[0][0];
        }

        double get_vsig(){
            double vsig = fabs(W[0][1]) + sqrt(get_SoundSpeed2());
            for(int j=1; j<=N_dust; j++){
                vsig += fabs(W[j][1]);
            }
            return vsig;
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
        double const_dt;
        double BC;
        int RiemannSolver;
        int N_cells;
        int N_dims;
        int N_dust;
        int N_vars;
        int DragIntegrator;
        double g0;
        std::vector<double> K;
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