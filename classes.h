#ifndef CELL_H
#define CELL_H

#include <vector>
#include <cmath>
#include <sstream>
#include "indices.h"


class Params{
    public:
        double CFL;
        double GAMMA;
        double sound_speed;
        double t_max;
        double dt_snap;
        double L;
        double Ly;
        double Lz;
        double dx;
        double dy;
        double dz;
        double const_dt;
        double BC;
        int ghost_in_input;
        int RiemannSolver;
        int N_cells;
        int N_cells_y;
        int N_cells_z;
        int N_ghost;
        int N_dims;
        int N_dust;
        int N_var_gas;
        int N_var_dust;
        int N_vars;
        int DragIntegrator;
        int apply_reconstruction;
        double g0;
        double Omega0;
        double q;
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


struct Grid {
    int Nx, Ny, Nz;   // number of active cells per dimension
    int Ng;            // number of ghost cells
    int N_dims;

    // Total cells including ghosts per dimension
    int Nx_tot, Ny_tot, Nz_tot;

    Grid() : Nx(0), Ny(1), Nz(1), Ng(0), N_dims(1),
             Nx_tot(0), Ny_tot(1), Nz_tot(1) {}  // safe zero-init

    Grid(const Params& p) {
        Ng     = p.N_ghost;
        N_dims = p.N_dims;
        Nx     = p.N_cells;
        Ny     = (N_dims >= 2) ? p.N_cells_y : 1;
        Nz     = (N_dims == 3) ? p.N_cells_z : 1;
        Nx_tot = Nx + 2 * Ng;
        Ny_tot = (N_dims >= 2) ? Ny + 2 * Ng : 1;
        Nz_tot = (N_dims == 3) ? Nz + 2 * Ng : 1;
    }

    // Convert (i, j, k) including ghosts to flat index
    inline int flat_idx(int i, int j = 0, int k = 0) const {
        int jj = (N_dims >= 2) ? j : 0;
        int kk = (N_dims == 3) ? k : 0;
        return i + Nx_tot * (jj + Ny_tot * kk);
    }

    // Total size of the flat array
    inline int size() const {
        return Nx_tot * Ny_tot * Nz_tot;
    }

    // Check if (i,j,k) is an active (non-ghost) cell
    inline bool is_active(int i, int j = 0, int k = 0) const {
        return i >= Ng && i < Nx + Ng &&
               j >= Ng && j < Ny + Ng &&
               k >= Ng && k < Nz + Ng;
    }
};




static constexpr int MAX_VARS = 5;   // rho, vx, vy, vz, P

struct VarArray {
    double data[MAX_VARS] = {};
    int    _size          = 0;

    void resize(int n)            { _size = n; }
    int  size()  const            { return _size; }

    double&       operator[](int i)       { return data[i]; }
    const double& operator[](int i) const { return data[i]; }

    // range-for support
    double*       begin()       { return data; }
    double*       end()         { return data + _size; }
    const double* begin() const { return data; }
    const double* end()   const { return data + _size; }
};

// ── Drop-in replacement for std::vector<VarArray> ────────────────────────────
static constexpr int MAX_DUST = 8;   // max dust species supported

struct SpeciesArray {
    VarArray data[MAX_DUST + 1] = {};
    int      _size              = 0;

    void resize(int n)                  { _size = n; }
    int  size()  const                  { return _size; }

    VarArray&       operator[](int i)       { return data[i]; }
    const VarArray& operator[](int i) const { return data[i]; }
};



class Cell{
    public:
        SpeciesArray W, U, Uold, dU, dW, F, FL, FR, K1, K2, K, Ln, Un;
        VarArray alpha;
        
        /*std::vector<std::vector<double>> W;
        std::vector<std::vector<double>> U;
        std::vector<std::vector<double>> Uold;
        std::vector<std::vector<double>> dU;
        std::vector<std::vector<double>> dW;
        std::vector<std::vector<double>> F;
        std::vector<std::vector<double>> FL;
        std::vector<std::vector<double>> FR;
        std::vector<std::vector<double>> K1;
        std::vector<std::vector<double>> K2;
        std::vector<std::vector<double>> K;
        std::vector<std::vector<double>> Ln;
        std::vector<std::vector<double>> Un;
        std::vector<double> alpha;*/
        double GAMMA;
        double sound_speed;
        int N_dust;
        int N_dims;
        int N_var_gas;
        int N_var_dust;
        double x_center;
        double y_center;
        double z_center;

        /*Cell& operator=(const Cell& other) {
            if (this != &other) {  // Prevent self-assignment
                // Copy all member variables
                W = other.W;
                U = other.U;
                F = other.F;
            }
            return *this;
        }*/

        void initialize() {
            alpha.resize(N_dust);

            W.resize(N_dust+1);
            U.resize(N_dust+1);
            Uold.resize(N_dust+1);
            dW.resize(N_dust+1);
            dU.resize(N_dust+1);
            F.resize(N_dust+1);
            FL.resize(N_dust+1);
            FR.resize(N_dust+1);
            K1.resize(N_dust+1);
            K2.resize(N_dust+1);
            K.resize(N_dust+1);
            Un.resize(N_dust+1);
            Ln.resize(N_dust+1);

            W[0].resize(N_var_gas);
            U[0].resize(N_var_gas);
            Uold[0].resize(N_var_gas);
            dW[0].resize(N_var_gas);
            dU[0].resize(N_var_gas);
            F[0].resize(N_var_gas);
            FL[0].resize(N_var_gas);
            FR[0].resize(N_var_gas);
            K1[0].resize(N_var_gas);
            K2[0].resize(N_var_gas);
            K[0].resize(N_var_gas);
            Un[0].resize(N_var_gas);
            Ln[0].resize(N_var_gas);

            for(int j=1; j<=N_dust; j++){
                W[j].resize(N_var_dust);
                U[j].resize(N_var_dust);  
                Uold[j].resize(N_var_dust);  
                dW[j].resize(N_var_dust);
                dU[j].resize(N_var_dust);   
                F[j].resize(N_var_dust);
                FL[j].resize(N_var_dust);
                FR[j].resize(N_var_dust);
                K1[j].resize(N_var_dust);
                K2[j].resize(N_var_dust);
                K[j].resize(N_var_dust);
                Un[j].resize(N_var_dust);
                Ln[j].resize(N_var_dust);
            }
        }


        void get_U_from_W(){
            // GAS
            double v2 = W[0][idx.vx] * W[0][idx.vx];
            U[0][idx.rho] = W[0][idx.rho];                                               
            U[0][idx.vx] = W[0][idx.rho] * W[0][idx.vx];  

            // DUST
            for(int j=1; j<=N_dust; j++){
                U[j][idx.rho] = W[j][idx.rho];
                U[j][idx.vx] = W[j][idx.rho] * W[j][idx.vx];              
            }

            if (N_dims >= 2){
                v2 += W[0][idx.vy] * W[0][idx.vy];
                U[0][idx.vy] = W[0][idx.rho] * W[0][idx.vy];
                for(int j=1; j<=N_dust; j++){
                    U[j][idx.vy] = W[j][idx.rho] * W[j][idx.vy];              
                }
            }

            if(N_dims == 3){
                v2 += W[0][idx.vz] * W[0][idx.vz];
                U[0][idx.vz] = W[0][idx.rho] * W[0][idx.vz];
                for(int j=1; j<=N_dust; j++){
                    U[j][idx.vz] = W[j][idx.rho] * W[j][idx.vz];              
                }
            }

            if(sound_speed < 0.){
                U[0][idx.P] = 0.5 * W[0][idx.rho] * v2 + W[0][idx.P] / (GAMMA-1.);  
            }else{
                U[0][idx.P] = 0.5 * W[0][idx.rho] * v2;
            }
        }

        void get_W_from_U(){
            // GAS
            double v2 = W[0][idx.vx]*W[0][idx.vx];
            W[0][idx.rho] = U[0][idx.rho];                    
            W[0][idx.vx] = U[0][idx.vx] / U[0][idx.rho];     
            // DUST
            for(int j=1; j<=N_dust; j++){
                W[j][idx.rho] = U[j][idx.rho];
                W[j][idx.vx] = U[j][idx.vx] / U[j][idx.rho];
            }

            if (N_dims >= 2){ 
                W[0][idx.vy] = U[0][idx.vy] / U[0][idx.rho];     
                v2 += W[0][idx.vy] * W[0][idx.vy];
                for(int j=1; j<=N_dust; j++){
                    W[j][idx.vy] = U[j][idx.vy] / U[j][idx.rho];            
                }
            }

            if (N_dims == 3){ 
                W[0][idx.vz] = U[0][idx.vz] / U[0][idx.rho];     
                v2 += W[0][idx.vz] * W[0][idx.vz];
                for(int j=1; j<=N_dust; j++){
                    W[j][idx.vz] = U[j][idx.vz] / U[j][idx.rho];            
                }
            }

            if(sound_speed < 0.){
                W[0][idx.P] = (U[0][idx.P] - 0.5 * W[0][idx.rho]*v2) * (GAMMA-1.);  
            }else{
                W[0][idx.P] = W[0][idx.rho] * pow(sound_speed,2);
            }
            
        }

        void get_F(){
            double v2 = W[0][idx.vx] * W[0][idx.vx];
            F[0][idx.rho] = W[0][idx.rho]*W[0][idx.vx];
            F[0][idx.vx] = W[0][idx.rho]*W[0][idx.vx]*W[0][idx.vx] + W[0][idx.P];

            if (N_dims >= 2){
                v2 += W[0][idx.vy] * W[0][idx.vy];
                F[0][idx.vy] = W[0][idx.rho]*W[0][idx.vx]*W[0][idx.vy];
            }

            if (N_dims == 3){
                v2 += W[0][idx.vz] * W[0][idx.vz];
                F[0][idx.vz] = W[0][idx.rho]*W[0][idx.vx]*W[0][idx.vz];
            }
            
            if(sound_speed < 0.){
                F[0][idx.P] = W[0][idx.vx] * (0.5 * W[0][idx.rho] * v2 + W[0][idx.P] * GAMMA / (GAMMA-1.)); 
            }else{
                F[0][idx.P] = W[0][idx.vx] * W[0][idx.rho] * (0.5 * v2 + pow(sound_speed,2)); 
            }
            
        }

        double get_SoundSpeed2() const {
            if(sound_speed < 0.){
                return GAMMA*W[0][idx.P]/W[0][idx.rho];
            }else{
                return sound_speed*sound_speed;
            }
        }

        double get_vsig() const {
            double v2 = W[0][idx.vx] * W[0][idx.vx];
            if (N_dims >= 2) v2 += W[0][idx.vy] * W[0][idx.vy];
            if (N_dims == 3) v2 += W[0][idx.vz] * W[0][idx.vz];

            double vsig = sqrt(v2) + sqrt(get_SoundSpeed2());

            for(int j=1; j<=N_dust; j++){
                v2 = W[j][idx.vx] * W[j][idx.vx];
                if (N_dims >= 2) v2 += W[j][idx.vy] * W[j][idx.vy];
                if (N_dims == 3) v2 += W[j][idx.vz] * W[j][idx.vz];
                vsig += sqrt(v2);
            }
            return vsig;
        }
};



#endif // CELL_H