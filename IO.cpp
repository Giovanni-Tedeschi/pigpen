#include "IO.h"
#include "indices.h"
#include "BoundaryConditions.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <unordered_map>

Indices idx; // Define global variable

Params read_param(std::string fname)
{
    Params p;

    // Read the param file
    std::unordered_map<std::string, std::string> configData;
    std::ifstream configFileStream(fname);
    for (std::string line{}; std::getline(configFileStream, line);)
    {
        // Remove comments starting with '#'
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Trim whitespace from the line
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        std::istringstream iss{line};
        if (std::string id{}, value{}; std::getline(std::getline(iss, id, ':'), value))
        {
            configData[id] = value;
        }
    }

    p.GAMMA = stof(configData["GAMMA"]);
    p.t_max = stof(configData["t_max"]);
    p.dt_snap = stof(configData["dt_snap"]);
    p.BC = stoi(configData["BoundaryCondition"]);
    p.DragIntegrator = stoi(configData["DragIntegrator"]);
    p.N_dust = stoi(configData["DustSpecies"]);
    p.N_dims = configData.count("N_dims") ? stoi(configData["N_dims"]) : 1;
    p.N_ghost = configData.count("N_ghost") ? stoi(configData["N_ghost"]) : 2;
    p.ghost_in_input = configData.count("ghost_in_input") ? stoi(configData["ghost_in_input"]) : 0;
    p.N_var_gas = 2 + p.N_dims;
    p.N_var_dust = 1 + p.N_dims;
    p.N_vars = p.N_var_gas + p.N_dust * p.N_var_dust;
    p.K.resize(p.N_dust);
    for(int j = 1; j<=p.N_dust; j++){
        p.K[j-1] = stof(configData["K_" + std::to_string(j)]);
    }

    p.input_file = configData["input_file"];
    p.input_file.erase(std::remove_if(p.input_file.begin(), p.input_file.end(), [](char c)
                                      { return c == '\n' || c == '\r'; }),
                       p.input_file.end());
    p.output_dir = configData["output_dir"];
    p.output_dir.erase(std::remove_if(p.output_dir.begin(), p.output_dir.end(), [](char c)
                                      { return c == '\n' || c == '\r'; }),
                       p.output_dir.end());

    p.L = configData.count("L") ? stof(configData["L"]) : 1.0;
    p.sound_speed = configData.count("sound_speed") ? stof(configData["sound_speed"]) : -1.0;
    p.RiemannSolver = configData.count("RiemannSolver") ? stoi(configData["RiemannSolver"]) : 1;
    p.apply_reconstruction = configData.count("apply_reconstruction") ? stoi(configData["apply_reconstruction"]) : 1;
    p.CFL = configData.count("CFL") ? stof(configData["CFL"]) : 0.1;
    p.const_dt = configData.count("dt") ? stof(configData["dt"]) : -1;
    p.g0 = configData.count("g0") ? stof(configData["g0"]) : 0.0;
    p.Omega0 = configData.count("Omega0") ? stof(configData["Omega0"]) : 0.0;
    p.q = configData.count("q") ? stof(configData["q"]) : 0.0;

    idx.rho = 0;
    idx.vx  = 1;
    if (p.N_dims >= 2) idx.vy = 2;
    if (p.N_dims == 3) idx.vz = 3;
    idx.P   = 1 + p.N_dims;

    return p;
}

std::vector<Cell> read_ic(Params &p)
{
    std::vector<Cell> c;

    // Read the initial condition file
    std::ifstream file;
    file.open(p.input_file);

    std::string line;
    std::vector<double> temp_W;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double value;
        while (iss >> value)
        {
            temp_W.push_back(value);
        }
    }

    file.close();

    // Determine N_cells
    p.N_cells = temp_W.size() / p.N_vars;
    if(p.ghost_in_input){
        p.N_cells -= 2*p.N_ghost;
    }
    p.dx = p.L / p.N_cells;
    c.resize(p.N_cells + 2*p.N_ghost);

    for (int i = 0; i < p.N_cells + 2*p.N_ghost; i++)
    {
        c[i].GAMMA = p.GAMMA;
        c[i].sound_speed = p.sound_speed;
        c[i].N_dust = p.N_dust;
        c[i].N_dims = p.N_dims;
        c[i].N_var_gas = p.N_var_gas;
        c[i].N_var_dust = p.N_var_dust;
        c[i].x_center = p.dx * (0.5 - p.N_ghost + i);
        c[i].initialize();
    }


    if(!p.ghost_in_input){
        for (int i = 0; i < p.N_cells; i++)
        {
            for(int k=0; k<p.N_var_gas; k++){
                c[i+p.N_ghost].W[0][k] = temp_W[i * p.N_vars + k];
            }
            
            for(int j=1; j<=p.N_dust; j++){
                for(int k=0; k<p.N_var_dust; k++){
                    c[i+p.N_ghost].W[j][k] = temp_W[i * p.N_vars + p.N_var_gas + p.N_var_dust*(j-1) + k];
                }
            }
            
            c[i + p.N_ghost].get_U_from_W();
        }

        apply_boundary_conditions(c, p);
        
    }else{
        for (int i = 0; i < p.N_cells + 2*p.N_ghost; i++)
        {
            for(int k=0; k<p.N_var_gas; k++){
                c[i].W[0][k] = temp_W[i * p.N_vars + k];
            }
            
            for(int j=1; j<=p.N_dust; j++){
                for(int k=0; k<p.N_var_dust; k++){
                    c[i].W[j][k] = temp_W[i * p.N_vars + p.N_var_gas + p.N_var_dust*(j-1) + k];
                }
            }
            
            c[i].get_U_from_W();
        }
    }

    return c;
}

void write_output(std::vector<Cell> c, Params p, Vars &v)
{
    if (v.t - v.k_snap * p.dt_snap >= 0)
    {
        printf("%lf %d\n", v.t, v.k_snap);
        std::string output_file = p.output_dir + std::to_string(v.k_snap) + ".txt";
        std::ofstream fp(output_file, std::ios::out);
        fp << std::scientific << std::setprecision(20); // 20 significant digits
        
        for (int i = p.N_ghost; i < p.N_cells + p.N_ghost; i++)
        {
            for(int k=0; k<p.N_var_gas; k++){
                fp << c[i].W[0][k] << " ";
            }

            for(int j = 1; j <= p.N_dust; j++){
                for(int k=0; k<p.N_var_dust; k++){
                    fp << c[i].W[j][k] << " ";    
                }
            }
            
            fp << v.t << " ";
            fp << "\n";
        }
        v.k_snap += 1;
        fp.close();
    }
}