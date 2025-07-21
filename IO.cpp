#include "IO.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <unordered_map>

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
    p.N_vars = 3 + 2 * p.N_dust;
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

    p.N_dims = configData.count("N_dims") ? stoi(configData["N_dims"]) : 1;
    p.L = configData.count("L") ? stof(configData["L"]) : 1.0;
    p.RiemannSolver = configData.count("RiemannSolver") ? stoi(configData["RiemannSolver"]) : 1;
    p.CFL = configData.count("CFL") ? stof(configData["CFL"]) : 0.1;
    p.const_dt = configData.count("dt") ? stof(configData["dt"]) : -1;
    p.g0 = configData.count("g0") ? stof(configData["g0"]) : 0.0;

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
    p.dx = p.L / p.N_cells;
    c.resize(p.N_cells + 2);

    // Resize vectors
    for (int i = 0; i <= p.N_cells + 1; i++)
    {
        c[i].GAMMA = p.GAMMA;
        c[i].N_dust = p.N_dust;
        c[i].initialize();
    }

    for (int i = 0; i < p.N_cells; i++)
    {
        c[i+1].W[0][0] = temp_W[i * p.N_vars + 0];
        c[i+1].W[0][1] = temp_W[i * p.N_vars + 1];
        c[i+1].W[0][2] = temp_W[i * p.N_vars + 2];
        
        for(int j=1; j<=p.N_dust; j++){
            c[i+1].W[j][0] = temp_W[i * p.N_vars + 3 + 2*(j-1)];
            c[i+1].W[j][1] = temp_W[i * p.N_vars + 4 + 2*(j-1)];
        }
        
        c[i + 1].get_U_from_W();
    }

    return c;
}

void write_output(std::vector<Cell> c, Params p, Vars &v)
{
    if (v.t - v.k_snap * p.dt_snap > 0)
    {
        printf("%lf %d\n", v.t, v.k_snap);
        std::string output_file = p.output_dir + std::to_string(v.k_snap) + ".txt";
        std::ofstream fp(output_file, std::ios::out);
        fp << std::scientific << std::setprecision(20); // 20 significant digits
        for (int i = 1; i <= p.N_cells; i++)
        {
            fp << c[i].W[0][0] << " " << c[i].W[0][1] << " " << c[i].W[0][2] << " ";

            for(int j = 1; j <= p.N_dust; j++){
                fp << c[i].W[j][0] << " " << c[i].W[j][1] << " ";
            }
            
            fp << v.t << " ";
            fp << "\n";
        }
        v.k_snap += 1;
        fp.close();
    }
}