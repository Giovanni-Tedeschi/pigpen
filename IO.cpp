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

    p.L  = configData.count("L")  ? stof(configData["L"])  : 1.0;
    p.Ly = configData.count("Ly") ? stof(configData["Ly"]) : p.L;
    p.Lz = configData.count("Lz") ? stof(configData["Lz"]) : p.L;
    p.N_cells   = configData.count("N_cells")   ? stoi(configData["N_cells"])   : -1;
    p.N_cells_y = configData.count("N_cells_y") ? stoi(configData["N_cells_y"]) : -1;
    p.N_cells_z = configData.count("N_cells_z") ? stoi(configData["N_cells_z"]) : -1;
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

std::vector<Cell> read_ic(Params &p, Grid& g)
{
    // Read the initial condition file
    std::ifstream file(p.input_file);

    std::vector<double> temp_W;
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double value;
        while (iss >> value)
            temp_W.push_back(value);
    }
    file.close();

    // Determine N_cells per dimension from the flat input
    int total_active_cells = temp_W.size() / p.N_vars;
    if (p.ghost_in_input)
        total_active_cells -= 2 * p.N_ghost * p.N_dims;

    if (p.N_cells < 0)
        p.N_cells = (p.N_dims == 1) ? total_active_cells
                  : (p.N_dims == 2) ? (int)round(sqrt((double)total_active_cells))
                                    : (int)round(cbrt((double)total_active_cells));

    if (p.N_cells_y < 0) p.N_cells_y = (p.N_dims >= 2) ? p.N_cells : 1;
    if (p.N_cells_z < 0) p.N_cells_z = (p.N_dims == 3) ? p.N_cells : 1;

    p.dx = p.L  / p.N_cells;
    p.dy = p.Ly / p.N_cells_y;
    p.dz = p.Lz / p.N_cells_z;

    // Rebuild the grid now that we know N_cells
    const_cast<Grid&>(g) = Grid(p);

    // Allocate flat cell array
    std::vector<Cell> c(g.size());

    // Initialize all cells (including ghosts)
    for (int k = 0; k < g.Nz_tot; k++)
    for (int j = 0; j < g.Ny_tot; j++)
    for (int i = 0; i < g.Nx_tot; i++)
    {
        int flat = g.flat_idx(i, j, k);
        c[flat].GAMMA       = p.GAMMA;
        c[flat].sound_speed = p.sound_speed;
        c[flat].N_dust      = p.N_dust;
        c[flat].N_dims      = p.N_dims;
        c[flat].N_var_gas   = p.N_var_gas;
        c[flat].N_var_dust  = p.N_var_dust;
        c[flat].x_center    = p.dx * (0.5 - p.N_ghost + i);
        c[flat].y_center    = p.dy * (0.5 - p.N_ghost + j);
        c[flat].z_center    = p.dz * (0.5 - p.N_ghost + k);
        c[flat].initialize();
    }

    // Fill active cells from temp_W
    // Input ordering: k outer -> j -> i inner (row-major, matching flat index)
    if (!p.ghost_in_input)
    {
        int input_cell = 0;
        for (int k = 0; k < g.Nz; k++)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

            for (int var = 0; var < p.N_var_gas; var++)
                c[flat].W[0][var] = temp_W[input_cell * p.N_vars + var];

            for (int s = 1; s <= p.N_dust; s++)
            for (int var = 0; var < p.N_var_dust; var++)
                c[flat].W[s][var] = temp_W[input_cell * p.N_vars + p.N_var_gas + p.N_var_dust * (s - 1) + var];

            c[flat].get_U_from_W();
            input_cell++;
        }

        apply_boundary_conditions(c, p, g);
    }
    else
    {
        int input_cell = 0;
        for (int k = 0; k < g.Nz_tot; k++)
        for (int j = 0; j < g.Ny_tot; j++)
        for (int i = 0; i < g.Nx_tot; i++)
        {
            int flat = g.flat_idx(i, j, k);

            for (int var = 0; var < p.N_var_gas; var++)
                c[flat].W[0][var] = temp_W[input_cell * p.N_vars + var];

            for (int s = 1; s <= p.N_dust; s++)
            for (int var = 0; var < p.N_var_dust; var++)
                c[flat].W[s][var] = temp_W[input_cell * p.N_vars + p.N_var_gas + p.N_var_dust * (s - 1) + var];

            c[flat].get_U_from_W();
            input_cell++;
        }
    }

    return c;
}


void write_output(const std::vector<Cell> &c, Params p, Vars &v, const Grid& g)
{
    if (v.t - v.k_snap * p.dt_snap >= 0)
    {
        printf("%lf %d\n", v.t, v.k_snap);
        std::string output_file = p.output_dir + std::to_string(v.k_snap) + ".txt";
        std::ofstream fp(output_file, std::ios::out);
        fp << std::scientific << std::setprecision(20);

        for (int k = 0; k < g.Nz; k++)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

            // Write gas variables
            for (int var = 0; var < p.N_var_gas; var++)
                fp << c[flat].W[0][var] << " ";

            // Write dust variables
            for (int s = 1; s <= p.N_dust; s++)
            for (int var = 0; var < p.N_var_dust; var++)
                fp << c[flat].W[s][var] << " ";

            // Write cell center coordinates
            fp << c[flat].x_center << " ";
            if (p.N_dims >= 2) fp << c[flat].y_center << " ";
            if (p.N_dims == 3) fp << c[flat].z_center << " ";

            fp << v.t << "\n";
        }

        v.k_snap += 1;
        fp.close();
    }
}