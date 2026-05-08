#include "IO.h"
#include "indices.h"
#include "BoundaryConditions.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <H5Cpp.h>

Indices idx; // Define global variable

namespace {
int parse_bc_code(const std::string& raw, int fallback)
{
    if (raw.empty()) return fallback;

    std::string s = raw;
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    if (s == "0" || s == "transmissive" || s == "trasmissive" || s == "outflow" || s == "open") return 0;
    if (s == "1" || s == "periodic") return 1;
    if (s == "2" || s == "reflective" || s == "reflecting" || s == "box" || s == "box-like") return 2;

    try {
        int v = std::stoi(s);
        return (v >= 0 && v <= 2) ? v : fallback;
    } catch (...) {
        return fallback;
    }
}

int read_face_bc(const std::unordered_map<std::string, std::string>& cfg,
                 const std::string& k1,
                 const std::string& k2 = "")
{
    auto it = cfg.find(k1);
    if (it != cfg.end()) return parse_bc_code(it->second, -1);

    if (!k2.empty()) {
        it = cfg.find(k2);
        if (it != cfg.end()) return parse_bc_code(it->second, -1);
    }

    return -1; // fallback to global BC
}
} // namespace

Params read_param(std::string fname)
{
    Params p;

    std::unordered_map<std::string, std::string> configData;
    std::ifstream configFileStream(fname);
    for (std::string line{}; std::getline(configFileStream, line);)
    {
        // skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        // find the colon separator
        size_t colon_pos = line.find(':');
        if (colon_pos == std::string::npos) continue;

        std::string key   = line.substr(0, colon_pos);
        std::string value = line.substr(colon_pos + 1);

        // trim leading/trailing whitespace and carriage returns from both
        auto trim = [](std::string& s) {
            s.erase(0, s.find_first_not_of(" \t\r\n"));
            s.erase(s.find_last_not_of(" \t\r\n") + 1);
        };
        trim(key);
        trim(value);

        if (!key.empty())
            configData[key] = value;
    }

    p.GAMMA      = configData.count("GAMMA")   ? stof(configData["GAMMA"])   : 1.4;
    p.t_max      = configData.count("t_max")   ? stof(configData["t_max"])   : 1.0;
    p.dt_snap    = configData.count("dt_snap") ? stof(configData["dt_snap"]) : 0.1;
    p.BC         = parse_bc_code(configData.count("BoundaryCondition") ? configData["BoundaryCondition"] : "periodic", 1);
    p.DragIntegrator = configData.count("DragIntegrator") ? stoi(configData["DragIntegrator"]) : 1;

    p.PTC    = configData.count("PTC")    ? stoi(configData["PTC"])    : 0;
    p.N_dims = configData.count("N_dims") ? stoi(configData["N_dims"]) : 1;
    p.N_ghost = configData.count("N_ghost") ? stoi(configData["N_ghost"]) : 2;
    p.ghost_in_input = configData.count("ghost_in_input") ? stoi(configData["ghost_in_input"]) : 0;

    p.input_file = configData["input_file"];
    p.input_file.erase(std::remove_if(p.input_file.begin(), p.input_file.end(),
                       [](char c){ return c == '\n' || c == '\r'; }), p.input_file.end());
    p.output_dir = configData["output_dir"];
    p.output_dir.erase(std::remove_if(p.output_dir.begin(), p.output_dir.end(),
                       [](char c){ return c == '\n' || c == '\r'; }), p.output_dir.end());

    p.L  = configData.count("L")  ? stof(configData["L"])  : 1.0;
    p.Ly = configData.count("Ly") ? stof(configData["Ly"]) : p.L;
    p.Lz = configData.count("Lz") ? stof(configData["Lz"]) : p.L;
    p.sound_speed = configData.count("sound_speed") ? stof(configData["sound_speed"]) : -1.0;
    p.RiemannSolver = configData.count("RiemannSolver") ? stoi(configData["RiemannSolver"]) : 1;
    p.apply_reconstruction = configData.count("apply_reconstruction") ? stoi(configData["apply_reconstruction"]) : 1;
    p.CFL      = configData.count("CFL") ? stof(configData["CFL"]) : 0.1;
    p.const_dt = configData.count("dt")  ? stof(configData["dt"])  : -1;
    p.g0       = configData.count("g0")       ? stof(configData["g0"])       : 0.0;
    p.Omega0   = configData.count("Omega0")   ? stof(configData["Omega0"])   : 0.0;
    p.q        = configData.count("q")        ? stof(configData["q"])        : 0.0;
    p.freeze_gas  = configData.count("freeze_gas")  ? stoi(configData["freeze_gas"])  : 0;
    p.inject_dust = configData.count("inject_dust") ? stoi(configData["inject_dust"]) : 0;

    p.BC_xmin = read_face_bc(configData, "BC_xmin", "BoundaryCondition_xmin");
    p.BC_xmax = read_face_bc(configData, "BC_xmax", "BoundaryCondition_xmax");
    if (p.N_dims >= 2) {
        p.BC_ymin = read_face_bc(configData, "BC_ymin", "BoundaryCondition_ymin");
        p.BC_ymax = read_face_bc(configData, "BC_ymax", "BoundaryCondition_ymax");
    }
    if (p.N_dims == 3) {
        p.BC_zmin = read_face_bc(configData, "BC_zmin", "BoundaryCondition_zmin");
        p.BC_zmax = read_face_bc(configData, "BC_zmax", "BoundaryCondition_zmax");
    }

    idx.rho = 0;
    idx.vx  = 1;
    if (p.N_dims >= 2) idx.vy = 2;
    if (p.N_dims == 3) idx.vz = 3;
    idx.P = 1 + p.N_dims;
    if (p.PTC == 1) {
        idx.s11 = 1 + p.N_dims;
        if (p.N_dims >= 2) { idx.s12 = idx.s11 + 1; idx.s22 = idx.s12 + 1; }
        if (p.N_dims == 3) { idx.s13 = idx.s22 + 1; idx.s23 = idx.s13 + 1; idx.s33 = idx.s23 + 1; }
    }

    return p;
}

// ---------------------------------------------------------------
// HDF5 helpers
// ---------------------------------------------------------------
namespace {

// Read a named 1D or ND dataset into a flat std::vector<double>
std::vector<double> h5_read_flat(H5::H5File& file, const std::string& path)
{
    H5::DataSet   ds   = file.openDataSet(path);
    H5::DataSpace sp   = ds.getSpace();
    hsize_t total = 1;
    int ndims = sp.getSimpleExtentNdims();
    std::vector<hsize_t> dims(ndims);
    sp.getSimpleExtentDims(dims.data());
    for (int d = 0; d < ndims; d++) total *= dims[d];
    std::vector<double> buf(total);
    ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
    return buf;
}

} // namespace

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
        c[flat].PTC         = p.PTC;
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

std::vector<Cell> read_ic_hdf5(Params& p, Grid& g)
{
    H5::H5File file(p.input_file, H5F_ACC_RDONLY);

    // --- detect N_dust from number of dust_rho_N datasets ---
    p.N_dust = 0;
    H5::Group fld = file.openGroup("fields");
    while (fld.nameExists("dust_rho_" + std::to_string(p.N_dust + 1)))
        p.N_dust++;

    // --- derive all variable counts from N_dust ---
    p.N_var_gas  = 2 + p.N_dims;
    p.N_var_dust = 1 + p.N_dims;
    if (p.PTC) p.N_var_dust += (int)(0.5 * p.N_dims * (p.N_dims + 1));
    p.N_vars = p.N_var_gas + p.N_dust * p.N_var_dust;

    // --- read K values from IC file ---
    p.K.resize(p.N_dust, 0.0);
    H5::Group drag = file.openGroup("drag");
    for (int s = 1; s <= p.N_dust; s++) {
        std::string key = "K_" + std::to_string(s);
        drag.openAttribute(key).read(H5::PredType::NATIVE_DOUBLE, &p.K[s-1]);
    }

    // --- gas fields ---
    auto gas_rho = h5_read_flat(file, "fields/gas_rho");
    int total_active = (int)gas_rho.size();

    if (p.N_cells < 0)
        p.N_cells = (p.N_dims == 1) ? total_active
                  : (p.N_dims == 2) ? (int)round(sqrt((double)total_active))
                                    : (int)round(cbrt((double)total_active));
    if (p.N_cells_y < 0) p.N_cells_y = (p.N_dims >= 2) ? p.N_cells : 1;
    if (p.N_cells_z < 0) p.N_cells_z = (p.N_dims == 3) ? p.N_cells : 1;

    p.dx = p.L  / p.N_cells;
    p.dy = p.Ly / p.N_cells_y;
    p.dz = p.Lz / p.N_cells_z;
    const_cast<Grid&>(g) = Grid(p);

    auto gas_vx = h5_read_flat(file, "fields/gas_vx");
    auto gas_P  = h5_read_flat(file, "fields/gas_P");
    std::vector<double> gas_vy, gas_vz;
    if (p.N_dims >= 2) gas_vy = h5_read_flat(file, "fields/gas_vy");
    if (p.N_dims == 3) gas_vz = h5_read_flat(file, "fields/gas_vz");

    // --- dust fields: one dataset per species ---
    std::vector<std::vector<double>> dust_rho(p.N_dust), dust_vx(p.N_dust);
    std::vector<std::vector<double>> dust_vy(p.N_dust), dust_vz(p.N_dust);
    for (int s = 1; s <= p.N_dust; s++) {
        std::string suf = "_" + std::to_string(s);
        dust_rho[s-1] = h5_read_flat(file, "fields/dust_rho" + suf);
        dust_vx[s-1]  = h5_read_flat(file, "fields/dust_vx"  + suf);
        if (p.N_dims >= 2) dust_vy[s-1] = h5_read_flat(file, "fields/dust_vy" + suf);
        if (p.N_dims == 3) dust_vz[s-1] = h5_read_flat(file, "fields/dust_vz" + suf);
    }

    // --- allocate & fill cells ---
    std::vector<Cell> c(g.size());
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
        c[flat].PTC         = p.PTC;
        c[flat].x_center    = p.dx * (0.5 - p.N_ghost + i);
        c[flat].y_center    = p.dy * (0.5 - p.N_ghost + j);
        c[flat].z_center    = p.dz * (0.5 - p.N_ghost + k);
        c[flat].initialize();
    }

    int ic = 0;
    for (int k = 0; k < g.Nz; k++)
    for (int j = 0; j < g.Ny; j++)
    for (int i = 0; i < g.Nx; i++)
    {
        int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);

        c[flat].W[0][idx.rho] = gas_rho[ic];
        c[flat].W[0][idx.vx]  = gas_vx[ic];
        if (p.N_dims >= 2) c[flat].W[0][idx.vy] = gas_vy[ic];
        if (p.N_dims == 3) c[flat].W[0][idx.vz] = gas_vz[ic];
        c[flat].W[0][idx.P]   = gas_P[ic];

        for (int s = 1; s <= p.N_dust; s++) {
            c[flat].W[s][idx.rho] = dust_rho[s-1][ic];
            c[flat].W[s][idx.vx]  = dust_vx[s-1][ic];
            if (p.N_dims >= 2) c[flat].W[s][idx.vy] = dust_vy[s-1][ic];
            if (p.N_dims == 3) c[flat].W[s][idx.vz] = dust_vz[s-1][ic];
        }

        c[flat].get_U_from_W();
        ic++;
    }

    apply_boundary_conditions(c, p, g);
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


void write_output_hdf5(const std::vector<Cell>& c, Params& p, Vars& v, const Grid& g)
{
    if (v.t - v.k_snap * p.dt_snap >= 0)
    {
        printf("%lf %d\n", v.t, v.k_snap);
        std::string output_file = p.output_dir + std::to_string(v.k_snap) + ".h5";
        H5::H5File file(output_file, H5F_ACC_TRUNC);

        int N_active = g.Nx * g.Ny * g.Nz;
        hsize_t dims1[1] = { (hsize_t)N_active };

        std::vector<double> gas_rho(N_active), gas_vx(N_active), gas_P(N_active);
        std::vector<double> gas_vy, gas_vz;
        if (p.N_dims >= 2) gas_vy.resize(N_active);
        if (p.N_dims == 3) gas_vz.resize(N_active);

        // one vector per species
        std::vector<std::vector<double>> dust_rho(p.N_dust, std::vector<double>(N_active));
        std::vector<std::vector<double>> dust_vx (p.N_dust, std::vector<double>(N_active));
        std::vector<std::vector<double>> dust_vy, dust_vz;
        if (p.N_dims >= 2) dust_vy.assign(p.N_dust, std::vector<double>(N_active));
        if (p.N_dims == 3) dust_vz.assign(p.N_dust, std::vector<double>(N_active));

        std::vector<double> x_coord(N_active), y_coord(N_active), z_coord(N_active);

        int ic = 0;
        for (int k = 0; k < g.Nz; k++)
        for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
        {
            int flat = g.flat_idx(i + p.N_ghost, j + p.N_ghost, k + p.N_ghost);
            gas_rho[ic] = c[flat].W[0][idx.rho];
            gas_vx[ic]  = c[flat].W[0][idx.vx];
            gas_P[ic]   = c[flat].W[0][idx.P];
            if (p.N_dims >= 2) gas_vy[ic] = c[flat].W[0][idx.vy];
            if (p.N_dims == 3) gas_vz[ic] = c[flat].W[0][idx.vz];

            for (int s = 1; s <= p.N_dust; s++) {
                dust_rho[s-1][ic] = c[flat].W[s][idx.rho];
                dust_vx[s-1][ic]  = c[flat].W[s][idx.vx];
                if (p.N_dims >= 2) dust_vy[s-1][ic] = c[flat].W[s][idx.vy];
                if (p.N_dims == 3) dust_vz[s-1][ic] = c[flat].W[s][idx.vz];
            }

            x_coord[ic] = c[flat].x_center;
            y_coord[ic] = c[flat].y_center;
            z_coord[ic] = c[flat].z_center;
            ic++;
        }

        H5::Group fld = file.createGroup("fields");
        auto write1D = [&](H5::Group& grp, const std::string& name, const std::vector<double>& buf) {
            H5::DataSpace sp(1, dims1);
            grp.createDataSet(name, H5::PredType::NATIVE_DOUBLE, sp)
               .write(buf.data(), H5::PredType::NATIVE_DOUBLE);
        };

        write1D(fld, "gas_rho", gas_rho);
        write1D(fld, "gas_vx",  gas_vx);
        write1D(fld, "gas_P",   gas_P);
        if (p.N_dims >= 2) write1D(fld, "gas_vy", gas_vy);
        if (p.N_dims == 3) write1D(fld, "gas_vz", gas_vz);

        for (int s = 1; s <= p.N_dust; s++) {
            std::string suf = "_" + std::to_string(s);
            write1D(fld, "dust_rho" + suf, dust_rho[s-1]);
            write1D(fld, "dust_vx"  + suf, dust_vx[s-1]);
            if (p.N_dims >= 2) write1D(fld, "dust_vy" + suf, dust_vy[s-1]);
            if (p.N_dims == 3) write1D(fld, "dust_vz" + suf, dust_vz[s-1]);
        }

        H5::Group coords = file.createGroup("coords");
        write1D(coords, "x", x_coord);
        if (p.N_dims >= 2) write1D(coords, "y", y_coord);
        if (p.N_dims == 3) write1D(coords, "z", z_coord);

        H5::Group meta = file.createGroup("params");
        H5::DataSpace scalar_space(H5S_SCALAR);
        meta.createAttribute("t",      H5::PredType::NATIVE_DOUBLE, scalar_space).write(H5::PredType::NATIVE_DOUBLE, &v.t);
        meta.createAttribute("k_snap", H5::PredType::NATIVE_INT,    scalar_space).write(H5::PredType::NATIVE_INT,    &v.k_snap);
        meta.createAttribute("N_dust", H5::PredType::NATIVE_INT,    scalar_space).write(H5::PredType::NATIVE_INT,    &p.N_dust);

        v.k_snap += 1;
    }
}