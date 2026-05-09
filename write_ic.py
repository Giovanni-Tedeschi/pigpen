import numpy as np
import h5py 

# ============================================================
# Generic helpers
# ============================================================

def write_ic(filepath, fields):
    """
    Write an IC file in row-major order (k outer, j middle, i inner).

    Parameters
    ----------
    filepath : str
    fields   : list of np.ndarray, all the same shape
               Order must match C++ variable ordering:
               [gas_rho, gas_vx, (gas_vy), (gas_vz), gas_P,
                dust1_rho, dust1_vx, (dust1_vy), (dust1_vz),
                dust2_rho, ...]
    """
    W = np.array([f.ravel(order='C') for f in fields]).T  # (N_total, N_vars)
    np.savetxt(filepath, W)

def bins_per_decade_from_grid(mass_grid):
    """
    Infer bins_per_decade from a log-spaced mass grid.
    Uses the ratio between the first two entries: m[1]/m[0] = 10^(1/bpd)
    """
    ratio = mass_grid[1] / mass_grid[0]
    return 1.0 / np.log10(ratio)

def write_ic_hdf5(filepath, fields, dims=1, N_dust=1, has_ptc=False, K=None, mass_grid=None):
    """
    Write an IC file in HDF5 format with auto-generated field names.

    Parameters
    ----------
    filepath  : str, should end in .h5
    fields    : list of np.ndarray in standard variable order
    dims      : number of spatial dimensions
    N_dust    : number of dust species
    has_ptc   : whether PTC stress tensor fields are included
    K         : list of drag coefficients, one per dust species (length N_dust)
    """
    if K is None:
        K = [0.0] * N_dust
    assert len(K) == N_dust, f"K must have length N_dust={N_dust}, got {len(K)}"

    # --- auto-build var_names ---
    var_names = ['gas_rho', 'gas_vx']
    if dims >= 2: var_names += ['gas_vy']
    if dims == 3: var_names += ['gas_vz']
    var_names += ['gas_P']
    for s in range(1, N_dust + 1):
        var_names += [f'dust_rho_{s}', f'dust_vx_{s}']
        if dims >= 2: var_names += [f'dust_vy_{s}']
        if dims == 3: var_names += [f'dust_vz_{s}']
        if has_ptc:
            var_names += [f'dust_s11_{s}', f'dust_s12_{s}', f'dust_s22_{s}']

    assert len(var_names) == len(fields), \
        f"Field count mismatch: {len(fields)} arrays but {len(var_names)} names.\n" \
        f"Expected: {var_names}"

    with h5py.File(filepath, 'w') as f:
        grp = f.create_group('fields')
        for name, arr in zip(var_names, fields):
            grp.create_dataset(name, data=np.asarray(arr, dtype=np.float64).ravel(order='C'))

        # --- drag coefficients, one per dust species ---
        drag = f.create_group('drag')
        for s in range(1, N_dust + 1):
            drag.attrs[f'K_{s}'] = float(K[s-1])

        # mass grid: only m and dm, no rho (that lives in the dust fluid fields)
        if mass_grid is not None:
            dm = mass_grid * np.log(10.0) / bins_per_decade_from_grid(mass_grid)
            coag = f.create_group('coagulation')
            coag.create_dataset('m',  data=np.asarray(mass_grid, dtype=np.float64))
            coag.create_dataset('dm', data=np.asarray(dm,        dtype=np.float64))



def make_grid(N, dims=1):
    """
    Returns (x, y, z) coordinate arrays of shape matching dims.
    In 1D: shape (N,)
    In 2D: shape (Ny, Nx) -- j outer, i inner
    In 3D: shape (Nz, Ny, Nx) -- k outer, j middle, i inner
    N can be a single int (uniform grid) or tuple (Nx, Ny, Nz).
    """
    if isinstance(N, int):
        N = (N,) * dims
    Nx = N[0]
    Ny = N[1] if dims >= 2 else 1
    Nz = N[2] if dims == 3 else 1

    x1d = np.linspace(0, 1, Nx, endpoint=False) + 0.5/Nx
    y1d = np.linspace(0, 1, Ny, endpoint=False) + 0.5/Ny
    z1d = np.linspace(0, 1, Nz, endpoint=False) + 0.5/Nz

    if dims == 1:
        return x1d, np.zeros_like(x1d), np.zeros_like(x1d)
    elif dims == 2:
        xx, yy = np.meshgrid(x1d, y1d, indexing='ij')   # shape (Nx, Ny)
        xx = xx.T                                         # shape (Ny, Nx) -- row-major
        yy = yy.T
        return xx, yy, np.zeros_like(xx)
    else:
        xx, yy, zz = np.meshgrid(x1d, y1d, z1d, indexing='ij')  # (Nx, Ny, Nz)
        xx = np.transpose(xx, (2, 1, 0))   # shape (Nz, Ny, Nx) -- row-major
        yy = np.transpose(yy, (2, 1, 0))
        zz = np.transpose(zz, (2, 1, 0))
        return xx, yy, zz


def _build_fields(x, dims, gas, dust_list, ptc_list=None):
    """
    Generic field list builder.

    Parameters
    ----------
    x         : x coordinate array (used only for shape via ones_like)
    dims      : number of spatial dimensions
    gas       : dict with keys 'rho', 'vx', 'vy', 'vz', 'P'
    dust_list : list of dicts, each with keys 'rho', 'vx', 'vy', 'vz'
    ptc_list  : list of dicts with keys 's11', 's12', 's22' (one per dust species, optional)
    """
    fields = [gas['rho'], gas['vx']]
    if dims >= 2: fields += [gas['vy']]
    if dims == 3: fields += [gas['vz']]
    fields += [gas['P']]
    for i, d in enumerate(dust_list):
        fields += [d['rho'], d['vx']]
        if dims >= 2: fields += [d['vy']]
        if dims == 3: fields += [d['vz']]
        if ptc_list is not None:
            fields += [ptc_list[i]['s11'], ptc_list[i]['s12'], ptc_list[i]['s22']]
    return fields


# ============================================================
# DUSTY BOX
# ============================================================

def write_box_A(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    gas = {'rho': ones,   'vx': ones,     'vy': ones*0, 'vz': ones*0, 'P': ones}
    d1  = {'rho': ones,   'vx': ones*2,   'vy': ones*0, 'vz': ones*0}
    d2  = {'rho': ones,   'vx': ones*0.5, 'vy': ones*0, 'vz': ones*0}
    fields = _build_fields(x, dims, gas, [d1, d2])
    write_ic(f"{folder}/DUSTYBOX/box_A.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYBOX/box_A.h5", fields, dims=dims, N_dust=2, K=[0.5, 1.0])

def write_box_B(N, folder, dims=1):
    write_box_A(N, folder, dims)

def write_box_C(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    gas = {'rho': ones,    'vx': ones,     'vy': ones*0, 'vz': ones*0, 'P': ones}
    d1  = {'rho': ones*10, 'vx': ones*2,   'vy': ones*0, 'vz': ones*0}
    d2  = {'rho': ones*100,'vx': ones*0.5, 'vy': ones*0, 'vz': ones*0}
    fields = _build_fields(x, dims, gas, [d1, d2])
    write_ic(f"{folder}/DUSTYBOX/box_C.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYBOX/box_C.h5", fields, dims=dims, N_dust=2, K=[5.0, 100.0])


# ============================================================
# SOUNDWAVE helpers
# ============================================================

def _wave_phase(x, y, z, direction):
    """Return the phase 2*pi*(coord along sweep direction)."""
    if direction == 'x': return 2 * np.pi * x
    if direction == 'y': return 2 * np.pi * y
    if direction == 'z': return 2 * np.pi * z
    raise ValueError(f"Unknown direction {direction}")


def _make_wave_A_fields(x, y, z, direction, dims):
    A = 1e-4
    cs = 1.0
    GAMMA = 1.00001
    phi = _wave_phase(x, y, z, direction)
    ones = np.ones_like(x)

    gasrho    = ones + A * np.cos(phi)
    vel_wave  = A * (-0.701960 * np.cos(phi) + 0.304924 * np.sin(phi))
    P         = gasrho * cs**2 / GAMMA

    dustrho_1 = ones * 2.24 + A * (0.165251 * np.cos(phi) + 1.247801 * np.sin(phi))
    dvel_wave = A * (-0.221645 * np.cos(phi) - 0.368534 * np.sin(phi))

    vx = vel_wave if direction == 'x' else ones*0
    vy = vel_wave if direction == 'y' else ones*0
    vz = vel_wave if direction == 'z' else ones*0
    dvx = dvel_wave if direction == 'x' else ones*0
    dvy = dvel_wave if direction == 'y' else ones*0
    dvz = dvel_wave if direction == 'z' else ones*0

    gas = {'rho': gasrho,   'vx': vx, 'vy': vy, 'vz': vz, 'P': P}
    d1  = {'rho': dustrho_1,'vx': dvx,'vy': dvy,'vz': dvz}
    return _build_fields(x, dims, gas, [d1])


def _make_wave_B_fields(x, y, z, direction, dims):
    A = 1e-4
    cs = 1.0
    GAMMA = 1.00001
    phi = _wave_phase(x, y, z, direction)
    ones = np.ones_like(x)

    gasrho   = ones + A * np.cos(phi)
    vel_wave = A * (-0.874365 * np.cos(phi) + 0.145215 * np.sin(phi))
    P        = gasrho * cs**2 / GAMMA

    dust_data = [
        (0.1,      A * (0.080588 * np.cos(phi) + 0.048719 * np.sin(phi)),
                   A * (-0.775380 * np.cos(phi) - 0.308952 * np.sin(phi))),
        (0.233333, A * (0.091607 * np.cos(phi) + 0.134955 * np.sin(phi)),
                   A * (-0.427268 * np.cos(phi) - 0.448704 * np.sin(phi))),
        (0.366667, A * (0.030927 * np.cos(phi) + 0.136799 * np.sin(phi)),
                   A * (-0.127928 * np.cos(phi) - 0.313967 * np.sin(phi))),
        (0.500000, A * (0.001451 * np.cos(phi) + 0.090989 * np.sin(phi)),
                   A * (-0.028963 * np.cos(phi) - 0.158693 * np.sin(phi))),
    ]

    vx = vel_wave if direction == 'x' else ones*0
    vy = vel_wave if direction == 'y' else ones*0
    vz = vel_wave if direction == 'z' else ones*0
    gas = {'rho': gasrho, 'vx': vx, 'vy': vy, 'vz': vz, 'P': P}

    dust_list = []
    for rho0, drho, dvel in dust_data:
        dvx = dvel if direction == 'x' else ones*0
        dvy = dvel if direction == 'y' else ones*0
        dvz = dvel if direction == 'z' else ones*0
        dust_list.append({'rho': ones*rho0 + drho, 'vx': dvx, 'vy': dvy, 'vz': dvz})

    return _build_fields(x, dims, gas, dust_list)


# ============================================================
# SOUNDWAVE_A
# ============================================================

def write_wave_A(N, folder, dims=1, direction='x'):
    x, y, z = make_grid(N, dims)
    fields = _make_wave_A_fields(x, y, z, direction, dims)
    write_ic(f"{folder}/DUSTYWAVE/wave_A_{dims}D_{direction}.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYWAVE/wave_A_{dims}D_{direction}.h5",
                  fields, dims=dims, N_dust=1, K=[5.6])


# ============================================================
# SOUNDWAVE_B
# ============================================================

def write_wave_B(N, folder, dims=1, direction='x'):
    x, y, z = make_grid(N, dims)
    fields = _make_wave_B_fields(x, y, z, direction, dims)
    write_ic(f"{folder}/DUSTYWAVE/wave_B_{dims}D_{direction}.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYWAVE/wave_B_{dims}D_{direction}.h5",
                  fields, dims=dims, N_dust=4, K=[1.0, 1.083038, 0.789959, 0.5])


# ============================================================
# EXT_FORCE
# ============================================================

def write_ext_force(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    cs = 1.0
    GAMMA = 1.00001
    gas = {'rho': ones,     'vx': ones*2.0,  'vy': ones*0, 'vz': ones*0, 'P': ones*cs**2/GAMMA}
    d1  = {'rho': ones*0.1, 'vx': ones*0.1,  'vy': ones*0, 'vz': ones*0}
    d2  = {'rho': ones*0.1, 'vx': ones*-0.5, 'vy': ones*0, 'vz': ones*0}
    fields = _build_fields(x, dims, gas, [d1, d2])
    write_ic(f"{folder}/EXT_FORCE/ext_force.inp", fields)
    write_ic_hdf5(f"{folder}/EXT_FORCE/ext_force.h5", fields, dims=dims, N_dust=2, K=[0.1, 0.075])


# ============================================================
# SHOCK_B
# ============================================================

def write_shock_B(N, folder):
    NL = int(N * 0.1)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    ones = np.ones(N)
    gas = {'rho': lr(1,16), 'vx': lr(2.0,0.125), 'vy': ones*0, 'vz': ones*0, 'P': lr(1,16)}
    dust_list = [{'rho': lr(1,16), 'vx': lr(2.0,0.125), 'vy': ones*0, 'vz': ones*0}
                 for _ in range(3)]
    fields = _build_fields(ones, 1, gas, dust_list)
    write_ic(f"{folder}/DUSTYSHOCK/shock_B.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYSHOCK/shock_B.h5", fields, dims=1, N_dust=3, K=[1.0, 3.0, 5.0])


# ============================================================
# SHOCK PTC
# ============================================================

def write_shock_PTC(N, folder):
    NL = int(N * 0.5)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    ones = np.ones(N)
    gas = {'rho': lr(1,1), 'vx': lr(0,0), 'vy': lr(0,0), 'vz': ones*0, 'P': lr(1,1)}
    d1  = {'rho': lr(1,0.125), 'vx': lr(0,0), 'vy': lr(0,0), 'vz': ones*0}
    ptc = [{'s11': lr(0.6, 0.2/0.125), 's12': lr(0.05, 0.1/0.125), 's22': lr(2.0, 0.2/0.125)}]
    fields = _build_fields(ones, 2, gas, [d1], ptc_list=ptc)
    write_ic(f"{folder}/PTC_SHOCK/shock_PTC.inp", fields)
    write_ic_hdf5(f"{folder}/PTC_SHOCK/shock_PTC.h5", fields, dims=2, N_dust=1, has_ptc=True, K=[0.0])

def write_shock_PTC_vacuum(N, folder):
    NL = int(N * 0.5)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    ones = np.ones(N)
    gas = {'rho': lr(1,1), 'vx': lr(0,0), 'vy': lr(0,0), 'vz': ones*0, 'P': lr(1,1)}
    d1  = {'rho': lr(1,0), 'vx': lr(0,0), 'vy': lr(0,0), 'vz': ones*0}
    ptc = [{'s11': lr(2.0,0), 's12': lr(0.05,0), 's22': lr(0.6,0)}]
    fields = _build_fields(ones, 2, gas, [d1], ptc_list=ptc)
    write_ic(f"{folder}/PTC_SHOCK/shock_PTC_vacuum.inp", fields)
    write_ic_hdf5(f"{folder}/PTC_SHOCK/shock_PTC_vacuum.h5", fields, dims=2, N_dust=1, has_ptc=True, K=[0.0])


# ============================================================
# STEADY_STATE_DRIFT  (2D, includes ghosts)
# ============================================================

def write_steady_state_drift(N, folder):
    N_ghost = 2
    cs = 0.05
    omega0 = 1
    chi0 = 0.005
    q = 1.5
    k2tilde = 2 * (2 - q)
    N_dust = 1
    ts = 0.1
    eps = 3.0

    A = 0; B = 1
    for i in range(N_dust):
        A += k2tilde * eps * ts / (1. + k2tilde * ts**2)
        B += eps / (1. + k2tilde * ts**2)
    psi = (A**2 + k2tilde * B**2) ** -1

    L = 1
    dx = L / N
    x = np.arange(dx*(0.5 - N_ghost), L + dx*N_ghost, dx)

    ones = np.ones_like(x)
    gasrho   = ones
    gasvelx  = ones * (A * chi0 * psi)
    gasvely  = ones * (-0.5 * k2tilde * B * chi0 * psi) - q * omega0 * x
    P        = ones * cs**2
    dustrho  = ones * eps
    dustvelx = (gasvelx + 2*ts*(gasvely + q*omega0*x)) / (1. + k2tilde * ts**2)
    dustvely = ((gasvely + q*omega0*x) - (2.-q)*ts*gasvelx) / (1. + k2tilde * ts**2) - q*omega0*x

    gas = {'rho': gasrho, 'vx': gasvelx, 'vy': gasvely, 'vz': ones*0, 'P': P}
    d1  = {'rho': dustrho,'vx': dustvelx,'vy': dustvely,'vz': ones*0}
    fields = _build_fields(x, 2, gas, [d1])
    write_ic(f"{folder}/STEADY_STATE_DRIFT/linA.inp", fields)
    write_ic_hdf5(f"{folder}/STEADY_STATE_DRIFT/linA.h5", fields, dims=2, N_dust=1, K=[30.0])


# ============================================================
# KELVIN-HELMHOLTZ
# ============================================================

def write_kelvin_helmholtz(N, folder, dims=2, dust_to_gas=0.01):
    assert dims >= 2, "KH instability requires at least 2D"
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)

    rho1, rho2 = 1.0, 2.0
    Dy   = 0.025
    Drho = 0.5*(rho2-rho1)
    vx1, vx2 = -0.5, 0.5
    A0 = 0.01
    P  = 2.5

    gasrho = rho1 + Drho * (np.tanh((y-0.25)/Dy) - np.tanh((y-0.75)/Dy) - 1.0)
    vx = 0.5*(vx2-vx1) * (np.tanh((y-0.25)/Dy) - np.tanh((y-0.75)/Dy) - 1.0) + vx2
    vy = A0 * np.sin(4 * np.pi * x)

    gas = {'rho': gasrho,           'vx': vx, 'vy': vy, 'vz': ones*0, 'P': ones*P}
    d1  = {'rho': gasrho*dust_to_gas,'vx': vx, 'vy': vy, 'vz': ones*0}
    fields = _build_fields(x, dims, gas, [d1])
    write_ic(f"{folder}/DUSTYKH/kh_{dims}D.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYKH/kh_{dims}D.h5", fields, dims=dims, N_dust=1, K=[1.0])


def write_kelvin_helmholtz_PTC(N, folder, dims=2, dust_to_gas=0.01):
    assert dims >= 2, "KH instability requires at least 2D"
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    zeros = np.zeros_like(x)

    rho1, rho2 = 1.0, 2.0
    Dy   = 0.025
    Drho = 0.5*(rho2-rho1)
    vx1, vx2 = -0.5, 0.5
    A0 = 0.01
    P  = 2.5

    gasrho = rho1 + Drho * (np.tanh((y-0.25)/Dy) - np.tanh((y-0.75)/Dy) - 1.0)
    vx = 0.5*(vx2-vx1) * (np.tanh((y-0.25)/Dy) - np.tanh((y-0.75)/Dy) - 1.0) + vx2
    vy = A0 * np.sin(4 * np.pi * x)

    gas = {'rho': gasrho,            'vx': vx, 'vy': vy, 'vz': ones*0, 'P': ones*P}
    d1  = {'rho': gasrho*dust_to_gas,'vx': vx, 'vy': vy, 'vz': ones*0}
    ptc = [{'s11': ones*1e-8, 's12': zeros, 's22': ones*1e-8}]
    fields = _build_fields(x, dims, gas, [d1], ptc_list=ptc)
    write_ic(f"{folder}/DUSTYKH/kh_{dims}D_ptc.inp", fields)
    write_ic_hdf5(f"{folder}/DUSTYKH/kh_{dims}D_ptc.h5", fields, dims=dims, N_dust=1, has_ptc=True, K=[1.0])


# ============================================================
# JET
# ============================================================

def write_jet(folder):
    N = (250, 100)
    dims = 2
    strain = 1
    Lx, Ly = 6, 2

    x, y, _ = make_grid(N, dims) * np.array([Lx, Ly, 0])[:,np.newaxis, np.newaxis]
    ones  = np.ones_like(x)
    zeros = np.zeros_like(x)

    gas = {'rho': ones, 'vx': ones*0.2, 'vy': (1-y)*strain, 'vz': zeros, 'P': ones}
    d1  = {'rho': ones*1e-8, 'vx': ones*0.2, 'vy': zeros, 'vz': zeros}
    fields = _build_fields(x, dims, gas, [d1])
    write_ic(f"{folder}/JET/jet.inp", fields)
    write_ic_hdf5(f"{folder}/JET/jet.h5", fields, dims=dims, N_dust=1, K=[1.0])

    ptc = [{'s11': ones*1e-8, 's12': ones*1e-8, 's22': ones*1e-8}]
    fields_ptc = _build_fields(x, dims, gas, [d1], ptc_list=ptc)
    write_ic(f"{folder}/JET/jet_ptc.inp", fields_ptc)
    write_ic_hdf5(f"{folder}/JET/jet_ptc.h5", fields_ptc, dims=dims, N_dust=1, has_ptc=True, K=[1.0])



def _coag_initial_conditions(m_min=1e-12, m_max=1e2, bins_per_decade=7, t0=1e-8):
    N_m = int(np.log10(m_max / m_min) * bins_per_decade) + 1
    m   = m_min * 10.0 ** (np.arange(N_m) / bins_per_decade)
    dm  = m * np.log(10.0) / bins_per_decade

    m0    = m[0]
    N0loc = 1.0 / m0
    alpha = 1.0
    tau0  = 2.0 / (alpha * N0loc * t0)

    # this IS the rho array that goes directly into dust_rho_1 ... dust_rho_Nm
    rho = (N0loc / m0) * tau0**2 * np.exp(tau0 * (1.0 - m / m0)) * m
    return m, dm, rho

def write_coagulation_test(N, folder, dims=2, bins_per_decade=7):
    import os
    os.makedirs(f"{folder}/COAGULATION", exist_ok=True)

    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)

    m, dm, rho_dist = _coag_initial_conditions()
    N_m = len(m)

    cs    = 1.0
    GAMMA = 1.4
    gas = {'rho': ones, 'vx': ones*0, 'vy': ones*0, 'vz': ones*0,
           'P':   ones * cs**2 / GAMMA}

    # one dust fluid per mass bin, uniform across all cells
    dust_list = [{'rho': ones * rho_dist[s], 'vx': ones*0, 'vy': ones*0, 'vz': ones*0}
                 for s in range(N_m)]

    fields = _build_fields(x, dims, gas, dust_list)

    write_ic_hdf5(
        f"{folder}/COAGULATION/coag_test_{dims}D.h5",
        fields,
        dims=dims,
        N_dust=N_m,
        K=[0.0] * N_m,   # drag not used in pure coagulation test
        mass_grid=m,
    )