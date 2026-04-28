import numpy as np

# ============================================================
# Generic helper
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


# ============================================================
# DUSTY BOX
# ============================================================

def write_box_A(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    fields = [ones*1, ones*1]                              # rho, vx
    if dims >= 2: fields += [ones*0]                       # vy
    if dims == 3: fields += [ones*0]                       # vz
    fields += [ones*1.0]                                   # P
    fields += [ones*1, ones*2]                             # dust1 rho, vx
    if dims >= 2: fields += [ones*0]                       # dust1 vy
    if dims == 3: fields += [ones*0]                       # dust1 vz
    fields += [ones*1, ones*0.5]                           # dust2 rho, vx
    if dims >= 2: fields += [ones*0]                       # dust2 vy
    if dims == 3: fields += [ones*0]                       # dust2 vz
    write_ic(f"{folder}/DUSTYBOX/box_A_{dims}D.inp", fields)

def write_box_B(N, folder, dims=1):
    write_box_A(N, folder, dims)   # same IC as A for now

def write_box_C(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    fields = [ones*1, ones*1]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    fields += [ones*1.0]
    fields += [ones*10,  ones*2]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    fields += [ones*100, ones*0.5]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    write_ic(f"{folder}/DUSTYBOX/box_C_{dims}D.inp", fields)


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
    """
    Constructs field list for SOUNDWAVE_A travelling in `direction`.
    Variable ordering: [rho, vx, (vy), (vz), P, dust_rho, dust_vx, (dust_vy), (dust_vz)]
    The wave velocity amplitude is placed in the component matching `direction`.
    """
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

    # Gas fields
    vx = vel_wave if direction == 'x' else ones*0
    vy = vel_wave if direction == 'y' else ones*0
    vz = vel_wave if direction == 'z' else ones*0
    fields = [gasrho, vx]
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [vz]
    fields += [P]

    # Dust fields
    dvx = dvel_wave if direction == 'x' else ones*0
    dvy = dvel_wave if direction == 'y' else ones*0
    dvz = dvel_wave if direction == 'z' else ones*0
    fields += [dustrho_1, dvx]
    if dims >= 2: fields += [dvy]
    if dims == 3: fields += [dvz]

    return fields


def _make_wave_B_fields(x, y, z, direction, dims):
    """
    Constructs field list for SOUNDWAVE_B travelling in `direction`.
    4 dust species.
    """
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
    fields = [gasrho, vx]
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [vz]
    fields += [P]

    for rho0, drho, dvel in dust_data:
        dvx = dvel if direction == 'x' else ones*0
        dvy = dvel if direction == 'y' else ones*0
        dvz = dvel if direction == 'z' else ones*0
        fields += [ones * rho0 + drho, dvx]
        if dims >= 2: fields += [dvy]
        if dims == 3: fields += [dvz]

    return fields


# ============================================================
# SOUNDWAVE_A  (1D, 2D, 3D  x  x/y/z direction)
# ============================================================

def write_wave_A(N, folder, dims=1, direction='x'):
    x, y, z = make_grid(N, dims)
    fields = _make_wave_A_fields(x, y, z, direction, dims)
    write_ic(f"{folder}/DUSTYWAVE/wave_A_{dims}D_{direction}.inp", fields)


# ============================================================
# SOUNDWAVE_B  (1D, 2D, 3D  x  x/y/z direction)
# ============================================================

def write_wave_B(N, folder, dims=1, direction='x'):
    x, y, z = make_grid(N, dims)
    fields = _make_wave_B_fields(x, y, z, direction, dims)
    write_ic(f"{folder}/DUSTYWAVE/wave_B_{dims}D_{direction}.inp", fields)


# ============================================================
# EXT_FORCE
# ============================================================

def write_ext_force(N, folder, dims=1):
    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    cs = 1.0
    GAMMA = 1.00001
    P = ones * cs**2 / GAMMA
    fields = [ones*1.0, ones*2.0]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    fields += [P]
    fields += [ones*0.1, ones*0.1]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    fields += [ones*0.1, ones*-0.5]
    if dims >= 2: fields += [ones*0]
    if dims == 3: fields += [ones*0]
    write_ic(f"{folder}/EXT_FORCE/ext_force_{dims}D.inp", fields)


# ============================================================
# SHOCK_B  (1D only -- shocks are inherently 1D)
# ============================================================

def write_shock_B(N, folder):
    NL = int(N * 0.1)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    gasrho = lr(1, 16);  gasvel = lr(2.0, 0.125);  P = lr(1, 16)
    fields = [gasrho, gasvel, P]
    for _ in range(3):
        fields += [lr(1, 16), lr(2.0, 0.125)]
    write_ic(f"{folder}/DUSTYSHOCK/shock_B.inp", fields)

# ============================================================
# SHOCK Particle Crossing Trajectory (PTC)
# ============================================================

def write_shock_PTC(N, folder):
    NL = int(N * 0.5)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    gasrho = lr(1, 1);  gasvelx = lr(0.0, 0.0); gasvely = lr(0.0, 0.0);  P = lr(1, 1)
    #dustrho = lr(1, 0.125); dustvelx = lr(0.0, 0.0); dustvely = lr(0.0, 0.0); dusts11 = lr(2.0, 0.2/0.125); dusts12 = lr(0.05, 0.1/0.125); dusts22 = lr(0.6, 0.2/0.125)
    dustrho = lr(1, 0.125); dustvelx = lr(0.0, 0.0); dustvely = lr(0.0, 0.0); dusts11 = lr(0.6, 0.2/0.125); dusts12 = lr(0.05, 0.1/0.125); dusts22 = lr(2.0, 0.2/0.125)
    fields = [gasrho, gasvelx, gasvely, P, dustrho, dustvelx, dustvely, dusts11, dusts12, dusts22]
    write_ic(f"{folder}/PTC_SHOCK/shock_PTC.inp", fields)


def write_shock_PTC_vacuum(N, folder):
    NL = int(N * 0.5)
    NR = N - NL
    def lr(vL, vR): return np.concatenate([np.ones(NL)*vL, np.ones(NR)*vR])
    gasrho = lr(1, 1);  gasvelx = lr(0.0, 0.0); gasvely = lr(0.0, 0.0);  P = lr(1, 1)
    dustrho = lr(1, 0.0); dustvelx = lr(0.0, 0.0); dustvely = lr(0.0, 0.0); dusts11 = lr(2.0, 0.0); dusts12 = lr(0.05, 0.0); dusts22 = lr(0.6, 0.0)
    fields = [gasrho, gasvelx, gasvely, P, dustrho, dustvelx, dustvely, dusts11, dusts12, dusts22]
    write_ic(f"{folder}/PTC_SHOCK/shock_PTC_vacuum.inp", fields)

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
    x = np.arange(dx*(0.5 - N_ghost), L + dx*N_ghost, dx)  # includes ghosts

    ones = np.ones_like(x)
    gasrho   = ones
    gasvelx  = ones * (A * chi0 * psi)
    gasvely  = ones * (-0.5 * k2tilde * B * chi0 * psi) - q * omega0 * x
    P        = ones * cs**2
    dustrho  = ones * eps
    dustvelx = (gasvelx + 2*ts*(gasvely + q*omega0*x)) / (1. + k2tilde * ts**2)
    dustvely = ((gasvely + q*omega0*x) - (2.-q)*ts*gasvelx) / (1. + k2tilde * ts**2) - q*omega0*x

    fields = [gasrho, gasvelx, gasvely, P, dustrho, dustvelx, dustvely]
    write_ic(f"{folder}/STEADY_STATE_DRIFT/linA.inp", fields)


# ============================================================
# KELVIN-HELMHOLTZ
# ============================================================

def write_kelvin_helmholtz(N, folder, dims=2, dust_to_gas=0.01):
    """
    Kelvin-Helmholtz instability IC matching the SPH setup:
    - Smooth tanh density and velocity profiles at y=0.25 and y=0.75
    - Two sinusoidal modes (sin(4*pi*x)) for vy perturbation
    - Dust comoving with gas, dust-to-gas ratio = dust_to_gas
    dims must be >= 2.
    """
    assert dims >= 2, "KH instability requires at least 2D"

    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)

    rho1, rho2 = 1.0, 2.0
    Dy   = 0.025          # shear layer width (matches SPH Dy)
    Drho = 0.5*(rho2-rho1)
    vx1, vx2 = -0.5, 0.5
    A0 = 0.01             # vy perturbation amplitude
    P  = 2.5              # uniform pressure

    # Smooth density profile (matches rho_KH in SPH code)
    gasrho = (rho1
              + Drho * (np.tanh((y - 0.25) / Dy) - np.tanh((y - 0.75) / Dy) - 1.0))

    # Smooth velocity profile (matches get_vel_x in SPH code)
    vx = (0.5*(vx2 - vx1)
          * (np.tanh((y - 0.25) / Dy) - np.tanh((y - 0.75) / Dy) - 1.0)
          + vx2)

    # 2-mode vy perturbation (matches sin(4*pi*x) in SPH code)
    vy = A0 * np.sin(4 * np.pi * x)

    dustrho = gasrho * dust_to_gas

    fields = [gasrho, vx]
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [ones * 0]
    fields += [ones * P]
    fields += [dustrho, vx]       # dust comoving with gas
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [ones * 0]

    write_ic(f"{folder}/DUSTYKH/kh_{dims}D.inp", fields)



def write_kelvin_helmholtz_PTC(N, folder, dims=2, dust_to_gas=0.01):
    """
    Kelvin-Helmholtz instability IC matching the SPH setup:
    - Smooth tanh density and velocity profiles at y=0.25 and y=0.75
    - Two sinusoidal modes (sin(4*pi*x)) for vy perturbation
    - Dust comoving with gas, dust-to-gas ratio = dust_to_gas
    dims must be >= 2.
    """
    assert dims >= 2, "KH instability requires at least 2D"

    x, y, z = make_grid(N, dims)
    ones = np.ones_like(x)
    zeros = np.zeros_like(x)

    rho1, rho2 = 1.0, 2.0
    Dy   = 0.025          # shear layer width (matches SPH Dy)
    Drho = 0.5*(rho2-rho1)
    vx1, vx2 = -0.5, 0.5
    A0 = 0.01             # vy perturbation amplitude
    P  = 2.5              # uniform pressure

    # Smooth density profile (matches rho_KH in SPH code)
    gasrho = (rho1
              + Drho * (np.tanh((y - 0.25) / Dy) - np.tanh((y - 0.75) / Dy) - 1.0))

    # Smooth velocity profile (matches get_vel_x in SPH code)
    vx = (0.5*(vx2 - vx1)
          * (np.tanh((y - 0.25) / Dy) - np.tanh((y - 0.75) / Dy) - 1.0)
          + vx2)

    # 2-mode vy perturbation (matches sin(4*pi*x) in SPH code)
    vy = A0 * np.sin(4 * np.pi * x)

    dustrho = gasrho * dust_to_gas

    fields = [gasrho, vx]
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [ones * 0]
    fields += [ones * P]
    fields += [dustrho, vx]       # dust comoving with gas
    if dims >= 2: fields += [vy]
    if dims == 3: fields += [ones * 0]

    fields += [ones * 1e-8, zeros, ones * 1e-8]

    write_ic(f"{folder}/DUSTYKH/kh_{dims}D_ptc.inp", fields)