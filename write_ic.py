import numpy as np

# BOX_multi_A
def write_box_A(N, folder):
    gasrho = np.ones((N)) * (1)
    gasvel = np.ones((N)) * (1)
    dustrho_1 = np.ones((N)) * (1)
    dustvel_1 = np.ones((N)) * (2) 
    dustrho_2 = np.ones((N)) * (1)
    dustvel_2 = np.ones((N)) * (0.5) 
    P = np.ones((N)) * 1.0
    W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYBOX/box_A.inp", W)

# BOX_multi_B
def write_box_B(N, folder):
    gasrho = np.ones((N)) * (1)
    gasvel = np.ones((N)) * (1)
    dustrho_1 = np.ones((N)) * (1)
    dustvel_1 = np.ones((N)) * (2) 
    dustrho_2 = np.ones((N)) * (1)
    dustvel_2 = np.ones((N)) * (0.5) 
    P = np.ones((N)) * 1.0
    W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYBOX/box_B.inp", W)

# BOX_multi_C
def write_box_C(N, folder):
    gasrho = np.ones((N)) * (1)
    gasvel = np.ones((N)) * (1)
    dustrho_1 = np.ones((N)) * (10)
    dustvel_1 = np.ones((N)) * (2) 
    dustrho_2 = np.ones((N)) * (100)
    dustvel_2 = np.ones((N)) * (0.5) 
    P = np.ones((N)) * 1.0
    W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYBOX/box_C.inp", W)


# SOUNDWAVE_A
def write_wave_A(N, folder):
    x = np.linspace(0,1,N)
    A = 1e-4
    gasrho = np.ones((N)) * 1.000000 + A * (1. * np.cos(2*np.pi*x) - 0. * np.sin(2*np.pi*x)) 
    gasvel = np.ones((N)) * A * (-0.701960 * np.cos(2*np.pi*x) + 0.304924 * np.sin(2*np.pi*x)) 
    cs = 1.
    GAMMA = 1.00001
    P = gasrho * cs * cs / GAMMA

    W = np.vstack((gasrho,gasvel,P))
    dustrho_1 = np.ones((N)) * 2.24 + A * (0.165251 * np.cos(2*np.pi*x) + 1.247801 * np.sin(2*np.pi*x)) 
    dustrho_1 = dustrho_1 
    dustvel_1 = np.ones((N)) * A * (-0.221645 * np.cos(2*np.pi*x) - 0.368534 * np.sin(2*np.pi*x)) 
    W = np.vstack((W,dustrho_1,dustvel_1))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYWAVE/wave_A.inp", W)


# SOUNDWAVE_B
def write_wave_B(N, folder):
    x = np.linspace(0,1,N)
    A = 1e-4
    gasrho = np.ones((N)) * 1.000000 + A * (1. * np.cos(2*np.pi*x) - 0. * np.sin(2*np.pi*x)) 
    gasvel = np.ones((N)) * A * (-0.874365 * np.cos(2*np.pi*x) + 0.145215 * np.sin(2*np.pi*x)) 

    dustrho_1 = np.ones((N)) * 0.1 + A * (0.080588 * np.cos(2*np.pi*x) + 0.048719 * np.sin(2*np.pi*x)) 
    dustvel_1 = np.ones((N)) * A * (-0.775380 * np.cos(2*np.pi*x) - 0.308952 * np.sin(2*np.pi*x)) 

    dustrho_2 = np.ones((N)) * 0.233333 + A * (0.091607 * np.cos(2*np.pi*x) + 0.134955 * np.sin(2*np.pi*x)) 
    dustvel_2 = np.ones((N)) * A * (-0.427268 * np.cos(2*np.pi*x) - 0.448704 * np.sin(2*np.pi*x)) 

    dustrho_3 = np.ones((N)) * 0.366667 + A * (0.030927 * np.cos(2*np.pi*x) + 0.136799 * np.sin(2*np.pi*x)) 
    dustvel_3 = np.ones((N)) * A * (-0.127928 * np.cos(2*np.pi*x) - 0.313967 * np.sin(2*np.pi*x)) 

    dustrho_4 = np.ones((N)) * 0.500000 + A * (0.001451 * np.cos(2*np.pi*x) + 0.090989 * np.sin(2*np.pi*x)) 
    dustvel_4 = np.ones((N)) * A * (-0.028963 * np.cos(2*np.pi*x) - 0.158693 * np.sin(2*np.pi*x)) 

    cs = 1.
    GAMMA = 1.00001
    P = gasrho * cs * cs / GAMMA
    W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2,dustrho_3,dustvel_3,dustrho_4,dustvel_4))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYWAVE/wave_B.inp", W)


# SOUNDWAVE_underrelax
def write_wave_underrelax(N, folder):
    x = np.linspace(0,1,N)
    A = 1e-4
    gasrho = np.ones((N)) * 1.000000 + A * (1. * np.cos(2*np.pi*x) - 0. * np.sin(2*np.pi*x)) 
    gasvel = np.ones((N)) * A * (0.9534626030322761 * np.cos(2*np.pi*x) + 2.5963568718457966e-05 * np.sin(2*np.pi*x)) 
    cs = 1.0
    GAMMA = 1.00001
    P = gasrho * cs * cs / GAMMA

    W = np.vstack((gasrho,gasvel,P))
    dustrho_1 = np.ones((N)) * 0.1 + A * (0.0999999657418787 * np.cos(2*np.pi*x) - 5.990780263718104e-05 * np.sin(2*np.pi*x)) 
    dustvel_1 = np.ones((N)) * A * (0.9534622919481045 * np.cos(2*np.pi*x) - 0.0005452349346200779 * np.sin(2*np.pi*x)) 
    W = np.vstack((W,dustrho_1,dustvel_1))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYWAVE/wave_underrelax.inp", W)


# SHOCK_multi_B
def write_shock_B(N, folder):
    N_part = int(N * 0.1)

    gasrhoL = np.ones((N_part))
    gasvelL = np.ones((N_part)) * 2.0

    dustrhoL_1 = np.ones((N_part))
    dustvelL_1 = np.ones((N_part)) * 2.0

    dustrhoL_2 = np.ones((N_part))
    dustvelL_2 = np.ones((N_part)) * 2.0

    dustrhoL_3 = np.ones((N_part))
    dustvelL_3 = np.ones((N_part)) * 2.0
    PL = gasrhoL

    N_part = int(N * 0.9)
    gasrhoR = np.ones((N_part)) * 16
    gasvelR = np.ones((N_part)) * 0.125

    dustrhoR_1 = np.ones((N_part)) * 16
    dustvelR_1 = np.ones((N_part)) * 0.125

    dustrhoR_2 = np.ones((N_part)) * 16
    dustvelR_2 = np.ones((N_part)) * 0.125

    dustrhoR_3 = np.ones((N_part)) * 16
    dustvelR_3 = np.ones((N_part)) * 0.125
    PR = gasrhoR

    gasrho = np.concatenate((gasrhoL, gasrhoR))
    gasvel = np.concatenate((gasvelL, gasvelR))

    dustrho_1 = np.concatenate((dustrhoL_1, dustrhoR_1))
    dustvel_1 = np.concatenate((dustvelL_1, dustvelR_1))

    dustrho_2 = np.concatenate((dustrhoL_2, dustrhoR_2))
    dustvel_2 = np.concatenate((dustvelL_2, dustvelR_2))

    dustrho_3 = np.concatenate((dustrhoL_3, dustrhoR_3))
    dustvel_3 = np.concatenate((dustvelL_3, dustvelR_3))
    P = np.concatenate((PL, PR))

    W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2,dustrho_3,dustvel_3))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/DUSTYSHOCK/shock_B.inp", W)

# EXT_FORCE
def write_ext_force(N, folder):
    x = np.linspace(0,1,N)
    gasrho = np.ones((N)) * 1.0
    gasvel = np.ones((N)) * 2.0
    cs = 1.0
    GAMMA = 1.00001
    P = gasrho * cs * cs / GAMMA

    W = np.vstack((gasrho,gasvel,P))
    dustrho_1 = np.ones((N)) * 0.1 
    dustvel_1 = np.ones((N)) * 0.1
    dustrho_2 = np.ones((N)) * 0.1
    dustvel_2 = np.ones((N)) * -0.5
    W = np.vstack((W,dustrho_1,dustvel_1, dustrho_2, dustvel_2))
    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/EXT_FORCE/ext_force.inp", W)


# STEADY_STATE_DRIFT
def write_steady_state_drift(N, folder):
    N_ghost = 2
    cs = 0.05
    omega0 = 1
    chi0 = 0.005
    q=1.5
    k2 = 2*(2-q) * pow(omega0,2) 
    k2tilde = 2*(2-q)

    N_dust = 1
    ts = 0.1
    eps = 3.0

    A = 0
    B = 1
    for i in range(0,N_dust):
        A += k2tilde * eps * ts / (1. + k2tilde * pow(ts,2))
        B += eps / (1. + k2tilde * pow(ts,2)) 

    psi = pow(pow(A,2) + k2tilde*pow(B,2), -1) 

    L = 1 
    dx = L / N 
    x = np.arange(dx*(1/2-N_ghost), L+(dx*N_ghost), dx)

    gasrho = np.ones((N+2*N_ghost)) * (1)
    gasvelx = np.ones((N+2*N_ghost)) * (A * chi0 * psi)
    gasvely = np.ones((N+2*N_ghost)) *  (-0.5*k2tilde*B*chi0*psi) - q*omega0*x 

    dustrho_i = np.ones((N+2*N_ghost)) * (1 * eps)
    dustvelx_i = (gasvelx + 2*ts*(gasvely + q*omega0*x)) / (1. + k2tilde*pow(ts,2)) 
    dustvely_i = np.ones((N+2*N_ghost)) * ((gasvely + q*omega0*x) - (2.-q)*ts*gasvelx) / (1. + k2tilde*pow(ts,2)) - q*omega0*x

    P = np.ones((N+2*N_ghost)) * 1.0 * pow(cs,2)
    W = np.vstack((gasrho,gasvelx,gasvely, P))
    for i in range(0,N_dust):
        W = np.vstack((W, dustrho_i, dustvelx_i, dustvely_i)) 


    W = np.swapaxes(W, 0, 1)
    np.savetxt(f"{folder}/STEADY_STATE_DRIFT/linA.inp", W) 