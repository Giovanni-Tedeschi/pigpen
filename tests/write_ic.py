import numpy as np

N = 512

# BOX_multi_A
gasrho = np.ones((N)) * (1)
gasvel = np.ones((N)) * (1)
dustrho_1 = np.ones((N)) * (1)
dustvel_1 = np.ones((N)) * (2) 
dustrho_2 = np.ones((N)) * (1)
dustvel_2 = np.ones((N)) * (0.5) 
P = np.ones((N)) * 1.0
W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYBOX/box_A.inp", W)

# BOX_multi_B
gasrho = np.ones((N)) * (1)
gasvel = np.ones((N)) * (1)
dustrho_1 = np.ones((N)) * (1)
dustvel_1 = np.ones((N)) * (2) 
dustrho_2 = np.ones((N)) * (1)
dustvel_2 = np.ones((N)) * (0.5) 
P = np.ones((N)) * 1.0
W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYBOX/box_B.inp", W)

# BOX_multi_C
gasrho = np.ones((N)) * (1)
gasvel = np.ones((N)) * (1)
dustrho_1 = np.ones((N)) * (10)
dustvel_1 = np.ones((N)) * (2) 
dustrho_2 = np.ones((N)) * (100)
dustvel_2 = np.ones((N)) * (0.5) 
P = np.ones((N)) * 1.0
W = np.vstack((gasrho,gasvel,P,dustrho_1,dustvel_1,dustrho_2,dustvel_2))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYBOX/box_C.inp", W)


# SOUNDWAVE
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
np.savetxt("DUSTYWAVE/wave_B.inp", W)

N=500

# SHOCK_multi_B
N_part = 40
gasrhoL = np.ones((N_part))
gasvelL = np.ones((N_part)) * 2.0

dustrhoL_1 = np.ones((N_part))
dustvelL_1 = np.ones((N_part)) * 2.0

dustrhoL_2 = np.ones((N_part))
dustvelL_2 = np.ones((N_part)) * 2.0

dustrhoL_3 = np.ones((N_part))
dustvelL_3 = np.ones((N_part)) * 2.0
PL = gasrhoL

N_part = 360
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
np.savetxt("DUSTYSHOCK/shock_B.inp", W)
