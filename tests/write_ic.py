import numpy as np

N = 1000

# SOUNDWAVE
x = np.linspace(0,2*np.pi,N)
A = 1e-4
gasrho = np.ones((N)) + A * np.sin(x) 
gasvel = A * np.sin(x)
dustrho = np.ones((N)) + A * np.sin(x) 
dustvel = A * np.sin(x)
cs = 1.
GAMMA = 1.00001
P = gasrho * cs * cs / GAMMA
W = np.vstack((gasrho,gasvel,P,dustrho,dustvel))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYWAVE/wave.inp", W)

# SHOCKS
# K = inf
N_part = 500
gasrhoL = np.ones((N_part))
dustrhoL = np.ones((N_part))
gasvelL = np.zeros((N_part))
dustvelL = np.zeros((N_part))
PL = np.ones((N_part))

gasrhoR = np.ones((N_part)) * 0.125
dustrhoR = np.ones((N_part)) * 0.125
gasvelR = np.zeros((N_part))
dustvelR = np.zeros((N_part))
PR = np.ones((N_part)) * 0.1

gasrho = np.concatenate((gasrhoL, gasrhoR))
dustrho = np.concatenate((dustrhoL, dustrhoR))
gasvel = np.concatenate((gasvelL, gasvelR))
dustvel = np.concatenate((dustvelL, dustvelR))
P = np.concatenate((PL, PR))

W = np.vstack((gasrho,gasvel,P,dustrho,dustvel))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYSHOCK/shock_Kinf.inp", W)

# K = 0
N_part = 500
gasrhoL = np.ones((N_part))
dustrhoL = np.ones((N_part)) * 0.125
gasvelL = np.zeros((N_part))
dustvelL = np.zeros((N_part))
PL = np.ones((N_part))

gasrhoR = np.ones((N_part)) * 0.125
dustrhoR = np.ones((N_part)) * 0.125
gasvelR = np.zeros((N_part))
dustvelR = np.zeros((N_part))
PR = np.ones((N_part)) * 0.1

gasrho = np.concatenate((gasrhoL, gasrhoR))
dustrho = np.concatenate((dustrhoL, dustrhoR))
gasvel = np.concatenate((gasvelL, gasvelR))
dustvel = np.concatenate((dustvelL, dustvelR))
P = np.concatenate((PL, PR))

W = np.vstack((gasrho,gasvel,P,dustrho,dustvel))
W = np.swapaxes(W, 0, 1)
np.savetxt("DUSTYSHOCK/shock_K0.inp", W)
