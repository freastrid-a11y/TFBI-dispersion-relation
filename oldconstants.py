import numpy as np

NR_MAX_ITER = 140
NR_MAX_BACKTRACK = 10
nu_i = 2.8e3                            # Ion-neutral collision freq (s^-1)
nu_e = 4e4                              # Electron-neutral collision freq (s^-1)
m_e = 9.1e-31                           # Electron mass (kg)
m_i = 31*1.6e-27                        # Ion mass (kg)
e = 1.6e-19                             # Elementary charge (C)
eps_0 = 8.854e-12                       # vacuum permittivity

B_0 = 0.25e-4                           # Static magnetic field (Tesla)

gamma = 1.0                             # Isothermal assumption for electrons
k_B = 1.380649e-23                      # Boltzmann constant (J/K)
T_i = 208.0 * k_B                       # Ion temperature in Joules (208 K * k_B)
T_e = T_i                               # Electron temperature in Joules (208 K * k_B)

C_s_squared = (T_i + gamma * T_e) / m_i
k_y = 0.0                               # perpendicular to drift

V_d = 600                               # drift velocity
V_th = np.sqrt(T_i/m_i)                 # thermal velocity
Omega_i = (e * B_0) / m_i               # ion cyclotron frequency
Omega_e = (e * B_0) / m_e               # electron cyclotron frequency
kappa = Omega_e / nu_e                  # frequency ratio
Psi_0 = (nu_e * nu_i) / (Omega_e * Omega_i)

n_0 = 1e11                              # mean ion density
omega_pe_squared = n_0*e**2/(eps_0*m_e) # electron plasma frequency
omega_pi_squared = n_0*e**2/(eps_0*m_i) # ion plasma frequency

E_0 = V_d * B_0

lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2) # Debye length squared

# Electron drift (V_d) is on the x-axis, so Electric Field (E_0) is on the y-axis.

a_i0_mag = (e * E_0) / m_i              # Acceleration vector (parallel to E_0)
a_i0_x = 0.0
a_i0_y = a_i0_mag

U0_mag = (e * E_0) / (m_i * nu_i)       # Primary Ion Pedersen drift
U0_x = 0
U0_y = U0_mag

beta_T = T_i / T_e                      # Temperature Ratio (beta_T)
V_Tn = np.sqrt(T_i / m_i)               # Neutral Thermal Velocity

# Anisotropic Thermal Velocities
# Heat is primarily generated along the Pedersen drift axis (y-axis)
T_eff = T_i + (m_i * U0_mag**2) / 3.0

anisotropy_ratio = 0.5
T_y = (3.0 * T_eff) / (1.0 + 2.0 * anisotropy_ratio)
T_x = anisotropy_ratio * T_y

V_Tx = np.sqrt(T_x / m_i)
V_Ty = np.sqrt(T_y / m_i)