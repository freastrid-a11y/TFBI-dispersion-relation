import numpy as np

NR_MAX_ITER = 140
NR_MAX_BACKTRACK = 10
nu_i = np.array([40, 75, 140, 261, 488, 912, 1705, 3188])         # Ion-neutral collision freq (s^-1)
nu_e = np.array([340, 635, 1187, 2219, 4148, 7755, 14497, 27101]) # Electron-neutral collision freq (s^-1)
m_i = 5e-26                                             # Ion mass
m_e = 9.1094e-31                                          # Electron mass
m_n_e = 1e-27
m_n_i = 5e-26
e = 1.6e-19                                             # Elementary charge (C)
eps_0 = 8.854e-12                                       # vacuum permittivity
delta_en = 3e-3                                         # electron energy loss fraction

B_0 = 50000e-9                                          # Static magnetic field (Tesla)
V_d = 2100                                              # drift velocity
E_0 = V_d * B_0

gamma = 1.0                                             # Isothermal assumption for electrons
k_B = 1.380649e-23                                      # Boltzmann constant (J/K)
T_n = 300 * k_B                                         # Ion temperature in Joules (300 K * k_B)
T_e = T_n + (2 * m_e * E_0**2) / (3 * delta_en * B_0**2)# Electron temperature in Joules (300 K * k_B)
C_s_squared = (T_n + gamma * T_e) / m_i
k_y = 0.0                                               # perpendicular to drift
V_th = np.sqrt(T_n/m_i)                                 # thermal velocity
Omega_i = (e * B_0) / m_i                               # ion cyclotron frequency
Omega_e = (e * B_0) / m_e                               # electron cyclotron frequency
kappa = Omega_e / nu_e
kappa_e = Omega_e / nu_e
kappa_i = np.array([4.01, 2.15, 1.15, 0.61, 0.33, 0.18, 0.09, 0.05])
Psi_0 = 1 / (kappa_i*kappa_e)

n_0 = 10e10#4e8                                               # mean ion density
omega_pe_squared = n_0*e**2/(eps_0*m_e)                 # electron plasma frequency
omega_pi_squared = n_0*e**2/(eps_0*m_i)                 # ion plasma frequency


lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)              # Debye length squared

# Electron drift (V_d) is on the x-axis, so Electric Field (E_0) is on the y-axis.

a_i0_mag = (e * E_0) / m_i                              # Acceleration vector (parallel to E_0)
a_i0_x = 0.0
a_i0_y = a_i0_mag

U0_mag = (e * E_0) / (m_i * nu_i)                       # Primary Ion Pedersen drift
U0_x = 0
U0_y = U0_mag

beta_T = T_n / T_e                                      # Temperature Ratio (beta_T)
V_Tn = np.sqrt(T_n / m_i)                               # Neutral Thermal Velocity

# Anisotropic Thermal Velocities
# Heat is primarily generated along the Pedersen drift axis (y-axis)
T_eff = T_n + (m_i * U0_mag**2) / 3.0

anisotropy_ratio = 0.5
T_y = (3.0 * T_eff) / (1.0 + 2.0 * anisotropy_ratio)
T_x = anisotropy_ratio * T_y

V_Tx = np.sqrt(T_x / m_i)
V_Ty = np.sqrt(T_y / m_i)
T_eff_i = m_i*(V_Tx**2+2*V_Ty**2)

alt_eff = [135, 129, 122, 116, 110, 104, 97, 91]        # effective altitudes


if __name__ == "__main__":
    print(T_e/T_n)