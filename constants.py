import numpy as np

NR_MAX_ITER = 140                                       # Maximum Newton-Raphson iterations (dimensionless)
NR_MAX_BACKTRACK = 10                                   # Maximum line search backtracks (dimensionless)

m_i = 5e-26                                             # Ion mass (kg)
m_e = 9.1094e-31                                        # Electron mass (kg)
m_n_e = 1e-27                                           # Neutral mass for electron collisions (kg)
m_n_i = 5e-26                                           # Neutral mass for ion collisions (kg)
e = 1.6e-19                                             # Elementary charge (C)
eps_0 = 8.854e-12                                       # Vacuum permittivity (F/m)
delta_en = 3e-3                                         # Electron energy loss fraction (dimensionless)

B_0 = 5e-05                                          # Static magnetic field (T)
V_d = 4200                                              # Drift velocity (m/s)
E_0 = V_d * B_0                                         # Background electric field (V/m)

gamma = 1.0                                             # Isothermal assumption for electrons / Polytropic index (dimensionless)
k_B = 1.380649e-23                                      # Boltzmann constant (J/K)
T_n = 300 * k_B                                         # Ion/Neutral temperature (J)
T_e = T_n + (2 * m_e * E_0**2) / (3 * delta_en * B_0**2)# Electron temperature (J)
print(T_e/T_n)
C_s_squared = (T_n + gamma * T_e) / m_i                 # Ion acoustic speed squared (m^2/s^2)
k_y = 0.0                                               # Wave vector perpendicular to drift (m^-1)
V_th = np.sqrt(T_n/m_i)                                 # Thermal velocity (m/s)
Omega_i = (e * B_0) / m_i                               # Ion cyclotron / gyro frequency (rad/s)
Omega_e = (e * B_0) / m_e                               # Electron cyclotron / gyro frequency (rad/s)

n_0 = 10e10                                             # Mean ion/plasma density (m^-3)
omega_pe_squared = n_0*e**2/(eps_0*m_e)                 # Electron plasma frequency squared (rad^2/s^2)
omega_pi_squared = n_0*e**2/(eps_0*m_i)                 # Ion plasma frequency squared (rad^2/s^2)

lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)              # Debye length squared (m^2)

# Altitude-dependent parameters (Analytic Fit)
alt_eff = np.array([135, 129, 122, 116, 110, 104, 97, 91]) # Effective altitudes (h in km)

# psi(h) fit
term1 = np.exp(-(alt_eff - 94.5) / 3.0) / (1.0 + np.exp((alt_eff - 120.0) / 15.0))
term2 = np.exp(-(alt_eff - 71.7) / 6.0) / (1.0 + np.exp(-(alt_eff - 120.0) / 15.0))
psi_h = term1 + term2                                   # Collision profile fit function (dimensionless)

Theta_0 = 1.35e-2                                       # Rounded mass-collision ratio constant (dimensionless)

kappa_i = Theta_0 / np.sqrt(psi_h)                      # Ion magnetization parameter array (dimensionless)
kappa_e = 1.0 / (Theta_0 * np.sqrt(psi_h))              # Electron magnetization parameter array (dimensionless)

nu_i = Omega_i / kappa_i                                # Ion-neutral collision frequency array (s^-1)
nu_e = Omega_e / kappa_e                                # Electron-neutral collision frequency array (s^-1)

kappa = kappa_e                                         # Alias for electron magnetization (dimensionless)
Psi_0 = 1 / (kappa_i * kappa_e)                         # Driving parameter Psi_0 array (dimensionless)

# Electron drift (V_d) is on the x-axis, so Electric Field (E_0) is on the y-axis.
a_i0_mag = 0#(e * E_0) / m_i                              # Acceleration vector magnitude (m/s^2)
a_i0_x = 0.0                                            # Acceleration vector x-component (m/s^2)
a_i0_y = a_i0_mag                                       # Acceleration vector y-component (m/s^2)

U0_mag = (e * E_0) / (m_i * nu_i)                       # Primary Ion Pedersen drift magnitude array (m/s)
U0_x = np.zeros_like(U0_mag)                            # Ion Pedersen drift x-component array (m/s)
U0_y = U0_mag                                           # Ion Pedersen drift y-component array (m/s)

beta_T = T_n / T_e                                      # Temperature ratio (dimensionless)
V_Tn = np.sqrt(T_n / m_i)                               # Neutral thermal velocity (m/s)

# Anisotropic Thermal Velocities
# Heat is primarily generated along the Pedersen drift axis (y-axis)
T_eff = T_n + (m_i * U0_mag**2) / 3.0                   # Effective temperature array (J)

anisotropy_ratio = 1                                    # Anisotropy ratio (dimensionless) assume an isotropic ion temperature initially
T_y = (3.0 * T_eff) / (1.0 + 2.0 * anisotropy_ratio)    # Effective temperature y-component array (J)
T_x = anisotropy_ratio * T_y                            # Effective temperature x-component array (J)

V_Tx = np.sqrt(T_x / m_i)                               # Thermal velocity x-component array (m/s)
V_Ty = np.sqrt(T_y / m_i)                               # Thermal velocity y-component array (m/s)
T_eff_i = m_i * (2 * V_Tx**2 + V_Ty**2) / 3.0           # Effective ion temperature array (J)