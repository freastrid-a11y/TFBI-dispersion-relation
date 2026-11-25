"""
Physical constants and plasma parameters used by every solver file.

Frame used throughout the repo:
    B along +z,  E_0 along +y,  E x B drift V_d along +x.

The electric field can be scaled with the E_SCALE environment variable
(E_SCALE=2 doubles it). Everything that depends on the field (T_e, drifts,
ion heating, ...) is computed below from E_0, so it all scales together.

Arrays hold one value per altitude in alt_eff; `run` picks the altitude
(run 4 = 110 km, run 5 = 104 km).
"""
import os
import numpy as np

# --- Particle masses and fundamental constants ---
m_i = 5e-26                                             # Ion mass (kg)
m_e = 9.1094e-31                                        # Electron mass (kg)
e = 1.6e-19                                             # Elementary charge (C)
eps_0 = 8.854e-12                                       # Vacuum permittivity (F/m)
delta_en = 3e-3                                         # Fraction of energy an electron loses per collision

# --- Fields, temperatures, gyrofrequencies ---
B_0 = 5e-05                                             # Magnetic field (T)
E_SCALE = float(os.environ.get("E_SCALE", "1.0"))       # Field multiplier, set by replot_variants.py
V_d = 2100 * E_SCALE                                    # E x B drift speed (m/s)
E_0 = V_d * B_0                                         # Electric field (V/m): 0.105 V/m when E_SCALE = 1
gamma = 1.0                                             # Electron pressure: 1 = isothermal
k_B = 1.380649e-23                                      # Boltzmann constant (J/K)
T_n = 300 * k_B                                         # Neutral (and background ion) temperature (J)
T_e = T_n + (2 * m_e * E_0**2) / (3 * delta_en * B_0**2)# Electron temperature, heated by friction (J)
Omega_i = (e * B_0) / m_i                               # Ion gyrofrequency (rad/s)
Omega_e = (e * B_0) / m_e                               # Electron gyrofrequency (rad/s)

# --- Density ---
n_0 = 10e10                                             # Plasma density (m^-3); note 10e10 = 1e11
omega_pe_squared = n_0*e**2/(eps_0*m_e)                 # Electron plasma frequency squared (rad^2/s^2)
omega_pi_squared = n_0*e**2/(eps_0*m_i)                 # Ion plasma frequency squared (rad^2/s^2)

lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)              # Debye length squared (m^2), uses the heated T_e

# --- Altitude profile ---
alt_eff = np.array([135, 129, 122, 116, 110, 104, 97, 91]) # Altitudes (km); `run` indexes this array

# psi(h): fitted collision profile. It sets how magnetized ions and electrons are.
term1 = np.exp(-(alt_eff - 94.5) / 3.0) / (1.0 + np.exp((alt_eff - 120.0) / 15.0))
term2 = np.exp(-(alt_eff - 71.7) / 6.0) / (1.0 + np.exp(-(alt_eff - 120.0) / 15.0))
psi_h = term1 + term2                                   # Collision profile (dimensionless)

Theta_0 = 1.35e-2                                       # Constant in the profile fit

# kappa = gyrofrequency / collision frequency. kappa >> 1: magnetized; kappa << 1: unmagnetized.
kappa_i = Theta_0 / np.sqrt(psi_h)                      # Ion kappa (per altitude)
kappa_e = 1.0 / (Theta_0 * np.sqrt(psi_h))              # Electron kappa (per altitude)

nu_i = Omega_i / kappa_i                                # Ion-neutral collision frequency (s^-1)
nu_e = Omega_e / kappa_e                                # Electron-neutral collision frequency <nu_en> (s^-1)

Psi_0 = 1 / (kappa_i * kappa_e)                         # Psi_0 = nu_e nu_i / (Omega_e Omega_i), the FB parameter

# --- Ion drift (E_0 is along +y) ---
U0_mag = (e * E_0) / (m_i * nu_i)                       # Ion Pedersen drift speed along E_0 (m/s)
U0_y = U0_mag                                           # Its y-component (all of it) (m/s)

V_Tn = np.sqrt(T_n / m_i)                               # Ion thermal speed at T_n (m/s)

# --- Ion heating by friction with neutrals ---
# The drifting ions are heated: T_eff = T_n + m_i U0^2 / 3 (average over directions).
T_eff = T_n + (m_i * U0_mag**2) / 3.0                   # Average heated ion temperature (J)

# How the heating is shared between directions: T_x = anisotropy_ratio * T_y.
# 1 means the same temperature in every direction (isotropic).
anisotropy_ratio = 1                                    # T_x / T_y
T_y = (3.0 * T_eff) / (1.0 + 2.0 * anisotropy_ratio)    # Ion temperature along E_0 (J)
T_x = anisotropy_ratio * T_y                            # Ion temperature across E_0 (J)

V_Tx = np.sqrt(T_x / m_i)                               # Ion thermal speed across E_0 (m/s)
V_Ty = np.sqrt(T_y / m_i)                               # Ion thermal speed along E_0 (m/s)
