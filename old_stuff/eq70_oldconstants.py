import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from oldroot import *
import matplotlib.pyplot as plt

def eq_70(omega, k_x, k_y):
    # Static electron drift only (V_d), matching eq 2.46's assumptions
    V_de_x = V_d * kappa**2 / (1 + kappa**2)
    V_de_y = -V_d * kappa / (1 + kappa**2)

    # No ion Pedersen drift and no ion acceleration (matches eq 2.46)
    U0_x_curr = 0.0
    U0_y_curr = 0.0
    a_i0_x_curr = 0.0
    a_i0_y_curr = 0.0

    # Isotropic ion/electron temperatures, no Joule heating (matches eq 2.46)
    V_Tx_curr = V_th
    V_Ty_curr = V_th

    k_sq = max(k_x**2 + k_y**2, 1e-12)
    growth_rate = np.imag(omega)

    k_dot_V_de = k_x * V_de_x + k_y * V_de_y
    term1_den = (-1j * omega * nu_e * (1 + kappa**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e * (1 + kappa**2) * k_dot_V_de)
    norm_I_electron = -omega_pe_squared / term1_den

    k_dot_a_i0 = k_x * a_i0_x_curr + k_y * a_i0_y_curr
    k_dot_U0 = k_x * U0_x_curr + k_y * U0_y_curr
    omega_i = omega + 1j * nu_i

    den1 = k_sq * V_th**2 + 1j * k_dot_a_i0
    sqrt_den1 = np.sqrt(2 * den1)
    z1 = omega_i / sqrt_den1

    den2 = k_x**2 * V_Tx_curr**2 + k_y**2 * V_Ty_curr**2 + 1j * k_dot_a_i0
    sqrt_den2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U0) / sqrt_den2

    norm_part1 = V_th**2 / (lambda_D_sq * beta_T * den1)
    part2 = 1.0 / (1.0 + (1j * nu_i / sqrt_den1) * Z_tilde(z1, growth_rate, nu_i))
    part3 = 1.0 + ((omega_i - k_dot_U0) / sqrt_den2) * Z_tilde(z2, growth_rate, nu_i)
    norm_I_ion = norm_part1 * part2 * part3

    # ion term must be subtracted to cancel against norm_I_electron in the k->0 (fluid) limit
    return 1.0 + norm_I_electron - norm_I_ion

# Same grids/levels/drift velocities as eq246.py, for a direct visual comparison.
kx_vals = np.linspace(-15.0, 15.0, 100)
ky_vals = np.linspace(-8.0, 8.0, 100)
custom_levels = np.linspace(0, 70, 60)
contour_plot(eq_70, kx_vals, ky_vals, custom_levels, "new_eq70_Vd_600")

kx_vals = np.linspace(-25.0, 25.0, 100)
ky_vals = np.linspace(-25.0, 25.0, 100)
custom_levels = np.linspace(0, 650, 60)
contour_plot(eq_70, kx_vals, ky_vals, custom_levels, "new_eq70_Vd_800", new_V_d=800)
