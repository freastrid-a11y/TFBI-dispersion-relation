from oldroot import *

def eq_246(omega, k_x, k_y):
    k_sq = k_x**2 + k_y**2
    k_mag = np.sqrt(k_sq)

    beta = (omega + 1j * nu_i) / (np.sqrt(2) * V_th * k_mag)

    term1_den = (-1j * omega * nu_e * (1 + kappa**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e * kappa * V_d * (kappa * k_x - k_y))
    term1 = omega_pe_squared / term1_den

    z_beta = Z(beta)
    z_prime_beta = Z_prime(beta)

    term2_den = (2 * V_th**2 * k_sq
                 + np.sqrt(2) * 1j * nu_i * V_th * k_mag * z_beta)

    term2 = (omega_pi_squared / term2_den) * z_prime_beta

    return 1.0 - term1 + term2

kx_vals = np.linspace(-15.0, 15.0, 100)
ky_vals = np.linspace(-8.0, 8.0, 100)
custom_levels = np.linspace(0, 70, 60)
contour_plot(eq_246, kx_vals, ky_vals, custom_levels, "new_params_eq246_Vd_600")

kx_vals = np.linspace(-25.0, 25.0, 100)
ky_vals = np.linspace(-25.0, 25.0, 100)
custom_levels = np.linspace(0, 650, 60)
contour_plot(eq_246, kx_vals, ky_vals, custom_levels, "new_params_eq246_Vd_800", new_V_d=800)