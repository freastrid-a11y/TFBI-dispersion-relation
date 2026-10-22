from root import *

def eq_73(omega, k_x, k_y, run):
    # Recalculate dynamic constants based on the current run
    V_de_x = (E_0 / B_0) * kappa_e[run]**2 / (1 + kappa_e[run]**2)
    V_de_y = -(E_0 / B_0) * kappa_e[run] / (1 + kappa_e[run]**2)

    # Heat is primarily generated along the Pedersen drift axis (y-axis)
    U0_mag_curr = (e * E_0) / (m_i * nu_i[run])
    U0_x_curr = 0
    U0_y_curr = U0_mag_curr

    T_eff_i_curr = T_n + (m_i * U0_mag_curr**2) / 3.0

    T_y_curr = (3.0 * T_eff_i_curr) / (1.0 + 2.0 * anisotropy_ratio)
    T_x_curr = anisotropy_ratio * T_y_curr

    V_Tx_curr = np.sqrt(T_x_curr / m_i)
    V_Ty_curr = np.sqrt(T_y_curr / m_i)

    a_i0_mag_curr = (e * E_0) / m_i
    a_i0_x_curr = 0.0
    a_i0_y_curr = a_i0_mag_curr

    lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)
    beta_T = T_n / T_e

    # Ensures k_sq is never strictly 0 to avoid division by zero
    k_sq = np.maximum(k_x**2 + k_y**2, 1e-12)

    # --- ELECTRON TERM (chi_e) ---
    k_dot_V_de = k_x * V_de_x + k_y * V_de_y
    term1_den = (-1j * omega * nu_e[run] * (1 + kappa_e[run]**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e[run] * (1 + kappa_e[run]**2) * k_dot_V_de)
    
    # This exactly matches the fluid electron susceptibility term
    chi_e = omega_pe_squared / term1_den

    # --- ION TERM (chi_i) ---
    k_dot_a_i0 = k_x * a_i0_x_curr + k_y * a_i0_y_curr
    k_dot_U0 = k_x * U0_x_curr + k_y * U0_y_curr
    omega_i = omega + 1j * nu_i[run]

    den1 = k_sq * V_Tn**2 + 1j * k_dot_a_i0
    sqrt_den1 = np.sqrt(2 * den1)
    z1 = omega_i / sqrt_den1

    den2 = k_x**2 * V_Tx_curr**2 + k_y**2 * V_Ty_curr**2 + 1j * k_dot_a_i0
    sqrt_den2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U0) / sqrt_den2

    # The denominator natively includes lambda_D_sq but no k^2 in the numerator, 
    # effectively evaluating the unnormalized ion susceptibility (chi_i)
    chi_i_part1 = V_Tn**2 / (lambda_D_sq * beta_T * den1)
    part2 = 1.0 / (1.0 + (1j * nu_i[run] / sqrt_den1) * Z_tilde(z1))
    part3 = 1.0 + ((omega_i - k_dot_U0) / sqrt_den2) * Z_tilde(z2)
    
    chi_i = chi_i_part1 * part2 * part3

    # --- NUMERICALLY SCALED DISPERSION RELATION ---
    # Multiplying by (k^2 * lambda_D^2) brings all extremely large susceptibility 
    # terms down to O(1), allowing the Newton solver to converge smoothly.
    scaling_factor = k_sq * lambda_D_sq
    
    return scaling_factor * (1.0 + chi_e + chi_i)

if __name__ == "__main__":
    for i in [4]:
        pts = 100
        kx_vals = np.linspace(-100, 100, pts)
        ky_vals = np.linspace(-100, 100, pts)

        contour_plot(
            lambda omega, kx, ky, run=i: eq_73(omega, kx, ky, run),
            kx_vals,
            ky_vals,
            alt_eff[i],
            run=i,
        )