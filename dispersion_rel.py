"""
Dispersion relation with fluid electrons and the full kinetic ion term:

    D = k^2 lambda_D^2 (1 + chi_e + chi_i)

    chi_e  fluid electrons (isothermal, constant collision rate nu_e)
    chi_i  kinetic ions, Eq. (73) including the i k.a_i0 terms

For kinetic electrons see electron_kinetic.py; for the "old terms" version
(ions without k.a_i0) see dispersion_rel_reduced.py.
"""
from root import *

def dispersion_rel(omega, k_x, k_y, run, a_i0_enabled=True):
    # Electron drift: E x B along x plus the Pedersen drift along -E_0 (-y)
    V_de_x = (E_0 / B_0)* kappa_e[run]**2 / (1 + kappa_e[run]**2)
    V_de_y = -(E_0 / B_0) * kappa_e[run] / (1 + kappa_e[run]**2)

    # Ion Pedersen drift U_0, along E_0 (+y)
    U0_mag_curr = (e * E_0) / (m_i * nu_i[run])
    U0_x_curr = 0
    U0_y_curr = U0_mag_curr

    # Ion temperatures heated by friction (see constants.py)
    T_eff_i_curr = T_n + (m_i * U0_mag_curr**2) / 3.0

    T_y_curr = (3.0 * T_eff_i_curr) / (1.0 + 2.0 * anisotropy_ratio)
    T_x_curr = anisotropy_ratio * T_y_curr

    V_Tx_curr = np.sqrt(T_x_curr / m_i)
    V_Ty_curr = np.sqrt(T_y_curr / m_i)

    # Ion acceleration by E_0 (set to 0 for the a_i0 = 0 variant)
    a_i0_mag_curr = (e * E_0) / m_i if a_i0_enabled else 0.0
    a_i0_x_curr = 0.0
    a_i0_y_curr = a_i0_mag_curr

    lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)
    beta_T = T_n / T_e

    # Keep k^2 above 0 to avoid dividing by zero
    k_sq = np.maximum(k_x**2 + k_y**2, 1e-12)

    # --- ELECTRON TERM (chi_e) ---
    # Fluid electrons: continuity + momentum without inertia, isothermal
    # pressure (gamma = 1), constant collision rate nu_e.
    k_dot_V_de = k_x * V_de_x + k_y * V_de_y
    term1_den = (-1j * omega * nu_e[run] * (1 + kappa_e[run]**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e[run] * (1 + kappa_e[run]**2) * k_dot_V_de)

    chi_e = omega_pe_squared / term1_den

    # --- ION TERM (chi_i), Eq. (73) ---
    k_dot_a_i0 = k_x * a_i0_x_curr + k_y * a_i0_y_curr
    k_dot_U0 = k_x * U0_x_curr + k_y * U0_y_curr
    omega_i = omega + 1j * nu_i[run]

    # den1 uses the neutral temperature, den2 the heated ion temperatures
    den1 = k_sq * V_Tn**2 + 1j * k_dot_a_i0
    sqrt_den1 = np.sqrt(2 * den1)
    z1 = omega_i / sqrt_den1

    den2 = k_x**2 * V_Tx_curr**2 + k_y**2 * V_Ty_curr**2 + 1j * k_dot_a_i0
    sqrt_den2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U0) / sqrt_den2

    # Prefactor of chi_i. It uses den2, as in Eq. (15); Eq. (73) as printed
    # has den1 here, which is a typo (the paper notes this correction).
    # chi_i = I_ion / (k^2 lambda_D^2), hence lambda_D^2 in the denominator.
    chi_i_part1 = V_Tn**2 / (lambda_D_sq * beta_T * den2)
    Z1 = Z_tilde(z1, omega_i)
    Z2 = Z_tilde(z2, omega_i)

    part2 = 1.0 / (
        1.0 + (1j * nu_i[run] / sqrt_den1) * Z1
    )

    part3 = 1.0 + z2 * Z2

    chi_i = chi_i_part1 * part2 * part3

    # --- MULTIPLY BY k^2 lambda_D^2 ---
    # chi_e and chi_i are huge; multiplying by k^2 lambda_D^2 brings D to
    # order 1, which helps the root finder.
    scaling_factor = k_sq * lambda_D_sq

    return scaling_factor * (1.0 + chi_e + chi_i)


def dispersion_rel_ai0zero(omega, k_x, k_y, run):
    """dispersion_rel with a_i0 = 0, to compare a_i0 on vs off."""
    return dispersion_rel(omega, k_x, k_y, run, a_i0_enabled=False)

if __name__ == "__main__":
    # Quick test: solve and plot one grid. See main.py for the full options.
    import main

    main.run(dispersion_rel, "growth_contour", "growth_contour_116km.png", run=4)
