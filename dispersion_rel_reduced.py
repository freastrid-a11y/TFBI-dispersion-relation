"""
"Old terms" dispersion relation:

    D = k^2 lambda_D^2 (1 + chi_e + chi_i)

    chi_e  fluid electrons, the same as in dispersion_rel.py
    chi_i  kinetic ions, Eq. (73) with the i k.a_i0 terms removed from
           den1 and den2

The field still enters the ions through the drift U_0 and the ion heating;
only the explicit i k.a_i0 terms are dropped.

Same inputs and outputs as dispersion_rel.dispersion_rel, so root.find_roots
and main.run work with it unchanged.
"""

from root import *


def dispersion_rel_reduced(omega, k_x, k_y, run):
    # Electron drift: E x B along x plus the Pedersen drift along -E_0 (-y)
    V_de_x = (E_0 / B_0) * kappa_e[run]**2 / (1 + kappa_e[run]**2)
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

    lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)
    beta_T = T_n / T_e

    # Keep k^2 above 0 to avoid dividing by zero
    k_sq = np.maximum(k_x**2 + k_y**2, 1e-12)

    # --- ELECTRON TERM (chi_e) ---
    # Same fluid electrons as dispersion_rel.py: isothermal, constant nu_e.
    k_dot_V_de = k_x * V_de_x + k_y * V_de_y
    term1_den = (-1j * omega * nu_e[run] * (1 + kappa_e[run]**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e[run] * (1 + kappa_e[run]**2) * k_dot_V_de)

    chi_e = omega_pe_squared / term1_den

    # --- ION TERM (chi_i) WITHOUT i k.a_i0 ---
    # Same as dispersion_rel.py, but den1 and den2 have no i k.a_i0 term.
    k_dot_U0 = k_x * U0_x_curr + k_y * U0_y_curr
    omega_i = omega + 1j * nu_i[run]

    den1 = k_sq * V_Tn**2
    sqrt_den1 = np.sqrt(2 * den1)
    z1 = omega_i / sqrt_den1

    den2 = k_x**2 * V_Tx_curr**2 + k_y**2 * V_Ty_curr**2
    sqrt_den2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U0) / sqrt_den2

    # Prefactor uses den2, as in Eq. (15); see dispersion_rel.py.
    chi_i_part1 = V_Tn**2 / (lambda_D_sq * beta_T * den2)
    Z1 = Z_tilde(z1, omega_i)
    Z2 = Z_tilde(z2, omega_i)

    part2 = 1.0 / (
        1.0 + (1j * nu_i[run] / sqrt_den1) * Z1
    )

    part3 = 1.0 + z2 * Z2

    chi_i = chi_i_part1 * part2 * part3

    # --- MULTIPLY BY k^2 lambda_D^2 ---
    # Brings D to order 1 for the root finder (see dispersion_rel.py).
    scaling_factor = k_sq * lambda_D_sq

    return scaling_factor * (1.0 + chi_e + chi_i)


if __name__ == "__main__":
    # Quick test: solve and plot one grid. See main.py for the full options.
    import main

    main.run(dispersion_rel_reduced, "growth_contour",
             "growth_contour_reduced_110km.png", run=4)
