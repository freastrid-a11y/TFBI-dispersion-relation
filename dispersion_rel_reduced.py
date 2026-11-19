"""
Reduced TFBI dispersion relation.

Same fluid electron susceptibility (chi_e) as dispersion_rel.py, but with the
ion term reduced by dropping the i k.a_i0 terms from the denominators
den1 and den2.  Since a_i0 enters the fluid-electron/kinetic-ion relation
only through those i k.a_i0 terms, this is the NTA (a_i0 -> 0) limit of
the ion susceptibility.

Interface matches dispersion_rel.dispersion_rel, so root.find_roots and
main.run work unchanged.
"""

from root import *


def dispersion_rel_reduced(omega, k_x, k_y, run):
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

    lambda_D_sq = (eps_0 * T_e) / (n_0 * e**2)
    beta_T = T_n / T_e

    # Ensures k_sq is never strictly 0 to avoid division by zero
    k_sq = np.maximum(k_x**2 + k_y**2, 1e-12)

    # --- ELECTRON TERM (chi_e) ---
    # Identical to the fluid electron term in dispersion_rel.py
    k_dot_V_de = k_x * V_de_x + k_y * V_de_y
    term1_den = (-1j * omega * nu_e[run] * (1 + kappa_e[run]**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e[run] * (1 + kappa_e[run]**2) * k_dot_V_de)

    chi_e = omega_pe_squared / term1_den

    # --- REDUCED ION TERM (chi_i) ---
    # Same as dispersion_rel.py but with the i k.a_i0 term dropped from
    # den1 and den2.
    k_dot_U0 = k_x * U0_x_curr + k_y * U0_y_curr
    omega_i = omega + 1j * nu_i[run]

    den1 = k_sq * V_Tn**2
    sqrt_den1 = np.sqrt(2 * den1)
    z1 = omega_i / sqrt_den1

    den2 = k_x**2 * V_Tx_curr**2 + k_y**2 * V_Ty_curr**2
    sqrt_den2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U0) / sqrt_den2

    # NOTE: prefactor uses den2 (bi-Maxwellian temperatures), as in manuscript
    # Eq. (15); see dispersion_rel.py for why Eq. (73)'s den1 is a typo.
    chi_i_part1 = V_Tn**2 / (lambda_D_sq * beta_T * den2)
    Z1 = Z_tilde(z1, omega_i)
    Z2 = Z_tilde(z2, omega_i)

    part2 = 1.0 / (
        1.0 + (1j * nu_i[run] / sqrt_den1) * Z1
    )

    part3 = 1.0 + z2 * Z2

    chi_i = chi_i_part1 * part2 * part3

    # --- NUMERICALLY SCALED DISPERSION RELATION ---
    scaling_factor = k_sq * lambda_D_sq

    return scaling_factor * (1.0 + chi_e + chi_i)


if __name__ == "__main__":
    # See main.py for running a method (e.g. growth_contour) against this
    # equation with cached CSV grid data.
    import main

    main.run(dispersion_rel_reduced, "growth_contour",
             "growth_contour_reduced_110km.png", run=4)
