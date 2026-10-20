from root import *

def eq_246(omega, k_x, k_y, run):
    # constants.py stores nu_i, nu_e, kappa_e as per-altitude-run arrays
    nu_i_r = nu_i[run]
    nu_e_r = nu_e[run]
    kappa_r = kappa_e[run]

    k_sq = k_x**2 + k_y**2
    k_mag = np.sqrt(k_sq)

    omega_i_combo = omega + 1j * nu_i_r
    beta = omega_i_combo / (np.sqrt(2) * V_th * k_mag)

    term1_den = (-1j * omega * nu_e_r * (1 + kappa_r**2)
                 + (gamma * T_e / m_e) * k_sq
                 + 1j * nu_e_r * kappa_r * V_d * (kappa_r * k_x - k_y))
    term1 = omega_pe_squared / term1_den

    # Route through the piecewise Z to keep the wofz argument in its stable half-plane.
    z_beta = Z_tilde(beta, omega_i_combo)
    z_prime_beta = Z_prime_tilde(beta, omega_i_combo)

    term2_den = (2 * V_th**2 * k_sq
                 + np.sqrt(2) * 1j * nu_i_r * V_th * k_mag * z_beta)

    term2 = (omega_pi_squared / term2_den) * z_prime_beta

    return 1.0 - term1 + term2

for i in [4, 5]:
    pts = 100
    kx_vals = np.linspace(-15.0, 15.0, pts)
    ky_vals = np.linspace(-8.0, 8.0, pts)
    # Growth rates with the new constants are much larger than the old fixed 0-70 scale, so auto-scale instead.
    contour_plot(lambda omega, kx, ky, run=i: eq_246(omega, kx, ky, run), kx_vals, ky_vals, alt_eff[i], pts, run=i, positive_only=True)
