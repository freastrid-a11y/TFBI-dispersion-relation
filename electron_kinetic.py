"""
"New terms" dispersion relation: kinetic ions + kinetic electrons,
manuscript Eq. (72):

    D = I_ion + I_electron + k^2 lambda_D^2 = 0

    I_ion       ion response, Eq. (73)
    I_electron  electron response, Eq. (75), with alpha_1, alpha_0, q from
                Eq. (67) as written on page 38 (lines 1086-1088)

Watch out: the paper's k_x and k_y mean different things in the ion and
electron parts.

* Ion part, Eq. (73): x is along E_0, and k is perpendicular to B.
* Electron part, Eqs. (67), (75): "k_y" is the part of k perpendicular to
  B (k_perp), and "k_x" is the part parallel to B (k_par).

This repo uses B along +z, E_0 (and a_i0) along +y, V_d along +x.
So in the electron part:

    paper "k_y^2"     ->  k_perp^2 = k_x^2 + k_y^2   (code k_x, k_y)
    paper "k_x"       ->  k_par (0 here: waves perpendicular to B)
    paper "k . a_e0"  ->  k_y * a_e0                 (a_e0 along +y)

Other definitions used (paper Eqs. 29-30, lines 1079-1090):
    a_e0 = e E_0 / m_e      (along E_0, the paper's sign convention)
    V_Te^2 = T_e / m_e
    nu_0 = <nu_en> / 3      (because nu_en(V) = nu_0 V^2 / V_Te^2)
    Delta = k . V_e0 - omega
    q_0 = the 4th root of q^4 with Re q > 0 and Im q > 0
    Z(x) = i sqrt(pi) w(x)
"""

import numpy as np
from scipy.special import wofz

from constants import (
    B_0, E_0, T_e, T_n, V_Tn, e, eps_0, kappa_e, m_e, n_0, nu_e, Omega_e,
)
from root import Z_tilde

SQRT_PI = np.sqrt(np.pi)


def Z(x):
    """Plasma dispersion function Z(x) = i sqrt(pi) w(x)."""
    return 1j * SQRT_PI * wofz(x)


def first_quadrant_quartic_root(q4):
    """
    Return the 4th root of q4 that lies in the first quadrant (Re > 0, Im > 0).

    The four 4th roots are q4**(1/4) times 1, i, -1, -i. We take the one
    whose smaller coordinate, min(Re q, Im q), is largest: that is the root
    furthest inside the first quadrant. This still gives a sensible answer
    when a root sits exactly on an axis.
    """
    q4 = np.asarray(q4, dtype=complex)
    base = q4 ** 0.25                                  # one of the 4 roots
    roots = np.stack([base * (1j) ** n for n in range(4)])
    score = np.minimum(roots.real, roots.imag)
    idx = np.argmax(score, axis=0)
    return np.take_along_axis(roots, idx[None, ...], axis=0)[0]


def electron_params(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True):
    """
    Compute the pieces of the electron term, Eq. (67): alpha_1, alpha_0, q^4.

    Returns a dict with D1, alpha0, alpha1, q4, nu_0, Delta, V_Te_sq, a_e0
    and k_dot_a. a_e0_enabled=False sets a_e0 = 0 (the paper's NTA limit).
    """
    kappa = kappa_e[run]
    nu_mean = nu_e[run]                    # <nu_en>
    nu_0 = nu_mean / 3.0                   # Eq. (30)
    V_Te_sq = T_e / m_e

    # Electron drift: E x B along x plus the Pedersen drift along -E_0 (-y)
    V_e0_x = (E_0 / B_0) * kappa**2 / (1.0 + kappa**2)
    V_e0_y = -(E_0 / B_0) * kappa / (1.0 + kappa**2)
    Delta = k_x * V_e0_x + k_y * V_e0_y - omega   # frequency seen by the drifting electrons

    a_e0 = (e * E_0 / m_e) if a_e0_enabled else 0.0
    k_dot_a = k_y * a_e0                   # a_e0 is along +y (E_0)

    k_perp_sq = k_x**2 + k_y**2            # the paper's "k_y^2" in the electron part
    D1 = k_perp_sq * V_Te_sq - 1j * k_dot_a              # k_perp^2 V_Te^2 - i k.a_e0

    # Page 38 (lines 1086-1088): only alpha_1 has an extra 2 in the
    # denominator. alpha_0 and q^4 are divided by 4 nu_0^2 D1.
    alpha1 = 5j * k_dot_a / (2.0 * D1)                                    # (67a)
    alpha0 = k_par**2 * V_Te_sq * Omega_e**2 / (4.0 * nu_0**2 * D1)       # (67b), 0 when k_par = 0
    q4 = -((3j * nu_0 * Delta + k_par**2 * V_Te_sq) * Omega_e**2
           / (4.0 * nu_0**2 * D1))                                        # (67c)

    return dict(D1=D1, alpha0=alpha0, alpha1=alpha1, q4=q4,
                nu_0=nu_0, Delta=Delta, V_Te_sq=V_Te_sq, a_e0=a_e0,
                k_dot_a=k_dot_a)


def I_electron(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True):
    """Kinetic electron term I_electron, Eq. (75)."""
    p = electron_params(omega, k_x, k_y, run, k_par, a_e0_enabled)
    q0 = first_quadrant_quartic_root(p["q4"])
    q4, a0, a1 = p["q4"], p["alpha0"], p["alpha1"]

    bracket = (1.0 + 2.0 * a1
               + (q4 + a0 + a1 * q0**2) / q0 * Z(q0)
               - 1j * (q4 + a0 - a1 * q0**2) / q0 * Z(1j * q0))
    return bracket


def I_ion(omega, k_x, k_y, run, a_i0_enabled=True):
    """
    Kinetic ion term I_ion, Eq. (73), with a_i0 and U_0 along +y.

    Same physics as chi_i in dispersion_rel.py, multiplied by
    k^2 lambda_D^2. a_i0_enabled=False sets a_i0 = 0 (the drift U_0 and
    the ion heating are kept).
    """
    from constants import (anisotropy_ratio, m_i, nu_i)

    # Ion drift and heated ion temperatures, see constants.py
    nu = nu_i[run]
    U0 = e * E_0 / (m_i * nu)
    T_eff = T_n + m_i * U0**2 / 3.0
    T_par = 3.0 * T_eff / (1.0 + 2.0 * anisotropy_ratio)  # along E_0 (y)
    T_perp = anisotropy_ratio * T_par                      # across E_0 (x)
    V_Tx_sq, V_Ty_sq = T_perp / m_i, T_par / m_i

    a_i0 = (e * E_0 / m_i) if a_i0_enabled else 0.0
    k_dot_a = k_y * a_i0
    k_dot_U0 = k_y * U0
    k_sq = k_x**2 + k_y**2
    beta_T = T_n / T_e
    omega_i = omega + 1j * nu

    # The two denominators of Eq. (73): den1 uses the neutral temperature,
    # den2 the heated ion temperatures.
    den1 = k_sq * V_Tn**2 + 1j * k_dot_a
    s1 = np.sqrt(2.0 * den1)
    den2 = k_x**2 * V_Tx_sq + k_y**2 * V_Ty_sq + 1j * k_dot_a
    s2 = np.sqrt(2.0 * den2)

    # The prefactor uses den2, as in Eq. (15). Eq. (73) as printed has den1
    # here, which is a typo (the paper notes this correction).
    first = k_sq * V_Tn**2 / (beta_T * den2)
    second = 1.0 / (1.0 + 1j * nu / s1 * Z_tilde(omega_i / s1, omega_i))
    z2 = (omega_i - k_dot_U0) / s2
    third = 1.0 + z2 * Z_tilde(z2, omega_i)
    return first * second * third


def eq_72_kinetic(omega, k_x, k_y, run, a_i0_enabled=True, a_e0_enabled=True,
                  k_par=0.0):
    """
    Full kinetic dispersion function, Eq. (72):
        D = I_ion + I_electron + k^2 lambda_D^2

    Has the same overall size as dispersion_rel (which returns
    k^2 lambda_D^2 (1 + chi_e + chi_i)), so root.find_roots can solve both.
    """
    # k^2 is kept above 1e-12 to avoid exactly 0
    k_sq = max(k_x**2 + k_y**2, 1e-12) if np.isscalar(k_x) else np.maximum(k_x**2 + k_y**2, 1e-12)
    lambda_D_sq = eps_0 * T_e / (n_0 * e**2)
    return (I_ion(omega, k_x, k_y, run, a_i0_enabled)
            + I_electron(omega, k_x, k_y, run, k_par, a_e0_enabled)
            + k_sq * lambda_D_sq)


def eq_72_kinetic_ae0zero(omega, k_x, k_y, run):
    """Eq. (72) with a_e0 = 0 (the paper's NTA electrons); a_i0 kept."""
    return eq_72_kinetic(omega, k_x, k_y, run, a_e0_enabled=False)


def eq_72_kinetic_ai0zero(omega, k_x, k_y, run):
    """Eq. (72) with a_i0 = 0; a_e0 kept."""
    return eq_72_kinetic(omega, k_x, k_y, run, a_i0_enabled=False)
