"""
Fully kinetic electron contribution to the TFBI dispersion relation,
manuscript Eqs. (67a-c), (70) and (75), and the full dispersion relation
Eq. (72):  I_ion + I_electron + k^2 lambda_D^2 = 0.

NOTATION WARNING (important)
----------------------------
The manuscript uses two different meanings of "k_x, k_y":

* Ion part, Eq. (73): x is along a_i0 (i.e. E_0), k lies in the plane
  perpendicular to B (k_z = 0).
* Electron part, Eqs. (23)-(25), (43)-(46), (67), (75): "k_y" is the
  magnitude of the component of k perpendicular to B (k_perp), while
  "k_x" is the component PARALLEL to B (k_par).  This follows from
  Eq. (23)/(25), where the k_x^2 term is the field-aligned diffusion
  V^2/nu_en, and from line 811: k_perp . a_e0 = k_y a_e0 cos(phi).

This repository uses the frame:  V_d (E x B drift) along +x,
E_0 (and a_i0) along +y, B along +z.  So in this module

    paper "k_y^2"      ->  k_perp^2 = k_x^2 + k_y^2      (code k_x, k_y)
    paper "k_x"        ->  k_par    (0 for strictly perpendicular waves)
    paper "k . a_e0"   ->  k_y * a_e0                    (a_e0 along +y = E_0)

Other definitions (paper lines 1058-1067, Eqs. 29-30):
    a_e0 = e E_0 / m_e   (the paper's sign convention: parallel to E_0)
    V_Te^2 = T_e / m_e
    nu_0 = <nu_en>/3   since nu_en(V) = nu_0 V^2/V_Te^2  and <nu_en> = 3 nu_0
    Delta_wk = k . V_e0 - omega
    q_0 = the quartic root of Eq. (67c) with Re q > 0 and Im q > 0
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
    """Ordinary plasma dispersion function Z(x) = i sqrt(pi) w(x)."""
    return 1j * SQRT_PI * wofz(x)


def first_quadrant_quartic_root(q4):
    """Return the 4th root of q4 lying in the open first quadrant (Re>0, Im>0).

    The four roots are q4**(1/4) * i**n.  Exactly one of them lies inside
    the first quadrant unless q4 is real and positive/negative-real with a
    root on an axis; we pick the root that is "most inside" the quadrant,
    i.e. maximises min(Re q, Im q), which is unique and continuous.
    """
    q4 = np.asarray(q4, dtype=complex)
    base = q4 ** 0.25                                  # principal root
    roots = np.stack([base * (1j) ** n for n in range(4)])
    score = np.minimum(roots.real, roots.imag)
    idx = np.argmax(score, axis=0)
    return np.take_along_axis(roots, idx[None, ...], axis=0)[0]


def electron_params(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True):
    """Return the building blocks D1, D2, alpha0, alpha1, q4 of Eqs. (67)."""
    kappa = kappa_e[run]
    nu_mean = nu_e[run]                    # <nu_en>
    nu_0 = nu_mean / 3.0                   # Eq. (30)
    V_Te_sq = T_e / m_e

    # Mean electron drift (Hall + Pedersen) in this repo's frame
    V_e0_x = (E_0 / B_0) * kappa**2 / (1.0 + kappa**2)
    V_e0_y = -(E_0 / B_0) * kappa / (1.0 + kappa**2)
    Delta = k_x * V_e0_x + k_y * V_e0_y - omega

    a_e0 = (e * E_0 / m_e) if a_e0_enabled else 0.0
    k_dot_a = k_y * a_e0                   # a_e0 along +y (E_0 direction)

    k_perp_sq = k_x**2 + k_y**2            # paper's "k_y^2" in the electron part
    D1 = k_perp_sq * V_Te_sq - 2j * k_dot_a              # k_y^2 V_Te^2 - 2i k.a_e0
    D2 = D1 - a_e0**2 / V_Te_sq                           # ... - a_e0^2/V_Te^2

    alpha1 = 5j * k_dot_a / D1                                            # (67a)
    alpha0 = k_par**2 * V_Te_sq * Omega_e**2 / (4.0 * nu_0**2 * D1)       # (67b)
    q4 = -((3j * nu_0 * Delta + k_par**2 * V_Te_sq) * Omega_e**2
           / (4.0 * nu_0**2 * D2))                                        # (67c)

    return dict(D1=D1, D2=D2, alpha0=alpha0, alpha1=alpha1, q4=q4,
                nu_0=nu_0, Delta=Delta, V_Te_sq=V_Te_sq, a_e0=a_e0,
                k_dot_a=k_dot_a)


def I_electron(omega, k_x, k_y, run, k_par=0.0, a_e0_enabled=True):
    """Kinetic electron contribution I_electron, manuscript Eq. (75)."""
    p = electron_params(omega, k_x, k_y, run, k_par, a_e0_enabled)
    q0 = first_quadrant_quartic_root(p["q4"])
    q4, a0, a1 = p["q4"], p["alpha0"], p["alpha1"]

    bracket = (1.0 + 2.0 * a1
               + (q4 + a0 + a1 * q0**2) / q0 * Z(q0)
               - 1j * (q4 + a0 - a1 * q0**2) / q0 * Z(1j * q0))
    return p["D1"] / p["D2"] * bracket


def I_ion(omega, k_x, k_y, run, a_i0_enabled=True):
    """Ion contribution I_ion, manuscript Eq. (73), in this repo's frame
    (a_i0 and U_0 along +y).  Identical physics to dispersion_rel's chi_i,
    multiplied by k^2 lambda_D^2."""
    from constants import (anisotropy_ratio, m_i, nu_i)

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

    den1 = k_sq * V_Tn**2 + 1j * k_dot_a
    s1 = np.sqrt(2.0 * den1)
    den2 = k_x**2 * V_Tx_sq + k_y**2 * V_Ty_sq + 1j * k_dot_a
    s2 = np.sqrt(2.0 * den2)

    first = k_sq * V_Tn**2 / (beta_T * den1)
    second = 1.0 / (1.0 + 1j * nu / s1 * Z_tilde(omega_i / s1, omega_i))
    z2 = (omega_i - k_dot_U0) / s2
    third = 1.0 + z2 * Z_tilde(z2, omega_i)
    return first * second * third


def eq_72_kinetic(omega, k_x, k_y, run, a_i0_enabled=True, a_e0_enabled=True,
                  k_par=0.0):
    """Full kinetic dispersion function, Eq. (72):
        I_ion + I_electron + k^2 lambda_D^2.

    Same overall normalisation as dispersion_rel (as dispersion_rel.dispersion_rel, which returns
    k^2 lambda_D^2 (1 + chi_e + chi_i)), so root.find_roots works unchanged.
    """
    k_sq = max(k_x**2 + k_y**2, 1e-12) if np.isscalar(k_x) else np.maximum(k_x**2 + k_y**2, 1e-12)
    lambda_D_sq = eps_0 * T_e / (n_0 * e**2)
    return (I_ion(omega, k_x, k_y, run, a_i0_enabled)
            + I_electron(omega, k_x, k_y, run, k_par, a_e0_enabled)
            + k_sq * lambda_D_sq)


def eq_72_kinetic_ae0zero(omega, k_x, k_y, run):
    """Eq. (72) in the NTA limit for electrons (a_e0 -> 0), a_i0 kept."""
    return eq_72_kinetic(omega, k_x, k_y, run, a_e0_enabled=False)


def eq_72_kinetic_ai0zero(omega, k_x, k_y, run):
    """Eq. (72) with a_i0 -> 0, kinetic electrons with a_e0 kept."""
    return eq_72_kinetic(omega, k_x, k_y, run, a_i0_enabled=False)