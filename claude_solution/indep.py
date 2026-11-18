"""
Independent implementation of Dimant et al. kinetic TFBI dispersion relation, Eq. (72),
written from the manuscript (not from the repo code), with the ion prefactor
denominator taken from Eq. (15) (bi-Maxwellian k_x^2 V_Tx^2 + k_y^2 V_Ty^2 + i k.a_i0).

Frame used here (same as the repo so grids can be compared point by point):
  E0 along +y, B along +z, E x B drift V0 along +x.
  Ion part: paper's "x" (E0 direction) -> our y.
  Electron part: paper's k_y -> k_perp, paper's k_x -> k_par (= 0 here).
"""
import numpy as np
from scipy.special import wofz

# ---- physical parameters (same plasma as the repo's constants.py) ----
mi, me, e, eps0, kB = 5e-26, 9.1094e-31, 1.6e-19, 8.854e-12, 1.380649e-23
delta = 3e-3
B = 5e-5
Tn = 300 * kB
n0 = 1e11
alts = np.array([135, 129, 122, 116, 110, 104, 97, 91])
psi = (np.exp(-(alts - 94.5) / 3) / (1 + np.exp((alts - 120) / 15))
       + np.exp(-(alts - 71.7) / 6) / (1 + np.exp(-(alts - 120) / 15)))
Th0 = 1.35e-2
kap_i = Th0 / np.sqrt(psi)
kap_e = 1 / (Th0 * np.sqrt(psi))
Om_i, Om_e = e * B / mi, e * B / me
nu_in = Om_i / kap_i
nu_en = Om_e / kap_e          # <nu_en>


def params(run, E_scale=1.0, ae0_sign=+1.0):
    Vd = 2100.0 * E_scale
    E0 = Vd * B
    Te = Tn + 2 * me * E0**2 / (3 * delta * B**2)
    ke = kap_e[run]
    U0 = e * E0 / (mi * nu_in[run])
    Teff = Tn + mi * U0**2 / 3          # isotropic bi-Maxwellian (anisotropy = 1)
    return dict(
        nu_i=nu_in[run], nu0=nu_en[run] / 3, Te=Te, VTe2=Te / me, VTn2=Tn / mi,
        VTx2=Teff / mi, VTy2=Teff / mi, U0=U0, ai0=e * E0 / mi,
        ae0=ae0_sign * e * E0 / me, betaT=Tn / Te, lD2=eps0 * Te / (n0 * e**2),
        Vex=Vd * ke**2 / (1 + ke**2), Vey=-Vd * ke / (1 + ke**2),
    )


def Z(x):
    return 1j * np.sqrt(np.pi) * wofz(x)


def Zt(x, wi):   # Eq. (74)
    return Z(x) if wi.imag > 0 else -Z(-x)


def q_first_quadrant(q4):
    r = q4 ** 0.25 * np.array([1, 1j, -1, -1j])
    return r[np.argmax(np.minimum(r.real, r.imag))]


def D(w, kx, ky, p):
    """I_ion + I_electron + k^2 lambda_D^2   (Eq. 72)."""
    k2 = kx * kx + ky * ky
    # ---------- ions, Eq. (15) with Z = i sqrt(pi) w ----------
    wi = w + 1j * p['nu_i']
    ka = ky * p['ai0']                       # a_i0 along E0 = +y
    d1 = k2 * p['VTn2'] + 1j * ka
    d2 = kx**2 * p['VTx2'] + ky**2 * p['VTy2'] + 1j * ka   # paper x <-> our y; isotropic anyway
    s1, s2 = np.sqrt(2 * d1), np.sqrt(2 * d2)
    z2 = (wi - ky * p['U0']) / s2
    Iion = (k2 * p['VTn2'] / (p['betaT'] * d2)
            * (1 + z2 * Zt(z2, wi)) / (1 + 1j * p['nu_i'] / s1 * Zt(wi / s1, wi)))
    # ---------- electrons, Eqs. (67), (69)/(75), k_par = 0 ----------
    kae = ky * p['ae0']
    VTe2, nu0 = p['VTe2'], p['nu0']
    A = k2 * VTe2 - 2j * kae
    Bd = A - p['ae0']**2 / VTe2
    Dw = kx * p['Vex'] + ky * p['Vey'] - w
    a1 = 5j * kae / A
    a0 = 0.0
    q4 = -(3j * nu0 * Dw) * Om_e**2 / (4 * nu0**2 * Bd)
    q = q_first_quadrant(q4)
    Iel = A / Bd * (1 + 2 * a1 + (q4 + a0 + a1 * q * q) / q * Z(q)
                    - 1j * (q4 + a0 - a1 * q * q) / q * Z(1j * q))
    return Iion + Iel + k2 * p['lD2']