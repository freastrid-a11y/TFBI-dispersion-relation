"""
Kinetic (eq. 70, simplified) vs fluid dispersion test bench.

Setup:  put this file, constants.py (the TEST one) and a COPY of root.py in the
same folder, then run  python eq70_test.py
Change the switches in constants.FLAGS (or the STAGES list below) one at a time.

v2 changes:
  * relative residual for root acceptance (eq_70 is divided by k^2 lambda_D^2, so
    its terms are huge and nearly cancel -- an absolute 1e-3 test is far too strict)
  * root tracking by continuation with linear extrapolation; fluid seed only at the start
  * finer k grid over the crossover region, fixed y-limits so the fluid curve is visible
  * optional printout of the tracked root and its residual
"""
import os
import numpy as np
import matplotlib.pyplot as plt

import constants as C
from root import Z_tilde, custom_newton_raphson_2d

RUN = 0
REL_TOL = 1e-8          # acceptance threshold on the relative residual


# ------------------------------------------------------------------ shared parameters
def params(run=RUN):
    F = C.FLAGS
    kE, kI = C.kappa_e[run], C.kappa_i[run]
    psi_e = 1 + kE**2

    V_de_x = (C.E_0 / C.B_0) * kE**2 / psi_e
    V_de_y = -(C.E_0 / C.B_0) * kE / psi_e
    V_di_x = (kI**2 * C.E_0) / (C.B_0 * (1 + kI**2))
    V_di_y = (kI * C.E_0) / (C.B_0 * (1 + kI**2))

    T_e = C.T_n + F["HEAT_ELECTRONS"] * (C.m_n_e / 3.0) * (V_de_x**2 + V_de_y**2)
    T_eff_i = C.T_n + F["HEAT_IONS"] * (C.m_n_i / 3.0) * (V_di_x**2 + V_di_y**2)

    ratio = C.anisotropy_ratio if F["ION_ANISOTROPY"] else 1.0
    T_y = 3.0 * T_eff_i / (1.0 + 2.0 * ratio)
    T_x = ratio * T_y

    if F["CONSISTENT_NORMALIZATION"]:
        lam_sq = C.eps0 * T_e / (C.e**2 * C.n0)
        beta_T = C.T_n / T_e
    else:
        lam_sq, beta_T = C.lambda_D_sq, C.beta_T

    return dict(
        psi_e=psi_e, V_de_x=V_de_x, V_de_y=V_de_y, T_e=T_e,
        V_Tx=np.sqrt(T_x / C.m_i), V_Ty=np.sqrt(T_y / C.m_i),
        U0_y=F["ION_DRIFT"] * C.e * C.E_0 / (C.m_i * C.nu_i[run]),
        a_i0_y=F["ION_FIELD_TERMS"] * C.e * C.E_0 / C.m_i,
        lam_sq=lam_sq, beta_T=beta_T,
        nu_i=C.nu_i[run], nu_e=C.nu_e[run],
    )


# ------------------------------------------------------------------ kinetic equation
def eq_70_test(omega, kx, ky, run=RUN, terms=False):
    p = params(run)
    k_sq = max(kx**2 + ky**2, 1e-12)
    growth = np.imag(omega)

    # electron (fluid) term
    k_dot_V = kx * p["V_de_x"] + ky * p["V_de_y"]
    D_e = (-1j * omega * p["nu_e"] * p["psi_e"]
           + (C.gamma * p["T_e"] / C.m_e) * k_sq
           + 1j * p["nu_e"] * p["psi_e"] * k_dot_V)
    I_el = C.FLAGS["ELECTRON_SIGN"] * C.omega_pe_squared / D_e

    # ion (kinetic) term
    k_dot_a = ky * p["a_i0_y"]
    k_dot_U = ky * p["U0_y"]
    omega_i = omega + 1j * p["nu_i"]

    den1 = k_sq * C.V_Tn**2 + 1j * k_dot_a
    s1 = np.sqrt(2 * den1)
    z1 = omega_i / s1

    den2 = kx**2 * p["V_Tx"]**2 + ky**2 * p["V_Ty"]**2 + 1j * k_dot_a
    s2 = np.sqrt(2 * den2)
    z2 = (omega_i - k_dot_U) / s2

    pre = C.V_Tn**2 / (p["lam_sq"] * p["beta_T"] * den1)
    part2 = 1.0 / (1.0 + (1j * p["nu_i"] / s1) * Z_tilde(z1, growth, p["nu_i"]))
    part3 = 1.0 + ((omega_i - k_dot_U) / s2) * Z_tilde(z2, growth, p["nu_i"])
    I_ion = pre * part2 * part3

    if terms:
        return I_el, I_ion
    # ion term must be subtracted to cancel against I_el in the k->0 (fluid) limit
    return 1.0 + I_el - I_ion


def rel_res(w, kx, ky, run=RUN):
    """Residual relative to the size of the terms that are cancelling."""
    el, ion = eq_70_test(w, kx, ky, run, terms=True)
    return abs(1.0 + el - ion) / (1.0 + abs(el) + abs(ion))


# ------------------------------------------------------------------ fluid dispersion (same parameters)
def fluid_roots(kx, ky, run=RUN):
    """
    Quasi-neutral fluid limit of the same model:
      w^2 + i w [nu_i + (m_e/m_i) nu_e psi]
          - k^2 (V_Tn^2 + gamma T_e/m_i) - i (m_e/m_i) nu_e psi k.V_de = 0
    (no ion heating, no a_i0, no U0 -- compare with those switches off)
    """
    p = params(run)
    k_sq = kx**2 + ky**2
    mr = C.m_e / C.m_i
    b = 1j * (p["nu_i"] + mr * p["nu_e"] * p["psi_e"])
    k_dot_V = kx * p["V_de_x"] + ky * p["V_de_y"]
    c = -k_sq * (C.V_Tn**2 + C.gamma * p["T_e"] / C.m_i) - 1j * mr * p["nu_e"] * p["psi_e"] * k_dot_V
    d = np.sqrt(b**2 - 4 * c)
    return ((-b + d) / 2, (-b - d) / 2)


# ------------------------------------------------------------------ root tracking along a ray
def scan_ray(theta_deg, ks, run=RUN, verbose=False):
    """
    Follow ONE kinetic root by continuation.
    Seeds: linear extrapolation of the last two roots, then the last root;
    the fluid root is used only for the first few points.
    A point is accepted if the RELATIVE residual < REL_TOL; otherwise it stays NaN
    and the tracker keeps its previous history (so it can try again at the next k).
    """
    th = np.deg2rad(theta_deg)
    fl = np.full(len(ks), np.nan, dtype=complex)
    kin = np.full(len(ks), np.nan, dtype=complex)
    f = lambda w, a, b: eq_70_test(w, a, b, run)

    hist = []
    for n, k in enumerate(ks):
        kx, ky = k * np.cos(th), k * np.sin(th)
        r1, r2 = fluid_roots(kx, ky, run)
        fr = r1 if r1.imag > r2.imag else r2
        fl[n] = fr

        seeds = []
        if len(hist) >= 2:
            seeds.append(2 * hist[-1] - hist[-2])
        if hist:
            seeds.append(hist[-1])
        if len(hist) < 3:
            seeds.append(fr)

        for s in seeds:
            root = custom_newton_raphson_2d(f, s, kx, ky, max_iter=C.NR_MAX_ITER,
                                            max_backtrack=C.NR_MAX_BACKTRACK)
            if np.isfinite(root) and rel_res(root, kx, ky, run) < REL_TOL:
                kin[n] = root
                hist.append(root)
                break

        if verbose:
            if np.isfinite(kin[n]):
                print(f"k={k:6.3f}  kinetic w={kin[n]:.2f}  fluid w={fl[n]:.2f}  "
                      f"rel_res={rel_res(kin[n], kx, ky, run):.1e}")
            else:
                print(f"k={k:6.3f}  NO ROOT ACCEPTED  (fluid w={fl[n]:.2f})")
    return fl, kin


def label_from_flags():
    F = C.FLAGS
    on = [k for k, v in F.items() if k not in ("ELECTRON_SIGN", "CONSISTENT_NORMALIZATION") and v]
    s = "+".join(on) if on else "all_off"
    return f"{s}  (el sign {F['ELECTRON_SIGN']:+d})"


def plot_scan(theta_deg=60, kmin=0.2, kmax=6.0, npts=150, stages=None,
              fname="ray_scan.png", ylim=(-3000, 3000), verbose=False):
    ks = np.linspace(kmin, kmax, npts)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    stages = stages or [("current FLAGS", {})]
    default_flags = dict(C.FLAGS)
    fluid_drawn = False

    for name, over in stages:
        C.FLAGS.update(default_flags)
        C.FLAGS.update(over)
        fl, kin = scan_ray(theta_deg, ks, verbose=verbose)
        if not fluid_drawn:
            ax[0].plot(ks, fl.imag, "k--", lw=2, label="fluid (this stage's T_e)")
            ax[1].plot(ks, fl.real, "k--", lw=2, label="fluid")
            fluid_drawn = True
        ax[0].plot(ks, kin.imag, label=f"kinetic: {name}")
        ax[1].plot(ks, kin.real, label=f"kinetic: {name}")
        n_fail = int(np.sum(~np.isfinite(kin)))
        print(f"[{name}] flags: {label_from_flags()}   points with no accepted root: {n_fail}/{npts}")
    C.FLAGS.update(default_flags)

    kc = C.nu_i[RUN] / C.V_Tn
    for a, yl in zip(ax, ("growth rate gamma [1/s]", "real frequency omega_r [1/s]")):
        a.axvline(kc, ls=":", color="gray", label=r"$kV_{Tn}=\nu_i$")
        a.axhline(0, color="gray", lw=0.5)
        a.set_xlabel("|k| [1/m]"); a.set_ylabel(yl)
        a.set_title(f"ray at {theta_deg} deg from kx axis")
    if ylim is not None:
        ax[0].set_ylim(*ylim)
    ax[0].legend(fontsize=7)
    plt.tight_layout()
    os.makedirs("plots", exist_ok=True)
    plt.savefig(os.path.join("plots", fname), dpi=200)
    plt.show()


# ------------------------------------------------------------------ sign diagnostic
def cancellation_check(k=0.2, theta_deg=60):
    """
    At a fluid root, I_el + I_ion should nearly cancel (quasi-neutrality).
    Prints |I_el + I_ion| / |I_el| for both electron signs.
    Run with all FLAGS off for a clean test.
    """
    th = np.deg2rad(theta_deg)
    kx, ky = k * np.cos(th), k * np.sin(th)
    saved = C.FLAGS["ELECTRON_SIGN"]
    for sgn in (+1, -1):
        C.FLAGS["ELECTRON_SIGN"] = sgn
        for w in fluid_roots(kx, ky):
            el, ion = eq_70_test(w, kx, ky, terms=True)
            print(f"sign {sgn:+d}  w={w:.1f}   |el-ion|/|el| = {abs(el - ion) / abs(el):.3e}")
    C.FLAGS["ELECTRON_SIGN"] = saved


# ------------------------------------------------------------------ what runs
# Each stage adds one piece of physics on top of the previous ones.
STAGES = [
    ("0 all off",            dict(HEAT_ELECTRONS=0, HEAT_IONS=0, ION_ANISOTROPY=0, ION_FIELD_TERMS=0, ION_DRIFT=0)),
    ("1 + U0",               dict(ION_DRIFT=1)),
    ("2 + ion heating",      dict(ION_DRIFT=1, HEAT_IONS=1)),
    ("3 + anisotropy",       dict(ION_DRIFT=1, HEAT_IONS=1, ION_ANISOTROPY=1)),
    ("4 + a_i0 terms",       dict(ION_DRIFT=1, HEAT_IONS=1, ION_ANISOTROPY=1, ION_FIELD_TERMS=1)),
    ("5 + electron heating", dict(ION_DRIFT=1, HEAT_IONS=1, ION_ANISOTROPY=1, ION_FIELD_TERMS=1, HEAT_ELECTRONS=1)),
]

if __name__ == "__main__":
    cancellation_check()                       # sign test, uses current FLAGS (default: all off)

    # Stage 0: kinetic must match fluid at small k, then roll over smoothly.
    # verbose=True prints the tracked root and residual at every k (set False once happy).
    plot_scan(theta_deg=-45, stages=[STAGES[0]], kmax=8, npts=150,
          fname="stage0_m45.png", verbose=True)

    # Once stage 0 looks right, uncomment to see where each added term breaks it:
    # plot_scan(stages=STAGES, kmax=20, npts=120, ylim=None, fname="all_stages.png")
