"""
Compare fluid electrons (dispersion_rel.dispersion_rel) with fully kinetic electrons
(electron_kinetic.eq_72_kinetic, manuscript Eq. 75) on the same k-grid.

    python compare_electrons.py --run 5            # solve (or load cached) and plot
    python compare_electrons.py --run 5 --solve-only --equation eq_72_kinetic
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from constants import alt_eff
from grid_cache import get_grid
from main import EQUATIONS

HERE = os.path.dirname(os.path.abspath(__file__))


def grid(eq_name, run, points, kmax):
    kx = np.linspace(-kmax, kmax, points)
    ky = np.linspace(-kmax, kmax, points)
    alt = int(alt_eff[run])
    return get_grid(EQUATIONS[eq_name], kx, ky, run,
                    csv_name=f"{alt}_{eq_name}_n{points}.csv")


def plot(run, points, kmax, save):
    kx, ky, g_f, w_f = grid("dispersion_rel", run, points, kmax)
    _, _, g_k, w_k = grid("eq_72_kinetic", run, points, kmax)
    alt = alt_eff[run]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), constrained_layout=True)
    vmax = np.nanmax(np.abs(np.concatenate([g_f.ravel(), g_k.ravel()])))
    vmin = min(np.nanmin(g_f), np.nanmin(g_k))
    norm = TwoSlopeNorm(vmin=min(vmin, -1), vcenter=0, vmax=max(np.nanmax(g_f), 1))
    for ax, data, title in [(axes[0], g_f, "Fluid electrons (dispersion_rel)"),
                            (axes[1], g_k, "Kinetic electrons (Eq. 75)")]:
        m = ax.pcolormesh(kx, ky, np.ma.masked_invalid(data), shading="auto",
                          cmap="RdBu_r", norm=norm)
        if np.nanmin(data) < 0 < np.nanmax(data):
            ax.contour(kx, ky, data, levels=[0], colors="k", linewidths=1.3)
        ax.set_title(title)
        fig.colorbar(m, ax=ax, label=r"Im$(\omega)$ [s$^{-1}$]")
    diff = g_k - g_f
    lim = np.nanmax(np.abs(diff))
    m = axes[2].pcolormesh(kx, ky, np.ma.masked_invalid(diff), shading="auto",
                           cmap="PuOr_r", vmin=-lim, vmax=lim)
    axes[2].set_title("Kinetic − fluid growth rate")
    fig.colorbar(m, ax=axes[2], label=r"$\Delta$Im$(\omega)$ [s$^{-1}$]")
    for ax in axes:
        ax.set_xlabel(r"$k_x$ (along $V_d$) [m$^{-1}$]")
        ax.set_ylabel(r"$k_y$ (along $E_0$) [m$^{-1}$]")
        ax.set_aspect("equal")
    fig.suptitle(f"Growth rate, effective altitude {alt:g} km — black line: γ = 0")
    fig.savefig(save, dpi=200)
    plt.close(fig)

    def stats(g):
        i = np.nanargmax(g)
        return np.nanmax(g), kx[i % len(kx)], ky[i // len(kx)], np.mean(g > 0)
    for name, g in [("fluid", g_f), ("kinetic", g_k)]:
        gm, kxm, kym, frac = stats(g)
        print(f"{alt:g} km {name:8s}: max γ = {gm:9.1f} s^-1 at k=({kxm:.1f},{kym:.1f}); "
              f"unstable fraction of grid = {frac:.1%}; NaN = {np.mean(np.isnan(g)):.1%}")
    print(f"Saved {save}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=int, default=5)
    ap.add_argument("--points", type=int, default=100)
    ap.add_argument("--kmax", type=float, default=100.0)
    ap.add_argument("--equation", default=None)
    ap.add_argument("--solve-only", action="store_true")
    a = ap.parse_args()
    if a.solve_only:
        grid(a.equation, a.run, a.points, a.kmax)
    else:
        plot(a.run, a.points, a.kmax,
             os.path.join(HERE, "plots", f"electrons_fluid_vs_kinetic_{int(alt_eff[a.run])}km.png"))