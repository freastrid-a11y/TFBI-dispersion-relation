"""
Four-panel growth-rate comparison: old (reduced: new a_i0/a_e0 terms dropped)
vs new (full kinetic Eq. 72), for E_0 and 2 E_0.

    python fourpanel.py --alt 110
    python fourpanel.py --alt 104
    python fourpanel.py --alt 110 --separate-scales   # one colour bar per panel

Reads the cached grids in data/ (regenerate them with replot_variants.py first):
    old, E   : data/<alt>_reduced.csv
    old, 2E  : data/<alt>_reduced_doubleE.csv
    new, E   : data/<alt>_kinetic.csv
    new, 2E  : data/<alt>_kinetic_doubleE.csv
By default all four panels share one colour scale so they can be compared directly.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PANELS = [  # (row, col, csv suffix, title)
    (0, 0, "_reduced", r"Old terms, $E_0$"),
    (0, 1, "_reduced_doubleE", r"Old terms, $2E_0$"),
    (1, 0, "_kinetic", r"New terms (Eq. 72), $E_0$"),
    (1, 1, "_kinetic_doubleE", r"New terms (Eq. 72), $2E_0$"),
]


def load(path):
    df = pd.read_csv(path).sort_values(["ky_m^-1", "kx_m^-1"])
    kx = np.sort(df["kx_m^-1"].unique()); ky = np.sort(df["ky_m^-1"].unique())
    g = df["growth_rate_s^-1"].to_numpy().reshape(len(ky), len(kx))
    return kx, ky, np.where(g > 0, g, np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt", type=int, default=110, choices=[104, 110])
    ap.add_argument("--separate-scales", action="store_true")
    ap.add_argument("--save", default=None, help="default: plots/four_panel_<alt>.png")
    a = ap.parse_args()

    data = {s: load(os.path.join(HERE, "data", f"{a.alt}{s}.csv")) for _, _, s, _ in PANELS}
    vmax = max(np.nanmax(g) for _, _, g in data.values())

    fig, axes = plt.subplots(2, 2, figsize=(11, 9.5), sharex=True, sharey=True, constrained_layout=True)
    mesh = None
    for r, c, s, title in PANELS:
        kx, ky, g = data[s]
        ax = axes[r, c]
        hi = np.nanmax(g) if a.separate_scales else vmax
        mesh = ax.pcolormesh(kx, ky, np.ma.masked_invalid(g), cmap="viridis", vmin=0, vmax=hi, shading="auto")
        if a.separate_scales:
            fig.colorbar(mesh, ax=ax, label=r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")
        i = np.unravel_index(np.nanargmax(g), g.shape)
        ax.set_title(f"{title}\n" + r"max $\gamma$ = %.0f s$^{-1}$ at $k$ = (%.0f, %.0f) m$^{-1}$" % (g[i], kx[i[1]], ky[i[0]]),
                     fontsize=10)
        ax.set_aspect("equal")
        if r == 1:
            ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
        if c == 0:
            ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    if not a.separate_scales:
        fig.colorbar(mesh, ax=axes, shrink=0.8, label=r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")
    fig.suptitle(f"Growth rate, altitude {a.alt} km (unstable roots only)")

    save = a.save or os.path.join(HERE, "plots", f"four_panel_{a.alt}.png")
    fig.savefig(save, dpi=200)
    print("saved", save)


if __name__ == "__main__":
    main()
