"""
Four-panel figure: "old terms" (top row) vs "new terms" (bottom row),
at E_0 (left) and 2 E_0 (right), for one altitude.

    old terms: fluid electrons, ions without i k.a_i0 (dispersion_rel_reduced)
    new terms: kinetic ions + kinetic electrons, Eq. 72 (eq_72_kinetic)

Plots the growth rate Im(omega), or the frequency Re(omega) with
--field frequency. Only unstable points (growth > 0) are shown.

    python fourpanel.py --alt 110
    python fourpanel.py --alt 104 --field frequency   # Re(omega) instead
    python fourpanel.py --alt 110 --separate-scales   # one colour bar per panel

Reads the previously saved grids in data/:
    old, E   : data/<alt>_reduced.csv
    old, 2E  : data/<alt>_reduced_doubleE.csv
    new, E   : data/<alt>kinetic.csv
    new, 2E  : data/<alt>_doubleEkinetic.csv
By default all four panels share one colour scale, so they can be compared directly.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PANELS = [  # (row, col, csv suffix, label, field multiplier)
    (0, 0, "_reduced", "Old terms", 1),
    (0, 1, "_reduced_doubleE", "Old terms", 2),
    # Must match the CSV names previously saved in data/ ("<alt>kinetic.csv"
    # and "<alt>_doubleEkinetic.csv"; the older "<alt>_kinetic*.csv" files
    # are stale copies).
    (1, 0, "kinetic", "New terms", 1),
    (1, 1, "_doubleEkinetic", "New terms", 2),
]


def panel_title(terms, e_mult):
    """Panel title with the field in mV/m, e.g. 'New terms, E = 105 mV/m'."""
    from constants import E_0
    return f"{terms}, $E$ = {round(e_mult * E_0 * 1e3)} mV/m"


def load(path):
    """Read a grid CSV; return kx, ky, growth and frequency, NaN where stable."""
    df = pd.read_csv(path).sort_values(["ky_m^-1", "kx_m^-1"])
    kx = np.sort(df["kx_m^-1"].unique()); ky = np.sort(df["ky_m^-1"].unique())
    g = df["growth_rate_s^-1"].to_numpy().reshape(len(ky), len(kx))
    f = df["frequency_rad_s^-1"].to_numpy().reshape(len(ky), len(kx))
    unstable = g > 0
    return kx, ky, np.where(unstable, g, np.nan), np.where(unstable, f, np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt", type=int, default=110, choices=[104, 110])
    ap.add_argument("--field", choices=["growth", "frequency"], default="growth",
                    help="growth: Im(omega); frequency: Re(omega) of the dominant unstable mode")
    ap.add_argument("--separate-scales", action="store_true")
    ap.add_argument("--save", default=None,
                    help="default: plots/four_panel_<alt>_growth.png (growth) or plots/four_panel_<alt>_freq.png (frequency)")
    a = ap.parse_args()

    data = {s: load(os.path.join(HERE, "data", f"{a.alt}{s}.csv")) for _, _, s, _, _ in PANELS}

    freq = a.field == "frequency"
    idx = 3 if freq else 2
    cmap = "RdBu_r" if freq else "viridis"
    label = (r"Real frequency $\mathrm{Re}(\omega)$ [rad s$^{-1}$]" if freq
             else r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")
    # Re(omega) can be negative, so use a range symmetric about 0 (keeps 0
    # in the middle of the red-blue colour map).
    panel_max = {s: np.nanmax(np.abs(d[idx])) if freq else np.nanmax(d[idx])
                 for s, d in data.items()}
    vmax = max(panel_max.values())

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.5), sharex=True, sharey=True, constrained_layout=True)
    mesh = None
    for r, c, s, terms, e_mult in PANELS:
        kx, ky = data[s][0], data[s][1]
        z = data[s][idx]
        ax = axes[r, c]
        hi = panel_max[s] if a.separate_scales else vmax
        lo = -hi if freq else 0
        mesh = ax.pcolormesh(kx, ky, np.ma.masked_invalid(z), cmap=cmap, vmin=lo, vmax=hi, shading="auto")
        if a.separate_scales:
            fig.colorbar(mesh, ax=ax, label=label)
        ax.set_title(panel_title(terms, e_mult), fontsize=10, pad=4)
        ax.set_aspect("equal")
        if r == 1:
            ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
        if c == 0:
            ax.set_ylabel(r"$k_y$ [m$^{-1}$]")
    if not a.separate_scales:
        fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.02, label=label)
    fig.suptitle(f"{'Real frequency' if freq else 'Growth rate'}, effective altitude {a.alt} km")

    save = a.save or os.path.join(HERE, "plots", f"four_panel_{a.alt}{'_freq' if freq else '_growth'}.png")
    fig.savefig(save, dpi=200)
    print("saved", save)


if __name__ == "__main__":
    main()
