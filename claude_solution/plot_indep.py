"""
Plot the independent Eq. 72 solution produced by solve.py.

    python plot_indep.py --alt 104                 # growth-rate map only
    python plot_indep.py --alt 104 --compare       # + repo CSV + relative difference

--compare looks for data/<alt>_kinetic.csv in this folder or any parent folder
(i.e. the repo's data/ directory); override with --csv.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))


def find_csv(alt):
    d = HERE
    while True:
        cand = os.path.join(d, "data", f"{alt}_kinetic.csv")
        if os.path.exists(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt", type=int, required=True, choices=[104, 110])
    ap.add_argument("--npz", default=None, help="default: out/indep_<alt>.npz")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--save", default=None, help="default: out/indep_<alt>[_compare].png")
    a = ap.parse_args()

    npz = a.npz or os.path.join(HERE, "out", f"indep_{a.alt}.npz")
    d = np.load(npz)
    om, k = d["omega"], d["k"]
    n = len(k)
    mine = np.where(np.isfinite(om) & (om.imag > 0), om.imag, np.nan)
    vmax = np.nanmax(mine)

    def growth_panel(ax, data, title):
        m = ax.pcolormesh(k, k, np.ma.masked_invalid(data), cmap="viridis",
                          vmin=0, vmax=vmax, shading="auto")
        plt.colorbar(m, ax=ax, label=r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")
        ax.set_title(title)

    if not a.compare:
        fig, ax = plt.subplots(figsize=(8, 6))
        growth_panel(ax, mine, f"Independent Eq. 72 growth rate — altitude {a.alt} km")
        axes = [ax]
    else:
        csv = a.csv or find_csv(a.alt)
        if csv is None:
            raise SystemExit(f"Could not find data/{a.alt}_kinetic.csv; pass --csv PATH")
        import pandas as pd
        df = pd.read_csv(csv)
        g = df["growth_rate_s^-1"].to_numpy().reshape(n, n)
        f = df["frequency_rad_s^-1"].to_numpy().reshape(n, n)
        repo = f + 1j * g
        theirs = np.where(g > 1.0, g, np.nan)   # same 1 s^-1 floor as solve.py's search box

        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
        growth_panel(axes[0], mine, f"{a.alt} km — independent solve")
        growth_panel(axes[1], theirs, f"{a.alt} km — {os.path.basename(csv)}")

        both = np.isfinite(mine) & np.isfinite(theirs)
        rel = np.where(both, np.abs(om - repo) / np.abs(om), np.nan)
        m = axes[2].pcolormesh(k, k, np.ma.masked_invalid(rel), cmap="Greys",
                               norm=LogNorm(1e-16, 1e-11), shading="auto")
        plt.colorbar(m, ax=axes[2], label=r"$|\omega_{\rm indep}-\omega_{\rm repo}|/|\omega|$")
        for mask, colour, lab in ((np.isfinite(mine) & ~np.isfinite(theirs), "#d62728", "unstable only in independent"),
                                  (~np.isfinite(mine) & np.isfinite(theirs), "#1f77b4", "unstable only in repo")):
            iy, ix = np.nonzero(mask)
            if len(ix):
                axes[2].scatter(k[ix], k[iy], s=40, facecolors="none", edgecolors=colour,
                                lw=1.5, label=f"{lab} ({mask.sum()} px)")
        if axes[2].get_legend_handles_labels()[0]:
            axes[2].legend(loc="upper right", fontsize=9)
        axes[2].set_title(f"{a.alt} km — relative difference")
        print(f"common unstable points: {both.sum()}, max rel diff {np.nanmax(rel):.2e}, "
              f"only independent: {(np.isfinite(mine) & ~np.isfinite(theirs)).sum()}, "
              f"only repo: {(~np.isfinite(mine) & np.isfinite(theirs)).sum()}")

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel(r"$k_x$ [m$^{-1}$]")
        ax.set_ylabel(r"$k_y$ [m$^{-1}$]")

    save = a.save or os.path.join(HERE, "out", f"indep_{a.alt}{'_compare' if a.compare else ''}.png")
    os.makedirs(os.path.dirname(save), exist_ok=True)
    fig.savefig(save, dpi=200, bbox_inches="tight")
    print("saved", save)


if __name__ == "__main__":
    main()