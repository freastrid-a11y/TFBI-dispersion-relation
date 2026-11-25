"""
Four-panel figure for the independent solution, laid out exactly like the
repo's fourpanel.py: old terms (reduced) vs new terms (Eq. 72), at E_0 and 2 E_0.

    python fourpanel_indep.py --alt 110
    python fourpanel_indep.py --alt 104 --field frequency
    python fourpanel_indep.py --alt 104 --compare      # relative difference vs the repo CSVs

Reads out/indep_<alt>{_reduced, _reduced_doubleE, , _doubleE}.npz, written by
    python solve.py --alt <alt> [--E-scale 2] [--model reduced] --out ...

Only unstable roots are plotted.  solve.py searches Im(omega) > 1 s^-1, so the
same 1 s^-1 floor is applied to the repo grids in --compare.

--compare draws |omega_indep - omega_repo| / |omega_indep| for each panel and
circles points that are unstable in only one of the two solutions.  It reads
the repo CSVs fourpanel.py uses:
    data/<alt>_reduced.csv, data/<alt>_reduced_doubleE.csv,
    data/<alt>kinetic.csv,  data/<alt>_doubleEkinetic.csv
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
E0_MV = 2100 * 5e-5 * 1e3          # V_d B_0 in mV/m, same plasma as indep.py

# (row, col, npz suffix, repo csv suffix, terms label, E multiplier)
PANELS = [
    (0, 0, "_reduced", "_reduced", "Old terms", 1),
    (0, 1, "_reduced_doubleE", "_reduced_doubleE", "Old terms", 2),
    (1, 0, "", "kinetic", "New terms", 1),
    (1, 1, "_doubleE", "_doubleEkinetic", "New terms", 2),
]
GMIN = 1.0


def load_indep(alt, suffix):
    d = np.load(os.path.join(HERE, "out", f"indep_{alt}{suffix}.npz"))
    om = d["omega"]
    return d["k"], np.where(np.isfinite(om) & (om.imag > GMIN), om, complex(np.nan, np.nan))


def load_repo(alt, suffix, n):
    df = pd.read_csv(os.path.join(REPO, "data", f"{alt}{suffix}.csv")).sort_values(["ky_m^-1", "kx_m^-1"])
    om = (df["frequency_rad_s^-1"].to_numpy() + 1j * df["growth_rate_s^-1"].to_numpy()).reshape(n, n)
    return np.where(np.isfinite(om) & (om.imag > GMIN), om, complex(np.nan, np.nan))


def title(terms, mult):
    return f"{terms}, $E$ = {round(mult * E0_MV)} mV/m"


def finish(fig, axes, suptitle, save):
    for r in range(2):
        axes[1, r].set_xlabel(r"$k_x$ [m$^{-1}$]")
        axes[r, 0].set_ylabel(r"$k_y$ [m$^{-1}$]")
    for ax in axes.ravel():
        ax.set_aspect("equal")
    fig.suptitle(suptitle)
    os.makedirs(os.path.dirname(save), exist_ok=True)
    fig.savefig(save, dpi=200)
    print("saved", save)


def plot_solution(alt, field, separate):
    freq = field == "frequency"
    data = {s: load_indep(alt, s) for _, _, s, _, _, _ in PANELS}
    val = {s: (om.real if freq else om.imag) for s, (_, om) in data.items()}
    pmax = {s: np.nanmax(np.abs(v)) for s, v in val.items()}
    vmax = max(pmax.values())
    cmap = plt.get_cmap("RdBu_r" if freq else "viridis").copy()
    cmap.set_bad("white")                     # stable / no root: blank, as in fourpanel.py
    label = (r"Real frequency $\mathrm{Re}(\omega)$ [rad s$^{-1}$]" if freq
             else r"Growth rate $\mathrm{Im}(\omega)$ [s$^{-1}$]")

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.5), sharex=True, sharey=True, constrained_layout=True)
    for r, c, s, _, terms, mult in PANELS:
        k = data[s][0]
        hi = pmax[s] if separate else vmax
        mesh = axes[r, c].pcolormesh(k, k, np.ma.masked_invalid(val[s]), cmap=cmap,
                                     vmin=-hi if freq else 0, vmax=hi, shading="auto")
        if separate:
            fig.colorbar(mesh, ax=axes[r, c], label=label)
        axes[r, c].set_title(title(terms, mult), fontsize=10, pad=4)
    if not separate:
        fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.02, label=label)
    save = os.path.join(HERE, "out", f"four_panel_indep_{alt}{'_freq' if freq else '_growth'}.png")
    finish(fig, axes, f"Independent solution: {'real frequency' if freq else 'growth rate'}, "
                      f"effective altitude {alt} km", save)


def plot_compare(alt):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9.5), sharex=True, sharey=True, constrained_layout=True)
    rels, stats = {}, []
    for r, c, s, cs, terms, mult in PANELS:
        k, mine = load_indep(alt, s)
        repo = load_repo(alt, cs, len(k))
        a, b = np.isfinite(mine), np.isfinite(repo)
        rel = np.where(a & b, np.abs(mine - repo) / np.abs(mine), np.nan)
        rels[(r, c)] = (k, rel, a & ~b, b & ~a)
        stats.append((title(terms, mult).replace("$", ""), cs, (a & b).sum(),
                      np.nanmedian(rel), np.nanmax(rel), (a & ~b).sum(), (b & ~a).sum()))
    finite = np.concatenate([v[1][np.isfinite(v[1]) & (v[1] > 0)] for v in rels.values()])
    norm = LogNorm(max(finite.min(), 1e-16), max(finite.max(), 1e-15))
    for (r, c), (k, rel, only_mine, only_repo) in rels.items():
        ax = axes[r, c]
        cm = plt.get_cmap("magma_r").copy(); cm.set_bad("white")
        mesh = ax.pcolormesh(k, k, np.ma.masked_invalid(rel), cmap=cm, norm=norm, shading="auto")
        for mask, colour, lab in ((only_mine, "#1f77b4", "unstable only in independent"),
                                  (only_repo, "#d62728", "unstable only in repo")):
            iy, ix = np.nonzero(mask)
            if len(ix):
                ax.scatter(k[ix], k[iy], s=6, c=colour, lw=0, label=f"{lab} ({mask.sum()})")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="lower left", fontsize=7, markerscale=2)
        st = next(x for x in stats if x[1] == PANELS[2 * r + c][3])
        ax.set_title(f"{PANELS[2 * r + c][4]}, $E$ = {round(PANELS[2 * r + c][5] * E0_MV)} mV/m\n"
                     f"median {st[3]:.1e}, max {st[4]:.1e}", fontsize=10, pad=4)
    fig.colorbar(mesh, ax=axes, shrink=0.8, pad=0.02,
                 label=r"$|\omega_{\rm indep}-\omega_{\rm repo}|\,/\,|\omega_{\rm indep}|$")
    finish(fig, axes, f"Independent vs repo four-panel grids, effective altitude {alt} km",
           os.path.join(HERE, "out", f"four_panel_compare_{alt}.png"))
    print(f"{'panel':28s} {'repo csv':18s} {'common':>6s} {'median rel':>10s} {'max rel':>9s} "
          f"{'only indep':>10s} {'only repo':>9s}")
    for t, cs, n, med, mx, oi, orp in stats:
        print(f"{t:28s} {cs:18s} {n:6d} {med:10.1e} {mx:9.1e} {oi:10d} {orp:9d}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt", type=int, default=110, choices=[104, 110])
    ap.add_argument("--field", choices=["growth", "frequency"], default="growth")
    ap.add_argument("--separate-scales", action="store_true")
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    if a.compare:
        plot_compare(a.alt)
    else:
        plot_solution(a.alt, a.field, a.separate_scales)


if __name__ == "__main__":
    main()
