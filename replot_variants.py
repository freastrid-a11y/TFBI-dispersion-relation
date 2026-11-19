"""
Regenerate the cached CSV grids and replot every variant:

  normal           dispersion_rel                         -> data/<alt>.csv,                 plots/<alt>.png
  ai0zero          dispersion_rel with a_i0 = 0           -> data/<alt>_ai0zero.csv,         plots/<alt>_ai0zero.png
  doubleE          dispersion_rel with E_0 doubled        -> data/<alt>_doubleE.csv,         plots/<alt>_double_E.png
  doubleE_ai0zero  E_0 doubled and a_i0 = 0      -> data/<alt>_doubleE_ai0zero.csv, plots/<alt>_double_E_ai0zero.png

Reduced-equation variants (dispersion_rel_reduced: fluid electrons, ion term
with the i k.a_i0 terms dropped; a_i0 no longer appears, so its "ai0zero"
variant is the same equation and is provided only for filename parity):

  reduced                  -> data/<alt>_reduced.csv,                 plots/<alt>_reduced.png
  reduced_ai0zero          -> data/<alt>_reduced_ai0zero.csv,         plots/<alt>_reduced_ai0zero.png
  reduced_doubleE          -> data/<alt>_reduced_doubleE.csv,         plots/<alt>_reduced_double_E.png

Each variant runs in its own subprocess so that doubling the electric field
(via the E_SCALE environment variable read by constants.py) consistently
rescales every derived constant (T_e, lambda_D, drifts, a_i0, ...), exactly
as if V_d had been edited in constants.py.

Example:
    python replot_variants.py                          # everything, both altitudes
    python replot_variants.py --variants doubleE       # just the doubled field
    python replot_variants.py --points 20              # coarse grid for a quick look
"""

import argparse
import inspect
import os
import subprocess
import sys

from methods import METHODS

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTS_DIR = os.path.join(HERE, "plots")

# variant name -> (equation in main.EQUATIONS, E-field scale, csv suffix, plot suffix)
VARIANTS = {
    "normal": ("dispersion_rel", 1.0, "", ""),
    "ai0zero": ("dispersion_rel_ai0zero", 1.0, "_ai0zero", "_ai0zero"),
    "doubleE": ("dispersion_rel", 2.0, "_doubleE", "_double_E"),
    "doubleE_ai0zero": ("dispersion_rel_ai0zero", 2.0, "_doubleE_ai0zero", "_double_E_ai0zero"),
    # dispersion_rel_reduced contains no a_i0 at all, so its ai0zero variant
    # uses the same equation (kept for filename parity with the full set).
    "reduced": ("dispersion_rel_reduced", 1.0, "_reduced", "_reduced"),
    "reduced_ai0zero": ("dispersion_rel_reduced", 1.0, "_reduced_ai0zero", "_reduced_ai0zero"),
    "reduced_doubleE": ("dispersion_rel_reduced", 2.0, "_reduced_doubleE", "_reduced_double_E"),
}

# alt_eff index -> altitude (km), used in CSV/plot filenames
RUNS = {4: 110, 5: 104}


def replot_all_variants(method="growth_contour", runs=None, variants=None,
                        points=100, k_max=100.0):
    """Regenerate the CSV grid and replot for every (run, variant) combination.

    Spawns one subprocess per combination so each one imports constants.py
    fresh with the variant's E_SCALE applied.
    """
    runs = list(RUNS) if runs is None else list(runs)
    variants = list(VARIANTS) if variants is None else list(variants)

    for run in runs:
        for name in variants:
            print(f"=== altitude {RUNS[run]} km (run {run}), variant {name} ===", flush=True)
            subprocess.run(
                [
                    sys.executable, os.path.abspath(__file__),
                    "--worker",
                    "--variants", name,
                    "--runs", str(run),
                    "--method", method,
                    "--points", str(points),
                    "--k-max", str(k_max),
                ],
                check=True,
                cwd=HERE,
            )


def _run_variant(variant, run, method, points, k_max):
    """Worker: regenerate one variant's CSV (force=True) and replot it.

    Must run in a fresh process: E_SCALE has to be set before constants.py
    (imported via main) computes its module-level derived values.
    """
    equation, e_scale, csv_suffix, plot_suffix = VARIANTS[variant]
    os.environ["E_SCALE"] = repr(e_scale)

    import main  # imported only now, after E_SCALE is set

    # Only plot positive growth rates (masking gamma <= 0) whenever the
    # chosen method supports it (e.g. growth_contour).
    positive_only = "positive_only" in inspect.signature(METHODS[method]).parameters
    if not positive_only:
        print(f"note: method '{method}' does not support positive-only masking; "
              "plotting all values.", flush=True)

    altitude = RUNS[run]
    os.makedirs(PLOTS_DIR, exist_ok=True)

    main.run(
        equation,
        method,
        os.path.join(PLOTS_DIR, f"{altitude}{plot_suffix}.png"),
        run=run,
        kx_min=-k_max,
        kx_max=k_max,
        ky_min=-k_max,
        ky_max=k_max,
        points=points,
        force=True,
        csv_name=f"{altitude}{csv_suffix}.csv",
        positive_only=positive_only,
    )


def main_cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", choices=sorted(METHODS), default="growth_contour")
    parser.add_argument("--runs", type=int, nargs="+", choices=sorted(RUNS),
                        default=sorted(RUNS), help="alt_eff indices (4=110 km, 5=104 km)")
    parser.add_argument("--variants", nargs="+", choices=sorted(VARIANTS),
                        default=sorted(VARIANTS))
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--k-max", type=float, default=100.0)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.worker:
        for run in args.runs:
            for variant in args.variants:
                _run_variant(variant, run, args.method, args.points, args.k_max)
    else:
        replot_all_variants(
            method=args.method,
            runs=args.runs,
            variants=args.variants,
            points=args.points,
            k_max=args.k_max,
        )


if __name__ == "__main__":
    main_cli()
