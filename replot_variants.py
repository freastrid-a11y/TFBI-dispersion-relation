"""
Regenerate the plot for every variant:

  normal           dispersion_rel                         -> plots/<alt>.png
  ai0zero          dispersion_rel with a_i0 = 0           -> plots/<alt>_ai0zero.png
  doubleE          dispersion_rel with E_0 doubled        -> plots/<alt>_double_E.png
  doubleE_ai0zero  E_0 doubled and a_i0 = 0               -> plots/<alt>_double_E_ai0zero.png

"Old terms" (dispersion_rel_reduced: fluid electrons, ions without i k.a_i0).
This equation has no a_i0 in it, so "reduced_ai0zero" gives exactly the same
result as "reduced"; it only exists so the filenames match the other sets:

  reduced                  -> plots/<alt>_reduced.png
  reduced_ai0zero          -> plots/<alt>_reduced_ai0zero.png
  reduced_doubleE          -> plots/<alt>_reduced_double_E.png

"New terms", kinetic electrons (eq_72_kinetic / eq_72_kinetic_ai0zero):

  normalkinetic            -> plots/<alt>kinetic.png
  ai0zerokinetic           -> plots/<alt>_ai0zerokinetic.png
  doubleEkinetic           -> plots/<alt>_double_Ekinetic.png
  doubleE_ai0zerokinetic   -> plots/<alt>_double_E_ai0zerokinetic.png

Each variant runs as a separate Python process. That way the E_SCALE
setting (2 for doubled field) is in place before constants.py is loaded,
so T_e, lambda_D, the drifts, a_i0, ... are all computed with the new field,
just as if V_d had been changed in constants.py.

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

# variant name -> (equation in main.EQUATIONS, E-field scale, plot suffix)
VARIANTS = {
    "normalkinetic": ("eq_72_kinetic", 1.0, "kinetic"),
    "ai0zerokinetic": ("eq_72_kinetic_ai0zero", 1.0, "_ai0zerokinetic"),
    "doubleEkinetic": ("eq_72_kinetic", 2.0, "_double_Ekinetic"),
    "doubleE_ai0zerokinetic": ("eq_72_kinetic_ai0zero", 2.0, "_double_E_ai0zerokinetic"),
    "normal": ("dispersion_rel", 1.0, ""),
    "ai0zero": ("dispersion_rel_ai0zero", 1.0, "_ai0zero"),
    "doubleE": ("dispersion_rel", 2.0, "_double_E"),
    "doubleE_ai0zero": ("dispersion_rel_ai0zero", 2.0, "_double_E_ai0zero"),
    # dispersion_rel_reduced has no a_i0, so reduced_ai0zero is the same
    # equation as reduced (kept only so the filenames match).
    "reduced": ("dispersion_rel_reduced", 1.0, "_reduced"),
    "reduced_ai0zero": ("dispersion_rel_reduced", 1.0, "_reduced_ai0zero"),
    "reduced_doubleE": ("dispersion_rel_reduced", 2.0, "_reduced_double_E"),
}

# alt_eff index -> altitude (km), used in plot filenames
RUNS = {4: 110, 5: 104}


def replot_all_variants(method="growth_contour", runs=None, variants=None,
                        points=100, k_max=100.0):
    """Regenerate the plot for every (run, variant) combination.

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
    """Worker: solve one variant's grid and plot it.

    Must run in a fresh process: E_SCALE has to be set before constants.py
    (imported via main) computes its module-level derived values.
    """
    equation, e_scale, plot_suffix = VARIANTS[variant]
    os.environ["E_SCALE"] = repr(e_scale)

    import main  # imported only now, after E_SCALE is set

    # Hide stable points (growth <= 0) if the plot method supports that
    # (e.g. growth_contour).
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
