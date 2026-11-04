"""
Main entry point: run a chosen method against a chosen dispersion equation,
using a cached CSV grid when one already exists for the given setup.

Example:
    python main.py --method growth_damping --save growth_damping_104km.png --run 5
"""

import argparse
import inspect

import numpy as np

from constants import alt_eff
from eq73 import eq_73, eq_73_ai0zero
from grid_cache import get_grid
from methods import METHODS

EQUATIONS = {
    "eq_73": eq_73,
    "eq_73_ai0zero": eq_73_ai0zero,
}


def run(equation, method, save_path, run=5, kx_min=-100.0, kx_max=100.0,
        ky_min=-100.0, ky_max=100.0, points=100, force=False, csv_name=None,
        positive_only=False):
    dispersion = EQUATIONS[equation] if isinstance(equation, str) else equation
    method_fn = METHODS[method]

    kwargs = {}
    if positive_only:
        if "positive_only" not in inspect.signature(method_fn).parameters:
            raise ValueError(
                f"Method '{method}' does not support --positive-only"
            )
        kwargs["positive_only"] = True

    kx_vals = np.linspace(kx_min, kx_max, points)
    ky_vals = np.linspace(ky_min, ky_max, points)

    kx_vals, ky_vals, growth, frequency = get_grid(
        dispersion, kx_vals, ky_vals, run, force=force, csv_name=csv_name
    )

    method_fn(kx_vals, ky_vals, growth, frequency, alt_eff[run], save_path, **kwargs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--equation", choices=sorted(EQUATIONS), default="eq_73")
    parser.add_argument("--method", choices=sorted(METHODS), required=True)
    parser.add_argument("--save", required=True, help="Output plot filename")
    parser.add_argument("--run", type=int, default=5, help="Index into alt_eff")
    parser.add_argument("--kx-min", type=float, default=-100.0)
    parser.add_argument("--kx-max", type=float, default=100.0)
    parser.add_argument("--ky-min", type=float, default=-100.0)
    parser.add_argument("--ky-max", type=float, default=100.0)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument(
        "--positive-only", action="store_true",
        help="Mask out non-positive growth rates (supported by growth_contour)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Recompute the grid even if a cached CSV already exists",
    )
    parser.add_argument(
        "--csv-name", default=None,
        help="Name for the cached grid CSV in data/ (default: auto-generated from equation/run/grid)",
    )
    args = parser.parse_args()

    run(
        args.equation,
        args.method,
        args.save,
        run=args.run,
        kx_min=args.kx_min,
        kx_max=args.kx_max,
        ky_min=args.ky_min,
        ky_max=args.ky_max,
        points=args.points,
        force=args.force,
        csv_name=args.csv_name,
        positive_only=args.positive_only,
    )


if __name__ == "__main__":
    main()
