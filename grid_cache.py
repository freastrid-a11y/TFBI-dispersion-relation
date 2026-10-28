"""
Shared grid-solving and CSV caching layer.

Each (equation, run, kx grid, ky grid) setup is solved once and cached as a
CSV under DATA_DIR. Re-running the same setup loads the CSV instead of
re-solving the dispersion relation. Methods in methods/methods.py consume
this cached data rather than each recomputing their own grid.

A cached CSV is automatically treated as stale (and recomputed) if
constants.py or the dispersion equation's source file has been modified
since the CSV was written, so editing constants.py doesn't silently leave
old results in place.
"""

import inspect
import os

import numpy as np
import pandas as pd

from root import solve_grid

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
CONSTANTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "constants.py")


def _grid_spec(kx_vals, ky_vals):
    return (
        f"kx{kx_vals[0]:g}_{kx_vals[-1]:g}_n{len(kx_vals)}"
        f"_ky{ky_vals[0]:g}_{ky_vals[-1]:g}_n{len(ky_vals)}"
    )


def csv_path(dispersion, kx_vals, ky_vals, run, data_dir=DATA_DIR):
    name = f"{dispersion.__name__}_run{run}_{_grid_spec(kx_vals, ky_vals)}.csv"
    return os.path.join(data_dir, name)


def _sources_mtime(dispersion):
    """Latest modification time among constants.py and the equation's source file."""
    paths = [CONSTANTS_PATH]

    try:
        paths.append(inspect.getfile(dispersion))
    except TypeError:
        pass

    return max(os.path.getmtime(p) for p in paths if os.path.exists(p))


def get_grid(dispersion, kx_vals, ky_vals, run, data_dir=DATA_DIR, force=False, csv_name=None):
    """
    Return (kx_vals, ky_vals, growth, frequency) for this exact setup.

    Loads the cached CSV if it already exists and is not stale; otherwise
    solves the grid and writes it out. growth and frequency are 2D arrays
    with shape (len(ky_vals), len(kx_vals)) holding Im(omega) and Re(omega).
    force=True always recomputes, regardless of staleness.
    csv_name overrides the auto-generated filename (still placed in data_dir).
    """
    kx_vals = np.asarray(kx_vals, dtype=float)
    ky_vals = np.asarray(ky_vals, dtype=float)
    shape = (len(ky_vals), len(kx_vals))

    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, csv_name) if csv_name else csv_path(dispersion, kx_vals, ky_vals, run, data_dir)
    exists = os.path.exists(path)

    stale = exists and _sources_mtime(dispersion) > os.path.getmtime(path)

    if stale:
        print(f"constants.py or {dispersion.__name__}'s source changed since {path} was written; recomputing.")

    if not force and not stale and exists:
        print(f"Loading cached grid: {path}")
        df = pd.read_csv(path)
        growth = df["growth_rate_s^-1"].to_numpy().reshape(shape)
        frequency = df["frequency_rad_s^-1"].to_numpy().reshape(shape)
        return kx_vals, ky_vals, growth, frequency

    omega_grid = solve_grid(dispersion, kx_vals, ky_vals, run)

    kx_grid, ky_grid = np.meshgrid(kx_vals, ky_vals)
    df = pd.DataFrame({
        "kx_m^-1": kx_grid.ravel(),
        "ky_m^-1": ky_grid.ravel(),
        "growth_rate_s^-1": omega_grid.imag.ravel(),
        "frequency_rad_s^-1": omega_grid.real.ravel(),
    })
    df.to_csv(path, index=False)
    print(f"Saved grid: {path}")

    return kx_vals, ky_vals, omega_grid.imag, omega_grid.real
