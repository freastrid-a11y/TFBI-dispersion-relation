# `compare_electrons.py` — Fluid vs. Kinetic Electron Comparison

## Overview

Script that compares **fluid electrons** ([`dispersion_rel.dispersion_rel`](dispersion_rel.md)) with **fully kinetic electrons** ([`electron_kinetic.eq_72_kinetic`](electron_kinetic.md), manuscript Eq. 75) on the same $k$-grid. It produces a three-panel figure:

1. growth rate $\mathrm{Im}(\omega)$ with fluid electrons,
2. growth rate with kinetic electrons,
3. the pointwise difference (kinetic − fluid),

with the $\gamma = 0$ contour drawn in black on the first two panels, and prints summary statistics for each model (max $\gamma$ and its location, unstable fraction of the grid, NaN fraction).

Grids are solved (or loaded) through [`grid_cache.get_grid`](grid_cache.md) with the cache-name convention `data/<altitude>_<equation>_n<points>.csv`, so repeat runs are instant. The figure is saved as `plots/electrons_fluid_vs_kinetic_<altitude>km.png`.

---

## API Reference

### `grid(eq_name, run, points, kmax)`

Solve (or load from cache) one equation's grid.

| Argument | Type | Description |
|---|---|---|
| `eq_name` | `str` | Key in `main.EQUATIONS` (e.g. `dispersion_rel`, `eq_72_kinetic`). |
| `run` | `int` | Altitude index into `alt_eff`. |
| `points` | `int` | Grid points per axis. |
| `kmax` | `float` | Grid spans ±`kmax` in both $k_x$ and $k_y$ (m⁻¹). |

**Returns:** `(kx_vals, ky_vals, growth, frequency)` as from `get_grid`, cached under `data/<altitude>_<eq_name>_n<points>.csv`.

### `plot(run, points, kmax, save)`

Build both grids (fluid and kinetic), render the three-panel comparison figure, print per-model statistics, and save to `save`.

### CLI

| Flag | Default | Description |
|---|---|---|
| `--run N` | `5` | Index into `alt_eff` (4 → 110 km, 5 → 104 km). |
| `--points N` | `100` | Grid resolution per axis. |
| `--kmax X` | `100.0` | Wave-number bound (m⁻¹). |
| `--equation NAME` | — | Only used with `--solve-only`: which equation to solve. |
| `--solve-only` | off | Just solve/cache the grid for `--equation`; no plotting. |

---

## Usage Examples

```bash
# Solve (or load cached) both models at 104 km and plot the comparison
python compare_electrons.py --run 5

# Pre-compute just the kinetic grid (e.g. in a separate long-running job)
python compare_electrons.py --run 5 --solve-only --equation eq_72_kinetic
```
