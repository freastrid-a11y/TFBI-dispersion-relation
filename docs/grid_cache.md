# `grid_cache.py` — Grid Solving & CSV Caching

## Overview

Solving the dispersion relation over a 100×100 wave-vector grid takes minutes; plotting takes seconds. This module decouples the two by caching each solved grid as a CSV under the repository's `data/` directory:

- Each **(equation, run, kx grid, ky grid)** setup maps to a deterministic file name, e.g.
  `eq_73_run5_kx-100_100_n100_ky-100_100_n100.csv`.
- On a repeat request the CSV is loaded instead of re-solving.
- A cached CSV is automatically treated as **stale** (and recomputed) if `constants.py` or the dispersion equation's source file has been modified since the CSV was written — so editing physics never silently reuses old results.
- `force=True` recomputes unconditionally; `csv_name` overrides the auto-generated name (useful for hand-named files like `104.csv` consumed by [`csvmethods.py`](csvmethods.md)).

The CSV is long-format with columns `kx_m^-1`, `ky_m^-1`, `growth_rate_s^-1` ($\mathrm{Im}\,\omega$), `frequency_rad_s^-1` ($\mathrm{Re}\,\omega$), one row per grid point.

---

## API Reference

### `get_grid(dispersion, kx_vals, ky_vals, run, data_dir=DATA_DIR, force=False, csv_name=None)`

Return the solved (or cached) grid for this exact setup.

| Argument | Type | Description |
|---|---|---|
| `dispersion` | `callable` | Dispersion residual, e.g. `eq_73`. Its `__name__` and source file participate in cache keying/staleness. |
| `kx_vals`, `ky_vals` | `array_like` of `float` | 1-D wave-number grids (m⁻¹). |
| `run` | `int` | Altitude index into `constants.py` arrays. |
| `data_dir` | `str` | Cache directory. Default: `data/` beside the module. |
| `force` | `bool` | Always recompute, even if a fresh cache exists. Default `False`. |
| `csv_name` | `str` or `None` | Override the auto-generated filename (still placed inside `data_dir`). |

**Returns:** `(kx_vals, ky_vals, growth, frequency)` — the 1-D grids as `ndarray[float]`, plus two 2-D `ndarray[float]` of shape `(len(ky_vals), len(kx_vals))` holding $\mathrm{Im}(\omega)$ and $\mathrm{Re}(\omega)$ (`nan` where no root was accepted).

### `csv_path(dispersion, kx_vals, ky_vals, run, data_dir=DATA_DIR)`

Compute the auto-generated cache path for a setup without solving anything.

| Argument | Type | Description |
|---|---|---|
| `dispersion` | `callable` | Equation whose `__name__` keys the file. |
| `kx_vals`, `ky_vals` | `array_like` | Grids (only endpoints and lengths are encoded). |
| `run` | `int` | Altitude index. |
| `data_dir` | `str` | Cache directory. |

**Returns:** `str` — e.g. `.../data/eq_73_run5_kx-100_100_n100_ky-100_100_n100.csv`.

### Module constants

| Name | Description |
|---|---|
| `DATA_DIR` | Absolute path to the `data/` cache directory next to the module. |
| `CONSTANTS_PATH` | Absolute path to `constants.py`, watched for staleness. |

---

## Usage Examples

Typical use goes through `main.py`, which calls `get_grid` for you:

```bash
# First run: solves and caches. Second run: loads the CSV instantly.
python main.py --method growth_contour --save growth_104km.png --run 5

# Recompute even though a cache exists
python main.py --method growth_contour --save growth_104km.png --run 5 --force

# Cache under a friendly name for later use with csvmethods.py
python main.py --method growth_contour --save growth_104km.png --run 5 --csv-name 104.csv
```

Direct use from Python:

```python
import numpy as np
from eq73 import eq_73
from grid_cache import get_grid, csv_path

kx = np.linspace(-100, 100, 100)
ky = np.linspace(-100, 100, 100)

print("cache file:", csv_path(eq_73, kx, ky, run=5))

kx, ky, growth, frequency = get_grid(eq_73, kx, ky, run=5)
print("max growth rate:", np.nanmax(growth), "1/s")
```
