# `csvmethods.py` — Post-Processing CSV Grids

## Overview

Analysis utilities that operate **directly on previously saved long-format CSV grids** in the repository's `data/` directory, without re-solving anything. They answer one physics question:

**How much does the ion acceleration term matter?** `plot_field_difference` takes the pointwise difference of a field (default: growth rate) between two CSVs sampled on the same grid — e.g. an `dispersion_rel` run vs. an `dispersion_rel_ai0zero` run — producing a map of $\Delta\,\mathrm{Im}(\omega)$ attributable to $\mathbf{a}_{i0}$.

CSV filenames processed in `__main__` follow the convention `data/<altitude>.csv` / `data/<altitude>_ai0zero.csv`, with columns `kx_m^-1`, `ky_m^-1`, `growth_rate_s^-1`, `frequency_rad_s^-1`.

---

## API Reference

### `load_field_grid(path, column="growth_rate_s^-1")`

Generic pivot of one CSV column into a 2-D grid.

**Returns:** `(kx_vals, ky_vals, grid)`.

### `plot_field_difference(path_nonzero, path_zero, save_path, column="growth_rate_s^-1", percentile=95)`

Contour of `path_nonzero` − `path_zero` for the chosen column, with a symmetric color scale clipped at `percentile` of $|{\rm diff}|$ so a few resonance outliers don't wash out the structure.

| Argument | Type | Description |
|---|---|---|
| `path_nonzero` | `str` | CSV of the full run (e.g. $a_{i0} \ne 0$). |
| `path_zero` | `str` | CSV of the reference run (e.g. $a_{i0} = 0$). Must share the exact same $k$ grid. |
| `save_path` | `str` | Output figure path. |
| `column` | `str` | CSV column to difference. Default growth rate. |
| `percentile` | `float` | Color-scale clip percentile of $\lvert\text{diff}\rvert$. Default `95`. |

**Returns:** `None`. **Raises** `ValueError` if the grids differ or no finite differences exist.

---

## Usage Examples

Run the built-in difference plots on existing CSVs (`data/104*.csv`, `data/110*.csv`):

```bash
python csvmethods.py
```

From Python:

```python
from csvmethods import plot_field_difference

# Effect of the ion acceleration term on the growth rate
plot_field_difference("data/104.csv", "data/104_ai0zero.csv", "plots/difference_104.png")
```
