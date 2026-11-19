# `fourpanel.py` — Four-Panel Old vs. New Growth-Rate Comparison

## Overview

Standalone plotting script producing a 2×2 growth-rate figure for one altitude, comparing the **old (reduced) relation** — [`dispersion_rel_reduced.py`](dispersion_rel_reduced.md), new $a_{i0}$/$a_{e0}$ terms dropped — against the **new fully kinetic relation** (Eq. 72, [`electron_kinetic.py`](electron_kinetic.md)), each at $E_0$ and $2E_0$:

| Panel | CSV read from `data/` |
|---|---|
| Old terms, $E_0$ | `<alt>_reduced.csv` |
| Old terms, $2E_0$ | `<alt>_reduced_doubleE.csv` |
| New terms (Eq. 72), $E_0$ | `<alt>_kinetic.csv` |
| New terms (Eq. 72), $2E_0$ | `<alt>_kinetic_doubleE.csv` |

It is a pure **post-processing** script: it only reads cached CSV grids (regenerate them with [`replot_variants.py`](replot_variants.md) first) and plots the unstable region ($\gamma > 0$; everything else is masked). Each panel's title reports the maximum growth rate and the wave vector where it occurs. By default all four panels share one colour scale so they can be compared directly.

The script makes **one figure per invocation** for a single altitude (`--alt`, default `110`); run it once per altitude to get both.

---

## CLI

| Flag | Default | Description |
|---|---|---|
| `--alt {104, 110}` | `110` | Altitude (km) whose cached grids to plot. |
| `--separate-scales` | off | Give each panel its own colour bar instead of one shared scale. |
| `--save PATH` | `plots/four_panel_<alt>.png` | Output figure filename. |

---

## Usage Examples

```bash
python fourpanel.py --alt 110
python fourpanel.py --alt 104
python fourpanel.py --alt 110 --separate-scales   # one colour bar per panel
```
