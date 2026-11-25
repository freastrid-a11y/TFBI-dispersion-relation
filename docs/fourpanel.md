# `fourpanel.py` — Four-Panel Old vs. New Growth-Rate Comparison

## Overview

Standalone plotting script producing a 2×2 growth-rate figure for one altitude, comparing the **old (reduced) relation** — [`dispersion_rel_reduced.py`](dispersion_rel_reduced.md), new $a_{i0}$/$a_{e0}$ terms dropped — against the **new fully kinetic relation** (Eq. 72, [`electron_kinetic.py`](electron_kinetic.md)), each at $E_0$ and $2E_0$. Panel titles show the actual field strength, computed from `E_0` in [`constants.py`](constants.md) and rounded to whole mV/m (e.g. "Old terms, $E$ = 105 mV/m"):

| Panel | CSV read from `data/` |
|---|---|
| Old terms, $E$ = 105 mV/m | `<alt>_reduced.csv` |
| Old terms, $E$ = 210 mV/m | `<alt>_reduced_doubleE.csv` |
| New terms, $E$ = 105 mV/m | `<alt>kinetic.csv` |
| New terms, $E$ = 210 mV/m | `<alt>_doubleEkinetic.csv` |

It is a pure **post-processing** script: it only reads previously saved CSV grids from `data/` and plots the unstable region ($\gamma > 0$; everything else is masked). It plots the **growth rate** $\mathrm{Im}(\omega)$ by default, or the **real frequency** $\mathrm{Re}(\omega)$ of the dominant unstable mode with `--field frequency` (signed, so it uses a diverging colour map symmetric about 0). By default all four panels share one colour scale so they can be compared directly.

The script makes **one figure per invocation** for a single altitude (`--alt`, default `110`); run it once per altitude to get both.

---

## CLI

| Flag | Default | Description |
|---|---|---|
| `--alt {104, 110}` | `110` | Altitude (km) whose saved grids to plot. |
| `--field {growth, frequency}` | `growth` | Plot $\mathrm{Im}(\omega)$ (growth) or $\mathrm{Re}(\omega)$ (frequency) of the dominant unstable mode. |
| `--separate-scales` | off | Give each panel its own colour bar instead of one shared scale. |
| `--save PATH` | `plots/four_panel_<alt>_growth.png`, or `plots/four_panel_<alt>_freq.png` for `--field frequency` | Output figure filename. |

---

## Usage Examples

```bash
python fourpanel.py --alt 110
python fourpanel.py --alt 104 --field frequency   # Re(omega) instead of growth
python fourpanel.py --alt 110 --separate-scales   # one colour bar per panel
```
