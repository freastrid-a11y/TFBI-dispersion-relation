# TFBI Solver Documentation

Welcome to the documentation for the Thermal Farley-Buneman Instability (TFBI) solver.

## Getting Started
To run the solver, use the main entry point:

```bash
python main.py --method growth_damping --save output.png --run 5
```

## Module Guide

| Module | Role |
|---|---|
| [dispersion_rel](dispersion_rel.md) | The TFBI dispersion relation with fluid electrons (manuscript Eq. 73) — the core math. |
| [dispersion_rel_reduced](dispersion_rel_reduced.md) | Reduced variant: same fluid electrons, ion term with the $i\,\mathbf{k}\cdot\mathbf{a}_{i0}$ terms dropped. |
| [electron_kinetic](electron_kinetic.md) | Fully kinetic electron contribution and the full kinetic dispersion relation (Eqs. 67, 70, 72, 75). |
| [root](root.md) | Modified plasma dispersion function, multistart root finder, grid sweep. |
| [constants](constants.md) | Physical constants and altitude-dependent parameter arrays. |
| [main](main.md) | Command-line entry point tying everything together. |
| [methods](methods.md) | Plotting methods (growth, frequency, ratio, phase velocity). |
| [csvmethods](csvmethods.md) | Post-processing of CSV grids in `data/` (a_i0 differences). |
| [replot_variants](replot_variants.md) | Batch script regenerating all variant plots (normal, a_i0 = 0, doubled E-field, reduced equation, kinetic electrons). |
| [fourpanel](fourpanel.md) | Four-panel figure comparing old (reduced) vs. new (kinetic Eq. 72) growth rates at E_0 and 2E_0. |