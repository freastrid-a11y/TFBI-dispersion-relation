# --- ISOLATED SMALL-K TEST ---
from root import custom_newton_raphson_2d, fluid_seed_pair_2d
from eq70 import eq_70
import numpy as np
run_idx = 4
k_vals = np.logspace(-6, -1, 6)

print("\n--- Testing True Small-k Limit ---")
for k in k_vals:
    # 1. Get the standard fluid roots to use as a seed
    seed1, seed2 = fluid_seed_pair_2d(k, 0.0, run=run_idx)
    
    # 2. Find the root using the newly normalized eq_70
    root = custom_newton_raphson_2d(
        func=lambda omega, kx, ky: eq_70(omega, kx, ky, run=run_idx),
        initial_guess=seed1, # Try the first fluid branch
        kx=k,
        ky=0.0,
        tol_residual=1e-2 # Relaxed tolerance to account for the massive 1/k^2 scaling
    )
    
    # 3. Print the growth rate
    growth_rate = np.imag(root)
    print(f"k = {k:.1e} | Root = {root.real:.3e} + {growth_rate:.3e}j | Growth Rate = {growth_rate:.3e}")