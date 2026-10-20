import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("csv/110_growth_rates_a0_zero.csv")

# Extract columns correctly
kx_all = df.iloc[:, 0].values
ky_all = df.iloc[:, 1].values
growth_all = df.iloc[:, 2].values

# Find the ky slice closest to 0 (since np.linspace often skips exact 0.0)
closest_ky = ky_all[np.argmin(np.abs(ky_all))]
mask = np.isclose(ky_all, closest_ky)

kx = kx_all[mask]
ky = ky_all[mask]
growth = growth_all[mask]

# Calculate k^2
k_squared = kx**2 + ky**2

# --- Plot 1: Full View (Showing the Singularity) ---
# Filter out the microscopic numerical noise directly at the origin
valid_full = k_squared > 1e-4

plt.figure(figsize=(8, 6))
plt.semilogy(kx[valid_full], (growth / k_squared)[valid_full], '.')
plt.grid(True)
plt.title(r"110 km: $\gamma / k^2$ at $k_y \approx 0$")
plt.xlabel(r"$k_x$ (1/m)")
plt.ylabel(r"$\gamma / k^2$")
plt.tight_layout()
plt.savefig("gamma_110_a0_zero.png", dpi=300)
plt.show()