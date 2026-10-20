from root import *

def inverted_ratio(kx, ky, run_index=5): 
    numerator = 1j * ky * e * E_0 / m_i
    V_Ty_scalar = V_Ty[run_index]
    denominator = (kx**2 + ky**2) * V_Ty_scalar**2 + 1e-15 
    return np.abs(numerator / denominator)

kx_vals = np.linspace(-12.0, 12.0, 200)
ky_vals = np.linspace(-12.0, 12.0, 200)
KX, KY = np.meshgrid(kx_vals, ky_vals)

ratio_magnitude = inverted_ratio(KX, KY, run_index=5)

fig, ax = plt.subplots(figsize=(8, 6))

# Use a linear scale, but cap the max level at 2.0. 
# This makes the 1.0 boundary (Drive = Thermal) jump out clearly.
max_val = min(ratio_magnitude.max(), 2.0)
levels = np.linspace(0, max_val, 60)

contours = ax.contourf(KX, KY, ratio_magnitude, levels=levels, cmap="viridis", extend="max")
colorbar = fig.colorbar(contours, ax=ax)
colorbar.set_label(r"$\left|\frac{i\,\mathbf{k}\cdot\mathbf{a}_{i0}}{k_x^2 V_{T_x}^2 + k_y^2 V_{T_y}^2}\right|$")

ax.set_xlabel("kx (1/m)")
ax.set_ylabel("ky (1/m)")
ax.set_title("Inverted ratio")
ax.grid(True, linestyle=":", alpha=0.6)

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "inverted_ratio_linear.png"), dpi=300)
plt.show()