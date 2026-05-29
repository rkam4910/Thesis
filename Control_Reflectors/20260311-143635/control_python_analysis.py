import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

df = pd.read_csv(r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Control_Reflectors\20260311-143635\pointcloud_20260311-143635_.csv')

r = df['range'].values
amp = df['amplitude'].values
pitch = df['pitch'].values
roll = df['roll'].values

# Cartesian conversion
x = r * np.sin(-pitch) * np.cos(roll)
y = r * np.sin(-pitch) * np.sin(roll)
z = -r * np.cos(-pitch)

# ── Parameters ────────────────────────────────────────────────────────────
CR1_RANGE = 19.79
CR2_RANGE = 23.77
CR_TOL    = 0.5
R_MAX     = 50
BIN_WIDTH = 2.0
MIN_POINTS = 5

# ── Clutter floor — exclude both CR clusters ──────────────────────────────
mask_scene = (r < R_MAX) & \
             ~((r >= CR1_RANGE - CR_TOL) & (r <= CR1_RANGE + CR_TOL)) & \
             ~((r >= CR2_RANGE - CR_TOL) & (r <= CR2_RANGE + CR_TOL))

clutter_floor = np.mean(amp[mask_scene])
clutter_std   = np.std(amp[mask_scene])

# ── Range-binned means (free space only, no CR zones) ─────────────────────
bins = np.arange(0, R_MAX + BIN_WIDTH, BIN_WIDTH)
bin_centers = bins[:-1] + BIN_WIDTH / 2
bin_means, bin_stds, bin_counts = [], [], []

for i in range(len(bins) - 1):
    in_bin = mask_scene & (r >= bins[i]) & (r < bins[i+1])
    if in_bin.sum() >= MIN_POINTS:
        bin_means.append(np.mean(amp[in_bin]))
        bin_stds.append(np.std(amp[in_bin]))
        bin_counts.append(in_bin.sum())
    else:
        bin_means.append(np.nan)
        bin_stds.append(np.nan)
        bin_counts.append(0)

bin_means = np.array(bin_means)
bin_stds  = np.array(bin_stds)
valid = ~np.isnan(bin_means)

slope, intercept, r_val, p_val, se = stats.linregress(bin_centers[valid], bin_means[valid])

# ── CR statistics ─────────────────────────────────────────────────────────
cr1_mask = (r >= CR1_RANGE - CR_TOL) & (r <= CR1_RANGE + CR_TOL)
cr2_mask = (r >= CR2_RANGE - CR_TOL) & (r <= CR2_RANGE + CR_TOL)

cr1_peak = amp[cr1_mask].max()
cr2_peak = amp[cr2_mask].max()
cr1_snr  = cr1_peak - clutter_floor
cr2_snr  = cr2_peak - clutter_floor
cr1_pts  = cr1_mask.sum()
cr2_pts  = cr2_mask.sum()

print("=" * 52)
print("CONTROL SCAN — FREE SPACE BASELINE")
print("=" * 52)
print(f"Clutter floor:         {clutter_floor:.1f} dB (σ = {clutter_std:.1f} dB)")
print(f"Attenuation slope:     {slope:.3f} ± {se:.3f} dB/m")
print(f"R²:                    {r_val**2:.3f}")
print(f"p-value:               {p_val:.4f}")
print()
print(f"CR1 @ {CR1_RANGE}m:  peak={cr1_peak:.1f} dB,  SNR={cr1_snr:.1f} dB,  pts={cr1_pts}")
print(f"CR2 @ {CR2_RANGE}m:  peak={cr2_peak:.1f} dB,  SNR={cr2_snr:.1f} dB,  pts={cr2_pts}")

# ── FIGURE ────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)

# Plot 1: Amplitude vs Range
ax1 = fig.add_subplot(gs[0, :])
ax1.scatter(r[r < R_MAX], amp[r < R_MAX], s=1, alpha=0.15, color='steelblue', label='All returns')
ax1.axvspan(CR1_RANGE - CR_TOL, CR1_RANGE + CR_TOL, alpha=0.15, color='red')
ax1.axvspan(CR2_RANGE - CR_TOL, CR2_RANGE + CR_TOL, alpha=0.15, color='green')
ax1.axhline(clutter_floor, color='orange', linestyle='--', lw=1.5,
            label=f'Clutter floor ({clutter_floor:.1f} dB)')
ax1.errorbar(bin_centers[valid], bin_means[valid], yerr=bin_stds[valid],
             fmt='ko', markersize=4, capsize=3, label='Bin mean ± 1σ')
fit_x = np.linspace(0, R_MAX, 100)
ax1.plot(fit_x, slope * fit_x + intercept, 'k--', lw=1.5,
         label=f'Linear fit: {slope:.3f} dB/m  (R²={r_val**2:.2f})')
ax1.annotate(f'CR1\n{cr1_peak:.0f} dB', xy=(CR1_RANGE, cr1_peak),
             xytext=(CR1_RANGE + 1.5, cr1_peak - 8), fontsize=8,
             arrowprops=dict(arrowstyle='->', color='red'), color='red')
ax1.annotate(f'CR2\n{cr2_peak:.0f} dB', xy=(CR2_RANGE, cr2_peak),
             xytext=(CR2_RANGE + 1.5, cr2_peak - 12), fontsize=8,
             arrowprops=dict(arrowstyle='->', color='green'), color='green')
ax1.set_xlabel('Range (m)', fontsize=11)
ax1.set_ylabel('Amplitude (dB)', fontsize=11)
ax1.set_title('Control Scan — Amplitude vs Range (no foliage)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=8, loc='lower right')
ax1.set_xlim(0, R_MAX)
ax1.set_ylim(-105, 75)
ax1.grid(True, alpha=0.3)

textstr = (f'Free-space slope: {slope:.3f} ± {se:.3f} dB/m\n'
           f'R² = {r_val**2:.3f},  p = {p_val:.4f}\n'
           f'Clutter floor: {clutter_floor:.1f} dB\n'
           f'CR1 SNR: {cr1_snr:.1f} dB  |  CR2 SNR: {cr2_snr:.1f} dB')
ax1.text(0.02, 0.05, textstr, transform=ax1.transAxes, fontsize=8,
         va='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

# Plot 2: Amplitude histogram
ax2 = fig.add_subplot(gs[1, 0])
ax2.hist(amp[r < R_MAX], bins=40, color='steelblue', edgecolor='white', alpha=0.8)
ax2.axvline(clutter_floor, color='orange', linestyle='--', lw=1.5,
            label=f'Clutter floor ({clutter_floor:.1f} dB)')
ax2.axvline(cr1_peak, color='red', linestyle=':', lw=1.5,
            label=f'CR1 peak ({cr1_peak:.0f} dB)')
ax2.axvline(cr2_peak, color='green', linestyle=':', lw=1.5,
            label=f'CR2 peak ({cr2_peak:.0f} dB)')
ax2.set_xlabel('Amplitude (dB)', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title('Amplitude Distribution', fontsize=11, fontweight='bold')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# Plot 3: XY point cloud coloured by amplitude
ax3 = fig.add_subplot(gs[1, 1])
mask_plot = r < R_MAX
sc = ax3.scatter(x[mask_plot], y[mask_plot], c=amp[mask_plot],
                 s=1, alpha=0.4, cmap='plasma', vmin=-90, vmax=70)
cr1_x = np.mean(x[cr1_mask]); cr1_y = np.mean(y[cr1_mask])
cr2_x = np.mean(x[cr2_mask]); cr2_y = np.mean(y[cr2_mask])
ax3.plot(cr1_x, cr1_y, 'r*', markersize=12, label='CR1')
ax3.plot(cr2_x, cr2_y, 'g*', markersize=12, label='CR2')
ax3.set_xlabel('X (m)', fontsize=11)
ax3.set_ylabel('Y (m)', fontsize=11)
ax3.set_title('XY Point Cloud (coloured by amplitude)', fontsize=11, fontweight='bold')
plt.colorbar(sc, ax=ax3, label='Amplitude (dB)')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)
ax3.set_aspect('equal')

plt.suptitle('Control Scan — Free Space Baseline (11 March 2026, No Foliage)',
             fontsize=13, fontweight='bold', y=1.01)
plt.savefig('control_scan_attenuation.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nFigure saved.")