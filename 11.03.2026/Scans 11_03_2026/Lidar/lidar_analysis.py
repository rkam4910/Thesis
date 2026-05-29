import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Output folder ─────────────────────────────────────────────────────────
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'lidar')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Saved: {path}")

# ── Path resolution helper ────────────────────────────────────────────────
def find_file(filename, search_dirs):
    """Search multiple directories for a file and return first match."""
    for d in search_dirs:
        candidate = os.path.join(d, filename)
        if os.path.isfile(candidate):
            return candidate
    return None

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIRS = [
    SCRIPT_DIR,
    os.path.join(SCRIPT_DIR, '..'),
    os.path.join(SCRIPT_DIR, '..', '..'),
    os.path.join(SCRIPT_DIR, '..', '..', '..'),
    r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis',
    r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\11.03.2026\Scans 11_03_2026\Lidar',
]

LIDAR_CSV_NAME  = 'Reeds_Scan.csv'
RADAR_CSV_NAME  = 'pointcloud_20260311-114835_.csv'

CSV_PATH = find_file(LIDAR_CSV_NAME, BASE_DIRS)
if CSV_PATH is None:
    raise FileNotFoundError(
        f"\n[ERROR] Could not find '{LIDAR_CSV_NAME}' in any of:\n"
        + "\n".join(f"  {d}" for d in BASE_DIRS)
        + f"\n\nPlease copy '{LIDAR_CSV_NAME}' to the same folder as this script:\n  {SCRIPT_DIR}"
    )
print(f"LiDAR CSV found: {CSV_PATH}")

df   = pd.read_csv(CSV_PATH)
df.columns = df.columns.str.strip()
x    = df['X'].values
y    = df['Y'].values
z    = df['Z'].values
refl = df['Reflectivity'].values
r    = np.sqrt(x**2 + y**2 + z**2)

# ── Known CR ranges (from radar analysis) ─────────────────────────────────
MARCH_CR_RANGE = 10.92
APRIL_CR1_RANGE = 13.88
REED_FRONT = 4.1       # front face of reeds (m)

# ── Vegetation attenuation values from radar ──────────────────────────────
RADAR_VEG_LOSS_MARCH = 58.3   # dB
RADAR_VEG_LOSS_APRIL = 19.2   # dB (lower bound, multipath affected)

# ── Print key statistics ──────────────────────────────────────────────────
print("=" * 55)
print("REEDS LiDAR — KEY STATISTICS")
print("=" * 55)
print(f"Total points:               {len(r):,}")
print(f"Reed front face (radar):    {REED_FRONT} m")
print(f"Points 0–4 m (near/front):  {(r < 4).sum():,}")
print(f"Points 4–5 m (reed face):   {((r>=4)&(r<5)).sum():,}")
print(f"Points 5–8 m (behind reeds):{((r>=5)&(r<8)).sum():,}")
print(f"Points 8–10 m:              {((r>=8)&(r<10)).sum():,}")
print(f"Points at March CR ({MARCH_CR_RANGE} m ±0.5): "
      f"{((r>=MARCH_CR_RANGE-0.5)&(r<=MARCH_CR_RANGE+0.5)).sum()}")
print(f"Points at April CR ({APRIL_CR1_RANGE} m ±0.5): "
      f"{((r>=APRIL_CR1_RANGE-0.5)&(r<=APRIL_CR1_RANGE+0.5)).sum()}")
print(f"Max reflectivity beyond 8 m: {refl[r>8].max() if (r>8).any() else 'N/A'}")


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Point count vs range (showing sharp cutoff at reed face)
# ══════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(10, 5))

# Filter out near-field noise (<1 m — sensor self-returns)
r_plot = r[(r >= 1) & (r < 16)]

bin_edges   = np.arange(1, 16.5, 0.5)
counts, _   = np.histogram(r_plot, bins=bin_edges)
bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2

ax1.bar(bin_centres, counts, width=0.45, color='steelblue',
        edgecolor='white', linewidth=0.4, label='LiDAR returns')

# Reed front face
ax1.axvline(REED_FRONT, color='green', linestyle='--', lw=1.8,
            label=f'Reed front face ({REED_FRONT} m)')

# CR range markers
ax1.axvline(MARCH_CR_RANGE, color='red', linestyle=':', lw=1.8,
            label=f'March CR ({MARCH_CR_RANGE} m) — 0 returns')
ax1.axvline(APRIL_CR1_RANGE, color='orange', linestyle=':', lw=1.8,
            label=f'April CR ({APRIL_CR1_RANGE} m) — 0 returns')

# Annotate the dropoff — point to the 6–6.5 m bin
drop_idx = np.argmin(np.abs(bin_centres - 6.25))
ax1.annotate('Sharp drop-off\nbeyond reed face',
             xy=(bin_centres[drop_idx], counts[drop_idx]),
             xytext=(8.0, counts.max() * 0.55),
             fontsize=9, color='darkred',
             arrowprops=dict(arrowstyle='->', color='darkred'))

ax1.set_xlabel('Range from LiDAR sensor (m)', fontsize=12)
ax1.set_ylabel('Number of returns', fontsize=12)
ax1.set_title('LiDAR Point Return Count vs Range — Reeds Scan\n'
              'Corner reflector positions marked at 10.92 m and 13.88 m',
              fontsize=12, fontweight='bold')
ax1.legend(fontsize=9, loc='upper right')
ax1.grid(True, alpha=0.3, axis='y')
ax1.set_xlim(1, 16)

# Text box — bottom right, clear of legend
ax1.text(0.63, 0.08,
         f'Returns at March CR ({MARCH_CR_RANGE} m ± 0.5 m):  0\n'
         f'Returns at April CR ({APRIL_CR1_RANGE} m ± 0.5 m):  0\n'
         f'LiDAR penetration depth:  < 8 m',
         transform=ax1.transAxes, fontsize=9,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
save('fig1_lidar_point_count_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — Reflectivity vs range
# ══════════════════════════════════════════════════════════════════════════
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(13, 5))
fig2.suptitle('LiDAR Reflectivity — Reeds Scan', fontsize=13, fontweight='bold')

# Scatter: reflectivity vs range
mask_plot = r < 16
ax2a.scatter(r[mask_plot], refl[mask_plot],
             s=0.5, alpha=0.15, color='steelblue')
ax2a.axvline(REED_FRONT, color='green', linestyle='--', lw=1.5,
             label=f'Reed front ({REED_FRONT} m)')
ax2a.axvline(MARCH_CR_RANGE, color='red', linestyle=':', lw=1.5,
             label=f'March CR ({MARCH_CR_RANGE} m)')
ax2a.axvline(APRIL_CR1_RANGE, color='orange', linestyle=':', lw=1.5,
             label=f'April CR ({APRIL_CR1_RANGE} m)')
ax2a.set_xlabel('Range (m)', fontsize=11)
ax2a.set_ylabel('Reflectivity (0–255)', fontsize=11)
ax2a.set_title('Reflectivity vs Range', fontsize=11, fontweight='bold')
ax2a.legend(fontsize=8)
ax2a.grid(True, alpha=0.3)
ax2a.set_xlim(0, 16)

# Mean reflectivity per range bin
bin_edges2   = np.arange(0, 16.5, 0.5)
bin_centres2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2
mean_refl, std_refl = [], []
for i in range(len(bin_edges2) - 1):
    mask = (r >= bin_edges2[i]) & (r < bin_edges2[i+1])
    if mask.sum() > 5:
        mean_refl.append(np.mean(refl[mask]))
        std_refl.append(np.std(refl[mask]))
    else:
        mean_refl.append(np.nan)
        std_refl.append(np.nan)

mean_refl = np.array(mean_refl)
std_refl  = np.array(std_refl)

ax2b.plot(bin_centres2, mean_refl, 'b-o', markersize=4, lw=1.5,
          label='Mean reflectivity per 0.5 m bin')
ax2b.fill_between(bin_centres2,
                  mean_refl - std_refl, mean_refl + std_refl,
                  alpha=0.2, color='blue', label='±1σ')
ax2b.axvline(REED_FRONT, color='green', linestyle='--', lw=1.5,
             label=f'Reed front ({REED_FRONT} m)')
ax2b.axvline(MARCH_CR_RANGE, color='red', linestyle=':', lw=1.5,
             label=f'March CR ({MARCH_CR_RANGE} m) — no returns')
ax2b.axvline(APRIL_CR1_RANGE, color='orange', linestyle=':', lw=1.5,
             label=f'April CR ({APRIL_CR1_RANGE} m) — no returns')
ax2b.set_xlabel('Range (m)', fontsize=11)
ax2b.set_ylabel('Mean Reflectivity', fontsize=11)
ax2b.set_title('Mean Reflectivity per Range Bin', fontsize=11, fontweight='bold')
ax2b.legend(fontsize=8)
ax2b.grid(True, alpha=0.3)
ax2b.set_xlim(0, 16)

plt.tight_layout()
save('fig2_lidar_reflectivity_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — Radar vs LiDAR attenuation comparison
# ══════════════════════════════════════════════════════════════════════════
fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle('Radar vs LiDAR — Attenuation Through Reeds',
              fontsize=13, fontweight='bold')

# LiDAR: normalised point density as function of range (within reeds region)
bin_edges3   = np.arange(3, 12, 0.5)
bin_centres3 = (bin_edges3[:-1] + bin_edges3[1:]) / 2
counts3, _   = np.histogram(r[(r >= 3) & (r < 12)], bins=bin_edges3)

# Normalise to the peak count (at reed face)
peak_count    = counts3.max()
norm_density  = counts3 / peak_count

ax3a.bar(bin_centres3, norm_density, width=0.45,
         color='steelblue', edgecolor='white', linewidth=0.4,
         label='LiDAR normalised point density')
ax3a.axvline(REED_FRONT, color='green', linestyle='--', lw=1.8,
             label=f'Reed front ({REED_FRONT} m)')
ax3a.axvline(MARCH_CR_RANGE, color='red', linestyle=':', lw=1.8,
             label=f'March CR ({MARCH_CR_RANGE} m)')
ax3a.axvline(APRIL_CR1_RANGE, color='orange', linestyle=':', lw=1.8,
             label=f'April CR ({APRIL_CR1_RANGE} m)')
ax3a.set_xlabel('Range (m)', fontsize=11)
ax3a.set_ylabel('Normalised point density', fontsize=11)
ax3a.set_title('LiDAR — Normalised Point Density vs Range\n'
               '(relative to peak at reed face)',
               fontsize=11, fontweight='bold')
ax3a.legend(fontsize=8)
ax3a.grid(True, alpha=0.3, axis='y')

# Radar vs LiDAR bar chart comparison
categories  = ['March\n(partial occlusion)', 'April CR1\n(heavy occlusion)']
radar_loss  = [RADAR_VEG_LOSS_MARCH, RADAR_VEG_LOSS_APRIL]
lidar_label = ['Total (0 returns\nbeyond reed face)', 'Total (0 returns\nbeyond reed face)']

x_pos = np.arange(len(categories))
width = 0.35

bars_radar = ax3b.bar(x_pos - width/2, radar_loss, width,
                      color='tomato', edgecolor='black', lw=0.8,
                      label='Radar vegetation loss (dB)')

# LiDAR: plot as effectively "total" — use a large representative value
# Since LiDAR has 0 returns at CR range, effective attenuation is unmeasurable
# We show this as a hatched bar indicating "complete attenuation"
ax3b.bar(x_pos + width/2, [100, 100], width,
         color='steelblue', edgecolor='black', lw=0.8,
         hatch='//', alpha=0.6, label='LiDAR — complete attenuation\n(0 returns at CR range)')

# Annotate radar bars
for bar, val in zip(bars_radar, radar_loss):
    ax3b.text(bar.get_x() + bar.get_width()/2, val + 1,
              f'{val} dB', ha='center', fontsize=9, fontweight='bold')

ax3b.text(x_pos[0] + width/2, 50, 'No returns\ndetected',
          ha='center', fontsize=8, color='white', fontweight='bold')
ax3b.text(x_pos[1] + width/2, 50, 'No returns\ndetected',
          ha='center', fontsize=8, color='white', fontweight='bold')

ax3b.set_xticks(x_pos)
ax3b.set_xticklabels(categories, fontsize=10)
ax3b.set_ylabel('Vegetation attenuation (dB)', fontsize=11)
ax3b.set_title('Radar vs LiDAR — Vegetation Attenuation\n'
               'at Corner Reflector Ranges',
               fontsize=11, fontweight='bold')
ax3b.legend(fontsize=8)
ax3b.grid(True, alpha=0.3, axis='y')
ax3b.set_ylim(0, 115)
ax3b.text(0.02, 0.97,
          'Note: April radar value (19.2 dB) is a lower bound\n'
          'due to multipath from gravel surface beneath CR.',
          transform=ax3b.transAxes, fontsize=7.5,
          va='top', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
save('fig3_radar_vs_lidar_attenuation')
plt.show()

# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — Normalised point density standalone
# ══════════════════════════════════════════════════════════════════════════
fig4, ax4 = plt.subplots(figsize=(8, 5))

bin_edges4   = np.arange(1, 16.5, 0.5)
bin_centres4 = (bin_edges4[:-1] + bin_edges4[1:]) / 2
counts4, _   = np.histogram(r[(r >= 1) & (r < 16)], bins=bin_edges4)
norm4        = counts4 / counts4.max()

ax4.bar(bin_centres4, norm4, width=0.45,
        color='steelblue', edgecolor='white', linewidth=0.4,
        label='LiDAR normalised point density')
ax4.axvline(REED_FRONT, color='green', linestyle='--', lw=1.8,
            label=f'Reed front face ({REED_FRONT} m)')
ax4.axvline(MARCH_CR_RANGE, color='red', linestyle=':', lw=1.8,
            label=f'March CR ({MARCH_CR_RANGE} m) — 0 returns')
ax4.axvline(APRIL_CR1_RANGE, color='orange', linestyle=':', lw=1.8,
            label=f'April CR ({APRIL_CR1_RANGE} m) — 0 returns')

ax4.set_xlabel('Range from LiDAR sensor (m)', fontsize=12)
ax4.set_ylabel('Normalised point density', fontsize=12)
ax4.set_title('LiDAR Normalised Point Density vs Range — Reeds Scan\n'
              '(relative to peak at reed face)',
              fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3, axis='y')
ax4.set_xlim(1, 16)
ax4.set_ylim(0, 1.1)

ax4.text(0.63, 0.75,
         f'Returns at March CR ({MARCH_CR_RANGE} m ± 0.5 m):  0\n'
         f'Returns at April CR ({APRIL_CR1_RANGE} m ± 0.5 m):  0\n'
         f'LiDAR penetration depth:  < 8 m',
         transform=ax4.transAxes, fontsize=9,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
save('fig4_lidar_normalised_density_standalone')
plt.show()

# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 5 — Attenuation through reeds: LiDAR vs Radar (line graph)
#  X axis: depth beyond reed front face (m)
#  Y axis: normalised signal level (dB), both referenced to 0 dB at reed face
# ══════════════════════════════════════════════════════════════════════════

# ── Load March radar data ─────────────────────────────────────────────────
MARCH_CSV = find_file(RADAR_CSV_NAME, BASE_DIRS)
if MARCH_CSV is None:
    raise FileNotFoundError(
        f"\n[ERROR] Could not find '{RADAR_CSV_NAME}' in any of:\n"
        + "\n".join(f"  {d}" for d in BASE_DIRS)
        + f"\n\nPlease copy '{RADAR_CSV_NAME}' to the same folder as this script:\n  {SCRIPT_DIR}"
    )
print(f"Radar CSV found: {MARCH_CSV}")
df_mar    = pd.read_csv(MARCH_CSV)
r_mar     = df_mar['range'].values
amp_mar   = df_mar['amplitude'].values

# ── Parameters ────────────────────────────────────────────────────────────
REED_FRONT   = 4.1    # reed front face for both sensors (m)
CR_DEPTH     = MARCH_CR_RANGE - REED_FRONT   # 6.82 m depth of March CR
BIN_WIDTH    = 0.5
MAX_DEPTH    = 8.0
depth_edges  = np.arange(0, MAX_DEPTH + BIN_WIDTH, BIN_WIDTH)
depth_centres = (depth_edges[:-1] + depth_edges[1:]) / 2

# ── LiDAR: point density per depth bin, normalised, converted to dB ───────
lidar_depth = r - REED_FRONT   # depth beyond reed face
lidar_counts = []
for i in range(len(depth_edges) - 1):
    mask = (lidar_depth >= depth_edges[i]) & (lidar_depth < depth_edges[i+1])
    lidar_counts.append(mask.sum())

lidar_counts  = np.array(lidar_counts, dtype=float)
peak_lidar    = lidar_counts[0]   # reference: count at reed face bin (depth 0–0.5 m)
lidar_norm    = lidar_counts / peak_lidar
# Avoid log(0) — clamp minimum to very small number
lidar_norm_safe = np.where(lidar_norm > 0, lidar_norm, 1e-6)
lidar_dB      = 10 * np.log10(lidar_norm_safe)

# ── Radar: mean amplitude per depth bin in vegetation zone only ────────────
# The radar has returns in two physically distinct regions:
#   (1) Reed body (0–~3m depth): vegetation scatter, attenuates with depth
#   (2) CR at 6.82m depth: strong specular return, isolated target
# Connecting these two regions with a line is misleading (the gap is empty
# space, not signal). Plot them separately.
CR_TOL_RAD  = 0.3
radar_depth = r_mar - REED_FRONT

# Vegetation zone: exclude CR returns, use mean amplitude per bin
cr_mask     = np.abs(r_mar - MARCH_CR_RANGE) < CR_TOL_RAD
radar_veg_amp = []
for i in range(len(depth_edges) - 1):
    bin_mask = (radar_depth >= depth_edges[i]) & (radar_depth < depth_edges[i+1])
    veg_mask = bin_mask & ~cr_mask
    if veg_mask.sum() >= 1:
        radar_veg_amp.append(np.mean(amp_mar[veg_mask]))
    else:
        radar_veg_amp.append(np.nan)

radar_veg_amp = np.array(radar_veg_amp)
ref_amp       = radar_veg_amp[0]                       # reed-face bin = 0 dB
radar_veg_dB  = radar_veg_amp - ref_amp

# Plot all vegetation bins up to (but not including) the CR bin
cr_bin_idx    = np.argmin(np.abs(depth_centres - CR_DEPTH))
veg_mask_plot = np.zeros(len(depth_centres), dtype=bool)
veg_mask_plot[:cr_bin_idx] = True    # everything before CR depth

# CR isolated point: amplitude relative to same reference
cr_range_mask   = np.abs(r_mar - MARCH_CR_RANGE) < CR_TOL_RAD
if cr_range_mask.sum() > 0:
    cr_amp_dB = np.max(amp_mar[cr_range_mask]) - ref_amp
else:
    cr_amp_dB = np.nan

# LiDAR complete attenuation depth
first_zero_idx = np.argmax(lidar_norm == 0)
zero_depth     = depth_centres[first_zero_idx] if first_zero_idx > 0 else MAX_DEPTH

Y_MIN, Y_MAX = -55, 10   # visible axis range

# ── Plot ──────────────────────────────────────────────────────────────────
fig5, ax5 = plt.subplots(figsize=(10, 6))

# LiDAR line — continuous, shows gradual attenuation then complete block
ax5.plot(depth_centres, lidar_dB, 'b^-', markersize=7, lw=2,
         label='LiDAR (905 nm) — normalised point density (dB)')

# Radar vegetation zone — line only through the reed body (no gap)
ax5.plot(depth_centres[veg_mask_plot], radar_veg_dB[veg_mask_plot],
         'rs-', markersize=7, lw=2,
         label='Radar (120 GHz) — mean vegetation clutter (dB)')

# CR: clipped to Y_MAX so it sits at the top of the chart.
# The actual value (+43.8 dB) is way above the axis range — annotated explicitly.
if not np.isnan(cr_amp_dB):
    CR_DISPLAY_Y = Y_MAX - 1          # pin the star just inside the top axis
    # Vertical dotted stem from bottom to top of chart at CR depth
    ax5.vlines(CR_DEPTH, Y_MIN, CR_DISPLAY_Y, colors='darkred',
               linestyles='dotted', lw=1.5, alpha=0.7)
    # Shaded column to highlight CR position
    ax5.axvspan(CR_DEPTH - 0.25, CR_DEPTH + 0.25,
                alpha=0.08, color='red')
    # Star at top of chart
    ax5.plot(CR_DEPTH, CR_DISPLAY_Y, 'r*', markersize=20, zorder=5,
             label=f'Radar: CR detected at {CR_DEPTH:.1f} m depth')
    # Upward arrow on the star to indicate it extends beyond the axis
    ax5.annotate('', xy=(CR_DEPTH, CR_DISPLAY_Y + 0.5),
                 xytext=(CR_DEPTH, CR_DISPLAY_Y - 3),
                 arrowprops=dict(arrowstyle='->', color='darkred', lw=2))
    # Text label at the CR position
    ax5.text(CR_DEPTH + 0.15, CR_DISPLAY_Y - 5,
             f'CR peak: +{cr_amp_dB:.1f} dB\n(above vegetation ref.)\nLiDAR: 0 returns here',
             fontsize=8.5, color='darkred', fontweight='bold', va='top')

# LiDAR complete attenuation vertical marker
ax5.axvline(zero_depth, color='steelblue', linestyle='--', lw=1.5, alpha=0.7,
            label=f'LiDAR: complete attenuation at ~{zero_depth:.1f} m depth')

# Shade the "LiDAR blind zone" beyond the attenuation limit
ax5.axvspan(zero_depth, MAX_DEPTH, alpha=0.06, color='blue',
            label='_nolegend_')
ax5.text(zero_depth + 0.1, Y_MIN + 4, 'LiDAR\nblind zone',
         fontsize=8, color='steelblue', alpha=0.8, style='italic')

ax5.set_xlabel('Depth beyond reed front face (m)', fontsize=12)
ax5.set_ylabel('Relative signal level (dB)', fontsize=12)
ax5.set_title('Foliage Penetration — LiDAR vs Radar\n'
              'Both normalised to 0 dB at reed front face',
              fontsize=12, fontweight='bold')
ax5.legend(fontsize=8.5, loc='lower left')
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, MAX_DEPTH)
ax5.set_ylim(Y_MIN, Y_MAX)

ax5.text(0.02, 0.98,
         f'LiDAR penetration limit: ~{zero_depth:.1f} m depth\n'
         f'Radar CR detected at: {CR_DEPTH:.1f} m depth\n'
         f'Two-way vegetation loss (radar): {RADAR_VEG_LOSS_MARCH} dB',
         transform=ax5.transAxes, fontsize=8.5, va='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.tight_layout()
save('fig5_attenuation_lidar_vs_radar')
plt.show()

print("\nAll figures saved to:", FIGURES_DIR)