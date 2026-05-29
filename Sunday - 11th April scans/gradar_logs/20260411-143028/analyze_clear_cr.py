import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONFIG
# ==========================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH   = os.path.join(SCRIPT_DIR, 'pointcloud_20260411-143028_.csv')

CR_X = 3.449
CR_Y = 1.560
CR_Z = 0.032

AMP_MIN         = -80
R_MAX           = 6.0
NEAR_FIELD_CUT  = 1.45
FENCE_ROLL_DEG  = 7.0
CONE_HALF_ANGLE = 8.0
BIN_WIDTH       = 0.15
MIN_PTS_BIN     = 3      # lower for density plot — we want to see drop-off
MIN_PTS_STATS   = 8      # higher for variance stats

SHRUB_FRONT  = 1.5
SHRUB_BACK   = 3.4
FOLIAGE_BACK = 2.5
CR_TOL       = 0.15
CR_RADIUS    = 0.18

CR_RANGE = np.sqrt(CR_X**2 + CR_Y**2 + CR_Z**2)

FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'dense_face_methods')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    plt.show()
    plt.close()

# ==========================================================
# LOAD + FILTER
# ==========================================================
df = pd.read_csv(CSV_PATH)
df = df[
    (df['amplitude'] >= AMP_MIN) &
    (df['range']     <= R_MAX)   &
    (df['range']     >= NEAR_FIELD_CUT) &
    (np.degrees(df['roll']) >= FENCE_ROLL_DEG)
].copy()

az = df['roll'].values
el = (np.pi / 2) + df['pitch'].values
r  = df['range'].values
amp = df['amplitude'].values

df['x'] = r * np.cos(el) * np.cos(az)
df['y'] = r * np.cos(el) * np.sin(az)
df['z'] = r * np.sin(el)

# ==========================================================
# CONE FILTER
# ==========================================================
cr_vec    = np.array([CR_X, CR_Y, CR_Z])
boresight = cr_vec / np.linalg.norm(cr_vec)

pts      = df[['x','y','z']].values
pts_norm = pts / np.linalg.norm(pts, axis=1, keepdims=True)
cos_theta = np.clip(pts_norm @ boresight, -1, 1)
theta_deg = np.degrees(np.arccos(cos_theta))

cone_mask = theta_deg <= CONE_HALF_ANGLE
cone      = df[cone_mask].copy()

r_cone   = cone['range'].values
amp_cone = cone['amplitude'].values
x_cone   = cone['x'].values
y_cone   = cone['y'].values
z_cone   = cone['z'].values

# ==========================================================
# CR MASK
# ==========================================================
xyz_dist = np.sqrt(
    (x_cone - CR_X)**2 +
    (y_cone - CR_Y)**2 +
    (z_cone - CR_Z)**2
)

cr_mask = (
    (r_cone >= CR_RANGE - CR_TOL) &
    (r_cone <= CR_RANGE + CR_TOL) &
    (xyz_dist <= CR_RADIUS)
)

bush_mask = (
    (r_cone >= SHRUB_FRONT) &
    (r_cone <= SHRUB_BACK)  &
    (~cr_mask)
)

cr_peak = float(np.max(amp_cone[cr_mask])) if cr_mask.sum() > 0 else np.nan

# ==========================================================
# METHOD 2 — Point density vs depth
# ==========================================================
edges   = np.arange(SHRUB_FRONT, R_MAX + BIN_WIDTH, BIN_WIDTH)
density_rows = []
for e in edges[:-1]:
    sel = (r_cone >= e) & (r_cone < e + BIN_WIDTH) & (~cr_mask)
    density_rows.append({
        'r_mid': e + BIN_WIDTH / 2,
        'count': int(sel.sum())
    })
density_df = pd.DataFrame(density_rows)
max_count  = density_df['count'].max()
density_df['normalised'] = density_df['count'] / max_count

print('\n' + '='*50)
print('METHOD 2 — Point density vs depth')
print('='*50)
print(f"{'Range (m)':<12} {'Count':<8} {'Normalised'}")
for _, row in density_df.iterrows():
    print(f"  {row['r_mid']:.2f}    |  {int(row['count']):3d}  |  {row['normalised']:.3f}")

# ==========================================================
# METHOD 3 — Amplitude variance vs depth
# ==========================================================
variance_rows = []
for e in edges[:-1]:
    sel = (r_cone >= e) & (r_cone < e + BIN_WIDTH) & bush_mask
    if sel.sum() >= MIN_PTS_STATS:
        variance_rows.append({
            'r_mid':    e + BIN_WIDTH / 2,
            'n':        int(sel.sum()),
            'mean':     float(np.mean(amp_cone[sel])),
            'std':      float(np.std(amp_cone[sel])),
            'median':   float(np.median(amp_cone[sel]))
        })
variance_df = pd.DataFrame(variance_rows)

print('\n' + '='*50)
print('METHOD 3 — Amplitude variance vs depth')
print('='*50)
print(f"{'Range mid (m)':<15} {'N pts':<8} {'Mean (dB)':<12} {'Std (dB)'}")
for _, row in variance_df.iterrows():
    print(f"    {row['r_mid']:.2f}       |  {int(row['n']):3d}  |  {row['mean']:+.1f}    |  {row['std']:.1f}")

# ==========================================================
# FIGURE 1 — Method 2: density vs depth
# ==========================================================
fig, ax = plt.subplots(figsize=(11, 5))

ax.bar(density_df['r_mid'], density_df['count'],
       width=BIN_WIDTH * 0.85,
       color='seagreen', edgecolor='darkgreen', alpha=0.8,
       label='Point count in cone')

ax.axvline(SHRUB_FRONT,  color='green',  lw=1.5, ls='--',
           label=f'Bush front ({SHRUB_FRONT} m)')
ax.axvline(FOLIAGE_BACK, color='blue',   lw=1.5, ls='--',
           label=f'Foliage end ({FOLIAGE_BACK} m)')
ax.axvline(SHRUB_BACK,   color='orange', lw=1.5, ls='--',
           label=f'Bush back ({SHRUB_BACK} m)')
ax.axvline(CR_RANGE,     color='red',    lw=1.8, ls='-',
           label=f'CR ({CR_RANGE:.2f} m)')

# Shade regions
ax.axvspan(SHRUB_FRONT,  FOLIAGE_BACK, alpha=0.10, color='green',
           label='Front foliage')
ax.axvspan(FOLIAGE_BACK, SHRUB_BACK,   alpha=0.06, color='orange',
           label='Open interior + woody core')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.20, color='red', label='CR window')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Point count in cone', fontsize=11)
ax.set_title('Cherry Ballart Dense Face — Cone Return Density vs Depth\n'
             'Point count per 0.25 m range bin (CR excluded)',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
save('method2_density_vs_depth')

# ==========================================================
# FIGURE 2 — Method 3: mean and variance vs depth
# ==========================================================
fig, ax1 = plt.subplots(figsize=(11, 5))

ax2 = ax1.twinx()

# Mean amplitude on left axis
ax1.plot(variance_df['r_mid'], variance_df['mean'],
         'b-o', lw=2, markersize=7, label='Mean amplitude', zorder=5)
ax1.fill_between(variance_df['r_mid'],
                 variance_df['mean'] - variance_df['std'],
                 variance_df['mean'] + variance_df['std'],
                 color='blue', alpha=0.12, label='±1σ envelope')
ax1.set_ylabel('Mean amplitude (dB)', fontsize=11, color='blue')
ax1.tick_params(axis='y', labelcolor='blue')

# Std dev on right axis
ax2.plot(variance_df['r_mid'], variance_df['std'],
         'r-s', lw=2, markersize=7, label='Std deviation', zorder=4)
ax2.set_ylabel('Standard deviation (dB)', fontsize=11, color='red')
ax2.tick_params(axis='y', labelcolor='red')

# Region markers
ax1.axvline(SHRUB_FRONT,  color='green',  lw=1.5, ls='--',
            label=f'Bush front ({SHRUB_FRONT} m)')
ax1.axvline(FOLIAGE_BACK, color='blue',   lw=1.5, ls='--',
            label=f'Foliage end ({FOLIAGE_BACK} m)')
ax1.axvline(SHRUB_BACK,   color='orange', lw=1.5, ls='--',
            label=f'Bush back ({SHRUB_BACK} m)')

ax1.axvspan(SHRUB_FRONT,  FOLIAGE_BACK, alpha=0.08, color='green')
ax1.axvspan(FOLIAGE_BACK, SHRUB_BACK,   alpha=0.05, color='orange')

ax1.set_xlabel('Range (m)', fontsize=11)
ax1.set_title('Cherry Ballart Dense Face — Amplitude Mean and Variance vs Depth\n'
              'Layered structure: dense foliage → open interior → woody stem core',
              fontsize=11, fontweight='bold')
ax1.grid(True, alpha=0.3)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='lower right')

plt.tight_layout()
save('method3_variance_vs_depth')

print(f'\nFigures saved to: {FIGURES_DIR}')