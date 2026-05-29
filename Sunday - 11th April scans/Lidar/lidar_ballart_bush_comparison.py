import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# CONFIG
# ==========================================================
LIDAR_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\Lidar\Bush_Scan_1.csv'

RADAR_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\gradar_logs\20260411-141145\pointcloud_20260411-141145_.csv'

FIGURES_DIR = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\Lidar\figures_python'
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    print(f"saved: {path}")
    plt.show()
    plt.close()

# Bush geometry
BUSH_FRONT = 1.5
BUSH_BACK  = 3.4
CR_RANGE   = 3.72
CR_TOL     = 0.3
R_MAX      = 6.0
BIN_WIDTH  = 0.25

# ==========================================================
# LOAD LIDAR
# ==========================================================
df_l = pd.read_csv(LIDAR_CSV)
df_l['X']     = pd.to_numeric(df_l.iloc[:, 8],  errors='coerce')
df_l['Y']     = pd.to_numeric(df_l.iloc[:, 9],  errors='coerce')
df_l['Z']     = pd.to_numeric(df_l.iloc[:, 10], errors='coerce')
df_l['Refl']  = pd.to_numeric(df_l.iloc[:, 11], errors='coerce')
df_l['range'] = np.sqrt(df_l['X']**2 + df_l['Y']**2 + df_l['Z']**2)
df_l = df_l[
    (df_l['Refl']  > 0)    &
    (df_l['Refl']  <= 130) &
    (df_l['range'] > 0)    &
    (df_l['range'] <= R_MAX)
].copy()

print(f"LiDAR points loaded: {len(df_l):,}")
print(f"LiDAR range:         {df_l['range'].min():.2f} to "
      f"{df_l['range'].max():.2f} m")

# ==========================================================
# LOAD RADAR
# ==========================================================
df_r = pd.read_csv(RADAR_CSV)
df_r = df_r[
    (df_r['amplitude'] >= -80) &
    (df_r['range']     <= R_MAX) &
    (df_r['range']     >= 1.45) &
    (np.degrees(df_r['roll']) >= 7.0)
].copy()

az = df_r['roll'].values
el = (np.pi / 2) + df_r['pitch'].values
r  = df_r['range'].values
df_r['x'] = r * np.cos(el) * np.cos(az)
df_r['y'] = r * np.cos(el) * np.sin(az)
df_r['z'] = r * np.sin(el)

print(f"Radar points loaded: {len(df_r):,}")

# ==========================================================
# PRINT SUMMARY
# ==========================================================
cr_sel   = (df_l['range'] >= CR_RANGE - CR_TOL) & \
           (df_l['range'] <= CR_RANGE + CR_TOL)
bush_sel = (df_l['range'] >= BUSH_FRONT) & \
           (df_l['range'] <= BUSH_BACK)
beyond   = df_l['range'] > BUSH_BACK

print('\n' + '='*50)
print('LIDAR CR WINDOW SUMMARY')
print('='*50)
print(f"Total points:              {len(df_l):,}")
print(f"Points in bush body:       {bush_sel.sum():,}")
print(f"Points beyond bush back:   {beyond.sum():,}")
print(f"Points in CR window:       {cr_sel.sum():,}")
if cr_sel.sum() > 0:
    print(f"CR window mean refl:       "
          f"{df_l.loc[cr_sel, 'Refl'].mean():.1f}")
    print(f"CR window max refl:        "
          f"{df_l.loc[cr_sel, 'Refl'].max():.1f}")
else:
    print("No returns in CR window")
print('='*50)

# ==========================================================
# FIGURE 1 — Plan view coloured by reflectivity
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 7))
sc = ax.scatter(df_l['X'], df_l['Y'],
                c=df_l['Refl'], s=3, cmap='jet',
                alpha=0.6, vmin=0, vmax=130)
plt.colorbar(sc, ax=ax, label='Reflectivity')
ax.axvline(BUSH_FRONT, color='green', lw=1.5, ls='--',
           label=f'Bush front ({BUSH_FRONT} m)')
ax.plot(CR_RANGE, 0, 'r*', markersize=15,
        label='CR position (3.72 m, 0 m)',
        zorder=6)
ax.axvline(BUSH_BACK,  color='green', lw=1.5, ls='-.',
           label=f'Bush back ({BUSH_BACK} m)')
ax.axvline(CR_RANGE,   color='red',   lw=1.5, ls='-',
           label=f'CR ({CR_RANGE} m)')
ax.set_xlabel('X (m)', fontsize=11)
ax.set_ylabel('Y (m)', fontsize=11)
ax.set_title('LiDAR Cherry Ballart Bush Scan — Plan View\n'
             'Coloured by Reflectivity',
             fontsize=11, fontweight='bold')
ax.set_aspect('equal')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
save('fig1_lidar_planview')

# ==========================================================
# FIGURE 2 — Side view coloured by reflectivity
# ==========================================================
fig, ax = plt.subplots(figsize=(11, 5))
sc = ax.scatter(df_l['X'], df_l['Z'],
                c=df_l['Refl'], s=3, cmap='jet',
                alpha=0.6, vmin=0, vmax=130)
plt.colorbar(sc, ax=ax, label='Reflectivity')
ax.axvline(BUSH_FRONT, color='green', lw=1.5, ls='--',
           label=f'Bush front ({BUSH_FRONT} m)')
ax.axvline(BUSH_BACK,  color='green', lw=1.5, ls='-.',
           label=f'Bush back ({BUSH_BACK} m)')
ax.axvline(CR_RANGE,   color='red',   lw=1.5, ls='-',
           label=f'CR ({CR_RANGE} m)')
ax.set_xlabel('X (m)', fontsize=11)
ax.set_ylabel('Z (m)', fontsize=11)
ax.set_title('LiDAR Cherry Ballart Bush Scan — Side View\n'
             'Coloured by Reflectivity',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
save('fig2_lidar_sideview')

# ==========================================================
# FIGURE 3 — Radar vs LiDAR normalised density vs range
# ==========================================================
edges = np.arange(0, R_MAX + BIN_WIDTH, BIN_WIDTH)
rmids = edges[:-1] + BIN_WIDTH / 2

# LiDAR density
counts_l = np.array([
    ((df_l['range'] >= e) & (df_l['range'] < e + BIN_WIDTH)).sum()
    for e in edges[:-1]], dtype=float)
counts_l_norm = counts_l / counts_l.max() \
    if counts_l.max() > 0 else counts_l

# Radar density
counts_r = np.array([
    ((df_r['range'] >= e) & (df_r['range'] < e + BIN_WIDTH)).sum()
    for e in edges[:-1]], dtype=float)
counts_r_norm = counts_r / counts_r.max() \
    if counts_r.max() > 0 else counts_r

fig, ax = plt.subplots(figsize=(11, 5))

ax.plot(rmids, counts_l_norm, color='steelblue', lw=2.5,
        label='LiDAR normalised density')
ax.plot(rmids, counts_r_norm, color='darkorange', lw=2.5,
        label='Radar normalised density')

ax.axvspan(BUSH_FRONT, BUSH_BACK, alpha=0.08, color='green',
           label=f'Bush body ({BUSH_FRONT}--{BUSH_BACK}\,m)')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.20, color='red',
           label=f'CR window ({CR_RANGE}\,m)')
ax.axvline(BUSH_FRONT, color='green', lw=1.5, ls='--')
ax.axvline(BUSH_BACK,  color='green', lw=1.5, ls='-.')
ax.axvline(CR_RANGE,   color='red',   lw=1.5, ls='-')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Normalised point density', fontsize=11)
ax.set_title('Radar vs LiDAR — Normalised Point Density vs Range\n'
             'Cherry Ballart Bush Scan',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, R_MAX)
plt.tight_layout()
save('fig3_radar_vs_lidar_density')