"""
Shrub Detection --- Improved version with cone filtering and cleaner plots
=========================================================================
Combines the strengths of:
  * reeds_detection_reference_to_control.py  (control-referenced SNR, vegetation loss)
  * bushAttenuationAnalysis.py               (cone extraction, range-binned profile)

What changed vs shrubDetection.py:
  1. Boresight is estimated from the *top-amplitude* cluster (the CR),
     not from the full-scan centroid.  The centroid-based boresight in
     bushAttenuationAnalysis.py gets pulled by near-field ground returns;
     this version points cleanly at the CR.
  2. All analysis happens inside a cone around that boresight
     (default 15 deg half-angle) so off-axis clutter doesn't contaminate
     the clutter-floor estimate or the amplitude-vs-range plot.
  3. The Amplitude vs Range plot now overlays a range-binned mean +- 1sigma
     profile on top of the raw cone points --- no more columns of
     hard-to-read scatter.
  4. Added a plan-view (X-Y) panel so you can see the cone filter is
     actually pointed at the CR, not at a random bit of ground.
  5. Added an attenuation-rate fit through the bush
     (shrub front -> CR front), so you get a dB/m number comparable to
     Cherry Ballart alongside the total vegetation loss.

Outputs go to figures/shrub_improved/ next to this script.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress


# ──────────────────────────────────────────────────────────────────────────
# CONFIG
# ──────────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# File paths
CONTROL_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Control_Reflectors\20260311-143635\pointcloud_20260311-143635_.csv'
SHRUB_CSV   = os.path.join(SCRIPT_DIR, 'pointcloud_20260411-130952_.csv')

# Control-scan CR ranges (same as the reeds analysis)
CR1_CTRL = 19.79
CR2_CTRL = 23.77
CR_TOL   = 0.3
R_MAX    = 50

# Shrub scan geometry
SHRUB_FRONT    = 1.5     # m  -- front face of shrub
SHRUB_BACK     = 3.4     # m  -- just before the CR
CR_RANGE_SHRUB = 3.72    # m  -- CR peak location (auto-verified below)
NEAR_FIELD_CUT = 1.45    # m  -- ignore < this (hardware / near-field)
R_MAX_SHRUB    = 5.0     # m  -- analysis cut-off

# Fence mask: drop returns at small azimuths (fence sits to one side of the
# scene in the shrub scan and dominates the top-amplitude list otherwise).
# Set FENCE_ROLL_CUT_DEG = 0.0 to disable.
FENCE_ROLL_CUT_DEG = 7.0  # drop points with roll < this (deg)

# Cone filter
CONE_HALF_ANGLE = 15.0   # deg
TOP_N_FOR_BORESIGHT = 20  # use mean direction of top-N strongest returns

# Manual boresight override (use when the CR is NOT detectable and the auto
# top-N detector would latch onto a random clutter pixel). Leave both as None
# to auto-detect. Values are the roll (azimuth) and pitch of the direction
# you actually pointed at the CR, in degrees.
MANUAL_BORESIGHT_ROLL_DEG  = None  # e.g. 25.0 for left-side CR
MANUAL_BORESIGHT_PITCH_DEG = None  # e.g. 0.0  for CR at radar height

# Range binning for the profile plot and attenuation fit
BIN_WIDTH      = 0.25    # m
MIN_POINTS_BIN = 5

# Output
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'shrub_improved')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  saved: {path}")


# ──────────────────────────────────────────────────────────────────────────
# 1. CONTROL BASELINE
# ──────────────────────────────────────────────────────────────────────────
df_ctrl = pd.read_csv(CONTROL_CSV)
r_c, amp_c = df_ctrl['range'].values, df_ctrl['amplitude'].values

mask_ctrl_bg = (
    (r_c < R_MAX)
    & ~((r_c >= CR1_CTRL - CR_TOL) & (r_c <= CR1_CTRL + CR_TOL))
    & ~((r_c >= CR2_CTRL - CR_TOL) & (r_c <= CR2_CTRL + CR_TOL))
)
CTRL_CLUTTER_FLOOR = float(np.mean(amp_c[mask_ctrl_bg]))
CTRL_CLUTTER_STD   = float(np.std(amp_c[mask_ctrl_bg]))
CTRL_CR1_PEAK      = float(amp_c[(r_c >= CR1_CTRL - CR_TOL)
                                 & (r_c <= CR1_CTRL + CR_TOL)].max())
CTRL_CR1_SNR       = CTRL_CR1_PEAK - CTRL_CLUTTER_FLOOR
FIXED_THRESHOLD    = CTRL_CLUTTER_FLOOR

print("=" * 66)
print("CONTROL BASELINE")
print("=" * 66)
print(f"  Clutter floor:     {CTRL_CLUTTER_FLOOR:+6.1f} dB  (sigma = {CTRL_CLUTTER_STD:.1f} dB)")
print(f"  CR1 peak @ {CR1_CTRL} m: {CTRL_CR1_PEAK:+6.1f} dB")
print(f"  CR1 SNR:           {CTRL_CR1_SNR:+6.1f} dB")
print(f"  Fixed threshold:   {FIXED_THRESHOLD:+6.1f} dB")


# ──────────────────────────────────────────────────────────────────────────
# 2. LOAD SHRUB SCAN + NEAR-FIELD/FAR-FIELD TRIM
# ──────────────────────────────────────────────────────────────────────────
df = pd.read_csv(SHRUB_CSV)
raw_n = len(df)
df = df[(df['range'] >= NEAR_FIELD_CUT) & (df['range'] <= R_MAX_SHRUB)].copy()
print(f"\nLoaded shrub scan: {raw_n:,} raw -> {len(df):,} points after trim "
      f"[{NEAR_FIELD_CUT} m, {R_MAX_SHRUB} m]")

# Apply the fence mask BEFORE computing the boresight, so the top-N
# amplitude cluster doesn't get pulled onto the fence.
if FENCE_ROLL_CUT_DEG > 0:
    before_fence = len(df)
    df = df[np.degrees(df['roll']) >= FENCE_ROLL_CUT_DEG].copy()
    print(f"Fence mask (roll >= {FENCE_ROLL_CUT_DEG} deg): "
          f"{before_fence:,} -> {len(df):,} points")


# ──────────────────────────────────────────────────────────────────────────
# 3. SPHERICAL -> CARTESIAN (same convention as bushAttenuationAnalysis.py)
#    az = roll,  el = pi/2 + pitch
# ──────────────────────────────────────────────────────────────────────────
az = df['roll'].values
el = (np.pi / 2) + df['pitch'].values
r  = df['range'].values
df['x'] = r * np.cos(el) * np.cos(az)
df['y'] = r * np.cos(el) * np.sin(az)
df['z'] = r * np.sin(el)


# ──────────────────────────────────────────────────────────────────────────
# 4. BORESIGHT (manual override OR top-N strongest returns)
# ──────────────────────────────────────────────────────────────────────────
if MANUAL_BORESIGHT_ROLL_DEG is not None and MANUAL_BORESIGHT_PITCH_DEG is not None:
    az_m = np.radians(MANUAL_BORESIGHT_ROLL_DEG)
    el_m = (np.pi / 2) + np.radians(MANUAL_BORESIGHT_PITCH_DEG)
    boresight = np.array([
        np.cos(el_m) * np.cos(az_m),
        np.cos(el_m) * np.sin(az_m),
        np.sin(el_m),
    ])
    cr_range_measured = CR_RANGE_SHRUB
    print(f"\nBoresight (MANUAL): {boresight.round(3)}   "
          f"(roll={MANUAL_BORESIGHT_ROLL_DEG} deg, pitch={MANUAL_BORESIGHT_PITCH_DEG} deg)")
else:
    top = df.nlargest(TOP_N_FOR_BORESIGHT, 'amplitude')
    mean_dir = top[['x', 'y', 'z']].mean().values
    boresight = mean_dir / np.linalg.norm(mean_dir)
    cr_range_measured = float(top['range'].median())
    top_peak_amp = float(top['amplitude'].max())

    print(f"\nBoresight (auto, from top-{TOP_N_FOR_BORESIGHT} strongest returns): "
          f"{boresight.round(3)}")
    print(f"  Median range of top-{TOP_N_FOR_BORESIGHT}: {cr_range_measured:.2f} m"
          f"   (config CR_RANGE_SHRUB = {CR_RANGE_SHRUB} m)")

    # Sanity check: if the strongest returns are weak, the CR probably isn't
    # really there and the boresight is being pulled to random clutter.
    if top_peak_amp < FIXED_THRESHOLD + 10:
        print(f"  WARNING: top amplitude ({top_peak_amp:+.1f} dB) is close to the "
              f"control threshold ({FIXED_THRESHOLD:+.1f} dB).")
        print( "           The CR is likely not detected. Consider setting "
               "MANUAL_BORESIGHT_ROLL_DEG / MANUAL_BORESIGHT_PITCH_DEG to the "
               "direction you pointed at the CR so the clutter stats are taken "
               "from the correct cone.")


# ──────────────────────────────────────────────────────────────────────────
# 5. CONE EXTRACTION
# ──────────────────────────────────────────────────────────────────────────
pts = df[['x', 'y', 'z']].values
pts_norm = pts / np.linalg.norm(pts, axis=1, keepdims=True)
cos_theta = np.clip(pts_norm @ boresight, -1, 1)
theta_deg = np.degrees(np.arccos(cos_theta))
cone_mask = theta_deg <= CONE_HALF_ANGLE
df['in_cone'] = cone_mask
cone = df[cone_mask].copy()

print(f"\nCone filter ({CONE_HALF_ANGLE} deg): "
      f"{cone_mask.sum():,} / {len(df):,} points kept")


# ──────────────────────────────────────────────────────────────────────────
# 6. DETECTION ANALYSIS ON CONE POINTS
# ──────────────────────────────────────────────────────────────────────────
r_cone   = cone['range'].values
amp_cone = cone['amplitude'].values

cr_mask      = (r_cone >= CR_RANGE_SHRUB - CR_TOL) & (r_cone <= CR_RANGE_SHRUB + CR_TOL)
behind_mask  = (r_cone >= SHRUB_FRONT) & (r_cone <= CR_RANGE_SHRUB + CR_TOL)
clutter_mask = behind_mask & ~cr_mask

cr_peak       = float(amp_cone[cr_mask].max())        if cr_mask.sum()      else np.nan
own_floor     = float(np.mean(amp_cone[clutter_mask])) if clutter_mask.sum() else np.nan
own_floor_std = float(np.std(amp_cone[clutter_mask]))  if clutter_mask.sum() else np.nan

snr_own   = cr_peak - own_floor
snr_ctrl  = cr_peak - CTRL_CLUTTER_FLOOR
cr_expect = CTRL_CR1_PEAK - 40 * np.log10(CR_RANGE_SHRUB / CR1_CTRL)
veg_loss  = cr_expect - cr_peak
cr_detected        = cr_peak > FIXED_THRESHOLD
meets_three_sigma  = snr_own > 3 * own_floor_std


# ──────────────────────────────────────────────────────────────────────────
# 7. RANGE-BINNED PROFILE (cone only) + ATTENUATION FIT THROUGH BUSH
# ──────────────────────────────────────────────────────────────────────────
def range_binned_stats(rr, aa, lo, hi, bin_width, min_pts):
    edges = np.arange(lo, hi + bin_width, bin_width)
    rows = []
    for e in edges[:-1]:
        sel = (rr >= e) & (rr < e + bin_width)
        if sel.sum() >= min_pts:
            rows.append({
                'r_mid':    e + bin_width / 2,
                'n':        int(sel.sum()),
                'amp_mean': float(np.mean(aa[sel])),
                'amp_std':  float(np.std(aa[sel])),
            })
    return pd.DataFrame(rows)

# Full profile across the analysis window
profile = range_binned_stats(r_cone, amp_cone,
                             NEAR_FIELD_CUT, R_MAX_SHRUB,
                             BIN_WIDTH, MIN_POINTS_BIN)

# Fit attenuation slope through the bush body only
# (shrub front to just before CR window)
fit_mask = (profile['r_mid'] >= SHRUB_FRONT) & \
           (profile['r_mid'] <= CR_RANGE_SHRUB - CR_TOL)
prof_fit = profile[fit_mask]

if len(prof_fit) >= 3:
    slope, intercept, rval, pval, se = linregress(prof_fit['r_mid'], prof_fit['amp_mean'])
    k_dB_per_m = slope
    fit_r2     = rval ** 2
    print(f"\nAttenuation fit through bush body "
          f"[{SHRUB_FRONT} -> {CR_RANGE_SHRUB - CR_TOL:.2f} m]:")
    print(f"  L(r) = {intercept:+.1f} + {k_dB_per_m:+.2f} * r  dB")
    print(f"  R^2 = {fit_r2:.3f},  p = {pval:.4f},  SE(slope) = {se:.2f}")
else:
    slope = intercept = rval = pval = se = np.nan
    k_dB_per_m = fit_r2 = np.nan
    print("\nAttenuation fit: not enough bins inside the bush body.")


# ──────────────────────────────────────────────────────────────────────────
# 8. PRINT DETECTION SUMMARY
# ──────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 66)
print(f"SHRUB SCAN 11/04/2026  --  cone +-{CONE_HALF_ANGLE} deg around CR")
print("=" * 66)
print(f"  Cone points:            {len(cone):,}")
print(f"  Behind-shrub cone pts:  {int(behind_mask.sum())}")
print(f"  CR-window pts:          {int(cr_mask.sum())}")
print(f"  Clutter pts:            {int(clutter_mask.sum())}")
print(f"  CR peak:                {cr_peak:+6.1f} dB  @ {CR_RANGE_SHRUB} m")
print(f"  Own clutter floor:      {own_floor:+6.1f} dB  (sigma = {own_floor_std:.1f} dB)")
print(f"  SNR (own):              {snr_own:+6.1f} dB")
print(f"  SNR (control):          {snr_ctrl:+6.1f} dB")
print(f"  Expected free-space:    {cr_expect:+6.1f} dB  @ {CR_RANGE_SHRUB} m")
print(f"  Vegetation loss:        {veg_loss:+6.1f} dB")
print(f"  k through bush body:    {k_dB_per_m:+.2f} dB/m  (R^2 = {fit_r2:.2f})"
      if not np.isnan(k_dB_per_m) else "  k through bush body:    n/a")
print(f"  Detected @ fixed thr:   {cr_detected}")
print(f"  Meets 3-sigma test:     {meets_three_sigma}  (3sigma = {3*own_floor_std:.1f} dB)")


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1 -- Amplitude vs Range (cone-filtered) with binned mean +- sigma
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(11, 6))

# Faded raw scatter (all in-cone points)
ax.scatter(r_cone, amp_cone, s=5, alpha=0.25, color='seagreen',
           label=f'Cone returns ({len(cone):,} pts)')

# Binned mean +- 1 sigma overlay
ax.errorbar(profile['r_mid'], profile['amp_mean'], yerr=profile['amp_std'],
            fmt='o', color='navy', capsize=3, markersize=5, lw=1.2,
            label=f'Binned mean +- 1sigma ({BIN_WIDTH} m bins)', zorder=4)

# Fit line through bush body
if not np.isnan(k_dB_per_m):
    r_line = np.linspace(SHRUB_FRONT, CR_RANGE_SHRUB - CR_TOL, 50)
    ax.plot(r_line, intercept + slope * r_line, 'r--', lw=1.8,
            label=f'Fit: k = {k_dB_per_m:+.2f} dB/m  (R^2 = {fit_r2:.2f})',
            zorder=5)

# Reference regions and thresholds
ax.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.10, color='orange',
           label='Shrub body')
ax.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.25, color='red', label=f'CR window ({CR_RANGE_SHRUB} m)')
ax.axhline(FIXED_THRESHOLD, color='black', lw=1.5,
           label=f'Control threshold ({FIXED_THRESHOLD:.1f} dB)')
ax.axhline(own_floor, color='orange', ls='--', lw=1.2,
           label=f'Own clutter floor ({own_floor:+.1f} dB)')
ax.axhline(cr_peak, color='red', ls=':', lw=1.2,
           label=f'CR peak ({cr_peak:+.1f} dB)')

ax.set_xlim(0, R_MAX_SHRUB + 0.25)
ax.set_ylim(-100, 70)
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Amplitude (dB)', fontsize=11)
ax.set_title('Shrub Scan 11/04/2026 -- Cone-Filtered Amplitude vs Range\n'
             f'Cone +-{CONE_HALF_ANGLE} deg around CR boresight, referenced to control baseline',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper left', ncol=2)
ax.grid(True, alpha=0.3)

ax.text(0.99, 0.04,
        f"CR peak: {cr_peak:+.1f} dB\n"
        f"SNR(own)  = {snr_own:+.1f} dB\n"
        f"SNR(ctrl) = {snr_ctrl:+.1f} dB\n"
        f"Veg loss  = {veg_loss:+.1f} dB\n"
        f"k through bush = {k_dB_per_m:+.2f} dB/m"
            if not np.isnan(k_dB_per_m) else
        f"CR peak: {cr_peak:+.1f} dB\n"
        f"SNR(own)  = {snr_own:+.1f} dB\n"
        f"SNR(ctrl) = {snr_ctrl:+.1f} dB\n"
        f"Veg loss  = {veg_loss:+.1f} dB",
        transform=ax.transAxes, fontsize=9, ha='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

plt.tight_layout()
save('fig1_cone_amplitude_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2 -- Plan view: cone filter sanity check
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 8))

# All points (grey)
ax.scatter(df['x'], df['y'], s=3, alpha=0.2, color='grey', label='All returns (trimmed)')
# In-cone points coloured by amplitude
sc = ax.scatter(cone['x'], cone['y'], c=cone['amplitude'], s=8, alpha=0.75,
                cmap='viridis', vmin=-80, vmax=40, label=f'Cone +-{CONE_HALF_ANGLE} deg')
cb = plt.colorbar(sc, ax=ax, label='Amplitude (dB)')

# CR location (projected boresight at CR range) and cone edges
cr_point = boresight * CR_RANGE_SHRUB
ax.plot(cr_point[0], cr_point[1], '*', color='red', markersize=22,
        markeredgecolor='black', label=f'CR boresight @ {CR_RANGE_SHRUB} m', zorder=6)

# Draw cone edges in XY (2D approximation)
cone_angles = np.array([-CONE_HALF_ANGLE, CONE_HALF_ANGLE])
bore_az = np.arctan2(boresight[1], boresight[0])
for a in cone_angles:
    edge_az = bore_az + np.radians(a)
    r_edge = np.linspace(0, R_MAX_SHRUB, 50)
    ax.plot(r_edge * np.cos(edge_az), r_edge * np.sin(edge_az),
            'r--', lw=1.0, alpha=0.7)

ax.plot(0, 0, 'ks', markersize=10, label='Radar')
ax.set_xlabel('X (m)', fontsize=11)
ax.set_ylabel('Y (m)', fontsize=11)
ax.set_title('Plan view -- cone filter sanity check',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')
plt.tight_layout()
save('fig2_plan_view')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3 -- Point density and amplitude histogram (cone)
# ══════════════════════════════════════════════════════════════════════════
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# Density per bin
axL.bar(profile['r_mid'], profile['n'], width=BIN_WIDTH * 0.85,
        color='seagreen', alpha=0.7, edgecolor='darkgreen', linewidth=0.5)
axL.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.15, color='orange',
            label='Shrub body')
axL.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
            alpha=0.25, color='red', label='CR window')
axL.set_xlabel('Range (m)', fontsize=11)
axL.set_ylabel('Point count (cone)', fontsize=11)
axL.set_title('Cone return density vs range', fontsize=11, fontweight='bold')
axL.legend(fontsize=9)
axL.grid(True, alpha=0.3, axis='y')

# Amplitude histogram (split into clutter vs CR-window)
axR.hist(amp_cone[clutter_mask], bins=30, alpha=0.65, color='seagreen',
         edgecolor='darkgreen', label=f'Clutter ({int(clutter_mask.sum())})')
axR.hist(amp_cone[cr_mask], bins=15, alpha=0.8, color='red',
         edgecolor='darkred', label=f'CR window ({int(cr_mask.sum())})')
axR.axvline(FIXED_THRESHOLD, color='black', lw=1.5, ls='-',
            label=f'Control threshold ({FIXED_THRESHOLD:+.1f} dB)')
axR.axvline(own_floor, color='orange', lw=1.2, ls='--',
            label=f'Own floor ({own_floor:+.1f} dB)')
axR.set_xlabel('Amplitude (dB)', fontsize=11)
axR.set_ylabel('Count', fontsize=11)
axR.set_title('Amplitude distribution (cone)', fontsize=11, fontweight='bold')
axR.legend(fontsize=9)
axR.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
save('fig3_density_and_histogram')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4 -- Summary table
# ══════════════════════════════════════════════════════════════════════════
fig, ax_t = plt.subplots(figsize=(11, 5))
ax_t.axis('off')

k_str = f"{k_dB_per_m:+.2f} (R^2 = {fit_r2:.2f})" if not np.isnan(k_dB_per_m) else '—'

rows = [
    ['Range to CR (m)',           f'{CR1_CTRL}',              f'{CR_RANGE_SHRUB}'],
    ['CR peak (dB)',              f'{CTRL_CR1_PEAK:+.1f}',    f'{cr_peak:+.1f}'],
    ['Expected free-space (dB)',  '—',                        f'{cr_expect:+.1f}'],
    ['Vegetation loss (dB)',      '0',                        f'{veg_loss:+.1f}'],
    ['k through bush (dB/m)',     '—',                        k_str],
    ['Clutter floor (dB)',        f'{CTRL_CLUTTER_FLOOR:+.1f}', f'{own_floor:+.1f}'],
    ['sigma clutter (dB)',        f'{CTRL_CLUTTER_STD:.1f}',  f'{own_floor_std:.1f}'],
    ['SNR vs own floor (dB)',     f'{CTRL_CR1_SNR:+.1f}',     f'{snr_own:+.1f}'],
    ['SNR vs control floor (dB)', f'{CTRL_CR1_SNR:+.1f}',     f'{snr_ctrl:+.1f}'],
    ['Cone points',               '—',                        f'{len(cone):,}'],
    ['3-sigma criterion met?',    'Yes',                      'Yes' if meets_three_sigma else 'No'],
    ['Detected above threshold?', 'Yes',                      'Yes' if cr_detected else 'No'],
]

tbl = ax_t.table(
    cellText=rows,
    colLabels=['Metric', 'Control', 'Shrub 11/04/2026 (cone-filtered)'],
    loc='center', cellLoc='center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.scale(1.15, 1.7)

for j in range(3):
    tbl[0, j].set_facecolor('#2c5f8a')
    tbl[0, j].set_text_props(color='white', fontweight='bold')
for i in range(1, len(rows) + 1):
    tbl[i, 1].set_facecolor('#f0f0f0')
    tbl[i, 2].set_facecolor('#e8f5e8')

ax_t.set_title('Shrub 11/04/2026 -- Cone-Filtered Detection Summary vs Control',
               fontsize=12, fontweight='bold', pad=18)
plt.tight_layout()
save('fig4_summary_table')
plt.show()

print("\n" + "=" * 66)
print(f"All figures saved to: {FIGURES_DIR}")
print("=" * 66)
