import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==========================================================
# CONFIG
# ==========================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(SCRIPT_DIR, 'pointcloud_20260411-141145_.csv')

# Known corner reflector centroid (Scan 2 — lower density face)
CR_X = 3.369311
CR_Y = 1.573796
CR_Z = 0.083744

# Analysis settings
AMP_MIN       = -80
R_MAX         = 6.0
NEAR_FIELD_CUT = 1.45

# Bush geometry
# SHRUB_FRONT: front face of the bush (where dense foliage begins)
# SHRUB_BACK: full bush back boundary (used for attenuation fit from front to back)
# FOLIAGE_BACK: end of the front foliage layer — kept for comparison
SHRUB_FRONT  = 1.5    # m — front face of dense foliage (from cone density plot)
SHRUB_BACK   = 3.4    # m — full bush back boundary (for full attenuation fit)
FOLIAGE_BACK = 2.5    # m — end of front foliage layer (for comparison)

CR_RANGE_SHRUB = np.sqrt(CR_X**2 + CR_Y**2 + CR_Z**2)

# Filters
FENCE_ROLL_CUT_DEG = 7.0
CONE_HALF_ANGLE    = 8.0
CR_TOL             = 0.15
BIN_WIDTH          = 0.25
MIN_POINTS_BIN     = 8

# Local clutter region: just before the CR, inside the bush
CLUTTER_START = 2.6
CLUTTER_END   = 3.4

# CR spatial gate
CR_RADIUS = 0.18

# Output folder
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'cherry_ballart_local_analysis')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    print(f"saved: {path}")
    plt.show()          # <-- show AFTER saving, before close
    plt.close()

def robust_sigma(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad

def range_binned_stats(rr, aa, lo, hi, bin_width, min_pts):
    edges = np.arange(lo, hi + bin_width, bin_width)
    rows = []
    for e in edges[:-1]:
        sel = (rr >= e) & (rr < e + bin_width)
        if sel.sum() >= min_pts:
            rows.append({
                'r_mid':      e + bin_width / 2,
                'n':          int(sel.sum()),
                'amp_mean':   float(np.mean(aa[sel])),
                'amp_std':    float(np.std(aa[sel])),
                'amp_median': float(np.median(aa[sel]))
            })
    return pd.DataFrame(rows)


# ==========================================================
# LOAD + FILTER
# ==========================================================
df = pd.read_csv(CSV_PATH)
df = df[
    (df['amplitude'] >= AMP_MIN) &
    (df['range']     <= R_MAX)   &
    (df['range']     >= NEAR_FIELD_CUT)
].copy()

if FENCE_ROLL_CUT_DEG > 0:
    df = df[np.degrees(df['roll']) >= FENCE_ROLL_CUT_DEG].copy()


# ==========================================================
# SPHERICAL -> CARTESIAN
# az = roll, el = pi/2 + pitch
# ==========================================================
az = df['roll'].values
el = (np.pi / 2) + df['pitch'].values
r  = df['range'].values

df['x'] = r * np.cos(el) * np.cos(az)
df['y'] = r * np.cos(el) * np.sin(az)
df['z'] = r * np.sin(el)


# ==========================================================
# CONE FILTER TOWARDS KNOWN CR DIRECTION
# ==========================================================
cr_vec    = np.array([CR_X, CR_Y, CR_Z], dtype=float)
boresight = cr_vec / np.linalg.norm(cr_vec)

pts       = df[['x', 'y', 'z']].values
pts_norm  = pts / np.linalg.norm(pts, axis=1, keepdims=True)
cos_theta = np.clip(pts_norm @ boresight, -1, 1)
theta_deg = np.degrees(np.arccos(cos_theta))

cone_mask = theta_deg <= CONE_HALF_ANGLE
cone      = df[cone_mask].copy()

r_cone   = cone['range'].values
amp_cone = cone['amplitude'].values


# ==========================================================
# CR / CLUTTER / FOLIAGE-ONLY MASKS
# ==========================================================
xyz_dist = np.sqrt(
    (cone['x'].values - CR_X)**2 +
    (cone['y'].values - CR_Y)**2 +
    (cone['z'].values - CR_Z)**2
)

# CR window: range gate + spatial radius around known centroid
cr_mask = (
    (r_cone >= CR_RANGE_SHRUB - CR_TOL) &
    (r_cone <= CR_RANGE_SHRUB + CR_TOL) &
    (xyz_dist <= CR_RADIUS)
)

# Local clutter: region just before the CR, excluding CR itself
clutter_mask = (
    (r_cone >= CLUTTER_START) &
    (r_cone <= CLUTTER_END)   &
    (~cr_mask)
)

# Front foliage only: used for attenuation fit
# Stops at FOLIAGE_BACK to exclude open interior and woody stem core
foliage_mask = (
    (r_cone >= SHRUB_FRONT)  &
    (r_cone <= FOLIAGE_BACK) &
    (~cr_mask)
)

# Full bush body: used for profile plot only (not fit)
bush_body_mask = (
    (r_cone >= SHRUB_FRONT) &
    (r_cone <= SHRUB_BACK)  &
    (~cr_mask)
)

cr_vals      = amp_cone[cr_mask]
clutter_vals = amp_cone[clutter_mask]
r_foliage    = r_cone[foliage_mask]
amp_foliage  = amp_cone[foliage_mask]
r_bush       = r_cone[bush_body_mask]
amp_bush     = amp_cone[bush_body_mask]


# ==========================================================
# DETECTION METRICS
# ==========================================================
cr_peak          = float(np.max(cr_vals))             if len(cr_vals)      else np.nan
cr_mean          = float(np.mean(cr_vals))            if len(cr_vals)      else np.nan
own_floor_mean   = float(np.mean(clutter_vals))       if len(clutter_vals) else np.nan
own_floor_median = float(np.median(clutter_vals))     if len(clutter_vals) else np.nan
own_floor_std    = float(np.std(clutter_vals))        if len(clutter_vals) else np.nan
own_floor_robust = float(robust_sigma(clutter_vals))  if len(clutter_vals) else np.nan

snr_meanfloor        = cr_peak - own_floor_mean   if not np.isnan(cr_peak) else np.nan
snr_medfloor         = cr_peak - own_floor_median if not np.isnan(cr_peak) else np.nan
local_threshold_3sig = (own_floor_median + 3 * own_floor_robust
                        if not np.isnan(own_floor_median) else np.nan)

detected_above_floor  = bool(cr_peak > own_floor_median)      if not np.isnan(cr_peak) else False
detected_above_3sigma = bool(cr_peak > local_threshold_3sig)  if not np.isnan(cr_peak) else False


# ==========================================================
# FULL BUSH ATTENUATION FIT
# Fit from SHRUB_FRONT to SHRUB_BACK — the complete bush depth
# This captures attenuation through the entire bush structure
# ==========================================================
profile_full_bush = range_binned_stats(
    r_bush, amp_bush,
    SHRUB_FRONT, SHRUB_BACK,
    BIN_WIDTH, MIN_POINTS_BIN
)

if len(profile_full_bush) >= 3:
    slope_fb, intercept_fb, rval_fb, pval_fb, se_fb = linregress(
        profile_full_bush['r_mid'], profile_full_bush['amp_mean']
    )
    k_full_bush = slope_fb
    r2_full_bush = rval_fb ** 2
else:
    slope_fb = intercept_fb = rval_fb = pval_fb = se_fb = np.nan
    k_full_bush = r2_full_bush = np.nan


# ==========================================================
# FRONT FOLIAGE ATTENUATION FIT (for comparison)
# Fit only SHRUB_FRONT to FOLIAGE_BACK — the region where
# amplitude genuinely decreases before the open interior
# ==========================================================
profile_foliage = range_binned_stats(
    r_foliage, amp_foliage,
    SHRUB_FRONT, FOLIAGE_BACK,
    BIN_WIDTH, MIN_POINTS_BIN
)

if len(profile_foliage) >= 3:
    slope_f, intercept_f, rval_f, pval_f, se_f = linregress(
        profile_foliage['r_mid'], profile_foliage['amp_mean']
    )
    k_foliage   = slope_f
    r2_foliage  = rval_f ** 2
else:
    slope_f = intercept_f = rval_f = pval_f = se_f = np.nan
    k_foliage = r2_foliage = np.nan


# ==========================================================
# FULL BUSH PROFILE (for plotting only — not fit)
# ==========================================================
profile_full = range_binned_stats(
    r_bush, amp_bush,
    SHRUB_FRONT, SHRUB_BACK,
    BIN_WIDTH, MIN_POINTS_BIN
)


# ==========================================================
# ROC CURVE
# Positive class = CR window points
# Negative class = local clutter points
# ==========================================================
all_scores = (np.concatenate([cr_vals, clutter_vals])
              if (len(cr_vals) + len(clutter_vals)) > 0 else np.array([]))
thresholds = (np.linspace(np.min(all_scores), np.max(all_scores), 200)
              if len(all_scores) > 0 else np.array([]))

tpr_list, fpr_list = [], []
for t in thresholds:
    tp = np.mean(cr_vals      >= t) if len(cr_vals)      else 0.0
    fp = np.mean(clutter_vals >= t) if len(clutter_vals) else 0.0
    tpr_list.append(tp)
    fpr_list.append(fp)

tpr_arr = np.array(tpr_list)
fpr_arr = np.array(fpr_list)

tpr_3sig = np.mean(cr_vals      >= local_threshold_3sig) if len(cr_vals)      else np.nan
fpr_3sig = np.mean(clutter_vals >= local_threshold_3sig) if len(clutter_vals) else np.nan

order = np.argsort(fpr_arr)
auc   = (np.trapz(tpr_arr[order], fpr_arr[order])
         if len(fpr_arr) > 1 else np.nan)


# ==========================================================
# PRINT SUMMARY
# ==========================================================
print("=" * 70)
print("CHERRY BALLART SHRUB ANALYSIS — Scan 2 (lower density face)")
print("=" * 70)
print(f"CR centroid:                ({CR_X:.6f}, {CR_Y:.6f}, {CR_Z:.6f}) m")
print(f"CR range:                   {CR_RANGE_SHRUB:.3f} m")
print(f"Cone half-angle:            {CONE_HALF_ANGLE} deg")
print(f"Cone points:                {len(cone):,}")
print(f"CR window points:           {int(cr_mask.sum())}")
print(f"Local clutter points:       {int(clutter_mask.sum())}")
print(f"Front foliage points:       {int(foliage_mask.sum())}")
print()
print(f"CR peak:                    {cr_peak:+.1f} dB")
print(f"CR mean:                    {cr_mean:+.1f} dB")
print(f"Local clutter median:       {own_floor_median:+.1f} dB")
print(f"Local clutter sigma:        {own_floor_std:.1f} dB")
print(f"Local robust sigma:         {own_floor_robust:.1f} dB")
print(f"SNR vs median floor:        {snr_medfloor:+.1f} dB")
print(f"Local 3σ threshold:         {local_threshold_3sig:+.1f} dB")
print(f"Detected above floor:       {detected_above_floor}")
print(f"Detected above 3σ:          {detected_above_3sigma}")
print()
print(f"Full bush fit region:       {SHRUB_FRONT} m to {SHRUB_BACK} m")
print(f"Full bush fit bins:         {len(profile_full_bush)}")
if not np.isnan(k_full_bush):
    print(f"Attenuation k (full bush):  {k_full_bush:+.2f} dB/m  "
          f"(R² = {r2_full_bush:.2f}, SE = {se_fb:.2f}, p = {pval_fb:.4f})")
else:
    print("Attenuation k (full bush):  n/a — not enough bins")
print()
print(f"Foliage fit region:         {SHRUB_FRONT} m to {FOLIAGE_BACK} m")
print(f"Foliage fit bins:           {len(profile_foliage)}")
if not np.isnan(k_foliage):
    print(f"Attenuation k (foliage):    {k_foliage:+.2f} dB/m  "
          f"(R² = {r2_foliage:.2f}, SE = {se_f:.2f}, p = {pval_f:.4f})")
else:
    print("Attenuation k (foliage):    n/a — not enough bins")
print()
print(f"ROC AUC:                    {auc:.3f}")
print(f"TPR at 3σ:                  {tpr_3sig:.3f}")
print(f"FPR at 3σ:                  {fpr_3sig:.3f}")
print("=" * 70)
print()
print("TIP: Full bush attenuation now covers SHRUB_FRONT to SHRUB_BACK")
print("     This provides complete bush depth attenuation analysis")
print("     Compare full bush vs foliage-only results for insights")
print("=" * 70)


# ==========================================================
# FIGURE 1 — Plan view with CR centroid
# ==========================================================
fig, ax = plt.subplots(figsize=(8, 7))
sc = ax.scatter(cone['x'], cone['y'], c=cone['amplitude'],
                s=10, cmap='viridis', alpha=0.75, vmin=-80, vmax=20)
ax.scatter(CR_X, CR_Y, s=200, c='red', marker='*',
           edgecolors='black', linewidths=1.2, label='CR centroid', zorder=6)
ax.text(CR_X + 0.08, CR_Y + 0.05,
        f"CR\n({CR_X:.3f}, {CR_Y:.3f})", color='red', fontsize=10)
ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)')
ax.set_title('Cherry Ballart — Plan View with Corner Reflector\nScan 2 (lower density face)')
ax.grid(True, alpha=0.3); ax.set_aspect('equal'); ax.legend()
plt.colorbar(sc, ax=ax, label='Amplitude (dB)')
plt.tight_layout()
save('fig1_plan_view_cr')


# ==========================================================
# FIGURE 2 — Full amplitude vs range with foliage fit
# ==========================================================
fig, ax = plt.subplots(figsize=(13, 6))

ax.scatter(r_cone, amp_cone, s=5, alpha=0.15, color='seagreen',
           label=f'All cone returns ({len(cone):,} pts)')
ax.scatter(r_foliage, amp_foliage, s=8, alpha=0.3, color='darkgreen',
           label=f'Front foliage returns ({len(r_foliage):,} pts)')

ax.errorbar(profile_foliage['r_mid'], profile_foliage['amp_mean'],
            yerr=profile_foliage['amp_std'],
            fmt='o', color='navy', capsize=3, markersize=6, lw=1.2,
            label=f'Foliage binned mean ± 1σ ({BIN_WIDTH} m bins)', zorder=5)

if not np.isnan(k_foliage):
    r_line = np.linspace(SHRUB_FRONT, FOLIAGE_BACK, 50)
    ax.plot(r_line, intercept_f + slope_f * r_line, 'r--', lw=2.0,
            label=f'Foliage fit: k = {k_foliage:+.2f} dB/m (R² = {r2_foliage:.2f})')

# Full bush profile (no fit line — just for context)
ax.errorbar(profile_full['r_mid'], profile_full['amp_mean'],
            yerr=profile_full['amp_std'],
            fmt='s', color='grey', capsize=2, markersize=4, lw=0.8, alpha=0.5,
            label=f'Full bush binned mean ± 1σ (context only)', zorder=3)

ax.axvspan(SHRUB_FRONT, FOLIAGE_BACK, alpha=0.12, color='green',
           label=f'Front foliage region ({SHRUB_FRONT}–{FOLIAGE_BACK} m)')
ax.axvspan(FOLIAGE_BACK, SHRUB_BACK, alpha=0.06, color='orange',
           label=f'Open interior + woody core ({FOLIAGE_BACK}–{SHRUB_BACK} m)')
ax.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.25, color='red',
           label=f'CR window ({CR_RANGE_SHRUB:.2f} m)')
ax.axhline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axhline(local_threshold_3sig, color='purple', lw=1.4, ls='--',
           label=f'Local 3σ threshold ({local_threshold_3sig:+.1f} dB)')
ax.axhline(cr_peak, color='red', lw=1.4, ls=':',
           label=f'CR peak ({cr_peak:+.1f} dB)')
ax.axvline(SHRUB_FRONT, color='green', lw=1.2, ls='--', alpha=0.7,
           label=f'Bush front face ({SHRUB_FRONT} m)')

ax.set_xlim(0, R_MAX + 0.1); ax.set_ylim(-100, 70)
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Amplitude vs Range\n'
             f'Attenuation fit through front foliage layer only ({SHRUB_FRONT}–{FOLIAGE_BACK} m)',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper left', ncol=2)
ax.grid(True, alpha=0.3)

summary = (
    f"CR peak: {cr_peak:+.1f} dB\n"
    f"Clutter median: {own_floor_median:+.1f} dB\n"
    f"SNR: {snr_medfloor:+.1f} dB\n"
    f"Foliage k = {k_foliage:+.2f} dB/m  (R²={r2_foliage:.2f})"
    if not np.isnan(k_foliage) else
    f"CR peak: {cr_peak:+.1f} dB\n"
    f"Clutter median: {own_floor_median:+.1f} dB\n"
    f"SNR: {snr_medfloor:+.1f} dB\n"
    f"Foliage k = n/a"
)
ax.text(0.985, 0.04, summary, transform=ax.transAxes, ha='right', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))

plt.tight_layout()
save('fig2_amplitude_vs_range_foliage_fit')


# ==========================================================
# FIGURE 3 — Front foliage profile (fit region only)
# ==========================================================
fig, ax = plt.subplots(figsize=(10, 5))

ax.plot(profile_foliage['r_mid'], profile_foliage['amp_mean'],
        'o-', color='navy', lw=2, markersize=7,
        label='Front foliage binned mean amplitude')
ax.fill_between(profile_foliage['r_mid'],
                profile_foliage['amp_mean'] - profile_foliage['amp_std'],
                profile_foliage['amp_mean'] + profile_foliage['amp_std'],
                color='skyblue', alpha=0.3, label='±1σ')

if not np.isnan(k_foliage):
    r_line = np.linspace(SHRUB_FRONT, FOLIAGE_BACK, 50)
    ax.plot(r_line, intercept_f + slope_f * r_line, 'r--', lw=2,
            label=f'Linear fit: k = {k_foliage:+.2f} dB/m  (R² = {r2_foliage:.2f})')

ax.axvspan(SHRUB_FRONT, FOLIAGE_BACK, alpha=0.10, color='green',
           label='Front foliage fit region')
ax.axhline(own_floor_median, color='black', lw=1.4, ls='-',
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Mean amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Front Foliage Attenuation Profile\n'
             f'Fit region: {SHRUB_FRONT}–{FOLIAGE_BACK} m  '
             f'(excludes open interior and woody stem core)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
plt.tight_layout()
save('fig3_foliage_attenuation_profile')


# ==========================================================
# FIGURE 4 — Full bush profile (for structural discussion)
# Shows layered structure: foliage / open interior / woody core
# ==========================================================
fig, ax = plt.subplots(figsize=(11, 5))

ax.plot(profile_full['r_mid'], profile_full['amp_mean'],
        'o-', color='darkgreen', lw=2, markersize=6,
        label='Full bush binned mean amplitude')
ax.fill_between(profile_full['r_mid'],
                profile_full['amp_mean'] - profile_full['amp_std'],
                profile_full['amp_mean'] + profile_full['amp_std'],
                color='lightgreen', alpha=0.3, label='±1σ')

if not np.isnan(k_foliage):
    r_line = np.linspace(SHRUB_FRONT, FOLIAGE_BACK, 50)
    ax.plot(r_line, intercept_f + slope_f * r_line, 'r--', lw=2,
            label=f'Front foliage fit: k = {k_foliage:+.2f} dB/m')

ax.axvspan(SHRUB_FRONT,  FOLIAGE_BACK, alpha=0.15, color='green',
           label=f'Front foliage ({SHRUB_FRONT}–{FOLIAGE_BACK} m)')
ax.axvspan(FOLIAGE_BACK, SHRUB_BACK,   alpha=0.10, color='orange',
           label=f'Open interior + woody core ({FOLIAGE_BACK}–{SHRUB_BACK} m)')
ax.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.25, color='red', label=f'CR window ({CR_RANGE_SHRUB:.2f} m)')
ax.axhline(own_floor_median, color='black', lw=1.4,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Mean amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Full Bush Amplitude Profile\n'
             'Layered structure: dense foliage → open interior → woody stem core',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
plt.tight_layout()
save('fig4_full_bush_profile')


# ==========================================================
# FIGURE 5 — Full Bush Attenuation Fit (Front to Back)
# ==========================================================
fig, ax = plt.subplots(figsize=(11, 6))

ax.errorbar(profile_full_bush['r_mid'], profile_full_bush['amp_mean'],
            yerr=profile_full_bush['amp_std'],
            fmt='o', color='darkblue', capsize=3, markersize=7, lw=1.5,
            label=f'Full bush binned mean ± 1σ ({BIN_WIDTH} m bins)', zorder=5)

if not np.isnan(k_full_bush):
    r_line = np.linspace(SHRUB_FRONT, SHRUB_BACK, 50)
    ax.plot(r_line, intercept_fb + slope_fb * r_line, 'r-', lw=2.5,
            label=f'Full bush fit: k = {k_full_bush:+.2f} dB/m (R² = {r2_full_bush:.2f})')

# Mark bush regions
ax.axvspan(SHRUB_FRONT, FOLIAGE_BACK, alpha=0.12, color='green',
           label=f'Front foliage ({SHRUB_FRONT}–{FOLIAGE_BACK} m)')
ax.axvspan(FOLIAGE_BACK, SHRUB_BACK, alpha=0.08, color='orange',
           label=f'Interior + woody core ({FOLIAGE_BACK}–{SHRUB_BACK} m)')

# Reference lines
ax.axhline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axhline(cr_peak, color='red', lw=1.4, ls=':',
           label=f'CR peak ({cr_peak:+.1f} dB)')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Mean amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Full Bush Attenuation Analysis\n'
             f'Complete depth fit: {SHRUB_FRONT}–{SHRUB_BACK} m (front face to back)',
             fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=9, loc='upper right')

summary_fb = (
    f"Full bush attenuation:\n"
    f"k = {k_full_bush:+.2f} dB/m\n"
    f"R² = {r2_full_bush:.2f}\n"
    f"Bins: {len(profile_full_bush)}\n"
    f"Range: {SHRUB_FRONT}–{SHRUB_BACK} m"
    if not np.isnan(k_full_bush) else
    f"Full bush attenuation:\n"
    f"k = n/a\n"
    f"Insufficient data"
)
ax.text(0.02, 0.98, summary_fb, transform=ax.transAxes, va='top', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9))

plt.tight_layout()
save('fig5_full_bush_attenuation_fit')


# ==========================================================
# FIGURE 6 — Histogram: CR vs local clutter
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 5))

ax.hist(clutter_vals, bins=20, alpha=0.65, color='seagreen',
        edgecolor='darkgreen', label=f'Local clutter ({int(clutter_mask.sum())} pts)')
ax.hist(cr_vals, bins=12, alpha=0.82, color='red',
        edgecolor='darkred', label=f'CR window ({int(cr_mask.sum())} pts)')
ax.axvline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axvline(local_threshold_3sig, color='purple', lw=1.4, ls='--',
           label=f'Local 3σ threshold ({local_threshold_3sig:+.1f} dB)')
ax.axvline(cr_peak, color='red', lw=1.4, ls=':',
           label=f'CR peak ({cr_peak:+.1f} dB)')
ax.set_xlabel('Amplitude (dB)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Amplitude Distribution: Corner Reflector vs Local Clutter\n'
             'Cherry Ballart — Scan 2 (lower density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y'); ax.legend(fontsize=9)
plt.tight_layout()
save('fig5_histogram_cr_vs_clutter')


# ==========================================================
# FIGURE 6 — Cone return density vs range
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 5))

density_profile = range_binned_stats(
    r_cone, amp_cone, NEAR_FIELD_CUT, R_MAX, BIN_WIDTH, 1
)
ax.bar(density_profile['r_mid'], density_profile['n'],
       width=BIN_WIDTH * 0.85, color='seagreen', alpha=0.7, edgecolor='darkgreen')
ax.axvspan(SHRUB_FRONT,  FOLIAGE_BACK, alpha=0.15, color='green',
           label=f'Front foliage ({SHRUB_FRONT}–{FOLIAGE_BACK} m)')
ax.axvspan(FOLIAGE_BACK, SHRUB_BACK,   alpha=0.08, color='orange',
           label=f'Open interior + woody core ({FOLIAGE_BACK}–{SHRUB_BACK} m)')
ax.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.25, color='red', label=f'CR window ({CR_RANGE_SHRUB:.2f} m)')
ax.axvline(SHRUB_FRONT, color='green', lw=1.5, ls='--',
           label=f'Bush front face ({SHRUB_FRONT} m)')
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Point count in cone', fontsize=11)
ax.set_title('Cherry Ballart — Cone Return Density vs Range\n'
             'Lower density face (CR visible)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y'); ax.legend(fontsize=9)
plt.tight_layout()
save('fig6_cone_density_vs_range')


# ==========================================================
# FIGURE 7 — ROC curve
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 6))

ax.plot(fpr_arr, tpr_arr, color='navy', lw=2,
        label=f'ROC curve (AUC = {auc:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Random classifier')
ax.plot(fpr_3sig, tpr_3sig, 'ro', markersize=10,
        label=f'3σ operating point\n(TPR={tpr_3sig:.2f}, FPR={fpr_3sig:.2f})')
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Detectability: Corner Reflector vs Local Clutter\n'
             'Cherry Ballart — Scan 2 (lower density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
plt.tight_layout()
save('fig7_roc_detectability')


# ==========================================================
# FIGURE 8 — Summary table
# ==========================================================
fig, ax_t = plt.subplots(figsize=(11, 5.5))
ax_t.axis('off')

k_str = (f"{k_foliage:+.2f} dB/m  (R² = {r2_foliage:.2f},  SE = {se_f:.2f})"
         if not np.isnan(k_foliage) else 'n/a')

rows = [
    ['CR centroid (x, y, z)',           f'({CR_X:.3f}, {CR_Y:.3f}, {CR_Z:.3f}) m'],
    ['Range to CR (m)',                  f'{CR_RANGE_SHRUB:.3f}'],
    ['CR peak (dB)',                     f'{cr_peak:+.1f}'],
    ['CR mean (dB)',                     f'{cr_mean:+.1f}'],
    ['Local clutter median (dB)',        f'{own_floor_median:+.1f}'],
    ['Local robust sigma (dB)',          f'{own_floor_robust:.1f}'],
    ['Local 3σ threshold (dB)',          f'{local_threshold_3sig:+.1f}'],
    ['SNR vs local median (dB)',         f'{snr_medfloor:+.1f}'],
    ['Foliage fit region (m)',           f'{SHRUB_FRONT} to {FOLIAGE_BACK}'],
    ['Attenuation k — front foliage',   k_str],
    ['Cone points',                      f'{len(cone):,}'],
    ['Detected above floor',             'Yes' if detected_above_floor  else 'No'],
    ['Detected above 3σ',               'Yes' if detected_above_3sigma else 'No'],
    ['ROC AUC',                         f'{auc:.3f}'],
    ['TPR at 3σ threshold',             f'{tpr_3sig:.3f}'],
    ['FPR at 3σ threshold',             f'{fpr_3sig:.3f}'],
]

tbl = ax_t.table(
    cellText=rows,
    colLabels=['Metric', 'Cherry Ballart — Scan 2 (lower density)'],
    loc='center', cellLoc='center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)
tbl.scale(1.2, 1.75)

for j in range(2):
    tbl[0, j].set_facecolor('#2c5f8a')
    tbl[0, j].set_text_props(color='white', fontweight='bold')
for i in range(1, len(rows) + 1):
    tbl[i, 1].set_facecolor('#e8f5e8')

ax_t.set_title('Cherry Ballart Shrub — Cone-Filtered Detection Summary\n'
               'Attenuation fit through front foliage layer only  |  Local clutter reference',
               fontsize=11, fontweight='bold', pad=20)
plt.tight_layout()
save('fig8_summary_table')

print(f"\nAll figures saved to: {FIGURES_DIR}")


# Add after bush_body_mask definition

# METHOD 3 — variance vs depth table
print('\nMETHOD 3 — Amplitude variance vs depth')
print(f"{'Range mid (m)':<15} {'N pts':<8} {'Mean (dB)':<12} {'Std (dB)'}")
for _, row in profile_full.iterrows():
    print(f"    {row['r_mid']:.2f}       |  {int(row['n']):3d}  |  {row['amp_mean']:+.1f}    |  {row['amp_std']:.1f}")