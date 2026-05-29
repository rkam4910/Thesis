import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==========================================================
# CONFIG
# ==========================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(SCRIPT_DIR, 'pointcloud_20260411-130952_.csv')

# For Scan 1 we do NOT know the exact CR centroid
# so we use the boresight auto-detection from the top-N strongest returns
# and search for the CR by range only
TOP_N_FOR_BORESIGHT = 20

# Analysis settings
AMP_MIN        = -80
R_MAX          = 6.0
NEAR_FIELD_CUT = 1.45

# Bush geometry — full bush depth, no truncation
# We fit the entire bush body since the layered structure
# means front-only fitting is not meaningful
SHRUB_FRONT  = 1.5
SHRUB_BACK   = 4.5
CR_RANGE     = 3.72   # approximate CR range from Scan 2 geometry (same setup)

# Filters
FENCE_ROLL_CUT_DEG = 7.0
CONE_HALF_ANGLE    = 15.0   # wider cone for Scan 1 since boresight is auto-detected
CR_TOL             = 0.3
BIN_WIDTH          = 0.25
MIN_POINTS_BIN     = 5

# Local clutter: region just before CR
CLUTTER_START = 2.6
CLUTTER_END   = 3.4

# Output folder
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'cherry_ballart_scan1')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=160, bbox_inches='tight')
    print(f"saved: {path}")
    plt.show()
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
# ==========================================================
az = df['roll'].values
el = (np.pi / 2) + df['pitch'].values
r  = df['range'].values

df['x'] = r * np.cos(el) * np.cos(az)
df['y'] = r * np.cos(el) * np.sin(az)
df['z'] = r * np.sin(el)


# ==========================================================
# AUTO BORESIGHT from top-N strongest returns
# ==========================================================
top       = df.nlargest(TOP_N_FOR_BORESIGHT, 'amplitude')
mean_dir  = top[['x', 'y', 'z']].mean().values
boresight = mean_dir / np.linalg.norm(mean_dir)

top_peak_amp   = float(top['amplitude'].max())
top_peak_range = float(top.loc[top['amplitude'].idxmax(), 'range'])

print("=" * 70)
print("CHERRY BALLART SHRUB ANALYSIS — Scan 1 (higher density face)")
print("=" * 70)
print(f"Auto boresight (top-{TOP_N_FOR_BORESIGHT}): {boresight.round(3)}")
print(f"Strongest return in top-{TOP_N_FOR_BORESIGHT}: "
      f"{top_peak_amp:+.1f} dB @ {top_peak_range:.2f} m")


# ==========================================================
# CONE FILTER
# ==========================================================
pts       = df[['x', 'y', 'z']].values
pts_norm  = pts / np.linalg.norm(pts, axis=1, keepdims=True)
cos_theta = np.clip(pts_norm @ boresight, -1, 1)
theta_deg = np.degrees(np.arccos(cos_theta))

cone_mask = theta_deg <= CONE_HALF_ANGLE
cone      = df[cone_mask].copy()

r_cone   = cone['range'].values
amp_cone = cone['amplitude'].values

print(f"Cone filter (+/-{CONE_HALF_ANGLE} deg): "
      f"{cone_mask.sum():,} / {len(df):,} points kept")


# ==========================================================
# CR / CLUTTER / BUSH MASKS
# ==========================================================
cr_mask = (r_cone >= CR_RANGE - CR_TOL) & (r_cone <= CR_RANGE + CR_TOL)

clutter_mask = (
    (r_cone >= CLUTTER_START) &
    (r_cone <= CLUTTER_END)   &
    (~cr_mask)
)

bush_mask = (
    (r_cone >= SHRUB_FRONT) &
    (r_cone <= SHRUB_BACK)  &
    (~cr_mask)
)

cr_vals      = amp_cone[cr_mask]
clutter_vals = amp_cone[clutter_mask]
r_bush       = r_cone[bush_mask]
amp_bush     = amp_cone[bush_mask]


# ==========================================================
# DETECTION METRICS
# ==========================================================
cr_peak          = float(np.max(cr_vals))            if len(cr_vals)      else np.nan
cr_mean          = float(np.mean(cr_vals))           if len(cr_vals)      else np.nan
own_floor_median = float(np.median(clutter_vals))    if len(clutter_vals) else np.nan
own_floor_std    = float(np.std(clutter_vals))       if len(clutter_vals) else np.nan
own_floor_robust = float(robust_sigma(clutter_vals)) if len(clutter_vals) else np.nan
snr_medfloor     = cr_peak - own_floor_median        if not np.isnan(cr_peak) else np.nan

local_threshold_3sig = (own_floor_median + 3 * own_floor_robust
                        if not np.isnan(own_floor_median) else np.nan)

detected_above_floor  = bool(cr_peak > own_floor_median)     if not np.isnan(cr_peak) else False
detected_above_3sigma = bool(cr_peak > local_threshold_3sig) if not np.isnan(cr_peak) else False


# ==========================================================
# FULL BUSH ATTENUATION FIT
# Fit entire bush body SHRUB_FRONT to SHRUB_BACK
# Positive slope = signal increases with depth (woody core dominates)
# Negative slope = genuine attenuation through dense canopy
# Both are physically meaningful — report and explain whichever you get
# ==========================================================
profile_full_bush = range_binned_stats(
    r_bush, amp_bush,
    SHRUB_FRONT, SHRUB_BACK,
    BIN_WIDTH, MIN_POINTS_BIN
)

if len(profile_full_bush) >= 3:
    slope, intercept, rval, pval, se = linregress(
        profile_full_bush['r_mid'], profile_full_bush['amp_mean']
    )
    k_full  = slope
    r2_full = rval ** 2
else:
    slope = intercept = rval = pval = se = np.nan
    k_full = r2_full = np.nan


# ==========================================================
# FOLIAGE-ONLY ATTENUATION FIT
# Fit front foliage region only (1.5m to 2.5m)
# This isolates dense canopy attenuation before woody core effects
# ==========================================================
profile_foliage = range_binned_stats(
    r_bush, amp_bush,
    1.5, 2.5,  # Front foliage region only
    BIN_WIDTH, MIN_POINTS_BIN
)

if len(profile_foliage) >= 3:
    slope, intercept, rval, pval, se = linregress(
        profile_foliage['r_mid'], profile_foliage['amp_mean']
    )
    k_foliage  = slope
    r2_foliage = rval ** 2
else:
    slope = intercept = rval = pval = se = np.nan
    k_foliage = r2_foliage = np.nan


# ==========================================================
# ROC CURVE
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
auc   = (np.trapezoid(tpr_arr[order], fpr_arr[order])
         if len(fpr_arr) > 1 else np.nan)


# ==========================================================
# PRINT SUMMARY
# ==========================================================
print()
print(f"CR window points:           {int(cr_mask.sum())}")
print(f"Local clutter points:       {int(clutter_mask.sum())}")
print(f"Bush body points:           {int(bush_mask.sum())}")
print(f"Bush profile bins:          {len(profile_full_bush)}")
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
if not np.isnan(k_full):
    print(f"Full bush attenuation k:    {k_full:+.2f} dB/m  "
          f"(R² = {r2_full:.2f},  SE = {se:.2f},  p = {pval:.4f})")
    if k_full > 0:
        print("  NOTE: Positive slope — woody stem core returns dominate over")
        print("        foliage absorption. Layered structure effect.")
    else:
        print("  Negative slope — genuine foliage attenuation detected.")
else:
    print("Full bush attenuation k:    n/a — not enough bins")

if not np.isnan(k_foliage):
    print(f"Foliage attenuation k:      {k_foliage:+.2f} dB/m  "
          f"(R² = {r2_foliage:.2f})")
else:
    print("Foliage attenuation k:      n/a — not enough bins")
print()
print(f"ROC AUC:                    {auc:.3f}")
print(f"TPR at 3σ:                  {tpr_3sig:.3f}")
print(f"FPR at 3σ:                  {fpr_3sig:.3f}")
print("=" * 70)


# ==========================================================
# FIGURE 1 — Plan view
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 7))
sc = ax.scatter(cone['x'], cone['y'], c=cone['amplitude'],
                s=8, cmap='viridis', alpha=0.7, vmin=-80, vmax=20)

# Mark boresight direction at CR range
cr_point = boresight * CR_RANGE
ax.plot(cr_point[0], cr_point[1], '*', color='red', markersize=20,
        markeredgecolor='black', label=f'Boresight @ {CR_RANGE} m', zorder=6)

# Draw cone edges
bore_az = np.arctan2(boresight[1], boresight[0])
for angle in [-CONE_HALF_ANGLE, CONE_HALF_ANGLE]:
    edge_az = bore_az + np.radians(angle)
    r_edge  = np.linspace(0, R_MAX, 50)
    ax.plot(r_edge * np.cos(edge_az), r_edge * np.sin(edge_az),
            'r--', lw=1.0, alpha=0.6)

ax.plot(0, 0, 'ks', markersize=10, label='Radar')
plt.colorbar(sc, ax=ax, label='Amplitude (dB)')
ax.set_xlabel('X (m)'); ax.set_ylabel('Y (m)')
ax.set_title('Cherry Ballart — Plan View\nScan 1 (higher density face), auto boresight')
ax.grid(True, alpha=0.3); ax.set_aspect('equal'); ax.legend(fontsize=9)
plt.tight_layout()
save('fig1_plan_view')


# ==========================================================
# FIGURE 2 — Amplitude vs range with full bush fit
# ==========================================================
fig, ax = plt.subplots(figsize=(13, 6))

ax.scatter(r_cone, amp_cone, s=5, alpha=0.15, color='steelblue',
           label=f'All cone returns ({len(cone):,} pts)')
ax.scatter(r_bush, amp_bush, s=8, alpha=0.3, color='darkblue',
           label=f'Bush body returns ({int(bush_mask.sum()):,} pts)')

ax.errorbar(profile_full_bush['r_mid'], profile_full_bush['amp_mean'],
            yerr=profile_full_bush['amp_std'],
            fmt='o', color='navy', capsize=3, markersize=6, lw=1.2,
            label=f'Binned mean ± 1σ ({BIN_WIDTH} m bins)', zorder=5)

if not np.isnan(k_full):
    r_line = np.linspace(SHRUB_FRONT, SHRUB_BACK, 100)
    ax.plot(r_line, intercept + slope * r_line, 'r--', lw=2.0,
            label=f'Full bush fit: k = {k_full:+.2f} dB/m (R² = {r2_full:.2f})')

ax.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.08, color='steelblue',
           label=f'Bush body ({SHRUB_FRONT}–{SHRUB_BACK} m)')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.25, color='red', label=f'CR window (~{CR_RANGE} m)')
ax.axhline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axhline(local_threshold_3sig, color='purple', lw=1.4, ls='--',
           label=f'Local 3σ threshold ({local_threshold_3sig:+.1f} dB)')
if not np.isnan(cr_peak):
    ax.axhline(cr_peak, color='red', lw=1.4, ls=':',
               label=f'CR peak ({cr_peak:+.1f} dB)')

ax.set_xlim(0, R_MAX + 0.1); ax.set_ylim(-100, 70)
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Amplitude vs Range\n'
             'Scan 1 (higher density face), full bush attenuation fit',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper left', ncol=2)
ax.grid(True, alpha=0.3)

k_str = f"{k_full:+.2f} dB/m" if not np.isnan(k_full) else "n/a"
r2_str = f"R²={r2_full:.2f}" if not np.isnan(r2_full) else ""
summary = (
    f"CR peak: {cr_peak:+.1f} dB\n"
    f"Clutter median: {own_floor_median:+.1f} dB\n"
    f"SNR: {snr_medfloor:+.1f} dB\n"
    f"Full bush k = {k_str}  {r2_str}"
)
ax.text(0.985, 0.04, summary, transform=ax.transAxes, ha='right', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
plt.tight_layout()
save('fig2_amplitude_vs_range')


# ==========================================================
# FIGURE 3 — Amplitude vs range with foliage attenuation fit
# ==========================================================
fig, ax = plt.subplots(figsize=(13, 6))

ax.scatter(r_cone, amp_cone, s=5, alpha=0.15, color='steelblue',
           label=f'All cone returns ({len(cone):,} pts)')
ax.scatter(r_bush, amp_bush, s=8, alpha=0.3, color='darkblue',
           label=f'Bush body returns ({int(bush_mask.sum()):,} pts)')

ax.errorbar(profile_foliage['r_mid'], profile_foliage['amp_mean'],
            yerr=profile_foliage['amp_std'],
            fmt='o', color='navy', capsize=3, markersize=6, lw=1.2,
            label=f'Binned mean ± 1σ ({BIN_WIDTH} m bins)', zorder=5)

if not np.isnan(k_foliage):
    r_line = np.linspace(1.5, 2.5, 100)
    ax.plot(r_line, intercept + slope * r_line, 'r--', lw=2.0,
            label=f'Foliage fit: k = {k_foliage:+.2f} dB/m (R² = {r2_foliage:.2f})')

ax.axvspan(1.5, 2.5, alpha=0.08, color='steelblue',
           label=f'Foliage region (1.5–2.5 m)')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.25, color='red', label=f'CR window (~{CR_RANGE} m)')
ax.axhline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axhline(local_threshold_3sig, color='purple', lw=1.4, ls='--',
           label=f'Local 3σ threshold ({local_threshold_3sig:+.1f} dB)')
if not np.isnan(cr_peak):
    ax.axhline(cr_peak, color='red', lw=1.4, ls=':',
               label=f'CR peak ({cr_peak:+.1f} dB)')

ax.set_xlim(0, R_MAX + 0.1); ax.set_ylim(-100, 70)
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Amplitude vs Range\n'
             'Scan 1 (higher density face), foliage attenuation fit',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper left', ncol=2)
ax.grid(True, alpha=0.3)

k_str = f"{k_foliage:+.2f} dB/m" if not np.isnan(k_foliage) else "n/a"
r2_str = f"R²={r2_foliage:.2f}" if not np.isnan(r2_foliage) else ""
summary = (
    f"CR peak: {cr_peak:+.1f} dB\n"
    f"Clutter median: {own_floor_median:+.1f} dB\n"
    f"SNR: {snr_medfloor:+.1f} dB\n"
    f"Foliage k = {k_str}  {r2_str}"
)
ax.text(0.985, 0.04, summary, transform=ax.transAxes, ha='right', fontsize=10,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
plt.tight_layout()
save('fig3_amplitude_vs_range_foliage_fit')


# ==========================================================
# FIGURE 4 — Bush profile showing layered structure
# ==========================================================
fig, ax = plt.subplots(figsize=(11, 5))

ax.plot(profile_full_bush['r_mid'], profile_full_bush['amp_mean'],
        'o-', color='steelblue', lw=2, markersize=7,
        label='Binned mean amplitude')
ax.fill_between(profile_full_bush['r_mid'],
                profile_full_bush['amp_mean'] - profile_full_bush['amp_std'],
                profile_full_bush['amp_mean'] + profile_full_bush['amp_std'],
                color='lightblue', alpha=0.3, label='±1σ')

if not np.isnan(k_full):
    r_line = np.linspace(SHRUB_FRONT, SHRUB_BACK, 100)
    ax.plot(r_line, intercept + slope * r_line, 'r--', lw=2,
            label=f'Linear fit: k = {k_full:+.2f} dB/m  (R² = {r2_full:.2f})')

ax.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.08, color='steelblue',
           label='Bush body')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.2, color='red', label=f'CR window (~{CR_RANGE} m)')
ax.axhline(own_floor_median, color='black', lw=1.4,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')

ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Mean amplitude (dB)', fontsize=11)
ax.set_title('Cherry Ballart — Full Bush Amplitude Profile\n'
             'Scan 1 (higher density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
plt.tight_layout()
save('fig3_bush_profile')


# ==========================================================
# FIGURE 5 — Histogram: CR vs local clutter
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 5))

ax.hist(clutter_vals, bins=20, alpha=0.65, color='steelblue',
        edgecolor='navy', label=f'Local clutter ({int(clutter_mask.sum())} pts)')
if len(cr_vals) > 0:
    ax.hist(cr_vals, bins=12, alpha=0.82, color='red',
            edgecolor='darkred', label=f'CR window ({int(cr_mask.sum())} pts)')
ax.axvline(own_floor_median, color='black', lw=1.8,
           label=f'Local clutter median ({own_floor_median:+.1f} dB)')
ax.axvline(local_threshold_3sig, color='purple', lw=1.4, ls='--',
           label=f'Local 3σ threshold ({local_threshold_3sig:+.1f} dB)')
if not np.isnan(cr_peak):
    ax.axvline(cr_peak, color='red', lw=1.4, ls=':',
               label=f'CR peak ({cr_peak:+.1f} dB)')

ax.set_xlabel('Amplitude (dB)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Amplitude Distribution: Corner Reflector vs Local Clutter\n'
             'Cherry Ballart — Scan 1 (higher density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y'); ax.legend(fontsize=9)
plt.tight_layout()
save('fig4_histogram')


# ==========================================================
# FIGURE 6 — Cone return density vs range
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 5))

density = range_binned_stats(r_cone, amp_cone, NEAR_FIELD_CUT, R_MAX, BIN_WIDTH, 1)
ax.bar(density['r_mid'], density['n'],
       width=BIN_WIDTH * 0.85, color='steelblue', alpha=0.7, edgecolor='navy')
ax.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.10, color='steelblue',
           label=f'Bush body ({SHRUB_FRONT}–{SHRUB_BACK} m)')
ax.axvspan(CR_RANGE - CR_TOL, CR_RANGE + CR_TOL,
           alpha=0.25, color='red', label=f'CR window (~{CR_RANGE} m)')
ax.set_xlabel('Range (m)', fontsize=11)
ax.set_ylabel('Point count in cone', fontsize=11)
ax.set_title('Cone Return Density vs Range\n'
             'Cherry Ballart — Scan 1 (higher density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y'); ax.legend(fontsize=9)
plt.tight_layout()
save('fig5_cone_density')


# ==========================================================
# FIGURE 7 — ROC curve
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 6))

ax.plot(fpr_arr, tpr_arr, color='steelblue', lw=2,
        label=f'ROC curve (AUC = {auc:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Random classifier')
if not np.isnan(tpr_3sig):
    ax.plot(fpr_3sig, tpr_3sig, 'ro', markersize=10,
            label=f'3σ operating point\n(TPR={tpr_3sig:.2f}, FPR={fpr_3sig:.2f})')

ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Detectability: Corner Reflector vs Local Clutter\n'
             'Cherry Ballart — Scan 1 (higher density face)',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
plt.tight_layout()
save('fig6_roc')


# ==========================================================
# FIGURE 8 — Summary table
# ==========================================================
fig, ax_t = plt.subplots(figsize=(11, 5.5))
ax_t.axis('off')

k_str  = (f"{k_full:+.2f} dB/m  (R² = {r2_full:.2f},  SE = {se:.2f},  p = {pval:.4f})"
          if not np.isnan(k_full) else 'n/a')
k_foliage_str = (f"{k_foliage:+.2f} dB/m  (R² = {r2_foliage:.2f})"
                 if not np.isnan(k_foliage) else 'n/a')
det_str = 'Yes' if detected_above_3sigma else 'No'

rows = [
    ['Scan',                         'Scan 1 — higher density face'],
    ['CR range (approx.) (m)',        f'{CR_RANGE:.2f}'],
    ['CR peak (dB)',                   f'{cr_peak:+.1f}' if not np.isnan(cr_peak) else 'n/a'],
    ['CR mean (dB)',                   f'{cr_mean:+.1f}' if not np.isnan(cr_mean) else 'n/a'],
    ['Local clutter median (dB)',      f'{own_floor_median:+.1f}'],
    ['Local robust sigma (dB)',        f'{own_floor_robust:.1f}'],
    ['Local 3σ threshold (dB)',        f'{local_threshold_3sig:+.1f}'],
    ['SNR vs local median (dB)',       f'{snr_medfloor:+.1f}' if not np.isnan(snr_medfloor) else 'n/a'],
    ['Full bush fit region (m)',       f'{SHRUB_FRONT} to {SHRUB_BACK}'],
    ['Attenuation k — full bush',      k_str],
    ['Attenuation k — foliage',        k_foliage_str],
    ['Cone points',                    f'{len(cone):,}'],
    ['Bush profile bins',              f'{len(profile_full_bush)}'],
    ['Detected above floor',           'Yes' if detected_above_floor  else 'No'],
    ['Detected above 3σ',             det_str],
    ['ROC AUC',                        f'{auc:.3f}'],
    ['TPR at 3σ threshold',            f'{tpr_3sig:.3f}'],
    ['FPR at 3σ threshold',            f'{fpr_3sig:.3f}'],
]

tbl = ax_t.table(
    cellText=rows,
    colLabels=['Metric', 'Cherry Ballart — Scan 1 (higher density)'],
    loc='center', cellLoc='center'
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)
tbl.scale(1.2, 1.65)

for j in range(2):
    tbl[0, j].set_facecolor('#2c5f8a')
    tbl[0, j].set_text_props(color='white', fontweight='bold')
for i in range(1, len(rows) + 1):
    tbl[i, 1].set_facecolor('#dce8f5')

ax_t.set_title('Cherry Ballart Shrub — Scan 1 Detection Summary\n'
               'Full bush attenuation fit  |  Local clutter reference only',
               fontsize=11, fontweight='bold', pad=20)
plt.tight_layout()
save('fig7_summary_table')

print(f"\nAll figures saved to: {FIGURES_DIR}")