"""
Shrub Scan 11/04/2026 -- Detection Analysis vs Control Baseline
================================================================
Applies the same SNR / detection / vegetation-loss pipeline used for the
reeds scans to the 11 April 2026 shrub scan, referenced against the
11 March 2026 control scan.

Parameters at the top of the file are the only things you should need to
change if the shrub front / CR range need tweaking.

Figures are written to a 'figures/shrub_11apr/' subfolder next to this script.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────
# CONFIG -- edit these if the defaults are wrong for your scan
# ──────────────────────────────────────────────────────────────────────────

# File paths -- update CONTROL_CSV if your local path differs
CONTROL_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Control_Reflectors\20260311-143635\pointcloud_20260311-143635_.csv'
# Shrub CSV: by default, look for the CSV next to this script.
# If your CSV lives somewhere else, replace this with the full path.
SHRUB_CSV   = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'pointcloud_20260411-130952_.csv',
)

# Control scan CR ranges -- same as reeds analysis
CR1_CTRL = 19.79
CR2_CTRL = 23.77
CR_TOL   = 0.3           # tight window around each CR peak
R_MAX    = 50            # control-scan range cap for clutter floor estimate

# Shrub scan geometry -- adjust if you know the real values
SHRUB_FRONT    = 1.5     # m  -- front face of shrub
SHRUB_BACK     = 3.5     # m  -- back face of shrub (just in front of CR)
CR_RANGE_SHRUB = 3.72    # m  -- CR peak location (from top-amp cluster)
NEAR_FIELD_CUT = 1.45    # m  -- ignore returns closer than this (hardware)
R_MAX_SHRUB    = 5.0     # m  -- cut-off for all shrub-scan processing

# Output
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, 'figures', 'shrub_11apr')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Saved: {path}")


# ──────────────────────────────────────────────────────────────────────────
# 1. CONTROL BASELINE (same as reeds analysis)
# ──────────────────────────────────────────────────────────────────────────
df_ctrl = pd.read_csv(CONTROL_CSV)
r_c   = df_ctrl['range'].values
amp_c = df_ctrl['amplitude'].values

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

print("=" * 60)
print("CONTROL BASELINE  (from 11 March 2026 control scan)")
print("=" * 60)
print(f"Clutter floor:     {CTRL_CLUTTER_FLOOR:+6.1f} dB   (σ = {CTRL_CLUTTER_STD:.1f} dB)")
print(f"CR1 peak @ {CR1_CTRL} m: {CTRL_CR1_PEAK:+6.1f} dB")
print(f"CR1 SNR:           {CTRL_CR1_SNR:+6.1f} dB")
print(f"Fixed threshold:   {FIXED_THRESHOLD:+6.1f} dB  (used as detection threshold)")


# ──────────────────────────────────────────────────────────────────────────
# 2. FREE-SPACE RANGE CORRECTION
# ──────────────────────────────────────────────────────────────────────────
def range_correct(amp_ref, r_ref, r_target):
    """Predicted free-space amplitude at r_target given a reference amp at r_ref,
    using the two-way free-space propagation law (-40 log10(r2/r1))."""
    return amp_ref - 40 * np.log10(r_target / r_ref)


# ──────────────────────────────────────────────────────────────────────────
# 3. SHRUB SCAN ANALYSIS
# ──────────────────────────────────────────────────────────────────────────
df_s = pd.read_csv(SHRUB_CSV)
r_s   = df_s['range'].values
amp_s = df_s['amplitude'].values

# Near-field and far-field trim
valid = (r_s >= NEAR_FIELD_CUT) & (r_s <= R_MAX_SHRUB)
r_s   = r_s[valid]
amp_s = amp_s[valid]

# Masks
cr_mask       = (r_s >= CR_RANGE_SHRUB - CR_TOL) & (r_s <= CR_RANGE_SHRUB + CR_TOL)
behind_mask   = (r_s >= SHRUB_FRONT) & (r_s <= CR_RANGE_SHRUB + CR_TOL)
clutter_mask  = behind_mask & ~cr_mask

cr_peak        = float(amp_s[cr_mask].max())         if cr_mask.sum()      > 0 else np.nan
own_floor      = float(np.mean(amp_s[clutter_mask])) if clutter_mask.sum() > 0 else np.nan
own_floor_std  = float(np.std(amp_s[clutter_mask]))  if clutter_mask.sum() > 0 else np.nan

snr_own  = cr_peak - own_floor
snr_ctrl = cr_peak - CTRL_CLUTTER_FLOOR

# Free-space expected amplitude at the CR range
cr_expected = range_correct(CTRL_CR1_PEAK, CR1_CTRL, CR_RANGE_SHRUB)
veg_loss    = cr_expected - cr_peak

# Detectability
cr_detected       = cr_peak > FIXED_THRESHOLD
three_sigma_meets = snr_own > 3 * own_floor_std

# Threshold sweep
thresholds = np.arange(-70, 70, 5)
cr_det_list, fp_list = [], []
tpr_list,    fpr_list = [], []

total_cr      = int(cr_mask.sum())
total_clutter = int(clutter_mask.sum())

for t in thresholds:
    d  = int(np.sum(amp_s[cr_mask]      >= t))
    fp = int(np.sum(amp_s[clutter_mask] >= t))
    cr_det_list.append(d)
    fp_list.append(fp)
    tpr_list.append(d  / total_cr      if total_cr      else 0.0)
    fpr_list.append(fp / total_clutter if total_clutter else 0.0)

cr_det_list = np.array(cr_det_list)
fp_list     = np.array(fp_list)
tpr_arr     = np.array(tpr_list)
fpr_arr     = np.array(fpr_list)

cr_det_fixed = int(np.sum(amp_s[cr_mask]      >= FIXED_THRESHOLD))
fp_fixed     = int(np.sum(amp_s[clutter_mask] >= FIXED_THRESHOLD))
tpr_fixed    = cr_det_fixed / total_cr      if total_cr      else 0.0
fpr_fixed    = fp_fixed     / total_clutter if total_clutter else 0.0

print()
print("=" * 60)
print("SHRUB SCAN  (11/04/2026)")
print("=" * 60)
print(f"Behind-shrub points:   {int(behind_mask.sum())}")
print(f"CR-window points:      {total_cr}")
print(f"Clutter points:        {total_clutter}")
print(f"CR peak @ {CR_RANGE_SHRUB} m:      {cr_peak:+6.1f} dB")
print(f"Own clutter floor:     {own_floor:+6.1f} dB   (σ = {own_floor_std:.1f} dB)")
print(f"SNR (own floor):       {snr_own:+6.1f} dB")
print(f"SNR (control floor):   {snr_ctrl:+6.1f} dB")
print(f"Expected free-space:   {cr_expected:+6.1f} dB  @ {CR_RANGE_SHRUB} m")
print(f"Vegetation loss:       {veg_loss:+6.1f} dB")
print(f"Detected @ fixed thr:  {cr_detected}  (thr = {FIXED_THRESHOLD:.1f} dB)")
print(f"Meets 3σ criterion:    {three_sigma_meets}  (3σ = {3*own_floor_std:.1f} dB)")
print(f"TPR @ fixed threshold: {tpr_fixed:.2f}")
print(f"FPR @ fixed threshold: {fpr_fixed:.2f}")


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 1 -- Amplitude vs Range
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(r_s, amp_s, s=3, alpha=0.25, color='seagreen', label='All returns')

ax.axvspan(SHRUB_FRONT, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.10, color='orange', label='Behind-shrub region')
ax.axvspan(CR_RANGE_SHRUB - CR_TOL, CR_RANGE_SHRUB + CR_TOL,
           alpha=0.30, color='red', label=f'CR window ({CR_RANGE_SHRUB} m)')

ax.axvline(SHRUB_FRONT, color='darkgreen', ls='--', lw=1.3,
           label=f'Shrub front ({SHRUB_FRONT} m)')
ax.axvline(SHRUB_BACK,  color='darkgreen', ls=':',  lw=1.0,
           label=f'Shrub back ({SHRUB_BACK} m)')
ax.axhline(FIXED_THRESHOLD, color='black', lw=1.5,
           label=f'Control threshold ({FIXED_THRESHOLD:.1f} dB)')
ax.axhline(own_floor, color='orange', ls='--', lw=1.2,
           label=f'Own clutter floor ({own_floor:.1f} dB)')
ax.axhline(cr_peak, color='red', ls=':', lw=1.2,
           label=f'CR peak ({cr_peak:.1f} dB)')

ax.set_xlim(0, R_MAX_SHRUB + 0.5)
ax.set_ylim(-100, 70)
ax.set_xlabel('Range (m)',    fontsize=11)
ax.set_ylabel('Amplitude (dB)', fontsize=11)
ax.set_title('Shrub Scan 11/04/2026 -- Amplitude vs Range\n'
             'Corner reflector behind shrub, referenced to control baseline',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=8, loc='upper right')
ax.grid(True, alpha=0.3)

ax.text(0.02, 0.04,
        f"CR peak: {cr_peak:.1f} dB  |  SNR(own) = {snr_own:.1f} dB  |  "
        f"SNR(ctrl) = {snr_ctrl:.1f} dB\n"
        f"Veg loss: {veg_loss:.1f} dB  |  Own floor: {own_floor:.1f} dB  "
        f"(σ = {own_floor_std:.1f} dB)",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

plt.tight_layout()
save('fig1_shrub_amplitude_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 2 -- ROC curve + fixed-threshold operating point
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(fpr_arr, tpr_arr, 'o-', color='seagreen', lw=1.8, markersize=5,
        label='Shrub scan 11/04/2026')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4, label='Random classifier')
ax.plot(fpr_fixed, tpr_fixed, '*', color='darkred', markersize=16, zorder=5,
        label=f'Fixed threshold  (TPR={tpr_fixed:.2f}, FPR={fpr_fixed:.2f})')
ax.set_xlabel('False Positive Rate', fontsize=11)
ax.set_ylabel('True Positive Rate',  fontsize=11)
ax.set_title(f'ROC Curve -- Shrub Scan 11/04/2026\n'
             f'Fixed control threshold = {FIXED_THRESHOLD:.1f} dB  (★ = operating point)',
             fontsize=12, fontweight='bold')
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
save('fig2_shrub_roc')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 3 -- CR points detected & false positives vs threshold
# ══════════════════════════════════════════════════════════════════════════
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Shrub Scan 11/04/2026 -- Threshold Sweep',
             fontsize=13, fontweight='bold')

axL.plot(thresholds, cr_det_list, 'o-', color='seagreen', lw=1.8, markersize=5,
         label='CR points detected')
axL.axvline(FIXED_THRESHOLD, color='black', lw=1.5,
            label=f'Fixed threshold ({FIXED_THRESHOLD:.1f} dB)')
axL.set_xlabel('Detection threshold (dB)', fontsize=11)
axL.set_ylabel('CR points detected',       fontsize=11)
axL.set_title('CR returns above threshold', fontsize=11, fontweight='bold')
axL.legend(fontsize=9); axL.grid(True, alpha=0.3); axL.set_xlim(-70, 70)

axR.plot(thresholds, fp_list, 's-', color='tomato', lw=1.8, markersize=5,
         label='False positives')
axR.axvline(FIXED_THRESHOLD, color='black', lw=1.5,
            label=f'Fixed threshold ({FIXED_THRESHOLD:.1f} dB)')
axR.set_xlabel('Detection threshold (dB)', fontsize=11)
axR.set_ylabel('False positive points',    fontsize=11)
axR.set_title('Clutter returns above threshold\n(behind-shrub region)',
              fontsize=11, fontweight='bold')
axR.legend(fontsize=9); axR.grid(True, alpha=0.3); axR.set_xlim(-70, 70)

plt.tight_layout()
save('fig3_shrub_threshold_sweep')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 4 -- SNR bar chart (shrub vs control reference)
# ══════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(7, 5))
labels   = ['Control\nCR1 (free-space)', 'Shrub 11/04\nCR']
snr_vals = [CTRL_CR1_SNR, snr_ctrl]
colors   = ['steelblue', 'seagreen']
bars = ax.bar(labels, snr_vals, color=colors, edgecolor='black',
              linewidth=0.8, width=0.45)
ax.axhline(0, color='black', lw=1.2)

for bar, val in zip(bars, snr_vals):
    ypos = val + 1.5 if val >= 0 else val - 5
    ax.text(bar.get_x() + bar.get_width() / 2, ypos,
            f'{val:.1f} dB', ha='center', va='bottom',
            fontsize=10, fontweight='bold')

if snr_ctrl < 0:
    ax.text(1, snr_ctrl - 9, 'UNDETECTED', ha='center',
            fontsize=9, color='darkred', fontweight='bold')

ax.set_ylabel('SNR above control clutter floor (dB)', fontsize=11)
ax.set_title('SNR vs Control Clutter Floor\n'
             'Shrub 11/04/2026 compared with free-space control',
             fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
save('fig4_shrub_snr_vs_control')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 5 -- Summary table
# ══════════════════════════════════════════════════════════════════════════
fig, ax_t = plt.subplots(figsize=(11, 4.5))
ax_t.axis('off')

col_labels = ['Metric', 'Control', 'Shrub 11/04/2026']
rows = [
    ['Range (m)',                 f'{CR1_CTRL}',              f'{CR_RANGE_SHRUB}'],
    ['CR peak (dB)',              f'{CTRL_CR1_PEAK:+.1f}',    f'{cr_peak:+.1f}'],
    ['Expected free-space (dB)',  '—',                        f'{cr_expected:+.1f}'],
    ['Vegetation loss (dB)',      '0',                        f'{veg_loss:+.1f}'],
    ['Clutter floor (dB)',        f'{CTRL_CLUTTER_FLOOR:+.1f}', f'{own_floor:+.1f}'],
    ['σ clutter (dB)',            f'{CTRL_CLUTTER_STD:.1f}',  f'{own_floor_std:.1f}'],
    ['SNR vs own floor (dB)',     f'{CTRL_CR1_SNR:+.1f}',     f'{snr_own:+.1f}'],
    ['SNR vs control floor (dB)', f'{CTRL_CR1_SNR:+.1f}',     f'{snr_ctrl:+.1f}'],
    ['TPR @ fixed threshold',     '—',                        f'{tpr_fixed:.2f}'],
    ['FPR @ fixed threshold',     '—',                        f'{fpr_fixed:.2f}'],
    ['3σ criterion met?',         'Yes',                      'Yes' if three_sigma_meets else 'No'],
    ['Detected?',                 'Yes',                      'Yes' if cr_detected else 'No'],
]

tbl = ax_t.table(cellText=rows, colLabels=col_labels,
                 loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.scale(1.15, 1.8)

header_color = '#2c5f8a'
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor(header_color)
    tbl[0, j].set_text_props(color='white', fontweight='bold')

for i in range(1, len(rows) + 1):
    tbl[i, 1].set_facecolor('#f0f0f0')
    tbl[i, 2].set_facecolor('#e8f5e8')

ax_t.set_title('Shrub 11/04/2026 -- Detection Summary vs Control Baseline',
               fontsize=12, fontweight='bold', pad=18)
plt.tight_layout()
save('fig5_shrub_summary_table')
plt.show()

print()
print("=" * 60)
print(f"All figures saved to: {FIGURES_DIR}")
print("=" * 60)
