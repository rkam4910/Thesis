import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ── Output folder — creates a 'figures' subfolder next to this script ─────
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Saved: {path}")

# ── File paths ────────────────────────────────────────────────────────────
CONTROL_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Control_Reflectors\20260311-143635\pointcloud_20260311-143635_.csv'
MARCH_CSV   = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260311-114835_.csv'
APRIL_CSV   = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260408-131953_.csv'

# ── Parameters ────────────────────────────────────────────────────────────
CR1_CTRL         = 19.79
CR2_CTRL         = 23.77
CR_TOL           = 0.3      # tight so April CR1/CR2 don't bleed into each other
R_MAX            = 50

MARCH_REED_FRONT = 4.1
APRIL_REED_FRONT = 3.0

MARCH_CR1_RANGE  = 10.92
MARCH_CR2_RANGE  = 10.92   # same range bin in March
APRIL_CR1_RANGE  = 13.88
APRIL_CR2_RANGE  = 14.01

# ── Control baseline ──────────────────────────────────────────────────────
df_ctrl = pd.read_csv(CONTROL_CSV)
r_c     = df_ctrl['range'].values
amp_c   = df_ctrl['amplitude'].values

mask_ctrl_bg = (
    (r_c < R_MAX) &
    ~((r_c >= CR1_CTRL - CR_TOL) & (r_c <= CR1_CTRL + CR_TOL)) &
    ~((r_c >= CR2_CTRL - CR_TOL) & (r_c <= CR2_CTRL + CR_TOL))
)

CTRL_CLUTTER_FLOOR = np.mean(amp_c[mask_ctrl_bg])
CTRL_CLUTTER_STD   = np.std(amp_c[mask_ctrl_bg])
CTRL_CR1_PEAK      = amp_c[(r_c >= CR1_CTRL - CR_TOL) & (r_c <= CR1_CTRL + CR_TOL)].max()
CTRL_CR1_SNR       = CTRL_CR1_PEAK - CTRL_CLUTTER_FLOOR
FIXED_THRESHOLD    = CTRL_CLUTTER_FLOOR

print("=" * 55)
print("CONTROL BASELINE")
print("=" * 55)
print(f"Clutter floor (fixed threshold): {CTRL_CLUTTER_FLOOR:.1f} dB  (σ = {CTRL_CLUTTER_STD:.1f} dB)")
print(f"CR1 peak:                        {CTRL_CR1_PEAK:.1f} dB at {CR1_CTRL} m")
print(f"CR1 SNR:                         {CTRL_CR1_SNR:.1f} dB")

def range_correct(amp_ref, r_ref, r_target):
    return amp_ref - 40 * np.log10(r_target / r_ref)

def fmt(val, decimals=1):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return '—'
    return f'{val:.{decimals}f}'

# ── Parse scan with separate CR1 / CR2 ───────────────────────────────────
def parse_scan(filepath, cr1_range, cr2_range, reed_front):
    df    = pd.read_csv(filepath)
    r     = df['range'].values
    amp   = df['amplitude'].values

    same_range   = abs(cr1_range - cr2_range) < 0.05
    far_cr       = max(cr1_range, cr2_range)
    mask_behind  = (r >= reed_front) & (r <= far_cr + CR_TOL)

    # If the two CRs are close enough that ±CR_TOL windows would overlap,
    # split at the midpoint so each CR gets its own exclusive range window
    if abs(cr1_range - cr2_range) < 2 * CR_TOL and not same_range:
        midpoint = (cr1_range + cr2_range) / 2
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <  midpoint)
        cr2_mask = (r >= midpoint)            & (r <= cr2_range + CR_TOL)
    else:
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <= cr1_range + CR_TOL)
        cr2_mask = (r >= cr2_range - CR_TOL) & (r <= cr2_range + CR_TOL)

    both_cr_mask = cr1_mask | cr2_mask
    clutter_mask = mask_behind & ~both_cr_mask

    cr1_peak = amp[cr1_mask].max() if cr1_mask.sum() > 0 else np.nan
    cr2_peak = amp[cr2_mask].max() if cr2_mask.sum() > 0 else np.nan

    own_floor     = np.mean(amp[clutter_mask]) if clutter_mask.sum() > 0 else np.nan
    own_floor_std = np.std(amp[clutter_mask])  if clutter_mask.sum() > 0 else np.nan

    cr1_snr_own  = cr1_peak - own_floor
    cr2_snr_own  = cr2_peak - own_floor
    cr1_snr_ctrl = cr1_peak - CTRL_CLUTTER_FLOOR
    cr2_snr_ctrl = cr2_peak - CTRL_CLUTTER_FLOOR

    cr1_expected  = range_correct(CTRL_CR1_PEAK, CR1_CTRL, cr1_range)
    cr2_expected  = range_correct(CTRL_CR1_PEAK, CR1_CTRL, cr2_range)
    cr1_veg_loss  = cr1_expected - cr1_peak
    cr2_veg_loss  = cr2_expected - cr2_peak if not np.isnan(cr2_peak) else np.nan

    thresholds = np.arange(-70, 70, 5)
    cr1_det_list, cr2_det_list, fp_list = [], [], []
    tpr1_list, tpr2_list, fpr_list      = [], [], []

    total_cr1     = cr1_mask.sum()
    total_cr2     = cr2_mask.sum()
    total_clutter = clutter_mask.sum()

    for t in thresholds:
        d1 = np.sum(amp[cr1_mask]      >= t)
        d2 = np.sum(amp[cr2_mask]      >= t)
        fp = np.sum(amp[clutter_mask]  >= t)
        cr1_det_list.append(d1)
        cr2_det_list.append(d2)
        fp_list.append(fp)
        tpr1_list.append(d1 / total_cr1      if total_cr1      > 0 else 0)
        tpr2_list.append(d2 / total_cr2      if total_cr2      > 0 else 0)
        fpr_list.append(fp  / total_clutter  if total_clutter  > 0 else 0)

    cr1_det_fixed = np.sum(amp[cr1_mask]     >= FIXED_THRESHOLD)
    cr2_det_fixed = np.sum(amp[cr2_mask]     >= FIXED_THRESHOLD)
    fp_fixed      = np.sum(amp[clutter_mask] >= FIXED_THRESHOLD)
    tpr1_fixed    = cr1_det_fixed / total_cr1     if total_cr1     > 0 else 0
    tpr2_fixed    = cr2_det_fixed / total_cr2     if total_cr2     > 0 else 0
    fpr_fixed     = fp_fixed      / total_clutter if total_clutter > 0 else 0

    return dict(
        r=r, amp=amp,
        mask_behind=mask_behind, cr1_mask=cr1_mask, cr2_mask=cr2_mask,
        clutter_mask=clutter_mask,
        cr1_peak=cr1_peak,       cr2_peak=cr2_peak,
        cr1_snr_own=cr1_snr_own, cr2_snr_own=cr2_snr_own,
        cr1_snr_ctrl=cr1_snr_ctrl, cr2_snr_ctrl=cr2_snr_ctrl,
        cr1_expected=cr1_expected, cr2_expected=cr2_expected,
        cr1_veg_loss=cr1_veg_loss, cr2_veg_loss=cr2_veg_loss,
        own_floor=own_floor, own_floor_std=own_floor_std,
        thresholds=thresholds,
        tpr1=np.array(tpr1_list), tpr2=np.array(tpr2_list),
        fpr=np.array(fpr_list),
        cr1_det=np.array(cr1_det_list), cr2_det=np.array(cr2_det_list),
        fp=np.array(fp_list),
        total_cr1=total_cr1, total_cr2=total_cr2, total_clutter=total_clutter,
        cr1_det_fixed=cr1_det_fixed, cr2_det_fixed=cr2_det_fixed,
        fp_fixed=fp_fixed,
        tpr1_fixed=tpr1_fixed, tpr2_fixed=tpr2_fixed, fpr_fixed=fpr_fixed,
        cr1_detected=(cr1_peak > FIXED_THRESHOLD),
        cr2_detected=(cr2_peak > FIXED_THRESHOLD),
        cr1_range=cr1_range, cr2_range=cr2_range, reed_front=reed_front,
    )

mar = parse_scan(MARCH_CSV, MARCH_CR1_RANGE, MARCH_CR2_RANGE, MARCH_REED_FRONT)
apr = parse_scan(APRIL_CSV, APRIL_CR1_RANGE, APRIL_CR2_RANGE, APRIL_REED_FRONT)

for label, s in [('MARCH', mar), ('APRIL', apr)]:
    print(f"\n{'='*55}\n{label} — CR1 vs CR2\n{'='*55}")
    print(f"Own clutter floor:  {s['own_floor']:.1f} dB  (σ = {s['own_floor_std']:.1f} dB)")
    print(f"CR1: peak={s['cr1_peak']:.1f} dB  veg_loss={s['cr1_veg_loss']:.1f} dB  detected={s['cr1_detected']}")
    print(f"CR2: peak={s['cr2_peak']:.1f} dB  veg_loss={fmt(s['cr2_veg_loss'])} dB  detected={s['cr2_detected']}")
    print(f"@ fixed threshold:  TPR_CR1={s['tpr1_fixed']:.2f}  TPR_CR2={s['tpr2_fixed']:.2f}  FPR={s['fpr_fixed']:.2f}")


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Amplitude vs Range (March and April, side by side)
# ══════════════════════════════════════════════════════════════════════════
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
fig1.suptitle('Reeds Scans — Amplitude vs Range\nMarch (partial occlusion) vs April (heavy occlusion)',
              fontsize=13, fontweight='bold')

for ax, s, title, dot_color in [
    (ax1, mar, 'March — CRs elevated, partial occlusion',   'steelblue'),
    (ax2, apr, 'April — CRs at ground level, heavy occlusion', 'tomato'),
]:
    ax.scatter(s['r'][s['r'] < 20], s['amp'][s['r'] < 20],
               s=1, alpha=0.2, color=dot_color, label='All returns')
    ax.axvspan(s['reed_front'], max(s['cr1_range'], s['cr2_range']) + CR_TOL,
               alpha=0.08, color='orange', label='Behind-foliage region')
    ax.axvspan(s['cr1_range'] - CR_TOL, s['cr1_range'] + CR_TOL,
               alpha=0.3, color='red', label=f'CR1 ({s["cr1_range"]} m)')
    same = abs(s['cr1_range'] - s['cr2_range']) < 0.05
    cr2_lbl = f'CR2 ({s["cr2_range"]} m) — same range as CR1' if same else f'CR2 ({s["cr2_range"]} m)'
    ax.axvspan(s['cr2_range'] - CR_TOL, s['cr2_range'] + CR_TOL,
               alpha=0.3, color='blue', label=cr2_lbl)
    ax.axvline(s['reed_front'], color='green', linestyle='--', lw=1.3,
               label=f"Reed front ({s['reed_front']} m)")
    ax.axhline(CTRL_CLUTTER_FLOOR, color='black', linestyle='-', lw=1.5,
               label=f'Control threshold ({CTRL_CLUTTER_FLOOR:.1f} dB)')
    ax.axhline(s['cr1_peak'], color='red',  linestyle=':', lw=1.2,
               label=f'CR1 peak ({s["cr1_peak"]:.1f} dB)')
    ax.axhline(s['cr2_peak'], color='blue', linestyle=':', lw=1.2,
               label=f'CR2 peak ({s["cr2_peak"]:.1f} dB)')
    ax.set_xlim(0, 20); ax.set_ylim(-100, 90)
    ax.set_xlabel('Range (m)', fontsize=11)
    ax.set_ylabel('Amplitude (dB)', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, alpha=0.3)
    cr2_label = f'{s["cr2_peak"]:.1f} dB  ✗ UNDETECTED' if not s['cr2_detected'] else f'{s["cr2_peak"]:.1f} dB'
    ax.text(0.03, 0.04,
            f"CR1: {s['cr1_peak']:.1f} dB  |  veg loss {s['cr1_veg_loss']:.1f} dB\n"
            f"CR2: {cr2_label}\n"
            f"Own clutter floor: {s['own_floor']:.1f} dB",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

plt.tight_layout()
save('fig1_reeds_amplitude_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — ROC curves (CR1 only — CR2 April undetectable)
# ══════════════════════════════════════════════════════════════════════════
fig2, ax3 = plt.subplots(figsize=(7, 6))
ax3.plot(mar['fpr'], mar['tpr1'], 'b^-', markersize=5, lw=1.8,
         label='March — CR1 (partial occlusion)')
ax3.plot(apr['fpr'], apr['tpr1'], 'rs-', markersize=5, lw=1.8,
         label='April — CR1 (heavy occlusion)')
ax3.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4, label='Random classifier')
ax3.plot(mar['fpr_fixed'], mar['tpr1_fixed'], 'b*', markersize=16, zorder=5,
         label=f'March @ fixed threshold  (TPR={mar["tpr1_fixed"]:.2f}, FPR={mar["fpr_fixed"]:.2f})')
ax3.plot(apr['fpr_fixed'], apr['tpr1_fixed'], 'r*', markersize=16, zorder=5,
         label=f'April @ fixed threshold  (TPR={apr["tpr1_fixed"]:.2f}, FPR={apr["fpr_fixed"]:.2f})')
ax3.text(0.52, 0.12,
         f'Note: April CR2 ({apr["cr2_peak"]:.1f} dB)\nfalls below clutter floor\n— excluded from ROC',
         transform=ax3.transAxes, fontsize=8,
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
ax3.set_xlabel('False Positive Rate', fontsize=12)
ax3.set_ylabel('True Positive Rate', fontsize=12)
ax3.set_title(f'ROC Curve — CR1 Detection\n'
              f'Fixed control threshold = {CTRL_CLUTTER_FLOOR:.1f} dB  (★ = operating point)',
              fontsize=12, fontweight='bold')
ax3.legend(fontsize=9); ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
plt.tight_layout()
save('fig2_reeds_roc_curves')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — CR points detected & false positives vs threshold
# ══════════════════════════════════════════════════════════════════════════
fig3, (ax4, ax5) = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle('Reeds Scans — Threshold Sweep\nCR1 and CR2 shown separately',
              fontsize=13, fontweight='bold')

# CR points detected
ax4.plot(mar['thresholds'], mar['cr1_det'], 'b^-',  markersize=5, lw=1.8,
         label='March CR1')
ax4.plot(mar['thresholds'], mar['cr2_det'], 'b^--', markersize=5, lw=1.2,
         alpha=0.6, label='March CR2')
ax4.plot(apr['thresholds'], apr['cr1_det'], 'rs-',  markersize=5, lw=1.8,
         label='April CR1')
ax4.plot(apr['thresholds'], apr['cr2_det'], 'rs--', markersize=5, lw=1.2,
         alpha=0.6, label='April CR2 (heavily occluded)')
ax4.axvline(CTRL_CLUTTER_FLOOR, color='black', linestyle='-', lw=1.5,
            label=f'Control threshold ({CTRL_CLUTTER_FLOOR:.1f} dB)')
ax4.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax4.set_ylabel('CR Points Detected', fontsize=11)
ax4.set_title('CR Points Detected vs Threshold', fontsize=11, fontweight='bold')
ax4.legend(fontsize=8); ax4.grid(True, alpha=0.3)
ax4.set_xlim(-70, 70)

# False positives
ax5.plot(mar['thresholds'], mar['fp'], 'b^-', markersize=5, lw=1.8, label='March')
ax5.plot(apr['thresholds'], apr['fp'], 'rs-', markersize=5, lw=1.8, label='April')
ax5.axvline(CTRL_CLUTTER_FLOOR, color='black', linestyle='-', lw=1.5,
            label=f'Control threshold ({CTRL_CLUTTER_FLOOR:.1f} dB)')
ax5.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax5.set_ylabel('False Positive Points', fontsize=11)
ax5.set_title('False Positives vs Threshold\n(Behind-Foliage Region)',
              fontsize=11, fontweight='bold')
ax5.legend(fontsize=8); ax5.grid(True, alpha=0.3)
ax5.set_xlim(-70, 70)

plt.tight_layout()
save('fig3_reeds_threshold_sweep')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — SNR bar chart (CR1 and CR2, all four reflectors)
# ══════════════════════════════════════════════════════════════════════════
fig4, ax6 = plt.subplots(figsize=(8, 5))
labels   = ['March\nCR1', 'March\nCR2', 'April\nCR1', 'April\nCR2']
snr_vals = [mar['cr1_snr_ctrl'], mar['cr2_snr_ctrl'],
            apr['cr1_snr_ctrl'], apr['cr2_snr_ctrl']]
colors   = ['steelblue', 'cornflowerblue', 'tomato', 'lightcoral']
bars     = ax6.bar(labels, snr_vals, color=colors, edgecolor='black', linewidth=0.8,
                   width=0.5)
ax6.axhline(0, color='black', lw=1.2)
for bar, val in zip(bars, snr_vals):
    ypos = val + 1.5 if val >= 0 else val - 5
    ax6.text(bar.get_x() + bar.get_width() / 2, ypos,
             f'{val:.1f} dB', ha='center', va='bottom', fontsize=9, fontweight='bold')
# Label undetected
undetected_idx = [i for i, v in enumerate(snr_vals) if v < 0]
for i in undetected_idx:
    ax6.text(i, snr_vals[i] - 9, 'UNDETECTED', ha='center',
             fontsize=8, color='darkred', fontweight='bold')
ax6.set_ylabel('SNR above Control Clutter Floor (dB)', fontsize=11)
ax6.set_title('SNR vs Control Clutter Floor\nCR1 and CR2 compared across March and April scans',
              fontsize=11, fontweight='bold')
ax6.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
save('fig4_reeds_snr_comparison')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 5 — Standalone detection summary table
# ══════════════════════════════════════════════════════════════════════════
fig5, ax_t = plt.subplots(figsize=(13, 5))
ax_t.axis('off')

col_labels = ['Metric', 'Control', 'March CR1', 'March CR2', 'April CR1', 'April CR2']
rows = [
    ['Range (m)',
     f'{CR1_CTRL}',
     f'{mar["cr1_range"]}',   f'{mar["cr2_range"]}',
     f'{apr["cr1_range"]}',   f'{apr["cr2_range"]}'],

    ['CR peak (dB)',
     f'{CTRL_CR1_PEAK:.1f}',
     f'{mar["cr1_peak"]:.1f}', fmt(mar["cr2_peak"]),
     f'{apr["cr1_peak"]:.1f}', f'{apr["cr2_peak"]:.1f}'],

    ['Expected free-space (dB)',
     '—',
     f'{mar["cr1_expected"]:.1f}', f'{mar["cr2_expected"]:.1f}',
     f'{apr["cr1_expected"]:.1f}', f'{apr["cr2_expected"]:.1f}'],

    ['Vegetation loss (dB)',
     '0',
     f'{mar["cr1_veg_loss"]:.1f}', fmt(mar["cr2_veg_loss"]),
     f'{apr["cr1_veg_loss"]:.1f}', '— (undetected)'],

    ['SNR vs ctrl floor (dB)',
     f'{CTRL_CR1_SNR:.1f}',
     f'{mar["cr1_snr_ctrl"]:.1f}', fmt(mar["cr2_snr_ctrl"]),
     f'{apr["cr1_snr_ctrl"]:.1f}', f'{apr["cr2_snr_ctrl"]:.1f}'],

    ['Own clutter floor (dB)',
     f'{CTRL_CLUTTER_FLOOR:.1f}',
     f'{mar["own_floor"]:.1f}', f'{mar["own_floor"]:.1f}',
     f'{apr["own_floor"]:.1f}', f'{apr["own_floor"]:.1f}'],

    ['TPR @ fixed threshold',
     '—',
     f'{mar["tpr1_fixed"]:.2f}', f'{mar["tpr2_fixed"]:.2f}',
     f'{apr["tpr1_fixed"]:.2f}', f'{apr["tpr2_fixed"]:.2f}'],

    ['FPR @ fixed threshold',
     '—',
     f'{mar["fpr_fixed"]:.2f}', f'{mar["fpr_fixed"]:.2f}',
     f'{apr["fpr_fixed"]:.2f}', f'{apr["fpr_fixed"]:.2f}'],

    ['Detected?',
     'Yes',
     'Yes', 'Yes',
     'Yes', 'No (below clutter floor)'],

    ['Occlusion level',
     'None',
     'Partial', 'Partial',
     'Heavy',   'Heavy'],
]

tbl = ax_t.table(cellText=rows, colLabels=col_labels, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.15, 1.9)

header_color = '#2c5f8a'
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor(header_color)
    tbl[0, j].set_text_props(color='white', fontweight='bold')

col_bg = {1: '#f0f0f0', 2: '#dce8f5', 3: '#eaf2fb', 4: '#fde8e8', 5: '#fff0f0'}
for i in range(1, len(rows) + 1):
    for j, color in col_bg.items():
        tbl[i, j].set_facecolor(color)

# Highlight undetected April CR2
for i in range(1, len(rows) + 1):
    tbl[i, 5].set_facecolor('#ffcccc')
    tbl[i, 5].set_text_props(color='darkred')

ax_t.set_title(
    'Control-Referenced Detection Summary\n'
    'CR1 and CR2 reported separately  —  March and April reeds scans',
    fontsize=12, fontweight='bold', pad=18
)
plt.tight_layout()
save('fig5_reeds_detection_summary_table')
plt.show()

print("\nAll figures saved to:", FIGURES_DIR)
