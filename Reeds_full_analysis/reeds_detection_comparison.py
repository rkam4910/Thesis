import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Output folder — saves to figures/scanspecific/ next to this script ────
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'scanspecific')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"Saved: {path}")

# ── File paths ────────────────────────────────────────────────────────────
MARCH_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260311-114835_.csv'
APRIL_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260408-131953_.csv'

# ── Parameters ────────────────────────────────────────────────────────────
MARCH_REED_FRONT = 4.1
APRIL_REED_FRONT = 3.0
CR_TOL           = 0.3      # tight tolerance so April CR1/CR2 don't bleed

MARCH_CR1_RANGE  = 10.92
MARCH_CR2_RANGE  = 10.92   # same range bin in March
APRIL_CR1_RANGE  = 13.88
APRIL_CR2_RANGE  = 14.01


# ── Parse scan — separate CR1 and CR2 ────────────────────────────────────
def parse_scan(filepath, cr1_range, cr2_range, reed_front):
    df  = pd.read_csv(filepath)
    r   = df['range'].values
    amp = df['amplitude'].values

    same_range  = abs(cr1_range - cr2_range) < 0.05
    far_cr      = max(cr1_range, cr2_range)
    mask_behind = (r >= reed_front) & (r <= far_cr + CR_TOL)

    # If CR1 and CR2 are so close their ±CR_TOL windows would overlap,
    # split at the midpoint so each gets an exclusive window
    if abs(cr1_range - cr2_range) < 2 * CR_TOL and not same_range:
        midpoint = (cr1_range + cr2_range) / 2
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <  midpoint)
        cr2_mask = (r >= midpoint)            & (r <= cr2_range + CR_TOL)
    else:
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <= cr1_range + CR_TOL)
        cr2_mask = (r >= cr2_range - CR_TOL) & (r <= cr2_range + CR_TOL)

    both_cr_mask = cr1_mask | cr2_mask
    clutter_mask = mask_behind & ~both_cr_mask

    cr1_peak      = amp[cr1_mask].max() if cr1_mask.sum() > 0 else np.nan
    cr2_peak      = amp[cr2_mask].max() if cr2_mask.sum() > 0 else np.nan
    own_floor     = np.mean(amp[clutter_mask]) if clutter_mask.sum() > 0 else np.nan
    own_floor_std = np.std(amp[clutter_mask])  if clutter_mask.sum() > 0 else np.nan

    cr1_snr = cr1_peak - own_floor
    cr2_snr = cr2_peak - own_floor

    thresholds    = np.arange(-70, 70, 5)
    cr1_det_list, cr2_det_list, fp_list           = [], [], []
    tpr1_list,    tpr2_list,    fpr_list          = [], [], []
    total_cr1     = cr1_mask.sum()
    total_cr2     = cr2_mask.sum()
    total_clutter = clutter_mask.sum()

    for t in thresholds:
        d1 = np.sum(amp[cr1_mask]     >= t)
        d2 = np.sum(amp[cr2_mask]     >= t)
        fp = np.sum(amp[clutter_mask] >= t)
        cr1_det_list.append(d1)
        cr2_det_list.append(d2)
        fp_list.append(fp)
        tpr1_list.append(d1 / total_cr1      if total_cr1      > 0 else 0)
        tpr2_list.append(d2 / total_cr2      if total_cr2      > 0 else 0)
        fpr_list.append(fp  / total_clutter  if total_clutter  > 0 else 0)

    return dict(
        r=r, amp=amp,
        mask_behind=mask_behind, cr1_mask=cr1_mask, cr2_mask=cr2_mask,
        clutter_mask=clutter_mask,
        cr1_peak=cr1_peak,   cr2_peak=cr2_peak,
        cr1_snr=cr1_snr,     cr2_snr=cr2_snr,
        own_floor=own_floor, own_floor_std=own_floor_std,
        thresholds=thresholds,
        tpr1=np.array(tpr1_list), tpr2=np.array(tpr2_list),
        fpr=np.array(fpr_list),
        cr1_det=np.array(cr1_det_list), cr2_det=np.array(cr2_det_list),
        fp=np.array(fp_list),
        total_cr1=total_cr1, total_cr2=total_cr2, total_clutter=total_clutter,
        cr1_detected=(cr1_peak > own_floor),
        cr2_detected=(not np.isnan(cr2_peak) and cr2_peak > own_floor),
        cr1_range=cr1_range, cr2_range=cr2_range,
        same_range=same_range, reed_front=reed_front,
    )


mar = parse_scan(MARCH_CSV, MARCH_CR1_RANGE, MARCH_CR2_RANGE, MARCH_REED_FRONT)
apr = parse_scan(APRIL_CSV, APRIL_CR1_RANGE, APRIL_CR2_RANGE, APRIL_REED_FRONT)

for label, s in [('MARCH', mar), ('APRIL', apr)]:
    print(f"\n{'='*52}\n{label}\n{'='*52}")
    print(f"Own clutter floor:  {s['own_floor']:.1f} dB  (σ = {s['own_floor_std']:.1f} dB)")
    print(f"CR1: peak={s['cr1_peak']:.1f} dB  SNR={s['cr1_snr']:.1f} dB  pts={s['total_cr1']}  detected={s['cr1_detected']}")
    print(f"CR2: peak={s['cr2_peak']:.1f} dB  SNR={s['cr2_snr']:.1f} dB  pts={s['total_cr2']}  detected={s['cr2_detected']}")
print(f"CR2: peak={s['cr2_peak']:.1f} dB  SNR={s['cr2_snr']:.1f} dB  pts={s['total_cr2']}  detected={s['cr2_detected']}")

# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Amplitude vs Range (March and April side by side)
# ══════════════════════════════════════════════════════════════════════════
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
fig1.suptitle('Scan-Specific Threshold — Amplitude vs Range\n'
              'March (partial occlusion) vs April (heavy occlusion)',
              fontsize=13, fontweight='bold')

for ax, s, title, dot_color in [
    (ax1, mar, 'March — CRs elevated, partial occlusion',       'steelblue'),
    (ax2, apr, 'April — CRs at ground level, heavy occlusion',  'tomato'),
]:
    ax.scatter(s['r'][s['r'] < 20], s['amp'][s['r'] < 20],
               s=1, alpha=0.2, color=dot_color, label='All returns')
    ax.axvspan(s['reed_front'], max(s['cr1_range'], s['cr2_range']) + CR_TOL,
               alpha=0.08, color='orange', label='Behind-foliage region')
    # CR1 band
    ax.axvspan(s['cr1_range'] - CR_TOL, s['cr1_range'] + CR_TOL,
               alpha=0.3, color='red', label=f'CR1 ({s["cr1_range"]} m)')
    # CR2 band — always show; label notes overlap if same range
    cr2_lbl = (f'CR2 ({s["cr2_range"]} m) — same range as CR1'
               if s['same_range'] else f'CR2 ({s["cr2_range"]} m)')
    ax.axvspan(s['cr2_range'] - CR_TOL, s['cr2_range'] + CR_TOL,
               alpha=0.3, color='blue', label=cr2_lbl)
    ax.axvline(s['reed_front'], color='green', linestyle='--', lw=1.3,
               label=f'Reed front ({s["reed_front"]} m)')
    ax.axhline(s['own_floor'], color='black', linestyle='-', lw=1.5,
               label=f'Clutter floor ({s["own_floor"]:.1f} dB)')
    ax.axhline(s['cr1_peak'], color='red', linestyle=':', lw=1.2,
               label=f'CR1 peak ({s["cr1_peak"]:.1f} dB)')
    cr2_peak_lbl = (f'CR2 peak ({s["cr2_peak"]:.1f} dB)  UNDETECTED'
                    if not s['cr2_detected'] else f'CR2 peak ({s["cr2_peak"]:.1f} dB)')
    ax.axhline(s['cr2_peak'], color='blue', linestyle=':', lw=1.2,
               label=cr2_peak_lbl)
    ax.set_xlim(0, 20); ax.set_ylim(-100, 90)
    ax.set_xlabel('Range (m)', fontsize=11)
    ax.set_ylabel('Amplitude (dB)', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, alpha=0.3)
    cr2_txt = (f'{s["cr2_peak"]:.1f} dB  UNDETECTED'
               if not s['cr2_detected'] else f'{s["cr2_peak"]:.1f} dB')
    ax.text(0.03, 0.04,
            f"CR1: {s['cr1_peak']:.1f} dB  |  SNR {s['cr1_snr']:.1f} dB\n"
            f"CR2: {cr2_txt}\n"
            f"Own clutter floor: {s['own_floor']:.1f} dB",
            transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

plt.tight_layout()
save('fig1_scanspec_amplitude_vs_range')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — ROC curves (CR1 and CR2 separately)
# ══════════════════════════════════════════════════════════════════════════
fig2, ax3 = plt.subplots(figsize=(7, 6))
ax3.plot(mar['fpr'], mar['tpr1'], 'b^-',  markersize=5, lw=1.8, label='March CR1')
ax3.plot(mar['fpr'], mar['tpr2'], 'b^--', markersize=5, lw=1.2, alpha=0.6, label='March CR2')
ax3.plot(apr['fpr'], apr['tpr1'], 'rs-',  markersize=5, lw=1.8, label='April CR1')
ax3.plot(apr['fpr'], apr['tpr2'], 'rs--', markersize=5, lw=1.2, alpha=0.6,
         label='April CR2 (heavily occluded)')
ax3.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4, label='Random classifier')

# Mark operating point at each scan's own clutter floor for CR1
for s, color, lbl in [(mar, 'blue', 'March'), (apr, 'red', 'April')]:
    idx = np.argmin(np.abs(s['thresholds'] - s['own_floor']))
    ax3.plot(s['fpr'][idx], s['tpr1'][idx], marker='*', markersize=14,
             color=color, zorder=5,
             label=f'{lbl} CR1 @ own floor  (TPR={s["tpr1"][idx]:.2f}, FPR={s["fpr"][idx]:.2f})')

if not apr['cr2_detected']:
    ax3.text(0.52, 0.08,
             f'April CR2 ({apr["cr2_peak"]:.1f} dB)\nfalls below clutter floor\n— excluded from ROC',
             transform=ax3.transAxes, fontsize=8,
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

ax3.set_xlabel('False Positive Rate', fontsize=12)
ax3.set_ylabel('True Positive Rate', fontsize=12)
ax3.set_title('ROC Curve — Scan-Specific Thresholds\nCR1 and CR2 shown separately\n'
              '(★ = operating point at own clutter floor)',
              fontsize=11, fontweight='bold')
ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
plt.tight_layout()
save('fig2_scanspec_roc_curves')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — Threshold sweep (CR1 and CR2 separately)
# ══════════════════════════════════════════════════════════════════════════
fig3, (ax4, ax5) = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle('Scan-Specific Threshold Sweep — CR1 and CR2 shown separately',
              fontsize=13, fontweight='bold')

# CR points detected
ax4.plot(mar['thresholds'], mar['cr1_det'], 'b^-',  markersize=5, lw=1.8, label='March CR1')
ax4.plot(mar['thresholds'], mar['cr2_det'], 'b^--', markersize=5, lw=1.2, alpha=0.6, label='March CR2')
ax4.plot(apr['thresholds'], apr['cr1_det'], 'rs-',  markersize=5, lw=1.8, label='April CR1')
ax4.plot(apr['thresholds'], apr['cr2_det'], 'rs--', markersize=5, lw=1.2, alpha=0.6,
         label='April CR2 (heavily occluded)')
ax4.axvline(mar['own_floor'], color='blue', linestyle='--', lw=1, alpha=0.6,
            label=f"March floor ({mar['own_floor']:.1f} dB)")
ax4.axvline(apr['own_floor'], color='red',  linestyle='--', lw=1, alpha=0.6,
            label=f"April floor ({apr['own_floor']:.1f} dB)")
ax4.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax4.set_ylabel('CR Points Detected', fontsize=11)
ax4.set_title('CR Points Detected vs Threshold', fontsize=11, fontweight='bold')
ax4.legend(fontsize=8); ax4.grid(True, alpha=0.3)
ax4.set_xlim(-70, 70)

# False positives
ax5.plot(mar['thresholds'], mar['fp'], 'b^-', markersize=5, lw=1.8, label='March')
ax5.plot(apr['thresholds'], apr['fp'], 'rs-', markersize=5, lw=1.8, label='April')
ax5.axvline(mar['own_floor'], color='blue', linestyle='--', lw=1, alpha=0.6,
            label=f"March floor ({mar['own_floor']:.1f} dB)")
ax5.axvline(apr['own_floor'], color='red',  linestyle='--', lw=1, alpha=0.6,
            label=f"April floor ({apr['own_floor']:.1f} dB)")
ax5.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax5.set_ylabel('False Positive Points', fontsize=11)
ax5.set_title('False Positives vs Threshold\n(Behind-Foliage Region)',
              fontsize=11, fontweight='bold')
ax5.legend(fontsize=8); ax5.grid(True, alpha=0.3)
ax5.set_xlim(-70, 70)

plt.tight_layout()
save('fig3_scanspec_threshold_sweep')
plt.show()


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — Summary table
# ══════════════════════════════════════════════════════════════════════════
fig4, ax_t = plt.subplots(figsize=(11, 4))
ax_t.axis('off')

col_labels = ['Metric', 'March CR1', 'March CR2', 'April CR1', 'April CR2']

def fmt(v):
    return f'{v:.1f}' if not np.isnan(v) else 'N/A'

rows = [
    ['Range (m)',
     f'{mar["cr1_range"]}',   f'{mar["cr2_range"]}',
     f'{apr["cr1_range"]}',   f'{apr["cr2_range"]}'],
    ['CR peak (dB)',
     fmt(mar['cr1_peak']),    fmt(mar['cr2_peak']),
     fmt(apr['cr1_peak']),    fmt(apr['cr2_peak'])],
    ['Own clutter floor (dB)',
     f'{mar["own_floor"]:.1f}', f'{mar["own_floor"]:.1f}',
     f'{apr["own_floor"]:.1f}', f'{apr["own_floor"]:.1f}'],
    ['SNR above own floor (dB)',
     fmt(mar['cr1_snr']),    fmt(mar['cr2_snr']),
     fmt(apr['cr1_snr']),    fmt(apr['cr2_snr'])],
    ['CR points detected',
     f'{mar["total_cr1"]}',  f'{mar["total_cr2"]}',
     f'{apr["total_cr1"]}',  f'{apr["total_cr2"]}'],
    ['Occlusion level',  'Partial', 'Partial', 'Heavy', 'Heavy'],
    ['Detected?',        'Yes',     'Yes',     'Yes',   'No (below floor)'],
]

tbl = ax_t.table(cellText=rows, colLabels=col_labels, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.15, 1.9)

header_color = '#2c5f8a'
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor(header_color)
    tbl[0, j].set_text_props(color='white', fontweight='bold')

col_bg = {1: '#dce8f5', 2: '#eaf2fb', 3: '#fde8e8', 4: '#ffcccc'}
for i in range(1, len(rows) + 1):
    for j, color in col_bg.items():
        tbl[i, j].set_facecolor(color)
    tbl[i, 4].set_text_props(color='darkred')

ax_t.set_title('Scan-Specific Detection Summary\nCR1 and CR2 reported separately',
               fontsize=12, fontweight='bold', pad=18)
plt.tight_layout()
save('fig4_scanspec_summary_table')
plt.show()

print("\nAll figures saved to:", FIGURES_DIR)