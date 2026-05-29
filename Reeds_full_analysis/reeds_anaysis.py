import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ── Output folder ─────────────────────────────────────────────────────────
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'reeds_full')
os.makedirs(FIGURES_DIR, exist_ok=True)

def save(name):
    path = os.path.join(FIGURES_DIR, name + '.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

# ── File paths ────────────────────────────────────────────────────────────
MARCH_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260311-114835_.csv'
APRIL_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Reeds_full_analysis\pointcloud_20260408-131953_.csv'

# ── Parameters ────────────────────────────────────────────────────────────
MARCH_REED_FRONT = 4.1
APRIL_REED_FRONT = 3.0
CR_TOL           = 0.3

MARCH_CR1_RANGE  = 10.92
MARCH_CR2_RANGE  = 10.92
APRIL_CR1_RANGE  = 13.88
APRIL_CR2_RANGE  = 14.01

BIN_WIDTH        = 0.5   # range bin width for attenuation profile (metres)


# ── Parse scan ────────────────────────────────────────────────────────────
def parse_scan(filepath, cr1_range, cr2_range, reed_front):
    df  = pd.read_csv(filepath)
    r   = df['range'].values
    amp = df['amplitude'].values

    # Cartesian conversion
    pitch = df['pitch'].values
    roll  = df['roll'].values
    el    = (np.pi / 2) + pitch
    x     = r * np.sin(-pitch) * np.cos(roll)
    y     = r * np.sin(-pitch) * np.sin(roll)
    z     = -r * np.cos(-pitch)

    same_range = abs(cr1_range - cr2_range) < 0.05
    far_cr     = max(cr1_range, cr2_range)
    mask_behind = (r >= reed_front) & (r <= far_cr + CR_TOL)

    if abs(cr1_range - cr2_range) < 2 * CR_TOL and not same_range:
        midpoint = (cr1_range + cr2_range) / 2
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <  midpoint)
        cr2_mask = (r >= midpoint)            & (r <= cr2_range + CR_TOL)
    else:
        cr1_mask = (r >= cr1_range - CR_TOL) & (r <= cr1_range + CR_TOL)
        cr2_mask = (r >= cr2_range - CR_TOL) & (r <= cr2_range + CR_TOL)

    both_cr_mask  = cr1_mask | cr2_mask
    clutter_mask  = mask_behind & ~both_cr_mask
    own_floor     = np.mean(amp[clutter_mask]) if clutter_mask.sum() > 0 else np.nan
    own_floor_std = np.std(amp[clutter_mask])  if clutter_mask.sum() > 0 else np.nan
    sigma3        = own_floor + 3 * own_floor_std

    cr1_peak = amp[cr1_mask].max() if cr1_mask.sum() > 0 else np.nan
    cr2_peak = amp[cr2_mask].max() if cr2_mask.sum() > 0 else np.nan
    cr1_snr  = cr1_peak - own_floor
    cr2_snr  = cr2_peak - own_floor

    cr1_detected = (not np.isnan(cr1_peak)) and (cr1_peak > sigma3)
    cr2_detected = (not np.isnan(cr2_peak)) and (cr2_peak > sigma3)

    thresholds    = np.arange(-70, 70, 5)
    cr1_det_list, cr2_det_list, fp_list  = [], [], []
    tpr1_list,    tpr2_list,    fpr_list = [], [], []
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
        tpr1_list.append(d1 / total_cr1     if total_cr1     > 0 else 0)
        tpr2_list.append(d2 / total_cr2     if total_cr2     > 0 else 0)
        fpr_list.append(fp  / total_clutter if total_clutter > 0 else 0)

    return dict(
        r=r, amp=amp, x=x, y=y, z=z,
        mask_behind=mask_behind, cr1_mask=cr1_mask, cr2_mask=cr2_mask,
        clutter_mask=clutter_mask, both_cr_mask=both_cr_mask,
        cr1_peak=cr1_peak,    cr2_peak=cr2_peak,
        cr1_snr=cr1_snr,      cr2_snr=cr2_snr,
        own_floor=own_floor,  own_floor_std=own_floor_std,
        sigma3=sigma3,
        cr1_detected=cr1_detected, cr2_detected=cr2_detected,
        thresholds=thresholds,
        tpr1=np.array(tpr1_list), tpr2=np.array(tpr2_list),
        fpr=np.array(fpr_list),
        cr1_det=np.array(cr1_det_list), cr2_det=np.array(cr2_det_list),
        fp=np.array(fp_list),
        total_cr1=total_cr1, total_cr2=total_cr2, total_clutter=total_clutter,
        cr1_range=cr1_range, cr2_range=cr2_range,
        same_range=same_range, reed_front=reed_front,
        far_cr=far_cr,
    )


mar = parse_scan(MARCH_CSV, MARCH_CR1_RANGE, MARCH_CR2_RANGE, MARCH_REED_FRONT)
apr = parse_scan(APRIL_CSV, APRIL_CR1_RANGE, APRIL_CR2_RANGE, APRIL_REED_FRONT)

# ── Print summary ─────────────────────────────────────────────────────────
for label, s in [('MARCH', mar), ('APRIL', apr)]:
    print(f"\n{'='*52}\n{label}\n{'='*52}")
    print(f"Behind-foliage region : {s['reed_front']} m to {s['far_cr'] + CR_TOL:.2f} m")
    print(f"Clutter returns       : {s['total_clutter']}")
    print(f"Own clutter floor     : {s['own_floor']:.1f} dB  (σ = {s['own_floor_std']:.1f} dB)")
    print(f"3σ threshold          : {s['sigma3']:.1f} dB")
    note = '  [same range bin]' if s['same_range'] else ''
    print(f"CR1: peak={s['cr1_peak']:.1f} dB  SNR={s['cr1_snr']:.1f} dB  pts={s['total_cr1']}  detected={s['cr1_detected']}{note}")
    print(f"CR2: peak={s['cr2_peak']:.1f} dB  SNR={s['cr2_snr']:.1f} dB  pts={s['total_cr2']}  detected={s['cr2_detected']}{note}")


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — 3D point clouds side by side, coloured by amplitude
# ══════════════════════════════════════════════════════════════════════════
fig1 = plt.figure(figsize=(16, 7))
fig1.suptitle('Reed Bed Radar Point Clouds — Coloured by Amplitude\n'
              'Reed Scan 1: March (partial occlusion) | Reed Scan 2: April (heavy occlusion)',
              fontsize=13, fontweight='bold')

for i, (s, title, cr1_lbl, cr2_lbl) in enumerate([
    (mar, 'Reed Scan 1 — March\nCRs elevated on stands',
          'CR1/CR2 combined\n(10.92 m)', None),
    (apr, 'Reed Scan 2 — April\nCRs at ground level',
          'CR1 (13.88 m)',  'CR2 (14.01 m)'),
]):
    ax = fig1.add_subplot(1, 2, i+1, projection='3d')

    # All background returns
    bg = ~s['both_cr_mask']
    sc = ax.scatter(s['x'][bg], s['y'][bg], s['z'][bg],
                    c=s['amp'][bg], cmap='jet', s=2, alpha=0.4,
                    vmin=-60, vmax=60)

    # CR1 cluster
    if s['cr1_mask'].sum() > 0:
        ax.scatter(s['x'][s['cr1_mask']], s['y'][s['cr1_mask']], s['z'][s['cr1_mask']],
                   c='red', s=20, alpha=1.0, label=cr1_lbl, zorder=5)

    # CR2 cluster (only April — March shares range bin)
    if not s['same_range'] and s['cr2_mask'].sum() > 0:
        ax.scatter(s['x'][s['cr2_mask']], s['y'][s['cr2_mask']], s['z'][s['cr2_mask']],
                   c='blue', s=20, alpha=1.0, label=cr2_lbl, zorder=5)

    cb = fig1.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
    cb.set_label('Amplitude (dB)', fontsize=9)
    ax.set_xlabel('X (m)', fontsize=9); ax.set_ylabel('Y (m)', fontsize=9)
    ax.set_zlabel('Z (m)', fontsize=9)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')
    ax.view_init(elev=25, azim=45)

plt.tight_layout()
save('fig1_3d_pointclouds')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — Amplitude vs Range (March and April side by side)
# ══════════════════════════════════════════════════════════════════════════
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
fig2.suptitle('Amplitude vs Range — Scan-Specific Clutter Floor\n'
              'Reed Scan 1: March (partial occlusion) | Reed Scan 2: April (heavy occlusion)',
              fontsize=13, fontweight='bold')

for ax, s, title, dot_color in [
    (ax1, mar, 'Reed Scan 1 — March: CRs elevated',          'steelblue'),
    (ax2, apr, 'Reed Scan 2 — April: CRs at ground level',   'tomato'),
]:
    ax.scatter(s['r'][s['r'] < 20], s['amp'][s['r'] < 20],
               s=1, alpha=0.2, color=dot_color, label='All returns')
    ax.axvspan(s['reed_front'], s['far_cr'] + CR_TOL,
               alpha=0.08, color='orange', label='Behind-foliage region')
    ax.axvspan(s['cr1_range'] - CR_TOL, s['cr1_range'] + CR_TOL,
               alpha=0.3, color='red',
               label=f'CR1 window ({s["cr1_range"]} m)' + (' [+ CR2]' if s['same_range'] else ''))
    if not s['same_range']:
        ax.axvspan(s['cr2_range'] - CR_TOL, s['cr2_range'] + CR_TOL,
                   alpha=0.3, color='blue', label=f'CR2 window ({s["cr2_range"]} m)')
    ax.axvline(s['reed_front'], color='green', linestyle='--', lw=1.3,
               label=f'Reed front ({s["reed_front"]} m)')
    ax.axhline(s['own_floor'], color='black', linestyle='-', lw=1.8,
               label=f'Clutter floor ({s["own_floor"]:.1f} dB)')
    ax.axhline(s['sigma3'], color='dimgrey', linestyle='-.', lw=1.3,
               label=f'3σ threshold ({s["sigma3"]:.1f} dB)')
    ax.axhline(s['cr1_peak'], color='red', linestyle=':', lw=1.5,
               label=f'CR1 peak ({s["cr1_peak"]:.1f} dB)')
    if not s['same_range']:
        ax.axhline(s['cr2_peak'], color='blue', linestyle=':', lw=1.5,
                   label=f'CR2 peak ({s["cr2_peak"]:.1f} dB)')
    ax.set_xlim(0, 20); ax.set_ylim(-100, 90)
    ax.set_xlabel('Range (m)', fontsize=11)
    ax.set_ylabel('Amplitude (dB)', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=7.5, loc='upper right')
    ax.grid(True, alpha=0.3)
    if s['same_range']:
        info = (f"CR (combined): {s['cr1_peak']:.1f} dB  |  SNR {s['cr1_snr']:.1f} dB\n"
                f"Floor: {s['own_floor']:.1f} dB (σ={s['own_floor_std']:.1f})  |  3σ: {s['sigma3']:.1f} dB\n"
                f"Detected (3σ): {'Yes' if s['cr1_detected'] else 'No'}")
    else:
        info = (f"CR1: {s['cr1_peak']:.1f} dB  SNR {s['cr1_snr']:.1f} dB  {'✓' if s['cr1_detected'] else '✗'}\n"
                f"CR2: {s['cr2_peak']:.1f} dB  SNR {s['cr2_snr']:.1f} dB  {'✓' if s['cr2_detected'] else '✗'}\n"
                f"Floor: {s['own_floor']:.1f} dB (σ={s['own_floor_std']:.1f})  |  3σ: {s['sigma3']:.1f} dB")
    ax.text(0.03, 0.04, info, transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))

plt.tight_layout()
save('fig2_amplitude_vs_range')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — Amplitude histograms (clutter distribution vs CR peaks)
# ══════════════════════════════════════════════════════════════════════════
fig3, axes = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle('Amplitude Distribution — Behind-Foliage Clutter vs Corner Reflector Returns\n'
              'Histogram shows clutter-only returns; CR peaks shown as vertical lines',
              fontsize=12, fontweight='bold')

for ax, s, title, hist_color in [
    (axes[0], mar, 'Reed Scan 1 — March', [0.2, 0.5, 0.8]),
    (axes[1], apr, 'Reed Scan 2 — April', [0.85, 0.25, 0.15]),
]:
    clutter_amp = s['amp'][s['clutter_mask']]
    ax.hist(clutter_amp, bins=40, color=hist_color, alpha=0.7,
            label=f'Behind-foliage clutter\n(n={len(clutter_amp)} returns)', edgecolor='none')
    ax.axvline(s['own_floor'], color='black', lw=2, linestyle='-',
               label=f'Mean floor ({s["own_floor"]:.1f} dB)')
    ax.axvline(s['sigma3'], color='dimgrey', lw=1.8, linestyle='-.',
               label=f'3σ threshold ({s["sigma3"]:.1f} dB)')
    ax.axvline(s['cr1_peak'], color='red', lw=2, linestyle='--',
               label=f'CR1 peak ({s["cr1_peak"]:.1f} dB)')
    if not s['same_range']:
        ax.axvline(s['cr2_peak'], color='blue', lw=2, linestyle='--',
                   label=f'CR2 peak ({s["cr2_peak"]:.1f} dB)')
    ax.set_xlabel('Amplitude (dB)', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-100, 80)

plt.tight_layout()
save('fig3_amplitude_histograms')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — Range-binned attenuation profile through the reeds
#  Shows mean amplitude in each range bin — how signal decays with depth
# ══════════════════════════════════════════════════════════════════════════
fig4, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
fig4.suptitle('Range-Binned Amplitude Profile Through Reed Canopy\n'
              'Mean ± 1σ amplitude per range bin — clutter returns only (CR windows excluded)',
              fontsize=12, fontweight='bold')

for ax, s, title, color in [
    (ax_l, mar, 'Reed Scan 1 — March', 'steelblue'),
    (ax_r, apr, 'Reed Scan 2 — April', 'tomato'),
]:
    # Use all returns (not just behind-foliage) to show full scene profile
    # but exclude CR windows so they don't skew the attenuation curve
    r_plot  = s['r'][~s['both_cr_mask']]
    a_plot  = s['amp'][~s['both_cr_mask']]

    r_max   = s['far_cr'] + 2.0
    bins    = np.arange(0, r_max + BIN_WIDTH, BIN_WIDTH)
    bin_centres = bins[:-1] + BIN_WIDTH / 2
    bin_means, bin_stds, bin_counts = [], [], []

    for b_lo, b_hi in zip(bins[:-1], bins[1:]):
        mask = (r_plot >= b_lo) & (r_plot < b_hi)
        if mask.sum() >= 3:
            bin_means.append(np.mean(a_plot[mask]))
            bin_stds.append(np.std(a_plot[mask]))
            bin_counts.append(mask.sum())
        else:
            bin_means.append(np.nan)
            bin_stds.append(np.nan)
            bin_counts.append(0)

    bin_means  = np.array(bin_means)
    bin_stds   = np.array(bin_stds)
    valid      = ~np.isnan(bin_means)

    ax.plot(bin_centres[valid], bin_means[valid], '-o', color=color,
            lw=2, markersize=5, label='Mean amplitude per bin')
    ax.fill_between(bin_centres[valid],
                    bin_means[valid] - bin_stds[valid],
                    bin_means[valid] + bin_stds[valid],
                    alpha=0.25, color=color, label='±1σ')

    # Reference lines
    ax.axvline(s['reed_front'], color='green', linestyle='--', lw=1.5,
               label=f'Reed front face ({s["reed_front"]} m)')
    ax.axhline(s['own_floor'], color='black', linestyle='-', lw=1.5,
               label=f'Clutter floor ({s["own_floor"]:.1f} dB)')
    ax.axhline(s['sigma3'], color='dimgrey', linestyle='-.', lw=1.2,
               label=f'3σ threshold ({s["sigma3"]:.1f} dB)')

    # Mark CR positions
    ax.axvline(s['cr1_range'], color='red', linestyle=':', lw=1.5,
               label=f'CR1 position ({s["cr1_range"]} m)')
    if not s['same_range']:
        ax.axvline(s['cr2_range'], color='blue', linestyle=':', lw=1.5,
                   label=f'CR2 position ({s["cr2_range"]} m)')

    # Shade reed region
    ax.axvspan(s['reed_front'], s['far_cr'] + CR_TOL,
               alpha=0.06, color='orange', label='Behind-foliage region')

    ax.set_xlabel('Range (m)', fontsize=11)
    ax.set_ylabel('Mean Amplitude (dB)', fontsize=11)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, r_max)
    ax.set_ylim(-80, 70)

plt.tight_layout()
save('fig4_range_binned_attenuation')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 5 — ROC curves (CR1 and CR2 separately, both scans)
# ══════════════════════════════════════════════════════════════════════════
fig5, ax5 = plt.subplots(figsize=(7, 6))

ax5.plot(mar['fpr'], mar['tpr1'], 'b^-',  markersize=5, lw=1.8,
         label='March CR1/CR2 (combined, partial occlusion)')
ax5.plot(apr['fpr'], apr['tpr1'], 'rs-',  markersize=5, lw=1.8,
         label='April CR1 (heavy occlusion)')
ax5.plot(apr['fpr'], apr['tpr2'], 'rs--', markersize=5, lw=1.2, alpha=0.65,
         label='April CR2 (heavy occlusion)')
ax5.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4, label='Random classifier')

# Operating points at clutter floor (★) and 3σ (◇)
for s, color, lbl, tpr_key in [
    (mar, 'blue',    'March',       'tpr1'),
    (apr, 'red',     'April CR1',   'tpr1'),
    (apr, 'darkred', 'April CR2',   'tpr2'),
]:
    idx_floor = np.argmin(np.abs(s['thresholds'] - s['own_floor']))
    idx_sig3  = np.argmin(np.abs(s['thresholds'] - s['sigma3']))
    ax5.plot(s['fpr'][idx_floor], s[tpr_key][idx_floor],
             marker='*', markersize=14, color=color, zorder=5,
             label=f'{lbl} @ floor (TPR={s[tpr_key][idx_floor]:.2f}, FPR={s["fpr"][idx_floor]:.2f})')
    ax5.plot(s['fpr'][idx_sig3], s[tpr_key][idx_sig3],
             marker='D', markersize=8, color=color, zorder=5, fillstyle='none',
             label=f'{lbl} @ 3σ (TPR={s[tpr_key][idx_sig3]:.2f}, FPR={s["fpr"][idx_sig3]:.2f})')

ax5.set_xlabel('False Positive Rate', fontsize=12)
ax5.set_ylabel('True Positive Rate', fontsize=12)
ax5.set_title('ROC Curves — Scan-Specific Thresholds\n'
              '★ = clutter floor operating point  ◇ = 3σ operating point',
              fontsize=11, fontweight='bold')
ax5.legend(fontsize=7.5, loc='lower right')
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, 1); ax5.set_ylim(0, 1)
plt.tight_layout()
save('fig5_roc_curves')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 6 — Threshold sweep (CR points detected and false positives)
# ══════════════════════════════════════════════════════════════════════════
fig6, (ax6l, ax6r) = plt.subplots(1, 2, figsize=(13, 5))
fig6.suptitle('Threshold Sweep — Behind-Foliage Region\n'
              'CR points detected and false positives as a function of detection threshold',
              fontsize=12, fontweight='bold')

ax6l.plot(mar['thresholds'], mar['cr1_det'], 'b^-',  markersize=5, lw=1.8,
          label='March CR1/CR2 (combined)')
ax6l.plot(apr['thresholds'], apr['cr1_det'], 'rs-',  markersize=5, lw=1.8,
          label='April CR1')
ax6l.plot(apr['thresholds'], apr['cr2_det'], 'rs--', markersize=5, lw=1.2, alpha=0.65,
          label='April CR2')
for s, color in [(mar, 'blue'), (apr, 'red')]:
    lbl = 'March' if color == 'blue' else 'April'
    ax6l.axvline(s['own_floor'], color=color, linestyle='--', lw=1.2, alpha=0.7,
                 label=f'{lbl} floor ({s["own_floor"]:.1f} dB)')
    ax6l.axvline(s['sigma3'],    color=color, linestyle=':',  lw=1.2, alpha=0.7,
                 label=f'{lbl} 3σ ({s["sigma3"]:.1f} dB)')
ax6l.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax6l.set_ylabel('CR Points Detected', fontsize=11)
ax6l.set_title('CR Points Detected vs Threshold', fontsize=11, fontweight='bold')
ax6l.legend(fontsize=8); ax6l.grid(True, alpha=0.3); ax6l.set_xlim(-70, 70)

ax6r.plot(mar['thresholds'], mar['fp'], 'b^-', markersize=5, lw=1.8, label='March')
ax6r.plot(apr['thresholds'], apr['fp'], 'rs-', markersize=5, lw=1.8, label='April')
for s, color in [(mar, 'blue'), (apr, 'red')]:
    lbl = 'March' if color == 'blue' else 'April'
    ax6r.axvline(s['own_floor'], color=color, linestyle='--', lw=1.2, alpha=0.7,
                 label=f'{lbl} floor ({s["own_floor"]:.1f} dB)')
    ax6r.axvline(s['sigma3'],    color=color, linestyle=':',  lw=1.2, alpha=0.7,
                 label=f'{lbl} 3σ ({s["sigma3"]:.1f} dB)')
ax6r.set_xlabel('Detection Threshold (dB)', fontsize=11)
ax6r.set_ylabel('False Positive Points', fontsize=11)
ax6r.set_title('False Positives vs Threshold\n(Behind-Foliage Region)',
               fontsize=11, fontweight='bold')
ax6r.legend(fontsize=8); ax6r.grid(True, alpha=0.3); ax6r.set_xlim(-70, 70)

plt.tight_layout()
save('fig6_threshold_sweep')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 7 — SNR bar chart comparison
# ══════════════════════════════════════════════════════════════════════════
fig7, ax7 = plt.subplots(figsize=(8, 5))

labels   = ['March\nCR1/CR2\n(combined)', 'April\nCR1', 'April\nCR2']
snr_vals = [mar['cr1_snr'], apr['cr1_snr'], apr['cr2_snr']]
colors   = ['steelblue', 'tomato', 'salmon']
detected = [mar['cr1_detected'], apr['cr1_detected'], apr['cr2_detected']]

bars = ax7.bar(labels, snr_vals, color=colors, edgecolor='black', linewidth=0.8, width=0.5)

# 3σ threshold lines — one per scan
ax7.axhline(mar['cr1_snr'] - mar['cr1_snr'] + (mar['sigma3'] - mar['own_floor']),
            color='steelblue', linestyle='--', lw=1.5, alpha=0.7,
            label=f'March 3σ margin ({mar["sigma3"] - mar["own_floor"]:.1f} dB above floor)')
ax7.axhline(apr['sigma3'] - apr['own_floor'],
            color='tomato', linestyle='--', lw=1.5, alpha=0.7,
            label=f'April 3σ margin ({apr["sigma3"] - apr["own_floor"]:.1f} dB above floor)')

# Annotate detected / not detected
for bar, det, snr in zip(bars, detected, snr_vals):
    label = '✓ Detected' if det else '✗ Not detected\n(below 3σ)'
    color = 'darkgreen' if det else 'darkred'
    ax7.text(bar.get_x() + bar.get_width()/2, snr + 1.5,
             label, ha='center', va='bottom', fontsize=9,
             fontweight='bold', color=color)

ax7.set_ylabel('SNR above clutter floor (dB)', fontsize=11)
ax7.set_title('Corner Reflector SNR — March vs April Reed Scans\n'
              'Dashed lines show 3σ detection threshold (SNR required for robust detection)',
              fontsize=11, fontweight='bold')
ax7.legend(fontsize=9)
ax7.grid(True, alpha=0.3, axis='y')
ax7.set_ylim(0, max(snr_vals) * 1.25)
plt.tight_layout()
save('fig7_snr_bar_chart')


# ══════════════════════════════════════════════════════════════════════════
#  FIGURE 8 — Summary table
# ══════════════════════════════════════════════════════════════════════════
fig8, ax8 = plt.subplots(figsize=(14, 5))
ax8.axis('off')

col_labels = ['Metric', 'March CR1', 'March CR2', 'April CR1', 'April CR2']

def fmt(v):
    return f'{v:.1f}' if not np.isnan(v) else 'N/A'

rows = [
    ['Range (m)',
     f'{mar["cr1_range"]}',          f'{mar["cr2_range"]} [shared bin]',
     f'{apr["cr1_range"]}',          f'{apr["cr2_range"]}'],
    ['Occlusion level',
     'Partial', 'Partial', 'Heavy', 'Heavy'],
    ['CR peak amplitude (dB)',
     fmt(mar['cr1_peak']),            fmt(mar['cr2_peak']),
     fmt(apr['cr1_peak']),            fmt(apr['cr2_peak'])],
    ['Clutter floor (dB)',
     f'{mar["own_floor"]:.1f}',       f'{mar["own_floor"]:.1f}',
     f'{apr["own_floor"]:.1f}',       f'{apr["own_floor"]:.1f}'],
    ['Clutter σ (dB)',
     f'{mar["own_floor_std"]:.1f}',   f'{mar["own_floor_std"]:.1f}',
     f'{apr["own_floor_std"]:.1f}',   f'{apr["own_floor_std"]:.1f}'],
    ['3σ threshold (dB)',
     f'{mar["sigma3"]:.1f}',          f'{mar["sigma3"]:.1f}',
     f'{apr["sigma3"]:.1f}',          f'{apr["sigma3"]:.1f}'],
    ['SNR above floor (dB)',
     fmt(mar['cr1_snr']),             fmt(mar['cr2_snr']),
     fmt(apr['cr1_snr']),             fmt(apr['cr2_snr'])],
    ['CR points in window',
     f'{mar["total_cr1"]}',           f'{mar["total_cr2"]}',
     f'{apr["total_cr1"]}',           f'{apr["total_cr2"]}'],
    ['Detected (3σ)?',
     'No — high clutter σ',           'No — same bin as CR1',
     'Yes',                           'No — SNR < 3σ margin'],
]

tbl = ax8.table(cellText=rows, colLabels=col_labels, loc='center', cellLoc='center')
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1.1, 2.0)

header_color = '#2c5f8a'
for j in range(len(col_labels)):
    tbl[0, j].set_facecolor(header_color)
    tbl[0, j].set_text_props(color='white', fontweight='bold')

col_bg = {1: '#dce8f5', 2: '#eaf2fb', 3: '#fde8e8', 4: '#ffcccc'}
for i in range(1, len(rows) + 1):
    for j, color in col_bg.items():
        tbl[i, j].set_facecolor(color)

# Highlight undetected cells
undetected_cols = [1, 2, 4]  # March CR1, March CR2, April CR2
for i in range(1, len(rows) + 1):
    for j in undetected_cols:
        tbl[i, j].set_text_props(color='darkred')
# April CR1 in green
for i in range(1, len(rows) + 1):
    tbl[i, 3].set_text_props(color='darkgreen')

ax8.set_title('Scan-Specific Detection Summary — Reed Scan 1 (March) vs Reed Scan 2 (April)\n'
              'Clutter floor estimated from behind-foliage region only  |  Detection criterion: peak > 3σ threshold',
              fontsize=11, fontweight='bold', pad=20)
plt.tight_layout()
save('fig8_summary_table')

print(f"\nAll 8 figures saved to: {FIGURES_DIR}")