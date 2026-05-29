import matplotlib.pyplot as plt
import numpy as np

# ── Data ─────────────────────────────────────────────────────────────────────
labels = [
    'Reed March\n(Partial)',
    'Reed April CR1\n(Heavy, multipath)',
    'Reed April CR2\n(Heavy)',
    'Cherry Ballart\nSparse',
    'Cherry Ballart\nDense',
]

cr_peak      = [18.0,  53.0, -11.0,  18.0, -41.0]
clutter_floor = [-31.8, -33.0, -33.0, -55.0, -58.0]
snr          = [49.8,  86.0,  22.0,  73.0,  17.0]
threshold_3s = [34.1,  22.5,  22.5,  -1.6, -22.4]  # 3sigma thresholds
detected     = [True,  True,  False,  True,  False]

x = np.arange(len(labels))
width = 0.30

# ── Colors ───────────────────────────────────────────────────────────────────
col_peak    = '#2E86AB'   # blue  — CR peak
col_clutter = '#A8A8A8'   # grey  — clutter floor
col_snr     = '#E84855'   # red   — SNR bar

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.subplots_adjust(wspace=0.35)

# ── Left plot: CR peak vs clutter floor ──────────────────────────────────────
ax1 = axes[0]

bars_peak    = ax1.bar(x - width/2, cr_peak,       width, label='CR peak (dB)',
                       color=col_peak,    edgecolor='white', linewidth=0.5)
bars_clutter = ax1.bar(x + width/2, clutter_floor, width, label='Clutter floor (dB)',
                       color=col_clutter, edgecolor='white', linewidth=0.5)

# 3-sigma threshold markers
for i, (xi, t3) in enumerate(zip(x, threshold_3s)):
    ax1.hlines(t3, xi - width, xi + width, colors='black',
               linewidths=1.2, linestyles='--', zorder=5)

# Detected / not detected annotation
for i, (xi, det) in enumerate(zip(x, detected)):
    symbol = '✓' if det else '✗'
    color  = 'green' if det else 'red'
    ax1.text(xi, max(cr_peak[i], clutter_floor[i]) + 3,
             symbol, ha='center', va='bottom', fontsize=13,
             color=color, fontweight='bold')

# LiDAR annotation strip
for xi in x:
    ax1.text(xi, -75, 'LiDAR\n✗', ha='center', va='bottom',
             fontsize=7.5, color='#888888', style='italic')

ax1.set_xticks(x)
ax1.set_xticklabels(labels, fontsize=8.5)
ax1.set_ylabel('Amplitude (dB)', fontsize=10)
ax1.set_title('CR Peak vs Clutter Floor', fontsize=11, fontweight='bold')
ax1.set_ylim(-85, 75)
ax1.axhline(0, color='black', linewidth=0.5, linestyle='-')
ax1.legend(fontsize=8.5, loc='upper right')
ax1.grid(axis='y', alpha=0.3, linewidth=0.5)

# dashed line legend entry
ax1.hlines([], [], [], colors='black', linewidths=1.2,
           linestyles='--', label='3σ threshold')
ax1.legend(fontsize=8.5, loc='upper right')

# ── Right plot: SNR above clutter floor ──────────────────────────────────────
ax2 = axes[1]

bar_colors = [col_snr if d else '#FFAAAA' for d in detected]
bars_snr = ax2.bar(x, snr, width=0.45, color=bar_colors,
                   edgecolor='white', linewidth=0.5)

# Annotate bar tops with detected/not detected
for i, (xi, s, det) in enumerate(zip(x, snr, detected)):
    symbol = '✓' if det else '✗'
    color  = 'green' if det else 'red'
    ax2.text(xi, s + 1.5, symbol, ha='center', va='bottom',
             fontsize=13, color=color, fontweight='bold')

# LiDAR annotation
for xi in x:
    ax2.text(xi, -6, 'LiDAR\n✗', ha='center', va='top',
             fontsize=7.5, color='#888888', style='italic')

ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=8.5)
ax2.set_ylabel('SNR above clutter floor (dB)', fontsize=10)
ax2.set_title('SNR Comparison Across Experiments', fontsize=11,
              fontweight='bold')
ax2.set_ylim(-10, 100)
ax2.axhline(0, color='black', linewidth=0.5)
ax2.grid(axis='y', alpha=0.3, linewidth=0.5)

# Legend for detected/not detected
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=col_snr,   label='Detected above 3σ'),
    Patch(facecolor='#FFAAAA', label='Not detected above 3σ'),
]
ax2.legend(handles=legend_elements, fontsize=8.5, loc='upper right')

# ── Shared caption info ───────────────────────────────────────────────────────
fig.suptitle(
    'Cross-Experiment Radar Detection Summary  |  '
    'Dashed lines = 3σ detection threshold  |  '
    'LiDAR: no detection in any foliage scan',
    fontsize=8.5, color='#444444', y=0.01
)

plt.savefig('/mnt/user-data/outputs/cross_experiment_summary.pdf',
            dpi=300, bbox_inches='tight')
plt.savefig('/mnt/user-data/outputs/cross_experiment_summary.png',
            dpi=300, bbox_inches='tight')
print("Saved.")