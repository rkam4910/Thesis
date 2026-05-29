import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==========================================================
# CONFIG
# ==========================================================
SPARSE_CSV = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\gradar_logs\20260411-141145\pointcloud_20260411-141145_.csv'
THICK_CSV  = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\gradar_logs\Undetected_CR\pointcloud_20260411-130952_.csv'

FIGURES_DIR = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\gradar_logs\comparison_figures'
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

# ==========================================================
# SPARSE SCAN PARAMS (lower density face)
# ==========================================================
CR_X_S = 3.369311
CR_Y_S = 1.573796
CR_Z_S = 0.083744
CONE_S = 8.0
CR_TOL_S = 0.15
CR_RADIUS_S = 0.18
CLUTTER_START_S = 2.6
CLUTTER_END_S   = 3.4

# ==========================================================
# THICK SCAN PARAMS (higher density face)
# ==========================================================
CR_RANGE_T      = 3.72
CONE_T          = 15.0
CR_TOL_T        = 0.3
CLUTTER_START_T = 2.6
CLUTTER_END_T   = 3.4

# ==========================================================
# HELPERS
# ==========================================================
def load_filter_convert(csv_path, amp_min=-80, r_max=6.0,
                         near_cut=1.45, fence_deg=7.0):
    df = pd.read_csv(csv_path)
    df = df[
        (df['amplitude'] >= amp_min) &
        (df['range']     <= r_max)   &
        (df['range']     >= near_cut) &
        (np.degrees(df['roll']) >= fence_deg)
    ].copy()
    az = df['roll'].values
    el = (np.pi / 2) + df['pitch'].values
    r  = df['range'].values
    df['x'] = r * np.cos(el) * np.cos(az)
    df['y'] = r * np.cos(el) * np.sin(az)
    df['z'] = r * np.sin(el)
    return df

def apply_cone(df, boresight, half_angle):
    pts      = df[['x','y','z']].values
    pts_norm = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    cos_t    = np.clip(pts_norm @ boresight, -1, 1)
    mask     = np.degrees(np.arccos(cos_t)) <= half_angle
    return df[mask].copy()

def get_cr_clutter(cone, r_cone, amp_cone,
                   cr_range, cr_tol, cr_radius,
                   cr_xyz, clutter_start, clutter_end,
                   use_spatial_gate=True):
    if use_spatial_gate:
        xyz_dist = np.sqrt(
            (cone['x'].values - cr_xyz[0])**2 +
            (cone['y'].values - cr_xyz[1])**2 +
            (cone['z'].values - cr_xyz[2])**2
        )
        cr_mask = (
            (r_cone >= cr_range - cr_tol) &
            (r_cone <= cr_range + cr_tol) &
            (xyz_dist <= cr_radius)
        )
    else:
        cr_mask = (
            (r_cone >= cr_range - cr_tol) &
            (r_cone <= cr_range + cr_tol)
        )
    clutter_mask = (
        (r_cone >= clutter_start) &
        (r_cone <= clutter_end)   &
        (~cr_mask)
    )
    return cr_mask, clutter_mask

def detection_metrics(cr_vals, clutter_vals):
    cr_peak    = float(np.max(cr_vals))         if len(cr_vals)      else np.nan
    floor_med  = float(np.median(clutter_vals)) if len(clutter_vals) else np.nan
    floor_rob  = float(robust_sigma(clutter_vals)) if len(clutter_vals) else np.nan
    thr_3sig   = floor_med + 3 * floor_rob      if not np.isnan(floor_med) else np.nan
    snr        = cr_peak - floor_med            if not np.isnan(cr_peak)  else np.nan
    det_floor  = bool(cr_peak > floor_med)
    det_3sig   = bool(cr_peak > thr_3sig)
    return cr_peak, floor_med, floor_rob, thr_3sig, snr, det_floor, det_3sig

def roc_curve(cr_vals, clutter_vals, thr_3sig):
    all_s = np.concatenate([cr_vals, clutter_vals])
    thrs  = np.linspace(np.min(all_s), np.max(all_s), 200)
    tpr_l, fpr_l = [], []
    for t in thrs:
        tpr_l.append(np.mean(cr_vals      >= t) if len(cr_vals)      else 0.0)
        fpr_l.append(np.mean(clutter_vals >= t) if len(clutter_vals) else 0.0)
    tpr = np.array(tpr_l)
    fpr = np.array(fpr_l)
    order = np.argsort(fpr)
    auc   = np.trapezoid(tpr[order], fpr[order]) if len(fpr) > 1 else np.nan
    tpr_3 = np.mean(cr_vals      >= thr_3sig) if len(cr_vals)      else np.nan
    fpr_3 = np.mean(clutter_vals >= thr_3sig) if len(clutter_vals) else np.nan
    return fpr, tpr, auc, tpr_3, fpr_3

# ==========================================================
# LOAD + PROCESS SPARSE
# ==========================================================
df_s   = load_filter_convert(SPARSE_CSV)
cr_vec_s   = np.array([CR_X_S, CR_Y_S, CR_Z_S])
bore_s     = cr_vec_s / np.linalg.norm(cr_vec_s)
cone_s     = apply_cone(df_s, bore_s, CONE_S)
r_s        = cone_s['range'].values
amp_s      = cone_s['amplitude'].values
cr_range_s = np.linalg.norm(cr_vec_s)

cr_mask_s, cl_mask_s = get_cr_clutter(
    cone_s, r_s, amp_s,
    cr_range_s, CR_TOL_S, CR_RADIUS_S,
    cr_vec_s, CLUTTER_START_S, CLUTTER_END_S,
    use_spatial_gate=True
)
cr_vals_s  = amp_s[cr_mask_s]
cl_vals_s  = amp_s[cl_mask_s]

(cr_peak_s, floor_s, rob_s,
 thr_s, snr_s, det_f_s, det_3_s) = detection_metrics(cr_vals_s, cl_vals_s)

fpr_s, tpr_s, auc_s, tpr3_s, fpr3_s = roc_curve(cr_vals_s, cl_vals_s, thr_s)

# ==========================================================
# LOAD + PROCESS THICK
# ==========================================================
df_t  = load_filter_convert(THICK_CSV)
top   = df_t.nlargest(20, 'amplitude')
bore_t = top[['x','y','z']].mean().values
bore_t = bore_t / np.linalg.norm(bore_t)
cone_t = apply_cone(df_t, bore_t, CONE_T)
r_t    = cone_t['range'].values
amp_t  = cone_t['amplitude'].values

cr_mask_t, cl_mask_t = get_cr_clutter(
    cone_t, r_t, amp_t,
    CR_RANGE_T, CR_TOL_T, None,
    None, CLUTTER_START_T, CLUTTER_END_T,
    use_spatial_gate=False
)
cr_vals_t = amp_t[cr_mask_t]
cl_vals_t = amp_t[cl_mask_t]

(cr_peak_t, floor_t, rob_t,
 thr_t, snr_t, det_f_t, det_3_t) = detection_metrics(cr_vals_t, cl_vals_t)

fpr_t, tpr_t, auc_t, tpr3_t, fpr3_t = roc_curve(cr_vals_t, cl_vals_t, thr_t)

# ==========================================================
# PRINT SUMMARY
# ==========================================================
print("\n" + "="*60)
print("COMPARISON SUMMARY")
print("="*60)
print(f"{'Metric':<30} {'Sparse':>10} {'Thick':>10}")
print("-"*60)
print(f"{'CR peak (dB)':<30} {cr_peak_s:>+10.1f} {cr_peak_t:>+10.1f}")
print(f"{'Clutter median (dB)':<30} {floor_s:>+10.1f} {floor_t:>+10.1f}")
print(f"{'SNR (dB)':<30} {snr_s:>+10.1f} {snr_t:>+10.1f}")
print(f"{'3sigma threshold (dB)':<30} {thr_s:>+10.1f} {thr_t:>+10.1f}")
print(f"{'Detected above floor':<30} {'Yes':>10} {'Yes' if det_f_t else 'No':>10}")
print(f"{'Detected above 3sigma':<30} {'Yes' if det_3_s else 'No':>10} {'Yes' if det_3_t else 'No':>10}")
print(f"{'ROC AUC':<30} {auc_s:>10.3f} {auc_t:>10.3f}")
print("="*60)

# ==========================================================
# FIGURE 1 — Side by side histograms
# ==========================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

for ax, cr_v, cl_v, cr_pk, fl, thr, label in [
    (axes[0], cr_vals_s, cl_vals_s, cr_peak_s, floor_s, thr_s,
     'Lower density face (CR visible)'),
    (axes[1], cr_vals_t, cl_vals_t, cr_peak_t, floor_t, thr_t,
     'Higher density face (CR occluded)')
]:
    ax.hist(cl_v, bins=20, alpha=0.65, color='seagreen',
            edgecolor='darkgreen',
            label=f'Local clutter ({len(cl_v)} pts)')
    ax.hist(cr_v, bins=12, alpha=0.82, color='red',
            edgecolor='darkred',
            label=f'CR window ({len(cr_v)} pts)')
    ax.axvline(fl,  color='black',  lw=1.8, ls='-',
               label=f'Clutter median ({fl:+.1f} dB)')
    ax.axvline(thr, color='purple', lw=1.4, ls='--',
               label=f'3σ threshold ({thr:+.1f} dB)')
    ax.axvline(cr_pk, color='red',  lw=1.4, ls=':',
               label=f'CR peak ({cr_pk:+.1f} dB)')
    ax.set_xlabel('Amplitude (dB)', fontsize=11)
    ax.set_title(label, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

axes[0].set_ylabel('Count', fontsize=11)
fig.suptitle('Amplitude Distribution Comparison\nCherry Ballart — Lower vs Higher Density Face',
             fontsize=12, fontweight='bold')
plt.tight_layout()
save('fig_comparison_histograms')

# ==========================================================
# FIGURE 2 — Side by side ROC curves
# ==========================================================
fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(fpr_s, tpr_s, color='green', lw=2,
        label=f'Lower density (AUC = {auc_s:.3f})')
ax.plot(fpr_t, tpr_t, color='red',   lw=2,
        label=f'Higher density (AUC = {auc_t:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5,
        label='Random classifier')
ax.plot(fpr3_s, tpr3_s, 'go', markersize=10,
        label=f'3σ lower density (TPR={tpr3_s:.2f}, FPR={fpr3_s:.2f})')
ax.plot(fpr3_t, tpr3_t, 'rs', markersize=10,
        label=f'3σ higher density (TPR={tpr3_t:.2f}, FPR={fpr3_t:.2f})')

ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curve Comparison\nCherry Ballart — Lower vs Higher Density Face',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
plt.tight_layout()
save('fig_comparison_roc')