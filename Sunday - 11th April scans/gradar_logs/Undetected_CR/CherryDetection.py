import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==========================================================
# FILE
# ==========================================================
CSV_PATH = r'C:\Users\raque\OneDrive - The University of Sydney (Students)\Documents\University of Sydney\Year 4\Thesis\2026 Scans\Sunday - 11th April scans\gradar_logs\Undetected_CR\pointcloud_20260411-130952_.csv'

# ==========================================================
# SETTINGS
# ==========================================================
AMP_MIN = -80
R_MAX = 9.0
NEAR_FIELD_CUT = 1.45
FENCE_ROLL_CUT_DEG = 7.0

# You moved the reflector slightly right, so search this region
ROI_X_MIN, ROI_X_MAX = 3.2, 4.6
ROI_Y_MIN, ROI_Y_MAX = 1.2, 2.0
ROI_Z_MIN, ROI_Z_MAX = -0.2, 0.4

# Local clutter region around expected CR depth
CLUTTER_X_MIN, CLUTTER_X_MAX = 2.7, 4.6
CLUTTER_Y_MIN, CLUTTER_Y_MAX = 1.0, 2.2
CLUTTER_Z_MIN, CLUTTER_Z_MAX = -0.4, 0.5

# Threshold sweep for bright-cluster candidate selection
CANDIDATE_THRESHOLDS = [-10, -5, 0, 5]

# Bush attenuation settings
SHRUB_FRONT = 1.5
SHRUB_BACK = 4.5
BIN_WIDTH = 0.25
MIN_POINTS_BIN = 8

# Folder to save figures
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(SCRIPT_DIR, "figures", "undetected_cr_analysis")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ==========================================================
# HELPERS
# ==========================================================
def save(name):
    path = os.path.join(FIGURES_DIR, f"{name}.png")
    plt.savefig(path, dpi=160, bbox_inches="tight")
    print(f"saved: {path}")

def robust_sigma(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return 1.4826 * mad

def classify_detection(peak, clutter_median, sigma_robust):
    if np.isnan(peak) or np.isnan(clutter_median) or np.isnan(sigma_robust):
        return "insufficient data"
    thr3 = clutter_median + 3 * sigma_robust
    if peak > thr3:
        return "likely detected"
    elif peak > clutter_median + 1.5 * sigma_robust:
        return "weak / ambiguous"
    else:
        return "not detectable"

def range_binned_stats(rr, aa, lo, hi, bin_width, min_pts):
    edges = np.arange(lo, hi + bin_width, bin_width)
    rows = []
    for e in edges[:-1]:
        sel = (rr >= e) & (rr < e + bin_width)
        if np.sum(sel) >= min_pts:
            rows.append({
                "r_mid": e + bin_width / 2,
                "n": int(np.sum(sel)),
                "amp_mean": float(np.mean(aa[sel])),
                "amp_std": float(np.std(aa[sel])),
                "amp_median": float(np.median(aa[sel])),
            })
    return pd.DataFrame(rows)

# ==========================================================
# LOAD
# ==========================================================
df = pd.read_csv(CSV_PATH)
df = df[
    (df["amplitude"] >= AMP_MIN) &
    (df["range"] <= R_MAX) &
    (df["range"] >= NEAR_FIELD_CUT)
].copy()

if FENCE_ROLL_CUT_DEG > 0:
    df = df[np.degrees(df["roll"]) >= FENCE_ROLL_CUT_DEG].copy()

# ==========================================================
# SPHERICAL -> CARTESIAN (same as MATLAB)
# ==========================================================
az = df["roll"].values
el = (np.pi / 2) + df["pitch"].values
r = df["range"].values

df["x"] = r * np.cos(el) * np.cos(az)
df["y"] = r * np.cos(el) * np.sin(az)
df["z"] = r * np.sin(el)

# ==========================================================
# ROI = possible moved CR region
# ==========================================================
roi = df[
    (df["x"] >= ROI_X_MIN) & (df["x"] <= ROI_X_MAX) &
    (df["y"] >= ROI_Y_MIN) & (df["y"] <= ROI_Y_MAX) &
    (df["z"] >= ROI_Z_MIN) & (df["z"] <= ROI_Z_MAX)
].copy()

# local clutter in surrounding bush region
clutter = df[
    (df["x"] >= CLUTTER_X_MIN) & (df["x"] <= CLUTTER_X_MAX) &
    (df["y"] >= CLUTTER_Y_MIN) & (df["y"] <= CLUTTER_Y_MAX) &
    (df["z"] >= CLUTTER_Z_MIN) & (df["z"] <= CLUTTER_Z_MAX)
].copy()

# remove ROI from clutter so clutter is background only
clutter = clutter[
    ~(
        (clutter["x"] >= ROI_X_MIN) & (clutter["x"] <= ROI_X_MAX) &
        (clutter["y"] >= ROI_Y_MIN) & (clutter["y"] <= ROI_Y_MAX) &
        (clutter["z"] >= ROI_Z_MIN) & (clutter["z"] <= ROI_Z_MAX)
    )
].copy()

# ==========================================================
# CLUTTER METRICS
# ==========================================================
if len(clutter) > 0:
    clutter_median = float(np.median(clutter["amplitude"]))
    clutter_mean = float(np.mean(clutter["amplitude"]))
    clutter_sigma = float(np.std(clutter["amplitude"]))
    clutter_robust = float(robust_sigma(clutter["amplitude"].values))
else:
    clutter_median = clutter_mean = clutter_sigma = clutter_robust = np.nan

threshold_3sigma = clutter_median + 3 * clutter_robust if not np.isnan(clutter_median) else np.nan

# ==========================================================
# CANDIDATE SEARCH
# ==========================================================
best_candidate = None

if len(roi) > 0:
    brightest = roi.loc[roi["amplitude"].idxmax()]
    best_candidate = {
        "method": "brightest point in ROI",
        "x": float(brightest["x"]),
        "y": float(brightest["y"]),
        "z": float(brightest["z"]),
        "range": float(brightest["range"]),
        "amp": float(brightest["amplitude"]),
        "npts": 1
    }

cluster_results = []
for thr in CANDIDATE_THRESHOLDS:
    pts = roi[roi["amplitude"] >= thr].copy()
    if len(pts) == 0:
        continue
    cluster_results.append({
        "threshold": thr,
        "npts": len(pts),
        "x": float(pts["x"].mean()),
        "y": float(pts["y"].mean()),
        "z": float(pts["z"].mean()),
        "amp_mean": float(pts["amplitude"].mean()),
        "amp_peak": float(pts["amplitude"].max())
    })

best_cluster = None
if len(cluster_results) > 0:
    best_cluster = max(cluster_results, key=lambda d: d["amp_peak"])

if best_candidate is not None:
    peak = best_candidate["amp"]
    margin = peak - clutter_median if not np.isnan(clutter_median) else np.nan
    decision = classify_detection(peak, clutter_median, clutter_robust)
else:
    peak = np.nan
    margin = np.nan
    decision = "no candidate in ROI"

# ==========================================================
# ROC USING ROI AS POSITIVES
# ==========================================================
roi_vals = roi["amplitude"].values
clutter_vals = clutter["amplitude"].values

all_scores = np.concatenate([roi_vals, clutter_vals]) if (len(roi_vals) + len(clutter_vals)) > 0 else np.array([])
thresholds = np.linspace(np.min(all_scores), np.max(all_scores), 200) if len(all_scores) > 0 else np.array([])

tpr, fpr = [], []
for t in thresholds:
    tpr.append(np.mean(roi_vals >= t) if len(roi_vals) else 0.0)
    fpr.append(np.mean(clutter_vals >= t) if len(clutter_vals) else 0.0)

tpr = np.array(tpr)
fpr = np.array(fpr)

auc = np.trapz(tpr[np.argsort(fpr)], fpr[np.argsort(fpr)]) if len(fpr) > 1 else np.nan
tpr_3sigma = np.mean(roi_vals >= threshold_3sigma) if len(roi_vals) and not np.isnan(threshold_3sigma) else np.nan
fpr_3sigma = np.mean(clutter_vals >= threshold_3sigma) if len(clutter_vals) and not np.isnan(threshold_3sigma) else np.nan

# ==========================================================
# BUSH-ONLY ATTENUATION ANALYSIS
# Exclude candidate region in range
# ==========================================================
if best_candidate is not None:
    candidate_range = best_candidate["range"]
else:
    candidate_range = 3.7

CR_TOL = 0.18

r_all = df["range"].values
amp_all = df["amplitude"].values

bush_only_mask = (
    (r_all >= SHRUB_FRONT) &
    (r_all <= SHRUB_BACK) &
    ~((r_all >= candidate_range - CR_TOL) & (r_all <= candidate_range + CR_TOL))
)

r_bush = r_all[bush_only_mask]
amp_bush = amp_all[bush_only_mask]

bush_profile = range_binned_stats(
    r_bush, amp_bush,
    SHRUB_FRONT, SHRUB_BACK,
    BIN_WIDTH, MIN_POINTS_BIN
)

if len(bush_profile) >= 3:
    slope_mean_bush, intercept_mean_bush, rval_mean_bush, pval_mean_bush, se_mean_bush = linregress(
        bush_profile["r_mid"], bush_profile["amp_mean"]
    )
    slope_median_bush, intercept_median_bush, rval_median_bush, pval_median_bush, se_median_bush = linregress(
        bush_profile["r_mid"], bush_profile["amp_median"]
    )

    k_mean_bush = slope_mean_bush
    k_median_bush = slope_median_bush
    r2_mean_bush = rval_mean_bush**2
    r2_median_bush = rval_median_bush**2

    r_line_bush = np.linspace(SHRUB_FRONT, SHRUB_BACK, 100)
else:
    slope_mean_bush = intercept_mean_bush = rval_mean_bush = pval_mean_bush = se_mean_bush = np.nan
    slope_median_bush = intercept_median_bush = rval_median_bush = pval_median_bush = se_median_bush = np.nan
    k_mean_bush = k_median_bush = r2_mean_bush = r2_median_bush = np.nan
    r_line_bush = np.array([])

# ==========================================================
# PRINT RESULTS
# ==========================================================
print("=" * 72)
print("UNDETECTED_CR SCAN ANALYSIS")
print("=" * 72)
print(f"Total filtered points:        {len(df):,}")
print(f"ROI points:                   {len(roi):,}")
print(f"Clutter points:               {len(clutter):,}")
print()
print(f"Clutter median:               {clutter_median:+.1f} dB")
print(f"Clutter mean:                 {clutter_mean:+.1f} dB")
print(f"Clutter sigma:                {clutter_sigma:.1f} dB")
print(f"Clutter robust sigma:         {clutter_robust:.1f} dB")
print(f"3σ threshold:                 {threshold_3sigma:+.1f} dB")
print()

if best_candidate is not None:
    print("Best candidate:")
    print(f"  X = {best_candidate['x']:.4f} m")
    print(f"  Y = {best_candidate['y']:.4f} m")
    print(f"  Z = {best_candidate['z']:.4f} m")
    print(f"  Range = {best_candidate['range']:.4f} m")
    print(f"  Peak amplitude = {best_candidate['amp']:+.1f} dB")
    print(f"  Margin above clutter median = {margin:+.1f} dB")
    print(f"  Decision: {decision}")
else:
    print("No candidate found in ROI.")

if best_cluster is not None:
    print()
    print("Best thresholded bright cluster:")
    print(f"  Threshold = {best_cluster['threshold']} dB")
    print(f"  Points = {best_cluster['npts']}")
    print(f"  Centroid = ({best_cluster['x']:.4f}, {best_cluster['y']:.4f}, {best_cluster['z']:.4f})")
    print(f"  Mean amplitude = {best_cluster['amp_mean']:+.1f} dB")
    print(f"  Peak amplitude = {best_cluster['amp_peak']:+.1f} dB")

print()
print("Bush-only attenuation:")
print(f"  Bush-only points = {len(r_bush):,}")
if not np.isnan(k_mean_bush):
    print(f"  Mean slope = {k_mean_bush:+.2f} dB/m (R² = {r2_mean_bush:.2f})")
    print(f"  Median slope = {k_median_bush:+.2f} dB/m (R² = {r2_median_bush:.2f})")
else:
    print("  Not enough bush-only bins for attenuation fit.")

print(f"ROC AUC = {auc:.3f}")
print("=" * 72)

# ==========================================================
# PLOT 1 -- plan view
# ==========================================================
fig, ax = plt.subplots(figsize=(8, 7))
sc = ax.scatter(df["x"], df["y"], c=df["amplitude"], s=10, cmap="viridis", alpha=0.65)

ax.add_patch(plt.Rectangle(
    (ROI_X_MIN, ROI_Y_MIN),
    ROI_X_MAX - ROI_X_MIN,
    ROI_Y_MAX - ROI_Y_MIN,
    fill=False, edgecolor="red", linewidth=2, label="Search ROI"
))

if best_candidate is not None:
    ax.scatter(best_candidate["x"], best_candidate["y"],
               s=180, c="yellow", edgecolors="black", linewidths=1.2,
               label="Best candidate")
    ax.text(best_candidate["x"] + 0.05, best_candidate["y"] + 0.03,
            f"Candidate\n({best_candidate['x']:.2f}, {best_candidate['y']:.2f})",
            color="red", fontsize=10)

if best_cluster is not None:
    ax.scatter(best_cluster["x"], best_cluster["y"],
               s=180, c="red", marker="*", edgecolors="black", linewidths=1.0,
               label="Bright cluster centroid")

ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_title("Undetected_CR Scan -- Plan View Candidate Search")
ax.grid(True, alpha=0.3)
ax.set_aspect("equal")
ax.legend()
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Amplitude (dB)")
plt.tight_layout()
save("fig1_plan_view_candidate_search")
plt.show()

# ==========================================================
# PLOT 2 -- histogram: ROI vs clutter
# ==========================================================
fig, ax = plt.subplots(figsize=(9, 5))

if len(clutter) > 0:
    ax.hist(clutter["amplitude"], bins=25, alpha=0.65, color="seagreen",
            edgecolor="darkgreen", label=f"Local clutter ({len(clutter)})")

if len(roi) > 0:
    ax.hist(roi["amplitude"], bins=20, alpha=0.75, color="red",
            edgecolor="darkred", label=f"Search ROI ({len(roi)})")

if not np.isnan(clutter_median):
    ax.axvline(clutter_median, color="black", lw=1.6, label=f"Clutter median ({clutter_median:+.1f} dB)")
if not np.isnan(threshold_3sigma):
    ax.axvline(threshold_3sigma, color="purple", lw=1.4, ls="--", label=f"3σ threshold ({threshold_3sigma:+.1f} dB)")
if best_candidate is not None:
    ax.axvline(best_candidate["amp"], color="red", lw=1.4, ls=":", label=f"Best candidate ({best_candidate['amp']:+.1f} dB)")

ax.set_xlabel("Amplitude (dB)")
ax.set_ylabel("Count")
ax.set_title("Undetected_CR Scan -- Candidate Detectability")
ax.grid(True, alpha=0.3, axis="y")
ax.legend()
plt.tight_layout()
save("fig2_histogram_candidate_detectability")
plt.show()

# ==========================================================
# PLOT 3 -- ROC using ROI as positives
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(fpr, tpr, color="navy", lw=2, label=f"ROI ROC (AUC = {auc:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5, label="Random classifier")
if not np.isnan(tpr_3sigma):
    ax.plot(fpr_3sigma, tpr_3sigma, "ro", markersize=9,
            label=f"3σ point (TPR={tpr_3sigma:.2f}, FPR={fpr_3sigma:.2f})")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("Undetected_CR Scan -- ROC Detectability")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
save("fig3_roc_detectability")
plt.show()

# ==========================================================
# PLOT 4 -- Bush-only attenuation
# ==========================================================
fig, ax = plt.subplots(figsize=(10, 5))

ax.scatter(r_bush, amp_bush, s=10, alpha=0.2,
           color="seagreen", label="Bush returns")

if len(bush_profile) > 0:
    ax.errorbar(
        bush_profile["r_mid"], bush_profile["amp_mean"],
        yerr=bush_profile["amp_std"],
        fmt="o", color="navy",
        capsize=3, label="Binned mean ±1σ"
    )

if not np.isnan(k_mean_bush):
    yy = intercept_mean_bush + slope_mean_bush * r_line_bush
    ax.plot(r_line_bush, yy, "r--",
            label=f"Fit = {k_mean_bush:+.2f} dB/m (R²={r2_mean_bush:.2f})")

ax.axvspan(candidate_range - CR_TOL, candidate_range + CR_TOL,
           color="red", alpha=0.18,
           label="Excluded candidate region")

ax.set_xlabel("Range (m)")
ax.set_ylabel("Amplitude (dB)")
ax.set_title("Bush-Only Attenuation (Candidate Excluded)")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
save("fig4_bush_only_attenuation")
plt.show()

# ==========================================================
# PLOT 5 -- Bush-only median profile
# ==========================================================
fig, ax = plt.subplots(figsize=(10, 5))

if len(bush_profile) > 0:
    ax.plot(
        bush_profile["r_mid"],
        bush_profile["amp_median"],
        "o-",
        color="darkgreen",
        label="Bush-only binned median amplitude"
    )

if not np.isnan(k_median_bush):
    ax.plot(
        r_line_bush,
        intercept_median_bush + slope_median_bush * r_line_bush,
        "r--",
        lw=2,
        label=f"Median fit: {k_median_bush:+.2f} dB/m (R² = {r2_median_bush:.2f})"
    )

ax.axvspan(SHRUB_FRONT, SHRUB_BACK, alpha=0.12, color="orange", label="Shrub body")
ax.axvspan(candidate_range - CR_TOL, candidate_range + CR_TOL, alpha=0.18, color="red", label="Excluded candidate region")
ax.set_xlabel("Range (m)")
ax.set_ylabel("Median amplitude (dB)")
ax.set_title("Bush-Only Median Amplitude Gradient")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
save("fig5_bush_only_median_gradient")
plt.show()

print(f"\nAll figures saved to: {FIGURES_DIR}")