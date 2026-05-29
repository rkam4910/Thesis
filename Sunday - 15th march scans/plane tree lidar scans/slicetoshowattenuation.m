% =========================================================================
% compare_radar_lidar.m
%
% Side-by-side comparison of radar and LiDAR scans of the same plane tree.
% Radar:  pointcloud_20260315-143039_.csv
% LiDAR:  2026-03-15 14-52-14_corner_small.csv
%
% Radar coordinate convention matches existing radar script (sph2cart with
% elevation offset). LiDAR coordinates are already Cartesian (metres).
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
file_radar = 'pointcloud_20260315-143039_.csv';
file_lidar = '2026-03-15 14-52-14_corner_small.csv';

RANGE_LIMIT  = 50;    % m — applied to both sensors
AMP_THRESH   = -80;   % dB — radar amplitude filter (matches existing script)
POINT_SIZE   = 4;     % scatter point size

% =========================================================================
fprintf('=== LOAD RADAR ===\n');
% =========================================================================
data = readmatrix(file_radar);

% Filter by amplitude and range (matches existing radar script)
data = data(data(:,5) >= AMP_THRESH, :);
data = data(data(:,2) <= RANGE_LIMIT, :);

az  = data(:,4);             % azimuth (rad)
el  = (pi/2) + data(:,3);   % elevation — matches existing script convention
r   = data(:,2);             % range (m)
amp = data(:,5);             % returned power (dB)

[xr, yr, zr] = sph2cart(az, el, r);

fprintf('  %d radar points loaded\n', numel(xr));
fprintf('  X: %.2f to %.2f m\n', min(xr), max(xr));
fprintf('  Z: %.2f to %.2f m\n', min(zr), max(zr));

% =========================================================================
fprintf('\n=== LOAD LIDAR ===\n');
% =========================================================================
T = readtable(file_lidar, 'VariableNamingRule', 'preserve');
get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};

xl   = get_col('X');
yl   = get_col('Y');
zl   = get_col('Z');
refl = get_col('Reflectivity');

% Remove invalid zero returns and apply range filter
valid   = ~(xl == 0 & yl == 0 & zl == 0);
xl = xl(valid); yl = yl(valid); zl = zl(valid); refl = refl(valid);

range_l = sqrt(xl.^2 + yl.^2 + zl.^2);
keep    = range_l <= RANGE_LIMIT;
xl = xl(keep); yl = yl(keep); zl = zl(keep); refl = refl(keep);

fprintf('  %d LiDAR points loaded\n', numel(xl));
fprintf('  X: %.2f to %.2f m\n', min(xl), max(xl));
fprintf('  Z: %.2f to %.2f m\n', min(zl), max(zl));

% =========================================================================
fprintf('\n=== PLOTS ===\n');
% =========================================================================

%% Figure 1: Side view (XZ) — key penetration comparison
figure('Name','Side view XZ — Radar vs LiDAR','NumberTitle','off', ...
    'Color','w','Units','normalized','Position',[0.02 0.3 0.96 0.5]);

ax1 = subplot(1,2,1);
scatter(ax1, xr, zr, POINT_SIZE, amp, 'filled');
colormap(ax1,'jet'); clim(ax1,[-100 100]);
cb1 = colorbar(ax1); cb1.Label.String = 'Amplitude (dB)';
axis(ax1,'equal'); grid(ax1,'on');
xlabel(ax1,'X — depth (m)'); ylabel(ax1,'Z — height (m)');
title(ax1,'Radar 120 GHz — side view (XZ)','FontSize',13,'FontWeight','bold');

ax2 = subplot(1,2,2);
scatter(ax2, xl, zl, POINT_SIZE, refl, 'filled');
colormap(ax2,'turbo'); clim(ax2,[0 255]);
cb2 = colorbar(ax2); cb2.Label.String = 'Reflectivity (0–255)';
axis(ax2,'equal'); grid(ax2,'on');
xlabel(ax2,'X — depth (m)'); ylabel(ax2,'Z — height (m)');
title(ax2,'LiDAR 905 nm — side view (XZ)','FontSize',13,'FontWeight','bold');

sgtitle('Side view — radar penetration vs LiDAR surface', ...
    'FontSize',15,'FontWeight','bold');

%% Figure 2: Top-down (XY) comparison
figure('Name','Top-down XY — Radar vs LiDAR','NumberTitle','off', ...
    'Color','w','Units','normalized','Position',[0.02 0.3 0.96 0.5]);

ax3 = subplot(1,2,1);
scatter(ax3, xr, yr, POINT_SIZE, amp, 'filled');
colormap(ax3,'jet'); clim(ax3,[-100 100]);
cb3 = colorbar(ax3); cb3.Label.String = 'Amplitude (dB)';
axis(ax3,'equal'); grid(ax3,'on');
xlabel(ax3,'X (m)'); ylabel(ax3,'Y (m)');
title(ax3,'Radar 120 GHz — top-down (XY)','FontSize',13,'FontWeight','bold');

ax4 = subplot(1,2,2);
scatter(ax4, xl, yl, POINT_SIZE, refl, 'filled');
colormap(ax4,'turbo'); clim(ax4,[0 255]);
cb4 = colorbar(ax4); cb4.Label.String = 'Reflectivity (0–255)';
axis(ax4,'equal'); grid(ax4,'on');
xlabel(ax4,'X (m)'); ylabel(ax4,'Y (m)');
title(ax4,'LiDAR 905 nm — top-down (XY)','FontSize',13,'FontWeight','bold');

sgtitle('Top-down comparison','FontSize',15,'FontWeight','bold');

%% Figure 3: 3D side-by-side
figure('Name','3D — Radar vs LiDAR','NumberTitle','off', ...
    'Color','w','Units','normalized','Position',[0.02 0.05 0.96 0.8]);

ax5 = subplot(1,2,1);
scatter3(ax5, xr, yr, zr, POINT_SIZE, amp, 'filled');
colormap(ax5,'jet'); clim(ax5,[-100 100]);
cb5 = colorbar(ax5); cb5.Label.String = 'Amplitude (dB)';
axis(ax5,'equal'); grid(ax5,'on'); view(ax5,45,30);
xlabel(ax5,'X (m)'); ylabel(ax5,'Y (m)'); zlabel(ax5,'Z (m)');
title(ax5,'Radar 120 GHz — 3D','FontSize',13,'FontWeight','bold');
rotate3d(ax5,'on');

ax6 = subplot(1,2,2);
scatter3(ax6, xl, yl, zl, POINT_SIZE, refl, 'filled');
colormap(ax6,'turbo'); clim(ax6,[0 255]);
cb6 = colorbar(ax6); cb6.Label.String = 'Reflectivity (0–255)';
axis(ax6,'equal'); grid(ax6,'on'); view(ax6,45,30);
xlabel(ax6,'X (m)'); ylabel(ax6,'Y (m)'); zlabel(ax6,'Z (m)');
title(ax6,'LiDAR 905 nm — 3D','FontSize',13,'FontWeight','bold');
rotate3d(ax6,'on');

sgtitle('3D point cloud comparison','FontSize',15,'FontWeight','bold');

%% Figure 4: GPS-aligned world frame Z slice (6-9m)
% Both sensors are first transformed into the common GPS world frame
% using the 2D rigid body transforms derived from corner reflector GPS.
% A Z slice is then valid because both clouds share the same vertical datum.

fprintf('\n=== SECTION 4: GPS ALIGNMENT + Z SLICE ===\n');

% ── GPS world frame (ENU, origin = Scan 1 radar GPS position) ──────────
% Values from gps_tree_scan_1.csv
olon = 151.19214356; olat = -33.89005759;

% Scan 1 reflector GPS -> ENU (metres)
lon1_small = 151.19207230; lat1_small = -33.88999004;
lon1_big   = 151.19211892; lat1_big   = -33.88994154;

gps2enu_xy = @(lon, lat) deal( ...
    (lon - olon) .* 111320 .* cos(deg2rad(olat)), ...
    (lat - olat) .* 110540);

[w1_small_e, w1_small_n] = gps2enu_xy(lon1_small, lat1_small);
[w1_big_e,   w1_big_n]   = gps2enu_xy(lon1_big,   lat1_big);

fprintf('  Scan 1 GPS world positions:\n');
fprintf('    Small reflector: E=%.3f  N=%.3f  (horiz dist=%.2fm)\n', ...
    w1_small_e, w1_small_n, norm([w1_small_e, w1_small_n]));
fprintf('    Big   reflector: E=%.3f  N=%.3f  (horiz dist=%.2fm)\n', ...
    w1_big_e,   w1_big_n,   norm([w1_big_e,   w1_big_n]));

% ── Radar reflector centroids (from high-amplitude cluster) ───────────
amp_full  = data(:,5);
r_full    = data(:,2);
az_full   = data(:,4);
el_full   = (pi/2) + data(:,3);
[xr_all, yr_all, ~] = sph2cart(az_full, el_full, r_full);

hi_mask = amp_full > 50;
r_hi    = r_full(hi_mask);
x_hi    = xr_all(hi_mask);
y_hi    = yr_all(hi_mask);

% Split by largest range gap
[r_sort, si] = sort(r_hi);
gaps = diff(r_sort);
[~, gi] = max(gaps);
near_idx = si(1:gi);
far_idx  = si(gi+1:end);

c1_near = [mean(x_hi(near_idx)), mean(y_hi(near_idx))];
c1_far  = [mean(x_hi(far_idx)),  mean(y_hi(far_idx))];

r1_near = mean(r_hi(near_idx));
r1_far  = mean(r_hi(far_idx));

% Match by range
if abs(r1_near - norm([w1_small_e, w1_small_n])) < abs(r1_near - norm([w1_big_e, w1_big_n]))
    s1_near_wf = [w1_small_e, w1_small_n];
    s1_far_wf  = [w1_big_e,   w1_big_n];
else
    s1_near_wf = [w1_big_e,   w1_big_n];
    s1_far_wf  = [w1_small_e, w1_small_n];
end

fprintf('  Radar centroid near: x=%.3f y=%.3f (range %.2fm)\n', c1_near, r1_near);
fprintf('  Radar centroid far:  x=%.3f y=%.3f (range %.2fm)\n', c1_far,  r1_far);

% ── 2D SVD rigid body transform for radar ────────────────────────────
sp1 = [c1_near; c1_far; 0, 0];
wp1 = [s1_near_wf; s1_far_wf; 0, 0];

mu_s = mean(sp1,1); mu_w = mean(wp1,1);
S    = (sp1 - mu_s)' * (wp1 - mu_w);
[U,~,V] = svd(S);
R2d  = V * diag([1, det(V*U')]) * U';
theta1 = atan2(R2d(2,1), R2d(1,1));
t1_2d  = mu_w - mu_s * R2d';

c1 = cos(theta1); s1_r = sin(theta1);
radar_xy_world = [xr, yr] * [c1 s1_r; -s1_r c1] + t1_2d;
xr_w = radar_xy_world(:,1);
yr_w = radar_xy_world(:,2);
zr_w = -zr;   % flip Z so up is positive

fprintf('  Radar yaw = %.2f deg  tx=%.3f  ty=%.3f\n', ...
    rad2deg(theta1), t1_2d(1), t1_2d(2));

% ── LiDAR reflector centroids ─────────────────────────────────────────
refl_hi = refl >= 200;
range_hi = sqrt(xl(refl_hi).^2 + yl(refl_hi).^2 + zl(refl_hi).^2);
x_lhi = xl(refl_hi); y_lhi = yl(refl_hi);

[r_lsort, lsi] = sort(range_hi);
lgaps = diff(r_lsort);
[~, lgi] = max(lgaps);
l_near_idx = lsi(1:lgi);
l_far_idx  = lsi(lgi+1:end);

cl_near = [mean(x_lhi(l_near_idx)), mean(y_lhi(l_near_idx))];
cl_far  = [mean(x_lhi(l_far_idx)),  mean(y_lhi(l_far_idx))];
rl_near = mean(r_lsort(1:lgi));
rl_far  = mean(r_lsort(lgi+1:end));

if abs(rl_near - norm([w1_small_e, w1_small_n])) < abs(rl_near - norm([w1_big_e, w1_big_n]))
    sl_near_wf = [w1_small_e, w1_small_n];
    sl_far_wf  = [w1_big_e,   w1_big_n];
else
    sl_near_wf = [w1_big_e,   w1_big_n];
    sl_far_wf  = [w1_small_e, w1_small_n];
end

fprintf('  LiDAR centroid near: x=%.3f y=%.3f (range %.2fm)\n', cl_near, rl_near);
fprintf('  LiDAR centroid far:  x=%.3f y=%.3f (range %.2fm)\n', cl_far,  rl_far);

% ── 2D SVD for LiDAR ──────────────────────────────────────────────────
sp2 = [cl_near; cl_far; 0, 0];
wp2 = [sl_near_wf; sl_far_wf; 0, 0];

mu_s2 = mean(sp2,1); mu_w2 = mean(wp2,1);
S2    = (sp2 - mu_s2)' * (wp2 - mu_w2);
[U2,~,V2] = svd(S2);
R2d2  = V2 * diag([1, det(V2*U2')]) * U2';
theta2 = atan2(R2d2(2,1), R2d2(1,1));
t2_2d  = mu_w2 - mu_s2 * R2d2';

c2 = cos(theta2); s2_r = sin(theta2);
lidar_xy_world = [xl, yl] * [c2 s2_r; -s2_r c2] + t2_2d;
xl_w = lidar_xy_world(:,1);
yl_w = lidar_xy_world(:,2);
zl_w = zl;   % LiDAR Z already upward positive

fprintf('  LiDAR yaw = %.2f deg  tx=%.3f  ty=%.3f\n', ...
    rad2deg(theta2), t2_2d(1), t2_2d(2));

% ── Z slice in world frame ─────────────────────────────────────────────
Z_LO = 6; Z_HI = 9;

radar_z_mask = zr_w >= Z_LO & zr_w <= Z_HI;
lidar_z_mask = zl_w >= Z_LO & zl_w <= Z_HI;

fprintf('\n  Radar points in world Z slice (%.0f–%.0fm): %d\n', Z_LO, Z_HI, sum(radar_z_mask));
fprintf('  LiDAR points in world Z slice (%.0f–%.0fm): %d\n', Z_LO, Z_HI, sum(lidar_z_mask));

figure('Name', sprintf('World frame Z slice %.0f–%.0fm — LiDAR only', Z_LO, Z_HI), ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [0.02 0.3 0.60 0.55]);

ax7 = axes;
plot(ax7, 0, 0, 'kx', 'MarkerSize', 14, 'LineWidth', 2); hold(ax7,'on');
scatter(ax7, xl_w(lidar_z_mask), yl_w(lidar_z_mask), POINT_SIZE+2, ...
    refl(lidar_z_mask), 'filled');
colormap(ax7,'turbo'); clim(ax7,[0 255]);
cb8 = colorbar(ax7); cb8.Label.String = 'Reflectivity (0–255)';
axis(ax7,'equal'); grid(ax7,'on');
xlabel(ax7,'East (m)'); ylabel(ax7,'North (m)');
title(ax7, sprintf('LiDAR 905 nm — horizontal slice Z = %.0f–%.0f m (GPS world frame)', ...
    Z_LO, Z_HI), 'FontSize', 13, 'FontWeight', 'bold');
legend(ax7, 'LiDAR origin', 'FontSize', 9);

% Annotation explaining radar exclusion
annotation('textbox', [0.02 0.01 0.96 0.08], ...
    'String', ['Note: Radar is excluded from this slice. The 120 GHz radar scans ' ...
               'downward at steep elevation angles (−30° to −99°), so a horizontal ' ...
               'Z slice does not correspond to a valid cross-section of its scan geometry. ' ...
               'The side view (XZ) in Figure 1 is used for radar penetration comparison instead.'], ...
    'FontSize', 9, 'EdgeColor', [0.6 0.6 0.6], 'BackgroundColor', [0.97 0.97 0.97], ...
    'FitBoxToText', 'off', 'HorizontalAlignment', 'left');

% ── Figure 5: Comparison at radar -70°, -75°, -80° vs LiDAR equivalents ─
COMPARE_ANGLES = [-70, -75, -80];
PITCH_TOL      = 2;
n              = numel(COMPARE_ANGLES);

pitch_deg_all = rad2deg(data(:,3));
range_l_all   = sqrt(xl.^2 + yl.^2 + zl.^2);
el_lidar_deg  = rad2deg(asin(zl ./ (range_l_all + eps)));

figure('Name', 'Radar vs LiDAR — 70, 75, 80 degree comparison', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', ...
    'Position', [0.01 0.05 0.98 0.85]);

for k = 1:n
    pd       = COMPARE_ANGLES(k);
    lidar_el = 90 + pd;   % equivalent LiDAR elevation angle

    r_mask = pitch_deg_all >= (pd       - PITCH_TOL) & pitch_deg_all <= (pd       + PITCH_TOL);
    l_mask = el_lidar_deg  >= (lidar_el - PITCH_TOL) & el_lidar_deg  <= (lidar_el + PITCH_TOL);

    % ── Radar (top row) ──────────────────────────────────────────────────
    ax_r = subplot(2, n, k);
    plot(ax_r, 0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2); hold(ax_r,'on');
    if any(r_mask)
        scatter(ax_r, xr(r_mask), yr(r_mask), POINT_SIZE, amp(r_mask), 'filled');
    end
    colormap(ax_r,'jet'); clim(ax_r,[-100 100]);
    cb = colorbar(ax_r); cb.Label.String = 'Amplitude (dB)';
    axis(ax_r,'equal'); grid(ax_r,'on');
    xlabel(ax_r,'X (m)'); ylabel(ax_r,'Y (m)');
    title(ax_r, sprintf('Radar 120 GHz — %.0f±%.0f°  (%d pts)', ...
        pd, PITCH_TOL, sum(r_mask)), 'FontSize', 11, 'FontWeight', 'bold');

    % ── LiDAR (bottom row) ───────────────────────────────────────────────
    ax_l = subplot(2, n, n + k);
    plot(ax_l, 0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2); hold(ax_l,'on');
    if any(l_mask)
        scatter(ax_l, xl(l_mask), yl(l_mask), POINT_SIZE, refl(l_mask), 'filled');
    end
    colormap(ax_l,'turbo'); clim(ax_l,[0 255]);
    cb = colorbar(ax_l); cb.Label.String = 'Reflectivity (0–255)';
    axis(ax_l,'equal'); grid(ax_l,'on');
    xlabel(ax_l,'X (m)'); ylabel(ax_l,'Y (m)');
    title(ax_l, sprintf('LiDAR 905 nm — %.0f±%.0f°  (%d pts)', ...
        lidar_el, PITCH_TOL, sum(l_mask)), 'FontSize', 11, 'FontWeight', 'bold');

    fprintf('  Radar %.0f°: %d pts  |  LiDAR %.0f°: %d pts\n', ...
        pd, sum(r_mask), lidar_el, sum(l_mask));
end

sgtitle('Matched angle slices — Radar (top) vs LiDAR (bottom) at −70°, −75°, −80°', ...
    'FontSize', 13, 'FontWeight', 'bold');

fprintf('\n  5 figures generated:\n');
fprintf('    1. Side view XZ — penetration comparison\n');
fprintf('    2. Top-down XY\n');
fprintf('    3. 3D side-by-side\n');
fprintf('    4. GPS world frame Z slice %.0f–%.0fm — LiDAR only\n', Z_LO, Z_HI);
fprintf('    5. Matched angle comparison at -70, -75, -80 deg\n');
fprintf('\nDone.\n');