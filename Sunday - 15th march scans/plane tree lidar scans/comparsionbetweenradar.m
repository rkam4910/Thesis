% =========================================================================
% plane_tree_analysis.m
%
% Three complementary comparisons of radar vs LiDAR on the plane tree:
%   1. Height profile (density vs Z)
%   2. Amplitude/reflectivity distribution
%   3. Side view point cloud (XZ)
%   4. Cone penetration comparison
%   5. Top down view (XY)
%
% Radar:  pointcloud_20260315-143039_.csv
% LiDAR:  2026-03-15 14-52-14_corner_small.csv
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
file_radar = 'pointcloud_20260315-143039_.csv';
file_lidar = '2026-03-15 14-52-14_corner_small.csv';

RANGE_LIMIT = 50;    % m
AMP_THRESH  = -80;   % dB
POINT_SIZE  = 4;
BIN_WIDTH   = 0.5;   % m

X_MIN = 8;
X_MAX = 40;
Z_MIN = 0.5;         % m — remove ground returns below this

HALF_ANGLE  = 15;    % degrees — radar cone half angle

% Radar tree centre for boresight
tree_centre_r = [22, 18, 8];

% -------------------------------------------------------------------------
% LOAD RADAR
% -------------------------------------------------------------------------
fprintf('=== LOAD RADAR ===\n');
data = readmatrix(file_radar);
data = data(data(:,5) >= AMP_THRESH, :);
data = data(data(:,2) <= RANGE_LIMIT, :);

az  = data(:,4);
el  = (pi/2) + data(:,3);
r   = data(:,2);
amp = data(:,5);

[xr, yr, zr] = sph2cart(az, el, r);

in_tree = xr >= X_MIN & xr <= X_MAX & zr >= Z_MIN;
xr  = xr(in_tree);
yr  = yr(in_tree);
zr  = zr(in_tree);
amp = amp(in_tree);

fprintf('  %d radar points after filtering\n', numel(xr));

% -------------------------------------------------------------------------
% LOAD LIDAR
% -------------------------------------------------------------------------
fprintf('\n=== LOAD LIDAR ===\n');
T = readtable(file_lidar, 'VariableNamingRule', 'preserve');
get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};

xl   = get_col('X');
yl   = get_col('Y');
zl   = get_col('Z');
refl = get_col('Reflectivity');

valid = ~(xl == 0 & yl == 0 & zl == 0);
xl = xl(valid); yl = yl(valid); zl = zl(valid); refl = refl(valid);

range_l = sqrt(xl.^2 + yl.^2 + zl.^2);
keep    = range_l <= RANGE_LIMIT;
xl = xl(keep); yl = yl(keep); zl = zl(keep); refl = refl(keep);

in_tree_l = xl >= X_MIN & xl <= X_MAX & zl >= Z_MIN;
xl   = xl(in_tree_l);
yl   = yl(in_tree_l);
zl   = zl(in_tree_l);
refl = refl(in_tree_l);

fprintf('  %d LiDAR points after filtering\n', numel(xl));

% =========================================================================
%% FIGURE 1: Height profile — density vs Z
% =========================================================================
z_edges  = Z_MIN : BIN_WIDTH : 20;
z_centres = z_edges(1:end-1) + BIN_WIDTH/2;

counts_r_z = histcounts(zr, z_edges);
counts_l_z = histcounts(zl, z_edges);

counts_radar_norm_z = counts_r_z / max(counts_r_z);
counts_lidar_norm_z = counts_l_z / max(counts_l_z);

figure('Color','w','Units','normalized','Position',[0.1 0.3 0.5 0.5]);

plot(counts_lidar_norm_z, z_centres, '-', ...
    'Color', [0.2 0.6 1.0], 'LineWidth', 2.0, ...
    'DisplayName', 'LiDAR 905 nm');
hold on;
plot(counts_radar_norm_z, z_centres, '-', ...
    'Color', [0.9 0.3 0.1], 'LineWidth', 2.0, ...
    'DisplayName', 'Radar 120 GHz');

xlabel('Normalised point density', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Height above ground - Z (m)', 'FontSize', 12, 'Interpreter', 'none');
title('Height profile - radar vs LiDAR point density', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
ylim([Z_MIN 20]);

exportgraphics(gcf, 'fig1_height_profile.png', 'Resolution', 300);
fprintf('Saved: fig1_height_profile.png\n');

% =========================================================================
%% FIGURE 3: Side view XZ — point cloud comparison
% =========================================================================
figure('Color','w','Units','normalized','Position',[0.02 0.05 0.96 0.55]);

ax1 = subplot(1,2,1);
scatter(ax1, xr, zr, POINT_SIZE, amp, 'filled');
colormap(ax1, 'jet');
clim(ax1, [-80 70]);
cb1 = colorbar(ax1);
cb1.Label.String = 'Amplitude (dB)';
cb1.Label.Interpreter = 'none';
axis(ax1, 'equal');
grid(ax1, 'on');
xlabel(ax1, 'X - depth (m)', 'Interpreter', 'none');
ylabel(ax1, 'Z - height (m)', 'Interpreter', 'none');
title(ax1, 'Radar 120 GHz - side view', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

ax2 = subplot(1,2,2);
scatter(ax2, xl, zl, POINT_SIZE, refl, 'filled');
colormap(ax2, 'turbo');
clim(ax2, [0 255]);
cb2 = colorbar(ax2);
cb2.Label.String = 'Reflectivity (0-255)';
cb2.Label.Interpreter = 'none';
axis(ax2, 'equal');
grid(ax2, 'on');
xlabel(ax2, 'X - depth (m)', 'Interpreter', 'none');
ylabel(ax2, 'Z - height (m)', 'Interpreter', 'none');
title(ax2, 'LiDAR 905 nm - side view', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

exportgraphics(gcf, 'fig3_side_view.png', 'Resolution', 300);
fprintf('Saved: fig3_side_view.png\n');

% =========================================================================
%% FIGURE 4: Top down view XY
% =========================================================================
figure('Color','w','Units','normalized','Position',[0.1 0.1 0.7 0.7]);

ax_top1 = subplot(1,2,1);
scatter(ax_top1, xr, yr, POINT_SIZE, amp, 'filled');
colormap(ax_top1, 'jet');
clim(ax_top1, [-80 70]);
cb3 = colorbar(ax_top1);
cb3.Label.String = 'Amplitude (dB)';
cb3.Label.Interpreter = 'none';
axis(ax_top1, 'equal');
grid(ax_top1, 'on');
xlabel(ax_top1, 'X (m)', 'Interpreter', 'none');
ylabel(ax_top1, 'Y (m)', 'Interpreter', 'none');
title(ax_top1, 'Radar - top down view (XY)', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

ax_top2 = subplot(1,2,2);
scatter(ax_top2, xl, yl, POINT_SIZE, refl, 'filled');
colormap(ax_top2, 'turbo');
clim(ax_top2, [0 255]);
cb4 = colorbar(ax_top2);
cb4.Label.String = 'Reflectivity (0-255)';
cb4.Label.Interpreter = 'none';
axis(ax_top2, 'equal');
grid(ax_top2, 'on');
xlabel(ax_top2, 'X (m)', 'Interpreter', 'none');
ylabel(ax_top2, 'Y (m)', 'Interpreter', 'none');
title(ax_top2, 'LiDAR - top down view (XY)', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

exportgraphics(gcf, 'fig4_top_down.png', 'Resolution', 300);
fprintf('Saved: fig4_top_down.png\n');

% =========================================================================
%% FIGURE 5: Cone penetration comparison
% =========================================================================

% ---- RADAR CONE ----
boresight_r = tree_centre_r / norm(tree_centre_r);

pts_r  = [xr, yr, zr];
norms_r = sqrt(sum(pts_r.^2, 2));
pts_r_unit = pts_r ./ norms_r;

cos_angle_r = pts_r_unit * boresight_r';
angle_r     = acosd(cos_angle_r);

in_cone_r    = angle_r <= HALF_ANGLE;
range_r_cone = norms_r(in_cone_r);

fprintf('Radar points in cone: %d\n', sum(in_cone_r));

% ---- LIDAR — all points within canopy range and height ----
% LiDAR uses rosette pattern covering full FOV, no cone extraction needed
range_l_all = sqrt(xl.^2 + yl.^2 + zl.^2);

in_canopy_l  = range_l_all >= 10 & range_l_all <= 45 & zl >= Z_MIN;
range_l_cone = range_l_all(in_canopy_l);

fprintf('LiDAR points in canopy volume: %d\n', sum(in_canopy_l));

% ---- BIN BY RANGE ----
r_edges   = 0 : BIN_WIDTH : 50;
r_centres = r_edges(1:end-1) + BIN_WIDTH/2;

counts_cone_r = histcounts(range_r_cone, r_edges);
counts_cone_l = histcounts(range_l_cone, r_edges);

counts_cone_r_norm = counts_cone_r / max(counts_cone_r);
counts_cone_l_norm = counts_cone_l / max(counts_cone_l);

% ---- PLOT ----
figure('Color','w','Units','normalized','Position',[0.1 0.3 0.6 0.45]);

plot(r_centres, counts_cone_l_norm, '-', ...
    'Color', [0.2 0.6 1.0], 'LineWidth', 2.0, ...
    'DisplayName', 'LiDAR 905 nm');
hold on;
plot(r_centres, counts_cone_r_norm, '-', ...
    'Color', [0.9 0.3 0.1], 'LineWidth', 2.0, ...
    'DisplayName', 'Radar 120 GHz');

xline(25, '--k', 'LineWidth', 1.2, 'Label', 'Canopy front', ...
    'LabelHorizontalAlignment', 'right', 'Interpreter', 'none', ...
    'HandleVisibility', 'off');
xline(37, '--k', 'LineWidth', 1.2, 'Label', 'Canopy back', ...
    'LabelHorizontalAlignment', 'left', 'Interpreter', 'none', ...
    'HandleVisibility', 'off');

xlabel('Range from sensor (m)', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Normalised point density', 'FontSize', 12, 'Interpreter', 'none');
title('Radar vs LiDAR point density through plane tree canopy', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
xlim([5 45]);
ylim([0 1.05]);

exportgraphics(gcf, 'fig5_cone_penetration.png', 'Resolution', 300);
fprintf('Saved: fig5_cone_penetration.png\n');

fprintf('\nAll figures saved.\n');