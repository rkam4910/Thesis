% =========================================================================
% align_lidar_scans.m
% Aligns 4 LiDAR scans (2 per day, 2 days) into a common GPS world frame
% using corner reflectors as tie points.
% Adapted from overlayofplanetree.m (working radar alignment script).
%
% Sunday (15 Mar): file_sun_big, file_sun_small  — GPS: file_gps1
% Monday (16 Mar): file_mon_big, file_mon_small  — GPS: file_gps2
%
% Both scans on the same day share the same LiDAR position, so one
% transform is solved per day and applied to both scans from that day.
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
file_sun_big   = '2026-03-15 14-50-46_corner_big.csv';
file_sun_small = '2026-03-15 14-52-14_corner_small.csv';
file_mon_big   = '2026-03-16_bigreflector_planetree.csv';
file_mon_small = '2026-03-16_smallreflector_planetree.csv';
file_gps1      = 'gps_tree scan 1.csv';
file_gps2      = 'tree scan 2.csv';


REFL_THRESH_SUN = 150;   % reflectivity threshold for Sunday scans
REFL_THRESH_MON = 150;   % reflectivity threshold for Monday scans
RANGE_LIMIT     = 60;    % m — discard far-field returns
FILTER_ZEROS    = true;
SUBSAMPLE_N     = 5;     % keep every Nth point for speed
POINT_SIZE      = 2;

% =========================================================================
fprintf('=== SECTION 1: LOAD LIDAR DATA ===\n');
% =========================================================================

function [x, y, z, refl] = load_lidar(fname, filter_zeros, range_limit, subsample_n)
    T = readtable(fname, 'VariableNamingRule', 'preserve');
    get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};
    x    = get_col('X');
    y    = get_col('Y');
    z    = get_col('Z');
    refl = get_col('Reflectivity');
    if filter_zeros
        v = ~(x==0 & y==0 & z==0);
        x=x(v); y=y(v); z=z(v); refl=refl(v);
    end
    rng  = sqrt(x.^2 + y.^2 + z.^2);
    keep = rng <= range_limit;
    x=x(keep); y=y(keep); z=z(keep); refl=refl(keep);
    idx = 1:subsample_n:numel(x);
    x=x(idx); y=y(idx); z=z(idx); refl=refl(idx);
end

[X_sb, Y_sb, Z_sb, R_sb] = load_lidar(file_sun_big,   FILTER_ZEROS, RANGE_LIMIT, SUBSAMPLE_N);
[X_ss, Y_ss, Z_ss, R_ss] = load_lidar(file_sun_small, FILTER_ZEROS, RANGE_LIMIT, SUBSAMPLE_N);
[X_mb, Y_mb, Z_mb, R_mb] = load_lidar(file_mon_big,   FILTER_ZEROS, RANGE_LIMIT, SUBSAMPLE_N);
[X_ms, Y_ms, Z_ms, R_ms] = load_lidar(file_mon_small, FILTER_ZEROS, RANGE_LIMIT, SUBSAMPLE_N);

fprintf('  Sun big:   %d pts  |  Sun small: %d pts\n', numel(X_sb), numel(X_ss));
fprintf('  Mon big:   %d pts  |  Mon small: %d pts\n', numel(X_mb), numel(X_ms));

% ── Reflectivity diagnostics ─────────────────────────────────────────────
fprintf('\n  Reflectivity ranges:\n');
fprintf('    Sun big:   min=%d  max=%d  mean=%.1f  pts>=200: %d  pts>=150: %d  pts>=100: %d\n', ...
    min(R_sb), max(R_sb), mean(R_sb), sum(R_sb>=200), sum(R_sb>=150), sum(R_sb>=100));
fprintf('    Sun small: min=%d  max=%d  mean=%.1f  pts>=200: %d  pts>=150: %d  pts>=100: %d\n', ...
    min(R_ss), max(R_ss), mean(R_ss), sum(R_ss>=200), sum(R_ss>=150), sum(R_ss>=100));
fprintf('    Mon big:   min=%d  max=%d  mean=%.1f  pts>=200: %d  pts>=150: %d  pts>=100: %d\n', ...
    min(R_mb), max(R_mb), mean(R_mb), sum(R_mb>=200), sum(R_mb>=150), sum(R_mb>=100));
fprintf('    Mon small: min=%d  max=%d  mean=%.1f  pts>=200: %d  pts>=150: %d  pts>=100: %d\n', ...
    min(R_ms), max(R_ms), mean(R_ms), sum(R_ms>=200), sum(R_ms>=150), sum(R_ms>=100));

% =========================================================================
fprintf('\n=== SECTION 2: EXTRACT REFLECTOR CENTROIDS ===\n');
% =========================================================================
% Use Sunday big scan for Sunday tie points (both Sunday scans same position)
% Use Monday big scan for Monday tie points (both Monday scans same position)

function [c_near, c_far] = get_centroids(x, y, z, refl, r_near_expected, r_far_expected)
    % Find two reflector centroids using tight range windows around GPS distances
    % Take only the top highest-reflectivity points within each window
    WINDOW  = 1.0;   % ± metres — tight window around expected GPS range
    TOP_PTS = 20;    % use only top N highest-reflectivity points per window

    rng = sqrt(x.^2 + y.^2);

    near_mask = abs(rng - r_near_expected) <= WINDOW;
    far_mask  = abs(rng - r_far_expected)  <= WINDOW;

    fprintf('    Near window %.1f±%.1fm: %d pts\n', r_near_expected, WINDOW, sum(near_mask));
    fprintf('    Far  window %.1f±%.1fm: %d pts\n', r_far_expected,  WINDOW, sum(far_mask));

    if sum(near_mask) < 1 || sum(far_mask) < 1
        error('No points found in reflector range window. Check expected ranges or increase WINDOW.');
    end

    function c = top_centroid(xw, yw, zw, rw)
        n = min(TOP_PTS, numel(xw));
        [~, si] = sort(rw, 'descend');
        si = si(1:n);
        c = [mean(xw(si)), mean(yw(si)), mean(zw(si))];
    end

    c_near = top_centroid(x(near_mask), y(near_mask), z(near_mask), refl(near_mask));
    c_far  = top_centroid(x(far_mask),  y(far_mask),  z(far_mask),  refl(far_mask));

    fprintf('    Near centroid: x=%.3f  y=%.3f  z=%.3f  (range=%.2fm)\n', c_near, norm(c_near(1:2)));
    fprintf('    Far  centroid: x=%.3f  y=%.3f  z=%.3f  (range=%.2fm)\n', c_far,  norm(c_far(1:2)));
end

% =========================================================================
fprintf('\n=== SECTION 2b: DIAGNOSTIC — raw reflector clusters (sensor frame) ===\n');
% =========================================================================
scan_names = {'Sun big', 'Sun small', 'Mon big', 'Mon small'};
Xs_all = {X_sb, X_ss, X_mb, X_ms};
Ys_all = {Y_sb, Y_ss, Y_mb, Y_ms};
Rs_all = {R_sb, R_ss, R_mb, R_ms};
thresh_all = [REFL_THRESH_SUN, REFL_THRESH_SUN, REFL_THRESH_MON, REFL_THRESH_MON];

figure('Name','Diagnostic — high-reflectivity clusters (sensor frame)', ...
    'Color','k','Units','normalized','Position',[0.02 0.05 0.96 0.85]);

for k = 1:4
    x = Xs_all{k}; y = Ys_all{k}; r = Rs_all{k};
    TOP_N     = 500;
    RANGE_CAP = 30;
    rng_all   = sqrt(x.^2 + y.^2);
    in_range  = rng_all <= RANGE_CAP;
    x = x(in_range); y = y(in_range); r = r(in_range);
    [~, si] = sort(r, 'descend');
    idx    = si(1:min(TOP_N, numel(x)));
    x_hi   = x(idx); y_hi = y(idx);
    rng_hi = sqrt(x_hi.^2 + y_hi.^2);
    [r_sort, si] = sort(rng_hi);
    [~, gi] = max(diff(r_sort));
    near_mask = false(size(x_hi)); near_mask(si(1:gi))    = true;
    far_mask  = false(size(x_hi)); far_mask(si(gi+1:end)) = true;

    % XY view
    ax = subplot(2, 4, k);
    set(ax,'Color','w','XColor','w','YColor','w');
    scatter(ax, x_hi(near_mask), y_hi(near_mask), 15, 'c', 'filled', 'DisplayName','Near');
    hold(ax,'on');
    scatter(ax, x_hi(far_mask),  y_hi(far_mask),  15, 'm', 'filled', 'DisplayName','Far');
    plot(ax, 0, 0, 'wx', 'MarkerSize', 12, 'LineWidth', 2);
    axis(ax,'equal'); grid(ax,'on'); ax.GridColor=[0.3 0.3 0.3];
    xlabel(ax,'X (m)','Color','w'); ylabel(ax,'Y (m)','Color','w');
    title(ax, sprintf('%s XY — near %.1fm (%d pts)  far %.1fm (%d pts)', ...
        scan_names{k}, mean(r_sort(1:gi)), sum(near_mask), ...
        mean(r_sort(gi+1:end)), sum(far_mask)), 'Color','w','FontSize',8);
    legend(ax,'show','TextColor','w','Color','k','FontSize',7);

    % Range histogram
    ax2 = subplot(2, 4, k+4);
    set(ax2,'Color','w','XColor','w','YColor','w');
    histogram(ax2, rng_hi, 40, 'FaceColor',[0.3 0.6 1], 'EdgeColor','none');
    xline(ax2, r_sort(gi), 'r--', 'LineWidth', 2);
    xlabel(ax2,'Range (m)','Color','w'); ylabel(ax2,'Count','Color','w');
    title(ax2, sprintf('%s range histogram — gap at %.1fm', scan_names{k}, r_sort(gi)), ...
        'Color','w','FontSize',8);
    grid(ax2,'on'); ax2.GridColor=[0.3 0.3 0.3];
end
sgtitle('High-reflectivity clusters — cyan=near  magenta=far  red=gap', ...
    'Color','w','FontSize',12);
fprintf('  Near/far ranges should match GPS horiz distances to reflectors\n');

fprintf('\n  Sunday (using big scan):\n');
% Sun GPS: small=9.96m, big=13.03m
[c_sun_near, c_sun_far] = get_centroids(X_sb, Y_sb, Z_sb, R_sb, 9.96, 13.03);

fprintf('  Monday (using big scan):\n');
% Mon GPS distances from Mon LiDAR: small=9.49m, big=13.63m
[c_mon_near, c_mon_far] = get_centroids(X_mb, Y_mb, Z_mb, R_mb, 9.49, 13.63);

% =========================================================================
fprintf('\n=== SECTION 3: LOAD GPS AND CONVERT TO ENU ===\n');
% =========================================================================
opts1 = detectImportOptions(file_gps1, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
G1    = readtable(file_gps1, opts1);
opts2 = detectImportOptions(file_gps2, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
G2    = readtable(file_gps2, opts2);

G1.Properties.VariableNames{1} = 'PointName';
G2.Properties.VariableNames{1} = 'PointName';
G1.PointName = strtrim(G1.PointName);
G2.PointName = strtrim(G2.PointName);

fprintf('  G1 point names: %s\n', strjoin(G1.PointName, ' | '));
fprintf('  G2 point names: %s\n', strjoin(G2.PointName, ' | '));

% World frame origin = Sunday LiDAR GPS position
idx  = strcmp(G1.PointName, 'radar/lidar');
olon = G1.Longitude(idx);
olat = G1.Latitude(idx);
oh   = G1.("Ellipsoidal height")(idx);
fprintf('  World origin: lon=%.8f  lat=%.8f  h=%.3fm\n', olon, olat, oh);

gps2enu = @(lon, lat, h) deal( ...
    (lon - olon) .* 111320 .* cos(deg2rad(olat)), ...
    (lat - olat) .* 110540, ...
     h - oh);

% Sunday GPS reflectors
[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector big')));
wf1_big = [e,n,u];

[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector small')));
wf1_small = [e,n,u];

% Monday GPS reflectors — use big1 (better reading)
[e,n,u] = gps2enu(G2.Longitude(strcmp(G2.PointName,'corner refelctor big1')), ...
                  G2.Latitude(strcmp(G2.PointName,'corner refelctor big1')), ...
                  G2.("Ellipsoidal height")(strcmp(G2.PointName,'corner refelctor big1')));
wf2_big = [e,n,u];

[e,n,u] = gps2enu(G2.Longitude(strcmp(G2.PointName,'corner refelctor small')), ...
                  G2.Latitude(strcmp(G2.PointName,'corner refelctor small')), ...
                  G2.("Ellipsoidal height")(strcmp(G2.PointName,'corner refelctor small')));
wf2_small = [e,n,u];

[e,n,u] = gps2enu(G2.Longitude(strcmp(G2.PointName,'Radar/Lidar')), ...
                  G2.Latitude(strcmp(G2.PointName,'Radar/Lidar')), ...
                  G2.("Ellipsoidal height")(strcmp(G2.PointName,'Radar/Lidar')));
wf2_lidar = [e,n,u];

fprintf('  Sun big   ENU: E=%7.3f  N=%7.3f  U=%7.3f  (horiz=%.2fm)\n', wf1_big,   norm(wf1_big(1:2)));
fprintf('  Sun small ENU: E=%7.3f  N=%7.3f  U=%7.3f  (horiz=%.2fm)\n', wf1_small, norm(wf1_small(1:2)));
fprintf('  Mon big   ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_big);
fprintf('  Mon small ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_small);
fprintf('  Mon LiDAR ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_lidar);

% =========================================================================
fprintf('\n=== SECTION 4: MATCH CENTROIDS TO GPS ===\n');
% =========================================================================

% Sunday — near window centred on small GPS (9.96m), far on big GPS (13.03m)
s1_near_wf = wf1_small;
s1_far_wf  = wf1_big;
fprintf('  Sun match: near->small (%.2fm)  far->big (%.2fm)\n', ...
    norm(wf1_small(1:2)), norm(wf1_big(1:2)));

% Monday — near window centred on small GPS (9.49m), far on big GPS (13.63m)
s2_near_wf = wf2_small;
s2_far_wf  = wf2_big;
fprintf('  Mon match: near->small (%.2fm)  far->big (%.2fm)\n', ...
    norm(wf2_small(1:2)-wf2_lidar(1:2)), norm(wf2_big(1:2)-wf2_lidar(1:2)));

% Tie points — 3 per day: near reflector, far reflector, LiDAR origin
% Sunday: origin maps to world origin (0,0) — ENU origin IS the Sunday LiDAR
% Monday: origin maps to wf2_lidar (Monday LiDAR GPS position in ENU)
sensor_pts_sun = [c_sun_near(1:2); c_sun_far(1:2); 0, 0];
world_pts_sun  = [s1_near_wf(1:2); s1_far_wf(1:2); 0, 0];

sensor_pts_mon = [c_mon_near(1:2); c_mon_far(1:2); 0, 0];
world_pts_mon  = [s2_near_wf(1:2); s2_far_wf(1:2); wf2_lidar(1:2)];

% =========================================================================
fprintf('\n=== SECTION 5: SOLVE 2D TRANSFORMS (SVD) ===\n');
% =========================================================================
function [R2d, t2d] = solve_2d_transform(sp, wp)
    % Solves 2D rigid body transform: p_world = R * p_sensor + t
    % Using SVD Procrustes method — more robust than atan2 single-vector approach.
    % det(V*U') term prevents reflection (mirror image) solutions.
    mu_s = mean(sp, 1);
    mu_w = mean(wp, 1);
    S    = (sp - mu_s)' * (wp - mu_w);
    [U, ~, V] = svd(S);
    R2d = V * diag([1, det(V*U')]) * U';
    t2d = (mu_w' - R2d * mu_s');   % column vector
end

[R_sun, t_sun] = solve_2d_transform(sensor_pts_sun, world_pts_sun);
[R_mon, t_mon] = solve_2d_transform(sensor_pts_mon, world_pts_mon);

theta_sun = atan2(R_sun(2,1), R_sun(1,1));
theta_mon = atan2(R_mon(2,1), R_mon(1,1));

fprintf('  Sunday: yaw=%.2f deg  t=[%.3f, %.3f] m\n', rad2deg(theta_sun), t_sun);
fprintf('  Monday: yaw=%.2f deg  t=[%.3f, %.3f] m\n', rad2deg(theta_mon), t_mon);

% =========================================================================
fprintf('\n=== SECTION 6: APPLY TRANSFORMS ===\n');
% =========================================================================
cs = cos(theta_sun); ss = sin(theta_sun);
cm = cos(theta_mon); sm = sin(theta_mon);

% Transposed rotation matrices for row-vector convention: p * R^T
Rsun = [cs ss; -ss cs];
Rmon = [cm sm; -sm cm];

xy_sb = [X_sb, Y_sb] * Rsun + t_sun';
xy_ss = [X_ss, Y_ss] * Rsun + t_sun';
xy_mb = [X_mb, Y_mb] * Rmon + t_mon';
xy_ms = [X_ms, Y_ms] * Rmon + t_mon';

% Z: no offset — heights unknown, each scan relative to its own sensor
% Apply Z height offset for Monday scans
% Monday LiDAR was 6.137m higher than Sunday based on GPS
z_offset_mon = wf2_lidar(3);   % = 6.137 m
z_offset_empirical = 0;        % adjust this if ground planes still don't align


Wsb = [xy_sb, Z_sb];
Wss = [xy_ss, Z_ss];
Wmb = [xy_mb, Z_mb];
Wms = [xy_ms, Z_ms];

all_e = [Wsb(:,1); Wss(:,1); Wmb(:,1); Wms(:,1)];
all_n = [Wsb(:,2); Wss(:,2); Wmb(:,2); Wms(:,2)];
fprintf('  World E range: %.2f to %.2f m\n', min(all_e), max(all_e));
fprintf('  World N range: %.2f to %.2f m\n', min(all_n), max(all_n));

% =========================================================================
fprintf('\n=== SECTION 7: PLOT ===\n');
% =========================================================================
figure('Name','LiDAR 4-scan overlay — GPS world frame', ...
    'Color','w','Units','normalized','Position',[0.02 0.05 0.7 0.85]);
ax = axes('Color','w','XColor','k','YColor','k','ZColor','k');
hold(ax,'on');

scatter3(ax, Wsb(:,1), Wsb(:,2), Wsb(:,3), POINT_SIZE, R_sb, 'o','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName','Sun big');
scatter3(ax, Wss(:,1), Wss(:,2), Wss(:,3), POINT_SIZE, R_ss, 'o','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName','Sun small');
scatter3(ax, Wmb(:,1), Wmb(:,2), Wmb(:,3), POINT_SIZE, R_mb, 's','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName','Mon big');
scatter3(ax, Wms(:,1), Wms(:,2), Wms(:,3), POINT_SIZE, R_ms, 's','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName','Mon small');

plot3(ax, wf1_big(1),   wf1_big(2),   wf1_big(3),   'r^','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Sun big refl GPS');
plot3(ax, wf1_small(1), wf1_small(2), wf1_small(3), 'g^','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Sun small refl GPS');
plot3(ax, wf2_big(1),   wf2_big(2),   wf2_big(3),   'rv','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Mon big refl GPS');
plot3(ax, wf2_small(1), wf2_small(2), wf2_small(3), 'gv','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Mon small refl GPS');
plot3(ax, wf2_lidar(1), wf2_lidar(2), wf2_lidar(3), 'bs','MarkerSize',12,'MarkerFaceColor','b','DisplayName','Mon LiDAR origin GPS');
plot3(ax, 0, 0, 0, 'kx','MarkerSize',14,'LineWidth',2,'DisplayName','Sun LiDAR origin (world origin)');

colormap(ax,'turbo');
cb = colorbar(ax); cb.Color='k'; cb.Label.String='Reflectivity (0-255)'; cb.Label.Color='k';
xlabel(ax,'East (m)','Color','k'); ylabel(ax,'North (m)','Color','k');
zlabel(ax,'Z (m)','Color','k');
title(ax,'4 LiDAR scans — GPS world frame (circles=Sunday, squares=Monday)', ...
    'Color','k','FontSize',12);
grid(ax,'on'); ax.GridColor=[0.8 0.8 0.8];
axis(ax,'equal'); view(ax,45,25); rotate3d(ax,'on');
legend(ax,'show','TextColor','k','Color','w','FontSize',8,'Location','northeast');

fprintf('  Plot complete.\n\nDone.\n');

% =========================================================================
%% LIDAR PENETRATION DEPTH — normalised density vs range
% =========================================================================

HALF_ANGLE_LIDAR = 15;   % degrees
BIN_WIDTH_LIDAR  = 0.5;  % m
tree_centre_world = [-11.9, 32.0, 9.0];

% ---- SUNDAY CONE (from world origin) ----
bs_sun = tree_centre_world(1:2) / norm(tree_centre_world(1:2));

pts_sun = [Wsb(:,1:2); Wss(:,1:2)];
z_sun   = [Wsb(:,3);   Wss(:,3)];
norms_sun = sqrt(sum(pts_sun.^2, 2));
pts_sun_unit = pts_sun ./ norms_sun;
angle_sun = acosd(pts_sun_unit * bs_sun');
in_cone_sun = angle_sun <= HALF_ANGLE_LIDAR & z_sun >= 0.5;
range_sun_cone = norms_sun(in_cone_sun);

fprintf('Sunday LiDAR points in cone: %d\n', sum(in_cone_sun));

% ---- MONDAY CONE (from monday lidar origin) ----
bs_mon_vec = tree_centre_world(1:2) - wf2_lidar(1:2);
bs_mon = bs_mon_vec / norm(bs_mon_vec);

pts_mon = [Wmb(:,1:2); Wms(:,1:2)];
z_mon   = [Wmb(:,3);   Wms(:,3)];
pts_mon_rel = pts_mon - wf2_lidar(1:2);
norms_mon = sqrt(sum(pts_mon_rel.^2, 2));
pts_mon_unit = pts_mon_rel ./ norms_mon;
angle_mon = acosd(pts_mon_unit * bs_mon');
in_cone_mon = angle_mon <= HALF_ANGLE_LIDAR & z_mon >= 0.5;
range_mon_cone = norms_mon(in_cone_mon);

fprintf('Monday LiDAR points in cone: %d\n', sum(in_cone_mon));

% ---- BIN AND PLOT ----
r_edges   = 0 : BIN_WIDTH_LIDAR : 60;
r_centres = r_edges(1:end-1) + BIN_WIDTH_LIDAR/2;

counts_sun = histcounts(range_sun_cone, r_edges);
counts_mon = histcounts(range_mon_cone, r_edges);
counts_merged = counts_sun + counts_mon;
counts_norm = counts_merged / max(counts_merged);

lidar_r_centres = r_centres;
lidar_counts_norm = counts_norm;
save('..\lidar_cone_density.mat', 'lidar_r_centres', 'lidar_counts_norm');
fprintf('Saved: lidar_cone_density.mat\n');

figure('Color','w','Units','normalized','Position',[0.1 0.3 0.6 0.45]);

plot(r_centres, counts_norm, '-', ...
    'Color', [0.2 0.6 1.0], 'LineWidth', 2.0, ...
    'DisplayName', 'LiDAR 905 nm (both days)');

xline(14, '--k', 'LineWidth', 1.2, 'Label', 'Canopy front', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');
xline(40, '--k', 'LineWidth', 1.2, 'Label', 'Canopy back', ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');

xlabel('Range from sensor (m)', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Normalised point density', 'FontSize', 12, 'Interpreter', 'none');
title('LiDAR point density through plane tree canopy', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
xlim([5 55]);
ylim([0 1.05]);

exportgraphics(gcf, 'fig_lidar_cone_density.png', 'Resolution', 300);
fprintf('Saved: fig_lidar_cone_density.png\n');