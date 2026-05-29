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
%
% Section 8 adds conical attenuation analysis in the GPS world frame,
% allowing Sunday and Monday cones to be directly compared since all
% clouds are in the same coordinate system.
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

REFL_THRESH_SUN = 150;
REFL_THRESH_MON = 150;
RANGE_LIMIT     = 60;
FILTER_ZEROS    = true;
SUBSAMPLE_N     = 5;
POINT_SIZE      = 2;

% Z offset: Monday LiDAR was 6.137m higher than Sunday LiDAR
DZ_MON_TO_SUN = 50.105 - 43.968;   % = 6.137 m

% ── Attenuation analysis settings (Section 8) ────────────────────────────
% Tree bounding box in ENU world frame — adjust after viewing Section 7 plot
ATT_TREE_E_MIN = -10;  ATT_TREE_E_MAX =  10;
ATT_TREE_N_MIN =   5;  ATT_TREE_N_MAX =  25;
ATT_TREE_Z_MIN =   0;  ATT_TREE_Z_MAX =  20;

CONE_HALF_ANGLE_DEG = 25;   % half-angle of conical section
BIN_START  =  2.0;           % m — start of first distance bin
BIN_END    = 25.0;           % m — end of last bin
BIN_WIDTH  =  1.0;           % m — bin width

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

function [c_near, c_far] = get_centroids(x, y, z, refl, r_near_expected, r_far_expected)
    WINDOW  = 1.0;
    TOP_PTS = 20;
    rng = sqrt(x.^2 + y.^2);
    near_mask = abs(rng - r_near_expected) <= WINDOW;
    far_mask  = abs(rng - r_far_expected)  <= WINDOW;
    fprintf('    Near window %.1f±%.1fm: %d pts\n', r_near_expected, WINDOW, sum(near_mask));
    fprintf('    Far  window %.1f±%.1fm: %d pts\n', r_far_expected,  WINDOW, sum(far_mask));
    if sum(near_mask) < 1 || sum(far_mask) < 1
        error('No points found in reflector range window.');
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
    TOP_N = 500; RANGE_CAP = 30;
    rng_all  = sqrt(x.^2 + y.^2);
    in_range = rng_all <= RANGE_CAP;
    x = x(in_range); y = y(in_range); r = r(in_range);
    [~, si] = sort(r, 'descend');
    idx    = si(1:min(TOP_N, numel(x)));
    x_hi   = x(idx); y_hi = y(idx);
    rng_hi = sqrt(x_hi.^2 + y_hi.^2);
    [r_sort, si] = sort(rng_hi);
    [~, gi] = max(diff(r_sort));
    near_mask = false(size(x_hi)); near_mask(si(1:gi))    = true;
    far_mask  = false(size(x_hi)); far_mask(si(gi+1:end)) = true;

    ax = subplot(2, 4, k);
    set(ax,'Color','k','XColor','w','YColor','w');
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

    ax2 = subplot(2, 4, k+4);
    set(ax2,'Color','k','XColor','w','YColor','w');
    histogram(ax2, rng_hi, 40, 'FaceColor',[0.3 0.6 1], 'EdgeColor','none');
    xline(ax2, r_sort(gi), 'r--', 'LineWidth', 2);
    xlabel(ax2,'Range (m)','Color','w'); ylabel(ax2,'Count','Color','w');
    title(ax2, sprintf('%s range histogram — gap at %.1fm', scan_names{k}, r_sort(gi)), ...
        'Color','w','FontSize',8);
    grid(ax2,'on'); ax2.GridColor=[0.3 0.3 0.3];
end
sgtitle('High-reflectivity clusters — cyan=near  magenta=far  red=gap', ...
    'Color','w','FontSize',12);

fprintf('\n  Sunday (using big scan):\n');
[c_sun_near, c_sun_far] = get_centroids(X_sb, Y_sb, Z_sb, R_sb, 9.96, 13.03);

fprintf('  Monday (using big scan):\n');
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

idx  = strcmp(G1.PointName, 'radar/lidar');
olon = G1.Longitude(idx);
olat = G1.Latitude(idx);
oh   = G1.("Ellipsoidal height")(idx);
fprintf('  World origin: lon=%.8f  lat=%.8f  h=%.3fm\n', olon, olat, oh);

gps2enu = @(lon, lat, h) deal( ...
    (lon - olon) .* 111320 .* cos(deg2rad(olat)), ...
    (lat - olat) .* 110540, ...
     h - oh);

[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector big')));
wf1_big = [e,n,u];

[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector small')));
wf1_small = [e,n,u];

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
s1_near_wf = wf1_small;
s1_far_wf  = wf1_big;
fprintf('  Sun match: near->small (%.2fm)  far->big (%.2fm)\n', ...
    norm(wf1_small(1:2)), norm(wf1_big(1:2)));

s2_near_wf = wf2_small;
s2_far_wf  = wf2_big;
fprintf('  Mon match: near->small (%.2fm)  far->big (%.2fm)\n', ...
    norm(wf2_small(1:2)-wf2_lidar(1:2)), norm(wf2_big(1:2)-wf2_lidar(1:2)));

sensor_pts_sun = [c_sun_near(1:2); c_sun_far(1:2); 0, 0];
world_pts_sun  = [s1_near_wf(1:2); s1_far_wf(1:2); 0, 0];

sensor_pts_mon = [c_mon_near(1:2); c_mon_far(1:2); 0, 0];
world_pts_mon  = [s2_near_wf(1:2); s2_far_wf(1:2); wf2_lidar(1:2)];

% =========================================================================
fprintf('\n=== SECTION 5: SOLVE 2D TRANSFORMS (SVD) ===\n');
% =========================================================================
function [R2d, t2d] = solve_2d_transform(sp, wp)
    mu_s = mean(sp, 1);
    mu_w = mean(wp, 1);
    S    = (sp - mu_s)' * (wp - mu_w);
    [U, ~, V] = svd(S);
    R2d = V * diag([1, det(V*U')]) * U';
    t2d = (mu_w' - R2d * mu_s');
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

Rsun = [cs ss; -ss cs];
Rmon = [cm sm; -sm cm];

xy_sb = [X_sb, Y_sb] * Rsun + t_sun';
xy_ss = [X_ss, Y_ss] * Rsun + t_sun';
xy_mb = [X_mb, Y_mb] * Rmon + t_mon';
xy_ms = [X_ms, Y_ms] * Rmon + t_mon';

% Z: subtract height offset from Monday clouds to align to Sunday datum
Wsb = [xy_sb, Z_sb];
Wss = [xy_ss, Z_ss];
Wmb = [xy_mb, Z_mb - DZ_MON_TO_SUN];
Wms = [xy_ms, Z_ms - DZ_MON_TO_SUN];

all_e = [Wsb(:,1); Wss(:,1); Wmb(:,1); Wms(:,1)];
all_n = [Wsb(:,2); Wss(:,2); Wmb(:,2); Wms(:,2)];
all_z = [Wsb(:,3); Wss(:,3); Wmb(:,3); Wms(:,3)];
fprintf('  World E range: %.2f to %.2f m\n', min(all_e), max(all_e));
fprintf('  World N range: %.2f to %.2f m\n', min(all_n), max(all_n));
fprintf('  World Z range: %.2f to %.2f m\n', min(all_z), max(all_z));

% =========================================================================
fprintf('\n=== SECTION 7: ALIGNMENT PLOT ===\n');
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
zlabel(ax,'Z relative to Sunday sensor (m)','Color','k');
title(ax,sprintf('4 LiDAR scans — GPS world frame\nZ offset Mon->Sun: -%.3fm applied', DZ_MON_TO_SUN), ...
    'Color','k','FontSize',12);
grid(ax,'on'); ax.GridColor=[0.8 0.8 0.8];
axis(ax,'equal'); view(ax,45,25); rotate3d(ax,'on');
legend(ax,'show','TextColor','k','Color','w','FontSize',8,'Location','northeast');

fprintf('  Alignment plot complete.\n');

% =========================================================================
fprintf('\n=== SECTION 8: CONICAL ATTENUATION ANALYSIS (world frame) ===\n');
% =========================================================================
% Now that all four clouds are in the GPS ENU world frame, the cone axes
% are defined by vectors from each LiDAR origin to each reflector GPS
% position — all in the same coordinate system.
%
% Sunday LiDAR origin in world frame = [0, 0, 0]
% Monday LiDAR origin in world frame = wf2_lidar (ENU)
%
% Four cones are defined:
%   1. Sun->small reflector  (from [0,0,0] toward wf1_small)
%   2. Sun->big reflector    (from [0,0,0] toward wf1_big)
%   3. Mon->small reflector  (from wf2_lidar toward wf2_small)
%   4. Mon->big reflector    (from wf2_lidar toward wf2_big)
%
% The pooled world-frame cloud is searched for each cone.
% Distance bins are in slant range from the relevant LiDAR origin.
% =========================================================================

% Pool all four aligned clouds (full arrays — Wsb etc. are subsampled,
% so reload using the transform directly on the full sensor arrays)
% For memory efficiency we reuse the already-subsampled Wsb/Wss/Wmb/Wms
% which are 1-in-5 subsamples — sufficient for density/reflectivity stats.
%
% Build unified pool with scan identity labels
pool_e    = [Wsb(:,1); Wss(:,1); Wmb(:,1); Wms(:,1)];
pool_n    = [Wsb(:,2); Wss(:,2); Wmb(:,2); Wms(:,2)];
pool_z    = [Wsb(:,3); Wss(:,3); Wmb(:,3); Wms(:,3)];
pool_refl = [R_sb;     R_ss;     R_mb;     R_ms];
pool_day  = [ones(size(R_sb)); ones(size(R_ss)); ...
             2*ones(size(R_mb)); 2*ones(size(R_ms))];  % 1=Sun 2=Mon

fprintf('  Pooled world-frame cloud: %d pts\n', numel(pool_e));
fprintf('  Tree bounding box: E=[%.1f %.1f]  N=[%.1f %.1f]  Z=[%.1f %.1f]\n', ...
    ATT_TREE_E_MIN, ATT_TREE_E_MAX, ATT_TREE_N_MIN, ATT_TREE_N_MAX, ...
    ATT_TREE_Z_MIN, ATT_TREE_Z_MAX);

% Define the four cones
% Each cone: origin (ENU), target GPS point (ENU), label, colour
cone_origins = {[0,0,0],       [0,0,0],      wf2_lidar,    wf2_lidar};
cone_targets = {wf1_small,     wf1_big,      wf2_small,    wf2_big};
cone_labels  = {'Sun->small',  'Sun->big',   'Mon->small', 'Mon->big'};
cone_colors  = {[1.0 0.9 0.1], [1.0 0.5 0.1], [0.2 0.5 1.0], [0.2 0.85 0.4]};

n_cones    = 4;
bin_edges  = BIN_START:BIN_WIDTH:BIN_END;
bin_centres= bin_edges(1:end-1) + BIN_WIDTH/2;
n_bins     = numel(bin_centres);
cone_rad   = deg2rad(CONE_HALF_ANGLE_DEG);

% Cone shell volumes for density normalisation
cone_vol = (2/3)*pi*(1-cos(cone_rad)) * ...
    ((bin_edges(2:end)).^3 - (bin_edges(1:end-1)).^3);

fprintf('\n  Cone half-angle: %.0f deg  bins: %.0f-%.0f m  width=%.1f m\n', ...
    CONE_HALF_ANGLE_DEG, BIN_START, BIN_END, BIN_WIDTH);

% Pre-allocate result arrays
bin_count_all = zeros(n_cones, n_bins);
bin_refl_all  = nan(n_cones, n_bins);
bin_dens_all  = zeros(n_cones, n_bins);

for c = 1:n_cones
    orig = cone_origins{c};
    tgt  = cone_targets{c};

    % Vector from origin to target — cone axis direction
    axis_vec = tgt - orig;
    axis_len = norm(axis_vec);
    axis_dir = axis_vec / axis_len;   % unit vector

    % Translate pool to cone origin
    dE = pool_e - orig(1);
    dN = pool_n - orig(2);
    dZ = pool_z - orig(3);

    % Slant range from this cone's origin
    slant = sqrt(dE.^2 + dN.^2 + dZ.^2);

    % Angular offset from cone axis using dot product
    % cos(ang) = (d . axis_dir) / |d|
    dot_prod = dE*axis_dir(1) + dN*axis_dir(2) + dZ*axis_dir(3);
    cos_ang  = dot_prod ./ max(slant, 1e-6);
    cos_ang  = min(max(cos_ang, -1), 1);   % clamp for acos
    ang_off  = acos(cos_ang);

    % Tree bounding box in world frame
    in_box = pool_e>=ATT_TREE_E_MIN & pool_e<=ATT_TREE_E_MAX & ...
             pool_n>=ATT_TREE_N_MIN & pool_n<=ATT_TREE_N_MAX & ...
             pool_z>=ATT_TREE_Z_MIN & pool_z<=ATT_TREE_Z_MAX;

    in_cone   = ang_off <= cone_rad;
    in_sample = in_cone & in_box;

    fprintf('\n  %s: %d pts in cone+box  (origin=[%.1f %.1f %.1f]  target=[%.1f %.1f %.1f])\n', ...
        cone_labels{c}, sum(in_sample), orig, tgt);

    for b = 1:n_bins
        in_bin = in_sample & slant>=bin_edges(b) & slant<bin_edges(b+1);
        n = sum(in_bin);
        bin_count_all(c,b) = n;
        if n > 0
            bin_refl_all(c,b) = mean(pool_refl(in_bin));
        end
    end

    bin_dens_all(c,:) = bin_count_all(c,:) ./ cone_vol;

    % Print table
    fprintf('    Range(m)  Count  Density  MeanRefl\n');
    for b = 1:n_bins
        if bin_count_all(c,b) > 0
            fprintf('    %5.1f  %6d  %7.2f  %6.1f\n', ...
                bin_centres(b), bin_count_all(c,b), ...
                bin_dens_all(c,b), bin_refl_all(c,b));
        end
    end
end

% ── Per-bin ratios and cumulative dB ─────────────────────────────────────
% Sun: small/big ratio (cone 1 / cone 2)
% Mon: small/big ratio (cone 3 / cone 4)
ratio_sun_wf = nan(1,n_bins);
ratio_mon_wf = nan(1,n_bins);

fprintf('\n  Sun small/big ratio per bin:\n');
fprintf('  Range  Refl_s  Refl_b  Ratio  Dens_s  Dens_b\n');
for b = 1:n_bins
    r1=bin_refl_all(1,b); r2=bin_refl_all(2,b);
    if ~isnan(r1)&&~isnan(r2)&&r2>0
        ratio_sun_wf(b)=r1/r2;
        fprintf('  %5.1f  %6.1f  %6.1f  %5.3f  %6.2f  %6.2f\n', ...
            bin_centres(b),r1,r2,ratio_sun_wf(b), ...
            bin_dens_all(1,b),bin_dens_all(2,b));
    end
end

fprintf('\n  Mon small/big ratio per bin:\n');
fprintf('  Range  Refl_s  Refl_b  Ratio  Dens_s  Dens_b\n');
for b = 1:n_bins
    r3=bin_refl_all(3,b); r4=bin_refl_all(4,b);
    if ~isnan(r3)&&~isnan(r4)&&r4>0
        ratio_mon_wf(b)=r3/r4;
        fprintf('  %5.1f  %6.1f  %6.1f  %5.3f  %6.2f  %6.2f\n', ...
            bin_centres(b),r3,r4,ratio_mon_wf(b), ...
            bin_dens_all(3,b),bin_dens_all(4,b));
    end
end

% Cumulative dB
fprintf('\n  Cumulative attenuation dB (Sun small vs big):\n');
cum_sun_wf = nan(1,n_bins); cum=0;
for b=1:n_bins
    if ~isnan(ratio_sun_wf(b))&&ratio_sun_wf(b)>0
        cum=cum+10*log10(ratio_sun_wf(b));
        cum_sun_wf(b)=cum;
        fprintf('  %5.1f  ratio=%5.3f  per-bin=%+6.2fdB  cum=%+7.2fdB\n', ...
            bin_centres(b),ratio_sun_wf(b),10*log10(ratio_sun_wf(b)),cum);
    end
end

fprintf('\n  Cumulative attenuation dB (Mon small vs big):\n');
cum_mon_wf = nan(1,n_bins); cum=0;
for b=1:n_bins
    if ~isnan(ratio_mon_wf(b))&&ratio_mon_wf(b)>0
        cum=cum+10*log10(ratio_mon_wf(b));
        cum_mon_wf(b)=cum;
        fprintf('  %5.1f  ratio=%5.3f  per-bin=%+6.2fdB  cum=%+7.2fdB\n', ...
            bin_centres(b),ratio_mon_wf(b),10*log10(ratio_mon_wf(b)),cum);
    end
end

% =========================================================================
fprintf('\n=== SECTION 9: ATTENUATION PLOTS (world frame) ===\n');
% =========================================================================
att_colors = cone_colors;

%% Figure: Point density vs distance
figure('Name','Point density vs distance — world frame','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.02 0.05 0.47 0.85]);
ax_d=axes('Color','k','XColor','w','YColor','w'); hold(ax_d,'on');
ls_all={'-','--','-','--'};
for c=1:n_cones
    plot(ax_d, bin_centres, bin_dens_all(c,:), ...
        [ls_all{c} 'o'],'Color',att_colors{c},'LineWidth',2,'MarkerSize',5, ...
        'DisplayName',cone_labels{c});
end
xlabel(ax_d,'Slant range from LiDAR origin (m)','Color','w');
ylabel(ax_d,'Point density (pts/m^3)','Color','w');
title(ax_d,sprintf('Point density — world frame (cone half-angle=%.0f deg)', ...
    CONE_HALF_ANGLE_DEG),'Color','w','FontSize',12);
grid(ax_d,'on'); ax_d.GridColor=[0.3 0.3 0.3]; ax_d.XLim=[BIN_START BIN_END];
legend(ax_d,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

%% Figure: Reflectivity, ratio, cumulative dB
figure('Name','Attenuation analysis — world frame','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.51 0.05 0.47 0.85]);

ax_r=subplot(3,1,1); set(ax_r,'Color','k','XColor','w','YColor','w'); hold(ax_r,'on');
for c=1:n_cones
    plot(ax_r, bin_centres, bin_refl_all(c,:), ...
        [ls_all{c} 'o'],'Color',att_colors{c},'LineWidth',2,'MarkerSize',5, ...
        'DisplayName',cone_labels{c});
end
ylabel(ax_r,'Mean reflectivity','Color','w');
title(ax_r,'Mean reflectivity vs distance — world frame','Color','w','FontSize',11);
grid(ax_r,'on'); ax_r.GridColor=[0.3 0.3 0.3];
ax_r.XLim=[BIN_START BIN_END]; ax_r.YLim=[0 100];
legend(ax_r,'show','TextColor','w','Color','k','FontSize',8,'Location','northwest');

ax_ratio=subplot(3,1,2); set(ax_ratio,'Color','k','XColor','w','YColor','w'); hold(ax_ratio,'on');
vs=~isnan(ratio_sun_wf); vm=~isnan(ratio_mon_wf);
if any(vs)
    plot(ax_ratio,bin_centres(vs),ratio_sun_wf(vs), ...
        'o-','Color',[1.0 0.7 0.1],'LineWidth',2,'MarkerSize',6,'DisplayName','Sun small/big');
end
if any(vm)
    plot(ax_ratio,bin_centres(vm),ratio_mon_wf(vm), ...
        's--','Color',[0.3 0.7 1.0],'LineWidth',2,'MarkerSize',6,'DisplayName','Mon small/big');
end
yline(ax_ratio,1,'w--','LineWidth',1.2,'DisplayName','Ratio=1 (equal)');
ylabel(ax_ratio,'Reflectivity ratio (small/big)','Color','w');
title(ax_ratio,'Small/big reflectivity ratio — Sun vs Mon','Color','w','FontSize',11);
grid(ax_ratio,'on'); ax_ratio.GridColor=[0.3 0.3 0.3]; ax_ratio.XLim=[BIN_START BIN_END];
legend(ax_ratio,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

ax_db=subplot(3,1,3); set(ax_db,'Color','k','XColor','w','YColor','w'); hold(ax_db,'on');
vs_db=~isnan(cum_sun_wf); vm_db=~isnan(cum_mon_wf);
if any(vs_db)
    plot(ax_db,bin_centres(vs_db),cum_sun_wf(vs_db), ...
        'o-','Color',[1.0 0.7 0.1],'LineWidth',2,'MarkerSize',6,'DisplayName','Sun cumulative dB');
end
if any(vm_db)
    plot(ax_db,bin_centres(vm_db),cum_mon_wf(vm_db), ...
        's--','Color',[0.3 0.7 1.0],'LineWidth',2,'MarkerSize',6,'DisplayName','Mon cumulative dB');
end
yline(ax_db,0,'w--','LineWidth',1,'DisplayName','0 dB (no loss)');
xlabel(ax_db,'Slant range from LiDAR origin (m)','Color','w');
ylabel(ax_db,'Cumulative attenuation (dB)','Color','w');
title(ax_db,'Cumulative attenuation — Sun vs Mon','Color','w','FontSize',11);
grid(ax_db,'on'); ax_db.GridColor=[0.3 0.3 0.3]; ax_db.XLim=[BIN_START BIN_END];
legend(ax_db,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

%% Figure: 3D world-frame cone visualisation
figure('Name','Cone sections — world frame 3D','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.02 0.05 0.65 0.85]);
ax_3d=axes('Color','k','XColor','w','YColor','w','ZColor','w'); hold(ax_3d,'on');

% Background cloud (faint)
scatter3(ax_3d,pool_e,pool_n,pool_z,1,[0.25 0.25 0.25],'filled', ...
    'MarkerFaceAlpha',0.08,'DisplayName','All pts');

% GPS reflector markers
plot3(ax_3d,wf1_big(1),   wf1_big(2),   wf1_big(3),   'r^','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Sun big refl');
plot3(ax_3d,wf1_small(1), wf1_small(2), wf1_small(3), 'g^','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Sun small refl');
plot3(ax_3d,wf2_big(1),   wf2_big(2),   wf2_big(3),   'rv','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Mon big refl');
plot3(ax_3d,wf2_small(1), wf2_small(2), wf2_small(3), 'gv','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Mon small refl');
plot3(ax_3d,wf2_lidar(1), wf2_lidar(2), wf2_lidar(3), 'bs','MarkerSize',12,'MarkerFaceColor','b','DisplayName','Mon LiDAR');
plot3(ax_3d,0,0,0,'wx','MarkerSize',14,'LineWidth',2,'DisplayName','Sun LiDAR');

% Draw cone axes and plot cone points
for c = 1:n_cones
    orig = cone_origins{c};
    tgt  = cone_targets{c};
    axis_vec = tgt - orig;
    axis_dir = axis_vec / norm(axis_vec);

    % Find cone points in pool
    dE=pool_e-orig(1); dN=pool_n-orig(2); dZ_c=pool_z-orig(3);
    slant_c=sqrt(dE.^2+dN.^2+dZ_c.^2);
    dot_p=dE*axis_dir(1)+dN*axis_dir(2)+dZ_c*axis_dir(3);
    cos_a=min(max(dot_p./max(slant_c,1e-6),-1),1);
    ang_c=acos(cos_a);
    in_box_c = pool_e>=ATT_TREE_E_MIN&pool_e<=ATT_TREE_E_MAX & ...
               pool_n>=ATT_TREE_N_MIN&pool_n<=ATT_TREE_N_MAX & ...
               pool_z>=ATT_TREE_Z_MIN&pool_z<=ATT_TREE_Z_MAX;
    in_c = ang_c<=cone_rad & in_box_c;

    if any(in_c)
        ss=1:4:sum(in_c);
        ep=pool_e(in_c); np=pool_n(in_c); zp=pool_z(in_c);
        scatter3(ax_3d,ep(ss),np(ss),zp(ss),5, ...
            repmat(reshape(att_colors{c},1,3),numel(ss),1),'filled', ...
            'MarkerFaceAlpha',0.6,'DisplayName',[cone_labels{c} sprintf(' (%d pts)',sum(in_c))]);
    end

    % Cone axis line (origin to 1.2× BIN_END along axis)
    ep2=orig(1)+axis_dir(1)*BIN_END*1.1;
    np2=orig(2)+axis_dir(2)*BIN_END*1.1;
    zp2=orig(3)+axis_dir(3)*BIN_END*1.1;
    plot3(ax_3d,[orig(1) ep2],[orig(2) np2],[orig(3) zp2],'--', ...
        'Color',att_colors{c},'LineWidth',2,'DisplayName',['Axis: ' cone_labels{c}]);
end

xlabel(ax_3d,'East (m)','Color','w'); ylabel(ax_3d,'North (m)','Color','w');
zlabel(ax_3d,'Z (m)','Color','w');
title(ax_3d,sprintf('World-frame cone sections (half-angle=%.0f deg)',CONE_HALF_ANGLE_DEG), ...
    'Color','w','FontSize',12);
axis(ax_3d,'equal'); grid(ax_3d,'on'); ax_3d.GridColor=[0.3 0.3 0.3];
view(ax_3d,45,25); rotate3d(ax_3d,'on');
legend(ax_3d,'show','TextColor','w','Color','k','FontSize',7,'Location','northeast');

fprintf('  3 attenuation figures generated.\n\nDone.\n');