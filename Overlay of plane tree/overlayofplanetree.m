% =========================================================================
% align_radar_scans.m
% Aligns two radar scans into a common GPS world frame using corner
% reflectors as tie points. Verification prints after every section.
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
file_scan1 = 'pointcloud_20260315-143039_.csv';
file_scan2 = 'pointcloud_20000101-140013_.csv';
file_gps1  = 'gps_tree scan 1.csv';
file_gps2  = 'tree scan 2.csv';

AMP_THRESH_SCAN1 = 50;   % dB — both reflectors strong in scan 1
AMP_THRESH_SCAN2 = 20;   % dB — weak reflector in scan 2 only reaches ~35 dB
RANGE_LIMIT      = 40;   % m  — discard far-field returns

% Rotation correction applied to Scan 2 after transform (degrees about Z)
% Change to -90 or 180 if overlay looks wrong
ROT_CORRECTION_DEG = 0;



% =========================================================================
fprintf('=== SECTION 1: LOAD RADAR DATA ===\n');
% =========================================================================
T1 = readtable(file_scan1);
T2 = readtable(file_scan2);

fprintf('  Scan 1 raw rows: %d  |  columns: %s\n', height(T1), strjoin(T1.Properties.VariableNames, ', '));
fprintf('  Scan 2 raw rows: %d\n', height(T2));
fprintf('  Scan 1 range: %.2f to %.2f m\n', min(T1.range), max(T1.range));
fprintf('  Scan 2 range: %.2f to %.2f m\n', min(T2.range), max(T2.range));
fprintf('  Scan 1 amplitude: %.1f to %.1f dB\n', min(T1.amplitude), max(T1.amplitude));
fprintf('  Scan 2 amplitude: %.1f to %.1f dB\n', min(T2.amplitude), max(T2.amplitude));

% =========================================================================
fprintf('\n=== SECTION 2: SPHERICAL -> CARTESIAN (sensor frame) ===\n');
% =========================================================================
% Scanner geometry: pitch is angle from nadir (straight down), negative values
% mean the scanner looks forward/down. Roll is azimuth scan angle.
%
%   horizontal_dist = r * sin(-pitch)
%   vertical (down)  = r * cos(-pitch)
%
%   X = r * sin(-pitch) * cos(roll)
%   Y = r * sin(-pitch) * sin(roll)
%   Z = -r * cos(-pitch)    (negative = below radar)
%
% Verified: small reflector at range=10.02m, pitch=-88deg gives
%   horizontal = 10.02 * sin(88deg) = 10.01m  (GPS = 9.96m) ✓

X1_all = T1.range .* sin(-T1.pitch) .* cos(T1.roll);
Y1_all = T1.range .* sin(-T1.pitch) .* sin(T1.roll);
Z1_all = -T1.range .* cos(-T1.pitch);
amp1_all = T1.amplitude;

X2_all = T2.range .* sin(-T2.pitch) .* cos(T2.roll);
Y2_all = T2.range .* sin(-T2.pitch) .* sin(T2.roll);
Z2_all = -T2.range .* cos(-T2.pitch);
amp2_all = T2.amplitude;

% Apply range filter
m1 = T1.range <= RANGE_LIMIT;
m2 = T2.range <= RANGE_LIMIT;

X1 = X1_all(m1); Y1 = Y1_all(m1); Z1 = Z1_all(m1); amp1 = amp1_all(m1);
X2 = X2_all(m2); Y2 = Y2_all(m2); Z2 = Z2_all(m2); amp2 = amp2_all(m2);

fprintf('  Scan 1: %d points after %.0fm range filter\n', sum(m1), RANGE_LIMIT);
fprintf('  Scan 2: %d points after %.0fm range filter\n', sum(m2), RANGE_LIMIT);
fprintf('  Scan 1 Cartesian X: %.2f to %.2f m\n', min(X1), max(X1));
fprintf('  Scan 1 Cartesian Y: %.2f to %.2f m\n', min(Y1), max(Y1));
fprintf('  Scan 1 Cartesian Z: %.2f to %.2f m\n', min(Z1), max(Z1));

% =========================================================================
fprintf('\n=== SECTION 3: EXTRACT REFLECTOR CENTROIDS ===\n');
% =========================================================================
% Reflectors appear as tight clusters at distinct ranges with high amplitude.
% Strategy: threshold amplitude, then split by largest range gap.

for scan_id = 1:2
    if scan_id == 1
        Xs = X1; Ys = Y1; Zs = Z1; amps = amp1;
        Rs = T1.range(m1); thresh = AMP_THRESH_SCAN1;
    else
        Xs = X2; Ys = Y2; Zs = Z2; amps = amp2;
        Rs = T2.range(m2); thresh = AMP_THRESH_SCAN2;
    end

    % Amplitude filter + remove far-field clutter (Z < -20m)
    mask  = (amps > thresh) & (Zs > -20);
    n_pts = sum(mask);
    fprintf('  Scan %d: %d points above %d dB threshold\n', scan_id, n_pts, thresh);

    if n_pts < 2
        error('Not enough high-amplitude points in Scan %d. Lower AMP_THRESH.', scan_id);
    end

    pts_r = Rs(mask);
    pts_x = Xs(mask); pts_y = Ys(mask); pts_z = Zs(mask);

    % Find largest range gap to split into two reflector clusters
    [pts_r_sorted, sort_idx] = sort(pts_r);
    gaps = diff(pts_r_sorted);
    [max_gap, gap_idx] = max(gaps);
    fprintf('  Scan %d: largest range gap = %.2f m (between %.2f m and %.2f m)\n', ...
        scan_id, max_gap, pts_r_sorted(gap_idx), pts_r_sorted(gap_idx+1));

    near_idx = sort_idx(1:gap_idx);
    far_idx  = sort_idx(gap_idx+1:end);

    c_near = [mean(pts_x(near_idx)), mean(pts_y(near_idx)), mean(pts_z(near_idx))];
    c_far  = [mean(pts_x(far_idx)),  mean(pts_y(far_idx)),  mean(pts_z(far_idx))];

    fprintf('  Scan %d near centroid: x=%6.3f  y=%6.3f  z=%6.3f  (range~%.1fm, %d pts, max amp=%.0f dB)\n', ...
        scan_id, c_near(1), c_near(2), c_near(3), mean(pts_r(near_idx)), length(near_idx), max(amps(mask & Rs<=pts_r_sorted(gap_idx))));
    fprintf('  Scan %d far  centroid: x=%6.3f  y=%6.3f  z=%6.3f  (range~%.1fm, %d pts, max amp=%.0f dB)\n', ...
        scan_id, c_far(1),  c_far(2),  c_far(3),  mean(pts_r(far_idx)),  length(far_idx),  max(amps(mask & Rs>pts_r_sorted(gap_idx))));
    fprintf('  Scan %d sensor-frame reflector separation: %.3f m\n', scan_id, norm(c_near - c_far));

    if scan_id == 1; c1_near = c_near; c1_far = c_far;
    else;            c2_near = c_near; c2_far = c_far; end
end

% =========================================================================
fprintf('\n=== SECTION 4: LOAD GPS AND CONVERT TO ENU ===\n');
% =========================================================================
opts1 = detectImportOptions(file_gps1, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
G1    = readtable(file_gps1, opts1);
opts2 = detectImportOptions(file_gps2, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
G2    = readtable(file_gps2, opts2);

G1.Properties.VariableNames{1} = 'PointName';
G2.Properties.VariableNames{1} = 'PointName';
G1.PointName = strtrim(G1.PointName);
G2.PointName = strtrim(G2.PointName);

fprintf('  G1 loaded: %d rows, %d cols\n', height(G1), width(G1));
fprintf('  G1 point names: %s\n', strjoin(G1.PointName, ' | '));
fprintf('  G2 loaded: %d rows, %d cols\n', height(G2), width(G2));
fprintf('  G2 point names: %s\n', strjoin(G2.PointName, ' | '));

% World frame origin = Scan 1 radar GPS position, Z=0 baseline
idx  = strcmp(G1.PointName, 'radar/lidar');
olon = G1.Longitude(idx);
olat = G1.Latitude(idx);
oh   = G1.("Ellipsoidal height")(idx);
fprintf('  World origin: lon=%.8f  lat=%.8f  h=%.3f m\n', olon, olat, oh);

gps2enu = @(lon, lat, h) deal( ...
    (lon - olon) .* 111320 .* cos(deg2rad(olat)), ...
    (lat - olat) .* 110540, ...
     h - oh);

% Scan 1 reflectors
[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector big')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector big')));
wf1_big = [e,n,u];

[e,n,u] = gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.Latitude(strcmp(G1.PointName,'corner reflector small')), ...
                  G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector small')));
wf1_small = [e,n,u];

% Scan 2 reflectors — use 'big1' (second, better reading) not average
idx_big2 = strcmp(G2.PointName, 'corner refelctor big1');
[e,n,u]  = gps2enu(G2.Longitude(idx_big2), G2.Latitude(idx_big2), ...
                   G2.("Ellipsoidal height")(idx_big2));
wf2_big = [e,n,u];

idx_s2   = strcmp(G2.PointName, 'corner refelctor small');  % typo in source data
[e,n,u]  = gps2enu(G2.Longitude(idx_s2), G2.Latitude(idx_s2), G2.("Ellipsoidal height")(idx_s2));
wf2_small = [e,n,u];

idx_r2   = strcmp(G2.PointName, 'Radar/Lidar');
[e,n,u]  = gps2enu(G2.Longitude(idx_r2), G2.Latitude(idx_r2), G2.("Ellipsoidal height")(idx_r2));
wf2_radar = [e,n,u];

fprintf('  Scan 1 big   reflector ENU: E=%7.3f  N=%7.3f  U=%7.3f  (horiz dist=%.2fm)\n', wf1_big,   norm(wf1_big(1:2)));
fprintf('  Scan 1 small reflector ENU: E=%7.3f  N=%7.3f  U=%7.3f  (horiz dist=%.2fm)\n', wf1_small, norm(wf1_small(1:2)));
fprintf('  Scan 2 big   reflector ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_big);
fprintf('  Scan 2 small reflector ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_small);
fprintf('  Scan 2 radar origin    ENU: E=%7.3f  N=%7.3f  U=%7.3f\n', wf2_radar);
fprintf('  Scan 1 GPS reflector separation: %.3f m\n', norm(wf1_big - wf1_small));
fprintf('  Scan 2 GPS reflector separation: %.3f m\n', norm(wf2_big  - wf2_small));

% =========================================================================
fprintf('\n=== SECTION 5: MATCH SENSOR CENTROIDS TO GPS POSITIONS ===\n');
% =========================================================================
% Match by comparing sensor range to GPS horizontal distance from each radar.
% The sensor near centroid should match the GPS reflector closer to that radar.

% Scan 1 — distances from Scan 1 radar (world origin = 0,0)
d1_big   = norm(wf1_big(1:2));
d1_small = norm(wf1_small(1:2));
fprintf('  Scan 1 GPS horiz distances from radar: big=%.2fm  small=%.2fm\n', d1_big, d1_small);
fprintf('  Scan 1 sensor ranges:                  near=%.2fm  far=%.2fm\n', norm(c1_near), norm(c1_far));

% Match sensor range to closest GPS distance
err_near_big   = abs(norm(c1_near) - d1_big);
err_near_small = abs(norm(c1_near) - d1_small);
if err_near_big < err_near_small
    s1_near_wf = wf1_big;   s1_near_name = 'big';
    s1_far_wf  = wf1_small; s1_far_name  = 'small';
else
    s1_near_wf = wf1_small; s1_near_name = 'small';
    s1_far_wf  = wf1_big;   s1_far_name  = 'big';
end
fprintf('  Scan 1 match: near->%s (err=%.2fm)  far->%s (err=%.2fm)\n', ...
    s1_near_name, min(err_near_big,err_near_small), s1_far_name, max(err_near_big,err_near_small));

% Scan 2 — distances from Scan 2 radar GPS position
d2_big   = norm(wf2_big(1:2)   - wf2_radar(1:2));
d2_small = norm(wf2_small(1:2) - wf2_radar(1:2));
fprintf('  Scan 2 GPS horiz distances from radar: big=%.2fm  small=%.2fm\n', d2_big, d2_small);
fprintf('  Scan 2 sensor ranges:                  near=%.2fm  far=%.2fm\n', norm(c2_near), norm(c2_far));

err_near_big2   = abs(norm(c2_near) - d2_big);
err_near_small2 = abs(norm(c2_near) - d2_small);
if err_near_big2 < err_near_small2
    s2_near_wf = wf2_big;   s2_near_name = 'big';
    s2_far_wf  = wf2_small; s2_far_name  = 'small';
else
    s2_near_wf = wf2_small; s2_near_name = 'small';
    s2_far_wf  = wf2_big;   s2_far_name  = 'big';
end
fprintf('  Scan 2 match: near->%s (err=%.2fm)  far->%s (err=%.2fm)\n', ...
    s2_near_name, min(err_near_big2,err_near_small2), s2_far_name, max(err_near_big2,err_near_small2));

% Third tie point: radar sensor origin (0,0,0) maps to its GPS world position.
% For Scan 1 the origin IS the world origin (0,0,0).
% For Scan 2 the origin maps to wf2_radar.
sensor_pts_1 = [c1_near; c1_far; 0, 0, 0];
world_pts_1  = [s1_near_wf; s1_far_wf; 0, 0, 0];

sensor_pts_2 = [c2_near; c2_far; 0, 0, 0];
world_pts_2  = [s2_near_wf; s2_far_wf; wf2_radar];

fprintf('  Tie points set (3 per scan: near reflector, far reflector, radar origin)\n');

% =========================================================================
fprintf('\n=== SECTION 6: SOLVE RIGID BODY TRANSFORMS (SVD) ===\n');
% =========================================================================
function [R, t] = solve_rigid_transform(sp, wp)
    mu_s = mean(sp, 1);
    mu_w = mean(wp, 1);
    S    = (sp - mu_s)' * (wp - mu_w);
    [U, ~, V] = svd(S);
    D = diag([1, 1, det(V * U')]);
    R = V * D * U';
    t = mu_w - mu_s * R';
end

[R1, t1] = solve_rigid_transform(sensor_pts_1, world_pts_1);
[R2, t2] = solve_rigid_transform(sensor_pts_2, world_pts_2);

fprintf('  R1 (rotation matrix Scan 1):\n');
fprintf('    [%6.3f %6.3f %6.3f]\n', R1(1,:));
fprintf('    [%6.3f %6.3f %6.3f]\n', R1(2,:));
fprintf('    [%6.3f %6.3f %6.3f]\n', R1(3,:));
fprintf('  t1 (translation Scan 1): [%.3f  %.3f  %.3f] m\n', t1);
fprintf('  R2 (rotation matrix Scan 2):\n');
fprintf('    [%6.3f %6.3f %6.3f]\n', R2(1,:));
fprintf('    [%6.3f %6.3f %6.3f]\n', R2(2,:));
fprintf('    [%6.3f %6.3f %6.3f]\n', R2(3,:));
fprintf('  t2 (translation Scan 2): [%.3f  %.3f  %.3f] m\n', t2);

fprintf('\n  Tie point residuals (should be < 1m ideally):\n');
for i = 1:3
    pw1 = sensor_pts_1(i,:) * R1' + t1;
    pw2 = sensor_pts_2(i,:) * R2' + t2;
    labels = {'near reflector', 'far reflector', 'vertical up'};
    fprintf('    Scan 1 %s: %.4f m\n', labels{i}, norm(pw1 - world_pts_1(i,:)));
    fprintf('    Scan 2 %s: %.4f m\n', labels{i}, norm(pw2 - world_pts_2(i,:)));
end

% =========================================================================
fprintf('\n=== SECTION 7: APPLY TRANSFORMS ===\n');
% =========================================================================
pts1_world = [X1, Y1, Z1] * R1' + t1;
pts2_world = [X2, Y2, Z2] * R2' + t2;

fprintf('  Scan 1 world-frame bounding box:\n');
fprintf('    E: %.2f to %.2f m\n', min(pts1_world(:,1)), max(pts1_world(:,1)));
fprintf('    N: %.2f to %.2f m\n', min(pts1_world(:,2)), max(pts1_world(:,2)));
fprintf('    U: %.2f to %.2f m\n', min(pts1_world(:,3)), max(pts1_world(:,3)));
fprintf('  Scan 2 world-frame bounding box (before rotation correction):\n');
fprintf('    E: %.2f to %.2f m\n', min(pts2_world(:,1)), max(pts2_world(:,1)));
fprintf('    N: %.2f to %.2f m\n', min(pts2_world(:,2)), max(pts2_world(:,2)));
fprintf('    U: %.2f to %.2f m\n', min(pts2_world(:,3)), max(pts2_world(:,3)));

% Apply rotation correction about Z axis centred on Scan 2 radar origin
theta  = deg2rad(ROT_CORRECTION_DEG);
Rz     = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
centre = wf2_radar;
pts2_world = (pts2_world - centre) * Rz' + centre;

fprintf('  Rotation correction: %+.0f degrees about Z (change ROT_CORRECTION_DEG if wrong)\n', ROT_CORRECTION_DEG);
fprintf('  Scan 2 world-frame bounding box (after rotation correction):\n');
fprintf('    E: %.2f to %.2f m\n', min(pts2_world(:,1)), max(pts2_world(:,1)));
fprintf('    N: %.2f to %.2f m\n', min(pts2_world(:,2)), max(pts2_world(:,2)));
fprintf('    U: %.2f to %.2f m\n', min(pts2_world(:,3)), max(pts2_world(:,3)));

% Sanity check: do the two clouds overlap?
overlap_E = [max(min(pts1_world(:,1)), min(pts2_world(:,1))), min(max(pts1_world(:,1)), max(pts2_world(:,1)))];
overlap_N = [max(min(pts1_world(:,2)), min(pts2_world(:,2))), min(max(pts1_world(:,2)), max(pts2_world(:,2)))];
if overlap_E(1) < overlap_E(2) && overlap_N(1) < overlap_N(2)
    fprintf('  SANITY CHECK PASSED: clouds overlap in E (%.1f to %.1f) and N (%.1f to %.1f)\n', overlap_E, overlap_N);
else
    fprintf('  WARNING: clouds do not overlap — check rotation correction or tie point matching\n');
end

% Z is ignored for alignment — GPS height uncertainty is too large and
% radar height above ground is unknown. Solve rotation/translation in
% XY only, then keep each scan's own Z relative to its radar height.

% 2D rigid body transform (rotation about Z + XY translation)
% From the SVD R matrices, extract just the yaw (Z rotation) component
function [theta, tx, ty] = solve_2d_transform(sp, wp)
    % sp, wp are Nx2 (XY only)
    mu_s = mean(sp, 1);
    mu_w = mean(wp, 1);
    S    = (sp - mu_s)' * (wp - mu_w);
    [U, ~, V] = svd(S);
    R2d  = V * diag([1, det(V*U')]) * U';
    theta = atan2(R2d(2,1), R2d(1,1));
    t2d  = mu_w - mu_s * R2d';
    tx   = t2d(1);
    ty   = t2d(2);
end

[theta1, tx1, ty1] = solve_2d_transform(...
    [sensor_pts_1(:,1), sensor_pts_1(:,2)], ...
    [world_pts_1(:,1),  world_pts_1(:,2)]);

[theta2, tx2, ty2] = solve_2d_transform(...
    [sensor_pts_2(:,1), sensor_pts_2(:,2)], ...
    [world_pts_2(:,1),  world_pts_2(:,2)]);

fprintf('  Scan 1: yaw=%.2f deg  tx=%.3f m  ty=%.3f m\n', rad2deg(theta1), tx1, ty1);
fprintf('  Scan 2: yaw=%.2f deg  tx=%.3f m  ty=%.3f m\n', rad2deg(theta2), tx2, ty2);

% Apply 2D transforms — XY rotated and translated
% Z: keep each scan relative to its own radar, then apply GPS height offset
% Option 1: shift Scan 2 Z down by wf2_radar(3) so both scans share
%           the same Z=0 baseline (Scan 1 radar height)
c1 = cos(theta1); s1 = sin(theta1);
c2 = cos(theta2); s2 = sin(theta2);

pts1_xy = [X1, Y1] * [c1 s1; -s1 c1] + [tx1, ty1];
pts2_xy = [X2, Y2] * [c2 s2; -s2 c2] + [tx2, ty2];

% Flip Z (scanner looks down so Z is negative upward — negate to make up = positive)
% Then offset Scan 2 by the GPS radar height difference so both share Z=0 baseline
z_offset = wf2_radar(3);   % = 6.137 m — Scan 2 radar was this much higher than Scan 1
z_offset_empirical = 4.3;   % metres — adjust if ground planes don't align
pts1_world = [pts1_xy, -Z1];
pts2_world = [pts2_xy, -Z2 + z_offset - z_offset_empirical];

fprintf('  Z Option 1: Scan 2 shifted by GPS radar height offset = %.3f m\n', z_offset);
fprintf('  Z: GPS offset=%.3fm  empirical correction=%.3fm  net shift=%.3fm\n', ...
    z_offset, z_offset_empirical, z_offset - z_offset_empirical);

% ── Statistical outlier removal ──────────────────────────────────────────
% For each point, compute the mean distance to its K nearest neighbours.
% Points whose mean neighbour distance exceeds MEAN_DIST + STD_MULTIPLIER
% standard deviations above the global mean are removed as outliers.

K              = 20;    % number of nearest neighbours to check
STD_MULTIPLIER = 1.5;   % how many std devs above mean = outlier (lower = more aggressive)

function pts_clean = remove_outliers(pts, k, std_mult)
    if isempty(pts); pts_clean = pts; return; end
    idx = knnsearch(pts, pts, 'K', k+1);   % +1 because point finds itself
    idx = idx(:, 2:end);                    % remove self
    dists = sqrt(sum((pts - pts(idx(:,1),:)).^2, 2));
    for i = 2:k
        dists = dists + sqrt(sum((pts - pts(idx(:,i),:)).^2, 2));
    end
    mean_dist = dists / k;
    threshold = mean(mean_dist) + std_mult * std(mean_dist);
    pts_clean = pts(mean_dist <= threshold, :);
end

n_before1 = size(pts1_world, 1);
n_before2 = size(pts2_world, 1);

pts1_xyz = pts1_world(:,1:3);
pts2_xyz = pts2_world(:,1:3);

pts1_clean_xyz = remove_outliers(pts1_xyz, K, STD_MULTIPLIER);
pts2_clean_xyz = remove_outliers(pts2_xyz, K, STD_MULTIPLIER);



% Re-match amplitude to cleaned points using nearest neighbour lookup
idx1 = knnsearch(pts1_xyz, pts1_clean_xyz, 'K', 1);
idx2 = knnsearch(pts2_xyz, pts2_clean_xyz, 'K', 1);
amp1 = amp1(idx1);
amp2 = amp2(idx2);
pts1_world = [pts1_clean_xyz, pts1_world(idx1, 4:end)];
pts2_world = [pts2_clean_xyz, pts2_world(idx2, 4:end)];

fprintf('  Outlier removal (K=%d, %.1f std): Scan1 %d->%d pts  Scan2 %d->%d pts\n', ...
    K, STD_MULTIPLIER, n_before1, size(pts1_world,1), n_before2, size(pts2_world,1));


% ── Exclusion zones — remove sidelobe noise regions ──────────────────────
% Two circular exclusion zones identified visually from the plot.
% Any point whose XY position falls within these zones is removed.

% Exclusion zone 1 — bottom-left sidelobe arc (Scan 1)
EX1_N = 12;   EX1_E = -6;   EX1_R = 10;   % centre North, centre East, radius (m)

% Exclusion zone 2 — top-right sidelobe arc (Scan 2)  
EX2_N = 44;   EX2_E = -24;  EX2_R = 8;

function [pts_out, amp_out] = exclude_zones(pts, amp, zones)
    keep = true(size(pts,1), 1);
    for z = 1:size(zones,1)
        n_c = zones(z,1); e_c = zones(z,2); r = zones(z,3);
        dist = sqrt((pts(:,2)-n_c).^2 + (pts(:,1)-e_c).^2);
        keep = keep & (dist > r);
    end
    pts_out = pts(keep,:);
    amp_out = amp(keep);
end

zones = [EX1_N, EX1_E, EX1_R;
         EX2_N, EX2_E, EX2_R];

[pts1_world, amp1] = exclude_zones(pts1_world, amp1, zones);
[pts2_world, amp2] = exclude_zones(pts2_world, amp2, zones);

fprintf('  After exclusion zones: Scan1=%d pts  Scan2=%d pts\n', ...
    size(pts1_world,1), size(pts2_world,1));
% =========================================================================
fprintf('\n=== SECTION 8: PLOT ===\n');
% =========================================================================
figure('Name','Radar Scan Overlay','Color','white','Position',[100 100 1200 800]);

scatter3(pts1_world(:,1), pts1_world(:,2), pts1_world(:,3), ...
    6, amp1, 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Scan 1');
hold on;
scatter3(pts2_world(:,1), pts2_world(:,2), pts2_world(:,3), ...
    6, amp2, 's', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Scan 2');

plot3(wf1_big(1),   wf1_big(2),   wf1_big(3),   'r^', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Scan 1 big reflector');
plot3(wf1_small(1), wf1_small(2), wf1_small(3), 'g^', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Scan 1 small reflector');
plot3(wf2_big(1),   wf2_big(2),   wf2_big(3),   'rv', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Scan 2 big reflector');
plot3(wf2_small(1), wf2_small(2), wf2_small(3), 'gv', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'g', 'DisplayName', 'Scan 2 small reflector');
plot3(wf2_radar(1), wf2_radar(2), wf2_radar(3), 'bs', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Scan 2 radar origin GPS');

colormap(jet); cb = colorbar; cb.Label.String = 'Amplitude (dB)';
clim([-100 100]);
xlabel('East (m)'); ylabel('North (m)'); zlabel('Up — rel. to Scan 1 radar Z=0 (m)');
title('Radar scan overlay — GPS world frame (ENU)');
legend('Location', 'northeast', 'FontSize', 9);
grid on; axis equal; view(45, 25);

fprintf('  Plot complete.\n');
fprintf('\n  If scans are still orthogonal, try ROT_CORRECTION_DEG = -90 or 180\n');
fprintf('  If scans overlap but are mirrored, swap s1_near_wf/s1_far_wf in Section 5\n');

% =========================================================================
%% CONE ATTENUATION — BOTH SCANS MERGED
% =========================================================================

HALF_ANGLE_TREE = 15;   % degrees — adjust if too few points
BIN_WIDTH_TREE  = 0.5;  % m
tree_centre_world = [-11.9, 32.0, 9.0];  % ENU world frame

% ---- SCAN 1 CONE ----
% Boresight from scan 1 radar origin (0,0,0) toward tree centre
bs1 = tree_centre_world / norm(tree_centre_world);

pts1 = pts1_world(:,1:3);
norms1 = sqrt(sum(pts1.^2, 2));
pts1_unit = pts1 ./ norms1;

angle1 = acosd(pts1_unit * bs1');
in_cone1 = angle1 <= HALF_ANGLE_TREE;
range1_cone = norms1(in_cone1);
amp1_cone   = amp1(in_cone1);

fprintf('Scan 1 points in cone: %d\n', sum(in_cone1));

% ---- SCAN 2 CONE ----
% Boresight from scan 2 radar origin toward tree centre
radar2_origin = wf2_radar(1:3);
bs2_vec = tree_centre_world - radar2_origin;
bs2 = bs2_vec / norm(bs2_vec);

pts2 = pts2_world(:,1:3);
% Range from scan 2 radar origin
pts2_rel = pts2 - radar2_origin;
norms2 = sqrt(sum(pts2_rel.^2, 2));
pts2_unit = pts2_rel ./ norms2;

angle2 = acosd(pts2_unit * bs2');
in_cone2 = angle2 <= HALF_ANGLE_TREE;
range2_cone = norms2(in_cone2);
amp2_cone   = amp2(in_cone2);

fprintf('Scan 2 points in cone: %d\n', sum(in_cone2));

% ---- BIN BY RANGE ----
r_edges   = 0 : BIN_WIDTH_TREE : 50;
r_centres = r_edges(1:end-1) + BIN_WIDTH_TREE/2;

counts1 = histcounts(range1_cone, r_edges);
counts2 = histcounts(range2_cone, r_edges);
counts_merged = counts1 + counts2;
counts_merged_norm = counts_merged / max(counts_merged);

radar_r_centres = r_centres;
radar_counts_norm = counts_merged_norm;
save('radar_cone_density.mat', 'radar_r_centres', 'radar_counts_norm')
fprintf('Saved: radar_cone_density.mat\n');

% Mean amplitude per bin — merged
amp_bin = zeros(1, numel(r_centres));
for b = 1:numel(r_centres)
    in_bin1 = range1_cone >= r_edges(b) & range1_cone < r_edges(b+1);
    in_bin2 = range2_cone >= r_edges(b) & range2_cone < r_edges(b+1);
    all_amp = [amp1_cone(in_bin1); amp2_cone(in_bin2)];
    if ~isempty(all_amp)
        amp_bin(b) = mean(all_amp);
    else
        amp_bin(b) = NaN;
    end
end

% ---- FIGURE 1: Normalised density ----
figure('Color','w','Units','normalized','Position',[0.1 0.3 0.6 0.45]);

plot(r_centres, counts_merged_norm, '-', ...
    'Color', [0.9 0.3 0.1], 'LineWidth', 2.0, ...
    'DisplayName', 'Radar 120 GHz (both scans)');

xline(14, '--k', 'LineWidth', 1.2, 'Label', 'Canopy front', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');
xline(40, '--k', 'LineWidth', 1.2, 'Label', 'Canopy back', ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');

xlabel('Range from sensor (m)', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Normalised point density', 'FontSize', 12, 'Interpreter', 'none');
title('Radar point density through plane tree canopy - both scans merged', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
xlim([5 50]);
ylim([0 1.05]);

exportgraphics(gcf, 'fig_tree_cone_density.png', 'Resolution', 300);
fprintf('Saved: fig_tree_cone_density.png\n');

% ---- FIGURE 2: Mean amplitude vs range ----
valid_bins = ~isnan(amp_bin) & counts_merged > 5;

figure('Color','w','Units','normalized','Position',[0.1 0.3 0.6 0.45]);

plot(r_centres(valid_bins), amp_bin(valid_bins), 'o-', ...
    'Color', [0.9 0.3 0.1], 'LineWidth', 2.0, 'MarkerSize', 5, ...
    'DisplayName', 'Mean amplitude per bin');

% Linear fit over canopy region — adjust bounds after seeing the plot
canopy_mask = valid_bins & r_centres >= 14 & r_centres <= 40;
if sum(canopy_mask) >= 3
    p = polyfit(r_centres(canopy_mask), amp_bin(canopy_mask), 1);
    fit_x = r_centres(canopy_mask);
    fit_y = polyval(p, fit_x);
    hold on;
    plot(fit_x, fit_y, '--k', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Linear fit: k=%.2f dB/m', p(1)));
    fprintf('Attenuation fit: k = %.2f dB/m  intercept = %.1f dB\n', p(1), p(2));
end
xline(14, '--k', 'LineWidth', 1.2, 'Label', 'Canopy front', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');
xline(40, '--k', 'LineWidth', 1.2, 'Label', 'Canopy back', ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');

xlabel('Range from sensor (m)', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Mean amplitude (dB)', 'FontSize', 12, 'Interpreter', 'none');
title('Mean radar amplitude vs range through plane tree canopy', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
xlim([5 50]);

exportgraphics(gcf, 'fig_tree_cone_amplitude.png', 'Resolution', 300);
fprintf('Saved: fig_tree_cone_amplitude.png\n');