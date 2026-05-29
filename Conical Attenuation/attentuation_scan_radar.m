% =========================================================================
% conical_attenuation_sunday.m
%
% Conical section attenuation analysis — Sunday plane tree scan
% Compares 120 GHz radar vs LiDAR 905 nm attenuation through canopy.
%
% No control scan required. Fits attenuation slope (dB/m) directly from
% the tree scan profiles for each sensor and each reflector cone.
%
% Files required (same directory as this script):
%   sundar_radar_scan.csv
%   sunday_lidar_corner_big.csv
%   sunday_lidar_corner_small.csv
%   gps_tree_scan_1.csv  (used for reference only — centroids pre-computed)
% =========================================================================

clear; clc; close all;

% =========================================================================
% SETTINGS
% =========================================================================

RADAR_FILE       = 'sundar_radar_scan.csv';
LIDAR_FILE_BIG   = 'sunday_lidar_corner_big.csv';
LIDAR_FILE_SMALL = 'sunday_lidar_corner_small.csv';

% Reflector centroids in RADAR sensor frame
% (computed from align_radar_scans.m Section 3 + GPS matching)
REF_BIG_XYZ_RADAR   = [11.615,  6.078,  0.100];   % big reflector,   range 13.11m
REF_SMALL_XYZ_RADAR = [ 5.657,  8.264, -0.256];   % small reflector, range 10.02m

% Reflector centroids in LIDAR sensor frame
% (detected from high-reflectivity returns in each LiDAR file)
REF_BIG_XYZ_LIDAR   = [12.892, -2.650, -0.480];   % big in LiDAR big file,   range 13.17m
REF_SMALL_XYZ_LIDAR = [ 9.640,  2.887, -0.157];   % small in LiDAR small file, range 10.07m

% Known GPS ranges (for reference lines on plots)
R_BIG   = 13.11;   % metres
R_SMALL = 10.02;   % metres

% Cone half-angle (degrees)
CONE_DEG = 10;

% Range bin width
BIN_W = 0.5;   % metres

% Minimum points per bin to include in profile
MIN_PTS = 5;

% Radar noise floor filter
RADAR_AMP_MIN   = -80;   % dB
RADAR_RANGE_MAX =  40;   % m

% LiDAR range limit
LIDAR_RANGE_MAX = 40;   % m

% ── Fit windows (chosen from inspecting profiles) ─────────────────────────
% Radar big:   attenuation visible from 10.75 to 12.75m (excludes reflector spike at 10.25m)
% Radar small: profile too noisy — fit 6.75 to 9.75m (before reflector spike)
% LiDAR big:   clear drop from 7.25 to 13.25m
% LiDAR small: drop from 8.75 to 10.25m
FIT_WIN = struct( ...
    'radar_big',   [10.75, 12.75], ...
    'radar_small', [ 6.75,  9.75], ...
    'lidar_big',   [ 7.25, 13.25], ...
    'lidar_small', [ 8.75, 10.25]);

% =========================================================================
% SECTION 1: LOAD DATA
% =========================================================================
fprintf('=== SECTION 1: LOAD DATA ===\n');

% ── Radar ─────────────────────────────────────────────────────────────────
T   = readtable(RADAR_FILE);
m   = T.range <= RADAR_RANGE_MAX & T.amplitude >= RADAR_AMP_MIN;
r_r = T.range(m);
amp = T.amplitude(m);
X_r =  r_r .* sin(-T.pitch(m)) .* cos(T.roll(m));
Y_r =  r_r .* sin(-T.pitch(m)) .* sin(T.roll(m));
Z_r = -r_r .* cos(-T.pitch(m));
pts_r = [X_r, Y_r, Z_r];
fprintf('  Radar: %d returns (within %dm, above %ddB)\n', size(pts_r,1), RADAR_RANGE_MAX, RADAR_AMP_MIN);

% ── LiDAR big ─────────────────────────────────────────────────────────────
[pts_lb, refl_lb] = load_lidar(LIDAR_FILE_BIG, LIDAR_RANGE_MAX);
fprintf('  LiDAR big:   %d returns\n', size(pts_lb,1));

% ── LiDAR small ───────────────────────────────────────────────────────────
[pts_ls, refl_ls] = load_lidar(LIDAR_FILE_SMALL, LIDAR_RANGE_MAX);
fprintf('  LiDAR small: %d returns\n', size(pts_ls,1));

% Convert LiDAR reflectivity to dB scale: 10*log10(refl+1)
% +1 avoids log(0). This is a relative scale — only the slope matters.
refl_lb_dB = 10 * log10(double(refl_lb) + 1);
refl_ls_dB = 10 * log10(double(refl_ls) + 1);

% =========================================================================
% SECTION 2: CONE PROFILES
% =========================================================================
fprintf('\n=== SECTION 2: CONICAL SECTION EXTRACTION (cone half-angle = %d deg) ===\n', CONE_DEG);

[ctr_rb, mn_rb, sd_rb, cnt_rb] = cone_profile(pts_r,  amp,         REF_BIG_XYZ_RADAR,   CONE_DEG, BIN_W, MIN_PTS);
[ctr_rs, mn_rs, sd_rs, cnt_rs] = cone_profile(pts_r,  amp,         REF_SMALL_XYZ_RADAR, CONE_DEG, BIN_W, MIN_PTS);
[ctr_lb, mn_lb, sd_lb, cnt_lb] = cone_profile(pts_lb, refl_lb_dB,  REF_BIG_XYZ_LIDAR,   CONE_DEG, BIN_W, MIN_PTS);
[ctr_ls, mn_ls, sd_ls, cnt_ls] = cone_profile(pts_ls, refl_ls_dB,  REF_SMALL_XYZ_LIDAR, CONE_DEG, BIN_W, MIN_PTS);

fprintf('  Radar big cone:    %d non-empty bins,  %d total points\n', sum(isfinite(mn_rb)), sum(cnt_rb));
fprintf('  Radar small cone:  %d non-empty bins,  %d total points\n', sum(isfinite(mn_rs)), sum(cnt_rs));
fprintf('  LiDAR big cone:    %d non-empty bins,  %d total points\n', sum(isfinite(mn_lb)), sum(cnt_lb));
fprintf('  LiDAR small cone:  %d non-empty bins,  %d total points\n', sum(isfinite(mn_ls)), sum(cnt_ls));

% =========================================================================
% SECTION 3: FIT ATTENUATION COEFFICIENTS
% =========================================================================
fprintf('\n=== SECTION 3: ATTENUATION FIT ===\n');

[alpha_rb, fx_rb, fy_rb, r2_rb] = fit_alpha(ctr_rb, mn_rb, FIT_WIN.radar_big(1),   FIT_WIN.radar_big(2));
[alpha_rs, fx_rs, fy_rs, r2_rs] = fit_alpha(ctr_rs, mn_rs, FIT_WIN.radar_small(1), FIT_WIN.radar_small(2));
[alpha_lb, fx_lb, fy_lb, r2_lb] = fit_alpha(ctr_lb, mn_lb, FIT_WIN.lidar_big(1),   FIT_WIN.lidar_big(2));
[alpha_ls, fx_ls, fy_ls, r2_ls] = fit_alpha(ctr_ls, mn_ls, FIT_WIN.lidar_small(1), FIT_WIN.lidar_small(2));

alpha_r_mean = mean([alpha_rb, alpha_rs], 'omitnan');
alpha_l_mean = mean([alpha_lb, alpha_ls], 'omitnan');

% =========================================================================
% SECTION 4: SUMMARY TABLE
% =========================================================================
fprintf('\n=== SECTION 4: RESULTS SUMMARY ===\n\n');
fprintf('  Sensor    Cone     Fit window (m)    alpha (dB/m)    R²\n');
fprintf('  ────────  ───────  ────────────────  ─────────────   ────\n');
fprintf('  Radar     Big      %4.2f – %4.2fm     %+.3f           %.3f\n', FIT_WIN.radar_big(1),   FIT_WIN.radar_big(2),   alpha_rb, r2_rb);
fprintf('  Radar     Small    %4.2f – %4.2fm     %+.3f           %.3f   [noisy — use with caution]\n', FIT_WIN.radar_small(1), FIT_WIN.radar_small(2), alpha_rs, r2_rs);
fprintf('  LiDAR     Big      %4.2f – %4.2fm     %+.3f           %.3f\n', FIT_WIN.lidar_big(1),   FIT_WIN.lidar_big(2),   alpha_lb, r2_lb);
fprintf('  LiDAR     Small    %4.2f – %4.2fm     %+.3f           %.3f\n', FIT_WIN.lidar_small(1), FIT_WIN.lidar_small(2), alpha_ls, r2_ls);
fprintf('\n');
fprintf('  Radar mean alpha (big cone only, reliable): %+.3f dB/m\n', alpha_rb);
fprintf('  LiDAR mean alpha (both cones):              %+.3f dB/m\n', alpha_l_mean);
fprintf('\n');
fprintf('  Literature (120 GHz, moderate canopy): -1.5 to -4.0 dB/m\n');
fprintf('  Literature (LiDAR 905 nm, dense canopy): -5 to -15 dB/m\n');
fprintf('  Note: lower than lit. values consistent with open plane tree canopy.\n');

% =========================================================================
% SECTION 5: FIGURE 1 — Side-by-side cone profiles
% =========================================================================
fprintf('\n=== SECTION 5: PLOTS ===\n');

fig1 = figure('Name','Conical attenuation profiles','Color','white', ...
    'Units','normalized','Position',[0.02 0.35 0.95 0.52]);
tl1  = tiledlayout(fig1, 1, 2, 'TileSpacing','compact','Padding','compact');
sgtitle('120 GHz Radar vs LiDAR 905 nm - Conical Section Attenuation Profiles (Sunday scan)', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

% ── Panel 1: Big reflector ────────────────────────────────────────────────
ax1 = nexttile(tl1);
plot_panel(ax1, ...
    ctr_rb, mn_rb, sd_rb, fx_rb, fy_rb, alpha_rb, r2_rb, ...
    ctr_lb, mn_lb, sd_lb, fx_lb, fy_lb, alpha_lb, r2_lb, ...
    R_BIG, FIT_WIN.radar_big(1), FIT_WIN.lidar_big(1), ...
    sprintf('Big reflector cone  (GPS range = %.2f m)', R_BIG));

% ── Panel 2: Small reflector ──────────────────────────────────────────────
ax2 = nexttile(tl1);
plot_panel(ax2, ...
    ctr_rs, mn_rs, sd_rs, fx_rs, fy_rs, alpha_rs, r2_rs, ...
    ctr_ls, mn_ls, sd_ls, fx_ls, fy_ls, alpha_ls, r2_ls, ...
    R_SMALL, FIT_WIN.radar_small(1), FIT_WIN.lidar_small(1), ...
    sprintf('Small reflector cone  (GPS range = %.2f m)', R_SMALL));

% =========================================================================
% SECTION 6: FIGURE 2 — Normalised comparison (big cone)
% =========================================================================
fig2 = figure('Name','Normalised attenuation comparison','Color','white', ...
    'Units','normalized','Position',[0.02 0.05 0.60 0.45]);
ax3  = axes(fig2);
hold(ax3, 'on');

% Normalise both profiles to 0 dB at their first valid bin in the canopy region
mn_rb_n = mn_rb - nanmean(mn_rb(ctr_rb >= 7.5  & ctr_rb <= 8.5));
mn_lb_n = mn_lb - nanmean(mn_lb(ctr_lb >= 7.0  & ctr_lb <= 7.75));

vr = isfinite(mn_rb_n);
vl = isfinite(mn_lb_n);

% Shaded std bands
fill(ax3, [ctr_rb(vr), fliplr(ctr_rb(vr))], ...
    [mn_rb_n(vr)-sd_rb(vr), fliplr(mn_rb_n(vr)+sd_rb(vr))], ...
    [0.2 0.5 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'DisplayName', 'Radar +/- 1 std');
fill(ax3, [ctr_lb(vl), fliplr(ctr_lb(vl))], ...
    [mn_lb_n(vl)-sd_lb(vl), fliplr(mn_lb_n(vl)+sd_lb(vl))], ...
    [0.8 0.2 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
    'DisplayName', 'LiDAR +/- 1 std');

% Profiles
plot(ax3, ctr_rb(vr), mn_rb_n(vr), 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', 'DisplayName', '120 GHz radar (big cone, normalised)');
plot(ax3, ctr_lb(vl), mn_lb_n(vl), 'r-s', 'LineWidth', 2, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'LiDAR 905 nm (big cone, normalised)');

% Fit lines (shifted to normalised scale)
if ~isempty(fx_rb)
    shift_r = nanmean(mn_rb_n(ctr_rb >= fx_rb(1)-0.4 & ctr_rb <= fx_rb(1)+0.4)) - fy_rb(1);
    plot(ax3, fx_rb, fy_rb + shift_r, 'b--', 'LineWidth', 2, ...
        'DisplayName', sprintf('Radar fit: alpha = %.2f dB/m  (R2 = %.2f)', alpha_rb, r2_rb));
end
if ~isempty(fx_lb)
    shift_l = nanmean(mn_lb_n(ctr_lb >= fx_lb(1)-0.4 & ctr_lb <= fx_lb(1)+0.4)) - fy_lb(1);
    plot(ax3, fx_lb, fy_lb + shift_l, 'r--', 'LineWidth', 2, ...
        'DisplayName', sprintf('LiDAR fit: alpha = %.2f dB/m  (R2 = %.2f)', alpha_lb, r2_lb));
end

xline(ax3, R_BIG,  'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Reflector range (%.2fm)', R_BIG));
yline(ax3, 0, 'Color', [0.6 0.6 0.6], 'LineStyle', ':', 'LineWidth', 0.8);

xlabel(ax3, 'Slant range from sensor (m)', 'FontSize', 11);
ylabel(ax3, 'Amplitude relative to canopy edge (dB)', 'FontSize', 11);
title(ax3, 'Radar vs LiDAR: normalised attenuation - big reflector cone', ...
    'Interpreter', 'none', ...
    'FontSize', 11, 'FontWeight', 'bold');
legend(ax3, 'Location', 'southwest', 'FontSize', 9, 'Interpreter', 'none');
grid(ax3, 'on'); xlim(ax3, [0, R_BIG * 1.35]);

% =========================================================================
% SECTION 7: FIGURE 3 — Point counts sanity check
% =========================================================================
fig3 = figure('Name','Point counts per bin','Color','white', ...
    'Units','normalized','Position',[0.63 0.05 0.35 0.45]);
tl3  = tiledlayout(fig3, 2, 2, 'TileSpacing','compact','Padding','compact');
title(tl3, 'Returns per range bin - sanity check', 'FontSize', 11, 'Interpreter', 'none');

datasets = { ctr_rb, cnt_rb, 'Radar - big cone',   [0.2 0.5 0.8];
             ctr_rs, cnt_rs, 'Radar - small cone',  [0.2 0.5 0.8];
             ctr_lb, cnt_lb, 'LiDAR - big cone',    [0.8 0.2 0.1];
             ctr_ls, cnt_ls, 'LiDAR - small cone',  [0.8 0.2 0.1] };

for k = 1:4
    ax = nexttile(tl3);
    bar(ax, datasets{k,1}, datasets{k,2}, 1.0, ...
        'FaceColor', datasets{k,4}, 'FaceAlpha', 0.65, 'EdgeColor', 'none');
    yline(ax, MIN_PTS, 'r:', 'LineWidth', 1.2);
    xlabel(ax, 'Range (m)', 'FontSize', 9);
    ylabel(ax, 'Returns', 'FontSize', 9);
    title(ax, datasets{k,3}, 'FontSize', 9, 'Interpreter', 'none');
    grid(ax, 'on');
end

fprintf('  3 figures generated.\n');
fprintf('\n=== DONE ===\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [pts, refl] = load_lidar(fname, rmax)
    T    = readtable(fname, 'VariableNamingRule', 'preserve');
    get  = @(n) T{:, find(strcmpi(T.Properties.VariableNames, n), 1)};
    X = get('X'); Y = get('Y'); Z = get('Z'); refl = get('Reflectivity');
    % Remove invalid zero returns
    valid = ~(X == 0 & Y == 0 & Z == 0);
    X=X(valid); Y=Y(valid); Z=Z(valid); refl=refl(valid);
    % Range filter
    r    = sqrt(X.^2 + Y.^2 + Z.^2);
    keep = r <= rmax;
    X=X(keep); Y=Y(keep); Z=Z(keep); refl=refl(keep);
    pts  = [X, Y, Z];
end

function [ctr, mn, sd, cnt] = cone_profile(pts, values, ref_xyz, cone_deg, bin_w, min_pts)
    % Compute axis unit vector toward reflector
    axis_vec = ref_xyz / norm(ref_xyz);
    % Compute angular offset of each point from cone axis
    pt_norms = vecnorm(pts, 2, 2);
    valid    = pt_norms > 0.1;   % exclude origin
    cos_ang  = min(1, max(-1, (pts(valid,:) * axis_vec') ./ pt_norms(valid)));
    in_cone  = rad2deg(acos(cos_ang)) <= cone_deg;
    % Extract cone points and values
    pts_c  = pts(valid,:);  pts_c  = pts_c(in_cone,:);
    val_c  = values(valid); val_c  = val_c(in_cone);
    r_c    = vecnorm(pts_c, 2, 2);
    % Bin by slant range
    r_max  = norm(ref_xyz) * 1.5;
    edges  = (0 : bin_w : r_max)';
    nb     = numel(edges) - 1;
    ctr    = edges(1:end-1) + bin_w/2;
    mn     = nan(nb, 1); sd = nan(nb, 1); cnt = zeros(nb, 1);
    for i = 1:nb
        m = r_c >= edges(i) & r_c < edges(i+1);
        if sum(m) >= min_pts
            mn(i)  = mean(val_c(m));
            sd(i)  = std(val_c(m));
            cnt(i) = sum(m);
        end
    end
    ctr = ctr'; mn = mn'; sd = sd'; cnt = cnt';
end

function [alpha, fx, fy, r2] = fit_alpha(ctr, mn, r_start, r_end)
    mask = isfinite(mn) & ctr >= r_start & ctr <= r_end;
    if sum(mask) < 3
        alpha = NaN; fx = []; fy = []; r2 = NaN;
        fprintf('    WARNING: fewer than 3 valid bins in fit window [%.2f, %.2f]\n', r_start, r_end);
        return;
    end
    x  = ctr(mask)';
    y  = mn(mask)';
    p  = polyfit(x, y, 1);
    alpha  = p(1);
    fx     = x;
    fy     = polyval(p, x);
    ss_res = sum((y - fy).^2);
    ss_tot = sum((y - mean(y)).^2);
    r2     = 1 - ss_res / max(ss_tot, eps);
end

function plot_panel(ax, ctr_r, mn_r, sd_r, fx_r, fy_r, alpha_r, r2_r, ...
                        ctr_l, mn_l, sd_l, fx_l, fy_l, alpha_l, r2_l, ...
                        r_ref, rs_radar, rs_lidar, ttl)
    hold(ax, 'on');
    % Std bands
    vr = isfinite(mn_r) & isfinite(sd_r);
    if any(vr)
        fill(ax, [ctr_r(vr), fliplr(ctr_r(vr))], ...
            [mn_r(vr)-sd_r(vr), fliplr(mn_r(vr)+sd_r(vr))], ...
            [0.2 0.5 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'DisplayName', 'Radar +/- 1 std');
    end
    vl = isfinite(mn_l) & isfinite(sd_l);
    if any(vl)
        fill(ax, [ctr_l(vl), fliplr(ctr_l(vl))], ...
            [mn_l(vl)-sd_l(vl), fliplr(mn_l(vl)+sd_l(vl))], ...
            [0.8 0.2 0.1], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
            'DisplayName', 'LiDAR +/- 1 std');
    end
    % Profiles
    vr2 = isfinite(mn_r); vl2 = isfinite(mn_l);
    plot(ax, ctr_r(vr2), mn_r(vr2), 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'b', 'DisplayName', '120 GHz radar (mean/bin)');
    plot(ax, ctr_l(vl2), mn_l(vl2), 'r-s', 'LineWidth', 1.8, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'LiDAR 905 nm (mean/bin)');
    % Fit lines
    if ~isempty(fx_r)
        plot(ax, fx_r, fy_r, 'b--', 'LineWidth', 2, ...
            'DisplayName', sprintf('Radar fit: alpha = %.2f dB/m  (R2 = %.2f)', alpha_r, r2_r));
    end
    if ~isempty(fx_l)
        plot(ax, fx_l, fy_l, 'r--', 'LineWidth', 2, ...
            'DisplayName', sprintf('LiDAR fit: alpha = %.2f dB/m  (R2 = %.2f)', alpha_l, r2_l));
    end
    % Reference lines
    xline(ax, r_ref,    'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Reflector range (%.2fm)', r_ref));
    xline(ax, rs_radar, 'Color', [0.2 0.5 0.2], 'LineStyle', '--', 'LineWidth', 0.9, ...
        'DisplayName', 'Radar fit window start');
    xline(ax, rs_lidar, 'Color', [0.8 0.3 0.1], 'LineStyle', '--', 'LineWidth', 0.9, ...
        'DisplayName', 'LiDAR fit window start');
    xlabel(ax, 'Slant range from sensor (m)', 'FontSize', 10);
    ylabel(ax, 'Mean amplitude / reflectivity (dB)', 'FontSize', 10);
    title(ax, ttl, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    legend(ax, 'Location', 'southwest', 'FontSize', 7.5, 'Interpreter', 'none');
    grid(ax, 'on');
    xlim(ax, [0, r_ref * 1.35]);
end