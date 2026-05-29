% =========================================================================
% analyse_lidar_scans.m
%
% Loads and analyses two LiDAR point cloud CSV files (Livox format).
% Compares point density, reflectivity, height distribution, and
% identifies corner reflector returns in each scan.
%
% Set CSV_FILES below or leave empty to use a file picker dialog.
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
CSV_FILES = {'2026-03-15 14-52-14_corner_small.csv', ...
             '2026-03-15 14-50-46_corner_big.csv'};

POINT_SIZE    = 2;       % Scatter point size
MAX_RANGE     = 40;      % Discard returns beyond this distance (m)
FILTER_ZEROS  = true;    % Remove X=Y=Z=0 invalid returns
REFL_THRESH   = 200;     % Reflectivity threshold for corner reflector detection
COLOR_BY      = 'reflectivity';   % 'reflectivity' | 'z' | 'file'

SCAN_LABELS   = {'Scan 1', 'Scan 2'};
SCAN_COLORS   = {[0.2 0.5 1.0], [0.2 0.85 0.4]};   % per-scan colours for file mode

% -------------------------------------------------------------------------

%% ── Load files ───────────────────────────────────────────────────────────
if numel(CSV_FILES) ~= 2
    error('Please specify exactly 2 CSV files in CSV_FILES.');
end

% Storage for each scan
scans = struct('x', {}, 'y', {}, 'z', {}, 'refl', {}, 'label', {});

fprintf('=== SECTION 1: LOAD DATA ===\n');
for k = 1:2
    fprintf('  Loading: %s\n', CSV_FILES{k});

    T = readtable(CSV_FILES{k}, 'VariableNamingRule', 'preserve');

    % Flexible column lookup — works regardless of column order or case
    get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};

    x    = get_col('X');
    y    = get_col('Y');
    z    = get_col('Z');
    refl = get_col('Reflectivity');

    % Remove invalid zero returns
    if FILTER_ZEROS
        valid = ~(x == 0 & y == 0 & z == 0);
        x = x(valid); y = y(valid); z = z(valid); refl = refl(valid);
    end

    % Range filter
    range = sqrt(x.^2 + y.^2 + z.^2);
    keep  = range <= MAX_RANGE;
    x = x(keep); y = y(keep); z = z(keep); refl = refl(keep);

    scans(k).x     = x;
    scans(k).y     = y;
    scans(k).z     = z;
    scans(k).refl  = refl;
    scans(k).label = SCAN_LABELS{k};

    fprintf('    %d points loaded (after filtering, within %.0fm)\n', numel(x), MAX_RANGE);
    fprintf('    X: %.3f to %.3f m\n', min(x), max(x));
    fprintf('    Y: %.3f to %.3f m\n', min(y), max(y));
    fprintf('    Z: %.3f to %.3f m\n', min(z), max(z));
    fprintf('    Reflectivity: %d to %d  (mean=%.1f)\n', min(refl), max(refl), mean(refl));
end

%% ── Corner reflector detection ───────────────────────────────────────────
fprintf('\n=== SECTION 2: CORNER REFLECTOR DETECTION ===\n');
for k = 1:2
    mask = scans(k).refl >= REFL_THRESH;
    n    = sum(mask);
    fprintf('  %s: %d points above reflectivity threshold %d\n', ...
        scans(k).label, n, REFL_THRESH);
    if n > 0
        fprintf('    Centroid: X=%.3f  Y=%.3f  Z=%.3f m\n', ...
            mean(scans(k).x(mask)), mean(scans(k).y(mask)), mean(scans(k).z(mask)));
        fprintf('    Max reflectivity: %d\n', max(scans(k).refl(mask)));
    end
    scans(k).refl_mask = mask;
end

%% ── Colour data ──────────────────────────────────────────────────────────
all_x    = [scans(1).x;    scans(2).x];
all_y    = [scans(1).y;    scans(2).y];
all_z    = [scans(1).z;    scans(2).z];
all_refl = [scans(1).refl; scans(2).refl];
all_file = [ones(numel(scans(1).x),1); 2*ones(numel(scans(2).x),1)];

switch lower(COLOR_BY)
    case {'reflectivity','intensity'}
        cdata = all_refl; clabel = 'Reflectivity';
    case 'z'
        cdata = all_z;    clabel = 'Z (m)';
    case 'file'
        cdata = all_file; clabel = 'Scan index';
    otherwise
        cdata = all_refl; clabel = 'Reflectivity';
end

%% ── Figure 1: 3D overlay ─────────────────────────────────────────────────
fprintf('\n=== SECTION 3: PLOTS ===\n');

fig1 = figure('Name', 'LiDAR — 3D overlay', 'NumberTitle', 'off', ...
    'Color', 'k', 'Units', 'normalized', 'Position', [0.05 0.1 0.55 0.8]);
ax1  = axes('Parent', fig1, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
hold(ax1, 'on');

scatter3(ax1, all_x, all_y, all_z, POINT_SIZE, cdata, 'filled');

% Mark corner reflectors
for k = 1:2
    mask = scans(k).refl_mask;
    if any(mask)
        scatter3(ax1, scans(k).x(mask), scans(k).y(mask), scans(k).z(mask), ...
            40, 'filled', 'MarkerFaceColor', SCAN_COLORS{k}, ...
            'DisplayName', [scans(k).label ' reflector']);
    end
end

colormap(ax1, 'turbo');
cb = colorbar(ax1); cb.Color = 'w';
cb.Label.String = clabel; cb.Label.Color = 'w';
xlabel(ax1,'X (m)','Color','w'); ylabel(ax1,'Y (m)','Color','w');
zlabel(ax1,'Z (m)','Color','w');
title(ax1, sprintf('3D LiDAR overlay — range \\leq %.0fm', MAX_RANGE), ...
    'Color','w','FontSize',14);
axis(ax1,'equal'); grid(ax1,'on');
ax1.GridColor = [0.3 0.3 0.3];
view(ax1, 45, 30); rotate3d(ax1, 'on');
legend(ax1, 'show', 'TextColor', 'w', 'Color', 'k');

%% ── Figure 2: Top-down view ──────────────────────────────────────────────
fig2 = figure('Name', 'LiDAR — top-down (XY)', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.35 0.55 0.3 0.4]);
ax2  = axes('Parent', fig2, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
scatter(ax2, all_x, all_y, POINT_SIZE, cdata, 'filled');
colormap(ax2,'turbo');
cb2 = colorbar(ax2); cb2.Color = 'w';
cb2.Label.String = clabel; cb2.Label.Color = 'w';
axis(ax2,'equal'); grid(ax2,'on');
ax2.GridColor = [0.3 0.3 0.3];
xlabel(ax2,'X (m)','Color','w'); ylabel(ax2,'Y (m)','Color','w');
title(ax2,'Top-down (XY)','Color','w','FontSize',13);

%% ── Figure 3: Side view (XZ) ─────────────────────────────────────────────
fig3 = figure('Name', 'LiDAR — side view (XZ)', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.65 0.55 0.3 0.4]);
ax3  = axes('Parent', fig3, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
scatter(ax3, all_x, all_z, POINT_SIZE, cdata, 'filled');
colormap(ax3,'turbo');
cb3 = colorbar(ax3); cb3.Color = 'w';
cb3.Label.String = clabel; cb3.Label.Color = 'w';
axis(ax3,'equal'); grid(ax3,'on');
ax3.GridColor = [0.3 0.3 0.3];
xlabel(ax3,'X (m)','Color','w'); ylabel(ax3,'Z (m)','Color','w');
title(ax3,'Side view (XZ)','Color','w','FontSize',13);

%% ── Figure 4: Reflectivity histograms side by side ───────────────────────
fig4 = figure('Name', 'Reflectivity distributions', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.55 0.28 0.4]);

for k = 1:2
    ax = subplot(1, 2, k);
    set(ax, 'Color', 'w', 'XColor', 'w', 'YColor', 'w');
    histogram(ax, scans(k).refl, 64, 'FaceColor', SCAN_COLORS{k}, 'EdgeColor', 'none');
    xline(REFL_THRESH, 'r--', 'LineWidth', 1.5);
    xlabel(ax,'Reflectivity','Color','w');
    ylabel(ax,'Count','Color','w');
    title(ax, scans(k).label, 'Color', 'w', 'FontSize', 12);
    grid(ax,'on'); ax.GridColor = [0.3 0.3 0.3];
end
sgtitle('Reflectivity distributions', 'Color', 'w', 'FontSize', 13);

%% ── Figure 5: Height (Z) histograms side by side ─────────────────────────
fig5 = figure('Name', 'Height distributions', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.35 0.1 0.28 0.4]);

for k = 1:2
    ax = subplot(1, 2, k);
    set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    histogram(ax, scans(k).z, 60, 'FaceColor', SCAN_COLORS{k}, 'EdgeColor', 'none');
    xlabel(ax,'Z (m)','Color','w');
    ylabel(ax,'Count','Color','w');
    title(ax, scans(k).label, 'Color', 'w', 'FontSize', 12);
    grid(ax,'on'); ax.GridColor = [0.3 0.3 0.3];
end
sgtitle('Height (Z) distributions', 'Color', 'w', 'FontSize', 13);

fprintf('  5 figures generated.\n');
fprintf('\nDone. Use mouse to rotate the 3D view.\n');