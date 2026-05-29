% =========================================================================
% align_radar_scans.m
% Aligns two radar scans into a common GPS world frame using corner
% reflectors as tie points. Verification prints after every section.
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
file_scan1 = 'pointcloud_20260525-125755_0_fft_cfar.csv';


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


fprintf('  Scan 1 raw rows: %d  |  columns: %s\n', height(T1), strjoin(T1.Properties.VariableNames, ', '));

fprintf('  Scan 1 range: %.2f to %.2f m\n', min(T1.range), max(T1.range));

fprintf('  Scan 1 amplitude: %.1f to %.1f dB\n', min(T1.peak), max(T1.peak));


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
amp1_all = T1.peak;



% Apply range filter
m1 = T1.range <= RANGE_LIMIT;


X1 = X1_all(m1); Y1 = Y1_all(m1); Z1 = Z1_all(m1); amp1 = amp1_all(m1);


fprintf('  Scan 1: %d points after %.0fm range filter\n', sum(m1), RANGE_LIMIT);

fprintf('  Scan 1 Cartesian X: %.2f to %.2f m\n', min(X1), max(X1));
fprintf('  Scan 1 Cartesian Y: %.2f to %.2f m\n', min(Y1), max(Y1));
fprintf('  Scan 1 Cartesian Z: %.2f to %.2f m\n', min(Z1), max(Z1));


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