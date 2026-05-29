% Cherry Ballart — Extended Analysis
% Methods 1, 2, 3
clear; close all;

% ==========================================================
% CONFIG
% ==========================================================
CSV_PATH = 'pointcloud_20260411-143028_.csv';

CR_X = 3.449; CR_Y = 1.560; CR_Z = 0.032;

AMP_MIN        = -80;
R_MAX          = 6.0;
NEAR_FIELD_CUT = 1.45;
FENCE_ROLL_DEG = 7.0;
CONE_HALF_ANGLE = 8.0;
BIN_WIDTH      = 0.25;
MIN_PTS_BIN    = 8;

SHRUB_FRONT  = 1.5;
SHRUB_BACK   = 3.4;
FOLIAGE_BACK = 2.5;
CR_TOL       = 0.15;
CR_RADIUS    = 0.18;

CR_RANGE = sqrt(CR_X^2 + CR_Y^2 + CR_Z^2);
VEG_DEPTH = SHRUB_BACK - SHRUB_FRONT;  % 1.9 m one-way

% ==========================================================
% LOAD + FILTER
% ==========================================================
data  = readmatrix(CSV_PATH);
r     = data(:,2);
pitch = data(:,3);
roll  = data(:,4);
amp   = data(:,5);

mask  = amp >= AMP_MIN & r <= R_MAX & r >= NEAR_FIELD_CUT & rad2deg(roll) >= FENCE_ROLL_DEG;
r     = r(mask); pitch = pitch(mask);
roll  = roll(mask); amp = amp(mask);

az = roll;
el = (pi/2) + pitch;
x  = r .* cos(el) .* cos(az);
y  = r .* cos(el) .* sin(az);
z  = r .* sin(el);

% ==========================================================
% CONE FILTER
% ==========================================================
cr_vec    = [CR_X, CR_Y, CR_Z];
boresight = cr_vec / norm(cr_vec);
pts       = [x, y, z];
pts_norm  = pts ./ vecnorm(pts, 2, 2);
cos_theta = min(max(pts_norm * boresight', -1), 1);
theta_deg = rad2deg(acos(cos_theta));
cone_mask = theta_deg <= CONE_HALF_ANGLE;

r_cone   = r(cone_mask);
amp_cone = amp(cone_mask);
x_cone   = x(cone_mask);
y_cone   = y(cone_mask);
z_cone   = z(cone_mask);

% ==========================================================
% MASKS
% ==========================================================
xyz_dist = sqrt((x_cone-CR_X).^2 + (y_cone-CR_Y).^2 + (z_cone-CR_Z).^2);

cr_mask = (r_cone >= CR_RANGE-CR_TOL) & ...
          (r_cone <= CR_RANGE+CR_TOL) & ...
          (xyz_dist <= CR_RADIUS);

near_mask = (r_cone >= SHRUB_FRONT) & ...
            (r_cone <= FOLIAGE_BACK) & ...
            (~cr_mask);

bush_mask = (r_cone >= SHRUB_FRONT) & ...
            (r_cone <= SHRUB_BACK)  & ...
            (~cr_mask);

cr_peak     = max(amp_cone(cr_mask));
near_mean   = mean(amp_cone(near_mask));
clutter_med = median(amp_cone(bush_mask));

% ==========================================================
% METHOD 1 — Near face mean vs CR peak (path loss)
% Two-way path = 2 * VEG_DEPTH
% ==========================================================
total_loss     = near_mean - cr_peak;
two_way_depth  = 2 * VEG_DEPTH;
loss_per_metre = total_loss / two_way_depth;

fprintf('\n========================================\n');
fprintf('METHOD 1 — Vegetation path loss\n');
fprintf('========================================\n');
fprintf('Near face mean (%.1f-%.1f m): %+.1f dB\n', SHRUB_FRONT, FOLIAGE_BACK, near_mean);
fprintf('CR peak:                      %+.1f dB\n', cr_peak);
fprintf('Total two-way loss:           %+.1f dB\n', total_loss);
fprintf('Two-way vegetation depth:     %.1f m\n', two_way_depth);
fprintf('Loss per metre:               %+.2f dB/m\n', loss_per_metre);

% ==========================================================
% HELPER — bin stats
% ==========================================================
function tbl = bin_stats(r_in, amp_in, r_lo, r_hi, bw, min_pts)
    edges = r_lo:bw:r_hi;
    tbl   = [];
    for i = 1:length(edges)-1
        sel = r_in >= edges(i) & r_in < edges(i)+bw;
        if sum(sel) >= min_pts
            tbl(end+1,:) = [edges(i)+bw/2, sum(sel), ...
                            mean(amp_in(sel)), std(amp_in(sel))];
        end
    end
end

% ==========================================================
% METHOD 2 — Point cloud density vs depth
% ==========================================================
% Use all cone returns, bin by range, count points per bin
% No minimum point threshold here — we want to see drop-off
edges_density = SHRUB_FRONT:BIN_WIDTH:R_MAX;
density_count = zeros(length(edges_density)-1, 1);
density_rmid  = zeros(length(edges_density)-1, 1);

for i = 1:length(edges_density)-1
    sel = r_cone >= edges_density(i) & r_cone < edges_density(i)+BIN_WIDTH & ~cr_mask;
    density_count(i) = sum(sel);
    density_rmid(i)  = edges_density(i) + BIN_WIDTH/2;
end

% Normalise by max count for easier comparison
density_norm = density_count / max(density_count);

fprintf('\n========================================\n');
fprintf('METHOD 2 — Point density vs depth\n');
fprintf('========================================\n');
fprintf('Range (m) | Count | Normalised\n');
for i = 1:length(density_rmid)
    fprintf('  %.2f    |  %3d  |  %.3f\n', density_rmid(i), density_count(i), density_norm(i));
end

% ==========================================================
% METHOD 3 — Clutter variance vs depth
% ==========================================================
prof = bin_stats(r_cone(bush_mask), amp_cone(bush_mask), ...
    SHRUB_FRONT, SHRUB_BACK+CR_TOL, BIN_WIDTH, 3);

fprintf('\n========================================\n');
fprintf('METHOD 3 — Amplitude variance vs depth\n');
fprintf('========================================\n');
fprintf('Range mid (m) | N pts | Mean (dB) | Std (dB)\n');
for i = 1:size(prof,1)
    fprintf('    %.2f       |  %3d  |  %+.1f    |  %.1f\n', ...
        prof(i,1), prof(i,2), prof(i,3), prof(i,4));
end

% ==========================================================
% FIGURES
% ==========================================================

% Figure 1 — Method 1: near face vs CR
figure;
bar_vals  = [near_mean, cr_peak];
bar_names = {'Near face mean (1.5-2.5m)', 'CR peak (3.785m)'};
b = bar(bar_vals, 'FaceColor', 'flat');
b.CData(1,:) = [0.2 0.6 0.2];
b.CData(2,:) = [0.8 0.2 0.2];
set(gca, 'XTickLabel', bar_names);
ylabel('Amplitude (dB)');
title(sprintf('Method 1: Vegetation Path Loss\nTotal loss = %.1f dB,  %.2f dB/m over %.1f m two-way', ...
    total_loss, loss_per_metre, two_way_depth));
grid on;
text(1, near_mean+0.5, sprintf('%+.1f dB', near_mean), 'HorizontalAlignment','center', 'FontSize', 11);
text(2, cr_peak+0.5,   sprintf('%+.1f dB', cr_peak),   'HorizontalAlignment','center', 'FontSize', 11);

% Figure 2 — Method 2: density vs depth
figure;
bar(density_rmid, density_count, BIN_WIDTH*0.85, 'FaceColor', [0.2 0.6 0.4], 'EdgeColor', 'k');
hold on;
xline(SHRUB_FRONT,  'g--', sprintf('Bush front (%.1f m)', SHRUB_FRONT),  'LineWidth', 1.5);
xline(FOLIAGE_BACK, 'b--', sprintf('Foliage end (%.1f m)', FOLIAGE_BACK),'LineWidth', 1.5);
xline(SHRUB_BACK,   'r--', sprintf('Bush back (%.1f m)', SHRUB_BACK),    'LineWidth', 1.5);
xline(CR_RANGE,     'r-',  sprintf('CR (%.2f m)', CR_RANGE),             'LineWidth', 1.5);
xlabel('Range (m)'); ylabel('Point count in cone');
title('Method 2: Cone Return Density vs Depth');
legend({'Point count', 'Bush front', 'Foliage end', 'Bush back', 'CR'}, 'Location', 'northeast');
grid on;

% Figure 3 — Method 3: variance vs depth
figure;
if ~isempty(prof)
    yyaxis left;
    plot(prof(:,1), prof(:,3), 'b-o', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', 'Mean amplitude');
    fill([prof(:,1); flipud(prof(:,1))], ...
         [prof(:,3)+prof(:,4); flipud(prof(:,3)-prof(:,4))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    ylabel('Mean amplitude (dB)');

    yyaxis right;
    plot(prof(:,1), prof(:,4), 'r-s', 'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', 'Std dev (variance)');
    ylabel('Standard deviation (dB)');

    hold on;
    xline(SHRUB_FRONT,  'g--', 'Bush front',  'LineWidth', 1.2);
    xline(FOLIAGE_BACK, 'b--', 'Foliage end', 'LineWidth', 1.2);
    xline(SHRUB_BACK,   'r--', 'Bush back',   'LineWidth', 1.2);

    xlabel('Range (m)');
    title('Method 3: Amplitude Mean and Variance vs Depth');
    legend({'Mean amplitude', '±1σ envelope', 'Std dev'}, 'Location', 'best');
    grid on;
end