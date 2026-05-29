clear; close all;

% ═══════════════════════════════════════════════════════════════════════════
%  REEDS SCAN ANALYSIS — April vs March
%  Compares radar detection of corner reflectors behind reed vegetation
% ═══════════════════════════════════════════════════════════════════════════

% ── LOAD DATA ─────────────────────────────────────────────────────────────
data_apr_raw = readmatrix('pointcloud_20260408-131953_.csv');
data_mar_raw = readmatrix('pointcloud_20260311-114835_.csv');

% ── FILTERS ───────────────────────────────────────────────────────────────
data_apr = data_apr_raw(data_apr_raw(:,5) >= -80, :);
data_apr = data_apr(data_apr(:,2) <= 22, :);

data_mar = data_mar_raw(data_mar_raw(:,5) >= -80, :);
data_mar = data_mar(data_mar(:,2) <= 15, :);

% ── POINT COUNT SUMMARY ───────────────────────────────────────────────────
fprintf('\n═══════════════════════════════════════════════\n');
fprintf(' REEDS SCAN ANALYSIS — April vs March\n');
fprintf('═══════════════════════════════════════════════\n\n');

fprintf('── Point Count Summary ─────────────────────────\n');
fprintf('  %-30s Raw: %5d   Filtered: %5d   Removed: %.1f%%\n', ...
    'April (20260408-131953)', size(data_apr_raw,1), size(data_apr,1), ...
    100*(size(data_apr_raw,1)-size(data_apr,1))/size(data_apr_raw,1));
fprintf('  %-30s Raw: %5d   Filtered: %5d   Removed: %.1f%%\n', ...
    'March (20260311-114835)', size(data_mar_raw,1), size(data_mar,1), ...
    100*(size(data_mar_raw,1)-size(data_mar,1))/size(data_mar_raw,1));
fprintf('═══════════════════════════════════════════════\n\n');

% ── CARTESIAN CONVERSION ──────────────────────────────────────────────────
az_apr = data_apr(:,4); el_apr = (pi/2) + data_apr(:,3); r_apr = data_apr(:,2);
[x_apr, y_apr, z_apr] = sph2cart(az_apr, el_apr, r_apr);

az_mar = data_mar(:,4); el_mar = (pi/2) + data_mar(:,3); r_mar = data_mar(:,2);
[x_mar, y_mar, z_mar] = sph2cart(az_mar, el_mar, r_mar);

% ── CORNER REFLECTOR MASKS ────────────────────────────────────────────────

% April CR1 — strong return at ~13.88m
apr_cr1_mask = data_apr(:,5) >= 20 & data_apr(:,2) >= 13 & data_apr(:,2) <= 15 & ...
               x_apr >= 12 & x_apr <= 14 & y_apr >= 4.5 & y_apr <= 6;

% April CR2 — weak attenuated return behind reeds
apr_cr2_mask = data_apr(:,5) >= -30 & data_apr(:,2) >= 12 & data_apr(:,2) <= 15 & ...
               x_apr >= 10 & x_apr <= 12 & y_apr >= 8 & y_apr <= 10;

% March CR1 — at ~10.9m, Y > 3.3
mar_cr1_mask = data_mar(:,5) >= 15 & data_mar(:,2) >= 10 & data_mar(:,2) <= 12 & ...
               y_mar >= 3.3 & y_mar <= 4.0;

% March CR2 — at ~10.9m, Y < 3.3
mar_cr2_mask = data_mar(:,5) >= 15 & data_mar(:,2) >= 10 & data_mar(:,2) <= 12 & ...
               y_mar >= 2.5 & y_mar <= 3.3;

% ── CLUTTER FLOOR ESTIMATE ────────────────────────────────────────────────
% Background clutter = all points excluding CR regions
apr_amp     = data_apr(:,5);
mar_amp     = data_mar(:,5);
apr_clutter = apr_amp(~apr_cr1_mask & ~apr_cr2_mask);
mar_clutter = mar_amp(~mar_cr1_mask & ~mar_cr2_mask);

apr_clutter_mean = mean(apr_clutter);
apr_clutter_std  = std(apr_clutter);
mar_clutter_mean = mean(mar_clutter);
mar_clutter_std  = std(mar_clutter);

% ── SNR CALCULATION ───────────────────────────────────────────────────────
apr_cr1_snr = max(data_apr(apr_cr1_mask,5)) - apr_clutter_mean;
apr_cr2_snr = max(data_apr(apr_cr2_mask,5)) - apr_clutter_mean;
mar_cr1_snr = max(data_mar(mar_cr1_mask,5)) - mar_clutter_mean;
mar_cr2_snr = max(data_mar(mar_cr2_mask,5)) - mar_clutter_mean;

% ── PRINT FULL ANALYSIS TABLE ─────────────────────────────────────────────
fprintf('── Corner Reflector Detection Summary ──────────\n\n');
fprintf('%-20s %-12s %-12s %-12s %-12s\n', '', 'April CR1', 'April CR2', 'March CR1', 'March CR2');
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%-20s %-12.2f %-12.2f %-12.2f %-12.2f\n', 'Range (m)', ...
    mean(r_apr(apr_cr1_mask)), mean(r_apr(apr_cr2_mask)), ...
    mean(r_mar(mar_cr1_mask)), mean(r_mar(mar_cr2_mask)));
fprintf('%-20s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'Peak amp (dB)', ...
    max(data_apr(apr_cr1_mask,5)), max(data_apr(apr_cr2_mask,5)), ...
    max(data_mar(mar_cr1_mask,5)), max(data_mar(mar_cr2_mask,5)));
fprintf('%-20s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'Mean amp (dB)', ...
    mean(data_apr(apr_cr1_mask,5)), mean(data_apr(apr_cr2_mask,5)), ...
    mean(data_mar(mar_cr1_mask,5)), mean(data_mar(mar_cr2_mask,5)));
fprintf('%-20s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'SNR above clutter', ...
    apr_cr1_snr, apr_cr2_snr, mar_cr1_snr, mar_cr2_snr);
fprintf('%-20s %-12d %-12d %-12d %-12d\n', 'Points in cluster', ...
    sum(apr_cr1_mask), sum(apr_cr2_mask), sum(mar_cr1_mask), sum(mar_cr2_mask));
fprintf('%-20s %-12s %-12s %-12s %-12s\n', 'Detected?', ...
    'YES', 'PARTIAL', 'YES', 'YES');

fprintf('\n── Clutter Floor ───────────────────────────────\n');
fprintf('  April : mean = %.1f dB,  std = %.1f dB\n', apr_clutter_mean, apr_clutter_std);
fprintf('  March : mean = %.1f dB,  std = %.1f dB\n', mar_clutter_mean, mar_clutter_std);
fprintf('═══════════════════════════════════════════════\n\n');

% ── FIGURE 1: April 3D coloured by amplitude ──────────────────────────────
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x_apr, y_apr, z_apr, 5, data_apr(:,5), 'filled'); hold on;
plot3(mean(x_apr(apr_cr1_mask)), mean(y_apr(apr_cr1_mask)), mean(z_apr(apr_cr1_mask)), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot3(mean(x_apr(apr_cr2_mask)), mean(y_apr(apr_cr2_mask)), mean(z_apr(apr_cr2_mask)), ...
      'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'LineWidth', 2);
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('March Reeds Scan — Coloured by Amplitude');
legend('Returns', 'CR1 (~13.88m)', 'CR2 (~13.88m)');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'april_reeds_amplitude.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 2: March 3D coloured by amplitude ──────────────────────────────
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x_mar, y_mar, z_mar, 5, data_mar(:,5), 'filled'); hold on;
plot3(mean(x_mar(mar_cr1_mask)), mean(y_mar(mar_cr1_mask)), mean(z_mar(mar_cr1_mask)), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot3(mean(x_mar(mar_cr2_mask)), mean(y_mar(mar_cr2_mask)), mean(z_mar(mar_cr2_mask)), ...
      'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'LineWidth', 2);
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('April Reeds Scan — Coloured by Amplitude');
legend('Returns', 'CR1 (~10.9m)', 'CR2 (~10.9m)');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'march_reeds_amplitude.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 3: Amplitude vs Range — both scans overlaid ────────────────────
figure;
set(gcf, 'Position', [100 100 900 450]);
scatter(r_apr, data_apr(:,5), 2, 'b', 'filled'); hold on;
scatter(r_mar, data_mar(:,5), 2, 'r', 'filled');
yline(apr_clutter_mean, 'b--', 'April clutter floor', 'LineWidth', 1.5);
yline(mar_clutter_mean, 'r--', 'March clutter floor', 'LineWidth', 1.5);
plot(mean(r_apr(apr_cr1_mask)), max(data_apr(apr_cr1_mask,5)), 'b^', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'April CR1');
plot(mean(r_apr(apr_cr2_mask)), max(data_apr(apr_cr2_mask,5)), 'bs', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'April CR2');
plot(mean(r_mar(mar_cr1_mask)), max(data_mar(mar_cr1_mask,5)), 'r^', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'March CR1');
plot(mean(r_mar(mar_cr2_mask)), max(data_mar(mar_cr2_mask,5)), 'rs', ...
     'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'March CR2');
xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('Amplitude vs Range — April vs March Reeds');
legend('April returns', 'March returns', 'April clutter floor', 'March clutter floor', ...
       'April CR1', 'April CR2', 'March CR1', 'March CR2', 'Location', 'southwest');
grid on;
exportgraphics(gcf, 'reeds_amp_vs_range.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 4: Amplitude histograms side by side ──────────────────────────
figure;
set(gcf, 'Position', [100 100 900 400]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile;
histogram(data_mar(:,5), 40, 'FaceColor', [0.9 0.3 0.2], 'EdgeColor', 'none');
hold on;
xline(mar_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 2);
xline(max(data_mar(mar_cr1_mask,5)), 'r-', 'CR1 peak', 'LineWidth', 2);
xline(max(data_mar(mar_cr2_mask,5)), 'b-', 'CR2 peak', 'LineWidth', 2);
xlabel('Amplitude (dB)'); ylabel('Count');
title('March Reeds — Amplitude Distribution');
grid on;

nexttile;
histogram(data_apr(:,5), 40, 'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'none');
hold on;
xline(apr_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 2);
xline(max(data_apr(apr_cr1_mask,5)), 'r-', 'CR1 peak', 'LineWidth', 2);
xline(max(data_apr(apr_cr2_mask,5)), 'b-', 'CR2 peak', 'LineWidth', 2);
xlabel('Amplitude (dB)'); ylabel('Count');
title('April Reeds — Amplitude Distribution');
grid on;

% ── FIGURE 5: SNR bar chart comparison ────────────────────────────────────
figure;
set(gcf, 'Position', [100 100 700 400]);
snr_vals = [mar_cr1_snr, mar_cr2_snr; apr_cr1_snr, apr_cr2_snr];
b = bar(snr_vals);
b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
set(gca, 'XTickLabel', {'March', 'April'});
ylabel('SNR above clutter floor (dB)');
title('Corner Reflector SNR — March vs April Reeds');
legend('CR1', 'CR2');
grid on;
exportgraphics(gcf, 'reeds_snr_comparison.pdf', 'ContentType', 'image', 'Resolution', 300);