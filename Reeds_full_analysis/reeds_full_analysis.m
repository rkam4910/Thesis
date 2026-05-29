clear; close all;

% ═══════════════════════════════════════════════════════════════════════════
%  REEDS SCAN ANALYSIS — April vs March
%  Produces all figures referenced in thesis LaTeX
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

% ── CARTESIAN CONVERSION ──────────────────────────────────────────────────
az_apr = data_apr(:,4); el_apr = (pi/2) + data_apr(:,3); r_apr = data_apr(:,2);
[x_apr, y_apr, z_apr] = sph2cart(az_apr, el_apr, r_apr);

az_mar = data_mar(:,4); el_mar = (pi/2) + data_mar(:,3); r_mar = data_mar(:,2);
[x_mar, y_mar, z_mar] = sph2cart(az_mar, el_mar, r_mar);

% ── CORNER REFLECTOR MASKS ────────────────────────────────────────────────
apr_cr1_mask = data_apr(:,5) >= 20 & data_apr(:,2) >= 13 & data_apr(:,2) <= 15 & ...
               x_apr >= 12 & x_apr <= 14 & y_apr >= 4.5 & y_apr <= 6;
apr_cr2_mask = data_apr(:,5) >= -30 & data_apr(:,2) >= 12 & data_apr(:,2) <= 15 & ...
               x_apr >= 10 & x_apr <= 12 & y_apr >= 8 & y_apr <= 10;
mar_cr1_mask = data_mar(:,5) >= 15 & data_mar(:,2) >= 10 & data_mar(:,2) <= 12 & ...
               y_mar >= 3.3 & y_mar <= 4.0;
mar_cr2_mask = data_mar(:,5) >= 15 & data_mar(:,2) >= 10 & data_mar(:,2) <= 12 & ...
               y_mar >= 2.5 & y_mar <= 3.3;

% ── CLUTTER FLOOR ─────────────────────────────────────────────────────────
apr_amp     = data_apr(:,5);
mar_amp     = data_mar(:,5);
apr_clutter = apr_amp(~apr_cr1_mask & ~apr_cr2_mask);
mar_clutter = mar_amp(~mar_cr1_mask & ~mar_cr2_mask);

apr_clutter_mean = mean(apr_clutter);
apr_clutter_std  = std(apr_clutter);
mar_clutter_mean = mean(mar_clutter);
mar_clutter_std  = std(mar_clutter);

% ── FIXED DETECTION THRESHOLD ─────────────────────────────────────────────
FIXED_THRESH = mean([apr_clutter_mean, mar_clutter_mean]);
fprintf('\nFixed detection threshold: %.1f dB (mean clutter floor)\n', FIXED_THRESH);

% ── SNR ───────────────────────────────────────────────────────────────────
apr_cr1_snr = max(data_apr(apr_cr1_mask,5)) - apr_clutter_mean;
apr_cr2_snr = max(data_apr(apr_cr2_mask,5)) - apr_clutter_mean;
mar_cr1_snr = max(data_mar(mar_cr1_mask,5)) - mar_clutter_mean;
mar_cr2_snr = max(data_mar(mar_cr2_mask,5)) - mar_clutter_mean;

% ── TWO-WAY ATTENUATION ESTIMATE ──────────────────────────────────────────
attenuation_cr1 = max(data_apr(apr_cr1_mask,5)) - max(data_mar(mar_cr1_mask,5));
fprintf('\n── Two-Way Attenuation Estimate ────────────────\n');
fprintf('  April CR1 peak : %.1f dB\n', max(data_apr(apr_cr1_mask,5)));
fprintf('  March CR1 peak : %.1f dB\n', max(data_mar(mar_cr1_mask,5)));
fprintf('  Two-way loss   : ~%.1f dB through reed layer\n', attenuation_cr1);
fprintf('  One-way loss   : ~%.1f dB through reed layer\n', attenuation_cr1/2);

% ── DETECTION TABLE ───────────────────────────────────────────────────────
fprintf('\n── Corner Reflector Detection Summary ──────────\n\n');
fprintf('%-22s %-12s %-12s %-12s %-12s\n', '', 'April CR1', 'April CR2', 'March CR1', 'March CR2');
fprintf('%s\n', repmat('-', 1, 72));
fprintf('%-22s %-12.2f %-12.2f %-12.2f %-12.2f\n', 'Range (m)', ...
    mean(r_apr(apr_cr1_mask)), mean(r_apr(apr_cr2_mask)), ...
    mean(r_mar(mar_cr1_mask)), mean(r_mar(mar_cr2_mask)));
fprintf('%-22s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'Peak amp (dB)', ...
    max(data_apr(apr_cr1_mask,5)), max(data_apr(apr_cr2_mask,5)), ...
    max(data_mar(mar_cr1_mask,5)), max(data_mar(mar_cr2_mask,5)));
fprintf('%-22s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'SNR above clutter (dB)', ...
    apr_cr1_snr, apr_cr2_snr, mar_cr1_snr, mar_cr2_snr);
fprintf('%-22s %-12.1f %-12.1f %-12.1f %-12.1f\n', 'Margin above thresh (dB)', ...
    max(data_apr(apr_cr1_mask,5)) - FIXED_THRESH, ...
    max(data_apr(apr_cr2_mask,5)) - FIXED_THRESH, ...
    max(data_mar(mar_cr1_mask,5)) - FIXED_THRESH, ...
    max(data_mar(mar_cr2_mask,5)) - FIXED_THRESH);
fprintf('%-22s %-12d %-12d %-12d %-12d\n', 'Points in cluster', ...
    sum(apr_cr1_mask), sum(apr_cr2_mask), sum(mar_cr1_mask), sum(mar_cr2_mask));
fprintf('%-22s %-12s %-12s %-12s %-12s\n', 'Detected?', ...
    'YES', 'YES (weak)', 'YES', 'YES');
fprintf('\n── Clutter Floor Reference ─────────────────────\n');
fprintf('  April : mean = %.1f dB,  std = %.1f dB\n', apr_clutter_mean, apr_clutter_std);
fprintf('  March : mean = %.1f dB,  std = %.1f dB\n', mar_clutter_mean, mar_clutter_std);
fprintf('  Fixed threshold : %.1f dB\n', FIXED_THRESH);
fprintf('═══════════════════════════════════════════════\n\n');

% ── RANGE PROFILES ────────────────────────────────────────────────────────
apr_slice = data_apr(:,4) >= 0.37 & data_apr(:,4) <= 0.40;
mar_slice = data_mar(:,4) >= 0.32 & data_mar(:,4) <= 0.35;

% ═══════════════════════════════════════════════════════════════════════════
%  FIGURES
% ═══════════════════════════════════════════════════════════════════════════

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
title('April Reeds Scan — Coloured by Amplitude');
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
title('March Reeds Scan — Coloured by Amplitude');
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
yline(FIXED_THRESH, 'k:', sprintf('Detection threshold (%.1f dB)', FIXED_THRESH), 'LineWidth', 2);
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
       'Detection threshold', 'April CR1', 'April CR2', 'March CR1', 'March CR2', ...
       'Location', 'southwest');
grid on;
exportgraphics(gcf, 'reeds_amp_vs_range.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 4: Amplitude histograms — March then April ─────────────────────
figure;
set(gcf, 'Position', [100 100 900 400]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile;
histogram(data_mar(:,5), 40, 'FaceColor', [0.9 0.3 0.2], 'EdgeColor', 'none');
hold on;
xline(mar_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 2);
xline(FIXED_THRESH, 'k:', 'Detection threshold', 'LineWidth', 2);
xline(max(data_mar(mar_cr1_mask,5)), 'r-', 'CR1 peak', 'LineWidth', 2);
xline(max(data_mar(mar_cr2_mask,5)), 'b-', 'CR2 peak', 'LineWidth', 2);
xlabel('Amplitude (dB)'); ylabel('Count');
title('March Reeds — Amplitude Distribution');
grid on;

nexttile;
histogram(data_apr(:,5), 40, 'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'none');
hold on;
xline(apr_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 2);
xline(FIXED_THRESH, 'k:', 'Detection threshold', 'LineWidth', 2);
xline(max(data_apr(apr_cr1_mask,5)), 'r-', 'CR1 peak', 'LineWidth', 2);
xline(max(data_apr(apr_cr2_mask,5)), 'b-', 'CR2 peak', 'LineWidth', 2);
xlabel('Amplitude (dB)'); ylabel('Count');
title('April Reeds — Amplitude Distribution');
grid on;

exportgraphics(gcf, 'reeds_amplitude_histograms.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 5: SNR bar chart — March then April ────────────────────────────
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

% ── FIGURE 6: Range profiles through CR azimuth slice ─────────────────────
figure;
set(gcf, 'Position', [100 100 900 400]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile;
r_slice_mar = data_mar(mar_slice, 2);
a_slice_mar = data_mar(mar_slice, 5);
scatter(r_slice_mar, a_slice_mar, 5, 'r', 'filled');
hold on;
yline(mar_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 1.5);
yline(FIXED_THRESH, 'k:', 'Detection threshold', 'LineWidth', 1.5);
xline(mean(r_mar(mar_cr1_mask)), 'r-', 'CR1', 'LineWidth', 2);
xline(mean(r_mar(mar_cr2_mask)), 'b-', 'CR2', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('March — Range Profile Through CR Azimuth');
grid on;

nexttile;
r_slice_apr = data_apr(apr_slice, 2);
a_slice_apr = data_apr(apr_slice, 5);
scatter(r_slice_apr, a_slice_apr, 5, 'b', 'filled');
hold on;
yline(apr_clutter_mean, 'k--', 'Clutter floor', 'LineWidth', 1.5);
yline(FIXED_THRESH, 'k:', 'Detection threshold', 'LineWidth', 1.5);
xline(mean(r_apr(apr_cr1_mask)), 'r-', 'CR1', 'LineWidth', 2);
xline(mean(r_apr(apr_cr2_mask)), 'b-', 'CR2', 'LineWidth', 2);
xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('April — Range Profile Through CR Azimuth');
grid on;

exportgraphics(gcf, 'reeds_range_profiles.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 7: Point cloud density comparison ──────────────────────────────
figure;
set(gcf, 'Position', [100 100 700 400]);
bar_data = [size(data_mar_raw,1), size(data_mar,1); ...
            size(data_apr_raw,1), size(data_apr,1)];
b3 = bar(bar_data);
b3(1).FaceColor = [0.4 0.4 0.4];
b3(2).FaceColor = [0.2 0.6 0.9];
set(gca, 'XTickLabel', {'March', 'April'});
ylabel('Number of points');
title('Point Cloud Density — Raw vs Filtered');
legend('Raw', 'Filtered');
grid on;
exportgraphics(gcf, 'reeds_point_density.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 8: Detection margin bar chart ──────────────────────────────────
figure;
set(gcf, 'Position', [100 100 700 400]);
margin_vals = [max(data_mar(mar_cr1_mask,5)) - FIXED_THRESH, ...
               max(data_mar(mar_cr2_mask,5)) - FIXED_THRESH; ...
               max(data_apr(apr_cr1_mask,5)) - FIXED_THRESH, ...
               max(data_apr(apr_cr2_mask,5)) - FIXED_THRESH];
b4 = bar(margin_vals);
b4(1).FaceColor = 'r';
b4(2).FaceColor = 'b';
hold on;
yline(0, 'k--', sprintf('Detection threshold (%.1f dB)', FIXED_THRESH), 'LineWidth', 2);
set(gca, 'XTickLabel', {'March', 'April'});
ylabel(sprintf('Margin above %.1f dB threshold (dB)', FIXED_THRESH));
title('CR Detection Margin — March vs April Reeds');
legend('CR1', 'CR2');
grid on;
exportgraphics(gcf, 'reeds_detection_margin.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 9: CR Peak Amplitude bar chart ─────────────────────────────────
figure;
set(gcf, 'Position', [100 100 700 400]);
cr_peaks = [max(data_mar(mar_cr1_mask,5)), max(data_mar(mar_cr2_mask,5)); ...
            max(data_apr(apr_cr1_mask,5)), max(data_apr(apr_cr2_mask,5))];
b5 = bar(cr_peaks);
b5(1).FaceColor = 'r';
b5(2).FaceColor = 'b';
hold on;
yline(FIXED_THRESH, 'k--', sprintf('Detection threshold (%.1f dB)', FIXED_THRESH), 'LineWidth', 2);
set(gca, 'XTickLabel', {'March', 'April'});
ylabel('Peak Amplitude (dB)');
title('Corner Reflector Peak Amplitude — March vs April Reeds');
legend('CR1', 'CR2');
grid on;
exportgraphics(gcf, 'reeds_cr_peak_amplitude.pdf', 'ContentType', 'image', 'Resolution', 300);

% ── FIGURE 10: Attenuation vs Range ───────────────────────────────────────
figure;
set(gcf, 'Position', [100 100 700 400]);
cr_ranges = [mean(r_mar(mar_cr1_mask)), mean(r_mar(mar_cr2_mask)), ...
             mean(r_apr(apr_cr1_mask)), mean(r_apr(apr_cr2_mask))];
cr_amps   = [max(data_mar(mar_cr1_mask,5)), max(data_mar(mar_cr2_mask,5)), ...
             max(data_apr(apr_cr1_mask,5)), max(data_apr(apr_cr2_mask,5))];
cr_labels = {'March CR1', 'March CR2', 'April CR1', 'April CR2'};
cr_colors = {'r', 'r', 'b', 'b'};
hold on;
for i = 1:4
    scatter(cr_ranges(i), cr_amps(i), 150, cr_colors{i}, 'filled');
    text(cr_ranges(i)+0.2, cr_amps(i), cr_labels{i}, 'FontSize', 9);
end
yline(FIXED_THRESH, 'k--', sprintf('Detection threshold (%.1f dB)', FIXED_THRESH), 'LineWidth', 2);
xlabel('Range (m)');
ylabel('Peak Amplitude (dB)');
title('CR Peak Amplitude vs Range — March vs April Reeds');
legend('March CRs', '', 'April CRs', 'Location', 'northwest');
grid on;
exportgraphics(gcf, 'reeds_attenuation_vs_range.pdf', 'ContentType', 'image', 'Resolution', 300);

fprintf('\nAll figures saved:\n');
fprintf('  1.  april_reeds_amplitude.pdf\n');
fprintf('  2.  march_reeds_amplitude.pdf\n');
fprintf('  3.  reeds_amp_vs_range.pdf\n');
fprintf('  4.  reeds_amplitude_histograms.pdf\n');
fprintf('  5.  reeds_snr_comparison.pdf\n');
fprintf('  6.  reeds_range_profiles.pdf\n');
fprintf('  7.  reeds_point_density.pdf\n');
fprintf('  8.  reeds_detection_margin.pdf\n');
fprintf('  9.  reeds_cr_peak_amplitude.pdf\n');
fprintf('  10. reeds_attenuation_vs_range.pdf\n');