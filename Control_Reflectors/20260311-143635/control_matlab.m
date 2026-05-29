clear; close all;

data = readmatrix('pointcloud_20260311-143635_.csv');

% Apply filters
data = data(data(:,5) >= -80, :);    % Amplitude filter
data = data(data(:,2) <= 50, :);     % Range filter

% Extract columns
% Col 1: timestamp, Col 2: range, Col 3: pitch, Col 4: roll, Col 5: amplitude
az = data(:, 4);
el = (pi/2) + data(:, 3);
r  = data(:, 2);

[x, y, z] = sph2cart(az, el, r);

% Figure 1: Coloured by amplitude
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 10, data(:,5), 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Radar Control Scan — Coloured by Amplitude');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'radar_control_amplitude.pdf', 'ContentType', 'vector');

% Figure 2: Coloured by height (Z)
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 10, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Radar Control Scan— Coloured by Height');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'radar_control_height.pdf', 'ContentType', 'vector');

% Figure 3: Top-down view coloured by range
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter(x, y, 10, r, 'filled');
colormap parula;
cb = colorbar;
cb.Label.String = 'Range (m)';
xlabel('X (m)'); ylabel('Y (m)');
title('Radar Control Scan— Top-Down View');
axis equal; grid on;
exportgraphics(gcf, 'radar_control_topdown.pdf', 'ContentType', 'vector');

% ── CORNER REFLECTOR ANALYSIS ─────────────────────────────────────────────
% CRs identified by range (frame alignment with LiDAR not assumed)
% CR1 at ~19.8 m, CR2 at ~23.8 m (confirmed from amplitude vs range plot)
% LiDAR ground truth separation: 4.82 m

cr_thresh = 20;                                % dB threshold
cr_mask   = data(:,5) >= cr_thresh & r <= 30;  % strong returns within 30 m

x_cr   = x(cr_mask);
y_cr   = y(cr_mask);
z_cr   = z(cr_mask);
r_cr   = r(cr_mask);
amp_cr = data(cr_mask, 5);
az_cr = az(cr_mask);

% Split by BOTH range AND azimuth (roll angle)
cr1_mask = r_cr <= 22  & az_cr >= 1.20 & az_cr <= 1.30;
cr2_mask = r_cr >  22  & r_cr <= 30 & az_cr >= 0.42 & az_cr <= 0.52;

% Compute radar separation
cr1_pos = [mean(x_cr(cr1_mask)), mean(y_cr(cr1_mask)), mean(z_cr(cr1_mask))];
cr2_pos = [mean(x_cr(cr2_mask)), mean(y_cr(cr2_mask)), mean(z_cr(cr2_mask))];
sep_radar = norm(cr1_pos - cr2_pos);
sep_lidar = 4.82;  % LiDAR ground truth

fprintf('\n── Corner Reflector Analysis ──────────────────\n');
fprintf('Total strong returns (>= %d dB, <= 30 m): %d points\n', cr_thresh, sum(cr_mask));

fprintf('\nReflector 1 (~19.8 m range):\n');
fprintf('  Points   : %d\n',        sum(cr1_mask));
fprintf('  Range    : %.2f m\n',    mean(r_cr(cr1_mask)));
fprintf('  Peak amp : %.1f dB\n',   max(amp_cr(cr1_mask)));
fprintf('  Mean amp : %.1f dB\n',   mean(amp_cr(cr1_mask)));

fprintf('\nReflector 2 (~23.8 m range):\n');
fprintf('  Points   : %d\n',        sum(cr2_mask));
fprintf('  Range    : %.2f m\n',    mean(r_cr(cr2_mask)));
fprintf('  Peak amp : %.1f dB\n',   max(amp_cr(cr2_mask)));
fprintf('  Mean amp : %.1f dB\n',   mean(amp_cr(cr2_mask)));

fprintf('\nRange difference between CRs:\n');
fprintf('  Radar  : %.2f m\n',      mean(r_cr(cr2_mask)) - mean(r_cr(cr1_mask)));
fprintf('  LiDAR  : %.2f m\n',      sep_lidar);
fprintf('═══════════════════════════════════════════════\n\n');

% Figure 4: Amplitude vs range
figure;
set(gcf, 'Position', [100 100 800 400]);
scatter(r, data(:,5), 2, 'b', 'filled'); hold on;
scatter(r_cr(cr1_mask), amp_cr(cr1_mask), 60, 'r', 'filled');
scatter(r_cr(cr2_mask), amp_cr(cr2_mask), 60, 'g', 'filled');
xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('Radar Control — Amplitude vs Range');
legend('All returns', 'Reflector 1 (~19.8 m)', 'Reflector 2 (~23.8 m)');
grid on;
exportgraphics(gcf, 'radar_control_amp_vs_range.pdf', 'ContentType', 'vector');

% Figure 5: Corner reflectors highlighted in 3D
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, [0.5 0.5 0.5], 'filled'); hold on;
scatter3(x_cr(cr1_mask), y_cr(cr1_mask), z_cr(cr1_mask), 100, 'r', 'filled');
scatter3(x_cr(cr2_mask), y_cr(cr2_mask), z_cr(cr2_mask), 100, 'b', 'filled');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Radar Control — Corner Reflectors Identified');
legend('Clutter', 'Reflector 1 (~19.8 m)', 'Reflector 2 (~23.8 m)');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'radar_control_reflectors.pdf', 'ContentType', 'vector');

% Draw boxes around corner reflectors in Figure 1
figure(1);
hold on;

% Box around CR1 (~19.8 m)
x1 = [4 8 8 4 4];
y1 = [17 17 20 20 17];
z1_bot = -2; z1_top = 0;
plot3(x1, y1, repmat(z1_bot,1,5), 'k-', 'LineWidth', 2);
plot3(x1, y1, repmat(z1_top,1,5), 'k-', 'LineWidth', 2);
for i = 1:4
    plot3([x1(i) x1(i)], [y1(i) y1(i)], [z1_bot z1_top], 'k-', 'LineWidth', 2);
end
text(6, 18.5, 1, 'CR1', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');

% Box around CR2 (~23.8 m)
x2 = [19 23 23 19 19];
y2 = [9 9 12 12 9];
z2_bot = -2; z2_top = 0;
plot3(x2, y2, repmat(z2_bot,1,5), 'k-', 'LineWidth', 2);
plot3(x2, y2, repmat(z2_top,1,5), 'k-', 'LineWidth', 2);
for i = 1:4
    plot3([x2(i) x2(i)], [y2(i) y2(i)], [z2_bot z2_top], 'k-', 'LineWidth', 2);
end
text(21, 10.5, 1, 'CR2', 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold');

exportgraphics(gcf, 'radar_control_amplitude_annotated.pdf', 'ContentType', 'image', 'Resolution', 300);