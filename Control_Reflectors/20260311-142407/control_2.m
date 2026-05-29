clear; close all;

data = readmatrix('pointcloud_20260311-142407_.csv');

% Apply filters
data = data(data(:,5) >= -90, :);    % Amplitude filter
data = data(data(:,2) <= 50, :);    % Range filter

% Extract columns
% Col 1: timestamp, Col 2: range, Col 3: pitch, Col 4: roll, Col 5: amplitude
az = data(:, 4);
el = (pi/2) + data(:, 3);
r  = data(:, 2);

[x, y, z] = sph2cart(az, el, r);

% Figure 1: Coloured by amplitude
figure;
scatter3(x, y, z, 10, data(:, 5), 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Scan 20260311-115610 — Coloured by Amplitude');
axis equal; grid on; view(45, 30);

% Figure 2: Coloured by height (z)
figure;
scatter3(x, y, z, 10, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Scan 20260311-115610 — Coloured by Height');
axis equal; grid on; view(45, 30);

% Figure 3: Top-down view coloured by range
figure;
scatter(x, y, 10, r, 'filled');
colormap parula;
cb = colorbar;
cb.Label.String = 'Range (m)';
xlabel('X (m)'); ylabel('Y (m)');
title('Scan 20260311-115610 — Top-Down View');
axis equal; grid on;