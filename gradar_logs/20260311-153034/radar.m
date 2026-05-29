clear; close all;
data = readmatrix('pointcloud_20260311-153034_.csv');

% Apply filters
data = data(data(:,5) >= -100, :);    % Amplitude filter
data = data(data(:,2) <= 30, :);     % Range filter

% Extract columns
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
title('Oldschool building — Coloured by Amplitude');
axis equal; grid on; view(45, 30);

% Figure 2: Coloured by height (z)
figure;
scatter3(x, y, z, 10, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Oldschool building — Coloured by Height');
axis equal; grid on; view(45, 30);