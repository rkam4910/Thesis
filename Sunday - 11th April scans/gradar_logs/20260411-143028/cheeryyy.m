% April Bush Scan — 11 Apr 2026, 13:09
clear; close all;

data = readmatrix('pointcloud_20260411-143028_.csv');

% Apply filters
data = data(data(:,5) >= -80, :);    % Amplitude filter
data = data(data(:,2) <= 6, :);     % Range filter

% Extract columns
az = data(:, 4);
el = (pi/2) + data(:, 3);
r  = data(:, 2);

[x, y, z] = sph2cart(az, el, r);

% Figure 1: Coloured by amplitude
figure;
scatter3(x, y, z, 5, data(:,5), 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Charry Ballart Shrub — Coloured by Amplitude');
axis equal; grid on; view(45, 30);

%11 Apr 2026, 13:09 

% Figure 2: Coloured by height (Z)
figure;
scatter3(x, y, z, 5, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('pointcloud\_20260411-130952\_ — Coloured by Height');
axis equal; grid on; view(45, 30);