clear; close all;

data = readmatrix('Reeds_Scan.csv');
% 
% Column indices (from CSV header):
%   9=X, 10=Y, 11=Z, 12=Reflectivity
x    = data(:, 9);
y    = data(:, 10);
z    = data(:, 11);
refl = data(:, 12);

% ── FILTERS ───────────────────────────────────────────────────────────────

% Compute 3D Euclidean distance from sensor to each point
dist = sqrt(x.^2 + y.^2 + z.^2);

% Filter 1: Remove invalid returns (reflectivity = 0 means no return)
valid_refl = refl > 0;

% Filter 2: Remove saturated/specular returns (very high reflectivity)
%max_refl = refl <= 130;

% Filter 3: Remove points at zero distance (sensor self-returns)
valid_dist = dist > 0;

% Filter 4: Keep only points within x m (focus on near-field scene)
max_dist = dist <= 10;

% Combine all filters
keep = valid_refl & valid_dist & max_dist;

% Apply to all data columns
x    = x(keep);
y    = y(keep);
z    = z(keep);
refl = refl(keep);

fprintf('Points after filtering: %d\n', sum(keep));

% Figure 1: Coloured by reflectivity
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, refl, 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Reflectivity';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Reeds Lidar Scan — Coloured by Reflectivity');
axis equal; grid on; view(45, 45);
%exportgraphics(gcf, 'lidar_control_reflectivity.pdf', 'ContentType', 'vector');

% Figure 2: Coloured by height (Z)
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Reeds Lidar Scan — Coloured by Height');
axis equal; grid on; view(45, 30);
%exportgraphics(gcf, 'lidar_control_height.pdf', 'ContentType', 'vector');

dist = sqrt(x.^2 + y.^2 + z.^2);

figure;
set(gcf, 'Position', [100 100 700 400]);
histogram(dist, 50, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
xlabel('Distance from Sensor (m)');
ylabel('Number of Point Returns');
title('Point Returns vs Distance — Reeds Scan');
grid on;