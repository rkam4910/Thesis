clear; close all;

data = readmatrix('lidar_control.csv');

% Column indices (from CSV header):
%   9=X, 10=Y, 11=Z, 12=Reflectivity
x    = data(:, 9);
y    = data(:, 10);
z    = data(:, 11);
refl = data(:, 12);

% Apply filters
dist = sqrt(x.^2 + y.^2 + z.^2);
keep = refl > 0 & refl <= 130 & dist > 0 & dist <= 20;
x    = x(keep);
y    = y(keep);
z    = z(keep);
refl = refl(keep);

% Figure 1: Coloured by reflectivity
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, refl, 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Reflectivity';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Control Scan — Coloured by Reflectivity');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'lidar_control_reflectivity.pdf', 'ContentType', 'vector');

% Figure 2: Coloured by height (Z)
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Control Scan — Coloured by Height');
axis equal; grid on; view(45, 30);
exportgraphics(gcf, 'lidar_control_height.pdf', 'ContentType', 'vector');

% ── CORNER REFLECTOR IDENTIFICATION BY HEIGHT SPIKE ───────────────────────

cr_lidar = z > 0.3;   % points above canopy level

fprintf('\n── LiDAR Corner Reflector Locations ───────────\n');
fprintf('Total CR candidates (Z > 0.3 m): %d points\n', sum(cr_lidar));

x_cr = x(cr_lidar);
y_cr = y(cr_lidar);
z_cr = z(cr_lidar);

cr1 = y_cr < 0;
cr2 = y_cr >= 0;

fprintf('\nReflector 1:\n');
fprintf('  Points      : %d\n',     sum(cr1));
fprintf('  Position XYZ: (%.2f, %.2f, %.2f) m\n', ...
        mean(x_cr(cr1)), mean(y_cr(cr1)), mean(z_cr(cr1)));

fprintf('\nReflector 2:\n');
fprintf('  Points      : %d\n',     sum(cr2));
fprintf('  Position XYZ: (%.2f, %.2f, %.2f) m\n', ...
        mean(x_cr(cr2)), mean(y_cr(cr2)), mean(z_cr(cr2)));

d = sqrt((mean(x_cr(cr1)) - mean(x_cr(cr2)))^2 + ...
         (mean(y_cr(cr1)) - mean(y_cr(cr2)))^2 + ...
         (mean(z_cr(cr1)) - mean(z_cr(cr2)))^2);
fprintf('\nSeparation between reflectors: %.2f m\n', d);
fprintf('═══════════════════════════════════════════════\n\n');



% Figure 3: Corner reflectors highlighted in 3D
figure;
set(gcf, 'Position', [100 100 800 500]);
scatter3(x, y, z, 2, [0.7 0.7 0.7], 'filled'); hold on;
scatter3(x_cr(cr1), y_cr(cr1), z_cr(cr1), 150, 'r', 'filled');
scatter3(x_cr(cr2), y_cr(cr2), z_cr(cr2), 150, 'b', 'filled');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('LiDAR Control — Corner Reflectors Identified');
legend('Canopy / Ground', 'Reflector 1', 'Reflector 2');
axis equal; grid on; view(45, 30);
%exportgraphics(gcf, 'lidar_control_reflectors.pdf', 'ContentType', 'vector');

fprintf('\nReflector 1 position: (%.2f, %.2f, %.2f) m\n', ...
        mean(x_cr(cr1)), mean(y_cr(cr1)), mean(z_cr(cr1)));
fprintf('Reflector 2 position: (%.2f, %.2f, %.2f) m\n', ...
        mean(x_cr(cr2)), mean(y_cr(cr2)), mean(z_cr(cr2)));
d = sqrt((mean(x_cr(cr1))-mean(x_cr(cr2)))^2 + ...
         (mean(y_cr(cr1))-mean(y_cr(cr2)))^2 + ...
         (mean(z_cr(cr1))-mean(z_cr(cr2)))^2);
fprintf('Separation: %.2f m\n', d);

% Figure 4: Reflectivity vs distance
dist_f = sqrt(x.^2 + y.^2 + z.^2);
figure;
set(gcf, 'Position', [100 100 800 400]);
scatter(dist_f, refl, 1, 'b', 'filled');
xlabel('Distance (m)'); ylabel('Reflectivity');
title('Reflectivity vs Distance — LiDAR Control');
grid on;
%exportgraphics(gcf, 'lidar_control_refl_vs_dist.pdf', 'ContentType', 'vector');

