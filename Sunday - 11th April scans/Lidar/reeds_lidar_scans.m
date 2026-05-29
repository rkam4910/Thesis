


data = readmatrix('small_corner_reflector_reeds.csv');

% Column indices (from CSV header):
%   9=X, 10=Y, 11=Z, 12=Reflectivity

x = data(:, 9);
y = data(:, 10);
z = data(:, 11);
refl = data(:, 12);

% Apply filters
dist = sqrt(x.^2 + y.^2 + z.^2);
keep = refl > 0 & refl <= 130 & dist > 0 & dist <= 15;
x    = x(keep);
y    = y(keep);
z    = z(keep);
refl = refl(keep);


% Figure 1: Coloured by reflectivity
figure;
scatter3(x, y, z, 2, refl, 'filled');
colormap jet;
cb = colorbar;
cb.Label.String = 'Reflectivity';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush\_Scan\_1 — Coloured by Reflectivity');
axis equal; grid on; view(45, 30);

% Figure 2: Coloured by height (Z)
figure;
scatter3(x, y, z, 2, z, 'filled');
colormap turbo;
cb = colorbar;
cb.Label.String = 'Height Z (m)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush\_Scan\_1 — Coloured by Height');
axis equal; grid on; view(45, 30);