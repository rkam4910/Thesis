clear; close all;

% ==========================================================
% CONFIG
% ==========================================================
BUSH_FRONT  = 1.5;
BUSH_BACK   = 3.4;
CR_RANGE    = 3.72;
CR_TOL      = 0.3;
BIN_WIDTH   = 0.25;
MIN_PTS_BIN = 5;

% ==========================================================
% LOAD SCAN 1
% ==========================================================
data1 = readmatrix('Bush_Scan_1.csv');
x1    = data1(:, 9);
y1    = data1(:, 10);
z1    = data1(:, 11);
refl1 = data1(:, 12);
dist1 = sqrt(x1.^2 + y1.^2 + z1.^2);
keep1 = refl1 > 0 & refl1 <= 130 & dist1 > 0 & dist1 <= 6;
x1    = x1(keep1);
y1    = y1(keep1);
z1    = z1(keep1);
refl1 = refl1(keep1);
dist1 = dist1(keep1);

% ==========================================================
% LOAD SCAN 2
% ==========================================================
data2 = readmatrix('Bush_Scan_2.csv');
x2    = data2(:, 9);
y2    = data2(:, 10);
z2    = data2(:, 11);
refl2 = data2(:, 12);
dist2 = sqrt(x2.^2 + y2.^2 + z2.^2);
keep2 = refl2 > 0 & refl2 <= 130 & dist2 > 0 & dist2 <= 6;
x2    = x2(keep2);
y2    = y2(keep2);
z2    = z2(keep2);
refl2 = refl2(keep2);
dist2 = dist2(keep2);

fprintf('Scan 1: %d points\n', sum(keep1));
fprintf('Scan 2: %d points\n', sum(keep2));

% ==========================================================
% HELPER — range binned stats
% ==========================================================
function tbl = bin_stats(dist, refl, r_lo, r_hi, bw, min_pts)
    edges = r_lo:bw:r_hi;
    tbl   = [];
    for i = 1:length(edges)-1
        sel = dist >= edges(i) & dist < edges(i) + bw;
        if sum(sel) >= min_pts
            tbl(end+1,:) = [edges(i)+bw/2, sum(sel), ...
                            mean(refl(sel)), std(refl(sel))];
        end
    end
end

% ==========================================================
% FIGURE 1 — Scan 1 coloured by reflectivity
% ==========================================================
figure;
scatter3(x1, y1, z1, 2, refl1, 'filled');
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush Scan 1 — Coloured by Reflectivity');
axis equal; grid on; view(45, 30);

% ==========================================================
% FIGURE 2 — Scan 1 coloured by height
% ==========================================================
figure;
scatter3(x1, y1, z1, 2, z1, 'filled');
colormap turbo; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush Scan 1 — Coloured by Height');
axis equal; grid on; view(45, 30);

% ==========================================================
% FIGURE 3 — Scan 2 coloured by reflectivity
% ==========================================================
figure;
scatter3(x2, y2, z2, 2, refl2, 'filled');
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush Scan 2 — Coloured by Reflectivity');
axis equal; grid on; view(45, 30);

% ==========================================================
% FIGURE 4 — Scan 2 coloured by height
% ==========================================================
figure;
scatter3(x2, y2, z2, 2, z2, 'filled');
colormap turbo; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Bush Scan 2 — Coloured by Height');
axis equal; grid on; view(45, 30);

% ==========================================================
% FIGURE 5 — Normalised density vs range (both scans)
% ==========================================================
edges = 0:BIN_WIDTH:6;
rmids = edges(1:end-1) + BIN_WIDTH/2;

counts1 = histcounts(dist1, edges);
counts2 = histcounts(dist2, edges);
counts1_norm = counts1 / max(counts1);
counts2_norm = counts2 / max(counts2);

figure; hold on;
plot(rmids, counts1_norm, 'b-', 'LineWidth', 2, 'DisplayName', 'Scan 1');
plot(rmids, counts2_norm, 'r-', 'LineWidth', 2, 'DisplayName', 'Scan 2');
xline(BUSH_FRONT, 'g--', sprintf('Bush front (%.1f m)', BUSH_FRONT), ...
      'LineWidth', 1.5);
xline(BUSH_BACK,  'g-.', sprintf('Bush back (%.1f m)', BUSH_BACK), ...
      'LineWidth', 1.5);
xline(CR_RANGE,   'r-',  sprintf('CR (%.2f m)', CR_RANGE), ...
      'LineWidth', 1.5);
xlabel('Range (m)', 'FontSize', 11);
ylabel('Normalised point density', 'FontSize', 11);
title('LiDAR Point Density vs Range — Cherry Ballart Bush Scans', ...
      'FontSize', 11, 'FontWeight', 'bold');
legend('FontSize', 9); grid on;

% ==========================================================
% FIGURE 6 — Mean reflectivity vs range (both scans)
% ==========================================================
prof1 = bin_stats(dist1, refl1, 0, 6, BIN_WIDTH, MIN_PTS_BIN);
prof2 = bin_stats(dist2, refl2, 0, 6, BIN_WIDTH, MIN_PTS_BIN);

figure; hold on;

if ~isempty(prof1)
    errorbar(prof1(:,1), prof1(:,3), prof1(:,4), ...
        'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', 'Scan 1 mean ±1σ');
end
if ~isempty(prof2)
    errorbar(prof2(:,1), prof2(:,3), prof2(:,4), ...
        'r-s', 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', 'Scan 2 mean ±1σ');
end

xline(BUSH_FRONT, 'g--', sprintf('Bush front (%.1f m)', BUSH_FRONT), ...
      'LineWidth', 1.5);
xline(BUSH_BACK,  'g-.', sprintf('Bush back (%.1f m)', BUSH_BACK), ...
      'LineWidth', 1.5);
xline(CR_RANGE,   'r-',  sprintf('CR (%.2f m)', CR_RANGE), ...
      'LineWidth', 1.5);

xlabel('Range (m)', 'FontSize', 11);
ylabel('Mean reflectivity', 'FontSize', 11);
title('LiDAR Mean Reflectivity vs Range — Cherry Ballart Bush Scans', ...
      'FontSize', 11, 'FontWeight', 'bold');
legend('FontSize', 9); grid on;

% ==========================================================
% PRINT SUMMARY
% ==========================================================
fprintf('\n========================================\n');
fprintf('LIDAR SUMMARY\n');
fprintf('========================================\n');
for s = 1:2
    if s == 1
        d = dist1; r = refl1; label = 'Scan 1';
    else
        d = dist2; r = refl2; label = 'Scan 2';
    end
    cr_sel = d >= CR_RANGE - CR_TOL & d <= CR_RANGE + CR_TOL;
    bush_sel = d >= BUSH_FRONT & d <= BUSH_BACK;
    fprintf('\n%s:\n', label);
    fprintf('  Total points:              %d\n', length(d));
    fprintf('  Points in CR window:       %d\n', sum(cr_sel));
    if sum(cr_sel) > 0
        fprintf('  CR window mean refl:       %.1f\n', mean(r(cr_sel)));
        fprintf('  CR window max refl:        %.1f\n', max(r(cr_sel)));
    else
        fprintf('  No points in CR window\n');
    end
    fprintf('  Bush body mean refl:       %.1f\n', mean(r(bush_sel)));
end
fprintf('========================================\n');