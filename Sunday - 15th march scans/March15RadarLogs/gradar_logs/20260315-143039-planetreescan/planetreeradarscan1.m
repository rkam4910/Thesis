% Clear workspace and close all figures
clear; close all;

data = readmatrix('pointcloud_20260315-143039_.csv');  % Initial 1 degree scan

%Data Structure:   
% (1) Time Stamp (unix?)
% (2) Range (m)
% (3) Elevation (Rad)
% (4) Azimuth (Rad)
% (5) Return Strength  

data = data(data(:,5) >= -80, :); %Filtering radar returns below a power of -80 (not sure of the units) to get rid of side lobe returns in free space (Arbitrary value chosen may need to adjust)
data2 = data(data(:,2) <= 50, :); %Filtering radar returns below range of 50 meters

az = data(:, 4);         % Azimuth in radians
el = (pi/2)+data(:, 3);  % Elevation in radians  !!Need to subtract 90 degrees (pi/2 rad) from each elivation because it thinks 90 degrees is straight up!! (Given already negative just adding pi/2 to each value)
r = data(:, 2);          % Range in meters

returned_power = data(:, 5);

% Convert spherical to Cartesian coordinates
[x, y, z] = sph2cart(az, el, r);

%%  Plot with return power determining colour
figure;

scatter3(x, y, z, 10, returned_power, 'filled');  % Scatter plot with color based on return power

% Plotting the figures next to eachother instead of on-top of one another
pos1 = get(gcf,'Position'); % get position of Figure(1) 
set(gcf,'Position', pos1 - [pos1(3)/2,0,0,0]) % Shift position of Figure(1) 

ax = gca;
ax.Clipping = 'off'; % Prevents the plot from dissapearing when moving it around or zooming in 

colormap jet;
c = colorbar;
c.Label.String = "Returned Power";

xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('Plane tree orthogonal, 0.5° sample seperation');
axis equal;

%% Plot with z axis determining colour
figure;

scatter3(x, y, z, 10, z, 'filled');  % Scatter plot with color based on return power

% Plotting the figures next to eachother instead of on-top of one another
set(gcf,'Position', get(gcf,'Position') + [0,0,150,0]); % When Figure(2) is not the same size as Figure(1)
pos2 = get(gcf,'Position');  % get position of Figure(2) 
set(gcf,'Position', pos2 + [pos1(3)/2,0,0,0]) % Shift position of Figure(2)

ax = gca;
ax.Clipping = 'off'; % Prevents the plot from dissapearing when moving it around or zooming in 

clim([-3,20]) %Limiting colours to between -3 and 20 meters z height relative to radar position
colormap jet;
c = colorbar;
c.Label.String = "Z Height (m)";

xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('Plane tree orthogonal, 0.5° sample seperation');
axis equal;

%% Plotting figure of tree with filtered range 

az2 = data2(:, 4);         % Azimuth in radians
el2 = (pi/2)+data2(:, 3);  % Elevation in radians  !!Need to subtract 90 degrees (pi/2 rad) from each elivation because it thinks 90 degrees is straight up!! (Given already negative just adding pi/2 to each value)
r2 = data2(:, 2);          % Range in meters

returned_power = data2(:, 5);

% Convert spherical to Cartesian coordinates
[x, y, z] = sph2cart(az2, el2, r2);
figure;

scatter3(x, y, z, 10, returned_power, 'filled'); 

colormap jet;
c = colorbar;
c.Label.String = "Returned Power (m)";

xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('Plane tree orthogonal, 0.5° sample seperation, 50 meter range limit');
axis equal;

%% Plot with z axis determining colour
figure;

scatter3(x, y, z, 10, z, 'filled');  % Scatter plot with color based on return power

ax = gca;
ax.Clipping = 'off'; % Prevents the plot from dissapearing when moving it around or zooming in 

clim([-3,20]) %Limiting colours to between -3 and 20 meters z height relative to radar position
colormap jet;
c = colorbar;
c.Label.String = "Z Height (m)";

xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Z (meters)');
title('Plane tree orthogonal, 0.5° sample seperation, 50 meter range limit');
axis equal;

% Top-down plots
figure;
scatter(x, y, 8, returned_power, 'filled'); view(2);
axis equal; grid on; colormap turbo; colorbar;
xlabel("x' (penetration, m)"); ylabel("y' (through centre, m)");
title("Top-down tree centred");

figure;
scatter3(x, y, z, 8, returned_power, 'filled');
axis equal; grid on; colormap turbo; colorbar; view([45 30]);
xlabel("x' (m)"); ylabel("y' (m)"); zlabel('z (m)');
title("3D view Top down tree");

%% === 3D rotation video ===
% Assumes x, y, z, A_dB already exist
%% === 3D rotation video: force constant 1920x1080 frames ===
fig = figure('Color','w','Units','pixels','Position',[100 100 1920 1080]);
scatter3(x, y, z, 8, returned_power, 'filled');
axis equal vis3d; grid on; colormap turbo; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Plane tree point cloud (3D rotate)');

xl = xlim; yl = ylim; zl = zlim; cax = caxis;

outFile = fullfile(pwd, 'tree_rotate.mp4');
v = VideoWriter(outFile,'MPEG-4'); v.FrameRate = 24; open(v);

% grab a first frame and use its size as the target
drawnow;
f1 = getframe(fig);              % <-- capture from the figure
targetH = 1080; targetW = 1920;  % enforce 1080p
f1.cdata = imresize(f1.cdata, [targetH targetW]);
writeVideo(v, f1);

N = 240; elv = 30;
for a = linspace(0,360,N)
    xlim(xl); ylim(yl); zlim(zl); caxis(cax);
    view(a, elv);
    drawnow;

    f = getframe(fig);                                % may vary a few px on Windows
    if size(f.cdata,1) ~= targetH || size(f.cdata,2) ~= targetW
        f.cdata = imresize(f.cdata, [targetH targetW]); % force exact size
    end
    writeVideo(v, f);
end

close(v);
disp(['Saved: ' outFile]);


%top down plot filtered height

z_min = 10;
z_max = 15;
z_roi = (z >= z_min) & (z <= z_max);

x_roi = (3 <= x) & (x <= 30);
y_roi = (6  <= y) & (y <= 28);
sel = x_roi & y_roi & z_roi;

idx = find(sel);

%% 
% Top-down plots
figure;
scatter(x(sel), y(sel), 8, returned_power(sel), 'filled'); view(2);
axis equal; grid on; colormap turbo; colorbar;
xlabel("x' (penetration, m)"); ylabel("y' (through centre, m)");
title("Top-down tree centred (filtered height 10-15m)");

hold on

z_min = 2;
z_max = 10;
z_roi = (z >= z_min) & (z <= z_max);

x_roi = (3 <= x) & (x <= 30);
y_roi = (6  <= y) & (y <= 28);
sel = x_roi & y_roi & z_roi;

idx = find(sel);

figure;
scatter(x(sel), y(sel), 8, returned_power(sel), 'filled'); view(2);
axis equal; grid on; colormap turbo; colorbar;
xlabel("x' (penetration, m)"); ylabel("y' (through centre, m)");
title("Top-down tree centred (filtered height 2-10m)");





%% regions selected (height limited to range )
z_min = 1;
z_max = 3;
z_roi = (z >= z_min) & (z <= z_max);

x_roi = (3 <= x) & (x <= 30);
y_roi = (6  <= y) & (y <= 28);
sel = x_roi & y_roi & z_roi;

idx = find(sel);

%% --- CONFIG
z_lo = 10; z_hi = 15;                 % height slice (m)
x_lo = 3;  x_hi = 30;                 % lateral window (m)
y_lo = 6;  y_hi = 28;                 % along-tree window (m)
amp_thr = -80;                        % amplitude threshold (dB)
az_bin_deg = 1;                       % azimuth bin width (degrees)
p_bin_w = 0.25;                       % penetration bin width (m)

%% --- INPUTS assumed from your script
% Variables: x, y, z, r, returned_power, az (radians)

%% --- Slice + base filters
base = (z >= z_lo & z <= z_hi) & ...
       (x >= x_lo & x <= x_hi) & ...
       (y >= y_lo & y <= y_hi) & ...
       (returned_power >= amp_thr);

x1 = x(base); y1 = y(base); z1 = z(base);
r1 = r(base); a1 = returned_power(base); az1 = az(base);

% Range-correct power ~ R^-4 (two-way spreading) -> add 40*log10(r)
a_corr = a1 + 40*log10(r1 + eps);

%% --- Bin by azimuth, get front face per ray, compute penetration p
az_deg = rad2deg(az1);
edges_az = floor(min(az_deg)):az_bin_deg:ceil(max(az_deg));
[idx_az, edges_az] = discretize(az_deg, edges_az);

p_all = []; a_all = [];             % penetration and corrected power for all points

for k = 1:numel(edges_az)-1
    ik = (idx_az == k);
    if nnz(ik) < 30, continue; end  % skip sparse bins

    r_k = r1(ik);
    a_k = a_corr(ik);

    % Robust "front" range for this azimuth bin (5th percentile ≈ near edge)
    r_front = prctile(r_k, 5);

    % Penetration distance along the LOS for this ray
    p_k = r_k - r_front;

    % Keep only points at or behind the front (p >= 0)
    keep = (p_k >= 0);
    p_all = [p_all; p_k(keep)];
    a_all = [a_all; a_k(keep)];
end

%% --- Clean limits for plotting
% (limit penetration range if you want)
p_max = prctile(p_all, 98);
in = isfinite(p_all) & isfinite(a_all) & p_all <= p_max;
p_all = p_all(in); a_all = a_all(in);

%% --- Plot scatter (light) + binned medians + linear fit (slope dB/m)
figure('Color','w');
scatter(p_all, a_all, 8, 'filled', 'MarkerFaceAlpha', 0.10); hold on;

% Bin by penetration for a smooth trend
p_edges = 0:p_bin_w:(ceil(p_max/p_bin_w)*p_bin_w);
[p_idx, p_edges] = discretize(p_all, p_edges);
p_med = accumarray(p_idx(~isnan(p_idx)), p_all(~isnan(p_idx)), [], @median, NaN);
a_med = accumarray(p_idx(~isnan(p_idx)), a_all(~isnan(p_idx)), [], @median, NaN);
p_ctr = p_edges(1:end-1) + diff(p_edges)/2;
valid = isfinite(p_med) & isfinite(a_med);

plot(p_ctr(valid), a_med(valid), 'LineWidth', 2);

% Linear fit on medians to estimate attenuation slope (dB/m)
pf = polyfit(p_ctr(valid), a_med(valid), 1);   % a ≈ pf(1)*p + pf(2)
slope_db_per_m = pf(1);
plot(p_ctr(valid), polyval(pf, p_ctr(valid)), '--', 'LineWidth', 1.5);

grid on; xlabel('Penetration distance p (m)');
ylabel('Range-corrected power A_{corr} (dB)');
title(sprintf('Attenuation in canopy (z = %g–%g m), slope = %.2f dB/m', z_lo, z_hi, slope_db_per_m));


%% --- CONFIG
z_lo = 10; z_hi = 15;                 % height slice (m)
x_lo = 3;  x_hi = 30;                 % lateral window (m)
y_lo = 6;  y_hi = 28;                 % along-tree window (m)
amp_thr = -80;                        % amplitude threshold (dB)
az_bin_deg = 0.5;                       % azimuth bin width (degrees)
p_bin_w = 0.25;                       % penetration bin width (m)

%% --- INPUTS assumed from your script
% Variables: x, y, z, r, returned_power, az (radians)

%% --- Slice + base filters
base = (z >= z_lo & z <= z_hi) & ...
       (x >= x_lo & x <= x_hi) & ...
       (y >= y_lo & y <= y_hi) & ...
       (returned_power >= amp_thr);

x1 = x(base); y1 = y(base); z1 = z(base);
r1 = r(base); a1 = returned_power(base); az1 = az(base);

% Range-correct power ~ R^-4 (two-way spreading) -> add 40*log10(r)
a_corr = a1 + 40*log10(r1 + eps);

%% --- Bin by azimuth, get front face per ray, compute penetration p
az_deg = rad2deg(az1);
edges_az = floor(min(az_deg)):az_bin_deg:ceil(max(az_deg));
[idx_az, edges_az] = discretize(az_deg, edges_az);

p_all = []; a_all = [];             % penetration and corrected power for all points

for k = 1:numel(edges_az)-1
    ik = (idx_az == k);
    if nnz(ik) < 30, continue; end  % skip sparse bins

    r_k = r1(ik);
    a_k = a_corr(ik);

    % Robust "front" range for this azimuth bin (5th percentile ≈ near edge)
    r_front = prctile(r_k, 5);

    % Penetration distance along the LOS for this ray
    p_k = r_k - r_front;

    % Keep only points at or behind the front (p >= 0)
    keep = (p_k >= 0);
    p_all = [p_all; p_k(keep)];
    a_all = [a_all; a_k(keep)];
end

%% --- Clean limits for plotting
% (limit penetration range if you want)
p_max = prctile(p_all, 98);
in = isfinite(p_all) & isfinite(a_all) & p_all <= p_max;
p_all = p_all(in); a_all = a_all(in);

%%
%% === Attenuation histogram: 10–15 m height ===
z_min = 10; 
z_max = 15;
z_roi = (z >= z_min) & (z <= z_max);

x_roi = (3 <= x) & (x <= 30);
y_roi = (6  <= y) & (y <= 28);
sel = z_roi & x_roi & y_roi;

% Select only those points
x_sel = x(sel);
y_sel = y(sel);
r_sel = r(sel);
p_sel = returned_power(sel);

% Sort points by range (penetration into tree)
[~, sortIdx] = sort(r_sel);
r_sorted = r_sel(sortIdx);
p_sorted = p_sel(sortIdx);

% Bin ranges every 0.5 m (adjust as needed)
binEdges = min(r_sorted):0.5:max(r_sorted);
[binCounts,~,binIdx] = histcounts(r_sorted, binEdges);

% Compute mean or median power per bin
binPower = accumarray(binIdx(~isnan(binIdx)), p_sorted(~isnan(binIdx)), [], @median, NaN);

binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% Normalize colors between strongest and weakest return
normPow = (binPower - min(binPower)) ./ (max(binPower) - min(binPower));

% Create histogram bars manually with color gradient
figure('Color','w');
hold on;
colormap turbo;
for i = 1:length(binCenters)
    if isnan(binPower(i)), continue; end
    color = interp1(linspace(0,1,256), colormap, normPow(i));
    bar(binCenters(i), binCounts(i), 1.0, 'FaceColor', color, 'EdgeColor', 'none');
end
hold off;

xlabel('Range from radar (m)');
ylabel('Number of returns');
title('Attenuation through tree canopy (10–15 m slice)');
colorbar('Ticks',[0 1],'TickLabels',{'Weak','Strong'},...
         'LabelString','Median returned power');
grid on;


%% === Attenuation curtain (z = 10–15 m) ===
% Assumes you already have: x, y, z, r, returned_power (all column vectors)

% -------- parameters (tweak as needed) --------
z_lo = 10; z_hi = 15;                  % height slice (m)
x_lo = 3;  x_hi = 30;                  % lateral x' window (m)
y_lo = -5; y_hi = 10;                  % transverse y' window (m)  (matches screenshot axes)
amp_thr = -80;                         % keep only reasonably strong points (dB)
front_pct = 5;                         % "front" = 5th percentile of x' in each y-strip

dy = 0.5;                              % y' bin width (m)
dd = 0.2;                              % depth bin width (m)

% -------- base mask and range correction --------
base = (z >= z_lo & z <= z_hi) & ...
       (x >= x_lo & x <= x_hi) & ...
       (y >= y_lo & y <= y_hi) & ...
       isfinite(r) & isfinite(returned_power) & ...
       (returned_power >= amp_thr);

x1 = x(base); y1 = y(base); r1 = r(base); a1 = returned_power(base);
% Two-way geometric spreading ~ r^-4 → add 40 log10 r (in dB)
a_corr = a1 + 40*log10(r1 + eps);

% -------- bin definitions --------
edges_y = y_lo:dy:y_hi;                         % transverse bins
cy = edges_y(1:end-1) + dy/2;                   % y' bin centers

% For depth (d): derive edges from data after alignment; but we need a guess first
% Use x-range to set a safe max depth; we’ll clamp later.
d_max_guess = (x_hi - x_lo);
edges_d = 0:dd:d_max_guess;
cd = edges_d(1:end-1) + dd/2;                   % depth bin centers

% -------- build curtain: rows = y' bins, cols = depth bins --------
ny = numel(edges_y)-1; nd = numel(edges_d)-1;
Curtain = nan(ny, nd);                           % median Acorr per (y, depth)

for iy = 1:ny
    in_y = (y1 >= edges_y(iy) & y1 < edges_y(iy+1));
    if ~any(in_y), continue; end

    xk = x1(in_y);
    ak = a_corr(in_y);

    % define the "front" (nearest x' within this y strip)
    x_front = prctile(xk, front_pct);           % robust against outliers

    % penetration depth along x' (clamp to >=0)
    dk = max(0, xk - x_front);

    % bin by depth
    [~, ~, idx_d] = histcounts(dk, edges_d);
    valid = idx_d > 0 & idx_d <= nd;

    if any(valid)
        % median per depth bin for this y row
        med_row = accumarray(idx_d(valid), ak(valid), [nd,1], @median, nan);
        Curtain(iy, :) = med_row.';             % row vector
    end
end

% Optional: trim empty leading/trailing depth columns
nonempty_cols = any(isfinite(Curtain), 1);
if any(nonempty_cols)
    Curtain = Curtain(:, nonempty_cols);
    cd = cd(nonempty_cols);
end

% -------- plot (like your screenshot) --------
figure('Color','w');
imagesc(cd, cy, Curtain); axis xy;
colormap turbo; cb = colorbar; grid on;
cb.Label.String = 'Median amplitude (dB, range-corrected)';
xlabel('Depth d along x'' (m)');
ylabel('y'' (m)');
title('Depth-aligned median amplitude (dB) — attenuation curtain');

% nice axes (optional)
xlim([min(cd) max(cd)]);
ylim([y_lo y_hi]);

% Optional: set a robust color scale (ignore extremes)
vals = Curtain(isfinite(Curtain));
if ~isempty(vals)
    caxis(prctile(vals, [2 98]));
end

