% Clear workspace and close all figures
clear; close all;

data = readmatrix('pointcloud_20000101-140013_.csv');  % Initial 1 degree scan

%Data Structure:   
% (1) Time Stamp (unix?)
% (2) Range (m)
% (3) Elevation (Rad)
% (4) Azimuth (Rad)
% (5) Return Strength  

data = data(data(:,5) >= -80, :); %Filtering radar returns below a power of -80 (not sure of the units) to get rid of side lobe returns in free space (Arbitrary value chosen may need to adjust)
data2 = data(data(:,2) <= 60, :); %Filtering radar returns below range of 50 meters

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
