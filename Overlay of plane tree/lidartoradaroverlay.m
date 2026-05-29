% =========================================================================
% overlay_radar_lidar.m
% Overlays aligned radar and LiDAR point clouds in a common GPS world frame.
% Runs both alignment pipelines then plots everything together.
%
% Radar:  2 scans from overlayofplanetree.m pipeline
% LiDAR:  4 scans from align_lidar_scans.m pipeline
% Both share the same ENU origin: Sunday radar/lidar GPS position
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
% Radar files
file_scan1 = 'pointcloud_20260315-143039_.csv';
file_scan2 = 'pointcloud_20000101-140013_.csv';

% LiDAR files
file_sun_big   = '2026-03-15 14-50-46_corner_big.csv';
file_sun_small = '2026-03-15 14-52-14_corner_small.csv';
file_mon_big   = '2026-03-16_bigreflector_planetree.csv';
file_mon_small = '2026-03-16_smallreflector_planetree.csv';

% GPS files (shared)
file_gps1 = 'gps_tree scan 1.csv';
file_gps2 = 'tree scan 2.csv';

% Radar settings
AMP_THRESH_SCAN1 = 50;
AMP_THRESH_SCAN2 = 20;
RADAR_RANGE_LIMIT = 40;
ROT_CORRECTION_DEG = 0;
z_offset_empirical = 4.3;

% LiDAR settings
RANGE_LIMIT  = 60;
FILTER_ZEROS = true;
SUBSAMPLE_N  = 5;

% Plot settings
RADAR_POINT_SIZE = 6;
LIDAR_POINT_SIZE = 2;

% Noise filters
K_KNN          = 20;
STD_MULTIPLIER = 1.5;
Z_MAX          = 25;
EX1_N=12; EX1_E=-6;  EX1_R=10;
EX2_N=44; EX2_E=-24; EX2_R=8;

% =========================================================================
fprintf('=== SECTION 1: LOAD AND ALIGN GPS ===\n');
% =========================================================================
opts1 = detectImportOptions(file_gps1,'Delimiter',',','VariableNamingRule','preserve');
G1    = readtable(file_gps1, opts1);
opts2 = detectImportOptions(file_gps2,'Delimiter',',','VariableNamingRule','preserve');
G2    = readtable(file_gps2, opts2);
G1.Properties.VariableNames{1} = 'PointName';
G2.Properties.VariableNames{1} = 'PointName';
G1.PointName = strtrim(G1.PointName);
G2.PointName = strtrim(G2.PointName);

% World origin = Sunday radar/lidar GPS
idx  = strcmp(G1.PointName,'radar/lidar');
olon = G1.Longitude(idx); olat = G1.Latitude(idx);
oh   = G1.("Ellipsoidal height")(idx);

gps2enu = @(lon,lat,h) deal(...
    (lon-olon).*111320.*cos(deg2rad(olat)), ...
    (lat-olat).*110540, h-oh);

[e,n,u]=gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector big')),G1.Latitude(strcmp(G1.PointName,'corner reflector big')),G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector big')));
wf1_big=[e,n,u];
[e,n,u]=gps2enu(G1.Longitude(strcmp(G1.PointName,'corner reflector small')),G1.Latitude(strcmp(G1.PointName,'corner reflector small')),G1.("Ellipsoidal height")(strcmp(G1.PointName,'corner reflector small')));
wf1_small=[e,n,u];
[e,n,u]=gps2enu(G2.Longitude(strcmp(G2.PointName,'corner refelctor big1')),G2.Latitude(strcmp(G2.PointName,'corner refelctor big1')),G2.("Ellipsoidal height")(strcmp(G2.PointName,'corner refelctor big1')));
wf2_big=[e,n,u];
[e,n,u]=gps2enu(G2.Longitude(strcmp(G2.PointName,'corner refelctor small')),G2.Latitude(strcmp(G2.PointName,'corner refelctor small')),G2.("Ellipsoidal height")(strcmp(G2.PointName,'corner refelctor small')));
wf2_small=[e,n,u];
[e,n,u]=gps2enu(G2.Longitude(strcmp(G2.PointName,'Radar/Lidar')),G2.Latitude(strcmp(G2.PointName,'Radar/Lidar')),G2.("Ellipsoidal height")(strcmp(G2.PointName,'Radar/Lidar')));
wf2_radar=[e,n,u];
[e,n,u]=gps2enu(G2.Longitude(strcmp(G2.PointName,'Radar/Lidar')),G2.Latitude(strcmp(G2.PointName,'Radar/Lidar')),G2.("Ellipsoidal height")(strcmp(G2.PointName,'Radar/Lidar')));
wf2_lidar=[e,n,u];

fprintf('  GPS loaded. World origin: lon=%.8f lat=%.8f h=%.3fm\n', olon, olat, oh);

% =========================================================================
fprintf('\n=== SECTION 2: LOAD AND ALIGN RADAR ===\n');
% =========================================================================
T1 = readtable(file_scan1);
T2 = readtable(file_scan2);

X1_all=T1.range.*sin(-T1.pitch).*cos(T1.roll); Y1_all=T1.range.*sin(-T1.pitch).*sin(T1.roll); Z1_all=-T1.range.*cos(-T1.pitch);
X2_all=T2.range.*sin(-T2.pitch).*cos(T2.roll); Y2_all=T2.range.*sin(-T2.pitch).*sin(T2.roll); Z2_all=-T2.range.*cos(-T2.pitch);
m1=T1.range<=RADAR_RANGE_LIMIT; m2=T2.range<=RADAR_RANGE_LIMIT;
X1=X1_all(m1);Y1=Y1_all(m1);Z1=Z1_all(m1);amp1=T1.amplitude(m1);
X2=X2_all(m2);Y2=Y2_all(m2);Z2=Z2_all(m2);amp2=T2.amplitude(m2);

% Reflector centroids
function [c_near,c_far]=radar_centroids(Xs,Ys,Zs,amps,Rs,thresh)
    mask=(amps>thresh)&(Zs>-20);
    pts_r=Rs(mask); pts_x=Xs(mask); pts_y=Ys(mask); pts_z=Zs(mask);
    [pts_r_sorted,sort_idx]=sort(pts_r);
    [~,gap_idx]=max(diff(pts_r_sorted));
    ni=sort_idx(1:gap_idx); fi=sort_idx(gap_idx+1:end);
    c_near=[mean(pts_x(ni)),mean(pts_y(ni)),mean(pts_z(ni))];
    c_far =[mean(pts_x(fi)),mean(pts_y(fi)),mean(pts_z(fi))];
end

[c1_near,c1_far]=radar_centroids(X1,Y1,Z1,amp1,T1.range(m1),AMP_THRESH_SCAN1);
[c2_near,c2_far]=radar_centroids(X2,Y2,Z2,amp2,T2.range(m2),AMP_THRESH_SCAN2);

% Match centroids to GPS
d1_big=norm(wf1_big(1:2)); d1_small=norm(wf1_small(1:2));
if abs(norm(c1_near)-d1_big)<abs(norm(c1_near)-d1_small); s1_near_wf=wf1_big; s1_far_wf=wf1_small;
else; s1_near_wf=wf1_small; s1_far_wf=wf1_big; end
d2_big=norm(wf2_big(1:2)-wf2_radar(1:2)); d2_small=norm(wf2_small(1:2)-wf2_radar(1:2));
if abs(norm(c2_near)-d2_big)<abs(norm(c2_near)-d2_small); s2_near_wf=wf2_big; s2_far_wf=wf2_small;
else; s2_near_wf=wf2_small; s2_far_wf=wf2_big; end

sensor_pts_1=[c1_near;c1_far;0,0,0]; world_pts_1=[s1_near_wf;s1_far_wf;0,0,0];
sensor_pts_2=[c2_near;c2_far;0,0,0]; world_pts_2=[s2_near_wf;s2_far_wf;wf2_radar];

function [theta,tx,ty]=solve_2d(sp,wp)
    mu_s=mean(sp,1); mu_w=mean(wp,1);
    S=(sp-mu_s)'*(wp-mu_w); [U,~,V]=svd(S);
    R2d=V*diag([1,det(V*U')])*U';
    theta=atan2(R2d(2,1),R2d(1,1));
    t=mu_w-mu_s*R2d'; tx=t(1); ty=t(2);
end

[theta1,tx1,ty1]=solve_2d(sensor_pts_1(:,1:2),world_pts_1(:,1:2));
[theta2,tx2,ty2]=solve_2d(sensor_pts_2(:,1:2),world_pts_2(:,1:2));

c1=cos(theta1);s1r=sin(theta1); c2=cos(theta2);s2r=sin(theta2);
pts1_xy=[X1,Y1]*[c1 s1r;-s1r c1]+[tx1,ty1];
pts2_xy=[X2,Y2]*[c2 s2r;-s2r c2]+[tx2,ty2];
z_offset=wf2_radar(3);
pts1_world=[pts1_xy,-Z1];
pts2_world=[pts2_xy,-Z2+z_offset-z_offset_empirical];

fprintf('  Radar Scan1: yaw=%.2f deg  Scan2: yaw=%.2f deg\n',rad2deg(theta1),rad2deg(theta2));

% =========================================================================
fprintf('\n=== SECTION 3: LOAD AND ALIGN LIDAR ===\n');
% =========================================================================
function [x,y,z,refl]=load_lidar(fname,filter_zeros,range_limit,subsample_n)
    T=readtable(fname,'VariableNamingRule','preserve');
    get_col=@(name) T{:,find(strcmpi(T.Properties.VariableNames,name),1)};
    x=get_col('X'); y=get_col('Y'); z=get_col('Z'); refl=get_col('Reflectivity');
    if filter_zeros; v=~(x==0&y==0&z==0); x=x(v);y=y(v);z=z(v);refl=refl(v); end
    rng=sqrt(x.^2+y.^2+z.^2); keep=rng<=range_limit;
    x=x(keep);y=y(keep);z=z(keep);refl=refl(keep);
    idx=1:subsample_n:numel(x); x=x(idx);y=y(idx);z=z(idx);refl=refl(idx);
end

[X_sb,Y_sb,Z_sb,R_sb]=load_lidar(file_sun_big,  FILTER_ZEROS,RANGE_LIMIT,SUBSAMPLE_N);
[X_ss,Y_ss,Z_ss,R_ss]=load_lidar(file_sun_small,FILTER_ZEROS,RANGE_LIMIT,SUBSAMPLE_N);
[X_mb,Y_mb,Z_mb,R_mb]=load_lidar(file_mon_big,  FILTER_ZEROS,RANGE_LIMIT,SUBSAMPLE_N);
[X_ms,Y_ms,Z_ms,R_ms]=load_lidar(file_mon_small,FILTER_ZEROS,RANGE_LIMIT,SUBSAMPLE_N);
fprintf('  LiDAR loaded: Sun=%d+%d pts  Mon=%d+%d pts\n',numel(X_sb),numel(X_ss),numel(X_mb),numel(X_ms));

% LiDAR reflector centroids using GPS-guided range windows
function [c_near,c_far]=lidar_centroids(x,y,z,refl,r_near,r_far)
    WINDOW=1.0; TOP_PTS=20;
    rng=sqrt(x.^2+y.^2);
    nm=abs(rng-r_near)<=WINDOW; fm=abs(rng-r_far)<=WINDOW;
    function c=top_c(xw,yw,zw,rw)
        n=min(TOP_PTS,numel(xw)); [~,si]=sort(rw,'descend'); si=si(1:n);
        c=[mean(xw(si)),mean(yw(si)),mean(zw(si))];
    end
    c_near=top_c(x(nm),y(nm),z(nm),refl(nm));
    c_far =top_c(x(fm),y(fm),z(fm),refl(fm));
end

[c_sun_near,c_sun_far]=lidar_centroids(X_sb,Y_sb,Z_sb,R_sb,9.96,13.03);
[c_mon_near,c_mon_far]=lidar_centroids(X_mb,Y_mb,Z_mb,R_mb,9.49,13.63);

sensor_pts_sun=[c_sun_near(1:2);c_sun_far(1:2);0,0];
world_pts_sun =[wf1_small(1:2); wf1_big(1:2);  0,0];
sensor_pts_mon=[c_mon_near(1:2);c_mon_far(1:2);0,0];
world_pts_mon =[wf2_small(1:2); wf2_big(1:2);  wf2_lidar(1:2)];

function [R2d,t2d]=solve_2d_svd(sp,wp)
    mu_s=mean(sp,1); mu_w=mean(wp,1);
    S=(sp-mu_s)'*(wp-mu_w); [U,~,V]=svd(S);
    R2d=V*diag([1,det(V*U')])*U';
    t2d=(mu_w'-R2d*mu_s');
end

[R_sun,t_sun]=solve_2d_svd(sensor_pts_sun,world_pts_sun);
[R_mon,t_mon]=solve_2d_svd(sensor_pts_mon,world_pts_mon);
th_sun=atan2(R_sun(2,1),R_sun(1,1)); th_mon=atan2(R_mon(2,1),R_mon(1,1));
cs=cos(th_sun);ss=sin(th_sun); cm=cos(th_mon);sm=sin(th_mon);
Rsun=[cs ss;-ss cs]; Rmon=[cm sm;-sm cm];

xy_sb=[X_sb,Y_sb]*Rsun+t_sun'; xy_ss=[X_ss,Y_ss]*Rsun+t_sun';
xy_mb=[X_mb,Y_mb]*Rmon+t_mon'; xy_ms=[X_ms,Y_ms]*Rmon+t_mon';
Wsb=[xy_sb,Z_sb]; Wss=[xy_ss,Z_ss];
Wmb=[xy_mb,Z_mb]; Wms=[xy_ms,Z_ms];
fprintf('  LiDAR Sun: yaw=%.2f deg  Mon: yaw=%.2f deg\n',rad2deg(th_sun),rad2deg(th_mon));

% =========================================================================
fprintf('\n=== SECTION 4: NOISE FILTERING (RADAR) ===\n');
% =========================================================================
% Z clip
keep1=pts1_world(:,3)<=Z_MAX; keep2=pts2_world(:,3)<=Z_MAX;
pts1_world=pts1_world(keep1,:); amp1=amp1(keep1);
pts2_world=pts2_world(keep2,:); amp2=amp2(keep2);

% kNN outlier removal
function pts_clean=remove_outliers(pts,k,std_mult)
    if isempty(pts); pts_clean=pts; return; end
    idx=knnsearch(pts,pts,'K',k+1); idx=idx(:,2:end);
    dists=zeros(size(pts,1),1);
    for i=1:k; dists=dists+sqrt(sum((pts-pts(idx(:,i),:)).^2,2)); end
    mean_dist=dists/k;
    threshold=mean(mean_dist)+std_mult*std(mean_dist);
    pts_clean=pts(mean_dist<=threshold,:);
end

pts1_clean=remove_outliers(pts1_world(:,1:3),K_KNN,STD_MULTIPLIER);
pts2_clean=remove_outliers(pts2_world(:,1:3),K_KNN,STD_MULTIPLIER);
idx1=knnsearch(pts1_world(:,1:3),pts1_clean,'K',1);
idx2=knnsearch(pts2_world(:,1:3),pts2_clean,'K',1);
amp1=amp1(idx1); pts1_world=pts1_world(idx1,:);
amp2=amp2(idx2); pts2_world=pts2_world(idx2,:);

% Exclusion zones
function [pts_out,amp_out]=exclude_zones(pts,amp,zones)
    keep=true(size(pts,1),1);
    for z=1:size(zones,1)
        dist=sqrt((pts(:,2)-zones(z,1)).^2+(pts(:,1)-zones(z,2)).^2);
        keep=keep&(dist>zones(z,3));
    end
    pts_out=pts(keep,:); amp_out=amp(keep);
end

zones=[EX1_N,EX1_E,EX1_R; EX2_N,EX2_E,EX2_R];
[pts1_world,amp1]=exclude_zones(pts1_world,amp1,zones);
[pts2_world,amp2]=exclude_zones(pts2_world,amp2,zones);
fprintf('  Radar after filtering: Scan1=%d  Scan2=%d pts\n',size(pts1_world,1),size(pts2_world,1));

% =========================================================================
fprintf('\n=== SECTION 5: COMBINED OVERLAY PLOT ===\n');
% =========================================================================
figure('Name','Radar + LiDAR overlay — GPS world frame','Color','w', ...
    'Units','normalized','Position',[0.01 0.05 0.97 0.88]);

ax = axes('Color','w','XColor','k','YColor','k','ZColor','k');
hold(ax,'on');

% LiDAR — plotted first (background), grey colourmap
scatter3(ax,Wsb(:,1),Wsb(:,2),Wsb(:,3),LIDAR_POINT_SIZE,[0.6 0.6 0.9],'o','filled','MarkerFaceAlpha',0.3,'DisplayName','LiDAR Sun big');
scatter3(ax,Wss(:,1),Wss(:,2),Wss(:,3),LIDAR_POINT_SIZE,[0.6 0.6 0.9],'o','filled','MarkerFaceAlpha',0.3,'DisplayName','LiDAR Sun small');
scatter3(ax,Wmb(:,1),Wmb(:,2),Wmb(:,3),LIDAR_POINT_SIZE,[0.4 0.8 0.6],'s','filled','MarkerFaceAlpha',0.3,'DisplayName','LiDAR Mon big');
scatter3(ax,Wms(:,1),Wms(:,2),Wms(:,3),LIDAR_POINT_SIZE,[0.4 0.8 0.6],'s','filled','MarkerFaceAlpha',0.3,'DisplayName','LiDAR Mon small');

% Radar — plotted on top, jet colourmap by amplitude
scatter3(ax,pts1_world(:,1),pts1_world(:,2),pts1_world(:,3),RADAR_POINT_SIZE,amp1,'o','filled','MarkerFaceAlpha',0.6,'DisplayName','Radar Scan 1');
scatter3(ax,pts2_world(:,1),pts2_world(:,2),pts2_world(:,3),RADAR_POINT_SIZE,amp2,'s','filled','MarkerFaceAlpha',0.6,'DisplayName','Radar Scan 2');

% GPS markers
plot3(ax,wf1_big(1),  wf1_big(2),  wf1_big(3),  'r^','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Big refl GPS (Sun)');
plot3(ax,wf1_small(1),wf1_small(2),wf1_small(3),'g^','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Small refl GPS (Sun)');
plot3(ax,wf2_big(1),  wf2_big(2),  wf2_big(3),  'rv','MarkerSize',12,'MarkerFaceColor','r','DisplayName','Big refl GPS (Mon)');
plot3(ax,wf2_small(1),wf2_small(2),wf2_small(3),'gv','MarkerSize',12,'MarkerFaceColor','g','DisplayName','Small refl GPS (Mon)');
plot3(ax,0,0,0,'kx','MarkerSize',14,'LineWidth',2,'DisplayName','World origin (Sun sensor)');

colormap(ax,'jet'); clim(ax,[-100 100]);
cb=colorbar(ax); cb.Color='k'; cb.Label.String='Radar amplitude (dB)'; cb.Label.Color='k';
xlabel(ax,'East (m)','Color','k'); ylabel(ax,'North (m)','Color','k'); zlabel(ax,'Z (m)','Color','k');
title(ax,'Radar + LiDAR overlay — GPS world frame (ENU)','Color','k','FontSize',13);
grid(ax,'on'); ax.GridColor=[0.8 0.8 0.8];
axis(ax,'equal'); view(ax,45,25); rotate3d(ax,'on');
legend(ax,'show','TextColor','k','Color','w','FontSize',8,'Location','northeast');
xlim(ax, [-30 20]);
ylim(ax, [-10 60]);
zlim(ax, [-5 25]);
fprintf('  Plot complete.\nDone.\n');