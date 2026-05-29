% =========================================================================
% analyse_lidar_scans.m
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
CSV_FILES  = {'2026-03-16_smallreflector_planetree.csv', ...
              '2026-03-16_bigreflector_planetree.csv', ...
              '2026-03-15 14-50-46_corner_big.csv', ...
              '2026-03-15 14-52-14_corner_small.csv'};

POINT_SIZE   = 2;
MAX_RANGE    = 30;
FILTER_ZEROS = true;
SUBSAMPLE_N  = 5;

% Labels and colours: Monday scans 1-2, Sunday scans 3-4
SCAN_LABELS = {'Mon small refl', 'Mon big refl', 'Sun big refl', 'Sun small refl'};
SCAN_COLORS = {[0.2 0.5 1.0], [0.2 0.85 0.4], [1.0 0.5 0.1], [1.0 0.9 0.1]};

% Reflector detection
TOP_REFL_N       = 30;
RANGE_WIN        = 1.5;
% Expected horizontal ranges [Mon small, Mon big, Sun big, Sun small]
R_EXPECTED       = [10.8, 13.6, 13.0, 9.96];

% Tree bounding box (sensor frame) — applied to all scans
TREE_X_MIN =  5.0;   TREE_X_MAX = 25.0;
TREE_Y_MIN = -8.0;   TREE_Y_MAX =  6.0;
TREE_Z_MIN =  0.5;   TREE_Z_MAX = 20.0;

% Conical section
CONE_HALF_ANGLE_DEG = 25;   % widened from 15 to capture more canopy

% Distance bins
BIN_START  =  2.0;
BIN_END    = 22.0;
BIN_WIDTH  =  1.0;

% =========================================================================
fprintf('=== SECTION 1: LOAD DATA ===\n');
% =========================================================================
n_scans = numel(CSV_FILES);
scans = struct('x',{},'y',{},'z',{},'refl',{},'label',{},...
               'xf',{},'yf',{},'zf',{},'reff',{});

for k = 1:n_scans
    fprintf('  Loading: %s\n', CSV_FILES{k});
    T = readtable(CSV_FILES{k}, 'VariableNamingRule', 'preserve');
    get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};

    x=get_col('X'); y=get_col('Y'); z=get_col('Z'); refl=get_col('Reflectivity');

    if FILTER_ZEROS
        v=~(x==0&y==0&z==0); x=x(v); y=y(v); z=z(v); refl=refl(v);
    end
    rng=sqrt(x.^2+y.^2+z.^2); keep=rng<=MAX_RANGE;
    x=x(keep); y=y(keep); z=z(keep); refl=refl(keep);

    scans(k).xf=x; scans(k).yf=y; scans(k).zf=z; scans(k).reff=refl;
    idx=1:SUBSAMPLE_N:numel(x);
    scans(k).x=x(idx); scans(k).y=y(idx);
    scans(k).z=z(idx); scans(k).refl=refl(idx);
    scans(k).label=SCAN_LABELS{k};

    fprintf('    Full-res: %d pts  |  Subsample: %d pts\n',numel(x),numel(idx));
    fprintf('    X: %.2f to %.2f  Y: %.2f to %.2f  Z: %.2f to %.2f m\n', ...
        min(x),max(x),min(y),max(y),min(z),max(z));
    fprintf('    Refl: mean=%.1f  max=%d  >=200:%d\n',mean(refl),max(refl),sum(refl>=200));
end

% =========================================================================
fprintf('\n=== SECTION 2: REFLECTOR DETECTION ===\n');
% =========================================================================
for k = 1:n_scans
    x=scans(k).xf; y=scans(k).yf; z=scans(k).zf; refl=scans(k).reff;
    hr=sqrt(x.^2+y.^2);
    win=abs(hr-R_EXPECTED(k))<=RANGE_WIN;

    x_w=x(win); y_w=y(win); z_w=z(win); r_w=refl(win);
    [rs,si]=sort(r_w,'descend'); n=min(TOP_REFL_N,numel(rs)); si=si(1:n);
    cx=mean(x_w(si)); cy=mean(y_w(si)); cz=mean(z_w(si));

    az=atan2(cy,cx); el=atan2(cz,sqrt(cx^2+cy^2));

    fprintf('  %s: max=%d median=%d | X=%.2f Y=%.2f Z=%.2f | az=%.1fdeg el=%.1fdeg\n', ...
        scans(k).label,rs(1),median(rs(1:n)),cx,cy,cz,rad2deg(az),rad2deg(el));

    scans(k).centroid    = [cx,cy,cz];
    scans(k).horiz_range = sqrt(cx^2+cy^2);
    scans(k).slant_range = sqrt(cx^2+cy^2+cz^2);
    scans(k).cone_az     = az;
    scans(k).cone_el     = el;
    scans(k).refl_median = median(rs(1:n));
    scans(k).refl_mean   = mean(rs(1:n));
    scans(k).refl_max    = rs(1);
    scans(k).refl_mask   = scans(k).refl >= prctile(scans(k).refl,99);
end

% =========================================================================
fprintf('\n=== SECTION 3: STANDARD PLOTS ===\n');
% =========================================================================

fig1=figure('Name','3D overlay — all scans','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.02 0.1 0.65 0.8]);
ax1=axes('Parent',fig1,'Color','k','XColor','w','YColor','w','ZColor','w');
hold(ax1,'on');
markers = {'o','s','^','d'};
for k=1:n_scans
    scatter3(ax1,scans(k).x,scans(k).y,scans(k).z,POINT_SIZE, ...
        repmat(reshape(SCAN_COLORS{k},1,3),numel(scans(k).x),1),'filled', ...
        'MarkerFaceAlpha',0.4,'DisplayName',scans(k).label);
end
for k=1:n_scans
    if isfield(scans(k),'centroid'), C=scans(k).centroid;
        scatter3(ax1,C(1),C(2),C(3),200,SCAN_COLORS{k},'^','filled', ...
            'DisplayName',[scans(k).label ' refl']); end
end
xlabel(ax1,'X (m)','Color','w'); ylabel(ax1,'Y (m)','Color','w'); zlabel(ax1,'Z (m)','Color','w');
title(ax1,'3D overlay — all 4 scans','Color','w','FontSize',12);
axis(ax1,'equal'); grid(ax1,'on'); ax1.GridColor=[0.3 0.3 0.3];
view(ax1,45,30); rotate3d(ax1,'on');
legend(ax1,'show','TextColor','w','Color','k','FontSize',8,'Location','northeast');

fig2=figure('Name','Top-down — all scans','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.35 0.55 0.3 0.4]);
ax2=axes('Parent',fig2,'Color','k','XColor','w','YColor','w'); hold(ax2,'on');
for k=1:n_scans
    scatter(ax2,scans(k).x,scans(k).y,POINT_SIZE, ...
        repmat(reshape(SCAN_COLORS{k},1,3),numel(scans(k).x),1),'filled', ...
        'MarkerFaceAlpha',0.4,'DisplayName',scans(k).label);
end
rectangle(ax2,'Position',[TREE_X_MIN TREE_Y_MIN ...
    TREE_X_MAX-TREE_X_MIN TREE_Y_MAX-TREE_Y_MIN], ...
    'EdgeColor','y','LineStyle','--','LineWidth',1.5);
for k=1:n_scans
    if ~isfield(scans(k),'cone_az'), continue; end
    theta=linspace(-deg2rad(CONE_HALF_ANGLE_DEG),deg2rad(CONE_HALF_ANGLE_DEG),60);
    cx_c=BIN_END*cos(scans(k).cone_az+theta);
    cy_c=BIN_END*sin(scans(k).cone_az+theta);
    plot(ax2,[0 cx_c 0],[0 cy_c 0],'--','Color',SCAN_COLORS{k},'LineWidth',1, ...
        'DisplayName',['Cone: ' scans(k).label]);
end
grid(ax2,'on'); ax2.GridColor=[0.3 0.3 0.3]; axis(ax2,'equal');
xlabel(ax2,'X (m)','Color','w'); ylabel(ax2,'Y (m)','Color','w');
title(ax2,'Top-down — cones & tree box','Color','w','FontSize',12);
legend(ax2,'show','TextColor','w','Color','k','FontSize',7,'Location','best');

fprintf('  2 overview figures generated.\n');

% =========================================================================
fprintf('\n=== SECTION 4: CONICAL ATTENUATION ANALYSIS ===\n');
% =========================================================================
bin_edges   = BIN_START:BIN_WIDTH:BIN_END;
bin_centres = bin_edges(1:end-1) + BIN_WIDTH/2;
n_bins      = numel(bin_centres);
cone_rad    = deg2rad(CONE_HALF_ANGLE_DEG);

% Cone shell volumes for density normalisation
cone_vol = (2/3)*pi*(1-cos(cone_rad)) * ...
    ((bin_edges(2:end)).^3 - (bin_edges(1:end-1)).^3);

fprintf('  Cone half-angle: %.0f deg  bins: %.0f-%.0f m  width=%.1f m\n', ...
    CONE_HALF_ANGLE_DEG, BIN_START, BIN_END, BIN_WIDTH);

for k = 1:n_scans
    if ~isfield(scans(k),'cone_az'), continue; end

    x=scans(k).xf; y=scans(k).yf; z=scans(k).zf; refl=scans(k).reff;
    slant=sqrt(x.^2+y.^2+z.^2);
    az_pt=atan2(y,x); el_pt=atan2(z,sqrt(x.^2+y.^2));

    daz=mod(az_pt-scans(k).cone_az+pi,2*pi)-pi;
    del=el_pt-scans(k).cone_el;
    ang_off=sqrt(daz.^2+del.^2);

    in_box=x>=TREE_X_MIN&x<=TREE_X_MAX & y>=TREE_Y_MIN&y<=TREE_Y_MAX & ...
           z>=TREE_Z_MIN&z<=TREE_Z_MAX;
    in_cone=ang_off<=cone_rad;
    in_sample=in_cone&in_box;

    fprintf('\n  %s: %d pts in cone+box\n',scans(k).label,sum(in_sample));

    bin_count=zeros(1,n_bins);
    bin_refl=nan(1,n_bins);
    bin_refl_hi=nan(1,n_bins);

    for b=1:n_bins
        in_bin=in_sample & slant>=bin_edges(b) & slant<bin_edges(b+1);
        n=sum(in_bin);
        bin_count(b)=n;
        if n>0
            r_bin=refl(in_bin);
            bin_refl(b)=mean(r_bin);
            thr=prctile(r_bin,90);
            hi=r_bin>=thr;
            if any(hi), bin_refl_hi(b)=mean(r_bin(hi)); end
        end
    end

    bin_density=bin_count./cone_vol;

    scans(k).bin_centres  = bin_centres;
    scans(k).bin_count    = bin_count;
    scans(k).bin_density  = bin_density;
    scans(k).bin_refl     = bin_refl;
    scans(k).bin_refl_hi  = bin_refl_hi;
    scans(k).in_cone_mask = in_sample;

    fprintf('    Range(m)  Count  Density  MeanRefl  Top10%%\n');
    for b=1:n_bins
        if bin_count(b)>0
            fprintf('    %5.1f  %6d  %7.2f  %6.1f  %6.1f\n', ...
                bin_centres(b),bin_count(b),bin_density(b),bin_refl(b),bin_refl_hi(b));
        end
    end
end

% ── Per-bin ratios: Mon small / Mon big, Sun big / Sun small ─────────────
fprintf('\n  Mon small/big ratio per bin:\n');
fprintf('  Range  Refl_s  Refl_b  Ratio  Dens_s  Dens_b\n');
ratio_mon = nan(1,n_bins);
for b=1:n_bins
    r1=scans(1).bin_refl(b); r2=scans(2).bin_refl(b);
    if ~isnan(r1)&&~isnan(r2)&&r2>0
        ratio_mon(b)=r1/r2;
        fprintf('  %5.1f  %6.1f  %6.1f  %5.3f  %6.2f  %6.2f\n', ...
            bin_centres(b),r1,r2,ratio_mon(b), ...
            scans(1).bin_density(b),scans(2).bin_density(b));
    end
end

fprintf('\n  Sun big/small ratio per bin:\n');
fprintf('  Range  Refl_b  Refl_s  Ratio  Dens_b  Dens_s\n');
ratio_sun = nan(1,n_bins);
for b=1:n_bins
    r3=scans(3).bin_refl(b); r4=scans(4).bin_refl(b);
    if ~isnan(r3)&&~isnan(r4)&&r4>0
        ratio_sun(b)=r3/r4;
        fprintf('  %5.1f  %6.1f  %6.1f  %5.3f  %6.2f  %6.2f\n', ...
            bin_centres(b),r3,r4,ratio_sun(b), ...
            scans(3).bin_density(b),scans(4).bin_density(b));
    end
end

% ── Cumulative attenuation in dB ─────────────────────────────────────────
fprintf('\n  Cumulative attenuation (dB) — Mon small relative to Mon big:\n');
fprintf('  (negative = small scan weaker = more attenuated)\n');
fprintf('  Range  Ratio  dB_per_bin  Cumulative_dB\n');
cum_db_mon = 0;
for b=1:n_bins
    if ~isnan(ratio_mon(b)) && ratio_mon(b)>0
        db_bin = 10*log10(ratio_mon(b));
        cum_db_mon = cum_db_mon + db_bin;
        fprintf('  %5.1f  %5.3f  %+8.2f dB  %+8.2f dB\n', ...
            bin_centres(b), ratio_mon(b), db_bin, cum_db_mon);
    end
end

fprintf('\n  Cumulative attenuation (dB) — Sun big relative to Sun small:\n');
cum_db_sun = 0;
for b=1:n_bins
    if ~isnan(ratio_sun(b)) && ratio_sun(b)>0
        db_bin = 10*log10(ratio_sun(b));
        cum_db_sun = cum_db_sun + db_bin;
        fprintf('  %5.1f  %5.3f  %+8.2f dB  %+8.2f dB\n', ...
            bin_centres(b), ratio_sun(b), db_bin, cum_db_sun);
    end
end

% =========================================================================
fprintf('\n=== SECTION 5: ATTENUATION PLOTS ===\n');
% =========================================================================

%% Figure 3: Point density vs distance — all 4 scans
fig3=figure('Name','Point density vs distance','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.02 0.05 0.47 0.85]);
ax3=axes('Parent',fig3,'Color','k','XColor','w','YColor','w');
hold(ax3,'on');
ls={'-','--','-','--'};
for k=1:n_scans
    if ~isfield(scans(k),'bin_density'), continue; end
    plot(ax3,scans(k).bin_centres,scans(k).bin_density, ...
        [ls{k} 'o'],'Color',SCAN_COLORS{k},'LineWidth',2,'MarkerSize',5, ...
        'DisplayName',scans(k).label);
end
xlabel(ax3,'Slant range (m)','Color','w');
ylabel(ax3,'Point density (pts/m^3)','Color','w');
title(ax3,sprintf('Point density in cone (half-angle=%.0f deg)',CONE_HALF_ANGLE_DEG), ...
    'Color','w','FontSize',12);
grid(ax3,'on'); ax3.GridColor=[0.3 0.3 0.3]; ax3.XLim=[BIN_START BIN_END];
legend(ax3,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

%% Figure 4: Mean reflectivity vs distance — all 4 scans
fig4=figure('Name','Reflectivity vs distance','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.51 0.05 0.47 0.85]);

ax4a=subplot(3,1,1,'Parent',fig4);
set(ax4a,'Color','k','XColor','w','YColor','w'); hold(ax4a,'on');
for k=1:n_scans
    if ~isfield(scans(k),'bin_refl'), continue; end
    plot(ax4a,scans(k).bin_centres,scans(k).bin_refl, ...
        [ls{k} 'o'],'Color',SCAN_COLORS{k},'LineWidth',2,'MarkerSize',5, ...
        'DisplayName',[scans(k).label ' mean']);
end
ylabel(ax4a,'Mean reflectivity','Color','w');
title(ax4a,'Mean reflectivity vs distance — all scans','Color','w','FontSize',11);
grid(ax4a,'on'); ax4a.GridColor=[0.3 0.3 0.3];
ax4a.XLim=[BIN_START BIN_END]; ax4a.YLim=[0 255];
legend(ax4a,'show','TextColor','w','Color','k','FontSize',8,'Location','northwest');

ax4b=subplot(3,1,2,'Parent',fig4);
set(ax4b,'Color','k','XColor','w','YColor','w'); hold(ax4b,'on');
valid_mon=~isnan(ratio_mon);
valid_sun=~isnan(ratio_sun);
if any(valid_mon)
    plot(ax4b,bin_centres(valid_mon),ratio_mon(valid_mon), ...
        'b-o','LineWidth',2,'MarkerSize',6,'DisplayName','Mon small/big');
end
if any(valid_sun)
    plot(ax4b,bin_centres(valid_sun),ratio_sun(valid_sun), ...
        'y--s','LineWidth',2,'MarkerSize',6,'DisplayName','Sun big/small');
end
yline(ax4b,1,'w--','LineWidth',1.2,'DisplayName','Ratio=1');
ylabel(ax4b,'Reflectivity ratio','Color','w');
title(ax4b,'Reflectivity ratio per bin (Mon: small/big  Sun: big/small)', ...
    'Color','w','FontSize',11);
grid(ax4b,'on'); ax4b.GridColor=[0.3 0.3 0.3]; ax4b.XLim=[BIN_START BIN_END];
legend(ax4b,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

ax4c=subplot(3,1,3,'Parent',fig4);
set(ax4c,'Color','k','XColor','w','YColor','w'); hold(ax4c,'on');
% Build cumulative dB arrays for plotting
cum_mon_arr = nan(1,n_bins); cum_sun_arr = nan(1,n_bins);
cum=0;
for b=1:n_bins
    if ~isnan(ratio_mon(b))&&ratio_mon(b)>0
        cum=cum+10*log10(ratio_mon(b));
        cum_mon_arr(b)=cum;
    end
end
cum=0;
for b=1:n_bins
    if ~isnan(ratio_sun(b))&&ratio_sun(b)>0
        cum=cum+10*log10(ratio_sun(b));
        cum_sun_arr(b)=cum;
    end
end
valid_cm=~isnan(cum_mon_arr); valid_cs=~isnan(cum_sun_arr);
if any(valid_cm)
    plot(ax4c,bin_centres(valid_cm),cum_mon_arr(valid_cm), ...
        'b-o','LineWidth',2,'MarkerSize',6,'DisplayName','Mon cumulative dB');
end
if any(valid_cs)
    plot(ax4c,bin_centres(valid_cs),cum_sun_arr(valid_cs), ...
        'y--s','LineWidth',2,'MarkerSize',6,'DisplayName','Sun cumulative dB');
end
yline(ax4c,0,'w--','LineWidth',1,'DisplayName','0 dB (no loss)');
xlabel(ax4c,'Slant range (m)','Color','w');
ylabel(ax4c,'Cumulative attenuation (dB)','Color','w');
title(ax4c,'Cumulative attenuation through canopy','Color','w','FontSize',11);
grid(ax4c,'on'); ax4c.GridColor=[0.3 0.3 0.3]; ax4c.XLim=[BIN_START BIN_END];
legend(ax4c,'show','TextColor','w','Color','k','FontSize',9,'Location','northwest');

%% Figure 5: 3D cone points coloured by distance bin
fig5=figure('Name','Cone sections — 3D by distance','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.02 0.05 0.65 0.85]);
ax5=axes('Parent',fig5,'Color','k','XColor','w','YColor','w','ZColor','w');
hold(ax5,'on');
% Background faint cloud from scan 1
scatter3(ax5,scans(1).x,scans(1).y,scans(1).z,1,[0.25 0.25 0.25],'filled', ...
    'MarkerFaceAlpha',0.1,'DisplayName','Background');
% For each scan, plot cone points coloured by slant range
for k=1:n_scans
    if ~isfield(scans(k),'in_cone_mask'), continue; end
    im=scans(k).in_cone_mask;
    if ~any(im), continue; end
    xp=scans(k).xf(im); yp=scans(k).yf(im); zp=scans(k).zf(im);
    slant_p=sqrt(xp.^2+yp.^2+zp.^2);
    % Subsample for speed
    ss=1:3:numel(xp);
    scatter3(ax5,xp(ss),yp(ss),zp(ss),4,slant_p(ss),'filled', ...
        'MarkerFaceAlpha',0.7,'DisplayName',[scans(k).label ' (by range)']);
    % Draw cone axis
    C=scans(k).centroid; Cv=C/norm(C)*BIN_END;
    plot3(ax5,[0 Cv(1)],[0 Cv(2)],[0 Cv(3)],'--','Color',SCAN_COLORS{k}, ...
        'LineWidth',2,'DisplayName',['Axis: ' scans(k).label]);
    scatter3(ax5,C(1),C(2),C(3),150,SCAN_COLORS{k},'^','filled', ...
        'DisplayName',['Reflector: ' scans(k).label]);
end
colormap(ax5,'turbo');
cb=colorbar(ax5); cb.Color='w';
cb.Label.String='Slant range (m)'; cb.Label.Color='w';
xlabel(ax5,'X (m)','Color','w'); ylabel(ax5,'Y (m)','Color','w');
zlabel(ax5,'Z (m)','Color','w');
title(ax5,sprintf('Cone sections — points coloured by distance\nhalf-angle=%.0f deg', ...
    CONE_HALF_ANGLE_DEG),'Color','w','FontSize',12);
axis(ax5,'equal'); grid(ax5,'on'); ax5.GridColor=[0.3 0.3 0.3];
view(ax5,45,30); rotate3d(ax5,'on');
legend(ax5,'show','TextColor','w','Color','k','FontSize',8,'Location','northeast');

fprintf('  5 figures generated.\n\nDone.\n');