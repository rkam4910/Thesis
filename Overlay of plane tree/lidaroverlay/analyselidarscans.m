% =========================================================================
% analyse_lidar_scans.m
% =========================================================================

clear; clc; close all;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------
CSV_FILES  = {'2026-03-16_smallreflector_planetree.csv', ...
              '2026-03-16_bigreflector_planetree.csv'};

POINT_SIZE   = 2;
MAX_RANGE    = 30;      % m — discard returns beyond this
FILTER_ZEROS = true;    % remove X=Y=Z=0 invalid returns
SUBSAMPLE_N  = 5;       % keep every Nth point to speed up plotting

SCAN_LABELS = {'Small reflector scan', 'Big reflector scan'};
SCAN_COLORS = {[0.2 0.5 1.0], [0.2 0.85 0.4]};

% Reflector detection — adaptive: take the top TOP_REFL_N points by
% reflectivity from a tight range window around the expected GPS distances.
% This avoids the threshold problem entirely.
TOP_REFL_N       = 30;    % number of brightest points to use for centroid
RANGE_WIN        = 1.5;   % ± m around expected range to search for reflector
R_EXPECTED_SMALL = 10.8;  % m — expected horizontal range to small reflector
R_EXPECTED_BIG   = 13.6;  % m — expected horizontal range to big reflector

% Attenuation settings
CYLINDER_RADIUS  = 0.5;   % m — cylinder radius around ray
CANOPY_Z_MIN     = 1.0;   % m — minimum Z to count as canopy (sensor-frame)
% Note: if reflector Z comes out negative (reflectors below sensor),
% the ray travels slightly downward. The script compensates by also
% searching points above the ray midpoint Z, not just above a fixed floor.

% =========================================================================
fprintf('=== SECTION 1: LOAD DATA ===\n');
% =========================================================================
scans = struct('x',{},'y',{},'z',{},'refl',{},'label',{},...
               'xf',{},'yf',{},'zf',{},'reff',{});

for k = 1:2
    fprintf('  Loading: %s\n', CSV_FILES{k});
    T = readtable(CSV_FILES{k}, 'VariableNamingRule', 'preserve');
    get_col = @(name) T{:, find(strcmpi(T.Properties.VariableNames, name), 1)};

    x = get_col('X'); y = get_col('Y'); z = get_col('Z'); refl = get_col('Reflectivity');

    if FILTER_ZEROS
        v = ~(x==0 & y==0 & z==0);
        x=x(v); y=y(v); z=z(v); refl=refl(v);
    end

    rng  = sqrt(x.^2 + y.^2 + z.^2);
    keep = rng <= MAX_RANGE;
    x=x(keep); y=y(keep); z=z(keep); refl=refl(keep);

    % Full-res for detection/attenuation
    scans(k).xf = x; scans(k).yf = y; scans(k).zf = z; scans(k).reff = refl;

    % Subsampled for plotting
    idx = 1:SUBSAMPLE_N:numel(x);
    scans(k).x = x(idx); scans(k).y = y(idx);
    scans(k).z = z(idx); scans(k).refl = refl(idx);
    scans(k).label = SCAN_LABELS{k};

    fprintf('    Full-res: %d pts  |  Plot subsample: %d pts\n', numel(x), numel(idx));
    fprintf('    X: %.3f to %.3f m\n', min(x), max(x));
    fprintf('    Y: %.3f to %.3f m\n', min(y), max(y));
    fprintf('    Z: %.3f to %.3f m\n', min(z), max(z));
    fprintf('    Reflectivity: %d to %d  (mean=%.1f)\n', min(refl), max(refl), mean(refl));
    fprintf('    Refl>=240: %d  >=220: %d  >=200: %d  >=150: %d\n', ...
        sum(refl>=240), sum(refl>=220), sum(refl>=200), sum(refl>=150));
end

% =========================================================================
fprintf('\n=== SECTION 2: CORNER REFLECTOR DETECTION ===\n');
% =========================================================================
% Strategy: use a range window around the known GPS distance to each
% reflector, then take the top TOP_REFL_N points by reflectivity within
% that window. This isolates the reflector without needing a fixed threshold.

r_expected = [R_EXPECTED_SMALL, R_EXPECTED_BIG];

for k = 1:2
    x    = scans(k).xf;
    y    = scans(k).yf;
    z    = scans(k).zf;
    refl = scans(k).reff;

    horiz_rng = sqrt(x.^2 + y.^2);
    win_mask  = abs(horiz_rng - r_expected(k)) <= RANGE_WIN;

    fprintf('  %s: %d pts in range window %.1f±%.1fm\n', ...
        scans(k).label, sum(win_mask), r_expected(k), RANGE_WIN);

    if sum(win_mask) < 1
        fprintf('    WARNING: no points in range window.\n');
        continue;
    end

    x_w = x(win_mask); y_w = y(win_mask); z_w = z(win_mask); r_w = refl(win_mask);

    % Sort by reflectivity descending, take top N
    [r_sort, si] = sort(r_w, 'descend');
    n   = min(TOP_REFL_N, numel(r_sort));
    si  = si(1:n);
    cx  = mean(x_w(si)); cy = mean(y_w(si)); cz = mean(z_w(si));
    hr  = sqrt(cx^2 + cy^2);
    sr  = sqrt(cx^2 + cy^2 + cz^2);

    fprintf('    Top-%d reflectivities: max=%d  median=%d  min=%d\n', ...
        n, r_sort(1), median(r_sort(1:n)), r_sort(n));
    fprintf('    Centroid: X=%.3f  Y=%.3f  Z=%.3f m\n', cx, cy, cz);
    fprintf('    Horizontal range: %.3f m  |  Slant range: %.3f m\n', hr, sr);

    scans(k).centroid    = [cx, cy, cz];
    scans(k).horiz_range = hr;
    scans(k).slant_range = sr;
    scans(k).refl_median = median(r_sort(1:n));
    scans(k).refl_mean   = mean(r_sort(1:n));
    scans(k).refl_max    = r_sort(1);
    scans(k).n_refl_pts  = n;

    % Reflector mask on subsampled data (for plotting) — top 1% of refl in window
    refl_thresh_plot = prctile(scans(k).refl, 99);
    scans(k).refl_mask = scans(k).refl >= refl_thresh_plot;
end

% =========================================================================
fprintf('\n=== SECTION 3: PLOTS ===\n');
% =========================================================================

fig1 = figure('Name','LiDAR — 3D overlay','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.05 0.1 0.65 0.8]);
ax1 = axes('Parent',fig1,'Color','k','XColor','w','YColor','w','ZColor','w');
hold(ax1,'on');
scatter3(ax1,scans(1).x,scans(1).y,scans(1).z,POINT_SIZE,scans(1).refl, ...
    'o','filled','MarkerFaceAlpha',0.6,'DisplayName',scans(1).label);
scatter3(ax1,scans(2).x,scans(2).y,scans(2).z,POINT_SIZE,scans(2).refl, ...
    's','filled','MarkerFaceAlpha',0.6,'DisplayName',scans(2).label);
for k=1:2
    if isfield(scans(k),'centroid')
        C = scans(k).centroid;
        scatter3(ax1,C(1),C(2),C(3),200,SCAN_COLORS{k},'^','filled', ...
            'DisplayName',[scans(k).label ' centroid']);
    end
end
colormap(fig1,'turbo');
cb=colorbar(ax1); cb.Color='w'; cb.Label.String='Reflectivity'; cb.Label.Color='w';
xlabel(ax1,'X (m)','Color','w'); ylabel(ax1,'Y (m)','Color','w'); zlabel(ax1,'Z (m)','Color','w');
title(ax1,'3D overlay','Color','w','FontSize',12);
axis(ax1,'equal'); grid(ax1,'on'); ax1.GridColor=[0.3 0.3 0.3];
view(ax1,45,30); rotate3d(ax1,'on');
legend(ax1,'show','TextColor','w','Color','k','FontSize',9);

fig2 = figure('Name','Top-down','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.35 0.55 0.3 0.4]);
ax2 = axes('Parent',fig2,'Color','k','XColor','w','YColor','w'); hold(ax2,'on');
scatter(ax2,scans(1).x,scans(1).y,POINT_SIZE,scans(1).refl,'o','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName',scans(1).label);
scatter(ax2,scans(2).x,scans(2).y,POINT_SIZE,scans(2).refl,'s','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName',scans(2).label);
colormap(fig2,'turbo'); grid(ax2,'on'); ax2.GridColor=[0.3 0.3 0.3]; axis(ax2,'equal');
xlabel(ax2,'X (m)','Color','w'); ylabel(ax2,'Y (m)','Color','w');
title(ax2,'Top-down (XY)','Color','w','FontSize',13);
legend(ax2,'show','TextColor','w','Color','k','FontSize',9);

fig3 = figure('Name','Side view','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.65 0.55 0.3 0.4]);
ax3 = axes('Parent',fig3,'Color','k','XColor','w','YColor','w'); hold(ax3,'on');
scatter(ax3,scans(1).x,scans(1).z,POINT_SIZE,scans(1).refl,'o','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName',scans(1).label);
scatter(ax3,scans(2).x,scans(2).z,POINT_SIZE,scans(2).refl,'s','filled', ...
    'MarkerFaceAlpha',0.5,'DisplayName',scans(2).label);
colormap(fig3,'turbo'); grid(ax3,'on'); ax3.GridColor=[0.3 0.3 0.3]; axis(ax3,'equal');
xlabel(ax3,'X (m)','Color','w'); ylabel(ax3,'Z (m)','Color','w');
title(ax3,'Side view (XZ)','Color','w','FontSize',13);
legend(ax3,'show','TextColor','w','Color','k','FontSize',9);

fig4 = figure('Name','Reflectivity distributions','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.05 0.55 0.28 0.4]);
for k=1:2
    ax=subplot(1,2,k); set(ax,'Color','k','XColor','w','YColor','w');
    histogram(ax,scans(k).refl,64,'FaceColor',SCAN_COLORS{k},'EdgeColor','none');
    if isfield(scans(k),'refl_median')
        xline(scans(k).refl_median,'r--','LineWidth',1.5);
    end
    xlabel(ax,'Reflectivity','Color','w'); ylabel(ax,'Count','Color','w');
    title(ax,scans(k).label,'Color','w','FontSize',11);
    grid(ax,'on'); ax.GridColor=[0.3 0.3 0.3];
end
sgtitle('Reflectivity distributions (red=reflector median)','Color','w','FontSize',12);

fig5 = figure('Name','Height distributions','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.35 0.1 0.28 0.4]);
for k=1:2
    ax=subplot(1,2,k); set(ax,'Color','k','XColor','w','YColor','w');
    histogram(ax,scans(k).z,60,'FaceColor',SCAN_COLORS{k},'EdgeColor','none');
    xline(0,'w--','LineWidth',1);
    xline(CANOPY_Z_MIN,'y--','LineWidth',1);
    if isfield(scans(k),'centroid')
        xline(scans(k).centroid(3),'r--','LineWidth',1.2);
    end
    xlabel(ax,'Z (m)','Color','w'); ylabel(ax,'Count','Color','w');
    title(ax,scans(k).label,'Color','w','FontSize',11);
    grid(ax,'on'); ax.GridColor=[0.3 0.3 0.3];
end
sgtitle('Z distributions  white=0  yellow=canopy min  red=reflector Z', ...
    'Color','w','FontSize',11);

fprintf('  5 figures generated.\n');

% =========================================================================
fprintf('\n=== SECTION 4: ATTENUATION ANALYSIS ===\n');
% =========================================================================
% The ray from sensor (0,0,0) to each reflector centroid travels mostly
% horizontally with a slight downward tilt (reflectors ~1m below sensor).
% We search the pooled full-res cloud for points that:
%   (a) fall within CYLINDER_RADIUS of the ray,
%   (b) lie between 5% and 95% along the ray (exclude ends),
%   (c) have Z > CANOPY_Z_MIN in sensor frame.
% These are canopy points the beam path passed through or near.
% We also report all-Z cylinder counts so you can see if the zero-intercept
% issue is a geometry problem vs a genuine clear path.

function [pd, tp] = ray_distances(qx, qy, qz, centroid)
    P  = centroid(:)'; Pm = norm(P); d = P/Pm;
    Q  = [qx(:), qy(:), qz(:)];
    proj = Q*d'; Q_perp = Q - proj*d;
    pd = sqrt(sum(Q_perp.^2,2)); tp = proj/Pm;
end

pool_x = [scans(1).xf; scans(2).xf];
pool_y = [scans(1).yf; scans(2).yf];
pool_z = [scans(1).zf; scans(2).zf];

fprintf('  Pooled cloud: %d pts\n', numel(pool_x));
fprintf('  Cylinder r=%.2fm  canopy Z>%.2fm\n\n', CYLINDER_RADIUS, CANOPY_Z_MIN);

for k = 1:2
    if ~isfield(scans(k),'centroid')
        fprintf('  %s: no centroid, skipping.\n', scans(k).label);
        continue;
    end
    C = scans(k).centroid;
    [pd, tp] = ray_distances(pool_x, pool_y, pool_z, C);

    in_cyl   = pd <= CYLINDER_RADIUS;
    in_seg   = tp >= 0.05 & tp <= 0.95;
    in_both  = in_cyl & in_seg;

    % Break down by Z band so we can diagnose geometry
    n_all    = sum(in_both);
    n_canopy = sum(in_both & pool_z > CANOPY_Z_MIN);
    n_ground = sum(in_both & pool_z <= CANOPY_Z_MIN);
    n_below  = sum(in_both & pool_z < 0);

    canopy_frac = n_canopy / max(n_all, 1);

    scans(k).n_intercept     = n_canopy;
    scans(k).n_cyl_total     = n_all;
    scans(k).canopy_fraction = canopy_frac;
    scans(k).intercept_mask  = in_both & pool_z > CANOPY_Z_MIN;

    fprintf('  %s:\n', scans(k).label);
    fprintf('    Centroid: X=%.3f  Y=%.3f  Z=%.3f  (slant=%.2fm)\n', C, scans(k).slant_range);
    fprintf('    Ray direction: horiz bearing=%.1f deg  elevation=%.2f deg\n', ...
        rad2deg(atan2(C(2),C(1))), rad2deg(atan2(C(3), sqrt(C(1)^2+C(2)^2))));
    fprintf('    Cylinder pts: %d total  |  %d canopy(Z>%.1f)  |  %d low(Z<=%.1f)  |  %d below 0\n', ...
        n_all, n_canopy, CANOPY_Z_MIN, n_ground, CANOPY_Z_MIN, n_below);
    fprintf('    Canopy fraction: %.3f\n', canopy_frac);
    fprintf('    Reflector: median=%d  mean=%.1f  max=%d  (from %d pts)\n', ...
        scans(k).refl_median, scans(k).refl_mean, scans(k).refl_max, scans(k).n_refl_pts);
end

% ── Reference selection and attenuation ──────────────────────────────────
if isfield(scans(1),'refl_median') && isfield(scans(2),'refl_median')
    if scans(1).refl_median >= scans(2).refl_median
        ref_idx=1; att_idx=2;
    else
        ref_idx=2; att_idx=1;
    end
    R_ref = scans(ref_idx).refl_median;
    R_att = scans(att_idx).refl_median;
    attenuation = 1 - R_att/R_ref;

    fprintf('\n  Reference: %s  (median=%d  canopy frac=%.3f)\n', ...
        scans(ref_idx).label, R_ref, scans(ref_idx).canopy_fraction);
    fprintf('  Attenuated: %s  (median=%d  canopy frac=%.3f)\n', ...
        scans(att_idx).label, R_att, scans(att_idx).canopy_fraction);
    fprintf('\n  Attenuation: %.1f%%\n', attenuation*100);
    fprintf('  (If canopy fraction is 0 for both, the cylinder radius may need\n');
    fprintf('   increasing, or the sensor may be scanning from the same side as\n');
    fprintf('   the reflectors with minimal canopy in between.)\n');
else
    attenuation = 0;
    ref_idx = 1; att_idx = 2;
end

% ── Figures 6 & 7: ray visualisation ─────────────────────────────────────
fig6 = figure('Name','Ray interception — top-down','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.05 0.05 0.55 0.85]);
ax6 = axes('Parent',fig6,'Color','k','XColor','w','YColor','w','ZColor','w');
hold(ax6,'on');
scatter(ax6,pool_x,pool_y,1,[0.35 0.35 0.35],'filled', ...
    'MarkerFaceAlpha',0.2,'DisplayName','All pts');
for k=1:2
    if ~isfield(scans(k),'intercept_mask'), continue; end
    im = scans(k).intercept_mask;
    if any(im)
        scatter(ax6,pool_x(im),pool_y(im),25,SCAN_COLORS{k},'filled', ...
            'MarkerFaceAlpha',0.9,'DisplayName', ...
            sprintf('%s intercepts (%d)', scans(k).label, sum(im)));
    end
    C = scans(k).centroid;
    tag = ''; if k==ref_idx, tag=' (REF)'; elseif k==att_idx, tag=' (ATT)'; end
    plot(ax6,[0 C(1)],[0 C(2)],'--','Color',SCAN_COLORS{k},'LineWidth',1.5, ...
        'DisplayName',['Ray: ' scans(k).label tag]);
    scatter(ax6,C(1),C(2),120,SCAN_COLORS{k},'^','filled', ...
        'DisplayName',['Reflector: ' scans(k).label]);
end
scatter(ax6,0,0,100,'w','x','LineWidth',2,'DisplayName','Sensor');
axis(ax6,'equal'); grid(ax6,'on'); ax6.GridColor=[0.3 0.3 0.3];
xlabel(ax6,'X (m)','Color','w'); ylabel(ax6,'Y (m)','Color','w');
title(ax6,sprintf('Top-down ray view — attenuation %.1f%%\ncylinder r=%.2fm', ...
    attenuation*100, CYLINDER_RADIUS),'Color','w','FontSize',11);
legend(ax6,'show','TextColor','w','Color','k','FontSize',8,'Location','best');

fig7 = figure('Name','Ray interception — side view','NumberTitle','off', ...
    'Color','k','Units','normalized','Position',[0.62 0.05 0.36 0.85]);
ax7 = axes('Parent',fig7,'Color','k','XColor','w','YColor','w');
hold(ax7,'on');
scatter(ax7,pool_x,pool_z,1,[0.35 0.35 0.35],'filled', ...
    'MarkerFaceAlpha',0.2,'DisplayName','All pts');
for k=1:2
    if ~isfield(scans(k),'intercept_mask'), continue; end
    im = scans(k).intercept_mask;
    if any(im)
        scatter(ax7,pool_x(im),pool_z(im),25,SCAN_COLORS{k},'filled', ...
            'MarkerFaceAlpha',0.9,'DisplayName',[scans(k).label ' intercepts']);
    end
    C = scans(k).centroid;
    plot(ax7,[0 C(1)],[0 C(3)],'--','Color',SCAN_COLORS{k},'LineWidth',1.5, ...
        'DisplayName',['Ray: ' scans(k).label]);
    scatter(ax7,C(1),C(3),120,SCAN_COLORS{k},'^','filled', ...
        'DisplayName',['Reflector: ' scans(k).label]);
end
yline(ax7,CANOPY_Z_MIN,'y--','LineWidth',1.2,'DisplayName',sprintf('Canopy min Z=%.1f',CANOPY_Z_MIN));
yline(ax7,0,'w--','LineWidth',0.8,'DisplayName','Sensor Z=0');
scatter(ax7,0,0,100,'w','x','LineWidth',2,'DisplayName','Sensor');
axis(ax7,'equal'); grid(ax7,'on'); ax7.GridColor=[0.3 0.3 0.3];
xlabel(ax7,'X (m)','Color','w'); ylabel(ax7,'Z (m)','Color','w');
title(ax7,'Side view XZ — ray vs canopy','Color','w','FontSize',11);
legend(ax7,'show','TextColor','w','Color','k','FontSize',8,'Location','best');

fprintf('\n  7 figures generated.\n\nDone.\n');