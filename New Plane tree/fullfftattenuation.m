%% Full CFAR Point Cloud + Peak-to-End Attenuation — All Beams
clear; close all;

%% ── LOAD DATA ────────────────────────────────────────────────────────────
fid       = fopen('pointcloud_20260525-125755_0_fft_cfar.csv', 'r');
raw_lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
raw_lines = raw_lines{1};
n_rows    = length(raw_lines);

raw_el      = zeros(n_rows, 1);
raw_az      = zeros(n_rows, 1);
fft_data    = cell(n_rows, 1);
fft_lengths = zeros(n_rows, 1);

for i = 1:n_rows
    cols = strsplit(raw_lines{i}, ',');
    if length(cols) >= 2
        raw_el(i) = str2double(cols{1});
        raw_az(i) = str2double(cols{2});
    end
    if length(cols) >= 5
        vals = str2double(strsplit(cols{5}, ';'));
        vals = vals(~isnan(vals));
        fft_data{i}    = vals;
        fft_lengths(i) = length(vals);
    end
end

%% ── PARAMETERS ───────────────────────────────────────────────────────────
bin_range     = 0.125;
canopy_start_physical = 12.7;
canopy_end_physical   = 35.2;
cut_lo        = 10.0;
cut_hi        = 45.0;
top_n         = 67;
smooth_win    = 8;
amp_thresh    = -82;
min_rsq       = 0.1;
display_floor = -90;

% CFAR
cfar_guard = 1;
cfar_train = 4;
cfar_pfa   = 0.3;
alpha      = cfar_train * (cfar_pfa^(-1/cfar_train) - 1);
fprintf('CFAR alpha: %.3f\n', alpha);

% Range axes
n_fft_bins = mode(fft_lengths);
range_full = (0:n_fft_bins-1) * bin_range + (2 * bin_range);
cut_mask   = range_full >= cut_lo & range_full <= cut_hi;
range_cut  = range_full(cut_mask);
n_cut      = length(range_cut);
in_canopy  = range_cut >= canopy_start_physical & range_cut <= canopy_end_physical;

valid_mask = fft_lengths == n_fft_bins;
valid_idx  = find(valid_mask);
n_valid    = length(valid_idx);
fprintf('Valid FFT rows: %d\n', n_valid);

%% ── MAIN LOOP ────────────────────────────────────────────────────────────
max_pts   = n_valid * top_n;
all_range = zeros(1, max_pts);
all_amp   = zeros(1, max_pts);
all_el    = zeros(1, max_pts);
all_az    = zeros(1, max_pts);
ptr       = 0;

atten_rate = nan(n_valid, 1);
atten_rsq  = nan(n_valid, 1);
atten_el   = zeros(n_valid, 1);
atten_az   = zeros(n_valid, 1);
atten_peak = nan(n_valid, 1);

fprintf('Processing %d beams...\n', n_valid);

for row = 1:n_valid

    if mod(row, 500) == 0
        fprintf('  Row %d / %d\n', row, n_valid);
    end

    ri       = valid_idx(row);
    fft_full = double(fft_data{ri});
    fft_cut  = fft_full(cut_mask);

    atten_el(row) = raw_el(ri);
    atten_az(row) = raw_az(ri);

    %% CFAR
    fft_lin  = 10.^(fft_cut / 10);
    cfar_det = false(1, n_cut);
    cfar_thr = nan(1, n_cut);

    for k = 1:n_cut
        lo    = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
        hi    = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
        train = [lo, hi];
        if isempty(train), continue; end
        noise_est   = mean(fft_lin(train));
        thresh_lin  = alpha * noise_est;
        cfar_thr(k) = 10*log10(thresh_lin);
        if fft_cut(k) > cfar_thr(k)
            cfar_det(k) = true;
        end
    end

    %% Top N detections
    det_bins = find(cfar_det);
    if ~isempty(det_bins)
        det_amps = fft_cut(det_bins);
        [~, ord] = sort(det_amps, 'descend');
        keep     = det_bins(ord(1:min(top_n, end)));
        n_keep   = length(keep);
        idx_out  = ptr+1 : ptr+n_keep;
        all_range(idx_out) = range_cut(keep);
        all_amp(idx_out)   = fft_cut(keep);
        all_el(idx_out)    = raw_el(ri);
        all_az(idx_out)    = raw_az(ri);
        ptr = ptr + n_keep;
    end

    %% Peak-to-end attenuation fit
    fft_smooth  = movmean(fft_cut, smooth_win);
    canopy_vals = fft_smooth(in_canopy);
    canopy_x    = range_cut(in_canopy);

    if sum(in_canopy) < 4, continue; end
    [~, pk_idx] = max(canopy_vals);
    if (length(canopy_vals) - pk_idx) < 4, continue; end

    x_decay = canopy_x(pk_idx:end);
    y_decay = canopy_vals(pk_idx:end);
    if (max(y_decay) - min(y_decay)) < 2, continue; end

    p      = polyfit(x_decay, y_decay, 1);
    fit_y  = polyval(p, x_decay);
    ss_res = sum((y_decay - fit_y).^2);
    ss_tot = sum((y_decay - mean(y_decay)).^2);
    if ss_tot == 0, continue; end

    atten_rate(row) = p(1);
    atten_rsq(row)  = 1 - ss_res / ss_tot;
    atten_peak(row) = canopy_x(pk_idx);
end

%% ── TRIM AND FILTER POINT CLOUD ──────────────────────────────────────────
all_range = all_range(1:ptr);
all_amp   = all_amp(1:ptr);
all_el    = all_el(1:ptr);
all_az    = all_az(1:ptr);

% Apply -90 dB filter to all arrays
amp_mask  = all_amp >= amp_thresh;
all_range = all_range(amp_mask);
all_amp   = all_amp(amp_mask);
all_el    = all_el(amp_mask);
all_az    = all_az(amp_mask);
fprintf('Point cloud: %d points after %.0f dB filter\n', length(all_amp), amp_thresh);

% Convert to Cartesian
el_sph       = (pi/2) + all_el;
[px, py, pz] = sph2cart(all_az, el_sph, all_range);

%% ── FILTER ATTENUATION RESULTS ───────────────────────────────────────────
valid_fit = ~isnan(atten_rate) & atten_rsq >= min_rsq & atten_rate < 0;
fprintf('\nBeams passing filter: %d / %d\n', sum(valid_fit), n_valid);

fit_rate = atten_rate(valid_fit);
fit_rsq  = atten_rsq(valid_fit);
fit_el   = rad2deg(atten_el(valid_fit));
fit_az   = rad2deg(atten_az(valid_fit));
fit_peak = atten_peak(valid_fit);

fprintf('Attenuation statistics:\n');
fprintf('  Mean:   %.3f dB/m\n', mean(fit_rate));
fprintf('  Median: %.3f dB/m\n', median(fit_rate));
fprintf('  Std:    %.3f dB/m\n', std(fit_rate));
fprintf('  Min:    %.3f dB/m\n', min(fit_rate));
fprintf('  Max:    %.3f dB/m\n', max(fit_rate));

%% ── FIGURE 1: 3D POINT CLOUD ─────────────────────────────────────────────
figure('Position', [50 50 900 700]);
scatter3(px, py, pz, 8, all_amp, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
clim([amp_thresh, 0]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('CFAR Point Cloud — All Beams, Top %d per beam (> %.0f dB)', ...
      top_n, amp_thresh));
axis equal; grid on;

%% ── FIGURE 2: TOP-DOWN VIEW ──────────────────────────────────────────────
figure('Position', [50 50 800 600]);
scatter(px, py, 6, all_amp, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
clim([amp_thresh, 0]);
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('CFAR Point Cloud — Top-down view (> %.0f dB)', amp_thresh));
axis equal; grid on;

%% ── FIGURE 3: ATTENUATION MAP (EL vs AZ) ────────────────────────────────
figure('Position', [50 50 900 600]);
scatter(fit_az, fit_el, 30, fit_rate, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'Attenuation (dB/m)';
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title(sprintf('Attenuation Rate Map — %d Valid Beams (R^2 \\geq %.2f)', ...
      sum(valid_fit), min_rsq));
grid on;

%% ── FIGURE 4: ATTENUATION vs ELEVATION ──────────────────────────────────
[el_sorted, el_order] = sort(fit_el);
rate_el_sorted = fit_rate(el_order);
rsq_el_sorted  = fit_rsq(el_order);

figure('Position', [50 50 900 500]);
hold on;
scatter(el_sorted, rate_el_sorted, 25, rsq_el_sorted, 'filled', ...
        'DisplayName', 'Beam attenuation');
colormap jet;
cb = colorbar; cb.Label.String = 'R^2';
clim([0.1, 0.9]);
plot(el_sorted, movmedian(rate_el_sorted, 20), 'r-', 'LineWidth', 2.0, ...
     'DisplayName', 'Running median (n=20)');
yline(mean(fit_rate), 'k--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('Mean: %.2f dB/m', mean(fit_rate)));
hold off;
grid on;
xlabel('Elevation (deg)');
ylabel('Attenuation rate (dB/m)');
title('Attenuation Rate vs Elevation Angle — All Valid Beams');
legend('Location', 'northeast');

%% ── FIGURE 5: ATTENUATION vs AZIMUTH ────────────────────────────────────
[az_sorted, az_order] = sort(fit_az);
rate_az_sorted = fit_rate(az_order);
rsq_az_sorted  = fit_rsq(az_order);

figure('Position', [50 50 900 500]);
hold on;
scatter(az_sorted, rate_az_sorted, 25, rsq_az_sorted, 'filled', ...
        'DisplayName', 'Beam attenuation');
colormap jet;
cb = colorbar; cb.Label.String = 'R^2';
clim([0.1, 0.9]);
plot(az_sorted, movmedian(rate_az_sorted, 20), 'r-', 'LineWidth', 2.0, ...
     'DisplayName', 'Running median (n=20)');
yline(mean(fit_rate), 'k--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('Mean: %.2f dB/m', mean(fit_rate)));
hold off;
grid on;
xlabel('Azimuth (deg)');
ylabel('Attenuation rate (dB/m)');
title('Attenuation Rate vs Azimuth Angle — All Valid Beams');
legend('Location', 'northeast');

%% ── FIGURE 6: HISTOGRAM ──────────────────────────────────────────────────
figure('Position', [50 50 700 450]);
histogram(fit_rate, 30, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on;
xline(mean(fit_rate),   'r-',  sprintf('Mean %.2f dB/m',   mean(fit_rate)),   'LineWidth', 2);
xline(median(fit_rate), 'g--', sprintf('Median %.2f dB/m', median(fit_rate)), 'LineWidth', 2);
hold off;
grid on;
xlabel('Attenuation Rate (dB/m)');
ylabel('Number of beams');
title('Distribution of Canopy Attenuation Rates — All Valid Beams');
legend('Beams', 'Mean', 'Median', 'Location', 'northwest');

%% ── FIGURE 7: ALL DECAY CURVES ───────────────────────────────────────────
figure('Position', [50 50 900 600]);
hold on;

cmap_lines = jet(sum(valid_fit));
valid_rows = find(valid_fit);
n_plotted  = 0;

for k = 1:length(valid_rows)
    row = valid_rows(k);
    ri  = valid_idx(row);

    fft_full     = double(fft_data{ri});
    fft_cut_k    = fft_full(cut_mask);
    fft_smooth_k = movmean(fft_cut_k, smooth_win);

    canopy_vals_k = fft_smooth_k(in_canopy);
    canopy_x_k    = range_cut(in_canopy);

    if isempty(canopy_vals_k), continue; end

    % Normalise to peak
    [pk_val, pk_i] = max(canopy_vals_k);
    if (length(canopy_vals_k) - pk_i) < 4, continue; end

    x_plot = canopy_x_k(pk_i:end);
    y_plot = canopy_vals_k(pk_i:end) - pk_val;

    % Only plot if decay is within display range
    if max(y_plot) < display_floor, continue; end

    plot(x_plot, y_plot, 'Color', [cmap_lines(k,:), 0.15], ...
         'LineWidth', 0.5, 'HandleVisibility', 'off');
    n_plotted = n_plotted + 1;
end

fprintf('Decay curves plotted: %d\n', n_plotted);

all_x    = linspace(canopy_start_physical, canopy_end_physical, 100);
mean_fit = mean(fit_rate) * (all_x - all_x(1));
plot(all_x, mean_fit, 'k-', 'LineWidth', 2.5, ...
     'DisplayName', sprintf('Mean fit: %.2f dB/m', mean(fit_rate)));
xline(canopy_end_physical, 'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
hold off;
grid on;
legend('Location', 'southwest');
xlabel('Range (m)');
ylabel('Relative amplitude (dB, normalised to peak)');
title(sprintf('All Valid Beam Decay Curves (n=%d) — Normalised to Peak', n_plotted));
xlim([canopy_start_physical, cut_hi]);
ylim([-100, 5]);

%% ── FIGURE 8: REFERENCE BEAM ─────────────────────────────────────────────
% Two canopy boundaries are used here:
%   canopy_start_physical = 12.7 m — physical canopy front from field measurement
%                                    used for CFAR window and all other experiments
%   signal_onset          = 17.75 m — range where FFT rises above noise for this beam
%                                     used only for the attenuation fit

target_el_rad = -1.385442;
target_az_rad =  0.536165;

dist_beam = abs(raw_el - target_el_rad) + 0.1 * abs(raw_az - target_az_rad);
dist_beam(fft_lengths ~= n_fft_bins) = inf;
[~, beam_row] = min(dist_beam);

fft_beam     = double(fft_data{beam_row});
fft_cut_b    = fft_beam(cut_mask);
fft_smooth_b = movmean(fft_cut_b, smooth_win);

% Clip to display floor
fft_cut_b_display = max(fft_cut_b, display_floor);

% CFAR on reference beam
fft_lin_b  = 10.^(fft_cut_b / 10);
cfar_det_b = false(1, n_cut);
cfar_thr_b = nan(1, n_cut);

for k = 1:n_cut
    lo    = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
    hi    = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
    train = [lo, hi];
    if isempty(train), continue; end
    noise_est     = mean(fft_lin_b(train));
    cfar_thr_b(k) = 10*log10(alpha * noise_est);
    if fft_cut_b(k) > cfar_thr_b(k)
        cfar_det_b(k) = true;
    end
end

det_bins_b = find(cfar_det_b);
[~, ord_b] = sort(fft_cut_b(det_bins_b), 'descend');
top_b      = det_bins_b(ord_b(1:min(top_n, end)));

% Beam-specific attenuation window
signal_onset  = 17.75;
in_canopy_b   = range_cut >= signal_onset & range_cut <= canopy_end_physical;
canopy_vals_b = fft_smooth_b(in_canopy_b);
canopy_x_b    = range_cut(in_canopy_b);

% Find peak within canopy window
[~, pk_b]  = max(canopy_vals_b);
peak_range = canopy_x_b(pk_b);
peak_amp   = canopy_vals_b(pk_b);

% Fit slope across full canopy window
x_decay_b = canopy_x_b;
y_decay_b = canopy_vals_b;
p_b       = polyfit(x_decay_b, y_decay_b, 1);

% Anchor fit line to peak amplitude at peak range
fit_b = p_b(1) * (x_decay_b - peak_range) + peak_amp;

rsq_b = 1 - sum((y_decay_b - fit_b).^2) / ...
            sum((y_decay_b - mean(y_decay_b)).^2);

% Clip to display range
x_fit_clip = x_decay_b(fit_b >= -95);
fit_b_clip = fit_b(fit_b >= -95);

fprintf('\nReference beam: El = %.2f deg, Az = %.2f deg\n', ...
        rad2deg(raw_el(beam_row)), rad2deg(raw_az(beam_row)));
fprintf('Physical canopy front: %.2f m\n', canopy_start_physical);
fprintf('Signal onset:          %.2f m\n', signal_onset);
fprintf('Peak at:               %.2f m,  amplitude: %.1f dB\n', peak_range, peak_amp);
fprintf('Attenuation:           %.3f dB/m,  R^2 = %.3f\n', p_b(1), rsq_b);

figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut_b_display, 1, 'FaceColor', [0.75 0.85 1.0], ...
        'EdgeColor', 'none', 'DisplayName', 'Raw FFT');
b.BaseValue = display_floor;
plot(range_cut, max(cfar_thr_b, display_floor), 'r-', 'LineWidth', 1.2, ...
     'DisplayName', 'CFAR threshold');
plot(range_cut, max(fft_smooth_b, display_floor), 'b-', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Smoothed (%d-bin MA)', smooth_win));
plot(x_fit_clip, fit_b_clip, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Attenuation fit: %.2f dB/m, R^2=%.2f', p_b(1), rsq_b));
scatter(range_cut(top_b), max(fft_cut_b(top_b), display_floor), 40, 'g', 'filled', ...
        'DisplayName', sprintf('Top %d detections', length(top_b)));
xline(peak_range,            'm:',  'Peak', ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, ...
      'HandleVisibility', 'off');
xline(canopy_start_physical, 'k--', 'Canopy front (field)', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(signal_onset,          'b--', 'Signal onsetCanopy front (field)', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(canopy_end_physical,   'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)');
ylabel('FFT Amplitude (dB)');
title(sprintf('Reference Beam — El = %.1f^\\circ, Az = %.1f^\\circ', ...
      rad2deg(raw_el(beam_row)), rad2deg(raw_az(beam_row))));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([-95, 0]);

%% ── FIGURE 9: MEAN ATTENUATION PROFILE ACROSS ALL BEAMS ─────────────────
% Bin all CFAR detection points by range and compute mean amplitude per bin
bin_width   = 0.5;
edges       = cut_lo:bin_width:cut_hi;
bin_ctrs    = edges(1:end-1) + bin_width/2;
n_bins_att  = length(bin_ctrs);

mean_amp_all = nan(1, n_bins_att);
std_amp_all  = nan(1, n_bins_att);

for i = 1:n_bins_att
    in_bin = all_range >= edges(i) & all_range < edges(i+1);
    if sum(in_bin) >= 5
        mean_amp_all(i) = mean(all_amp(in_bin));
        std_amp_all(i)  = std(all_amp(in_bin));
    end
end

% Fit mean attenuation from signal onset to canopy end
fit_mask_all = bin_ctrs >= signal_onset & bin_ctrs <= canopy_end_physical ...
               & ~isnan(mean_amp_all);
p_all   = polyfit(bin_ctrs(fit_mask_all), mean_amp_all(fit_mask_all), 1);
fit_all = polyval(p_all, bin_ctrs(fit_mask_all));
rsq_all = 1 - sum((mean_amp_all(fit_mask_all) - fit_all).^2) / ...
              sum((mean_amp_all(fit_mask_all) - mean(mean_amp_all(fit_mask_all))).^2);

fprintf('\nMean attenuation across all beams: %.3f dB/m,  R^2 = %.3f\n', ...
        p_all(1), rsq_all);

figure('Position', [50 50 900 500]);
hold on;
fill([bin_ctrs(~isnan(mean_amp_all)), fliplr(bin_ctrs(~isnan(mean_amp_all)))], ...
     [mean_amp_all(~isnan(mean_amp_all)) + std_amp_all(~isnan(mean_amp_all)), ...
      fliplr(mean_amp_all(~isnan(mean_amp_all)) - std_amp_all(~isnan(mean_amp_all)))], ...
     [0.8 0.8 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'plus minus 1 sigma');
plot(bin_ctrs, mean_amp_all, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10, ...
     'DisplayName', 'Mean amplitude');
plot(bin_ctrs(fit_mask_all), fit_all, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Attenuation fit: %.2f dB/m, R^2=%.2f', p_all(1), rsq_all));
xline(canopy_start_physical, 'k--', 'Signal onset', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(signal_onset,          'b--', 'Canopy front', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(canopy_end_physical,   'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)');
ylabel('Mean FFT Amplitude (dB)');
title(sprintf('Mean Attenuation Profile — All CFAR Detections\n%.2f dB/m from signal onset to canopy back', p_all(1)));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([-95, 0]);

%% ── FIGURE 10: HORIZONTAL SLICE HEAT MAP ────────────────────────────────
grid_res_10 = 0.3;
mid_el_deg  = median(rad2deg(all_el));
el_tol_deg  = 2.0;

slice_10     = abs(rad2deg(all_el) - mid_el_deg) <= el_tol_deg;
px_10        = px(slice_10);
py_10        = py(slice_10);
amp_10       = all_amp(slice_10);

fprintf('\nHorizontal slice at El = %.1f deg (tol +/- %.1f deg): %d points\n', ...
        mid_el_deg, el_tol_deg, sum(slice_10));

x_edges_10  = min(px_10)-1 : grid_res_10 : max(px_10)+1;
y_edges_10  = min(py_10)-1 : grid_res_10 : max(py_10)+1;
n_x10       = length(x_edges_10) - 1;
n_y10       = length(y_edges_10) - 1;
x_ctrs_10   = x_edges_10(1:end-1) + grid_res_10/2;
y_ctrs_10   = y_edges_10(1:end-1) + grid_res_10/2;

hmap_10 = nan(n_y10, n_x10);
for ix = 1:n_x10
    for iy = 1:n_y10
        in_cell = px_10 >= x_edges_10(ix) & px_10 < x_edges_10(ix+1) & ...
                  py_10 >= y_edges_10(iy) & py_10 < y_edges_10(iy+1);
        if sum(in_cell) >= 2
            hmap_10(iy, ix) = mean(amp_10(in_cell));
        end
    end
end

figure('Position', [50 50 800 700]);
imagesc(x_ctrs_10, y_ctrs_10, hmap_10);
set(gca, 'YDir', 'normal');
colormap jet;
cb = colorbar; cb.Label.String = 'Mean FFT Amplitude (dB)';
clim([-90, -30]);
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Horizontal Slice Heat Map — El = %.1f^\\circ (\\pm%.1f^\\circ)', ...
      mid_el_deg, el_tol_deg));
axis equal tight; grid on;

%% ── FIGURE 11: VERTICAL SLICE COMPARISON — SINGLE RETURN vs FFT CFAR ────
grid_res_11 = 0.3;
mid_az_11   = median(rad2deg(all_az));
az_tol_11   = 2.0;

%% FFT CFAR slice
slice_fft_11  = abs(rad2deg(all_az) - mid_az_11) <= az_tol_11;
px_fft_11     = px(slice_fft_11);
pz_fft_11     = pz(slice_fft_11);
amp_fft_11    = all_amp(slice_fft_11);
px_fft_11     = px_fft_11(:);
pz_fft_11     = pz_fft_11(:);
amp_fft_11    = amp_fft_11(:);
fprintf('\nFFT CFAR vertical slice: %d points\n', length(px_fft_11));

%% Single strongest return per beam
% Load original data freshly to avoid any variable conflicts
data_s        = readmatrix('pointcloud_20260525-125755_0_fft_cfar.csv');
data_s        = data_s(~(data_s(:,3) == 0 | data_s(:,3) == 1000), :);
data_s        = data_s(data_s(:,3) >= -85, :);

az_s          = data_s(:,2);
el_s          = data_s(:,1);
r_s           = data_s(:,4);
amp_s         = data_s(:,3);
el_sph_s      = (pi/2) + el_s;
[px_s, py_s, pz_s] = sph2cart(az_s, el_sph_s, r_s);

az_s_deg      = rad2deg(az_s);
slice_s       = abs(az_s_deg - mid_az_11) <= az_tol_11;
px_s_sl       = px_s(slice_s);
pz_s_sl       = pz_s(slice_s);
amp_s_sl      = amp_s(slice_s);
el_s_sl       = el_s(slice_s);

% Keep strongest return per unique elevation
unique_el_s   = unique(el_s_sl);
n_uel         = length(unique_el_s);
px_sing_11    = zeros(n_uel, 1);
pz_sing_11    = zeros(n_uel, 1);
amp_sing_11   = zeros(n_uel, 1);

for i = 1:n_uel
    mask_i    = el_s_sl == unique_el_s(i);
    amps_i    = amp_s_sl(mask_i);
    [~, best] = max(amps_i);
    px_i      = px_s_sl(mask_i);
    pz_i      = pz_s_sl(mask_i);
    px_sing_11(i)  = px_i(best);
    pz_sing_11(i)  = pz_i(best);
    amp_sing_11(i) = amps_i(best);
end

fprintf('Single strongest return per beam: %d points\n', n_uel);

%% Shared grid
x_min_11   = min([px_fft_11; px_sing_11]) - 1;
x_max_11   = max([px_fft_11; px_sing_11]) + 1;
z_min_11   = min([pz_fft_11; pz_sing_11]) - 1;
z_max_11   = max([pz_fft_11; pz_sing_11]) + 1;

x_edges_11 = x_min_11 : grid_res_11 : x_max_11;
z_edges_11 = z_min_11 : grid_res_11 : z_max_11;
n_x11      = length(x_edges_11) - 1;
n_z11      = length(z_edges_11) - 1;
x_ctrs_11  = x_edges_11(1:end-1) + grid_res_11/2;
z_ctrs_11  = z_edges_11(1:end-1) + grid_res_11/2;

fprintf('Grid: %d x %d  |  x: %.1f to %.1f  |  z: %.1f to %.1f\n', ...
        n_x11, n_z11, x_min_11, x_max_11, z_min_11, z_max_11);

%% FFT CFAR heat map
hmap_fft_11 = nan(n_z11, n_x11);
for ix = 1:n_x11
    for iz = 1:n_z11
        in_cell = px_fft_11 >= x_edges_11(ix) & px_fft_11 < x_edges_11(ix+1) & ...
                  pz_fft_11 >= z_edges_11(iz) & pz_fft_11 < z_edges_11(iz+1);
        if sum(in_cell) >= 1
            hmap_fft_11(iz, ix) = mean(amp_fft_11(in_cell));
        end
    end
end

%% Single return heat map
hmap_sing_11 = nan(n_z11, n_x11);
for ix = 1:n_x11
    for iz = 1:n_z11
        in_cell = px_sing_11 >= x_edges_11(ix) & px_sing_11 < x_edges_11(ix+1) & ...
                  pz_sing_11 >= z_edges_11(iz) & pz_sing_11 < z_edges_11(iz+1);
        if sum(in_cell) >= 1
            hmap_sing_11(iz, ix) = mean(amp_sing_11(in_cell));
        end
    end
end

fprintf('FFT hmap non-nan:    %d cells\n', sum(~isnan(hmap_fft_11(:))));
fprintf('Single hmap non-nan: %d cells\n', sum(~isnan(hmap_sing_11(:))));

%% Plot
figure('Position', [50 50 1400 650]);

subplot(1,2,1);
imagesc(x_ctrs_11, z_ctrs_11, hmap_sing_11);
set(gca, 'YDir', 'normal');
colormap jet;
cb = colorbar; cb.Label.String = 'Amplitude (dB)';
clim([-90, -30]);
xlabel('X (m) — Range');
ylabel('Z (m) — Height');
title(sprintf('Single Strongest Return Per Beam\nAz = %.1f^\\circ (\\pm%.1f^\\circ),  %d points', ...
      mid_az_11, az_tol_11, n_uel));
xlim([10, 41]); ylim([-5, 22]);   % <-- zoomed
grid on;

subplot(1,2,2);
imagesc(x_ctrs_11, z_ctrs_11, hmap_fft_11);
set(gca, 'YDir', 'normal');
colormap jet;
cb = colorbar; cb.Label.String = 'Mean FFT Amplitude (dB)';
clim([-90, -30]);
xlabel('X (m) — Range');
ylabel('Z (m) — Height');
title(sprintf('FFT CFAR Top %d Per Beam\nAz = %.1f^\\circ (\\pm%.1f^\\circ),  %d points', ...
      top_n, mid_az_11, az_tol_11, length(px_fft_11)));
xlim([10, 41]); ylim([-5, 22]);   % <-- zoomed
 grid on;

sgtitle(sprintf('Vertical Slice — Single Return vs FFT CFAR  |  Az = %.1f^\\circ', mid_az_11));
%% ── SAVE ALL FIGURES ─────────────────────────────────────────────────────
save_dir = fullfile(pwd, 'figures');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Saving to: %s\n', save_dir);

fig_names = {
    'fig01_3d_point_cloud',
    'fig02_topdown_view',
    'fig03_attenuation_map',
    'fig04_attenuation_vs_elevation',
    'fig05_attenuation_vs_azimuth',
    'fig06_attenuation_histogram',
    'fig07_all_decay_curves',
    'fig08_reference_beam',
    'fig09_mean_attenuation_profile',
    'fig10_horizontal_slice_heatmap',
    'fig11_vertical_slice_comparison'
};

all_figs = findall(0, 'Type', 'figure');
all_figs = sort(all_figs);   % sort by creation order

for k = 1:min(length(all_figs), length(fig_names))
    fname = fullfile(save_dir, fig_names{k});
    exportgraphics(all_figs(k), [fname '.png'], 'Resolution', 200);
    fprintf('Saved: %s.png\n', fname);
end
fprintf('Done. %d figures saved to %s/\n', length(fig_names), save_dir);