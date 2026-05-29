%% Single Beam CFAR Analysis — Step by Step
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

% Find the target beam
target_el_rad = -1.385442;
target_az_rad =  0.536165;
bin_range     = 0.125;
canopy_start  = 12.7;
canopy_end    = 35.2;

dist_beam  = abs(raw_el - target_el_rad) + 0.1 * abs(raw_az - target_az_rad);
dist_beam(fft_lengths ~= mode(fft_lengths)) = inf;
[~, row_idx] = min(dist_beam);

fft_full   = double(fft_data{row_idx});
n_bins     = length(fft_full);
range_full = (0:n_bins-1) * bin_range + (2 * bin_range);

fprintf('Beam: El = %.2f deg, Az = %.2f deg\n', ...
        rad2deg(raw_el(row_idx)), rad2deg(raw_az(row_idx)));

%% ── STEP 1: CUT TO CANOPY WINDOW (10–45 m) ──────────────────────────────
cut_lo   = 10.0;
cut_hi   = 45.0;
cut_mask = range_full >= cut_lo & range_full <= cut_hi;

range_cut = range_full(cut_mask);
fft_cut   = fft_full(cut_mask);
n_cut     = length(range_cut);

figure('Position', [50 50 900 500]);
b = bar(range_cut, fft_cut, 1, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
b.BaseValue = min(fft_cut) - 5;
hold on;
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title('Step 1 — Raw FFT cut to canopy window (10–45 m)');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut) - 5, 0]);

%% ── STEP 2: RUN CA-CFAR ──────────────────────────────────────────────────
cfar_guard = 1;
cfar_train = 4;
cfar_pfa   = 1e-1;
alpha      = cfar_train * (cfar_pfa^(-1/cfar_train) - 1);

fft_linear  = 10.^(fft_cut / 10);
cfar_thresh_db  = nan(1, n_cut);
cfar_det        = false(1, n_cut);

for k = 1:n_cut
    lo = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
    hi = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
    train = [lo, hi];
    if isempty(train), continue; end
    
    noise_est         = mean(fft_linear(train));
    thresh_lin        = alpha * noise_est;
    cfar_thresh_db(k) = 10*log10(thresh_lin);
    
    % Fix: compare in dB directly rather than linear
    if fft_cut(k) > cfar_thresh_db(k)
        cfar_det(k) = true;
    end
end

fprintf('CFAR detections: %d / %d bins\n', sum(cfar_det), n_cut);

figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut, 1, 'FaceColor', [0.75 0.85 1.0], 'EdgeColor', 'none', ...
        'DisplayName', 'Raw FFT');
b.BaseValue = min(fft_cut) - 5;
plot(range_cut, cfar_thresh_db, 'r-', 'LineWidth', 1.5, 'DisplayName', 'CFAR threshold');
scatter(range_cut(cfar_det), fft_cut(cfar_det), 40, 'g', 'filled', ...
        'DisplayName', 'Detections');
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title('Step 2 — CA-CFAR applied to canopy window');
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut) - 5, 0]);

%% ── STEP 3: PICK TOP 60 DETECTIONS ──────────────────────────────────────
det_bins = find(cfar_det);
det_amps = fft_cut(det_bins);
[~, ord] = sort(det_amps, 'descend');
top_n    = min(67, length(det_bins));
top_bins = det_bins(ord(1:top_n));
top_amps = fft_cut(top_bins);
top_ranges = range_cut(top_bins);

fprintf('Top %d detections: range %.2f–%.2f m, amp %.1f–%.1f dB\n', ...
        top_n, min(top_ranges), max(top_ranges), min(top_amps), max(top_amps));

figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut, 1, 'FaceColor', [0.75 0.85 1.0], 'EdgeColor', 'none', ...
        'DisplayName', 'Raw FFT');
b.BaseValue = min(fft_cut) - 5;
plot(range_cut, cfar_thresh_db, 'r-', 'LineWidth', 1.5, 'DisplayName', 'CFAR threshold');
scatter(top_ranges, top_amps, 60, 'r', 'filled', 'DisplayName', ...
        sprintf('Top %d detections', top_n));
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title(sprintf('Step 3 — Top %d CFAR detections highlighted', top_n));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut) - 5, 0]);

%% ── STEP 4: BUILD POINT CLOUD FROM TOP DETECTIONS ───────────────────────
beam_el  = raw_el(row_idx);
beam_az  = raw_az(row_idx);
el_sph   = (pi/2) + beam_el;

[px, py, pz] = sph2cart(repmat(beam_az, 1, top_n), ...
                         repmat(el_sph,  1, top_n), ...
                         top_ranges);

figure('Position', [50 50 800 600]);
scatter3(px, py, pz, 60, top_amps, 'filled');
colormap jet; 
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Step 4 — Point cloud from top %d CFAR detections\nEl = %.1f°, Az = %.1f°', ...
      top_n, rad2deg(beam_el), rad2deg(beam_az)));
grid on; axis equal;

%% ── STEPS 5 & 6: SINGLE BEAM ATTENUATION THROUGH CANOPY ─────────────────
% Smooth the raw FFT and fit linear attenuation inside canopy only
smooth_win  = 8;
fft_smooth  = movmean(fft_cut, smooth_win);
in_canopy   = range_cut >= canopy_start & range_cut <= canopy_end;

p      = polyfit(range_cut(in_canopy), fft_smooth(in_canopy), 1);
fit_y  = polyval(p, range_cut(in_canopy));
r_sq   = 1 - sum((fft_smooth(in_canopy) - fit_y).^2) / ...
             sum((fft_smooth(in_canopy) - mean(fft_smooth(in_canopy))).^2);

fprintf('\nSingle beam attenuation: %.3f dB/m,  R^2 = %.3f\n', p(1), r_sq);

figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut, 1, 'FaceColor', [0.75 0.85 1.0], 'EdgeColor', 'none', ...
        'DisplayName', 'Raw FFT');
b.BaseValue = min(fft_cut) - 5;
plot(range_cut, fft_smooth, 'b-', 'LineWidth', 1.8, 'DisplayName', ...
     sprintf('Smoothed (%d-bin MA)', smooth_win));
plot(range_cut(in_canopy), fit_y, 'r--', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Linear fit: %.2f dB/m, R^2=%.2f', p(1), r_sq));
scatter(top_ranges, top_amps, 40, 'r', 'filled', ...
        'DisplayName', sprintf('Top %d CFAR detections', top_n));
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title(sprintf('Steps 5 & 6 — Single beam attenuation\nEl = %.1f°, Az = %.1f°', ...
      rad2deg(beam_el), rad2deg(beam_az)));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut) - 5, 0]);


%% ── COMBINED: CFAR + POINT CLOUD + ATTENUATION ───────────────────────────

%% CFAR
cfar_guard = 1;
cfar_train = 4;
cfar_pfa   = 0.3;
alpha      = cfar_train * (cfar_pfa^(-1/cfar_train) - 1);

fft_linear     = 10.^(fft_cut / 10);
cfar_thresh_db = nan(1, n_cut);
cfar_det       = false(1, n_cut);

for k = 1:n_cut
    lo = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
    hi = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
    train = [lo, hi];
    if isempty(train), continue; end
    noise_est         = mean(fft_linear(train));
    thresh_lin        = alpha * noise_est;
    cfar_thresh_db(k) = 10*log10(thresh_lin);
    if fft_cut(k) > cfar_thresh_db(k)
        cfar_det(k) = true;
    end
end

%% Take top N — either CFAR detections or raw peaks
use_cfar = true;   % set false to skip CFAR and just take top N by amplitude

if use_cfar
    det_bins = find(cfar_det);
else
    [~, det_bins] = sort(fft_cut, 'descend');
end

top_n      = min(67, length(det_bins));
det_amps   = fft_cut(det_bins);
[~, ord]   = sort(det_amps, 'descend');
top_bins   = det_bins(ord(1:top_n));
top_amps   = fft_cut(top_bins);
top_ranges = range_cut(top_bins);

fprintf('Using %d detections\n', top_n);

%% Point cloud from top detections
beam_el = raw_el(row_idx);
beam_az = raw_az(row_idx);
el_sph  = (pi/2) + beam_el;

[px, py, pz] = sph2cart(repmat(beam_az, 1, top_n), ...
                         repmat(el_sph,  1, top_n), ...
                         top_ranges);

figure('Position', [50 50 800 600]);
scatter3(px, py, pz, 60, top_amps, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Point Cloud — Top %d CFAR Detections\nEl = %.1f°, Az = %.1f°', ...
      top_n, rad2deg(beam_el), rad2deg(beam_az)));
grid on; axis equal;

%% Attenuation — peak to canopy end fit
smooth_win  = 8;
fft_smooth  = movmean(fft_cut, smooth_win);
in_canopy   = range_cut >= canopy_start & range_cut <= canopy_end;
canopy_vals = fft_smooth(in_canopy);
canopy_x    = range_cut(in_canopy);

[~, pk_idx] = max(canopy_vals);
fprintf('Canopy peak at: %.2f m,  amplitude: %.1f dB\n', ...
        canopy_x(pk_idx), canopy_vals(pk_idx));

if (length(canopy_vals) - pk_idx) >= 4
    x_decay = canopy_x(pk_idx:end);
    y_decay = canopy_vals(pk_idx:end);

    p     = polyfit(x_decay, y_decay, 1);
    fit_y = polyval(p, x_decay);
    r_sq  = 1 - sum((y_decay - fit_y).^2) / ...
                sum((y_decay - mean(y_decay)).^2);

    fprintf('Peak-to-end attenuation: %.3f dB/m,  R^2 = %.3f\n', p(1), r_sq);
else
    warning('Not enough points after peak for fit');
    x_decay = []; fit_y = []; p = [0 0]; r_sq = 0;
end

%% Combined attenuation figure
figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut, 1, 'FaceColor', [0.75 0.85 1.0], ...
        'EdgeColor', 'none', 'DisplayName', 'Raw FFT');
b.BaseValue = min(fft_cut) - 5;
plot(range_cut, cfar_thresh_db, 'r-', 'LineWidth', 1.2, ...
     'DisplayName', 'CFAR threshold');
plot(range_cut, fft_smooth, 'b-', 'LineWidth', 1.8, ...
     'DisplayName', sprintf('Smoothed (%d-bin MA)', smooth_win));
if (length(canopy_vals) - pk_idx) >= 4
    plot(x_decay, fit_y, 'k--', 'LineWidth', 2.0, ...
         'DisplayName', sprintf('Peak-to-end fit: %.2f dB/m, R^2=%.2f', p(1), r_sq));
    xline(canopy_x(pk_idx), 'm:', 'Peak', 'LineWidth', 1.5, ...
          'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
scatter(top_ranges, top_amps, 40, 'g', 'filled', ...
        'DisplayName', sprintf('Top %d detections', top_n));
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title(sprintf('Single Beam — CFAR Detections and Attenuation Fit\nEl = %.1f°, Az = %.1f°', ...
      rad2deg(beam_el), rad2deg(beam_az)));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut) - 5, 0]);

%% Attenuation vs range — detections only
figure('Position', [50 50 900 400]);
hold on;
scatter(top_ranges, top_amps, 40, top_amps, 'filled', 'DisplayName', 'CFAR detections');
colormap jet; cb = colorbar; cb.Label.String = 'Amplitude (dB)';
if (length(canopy_vals) - pk_idx) >= 4
    plot(x_decay, fit_y, 'k--', 'LineWidth', 2.0, ...
         'DisplayName', sprintf('Peak-to-end fit: %.2f dB/m', p(1)));
end
xline(canopy_start, 'k--', 'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title('Attenuation Profile — CFAR Detection Points Only');
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);