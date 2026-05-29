%% All Beams CFAR Point Cloud and Attenuation Analysis
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

fprintf('Loading %d rows...\n', n_rows);
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
bin_range    = 0.125;
canopy_start = 12.7;
canopy_end   = 35.2;
cut_lo       = 10.0;
cut_hi       = 45.0;
top_n        = 60;

% CFAR parameters
cfar_guard = 2;
cfar_train = 8;
cfar_pfa   = 1e-2;
alpha      = cfar_train * (cfar_pfa^(-1/cfar_train) - 1);
fprintf('CFAR alpha: %.3f\n', alpha);

% Build range axes
n_fft_bins = mode(fft_lengths);
range_full = (0:n_fft_bins-1) * bin_range + (2 * bin_range);
cut_mask   = range_full >= cut_lo & range_full <= cut_hi;
range_cut  = range_full(cut_mask);
n_cut      = length(range_cut);

% Valid rows only
valid_mask = fft_lengths == n_fft_bins;
valid_idx  = find(valid_mask);
n_valid    = length(valid_idx);
fprintf('Valid FFT rows: %d\n', n_valid);

%% ── STEPS 1-5: CFAR ON ALL BEAMS → POINT CLOUD ──────────────────────────
% Pre-allocate output arrays (generous upper bound)
max_pts   = n_valid * top_n;
all_range = zeros(1, max_pts);
all_amp   = zeros(1, max_pts);
all_el    = zeros(1, max_pts);
all_az    = zeros(1, max_pts);
ptr       = 0;

fprintf('Running CFAR on %d beams...\n', n_valid);
for row = 1:n_valid

    if mod(row, 500) == 0
        fprintf('  Processing row %d / %d\n', row, n_valid);
    end

    ri = valid_idx(row);

    %% Step 1: cut to canopy window
    fft_full = double(fft_data{ri});
    fft_cut  = fft_full(cut_mask);

    %% Step 2: CA-CFAR in linear domain
    fft_lin     = 10.^(fft_cut / 10);
    cfar_det    = false(1, n_cut);

    for k = 1:n_cut
        lo    = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
        hi    = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
        train = [lo, hi];
        if isempty(train), continue; end
        if fft_lin(k) > alpha * mean(fft_lin(train))
            cfar_det(k) = true;
        end
    end

    %% Step 3: pick top N by amplitude
    det_bins = find(cfar_det);
    if isempty(det_bins), continue; end

    det_amps = fft_cut(det_bins);
    [~, ord] = sort(det_amps, 'descend');
    keep     = det_bins(ord(1:min(top_n, end)));
    n_keep   = length(keep);

    %% Step 4: store detections
    idx_range        = ptr+1 : ptr+n_keep;
    all_range(idx_range) = range_cut(keep);
    all_amp(idx_range)   = fft_cut(keep);
    all_el(idx_range)    = raw_el(ri);
    all_az(idx_range)    = raw_az(ri);
    ptr = ptr + n_keep;
end

% Filter out sidelobe noise below -90 dB
amp_threshold = -80;
keep_mask = all_amp >= amp_threshold;

all_range = all_range(keep_mask);
all_amp   = all_amp(keep_mask);
all_el    = all_el(keep_mask);
all_az    = all_az(keep_mask);

fprintf('Points after -90 dB filter: %d (removed %d)\n', ...
        sum(keep_mask), sum(~keep_mask));

% Recompute Cartesian
el_sph       = (pi/2) + all_el;
[px, py, pz] = sph2cart(all_az, el_sph, all_range);

% Trim to actual size
all_range = all_range(1:ptr);
all_amp   = all_amp(1:ptr);
all_el    = all_el(1:ptr);
all_az    = all_az(1:ptr);
fprintf('Total detections across all beams: %d\n', ptr);

%% Step 5: convert to Cartesian point cloud
el_sph       = (pi/2) + all_el;
[px, py, pz] = sph2cart(all_az, el_sph, all_range);

% ── Figure 1: 3D point cloud coloured by amplitude
figure('Position', [50 50 900 700]);
scatter3(px, py, pz, 8, all_amp, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('CFAR Point Cloud — All Beams, Top %d per beam', top_n));
axis equal; grid on;

% ── Figure 2: top-down view (X-Y)
figure('Position', [50 50 800 600]);
scatter(px, py, 6, all_amp, 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'FFT Amplitude (dB)';
xlabel('X (m)'); ylabel('Y (m)');
title('CFAR Point Cloud — Top-down view');
axis equal; grid on;

% ── Figure 3: range vs amplitude coloured by elevation
figure('Position', [50 50 900 500]);
scatter(all_range, all_amp, 6, rad2deg(all_el), 'filled');
colormap jet;
cb = colorbar; cb.Label.String = 'Elevation (deg)';
xlabel('Range (m)'); ylabel('Amplitude (dB)');
title('All CFAR Detections — Amplitude vs Range');
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom');
grid on;

%% ── STEPS 6 & 7: SINGLE BEAM ATTENUATION ────────────────────────────────
target_el_rad = -1.385442;
target_az_rad =  0.536165;

dist_beam = abs(raw_el - target_el_rad) + 0.1 * abs(raw_az - target_az_rad);
dist_beam(fft_lengths ~= n_fft_bins) = inf;
[~, beam_row] = min(dist_beam);

fft_beam  = double(fft_data{beam_row});
fft_cut_b = fft_beam(cut_mask);

% Re-run CFAR on this beam
fft_lin_b      = 10.^(fft_cut_b / 10);
cfar_det_b     = false(1, n_cut);
cfar_thresh_b  = nan(1, n_cut);

for k = 1:n_cut
    lo    = max(1,     k - cfar_guard - cfar_train) : max(1,    k - cfar_guard - 1);
    hi    = min(n_cut, k + cfar_guard + 1)          : min(n_cut, k + cfar_guard + cfar_train);
    train = [lo, hi];
    if isempty(train), continue; end
    noise_est        = mean(fft_lin_b(train));
    cfar_thresh_b(k) = 10*log10(alpha * noise_est);
    if fft_lin_b(k) > alpha * noise_est
        cfar_det_b(k) = true;
    end
end

% Top detections on this beam
det_bins_b = find(cfar_det_b);
[~, ord_b] = sort(fft_cut_b(det_bins_b), 'descend');
top_b      = det_bins_b(ord_b(1:min(top_n, end)));

% Smooth and fit attenuation
smooth_win = 8;
fft_smooth = movmean(fft_cut_b, smooth_win);
in_canopy  = range_cut >= canopy_start & range_cut <= canopy_end;
p          = polyfit(range_cut(in_canopy), fft_smooth(in_canopy), 1);
fit_y      = polyval(p, range_cut(in_canopy));
r_sq       = 1 - sum((fft_smooth(in_canopy) - fit_y).^2) / ...
                 sum((fft_smooth(in_canopy) - mean(fft_smooth(in_canopy))).^2);

fprintf('\nSingle beam — El = %.2f deg, Az = %.2f deg\n', ...
        rad2deg(raw_el(beam_row)), rad2deg(raw_az(beam_row)));
fprintf('Attenuation: %.3f dB/m,  R^2 = %.3f\n', p(1), r_sq);

% ── Figure 4: single beam attenuation
figure('Position', [50 50 900 500]);
hold on;
b = bar(range_cut, fft_cut_b, 1, 'FaceColor', [0.75 0.85 1.0], ...
        'EdgeColor', 'none', 'DisplayName', 'Raw FFT');
b.BaseValue = min(fft_cut_b) - 5;
plot(range_cut, cfar_thresh_b, 'r-',  'LineWidth', 1.5, ...
     'DisplayName', 'CFAR threshold');
plot(range_cut, fft_smooth,    'b-',  'LineWidth', 1.8, ...
     'DisplayName', sprintf('Smoothed (%d-bin MA)', smooth_win));
plot(range_cut(in_canopy), fit_y, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', sprintf('Fit: %.2f dB/m, R^2=%.2f', p(1), r_sq));
scatter(range_cut(top_b), fft_cut_b(top_b), 50, 'r', 'filled', ...
        'DisplayName', sprintf('Top %d detections', min(top_n, length(top_b))));
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)'); ylabel('FFT Amplitude (dB)');
title(sprintf('Single Beam Attenuation — El = %.1f°, Az = %.1f°', ...
      rad2deg(raw_el(beam_row)), rad2deg(raw_az(beam_row))));
legend('Location', 'northeast');
xlim([cut_lo, cut_hi]);
ylim([min(fft_cut_b) - 5, 0]);