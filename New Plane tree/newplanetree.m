% Clear workspace and close all figures
clear; close all;

data = readmatrix('pointcloud_20260525-125755_0_fft_cfar.csv');

%Data Structure:   
% (1) Elevation (Rad)
% (2) Azimuth (Rad)
% (3) Return Strength (dB)
% (4) Range (m)
% (5) FFT values (semicolon separated)
% (6) CFAR values (semicolon separated)

% Remove invalid peak values
data = data(~(data(:,3) == 0 | data(:,3) == 1000), :);

% Filter radar returns below -80 dB
data = data(data(:,3) >= -80, :);

az = data(:, 2);
el = (pi/2) + data(:, 1);
r  = data(:, 4);
returned_power = data(:, 3);

[x, y, z] = sph2cart(az, el, r);

%% Read raw file to extract FFT (col 5) and CFAR (col 6)
fid = fopen('pointcloud_20260525-125755_0_fft_cfar.csv', 'r');
raw_lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
raw_lines = raw_lines{1};

n_rows = length(raw_lines);
fft_data     = cell(n_rows, 1);
cfar_data    = cell(n_rows, 1);
fft_lengths  = zeros(n_rows, 1);
cfar_lengths = zeros(n_rows, 1);

for i = 1:n_rows
    cols = strsplit(raw_lines{i}, ',');
    if length(cols) >= 5
        fft_vals = strsplit(cols{5}, ';');
        fft_vals = str2double(fft_vals);
        fft_vals = fft_vals(~isnan(fft_vals));
        fft_data{i} = fft_vals;
        fft_lengths(i) = length(fft_vals);
    end
    if length(cols) >= 6
        cfar_vals = strsplit(cols{6}, ';');
        cfar_vals = str2double(cfar_vals);
        cfar_vals = cfar_vals(~isnan(cfar_vals));
        cfar_data{i} = cfar_vals;
        cfar_lengths(i) = length(cfar_vals);
    end
end

fprintf('--- FFT ---\n');
fprintf('Min length: %d\n', min(fft_lengths));
fprintf('Max length: %d\n', max(fft_lengths));
fprintf('Most common length: %d\n', mode(fft_lengths));
fprintf('--- CFAR ---\n');
fprintf('Min length: %d\n', min(cfar_lengths));
fprintf('Max length: %d\n', max(cfar_lengths));
fprintf('Most common length: %d\n', mode(cfar_lengths));

%% Canopy boundaries
canopy_start = 12.7;
canopy_end   = 35.2;

%% Figure 1: 3D point cloud coloured by return power
figure;
scatter3(x, y, z, 10, returned_power, 'filled');
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Plane Tree Radar Point Cloud - Return Power');
axis equal;

%% Figure 2: 3D point cloud coloured by height
figure;
scatter3(x, y, z, 10, z, 'filled');
clim([-3, 20]);
colormap jet; colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Plane Tree Radar Point Cloud - Height');
axis equal;

%% Figure 3: Attenuation vs Range (bin averaged peak amplitude)

bin_width   = 0.5;
edges       = 0:bin_width:50;
bin_centres = edges(1:end-1) + bin_width/2;
n_bins      = length(bin_centres);

mean_amp = nan(1, n_bins);
std_amp  = nan(1, n_bins);

r_all              = data(:, 4);
returned_power_all = data(:, 3);

for i = 1:n_bins
    in_bin = r_all >= edges(i) & r_all < edges(i+1);
    if sum(in_bin) >= 5
        mean_amp(i) = mean(returned_power_all(in_bin));
        std_amp(i)  = std(returned_power_all(in_bin));
    end
end

fit_mask  = bin_centres >= canopy_start & bin_centres <= canopy_end & ~isnan(mean_amp);
p         = polyfit(bin_centres(fit_mask), mean_amp(fit_mask), 1);
fit_line  = polyval(p, bin_centres(fit_mask));
r_squared = 1 - sum((mean_amp(fit_mask) - fit_line).^2) / ...
                sum((mean_amp(fit_mask) - mean(mean_amp(fit_mask))).^2);

fprintf('\nAttenuation rate: %.3f dB/m, R^2 = %.2f\n', p(1), r_squared);

figure;
hold on;
fill([bin_centres(fit_mask), fliplr(bin_centres(fit_mask))], ...
     [mean_amp(fit_mask) + std_amp(fit_mask), ...
      fliplr(mean_amp(fit_mask) - std_amp(fit_mask))], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(bin_centres, mean_amp, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 10);
plot(bin_centres(fit_mask), fit_line, 'r--', 'LineWidth', 1.5);
xline(canopy_start, 'k--', 'Canopy front', 'LabelVerticalAlignment', 'bottom');
xline(canopy_end,   'k--', 'Canopy back',  'LabelVerticalAlignment', 'bottom');
hold off;
grid on;
xlabel('Range (m)');
ylabel('Mean Amplitude (dB)');
title('Attenuation vs Range through Plane Tree Canopy');
legend('±1\sigma', 'Mean amplitude', ...
       sprintf('Linear fit: %.3f dB/m, R^2=%.2f', p(1), r_squared), ...
       'Location', 'northeast');
xlim([5, 50]);

%% Figure 4: CFAR threshold vs FFT amplitude vs range

bin_range  = 0.125;
n_bins_fft = mode(fft_lengths);
range_axis = (0:n_bins_fft-1) * bin_range + (2 * bin_range);

% Build FFT matrix and normalise to dB relative to max
valid_fft   = fft_lengths == n_bins_fft;
fft_matrix  = cell2mat(fft_data(valid_fft));
fft_linear  = abs(fft_matrix);
fft_norm    = fft_linear ./ max(fft_linear(:));
fft_db      = 20*log10(fft_norm);

% Build CFAR matrix - already in dB
cfar_mode_len = mode(cfar_lengths(cfar_lengths > 0));
valid_cfar    = cfar_lengths == cfar_mode_len;
cfar_matrix   = cell2mat(cfar_data(valid_cfar));
if cfar_mode_len < n_bins_fft
    cfar_matrix = [cfar_matrix, ...
                   nan(size(cfar_matrix,1), n_bins_fft - cfar_mode_len)];
end

% Shift FFT dB scale to align with CFAR dB scale
mean_cfar_db = mean(cfar_matrix, 1, 'omitnan');
std_cfar_db  = std(cfar_matrix,  0, 1, 'omitnan');

mean_fft_db  = mean(fft_db, 1, 'omitnan');
std_fft_db   = std(fft_db,  0, 1, 'omitnan');

cfar_max     = max(mean_cfar_db, [], 'omitnan');
fft_max      = max(mean_fft_db,  [], 'omitnan');
shift        = cfar_max - fft_max;
mean_fft_db  = mean_fft_db + shift;
std_fft_db   = std_fft_db  + shift;

% Detections where FFT exceeds CFAR
detection_mask   = mean_fft_db > mean_cfar_db;
detection_ranges = range_axis(detection_mask);
detection_amps   = mean_fft_db(detection_mask);

fprintf('\nNumber of range bins with detections: %d / %d\n', ...
        sum(detection_mask), n_bins_fft);

figure;
hold on;
fill([range_axis, fliplr(range_axis)], ...
     [mean_fft_db + std_fft_db, fliplr(mean_fft_db - std_fft_db)], ...
     [0.8 0.8 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
fill([range_axis, fliplr(range_axis)], ...
     [mean_cfar_db + std_cfar_db, fliplr(mean_cfar_db - std_cfar_db)], ...
     [1.0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(range_axis, mean_fft_db,  'b-', 'LineWidth', 1.2, ...
     'DisplayName', 'Mean FFT amplitude (normalised)');
plot(range_axis, mean_cfar_db, 'r-', 'LineWidth', 1.2, ...
     'DisplayName', 'Mean CFAR threshold');
scatter(detection_ranges, detection_amps, 20, 'g', 'filled', ...
        'DisplayName', 'Detection (FFT > CFAR)');
xline(canopy_start, 'k--', 'Canopy front', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)');
ylabel('Amplitude (dB)');
title('FFT Amplitude and CFAR Threshold vs Range');
legend('FFT ±1\sigma', 'CFAR ±1\sigma', ...
       'Mean FFT amplitude (normalised)', ...
       'Mean CFAR threshold', ...
       'Detection (FFT > CFAR)', ...
       'Location', 'northeast');
xlim([5, 60]);

% CFAR offset in bins
cfar_offset = round((n_bins_fft - cfar_mode_len) / 2);
cfar_range_axis = range_axis(cfar_offset+1 : cfar_offset+cfar_mode_len);

% Plot CFAR on its own range axis
figure;
hold on;
plot(range_axis,      mean_fft_db,              'b-', 'LineWidth', 1.2, ...
     'DisplayName', 'Mean FFT (normalised)');
plot(cfar_range_axis, mean_cfar_db(1:cfar_mode_len), 'r-', 'LineWidth', 1.2, ...
     'DisplayName', 'Mean CFAR threshold');
xline(canopy_start, 'k--', 'Canopy front', 'HandleVisibility', 'off');
xline(canopy_end,   'k--', 'Canopy back',  'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Range (m)');
ylabel('Amplitude (dB)');
title('FFT Amplitude and CFAR Threshold vs Range (aligned)');
legend('Location', 'northeast');
xlim([5, 60]);

%% Figure 5a: FFT amplitude histogram with canopy range indicated

% Flatten FFT matrix and corresponding range axis
fft_db_flat    = fft_db(:);
range_repeated = repmat(range_axis, size(fft_db, 1), 1);
range_flat     = range_repeated(:);

% Split into three regions
in_pre_canopy  = range_flat < canopy_start;
in_canopy      = range_flat >= canopy_start & range_flat <= canopy_end;
in_post_canopy = range_flat > canopy_end;

figure;
hold on;
histogram(fft_db_flat(in_pre_canopy),  'BinWidth', 1, ...
          'FaceColor', [0.2 0.6 1.0], 'FaceAlpha', 0.6, ...
          'DisplayName', sprintf('Pre-canopy (< %.1f m)', canopy_start));
histogram(fft_db_flat(in_canopy),      'BinWidth', 1, ...
          'FaceColor', [0.1 0.7 0.1], 'FaceAlpha', 0.6, ...
          'DisplayName', sprintf('In-canopy (%.1f-%.1f m)', canopy_start, canopy_end));
histogram(fft_db_flat(in_post_canopy), 'BinWidth', 1, ...
          'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.6, ...
          'DisplayName', sprintf('Post-canopy (> %.1f m)', canopy_end));
hold off;
grid on;
xlabel('FFT Amplitude (dB)');
ylabel('Count');
title('FFT Amplitude Distribution by Canopy Region');
legend('Location', 'northwest');

%% Figure 5b: Raw FFT waterfall - all scan rows vs range

%% Filtered waterfall — clip colour axis for clarity
figure;
imagesc(range_axis, 1:size(fft_db,1), fft_db);
colormap jet;
c = colorbar;
c.Label.String = 'FFT Amplitude (dB)';
xlabel('Range (m)');
ylabel('Scan row');
title('Raw FFT Amplitude - All Scan Rows vs Range (clipped)');
xline(canopy_start, 'w--', 'Canopy front', 'LineWidth', 1.5, ...
      'LabelVerticalAlignment', 'bottom');
xline(canopy_end,   'w--', 'Canopy back',  'LineWidth', 1.5, ...
      'LabelVerticalAlignment', 'bottom');
xlim([5, 60]);
clim([-30, 0]);  % clip: ignore noise floor below -120, saturate strong returns above -40
% Replace -Inf with a floor value before plotting
fft_db_clean = fft_db;
fft_db_clean(isinf(fft_db_clean)) = -220;

fprintf('FFT dB min: %.1f, max: %.1f, mean: %.1f\n', ...
        min(fft_db_clean(:)), max(fft_db_clean(:)), mean(fft_db_clean(:)));
%% Figure 6: Range vs Amplitude for middle elevation and azimuth row

% Get unique sorted elevation and azimuth values
unique_el = unique(data(:,1));
unique_az = unique(data(:,2));

% Pick the middle value in each
mid_el = unique_el(round(end/2));
mid_az = unique_az(round(end/2));

fprintf('\nMiddle elevation: %.4f rad (%.2f deg)\n', mid_el, rad2deg(mid_el));
fprintf('Middle azimuth:   %.4f rad (%.2f deg)\n', mid_az, rad2deg(mid_az));

% Snap each point to its nearest unique elevation and azimuth
[~, el_idx] = min(abs(unique_el - mid_el));
[~, az_idx] = min(abs(unique_az - mid_az));

snapped_el = unique_el(el_idx);
snapped_az = unique_az(az_idx);

fprintf('Snapped elevation: %.4f rad (%.2f deg)\n', snapped_el, rad2deg(snapped_el));
fprintf('Snapped azimuth:   %.4f rad (%.2f deg)\n', snapped_az, rad2deg(snapped_az));

% Use exact match on snapped values
row_mask = data(:,1) == snapped_el & data(:,2) == snapped_az;
fprintf('Points in selected row (exact snap): %d\n', sum(row_mask));

% If still empty, relax to just the middle elevation beam across all azimuths
if sum(row_mask) < 2
    fprintf('Falling back to middle elevation only...\n');
    row_mask = data(:,1) == snapped_el;
    fprintf('Points at middle elevation (all azimuths): %d\n', sum(row_mask));
end

row_range = data(row_mask, 4);
row_amp   = data(row_mask, 3);

% Sort by range
[row_range, sort_idx] = sort(row_range);
row_amp = row_amp(sort_idx);

% Bin into 0.5 m range bins, take max amplitude per bin
bin_width = 0.5;
edges     = 0:bin_width:ceil(max(row_range));
bin_ctrs  = edges(1:end-1) + bin_width/2;
bin_amp   = nan(1, length(bin_ctrs));

for i = 1:length(bin_ctrs)
    in_bin = row_range >= edges(i) & row_range < edges(i+1);
    if any(in_bin)
        bin_amp(i) = max(row_amp(in_bin));
    end
end

valid = ~isnan(bin_amp);

figure;
bar(bin_ctrs(valid), bin_amp(valid), 1, ...
    'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
hold on;
xline(canopy_start, 'k--', 'Canopy front', ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
xline(canopy_end,   'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('Range (m)');
ylabel('Peak Return Amplitude (dB)');
title(sprintf('Range vs Amplitude — Middle Ray\n(El = %.2f°, Az = %.2f°)', ...
      rad2deg(snapped_el), rad2deg(snapped_az)));
xlim([0, max(bin_ctrs(valid)) + 1]);
%% Figure 7: Raw FFT profile for specific row (El=-1.385442, Az=0.536165)

% Parse elevation and azimuth from raw lines (needed to index fft_data)
raw_el = zeros(n_rows, 1);
raw_az = zeros(n_rows, 1);

for i = 1:n_rows
    cols = strsplit(raw_lines{i}, ',');
    if length(cols) >= 2
        raw_el(i) = str2double(cols{1});
        raw_az(i) = str2double(cols{2});
    end
end

% Find the raw row matching target elevation and azimuth
target_el = -1.385442;
target_az =  0.536165;

dist = abs(raw_el - target_el) + abs(raw_az - target_az);
dist(fft_lengths ~= mode(fft_lengths)) = inf;
[~, row_idx] = min(dist);

fprintf('Selected raw row index: %d\n', row_idx);
fprintf('El = %.4f rad (%.2f deg), Az = %.4f rad (%.2f deg)\n', ...
        raw_el(row_idx), rad2deg(raw_el(row_idx)), ...
        raw_az(row_idx), rad2deg(raw_az(row_idx)));

fft_vals   = fft_data{row_idx};
fft_db_row = double(fft_vals);

bin_range      = 0.125;
n_bins_row     = length(fft_vals);
range_axis_row = (0:n_bins_row-1) * bin_range + (2 * bin_range);

figure;
b = bar(range_axis_row, fft_db_row, 1, ...
    'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
b.BaseValue = min(fft_db_row) - 5;
hold on;
xline(canopy_start, 'k--', 'Canopy front', ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
xline(canopy_end,   'k--', 'Canopy back', ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('Range (m)');
ylabel('FFT Amplitude (dB)');
title(sprintf('Raw FFT Profile — El = %.2f°, Az = %.2f°', ...
      rad2deg(raw_el(row_idx)), rad2deg(raw_az(row_idx))));
xlim([0, max(range_axis_row) + 1]);
ylim([min(fft_db_row) - 5, 0]);

%% Count points at target elevation

%% 
target_el = -1.385442;
el_tol    = 1e-4;   % tight tolerance for exact match

n_at_el = sum(abs(raw_el - target_el) <= el_tol);
fprintf('Number of raw rows at elevation %.4f rad (%.2f deg): %d\n', ...
        target_el, rad2deg(target_el), n_at_el);

% Also show the unique azimuths at this elevation
az_at_el = raw_az(abs(raw_el - target_el) <= el_tol);
fprintf('Unique azimuths at this elevation: %d\n', length(unique(az_at_el)));
fprintf('Azimuth range: %.2f deg to %.2f deg\n', ...
        rad2deg(min(az_at_el)), rad2deg(max(az_at_el)));