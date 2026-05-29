% =========================================================================
% radar_vs_lidar_plane_tree.m
% Overlays radar and LiDAR normalised point density through plane tree
% =========================================================================

clear; clc; close all;

load('radar_cone_density.mat');
load('lidar_cone_density.mat');

figure('Color','w','Units','normalized','Position',[0.1 0.3 0.6 0.45]);

plot(lidar_r_centres, lidar_counts_norm, '-', ...
    'Color', [0.2 0.6 1.0], 'LineWidth', 2.0, ...
    'DisplayName', 'LiDAR 905 nm');
hold on;
plot(radar_r_centres, radar_counts_norm, '-', ...
    'Color', [0.9 0.3 0.1], 'LineWidth', 2.0, ...
    'DisplayName', 'Radar 120 GHz');

xline(14, '--k', 'LineWidth', 1.2, 'Label', 'Canopy front', ...
    'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');
xline(40, '--k', 'LineWidth', 1.2, 'Label', 'Canopy back', ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'none', 'HandleVisibility', 'off');

xlabel('Range from sensor (m)', 'FontSize', 12, 'Interpreter', 'none');
ylabel('Normalised point density', 'FontSize', 12, 'Interpreter', 'none');
title('Radar vs LiDAR point density through plane tree canopy', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 11, 'Interpreter', 'none');
grid on;
xlim([5 50]);
ylim([0 1.05]);

exportgraphics(gcf, 'fig_radar_vs_lidar_plane_tree.png', 'Resolution', 300);
fprintf('Saved: fig_radar_vs_lidar_plane_tree.png\n');