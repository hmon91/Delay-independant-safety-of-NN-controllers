%EXAMPLE_SECTOR_BOUND  Reproduce Fig. 3: NN mapping vs. its local sector bound.
%   Computes the network-level sector [gamma1, gamma2] over the output range
%   y in [0, 4.5] and overlays it on the true NN mapping Phi(y). Exports
%   nn_sector_plot.pdf into this folder.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

net  = load_nn();
y_lb = 0;  y_ub = 4.5;

[gamma1, gamma2] = sector_bound_ffnn(net, y_lb, y_ub);
fprintf('Local sector over y in [%.2f, %.2f]:  gamma1 = %.2f,  gamma2 = %.2f\n', ...
        y_lb, y_ub, gamma1, gamma2);

% --- true NN mapping and the two sector lines on a grid ---
y  = linspace(y_lb, y_ub, 200);
NN = arrayfun(@(v) nn_forward(v, net), y);
lo = gamma1 * y;
hi = gamma2 * y;

% --- publication-style plot ---
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
figure('Color', 'w'); hold on; box on; grid on; grid minor;
p1 = plot(y, NN, 'k',  'LineWidth', 1.6);
p2 = plot(y, lo, '--', 'Color', [0.00 0.35 0.90], 'LineWidth', 1.5);
p3 = plot(y, hi, '--', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.5);
fill([y, fliplr(y)], [lo, fliplr(hi)], [1 0.85 0.25], 'FaceAlpha', 0.22, 'EdgeColor', 'none');
xlabel('$y$', 'FontSize', 20);
ylabel('$\mathrm{Output}$', 'FontSize', 20);
xlim([y_lb, y_ub]);
set(gca, 'LineWidth', 1.1, 'FontSize', 14, 'TickDir', 'out');
legend([p1 p2 p3], {'$\mathrm{NN}(y)$', '$\gamma_1 y$', '$\gamma_2 y$'}, ...
       'Location', 'northeast', 'FontSize', 16);
title(sprintf('$y\\in[%.1f,%.1f],\\;\\gamma_1=%.2f,\\;\\gamma_2=%.2f$', ...
      y_lb, y_ub, gamma1, gamma2), 'FontSize', 16);

set(gcf, 'Renderer', 'painters');
exportgraphics(gca, fullfile(here, 'nn_sector_plot.pdf'), 'ContentType', 'vector');
