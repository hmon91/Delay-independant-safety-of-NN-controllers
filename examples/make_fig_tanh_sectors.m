%MAKE_FIG_TANH_SECTORS  Reproduce Fig. 1: tanh sector relaxations (3 cases).
%   Self-contained (no network weights needed). Draws the three sector cases
%   used by SECTOR_BOUND_FFNN / TANH_SECTOR_SLOPES and exports
%   tanh_sector_triptych.pdf.
%     (i)   both endpoints >= 0
%     (ii)  both endpoints <= 0
%     (iii) sign-crossing interval (slope 1, absolute-value envelope)
close all; clear;
here = fileparts(mfilename('fullpath'));

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultLegendInterpreter',        'latex');

nu      = linspace(-3, 3, 500);
tanh_nu = tanh(nu);
xlimv   = [-3 3];
ylimv   = [-1.2 1.2];

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 10.5 3.1]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------- Tile 1: both endpoints positive ----------
ax1 = nexttile; hold(ax1, 'on');
cfg1.nu_bar = 1.5;  cfg1.nu_underline = 0.8;  cfg1.mode = "pos";
h1 = draw_sector_tile(ax1, nu, tanh_nu, cfg1, xlimv, ylimv);
title(ax1, '$(i)\; \underline{\nu},\bar{\nu} \ge 0$', 'FontSize', 20);
legend(ax1, [h1.hbeta, h1.halpha], {'$\beta \nu$', '$\alpha \nu$'}, 'Location', 'best', 'FontSize', 20);

% ---------- Tile 2: both endpoints negative ----------
ax2 = nexttile; hold(ax2, 'on');
cfg2.nu_bar = -1.5;  cfg2.nu_underline = -0.8;  cfg2.mode = "neg";
h2 = draw_sector_tile(ax2, nu, tanh_nu, cfg2, xlimv, ylimv);
title(ax2, '$(ii)\; \underline{\nu},\bar{\nu} \le 0$', 'FontSize', 20);
legend(ax2, [h2.hbeta, h2.halpha], {'$\beta \nu$', '$\alpha \nu$'}, 'Location', 'best', 'FontSize', 20);

% ---------- Tile 3: mixed (sign-crossing) ----------
ax3 = nexttile; hold(ax3, 'on');
cfg3.nu_bar = 1.5;  cfg3.nu_underline = -1;  cfg3.mode = "mixed";
h3 = draw_sector_tile(ax3, nu, tanh_nu, cfg3, xlimv, ylimv);
title(ax3, '$(iii)\; \underline{\nu} < 0,\bar{\nu} > 0$', 'FontSize', 20);
legend(ax3, [h3.hbeta1, h3.halpha2], {'$\beta \nu$', '$\alpha \nu$'}, ...
       'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 20);

% Horizontal y-axis labels near the top of each axis
for ax = [ax1 ax2 ax3]
    yl = ax.YLabel; yl.Interpreter = 'latex'; yl.String = '$f(\nu)$'; yl.Units = 'data';
    xl = xlim(ax); ylm = ylim(ax); dx = diff(xl); dy = diff(ylm);
    yl.Position(1) = 0 - 0.06*dx;  yl.Position(2) = ylm(2) - 0.04*dy;
    yl.Rotation = 0; yl.HorizontalAlignment = 'right'; yl.VerticalAlignment = 'top';
end

set(fig, 'Renderer', 'opengl');
exportgraphics(fig, fullfile(here, 'tanh_sector_triptych.pdf'), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

% ============================ helper ============================
function H = draw_sector_tile(ax, nu, tanh_nu, cfg, xlimv, ylimv)
H.htanh = plot(ax, nu, tanh_nu, 'b', 'LineWidth', 2);
xlim(ax, xlimv); ylim(ax, ylimv); grid(ax, 'on'); box(ax, 'off');
ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 1.2; ax.TickLength = [0 0];
ax.XTickLabel = []; ax.YTickLabel = [];
xlabel(ax, '$\nu$', 'FontSize', 20); ylabel(ax, '$f(\nu)$', 'FontSize', 20);
plot(ax, xlimv, [0 0], 'k', 'LineWidth', 1.5);
plot(ax, [0 0], ylimv, 'k', 'LineWidth', 1.5);
fill(ax, [2.9, 3, 2.9], [-0.05, 0, 0.05], 'k', 'EdgeColor', 'k');
fill(ax, [-0.05, 0, 0.05], [1.1, 1.2, 1.1], 'k', 'EdgeColor', 'k');

switch cfg.mode
    case "pos"
        nb = cfg.nu_bar; nu_ = cfg.nu_underline;
        sb = tanh(nb)/nb; su = tanh(nu_)/nu_;
        H.hbeta  = plot(ax, nu, sb*nu, '--r', 'LineWidth', 1.5);
        H.halpha = plot(ax, nu, su*nu, '--g', 'LineWidth', 1.5);
        plot(ax, nb, tanh(nb), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        plot(ax, nu_, tanh(nu_), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        plot(ax, [nb nb], [0 tanh(nb)], '--k', 'LineWidth', 1.2);
        plot(ax, [nu_ nu_], [0 tanh(nu_)], '--k', 'LineWidth', 1.2);
        nf = linspace(nu_, nb, 200);
        fill(ax, [nf, fliplr(nf)], [sb*nf, fliplr(su*nf)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        text(ax, nb, -0.18, '$\overline{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
        text(ax, nu_, -0.18, '$\underline{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
    case "neg"
        nb = cfg.nu_bar; nu_ = cfg.nu_underline;
        sb = tanh(nb)/nb; su = tanh(nu_)/nu_;
        H.halpha = plot(ax, nu, sb*nu, '--g', 'LineWidth', 1.5);
        H.hbeta  = plot(ax, nu, su*nu, '--r', 'LineWidth', 1.5);
        plot(ax, nb, tanh(nb), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        plot(ax, nu_, tanh(nu_), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        plot(ax, [nb nb], [0 tanh(nb)], '--k', 'LineWidth', 1.2);
        plot(ax, [nu_ nu_], [0 tanh(nu_)], '--k', 'LineWidth', 1.2);
        nf = linspace(nu_, nb, 200);
        fill(ax, [nf, fliplr(nf)], [sb*nf, fliplr(su*nf)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        text(ax, nb, 0.18, '$\underline{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
        text(ax, nu_, 0.18, '$\overline{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
    case "mixed"
        nb = cfg.nu_bar; nu_ = cfg.nu_underline;
        beta2 = tanh(nb)/nb; beta1 = tanh(nu_)/nu_; alpha2 = 1; alpha1 = 1;
        nu_pos = linspace(0, nb, 250); nu_neg = linspace(nu_, 0, 250);
        nu_pos1 = linspace(0, 3, 250);  nu_neg1 = linspace(-3, 0, 250);
        ya = beta1*nu_neg; yb = alpha1*nu_neg;
        H.halpha1 = plot(ax, nu_neg1, beta1*nu_neg1, ':g', 'LineWidth', 2.5);
        H.hbeta1  = plot(ax, nu_neg1, alpha1*nu_neg1, ':r', 'LineWidth', 2.5);
        fill(ax, [nu_neg, fliplr(nu_neg)], [max(ya,yb), fliplr(min(ya,yb))], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        H.halpha2 = plot(ax, nu_pos1, alpha2*nu_pos1, '--g', 'LineWidth', 1.5);
        H.hbeta2  = plot(ax, nu_pos1, beta2*nu_pos1, '--r', 'LineWidth', 1.5);
        fill(ax, [nu_pos, fliplr(nu_pos)], [beta2*nu_pos, fliplr(alpha2*nu_pos)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(ax, nb, tanh(nb), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        plot(ax, nu_, tanh(nu_), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        plot(ax, [nb nb], [0 tanh(nb)], '--k', 'LineWidth', 1.2);
        plot(ax, [nu_ nu_], [0 tanh(nu_)], '--k', 'LineWidth', 1.2);
        text(ax, nb, -0.18, '$\bar{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
        text(ax, nu_, 0.18, '$\underline{\nu}$', 'FontSize', 20, 'HorizontalAlignment', 'center');
end
end
