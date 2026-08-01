%EXAMPLE_C1_INTERVAL  Risk configuration C1: interval uncertainty, no delay.
%   1) Computes the network sector [gamma1, gamma2] and the delay-independent
%      positivity certificate (Theorem 2) for the whole interval box.
%   2) Monte-Carlo check: samples 4 random plants from the interval and, for
%      each, simulates 100 random initial conditions of
%           xdot = (A0 + A1) x + B * Phi(C x).
%   Reproduces Fig. 4 and exports mc_interval_nn_no_delay.pdf.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

net = load_nn();
phi = @(y) nn_forward(y, net);

% ----- Interval plant (no delay) -----
A0L = [-8   2   1;   3  -10   2;   1   2  -8];
A0U = [-7.5 2.5 1.5; 4.5 -9.5  2.5; 1.5 2.5 -7.5];
A1L = 3.0*ones(3);
A1U = 3.5*ones(3);
B   = [1;1;1];
C   = [1 1 1];

% ----- Certificate over the interval box (Theorem 2) -----
[g1, g2] = sector_bound_ffnn(net, 0, 4.5);
[pass, info] = certify_positive_lure(A0L+A1L, A0U+A1U, B, C, g1, g2);
fprintf('C1 sector: gamma1=%.3f, gamma2=%.3f\n', g1, g2);
fprintf('   Metzler margin (min off-diag of A_lo+B*g1*C) = %.4f  -> %s\n', ...
        info.metzler_margin, ternary(info.metzler_ok,'OK','FAIL'));
fprintf('   Hurwitz max Re eig (A_up+B*g2*C)             = %.4f  -> %s\n', ...
        info.max_real_eig, ternary(info.hurwitz_ok,'OK','FAIL'));
fprintf('   CERTIFICATE: %s\n\n', ternary(pass,'PASSED (locally exp. stable)','NOT satisfied'));

% ----- Monte-Carlo over sampled interval systems -----
rng(11);
numSystems = 4;  numIC = 100;  Tend = 5;  x0_box = 1.5;
tol_final = 1e-4;  blow_guard = 1e4;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
              'Events', @(t,x) blowup_event(t,x,blow_guard));

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontWeight','bold');
LW_thin=0.6; LW_bold=1.6; FS_ax=12; FS_lbl=18; FS_ttl=18;
COL_thin=[0.65 0.65 0.65]; COL_bold=[0 0 0];

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 13.5]);
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax_list = gobjects(numSystems,1); ymax_list = zeros(numSystems,1);

fprintf('--- C1 Monte-Carlo (per sampled system) ---\n');
for k = 1:numSystems
    A0 = A0L + (A0U-A0L).*rand(3);
    A1 = A1L + (A1U-A1L).*rand(3);
    A  = A0 + A1;

    nexttile(tl,k); ax=gca; ax_list(k)=ax; hold(ax,'on'); grid(ax,'on'); ax.FontSize=FS_ax;
    norms_end=zeros(numIC,1); t_cell=cell(numIC,1); xn_cell=cell(numIC,1); stable=0;

    for i = 1:numIC
        x0 = -x0_box + 2*x0_box*rand(3,1);
        odefun = @(t,x) A*x + B*phi(C*x);
        tspan  = linspace(0, Tend, 300);
        [t,x]  = ode45(odefun, tspan, x0, opts);
        xn = vecnorm(x,2,2);
        plot(ax, t, xn, 'Color', COL_thin, 'LineWidth', LW_thin);
        t_cell{i}=t; xn_cell{i}=xn; norms_end(i)=xn(end);
        if norms_end(i) <= tol_final && all(~isnan(xn)), stable=stable+1; end
    end
    [~,im] = min(abs(norms_end - median(norms_end)));
    plot(ax, t_cell{im}, xn_cell{im}, 'Color', COL_bold, 'LineWidth', LW_bold);
    xlabel(ax,'$t$ (s)','FontSize',FS_lbl); ylabel(ax,'$\|x(t)\|_2$','FontSize',FS_lbl);
    title(ax, sprintf('Sample \\#%d: stable %d/%d', k, stable, numIC), 'FontSize', FS_ttl);
    ymax_list(k) = max(cellfun(@(v) max(v), xn_cell));
    fprintf('  #%d: stable %d/%d\n', k, stable, numIC);
end
for k=1:numSystems, xlim(ax_list(k),[0 Tend]); ylim(ax_list(k),[0 max(max(ymax_list),1e-3)]); end
exportgraphics(fig, fullfile(here,'mc_interval_nn_no_delay.pdf'), 'ContentType','vector');

% ----- local helpers -----
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
function [v,isterminal,direction] = blowup_event(~,x,guard)
v = guard - norm(x,2); isterminal = 1; direction = -1;
end
