%EXAMPLE_C3_COMBINED  Risk configuration C3: simultaneous delay + interval.
%   1) Computes the network sector [gamma1, gamma2] and the delay-independent
%      positivity certificate (Theorem 4) for the interval box.
%   2) Monte-Carlo check: draws 5 random plants (A0, A1) from the interval box
%      AND a random delay tau per plant, then simulates 10 componentwise-
%      positive initial histories each, plotting the output y(t) = C x(t).
%   Reproduces Fig. 6 and exports mc_delay_interval_outputs.pdf.
%   (No pre-rendered image ships with the repo; run this script to generate it.)
clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

net = load_nn();
phi = @(y) nn_forward(y, net);

% ----- Interval plant & IO -----
A0L = [-8   2   1;   3  -10   2;   1   2  -8];
A0U = [-7.5 2.5 1.5; 4.5 -9.5  2.5; 1.5 2.5 -7.5];
A1L = 3.0*ones(3);  A1U = 3.5*ones(3);
B = [1;1;1];  C = [1 1 1];

% ----- Delay-independent certificate over the box (Theorem 4) -----
[g1, g2] = sector_bound_ffnn(net, 0, 4.5);
[pass, info] = certify_positive_lure(A0L+A1L, A0U+A1U, B, C, g1, g2);
fprintf('C3 sector: gamma1=%.3f, gamma2=%.3f\n', g1, g2);
fprintf('   Metzler margin = %.4f -> %s;  Hurwitz max Re = %.4f -> %s\n', ...
        info.metzler_margin, ternary(info.metzler_ok,'OK','FAIL'), ...
        info.max_real_eig,   ternary(info.hurwitz_ok,'OK','FAIL'));
fprintf('   CERTIFICATE: %s (all interval realizations, every delay tau >= 0)\n\n', ...
        ternary(pass,'PASSED','NOT satisfied'));

% ----- Monte-Carlo (output only) -----
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
LW_thin=0.8; FS_ax=13; FS_lbl=18;

rng(21);
numSystems=5; numIC=10; TendSim=100; TendPlot=20;
tau_min=0.2; tau_max=10; center_box=1.5; r_min=0.3; r_max=0.8;
colors = lines(numSystems);
S(numSystems) = struct('A0',[],'A1',[],'tau',[],'c',[],'r',[]);
for i=1:numSystems
    r = r_min + (r_max-r_min)*rand;
    c = r + (center_box-r)*rand(3,1);        % ensures the IC ball stays >= 0
    S(i).A0=A0L+(A0U-A0L).*rand(3);
    S(i).A1=A1L+(A1U-A1L).*rand(3);
    S(i).tau=tau_min+(tau_max-tau_min)*rand;
    S(i).c=c; S(i).r=r;
end

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 7.5]);
ax  = axes(fig); hold(ax,'on'); grid(ax,'on'); ax.FontSize=FS_ax;
for i=1:numSystems
    A0=S(i).A0; A1=S(i).A1; tau=S(i).tau; col=colors(i,:);
    ddefun = @(t,x,Z) A0*x + A1*Z(:,1) + B*phi(C*Z(:,1));
    for j=1:numIC
        x0 = sample_ball(S(i).c, S(i).r, 3);
        sol = dde23(ddefun, tau, @(t) x0, [0 TendSim], ddeset('RelTol',1e-8,'AbsTol',1e-10));
        t = sol.x(:); y = (C*sol.y).';
        plot(ax, t(t<=TendPlot), y(t<=TendPlot), 'Color', col, 'LineWidth', LW_thin);
    end
end
xlabel(ax,'$t$ (s)','FontSize',FS_lbl); ylabel(ax,'$y(t)$','FontSize',FS_lbl);
h = gobjects(numSystems,1); labels = cell(numSystems,1);
for i=1:numSystems
    labels{i} = sprintf('sys %d: $\\tau{=}%.2f\\,\\mathrm{s}$', i, S(i).tau);
    h(i) = plot(ax, nan, nan, '-', 'Color', colors(i,:), 'LineWidth', 2);
end
hL = legend(ax, h, labels, 'FontSize', 15); set(hL,'Interpreter','latex');
xlim(ax,[0 TendPlot]);
exportgraphics(fig, fullfile(here,'mc_delay_interval_outputs.pdf'), 'ContentType','vector');

% ----- local helpers -----
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
function x = sample_ball(c, r, n)
v = randn(n,1); v = v/norm(v); rho = r*rand^(1/n); x = c + rho*v;
end
