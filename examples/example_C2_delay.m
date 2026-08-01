%EXAMPLE_C2_DELAY  Risk configuration C2: state/input delay, no interval.
%   1) Computes the network sector [gamma1, gamma2] and the DELAY-INDEPENDENT
%      positivity certificate (Theorem 3) for the delayed Lur'e system
%           xdot = A0 x(t) + A1 x(t-tau) + B * Phi(C x(t-tau)).
%   2) Monte-Carlo check for four delays tau in {0.2, 2, 8, 16} s, each with
%      100 random constant histories. The certificate holds for ALL tau >= 0.
%   Reproduces Fig. 5 and exports mc_delay_nn.pdf.
clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

net = load_nn();
phi = @(y) nn_forward(y, net);

% ----- Delayed plant -----
A0 = [-8 2 1; 3 -10 2; 1 2 -8];
A1 = 3*ones(3);
B  = [1;1;1];
C  = [1 1 1];

% ----- Delay-independent certificate (Theorem 3) -----
[g1, g2] = sector_bound_ffnn(net, 0, 4.5);
[pass, info] = certify_positive_lure(A0+A1, A0+A1, B, C, g1, g2);
fprintf('C2 sector: gamma1=%.3f, gamma2=%.3f\n', g1, g2);
fprintf('   Metzler margin = %.4f -> %s;  Hurwitz max Re = %.4f -> %s\n', ...
        info.metzler_margin, ternary(info.metzler_ok,'OK','FAIL'), ...
        info.max_real_eig,   ternary(info.hurwitz_ok,'OK','FAIL'));
fprintf('   CERTIFICATE: %s (valid for every delay tau >= 0)\n\n', ...
        ternary(pass,'PASSED','NOT satisfied'));

% ----- Monte-Carlo -----
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
LW_thin=0.6; LW_bold=1.6; FS_ax=13; FS_lbl=18; FS_ttl=18;
COL_thin=[0.70 0.70 0.70]; COL_bold=[0 0 0];
NPTS_thin=300; NPTS_bold=1200;

rng(11);
tau_list=[0.2 2 8 16]; numIC=100;
TendSim=20; TendPlot=20; x0_box=1.5; tol_final=1;
opts = ddeset('RelTol',1e-8,'AbsTol',1e-10);

fig = figure('Color','w','Units','centimeters','Position',[2 2 18 13.5]);
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax1 = gobjects(numel(tau_list),1); ymax1 = zeros(numel(tau_list),1);

fprintf('--- C2 delay-only Monte-Carlo ---\n');
for k = 1:numel(tau_list)
    tau = tau_list(k);
    nexttile(tl,k); ax=gca; ax1(k)=ax; hold(ax,'on'); grid(ax,'on'); ax.FontSize=FS_ax;
    norms_end=zeros(numIC,1); t_cell=cell(numIC,1); xn_cell=cell(numIC,1); x0_list=cell(numIC,1);

    ddefun = @(t,x,Z) A0*x + A1*Z(:,1) + B*phi(C*Z(:,1));
    for i = 1:numIC
        x0 = -x0_box + 2*x0_box*rand(3,1);
        sol = dde23(ddefun, tau, @(t) x0, [0 TendSim], opts);
        tp  = linspace(0, TendPlot, NPTS_thin);
        yp  = deval(sol, tp);
        xn  = vecnorm(yp,2,1)';
        plot(ax, tp, xn, 'Color', COL_thin, 'LineWidth', LW_thin);
        t_cell{i}=tp; xn_cell{i}=xn;
        norms_end(i) = norm(deval(sol,TendSim),2);
        x0_list{i}   = x0;
    end
    [~,im] = min(abs(norms_end - median(norms_end)));
    sol_m = dde23(ddefun, tau, @(t) x0_list{im}, [0 TendSim], opts);
    tb = linspace(0, TendPlot, NPTS_bold);
    plot(ax, tb, vecnorm(deval(sol_m,tb),2,1)', 'Color', COL_bold, 'LineWidth', LW_bold);

    bad = find(isnan(norms_end)|isinf(norms_end)|(norms_end>tol_final));
    xlabel(ax,'$t$ (s)','FontSize',FS_lbl); ylabel(ax,'$\|x(t)\|_2$','FontSize',FS_lbl);
    title(ax, sprintf('$\\tau=%.1f$ s: stable %d/%d', tau, numIC-numel(bad), numIC), 'FontSize', FS_ttl);
    ymax1(k) = max(cellfun(@(v) max(v), xn_cell));
    fprintf('  tau=%.1f: stable %d/%d\n', tau, numIC-numel(bad), numIC);
end
for k=1:numel(tau_list), xlim(ax1(k),[0 TendPlot]); ylim(ax1(k),[0 max(max(ymax1),1e-3)]); end
exportgraphics(fig, fullfile(here,'mc_delay_nn.pdf'), 'ContentType','vector');

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
