function [gamma1, gamma2] = sector_bound_ffnn(net, ylb, yub)
%SECTOR_BOUND_FFNN  Network-level local sector bound for a FFNN (Theorem 1).
%   [gamma1, gamma2] = SECTOR_BOUND_FFNN(net, ylb, yub) computes matrices such
%   that the network output satisfies the elementwise sandwich
%
%           gamma1 * y  <=  Phi(y)  <=  gamma2 * y      for all  ylb <= y <= yub.
%
%   Unlike elementwise relaxations (CROWN, IBP), this returns a single
%   *network-level* [gamma1, gamma2] sector suitable for Lur'e/positive-system
%   analysis. Only 'tanh' and 'linear' layers are supported.
%
%   The construction follows the paper:
%     Step A  - propagate a linear-in-y enclosure through each weight using the
%               W = W+ + W- split (see PREACTIVATION_BOUNDS);
%     Step B  - relax each tanh with per-neuron sector slopes, handling the
%               sign-crossing case by row replacement (see TANH_SECTOR_SLOPES).
%
%   Example (reproduces Fig. 3, gamma1 ~ -3.01, gamma2 ~ -2.42):
%       net = load_nn();
%       [g1, g2] = sector_bound_ffnn(net, 0, 4.5);
%
%   See also PREACTIVATION_BOUNDS, TANH_SECTOR_SLOPES, CERTIFY_POSITIVE_LURE.
W   = net.W;
act = net.act;
q   = numel(W);
ylb = ylb(:);  yub = yub(:);

% ---- (1) Pre-activation ranges per layer (numeric; used to pick tanh slopes) ----
Lv = cell(q,1);  Uv = cell(q,1);
[Lv{1}, Uv{1}] = preactivation_bounds(W{1}, ylb, yub);
for i = 2:q
    [Lw, Uw] = act_range(act{i-1}, Lv{i-1}, Uv{i-1});
    [Lv{i}, Uv{i}] = preactivation_bounds(W{i}, Lw, Uw);
end

% ---- (2) Linear-in-y sector matrices, propagated layer by layer ----
Mup = [];  Mlo = [];
for i = 1:q
    % --- weight step (pre-activation enclosure) ---
    if i == 1
        Uvm = W{1};   Lvm = W{1};              % first pre-activation is exact
    else
        Wp = max(W{i}, 0);  Wm = min(W{i}, 0);
        Uvm = Wp*Mup + Wm*Mlo;
        Lvm = Wp*Mlo + Wm*Mup;
    end
    % --- activation step ---
    switch lower(act{i})
        case 'tanh'
            [A, B, F] = tanh_sector_slopes(Lv{i}, Uv{i});
            Mup = replace_rows(F, A*Uvm, +1);
            Mlo = replace_rows(F, B*Lvm, -1);
        case {'linear', 'none', ''}
            Mup = Uvm;  Mlo = Lvm;
        otherwise
            error('sector_bound_ffnn:activation', ...
                  'Activation "%s" not supported (use tanh or linear).', act{i});
    end
end

gamma1 = Mlo;
gamma2 = Mup;
end

% ------------------------------------------------------------------------------
function [Lw, Uw] = act_range(actName, Lv, Uv)
%ACT_RANGE  Propagate a numeric interval through a monotone activation.
switch lower(actName)
    case 'tanh',                Lw = tanh(Lv);  Uw = tanh(Uv);
    case {'linear','none',''},  Lw = Lv;        Uw = Uv;
    otherwise
        error('sector_bound_ffnn:activation', 'Unknown activation "%s".', actName);
end
end

function M = replace_rows(flag, M, s)
%REPLACE_ROWS  Rows flagged (sign-crossing) are replaced by s*abs(row).
idx = find(flag ~= 0);
M(idx,:) = s * abs(M(idx,:));
end
