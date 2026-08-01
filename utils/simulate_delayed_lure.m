function out = simulate_delayed_lure(A0, A1, B, C, tau, pi_fun, Tspan, phi, x0, dt)
%SIMULATE_DELAYED_LURE  Simulate one trajectory of a delayed Lur'e system.
%       xdot = A0 x(t) + A1 x(t-tau) + B * pi_fun(C x(t-tau), t)
%
%   Inputs
%     A0,A1,B,C : system matrices
%     tau       : scalar delay >= 0
%     pi_fun    : handle u = pi_fun(y_delayed, t), y_delayed = C*x(t-tau)
%     Tspan     : [t0 tf]
%     phi       : history handle x = phi(t) for t <= t0  ([] => constant x0)
%     x0        : initial state (used only when phi is [])
%     dt        : output sampling step (optional)
%   Output struct OUT with fields t, x, y, u, sol.
%
%   Example (constant NN controller from the repo):
%       net = load_nn();  pifun = @(y,t) nn_forward(y, net);
%       out = simulate_delayed_lure([-8 2 1;3 -10 2;1 2 -8], 3*ones(3), ...
%             [1;1;1], [1 1 1], 2, pifun, [0 20], [], [2;2;2], 1e-2);
%       plot(out.t, out.y);
if nargin < 10 || isempty(dt), dt = 1e-3*(Tspan(2)-Tspan(1)+eps) + 1e-3; end
if isempty(phi), phi = @(t) x0; end

dde_rhs = @(t,x,Z) A0*x + A1*Z(:,1) + B*pi_fun(C*Z(:,1), t);
sol = dde23(dde_rhs, tau, phi, Tspan, ddeset('RelTol',1e-7,'AbsTol',1e-9));

t = (Tspan(1):dt:Tspan(2)).';
x = deval(sol, t).';
T0 = Tspan(1);  Nt = numel(t);
y = zeros(Nt, size(C,1));  u = zeros(Nt, size(B,2));
for k = 1:Nt
    tlag = t(k) - tau;
    if tlag < T0 - 1e-12, xlag = phi(tlag); else, xlag = deval(sol, tlag); end
    y(k,:) = (C*xlag).';
    u(k,:) = (pi_fun(C*xlag, t(k))).';
end
out.t = t;  out.x = x;  out.y = y;  out.u = u;  out.sol = sol;
end
