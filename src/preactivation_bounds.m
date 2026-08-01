function [L, U] = preactivation_bounds(W, xmin, xmax)
%PREACTIVATION_BOUNDS  Elementwise bounds of W*x over a box xmin <= x <= xmax.
%   [L, U] = PREACTIVATION_BOUNDS(W, xmin, xmax) returns the tightest
%   elementwise lower/upper bounds of the linear map W*x when each coordinate
%   x_i lies in [xmin(i), xmax(i)]. This is the W = W+ + W- interval split used
%   in Step A of the network-level sector bound.
%
%   Scalars are accepted for xmin/xmax (broadcast to all coordinates).
xmin = xmin(:);  xmax = xmax(:);
Wpos = (W >= 0);
Wneg = (W <  0);
% Lower bound: use xmin where the weight is >= 0, xmax where it is < 0.
L = sum(Wpos .* (W .* xmin.') + Wneg .* (W .* xmax.'), 2);
% Upper bound: use xmax where the weight is >= 0, xmin where it is < 0.
U = sum(Wpos .* (W .* xmax.') + Wneg .* (W .* xmin.'), 2);
end
