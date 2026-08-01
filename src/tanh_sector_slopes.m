function [Alpha, Beta, flag] = tanh_sector_slopes(vlow, vup)
%TANH_SECTOR_SLOPES  Per-neuron local sector slopes for a tanh activation.
%   [Alpha, Beta, flag] = TANH_SECTOR_SLOPES(vlow, vup) returns diagonal slope
%   matrices such that, on each interval [vlow(i), vup(i)],
%
%           Beta(i,i)*v <= tanh(v) <= Alpha(i,i)*v .
%
%   Three cases (paper Fig. 1):
%     (i)   v in [vlow, vup] with vlow, vup >= 0
%     (ii)  v in [vlow, vup] with vlow, vup <= 0
%     (iii) sign-crossing interval (vlow < 0 < vup): the slope is 1 and the
%           neuron is marked in FLAG so the caller can take absolute values.
vlow = vlow(:);  vup = vup(:);
n = numel(vlow);
alpha = zeros(n,1);  beta = zeros(n,1);  flag = zeros(n,1);
for i = 1:n
    if vlow(i)*vup(i) < 0            % (iii) sign-crossing
        alpha(i) = 1;  beta(i) = 1;  flag(i) = 1;
    elseif vup(i) > 0               % (i) both non-negative
        alpha(i) = 1;
        beta(i)  = tanh(vup(i))/vup(i);
    elseif vlow(i) < 0              % (ii) both non-positive
        alpha(i) = tanh(vlow(i))/vlow(i);
        beta(i)  = 1;
    end
end
Alpha = diag(alpha);
Beta  = diag(beta);
end
