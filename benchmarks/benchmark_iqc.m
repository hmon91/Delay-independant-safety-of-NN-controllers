%BENCHMARK_IQC  Runtime of the IQC/SDP baseline (Table I).
%   Runs VERIFY_IQC for the two discretization/delay configurations reported in
%   the paper and averages the solve time. The IQC test is configuration
%   dependent: at the coarse step it typically CANNOT certify, and at the fine
%   step it certifies but is orders of magnitude slower than the positivity test.
%
%   REQUIRES: CVX + SeDuMi.
clc; clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'), fullfile(here, '..', 'verification'));

configs = [0.10, 0.70;    % dt, tau   (coarse)  -> typically "cannot certify"
           0.01, 0.07];   % dt, tau   (fine)     -> certifies, ~2 s
Nrep = 5;

for c = 1:size(configs,1)
    dt = configs(c,1);  tau = configs(c,2);
    rt = zeros(Nrep,1);  status = 'n/a';
    for k = 1:Nrep
        if exist('cvx_clear','file'), cvx_clear; end
        t0 = tic;
        status = verify_iqc(dt, tau);
        rt(k) = toc(t0);
    end
    fprintf('IQC  dt=%.2f tau=%.2f :  status=%-12s  mean runtime = %.3f s\n', ...
            dt, tau, status, mean(rt));
end
