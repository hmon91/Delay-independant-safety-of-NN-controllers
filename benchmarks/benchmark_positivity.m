%BENCHMARK_POSITIVITY  Runtime of the positivity certificate (Table I).
%   Times the full positivity pipeline (network sector bound + Metzler/Hurwitz
%   test) over many repetitions and reports the average wall-clock time. On the
%   paper's example this is ~5-6e-5 s and is delay/discretization INDEPENDENT.
clc; clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

net = load_nn();
A0 = [-8 2 1; 3 -10 2; 1 2 -8];  A1 = 3*ones(3);
B  = [1;1;1];  C = [1 1 1];
y_lb = 0;  y_ub = 4.5;

Nrep = 10000;
runtimes = zeros(Nrep,1);
for k = 1:Nrep
    t0 = tic;
    [g1, g2] = sector_bound_ffnn(net, y_lb, y_ub);
    certify_positive_lure(A0+A1, A0+A1, B, C, g1, g2);
    runtimes(k) = toc(t0);
end
fprintf('Positivity certificate: mean runtime over %d runs = %.3e s\n', Nrep, mean(runtimes));
fprintf('(delay- and discretization-independent)\n');
