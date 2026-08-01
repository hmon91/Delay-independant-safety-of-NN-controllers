function [status, radii, info] = verify_iqc(dt, tau, deltav1)
%VERIFY_IQC  IQC/SDP baseline for the delayed + interval NN feedback loop.
%   [status, radii, info] = VERIFY_IQC(dt, tau, deltav1) runs the
%   state-of-the-art IQC pipeline (adapted from Yin-Seiler-Arcak, ref. [6])
%   tailored to this paper's risk model: a single integer delay d = round(tau/dt)
%   plus entrywise interval uncertainty on A0, A1, with the tanh network handled
%   by a local-sector IQC. Solves one LMI to certify a region of attraction.
%
%   Inputs (all optional):
%     dt       - discretization step        (default 0.01)
%     tau      - physical delay in seconds   (default 0.07)
%     deltav1  - local box half-width on v1  (default 4.5)
%   Outputs:
%     status   - CVX solver status string ('Solved', 'Infeasible', ...)
%     radii    - principal radii of the ROA x-ellipsoid ([] if infeasible)
%     info     - struct with P and the parsed feasibility flag
%
%   REQUIRES: CVX (http://cvxr.com) with the SeDuMi solver. This is the paper's
%   IQC comparison; it is the original formulation reorganized as a function.
%   Contrast with the positivity certificate in CERTIFY_POSITIVE_LURE, which is
%   an LMI-free, delay-independent linear test.
if nargin < 1 || isempty(dt),      dt = 0.01;  end
if nargin < 2 || isempty(tau),     tau = 0.07; end
if nargin < 3 || isempty(deltav1), deltav1 = 4.5; end

here = fileparts(mfilename('fullpath'));
dataDir = fullfile(here, '..', 'data');

%% Continuous-time interval data (same box as example C1)
A0L = [-8   2   1;   3 -10  2; 1  2 -8];
A0U = [-7.5 2.5 1.5; 4.5 -9.5 2.5; 1.5 2.5 -7.5];
A1L = 3.0*ones(3);
A1U = 3.5*ones(3);
Bc  = [1;1;1];
A0c = 0.5*(A0L + A0U);  A1c = 0.5*(A1L + A1U);
R0  = 0.5*(A0U - A0L);  R1  = 0.5*(A1U - A1L);
nG  = size(A0c,1);  nu = 1;
C   = [1 1 1];

%% NN weights (network sees y = C x at the first layer)
W1 = readmatrix(fullfile(dataDir,'W1.csv'));
W2 = readmatrix(fullfile(dataDir,'W2.csv'));
W3 = readmatrix(fullfile(dataDir,'W3.csv'));
W1 = W1*C;
n1 = size(W1,1); n2 = size(W2,1); n3 = size(W3,1);
nphi = n1 + n2;
b1 = zeros(n1,1); b2 = zeros(n2,1); b3 = zeros(n3,1);

%% Discretization + single delay port on [x; u]
d   = max(1, round(tau/dt));
A0d = eye(nG) + dt*A0c;  A1d = dt*A1c;  Bd = dt*Bc;
AG  = A0d;  BG1 = [A1d, Bd];  BG2 = zeros(nG,nu);  nq = nG + nu;
CG  = [eye(nG); zeros(nu,nG)];  DG1 = zeros(nq,nq);  DG2 = [zeros(nG,nu); eye(nu)];

%% NN local-sector bounds
xeq = zeros(nG,1);
v1eq = W1*xeq + b1;  w1eq = tanh(v1eq);
v2eq = W2*w1eq + b2; w2eq = tanh(v2eq);
v3eq = W3*w2eq + b3;
v1up = v1eq + deltav1;  v1lb = v1eq - deltav1;
alpha1 = min((tanh(v1up)-tanh(v1eq))./(v1up-v1eq), (tanh(v1eq)-tanh(v1lb))./(v1eq-v1lb));
beta = 1;
v2up = zeros(n2,1); v2lb = zeros(n2,1);
for i = 1:n2
    v2up(i) = W2(i,:)*0.5*(tanh(v1up)+tanh(v1lb)) + b2(i) + abs(W2(i,:))*0.5*abs(tanh(v1up)-tanh(v1lb));
    v2lb(i) = W2(i,:)*0.5*(tanh(v1up)+tanh(v1lb)) + b2(i) - abs(W2(i,:))*0.5*abs(tanh(v1up)-tanh(v1lb));
end
alpha2 = min((tanh(v2up)-tanh(v2eq))./(v2up-v2eq), (tanh(v2eq)-tanh(v2lb))./(v2eq-v2lb));
Alpha = blkdiag(diag(alpha1), diag(alpha2));
Beta  = beta*eye(nphi);
Psi_phi = [Beta, -eye(nphi); -Alpha, eye(nphi)];

N   = blkdiag(W1,W2,W3);
Nux = N(n1+n2+1:end, 1:nG);   Nuw = N(n1+n2+1:end, nG+1:end);
Nvx = N(1:n1+n2, 1:nG);       Nvw = N(1:n1+n2, nG+1:end);

%% Delay IQC filter: r = q - p_{k-d}
n_p = nq;  npsi = d*n_p;
J = diag(ones(d-1,1), -1);
Apsi = kron(J, eye(n_p));
Bpsi1 = [eye(n_p); zeros((d-1)*n_p, n_p)];  Bpsi2 = zeros(npsi, n_p);
Cpsi  = [zeros(n_p, (d-1)*n_p), -eye(n_p)];
Dpsi1 = zeros(n_p, n_p);  Dpsi2 = eye(n_p);
M_delay = eye(n_p);

%% Augmented plant o delay-filter
Abar = [AG, zeros(nG,npsi); Bpsi1*CG, Apsi];
Bbar = [BG1, BG2; Bpsi1*DG1+Bpsi2, Bpsi1*DG2];
Cbar = [Dpsi1*CG, Cpsi];
Dbar = [Dpsi1*DG1+Dpsi2, Dpsi1*DG2];
nzeta = nG + npsi;

%% Interval LFT: per-entry boxes for A0, A1
[idx_r, idx_c] = ndgrid(1:nG, 1:nG);  idx_r = idx_r(:); idx_c = idx_c(:);
s0 = nG*nG; s1 = nG*nG; s_theta = s0 + s1;
rho0d = dt*R0(:);  rho1d = dt*R1(:);  rho_d = [rho0d; rho1d];
E0 = zeros(nG,s0); for j=1:s0, E0(idx_r(j),j)=1; end
E1 = zeros(nG,s1); for j=1:s1, E1(idx_r(j),j)=1; end
Etheta = [E0, E1];
Bbar_ext = [Bbar, [Etheta; zeros(npsi, s_theta)]];

%% Selection operators (master stack S_vec = [zeta; w_phi; q; w_theta])
Rzeta   = [eye(nzeta), zeros(nzeta, nphi+nq+s_theta)];
Rq      = [zeros(nq, nzeta+nphi), eye(nq), zeros(nq, s_theta)];
Rwtheta = [zeros(s_theta, nzeta+nphi+nq), eye(s_theta)];
Ru = [Nux, zeros(nu,npsi), Nuw, zeros(nu,nq), zeros(nu,s_theta)];
RV = [Rzeta; Rq; Ru; Rwtheta];
Nvzeta = [Nvx, zeros(nphi,npsi)];
Rphi = [Nvzeta, Nvw, zeros(nphi, nq+s_theta);
        zeros(nphi, nzeta), eye(nphi), zeros(nphi, nq+s_theta)];
stack_dim = nzeta + nphi + nq + s_theta;
Zsel = zeros(s_theta, stack_dim);  Wsel = zeros(s_theta, stack_dim);
for j = 1:s0
    c = idx_c(j);
    Zsel(j,c) = 1;  Wsel(j, nzeta+nphi+nq+j) = 1;
end
for j = 1:s1
    c = idx_c(j);  row = s0 + j;
    Zsel(row, nzeta+nphi+c) = 1;  Wsel(row, nzeta+nphi+nq+row) = 1;
end

%% CVX: main LMI
cvx_begin sdp quiet
    cvx_solver sedumi
    variable P(nzeta, nzeta) symmetric
    variable Tphi(nphi, nphi) diagonal
    variable lambda_tau
    variable mu_theta(s_theta)
    P >= 1e-8*eye(nzeta);
    Tphi(:) >= 0;  lambda_tau >= 0;  mu_theta >= 0;

    Mphi = [zeros(nphi), Tphi; Tphi, zeros(nphi)];
    Qphi = Rphi' * Psi_phi' * Mphi * Psi_phi * Rphi;

    S11 = Abar'*P*Abar - P;   S12 = Abar'*P*Bbar_ext;   S22 = Bbar_ext'*P*Bbar_ext;
    S   = [S11, S12; S12', S22];
    QV  = RV' * S * RV;

    Ddelay_ext = [Dbar, zeros(size(Dbar,1), s_theta)];
    Qdelay = lambda_tau * RV' * [Cbar, Ddelay_ext]' * M_delay * [Cbar, Ddelay_ext] * RV;

    Qtheta = 0;
    for j = 1:s_theta
        Mbox_j = diag([rho_d(j)^2, -1]);
        Rj = [Zsel(j,:); Wsel(j,:)];
        Qtheta = Qtheta + mu_theta(j) * (Rj' * Mbox_j * Rj);
    end

    Qtotal = QV + Qphi + Qdelay + Qtheta;
    Qtotal <= -1e-10*eye(stack_dim);

    for i = 1:n1
        [deltav1^2, [W1(i,:) zeros(1,npsi)]; [W1(i,:) zeros(1,npsi)]', P] >= 0;
    end
    minimize( trace(P(1:nG,1:nG)) )
cvx_end

status = cvx_status;
info.P = P;
info.feasible = ~isempty(strfind(lower(status),'solved'));
Px = P(1:nG,1:nG);
if all(isfinite(Px(:))) && all(eig(Px) > 0)
    radii = 1./sqrt(eig(Px));
    radii = radii(:).';
else
    radii = [];
end
end
