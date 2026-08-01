function [pass, info] = certify_positive_lure(A_lo, A_up, B, C, gamma1, gamma2, tol)
%CERTIFY_POSITIVE_LURE  Delay-independent positivity certificate (Thm 2-4).
%   [pass, info] = CERTIFY_POSITIVE_LURE(A_lo, A_up, B, C, gamma1, gamma2)
%   checks the linear, delay-INDEPENDENT stability certificate for an
%   NN-in-the-loop Lur'e system whose controller is sector bounded in
%   [gamma1, gamma2] (from SECTOR_BOUND_FFNN).
%
%   Inputs
%     A_lo  - sum of the LOWER interval endpoint matrices, A0 + sum_i Ai
%     A_up  - sum of the UPPER interval endpoint matrices, A0bar + sum_i Aibar
%     B, C  - input / output matrices (assumed elementwise nonnegative)
%     gamma1, gamma2 - lower / upper network sector slopes
%     tol   - numerical tolerance (default 1e-9)
%
%   Certificate (holds for EVERY delay tau >= 0 and every admissible interval):
%     (Metzler, lower sector)  M = A_lo + B*gamma1*C  has nonnegative off-diagonals
%     (Hurwitz, upper sector)  H = A_up + B*gamma2*C  has all eigenvalues Re < 0
%
%   For the interval-only case (C1) this is exactly Theorem 2. For the delay
%   cases (C2/C3) positivity additionally requires A0 Metzler and each delayed
%   block Ai + B*gamma1*C >= 0, which hold by construction for the shipped
%   examples; the operative pass/fail test is the Metzler + Hurwitz pair below.
%
%   Output INFO contains the margins (metzler_margin, max_real_eig), the boolean
%   sub-tests, the eigenvalues of H, and the matrices M and H.
if nargin < 7 || isempty(tol), tol = 1e-9; end

M = A_lo + B*gamma1*C;
H = A_up + B*gamma2*C;

Moff = M;
Moff(1:size(M,1)+1:end) = +Inf;                      % ignore the diagonal
info.metzler_margin = min(Moff(:));                   % true min off-diagonal; >=0 => Metzler
info.metzler_ok     = info.metzler_margin >= -tol;

info.eigs         = eig(H);
info.max_real_eig = max(real(info.eigs));            % <  0  => Hurwitz
info.hurwitz_ok   = info.max_real_eig < -tol;

info.M = M;  info.H = H;
pass = info.metzler_ok && info.hurwitz_ok;
end
