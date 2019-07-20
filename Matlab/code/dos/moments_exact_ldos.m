% Local density of states using exact three-term recurrence.
%
% Input:
%   A: Graph matrix.
%   N: Number of moments to compute.
%   
% Output:
%   c: Local density of states up to N-th moment.
%
function c = moments_exact_ldos(varargin)

defaults = {'Afun', NaN, 'n', NaN, 'N', 10};
[Afun, n, N] = mfuncify(defaults, varargin{:});

c = zeros(N, n);
Tk = eye(n);
Tj = Afun(eye(n));

c(1, :) = 1;
c(2, :) = diag(Tj)';
for i = 3:N
    Ti = 2 * Afun(Tj) - Tk;
    c(i, :) = diag(Ti)';
    Tk = Tj;
    Tj = Ti;
end