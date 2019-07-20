% Density of states using exact three-term recurrence.
%
% Input:
%   A: Graph matrix.
%   N: Number of moments to compute.
%   
% Output:
%   c: Density of states up to N-th moment.
%
function c = moments_exact_dos(varargin)

defaults = {'Afun', NaN, 'n', NaN, 'N', 10};
[Afun, n, N] = mfuncify(defaults, varargin{:});

c = zeros(N, 1);
Tk = eye(n);
Tj = Afun(eye(n));
c(1) = n;
c(2) = trace(Tj);
for i = 3:N
    Ti = 2 * Afun(Tj) - Tk;
    c(i) = trace(Ti);
    Tk = Tj;
    Tj = Ti;
end