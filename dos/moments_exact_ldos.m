% Local density of states using exact three-term recurrence.
%
% Input:
%   A: Graph matrix.
%   N: Number of moments to compute.
%   
% Output:
%   c: Local density of states up to N-th moment.

function c = moments_exact_ldos(A, N)
n = length(A);
c = zeros(N, n);
Tk = eye(n);
Tj = A;
c(1, :) = 1;
c(2, :) = diag(A)';
for i = 3:N
    Ti = 2*A*Tj - Tk;
    c(i, :) = diag(Ti)';
    Tk = Tj;
    Tj = Ti;
end
end