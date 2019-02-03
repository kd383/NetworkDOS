% Density of states using exact three-term recurrence.
%
% Input:
%   A: Graph matrix.
%   N: Number of moments to compute.
%   
% Output:
%   c: Density of states up to N-th moment.

function c = moments_exact_dos(A, N)
n = length(A);
c = zeros(N, 1);
Tk = eye(n);
Tj = A;
c(1) = n;
c(2) = trace(A);
for i = 3:N
    Ti = 2*A*Tj - Tk;
    c(i) = trace(Ti);
    Tk = Tj;
    Tj = Ti;
end
end