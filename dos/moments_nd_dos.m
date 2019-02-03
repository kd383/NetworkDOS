% Density of states using nested dissection method.
%
% Input:
%   A: Graph matrix.
%   N: Number of moments to compute.
%   perm: Nested dissection ordering.
%   tree: Nested dissection partition tree.
%   
% Output:
%   c: Density of states up to N-th moment.

function c = moments_nd_dos(A, N, perm, tree)
if nargin < 4
    [perm, tree] = ndmetis(A);
end
Ap = A(perm, perm);
assert(CheckValidTree(tree, Ap),'The nested dissection tree is not valid.');

n = length(A);
c = zeros(N, 1);
Tk = MatrixPtr(speye(n));
Tj = MatrixPtr(Ap);
c(1) = n;
c(2) = trace(Tj);
for i = 3:N
    Ti = MatrixPtr(zeros(n));
    BuildRootColumns(tree, Ap, Ti, Tj, Tk);
    BuildLeafBlocks(tree, Ap, Ti, Tj, Tk);
    c(i) = trace(Ti);
    Tk = Tj;
    Tj = Ti;
end
end