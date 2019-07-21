% Q = zero_filter(dname,option)
% Construct the filter for node duplicate which leads to eigenvalue 0.
% It is the most common motif.
%
% Input:
%   A: The adjacency matrix.
%
% Output:
%   Q: Nnode x Nfil filter matrix (eigenvectors/null vectors)
%
function Q = zero_filter(A)

vinfo = vertex_hashing(A, 1);
sqd = sqrt(sum(A));
n = size(A,1);
Q = filter_construct(vinfo, sqd, n);