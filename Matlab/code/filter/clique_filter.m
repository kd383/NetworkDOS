% Q = clique_filter(A, nclique)
% Create filters for cliques of specific size. A clique of size nnode
% produces an eigenvalue of -1/nnode.
% 
% Input:
%   A: The adjacency matrix.
%   nclique: Size of the clique.
%
% Output:
%   Q: Nnode x Nfil filter matrix (eigenvectors of normlized adjacency).
%
function Q = clique_filter(A, nclique)

vinfo = vertex_hashing(A, 2);
deg = sum(A);
sqd = sqrt(deg);
cdeg = zeros(size(vinfo, 1), 1);
for i = 1:length(vinfo)
    cdeg(i) = deg(vinfo{i}(1));
end

n = size(A,1);
Q = filter_construct(vinfo(cdeg==nclique), sqd, n);