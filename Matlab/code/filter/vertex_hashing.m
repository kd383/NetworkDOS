% vinfo = vertex_hashing(A)
% Use random hashing on the neighbor list of nodes. Find the ones
% satisfying certin motif patterns.
%
% Input:
%   G: Graph matrix, usually (normalized) adjacency.
%   pat: Hashing pattern. 1: x->G*x. 2: x->(G+I)x.
%
% Output:
%   vinfo: Sets of vertices that satisfies the given patter.
%
% Note: This function currently finds all nodes with identical hash weight.
% It may get updated in the feature in order to caputer other patterns,
% such as the sum of two weights equal to a third.
%
function vinfo = vertex_hashing(G, pat)

% Determine hashing function
switch pat
    case 1
        hashfun = @(x) G*x;
    case 2
        hashfun = @(x) G*x + x;
    otherwise
        error('Pattern %d not yet defined.', pat);
end

% Hashing with random vector.
n = size(G, 1);
w = hashfun(randn(n,1));
[~, ~, IC] = uniquetol(w);
k = max(IC);
vinfo = cell(k, 1);
kt = 0;

% Find all sets of nodes with idential hashing weights
for i = 1:k
    idx = find(IC == i);
    if length(idx) > 1
        kt = kt+1;
        vinfo{kt} = idx;
    end
end
vinfo(kt+1:end) = [];