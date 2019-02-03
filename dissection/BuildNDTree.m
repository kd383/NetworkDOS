% DONT USE: Heuristic for retrieving a nested dissection tree from its
% ordering.
function tree = BuildNDTree(A, index, nLevel, parent)
if nargin < 4,  parent = []; end
if nargin < 3,  nLevel = 5; end
if nargin < 2,  index = 1:length(A); end
tree = TreeNode(index, parent);
if nLevel == 0 || ~nnz(A)
    return
end
n = length(index);
for i = ceil(0.9*n):n
    if find(A(:,i), 1, 'first') < n/2.2
        break;
    end
end 
j = find(sum(A(1:ceil(0.45*(i-1)), 1:i-1)), 1, 'last');
k = find(sum(A(1:j, j+1:i-1)), 1, 'first');
if ~isempty(k), i = k;  end
tree.root = index(i:end);
tree.left = BuildNDTree(A(1:j, 1:j), index(1:j), nLevel-1, tree);
tree.right = BuildNDTree(A(j+1:i-1, j+1:i-1), index(j+1:i-1), nLevel-1, tree);
end