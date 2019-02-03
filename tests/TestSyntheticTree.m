nLevel = 4; % Consists of (nLevel-1) levels of separator and one level of leaf blocks
nLeafNodes = 250;
nRootNodes = 50;
n = nLeafNodes * 2^(nLevel-1) + nRootNodes * (2^(nLevel-1)-1);
index = 1:n;
% Build the nested dissection index tree.
tree = BuildTree(index, nLevel-1, nRootNodes,[]);
% Create a matrix based on the tree.
A = BuildTreeGraph(tree);
sqrt_degree = sqrt(sum(A));
% Normalized it to have spectral radius 1.
NA = A ./ sqrt_degree ./ sqrt_degree';

T{1} = MatrixPtr(eye(n));
T{2} = MatrixPtr(NA);
Ttest = T;

% Base method of straight-forward recursion.
tic;
for i = 3:10
    Ttest{i} = MatrixPtr(2*NA*Ttest{i-1}.data - Ttest{i-2}.data);
end
Tbase = toc;

% Nested Dissection Method.
% The computation of Ti(NA) is broken into two steps.
% 1. We calculuate the (incomplete) columns for the root indices (separators).
% 2. We fill in the leaf blocks with three-term recurrence plus separator
%    correction.
tic;
for i = 3:10
    T{i} = MatrixPtr(zeros(n));
    BuildRootColumns(tree, NA, T{i}, T{i-1}, T{i-2});
    BuildLeafBlocks(tree, NA, T{i}, T{i-1}, T{i-2});
end
Tnd = toc;

error = 0;
for i = 1:10
    error = max(error, norm(diag(Ttest{i}) - diag(T{i})));
end

fprintf('---------------------------\n');
fprintf('Diagonal error is %.6f.\n', error);
fprintf('---------------------------\n');
fprintf('Base method takes %.6f seconds.\n', Tbase);
fprintf('Nested Dissection method takes %.6f seconds.\n', Tnd);


% Auxiliary function for building a tree of indicies.
% Each level has @nRootNodes separators, then the rest of the indicies are
% divided evenly between left and right sub-trees.
function tree = BuildTree(index, nLevel, nRootNodes, parent)
if nLevel == 0
    tree = TreeNode();
    tree.index = index;
    tree.parent = parent;
    return
end
tree = TreeNode();
tree.index = index;
right_end = length(index) - nRootNodes;
left_end = right_end/2;
tree.parent = parent;
tree.root = index(right_end+1:end);
tree.left = BuildTree(index(1:left_end), nLevel-1, nRootNodes, tree);
tree.right = BuildTree(index(left_end+1:right_end), nLevel-1, nRootNodes, tree);
end

% Build a graph(matrix) that satisfies the given nested dissection tree.
function graph = BuildTreeGraph(tree)
n = length(tree.index);
if isempty(tree.left) && isempty(tree.right)
    graph = rand(n);
    return
end
graph = zeros(n);
nR = length(tree.root);
nl = length(tree.left.index);
nr = length(tree.right.index);
graph(:, end-nR+1:end) = rand(n, nR);
graph(1:nl, 1:nl) = BuildTreeGraph(tree.left);
graph(nl+1:nl+nr, nl+1:nl+nr) = BuildTreeGraph(tree.right);
if isempty(tree.parent)
    graph(graph>0.1) = 0;
    graph = (graph + graph')/2;
    graph = sparse(graph);
end
end