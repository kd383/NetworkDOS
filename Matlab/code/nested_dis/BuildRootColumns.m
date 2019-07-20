% Three-term recurrence for the separators. Each separator gets corrections
% from separators above along the path to the root.
function BuildRootColumns(tree, A, Ti, Tj, Tk)
nRoot = length(tree.root);
nZ = max(ceil(nRoot), min(nRoot, 1));
Ti.data(tree.index, tree.root(1:nZ)) = ...
    2*A(tree.index,tree.index)*Tj.data(tree.index,tree.root(1:nZ)) ...
    - Tk.data(tree.index,tree.root(1:nZ));
parent = tree.parent;
while ~isempty(parent)
    nPZ = max(min(length(parent.root), 1), ceil(length(parent.root)));
    Ti.data(tree.index, tree.root(1:nZ)) = Ti.data(tree.index, tree.root(1:nZ)) + ...
        2*A(tree.index, parent.root(1:nPZ))*Tj.data(tree.root(1:nZ), parent.root(1:nPZ))';
    parent = parent.parent;
end
if ~isempty(tree.left)
    BuildRootColumns(tree.left, A, Ti, Tj, Tk);
end
if ~isempty(tree.right)
    BuildRootColumns(tree.right, A, Ti, Tj, Tk);
end
end