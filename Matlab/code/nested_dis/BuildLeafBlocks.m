% Three-term recurrence for the leaf blocks. Each block gets corrections
% from separators along the path to the root.
function BuildLeafBlocks(tree, A, Ti, Tj, Tk)
if ~isempty(tree.left)
    BuildLeafBlocks(tree.left, A, Ti, Tj, Tk);
end
if ~isempty(tree.right)
    BuildLeafBlocks(tree.right, A, Ti, Tj, Tk);
end
if isempty(tree.left) && isempty(tree.right)
    Ti.data(tree.index, tree.index) = ...
        2*A(tree.index, tree.index)*Tj.data(tree.index, tree.index) ...
        - Tk.data(tree.index, tree.index);
    parent = tree.parent;
    while ~isempty(parent)
        nPZ = max(min(length(parent.root),1), ceil(length(parent.root)));
        Ti.data(tree.index, tree.index) = Ti.data(tree.index, tree.index) + ...
            2*A(tree.index, parent.root(1:nPZ))*Tj.data(tree.index, parent.root(1:nPZ))';
        parent = parent.parent;
    end  
end
end