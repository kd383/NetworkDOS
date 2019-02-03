% A tree node structure to represent the nested dissection.
classdef TreeNode < handle
    properties
        index % Nodes in this partition.
        parent % Tree node at previous level.
        root % Separator index.
        left % Left partition.
        right % Right partition.
    end
    methods
        function obj = TreeNode(index, parent)
            if nargin < 1,  return; end
            if nargin < 2, parent = []; end
            obj.index = index;
            obj.parent = parent;
        end
        
        % Helper function check if tree represents a valid nested
        % dissection.
        function valid = CheckValidTree(obj, A)
            valid_left = isempty(obj.left) || CheckValidTree(obj.left, A);
            valid_right = isempty(obj.right) || CheckValidTree(obj.right, A);
            valid_mix = isempty(obj.left) || isempty(obj.right) ||  ...
                ~nnz(A(obj.left.index, obj.right.index));
            valid = valid_left && valid_right && valid_mix;
        end
    end
end

