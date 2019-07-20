% A temporary method for calling ndmetis.c routine.
% It will be replaced by a mex function in the future.
%
% Input:
%   A: The adjacency matrix.
%   graphname: A string for the name of the graph.
%
% Output:
%   perm: The permutation vector for nested dissection ordering with
%         separators at the end.
%   tree: A tree that holds the separators and left/right partitions.

function [perm, tree] = ndmetis(A, graphname)
if nargin < 2, graphname = 'test'; end
% Create the path of input and output.
path = fileparts(mfilename('fullpath'));
graph_path = [path, '/', graphname, '.graph'];
% Write the adjacency matrix into METIS format.
write_metis_graph(A, graph_path);
cmd_ndmetis = [path, '/ndmetis ', graph_path];
cmd_cleanup = ['rm ', graph_path, ' ', graph_path, '.iperm ', './ndtree.txt'];

% Call ndmetis.
[status, result] = system(cmd_ndmetis);

% Load the permutation vector.
iperm = load([graph_path, '.iperm']) + 1;
perm(iperm) = 1:length(iperm);

% Find the separators and create the tree.
tree = TreeNode;
f = fopen('./ndtree.txt', 'r');
s = fgetl(f);
while ischar(s)
    if contains(s, 'Index')
        tree.index = iperm(str2num(s(7:end)) + 1);
    elseif contains(s, 'Root')
        tree.root = iperm(str2num(s(6:end)) + 1);
    elseif contains(s, 'Left')
        tree.left = TreeNode([], tree);
        tree = tree.left;
    elseif contains(s, 'Right')
        tree.right = TreeNode([], tree);
        tree = tree.right;
    elseif contains(s, 'end')
        tree = tree.parent;
    else
        fprintf('Unexpected line.\n');
    end
    s = fgetl(f);
end
fclose(f);
system(cmd_cleanup);
end