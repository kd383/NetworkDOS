% Helper function that convert an adjacency matrix into a METIS graph file.
%
% Input:
%   A: Adjacency matrix.
%   graph_path: The directory for the output graph.

function write_metis_graph(A, graph_path)
if nargin < 2, graph_path = './test.graph'; end
f = fopen(graph_path, 'w');

% IMPORTANT: Self-edges in graph is ignored because they don't affect
% nested dissection.
n = length(A);
A(1:n+1:end) = 0;
fprintf(f, '%d %d\n', n, nnz(A)/2);
for i = 1:length(A)
    neighbor = find(A(:,i));
    for j = 1:length(neighbor)-1
        fprintf(f, '%d ', neighbor(j));
    end
    fprintf(f, '%d', neighbor(end));
    if i < length(A), fprintf(f, '\n'); end
end
fclose(f);
end