clear all
close all
clc

graph_name = 'pgp-cc'; % The test graph name.
A = load_graph(graph_name);
N = matrix_normalize(A);

% Create the nested dissection.
[perm, tree] = ndmetis(A, graph_name);
Np = N(perm, perm);
assert(CheckValidTree(tree, Np),'The nested dissection tree is not valid.')

n = length(perm);
T{1} = MatrixPtr(eye(n));
T{2} = MatrixPtr(Np);
Tbase = T;

% Naive three-term recurrence as baseline.
tic;
for i = 3:10
    Tbase{i} = MatrixPtr(2*Np*Tbase{i-1}.data - Tbase{i-2}.data);
end
time_base = toc;

% Nested dissection three-term recurrence.
tic;
for i = 3:10
    T{i} = MatrixPtr(zeros(n));
    BuildRootColumns(tree, Np, T{i}, T{i-1}, T{i-2});
    BuildLeafBlocks(tree, Np, T{i}, T{i-1}, T{i-2});
end
time_nd = toc;

err = 0;
for i = 1:10
    err = max(err, norm(diag(Tbase{i}) - diag(T{i})));
end

fprintf('---------------------------\n');
fprintf('Diagonal error is %.6f.\n', err);
fprintf('---------------------------\n');
fprintf('Base method takes %.6f seconds.\n', time_base);
fprintf('Nested Dissection method takes %.6f seconds.\n', time_nd);