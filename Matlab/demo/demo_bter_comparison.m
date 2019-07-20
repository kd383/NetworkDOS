clear all;
close all;
clc

graph_name = 'Erdos02-cc';
A = load_graph(graph_name);
% Generate the Bter graph
A_bter = bter_model(A);
N = matrix_normalize(A);
N_bter = matrix_normalize(A_bter);

demo_dos('cheb', A);
axis([-1 1 0 500]);
demo_dos('cheb', A_bter);
axis([-1 1 0 500]);
demo_ldos('cheb', A);
demo_ldos('cheb', A_bter);


% Bter graph construction
% Check the instruction at the link:
% http://www.sandia.gov/~tgkolda/feastpack/doc_bter_match.html
function A_bter = bter_model(A)
nd = accumarray(nonzeros(sum(A, 2)), 1);
[ccd, ~] = ccperdeg(A);
[E1,E2] = bter(nd,ccd);
A_bter = bter_edges2graph(E1,E2);
end