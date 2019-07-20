% This is a test script for a comparison between all our methods on rodger
% dataset
clear all
close all
clc

% Shortlist: Use this if you want a quick test.
% graph_names = {'as19991115', 'Erdos02-cc', 'yeast-cc'};

graph_names = {
    'as19991115', ...
    'as-caida20060911', ...
    'Erdos02-cc', ...
    'homo-cc', ...
    'marvel-chars-cc', ...
    'marvel-comics-cc', ...
    'musm-cc', ...
    'pgp-cc', ...
    'yeast-cc', ...
    };

path = fileparts(mfilename('fullpath'));
f_dos = fopen([path '/dos_result.txt'], 'w');
f_ldos = fopen([path '/ldos_result.txt'], 'w');
fprintf(f_dos, '%15s%9s%20s%20s%20s\n', 'Graph Name', 'Exact', ...
    'Nest Dissection','Chebyshev','Lanczos');
fprintf(f_ldos, '%15s%9s%20s%20s\n', 'Graph Name', 'Exact', 'Nest Dissection',...
    'Chebyshev');
for i = 1:length(graph_names)
    graph_name = graph_names{i};
    fprintf(f_dos, '%15s', graph_name);
    fprintf(f_ldos, '%15s', graph_name);
    A = load_graph(graph_name);
    N = matrix_normalize(A);
    n = length(A);
    [perm, tree] = ndmetis(A, graph_name);
    Np = N(perm, perm);
    assert(CheckValidTree(tree, Np),'The nested dissection tree is not valid.');

    tic;
    c_exact = moments_exact_dos(N, 20);
    time = toc;
    fprintf(f_dos, '%8.3fs', time);

    tic;
    c_nd = moments_nd_dos(N, 20, perm, tree);
    time = toc;
    err = norm(c_nd - c_exact) / norm(c_exact);
    fprintf(f_dos, '%8.3fs %10.3e', time, err);

    tic;
    c_cheb = moments_cheb_dos(N, ceil(n/4), 20);
    time = toc;
    err = norm(c_cheb - c_exact) / norm(c_exact);
    fprintf(f_dos, '%8.3fs %10.3e', time, err);

    tic;
    c_lanz = moments_lan_dos(N, ceil(n/4), 20);
    time = toc;
    err = norm(c_lanz - c_exact) / norm(c_exact);
    fprintf(f_dos, '%8.3fs %10.3e', time, err);
    
    tic;
    c_exact = moments_exact_ldos(N, 20);
    time = toc;
    fprintf(f_ldos, '%8.3fs', time);
    
    tic;
    c_nd = moments_nd_ldos(N, 20, perm, tree);
    time = toc;
    err = mean(sqrt(sum((c_nd-c_exact).^2)./sum(c_exact.^2)));
    fprintf(f_ldos, '%8.3fs %10.3e', time, err);
    
    tic;
    c_cheb = moments_cheb_ldos(N, ceil(n/4), 20);
    time = toc;
    err = mean(sqrt(sum((c_cheb-c_exact).^2)./sum(c_exact.^2)));
    fprintf(f_ldos, '%8.3fs %10.3e', time, err);
    
    if i~=length(graph_names)
        fprintf(f_dos, '\n');
        fprintf(f_ldos, '\n');
    end
end
fclose(f_dos);
fclose(f_ldos);