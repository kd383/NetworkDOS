clear all
close all
clc

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

f_dos = fopen('dos_result.txt', 'w');
f_ldos = fopen('ldos_result.txt', 'w');
for i = 1:length(graph_names)
    graph_name = graph_names{i};
    fprintf(f_dos, [graph_name ' ']);
    fprintf(f_ldos, [graph_name ' ']);
    A = load_graph(graph_name);
    N = matrix_normalize(A);
    n = length(A);
    [perm, tree] = ndmetis(A, graph_name);
    Np = N(perm, perm);
    assert(CheckValidTree(tree, Np),'The nested dissection tree is not valid.');

    tic;
    c_exact = moments_exact_dos(N, 20);
    time = toc;
    fprintf(f_dos, '%.3f ', time);

    tic;
    c_nd = moments_nd_dos(N, 20, perm, tree);
    time = toc;
    err = norm(c_nd - c_exact) / norm(c_exact);
    fprintf(f_dos, '%.3f %.3e ', time, err);

    tic;
    c_cheb = moments_cheb_dos(N, ceil(n/4), 20);
    time = toc;
    err = norm(c_cheb - c_exact) / norm(c_exact);
    fprintf(f_dos, '%.3f %.3e ', time, err);

    tic;
    c_lanz = moments_lan_dos(N, ceil(n/4), 20);
    time = toc;
    err = norm(c_lanz - c_exact) / norm(c_exact);
    fprintf(f_dos, '%.3f %.3e ', time, err);
    
    tic;
    c_exact = moments_exact_ldos(N, 20);
    time = toc;
    fprintf(f_ldos, '%.3f ', time);
    
    tic;
    c_nd = moments_nd_ldos(N, 20, perm, tree);
    time = toc;
    err = mean(sqrt(sum((c_nd-c_exact).^2)./sum(c_exact.^2)));
    fprintf(f_ldos, '%.3f %.3e ', time, err);
    
    tic;
    c_cheb = moments_cheb_ldos(N, ceil(n/4), 20);
    time = toc;
    err = mean(sqrt(sum((c_cheb-c_exact).^2)./sum(c_exact.^2)));
    fprintf(f_ldos, '%.3f %.3e ', time, err);
    
    if i~=length(graph_names)
        fprintf(f_dos, '\n');
        fprintf(f_ldos, '\n');
    end
end
fclose(f_dos);
fclose(f_ldos);
