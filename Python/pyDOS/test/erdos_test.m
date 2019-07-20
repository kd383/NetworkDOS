%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is a testing script for pyDOS on Erdos matrix
%
%   Author: Kun Dong, kd383@cornell.edu
%
%   11/28/2016 - Created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the path for graph-dos repository
% change to local directory
GRAPH_DOS_PATH = '../../../graph-dos/src';
addpath(GRAPH_DOS_PATH)

% load the Erdos Matrix
MA = load('../data/erdos02-cc.mat');
MA = MA.A;

% Run erdos_test.py before proceeding
if isempty(dir('*.mat'))
    system('python erdos_test.py');
end
flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check matrix_generation.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_mat_gen.mat');
error = zeros(4,1);
error(1) = norm(py.A1-MA,'fro');
error(2) = norm(py.A2-matrix_laplacian(MA),'fro');
error(3) = norm(py.A3-matrix_normalize(MA),'fro');
error(4) = norm(py.A4-matrix_slaplacian(MA),'fro');
if max(error) < 1e-8
    fprintf('matrix_generation.py working!\n');
else
    flag = 1;
end
clear('py','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check mfunc_generation.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_mfc_gen.mat');
error = zeros(5,1);
Mf.M1 = mfunc_laplacian(MA);
error(1) = norm(py.Afun1x - Mf.M1(py.x),'fro');
Mf.M2 = mfunc_modularity(MA);
error(2) = norm(py.Afun2x - Mf.M2(py.x),'fro');
Mf.M3 = mfunc_normalize(MA);
error(3) = norm(py.Afun3x - Mf.M3(py.x),'fro');
Mf.M4 = mfunc_slaplacian(MA);
error(4) = norm(py.Afun4x - Mf.M4(py.x),'fro');
Mf.M5 = mfunc_cheb_poly(py.c,MA);
error(5) = norm(py.Afun5x - Mf.M5(py.x),'fro');
if max(error) < 1e-8
    fprintf('mfunc_generation.py working!\n');
else
    flag = 1;
end
clear('py','Mf','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check eig_rand1.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_eig_rand.mat');
error = zeros(3,1);
M.Afun = @(x)py.U*(py.Acore*(py.U'*x));
[M.V,M.D] = eig_rand1(py.Z,M.Afun);
M.V = bsxfun(@times,M.V,sign(M.V(1,:).*py.V(1,:)));
error(1) = norm(py.V - M.V,'fro');
error(2) = norm(diag(py.D) - M.D,'fro');
error(3) = norm(M.Afun(py.V)-py.V*diag(py.D));
if max(error) < 1e-8
    fprintf('eig_rand1.py working!\n');
else
    flag = 1;
end
clear('py','M','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check moment_filter.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_mmt_fil.mat');
error = zeros(2,1);
M.cj = filter_jackson(py.c);
M.cl = filter_lorentz(py.c);
error(1) = norm(py.cj - M.cj,'fro');
error(2) = norm(py.cl - M.cl,'fro');
if max(error) < 1e-8
    fprintf('moment_filter.py working!\n');
else
    flag = 1;
end
clear('py','M','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check rescale_matrix.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_rsc_mat.mat');
error = zeros(4,1);
[M.As,M.ab] = rescale_matrix(MA);
[M.Ls,M.cd] = rescale_mfunc(matrix_normalize(MA));
error(1) = norm(py.As - M.As,'fro');
error(2) = norm(py.ab - M.ab,'fro');
error(3) = norm(py.Lx - M.Ls(py.x),'fro');
error(4) = norm(py.cd - M.cd,'fro');
if max(error) < 1e-8
    fprintf('rescale_matrix.py working!\n');
else
    flag = 1;
end
clear('py','M','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check moment_comp.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_mmt_cmp.mat');
error = zeros(9,1);
M.L = matrix_normalize(MA);
M.Lfun = mfunc_normalize(MA);
M.n = length(M.L);
[M.c,M.cs] = moments_cheb_dos(M.L,py.Z,20);
error(1:2) = [norm(py.c - M.c,'fro');norm(py.cs - M.cs,'fro')];
[M.cl,M.csl] = moments_cheb_ldos(M.L,py.Z,20);
error(3:4) = [norm(py.cl - M.cl,'fro');norm(py.csl - M.csl,'fro')];
[M.cfun,M.csfun] = moments_cheb_dos(M.Lfun,M.n,py.Z,20);
error(5:6) = [norm(py.cfun - M.cfun,'fro');norm(py.csfun - M.csfun,'fro')];
[M.cfunl,M.csfunl] = moments_cheb_ldos(M.Lfun,M.n,py.Z,20);
error(7:8) = [norm(py.cfunl - M.cfunl,'fro');norm(py.csfunl - M.csfunl,'fro')];
M.cd = moments_delta(py.x,20);
error(9) = norm((py.cd- M.cd)./py.cd,'fro');
if max(error) < 1e-8
    fprintf('moment_comp.py working!\n');
else
    flag = 1;
end
clear('py','M','error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check plot_cheb.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py = load('test_plt_chb.mat');
error = zeros(6,1);
M.y = plot_cheb(py.c);
error(1) = norm(py.y - M.y,'fro');
[M.yl,M.idx] = plot_cheb_ldos(py.cl);
error(2) = norm(py.yl - M.yl,'fro');
error(3) = min(norm(M.idx-1-py.idx,'fro'),norm(M.idx-1-flipud(py.idx),'fro'));
M.yh = plot_chebhist(py.c);
error(4) = norm(py.yh - M.yh,'fro');
M.yi = plot_chebint(py.c);
error(5) = norm(py.yi - M.yi,'fro');
M.yp = plot_chebp(py.c);
error(6) = norm(py.yp - M.yp,'fro');
if max(error) < 1e-8
    fprintf('plot_cheb.py working!\n');
else
    flag = 1;
end
clear('py','M','error');


if flag 
    fprintf('Something is wrong!\n');
else
    fprintf('Hooray!\n');
    delete('*.mat');
end