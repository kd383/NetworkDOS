% demo_ldos(method, dname, Nprobe, Nmoment, Nbin)
%
% Loads the graph in the rodger collection (with reference eigenvalues)
% and plots histograms computed with true eigenvalue distribution and with
% the  estimate.
%
% Input:
%   method: One of the approximation methods.
%           Possible input value: 'cheb'(default), 'lan', 'exact', 'nd'.
%   dname:  Data set name, or an adjacency matrix A.
%   Nprobe: Number of probe vectors for moment estimation. Default: 20.
%   Nmoment:  Number of Chebyshev moments. Default: 1000.
%   Nbin:   Number of histogram bins. Default: 50.
%
function demo_ldos(method, dname, Nprobe, Nmoment, Nbin)

if nargin < 1, method = 'cheb'; end
if nargin < 2, dname = 'Erdos02-cc'; end
if nargin < 3, Nprobe = 20;  end
if nargin < 4, Nmoment = 1000; end
if nargin < 5, Nbin  = 50;   end

% Load data and precompute.
tic;
if ischar(dname)
    [A,lambda] = load_graph(dname);
else
    A = dname;
    if length(A) > 1e4
        warning('A might be too large for exact computation.');
    end
end

  
N = matrix_normalize(A, 's');
fprintf('Time to load and convert: %g\n', toc);

% Compute spectral density as Chebyshev moments.
tic;
switch method
    case 'cheb'
        cl = moments_cheb_ldos(N, Nprobe, Nmoment);
    case 'lan'
        node_ind = ceil(linspace(1, length(N),100));
        cl = moments_lan_ldos(N, node_ind, Nprobe, Nmoment);
    case 'exact'
        cl = moments_exact_ldos(N, Nmoment);
    case 'nd'
        cl = moments_nd_ldos(N, Nmoment);
    otherwise
        error('Method unrecognized.');
end

% Apply Jackson filter to remove Gibbs phenomenon.
cf = filter_jackson(cl);
fprintf('Time to compute Chebyshev moments: %g\n', toc);

% Plot
plot_cheb_ldos(Nbin, cf);