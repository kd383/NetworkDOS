% demo_dos(method, dname, Nprobe, Nmoment, Nbin)
%
% Loads the graph in the rodger collection (with reference eigenvalues)
% and plots histograms computed with true eigenvalue distribution and with
% the estimate.
%
% Input:
%   method: One of the approximation methods.
%           Possible input value: 'cheb'(default), 'lan', 'exact', 'nd'.
%   dname:  Data set name, or an adjacency matrix A.
%   Nprobe: Number of probe vectors for moment estimation. Default: 20.
%   Nmoment:  Number of Chebyshev moments. Default: 1000.
%   Nbin:   Number of histogram bins. Default: 50.
%
function demo_dos(method, dname, Nprobe, Nmoment, Nbin)

if nargin < 1, method = 'cheb'; end
if nargin < 2, dname = 'Erdos02-cc'; end
if nargin < 3, Nprobe = 20;  end
if nargin < 4, Nmoment = 1000; end
if nargin < 5, Nbin  = 50;   end

% Load data and precompute.
tic;
if ischar(dname)
    [A,lambda] = load_graph(dname);
    lambda = 1-lambda;
else
    A = dname;
end
if length(A) > 1e4
    warning('A might be too large for exact computation.');
end
  
N = matrix_normalize(A, 's');
% Compute eigenvalues if not possible.
if ~exist('lambda','var')
    lambda = eig(full(N));
end
fprintf('Time to load and convert: %g\n', toc);

% Compute spectral density as Chebyshev moments.
tic;
switch method
    case 'cheb'
        c = moments_cheb_dos(N, Nprobe, Nmoment);
    case 'lan'
        c = moments_lan_dos(N, Nprobe, Nmoment);
    case 'exact'
        c = moments_exact_dos(N, Nmoment);
    case 'nd'
        c = moments_nd_dos(N, Nmoment);
    otherwise
        error('Method unrecognized.');
end

% Apply Jackson filter to remove Gibbs phenomenon.
cf = filter_jackson(c);
fprintf('Time to compute Chebyshev moments: %g\n', toc);

% Plot
figure('outerposition',[0 0 1050 900]);
hold on
box on
lmin = max(min(lambda),-1);
lmax = min(max(lambda), 1);
x = linspace(lmin,lmax,Nbin+1);
xm = (x(1:end-1)+x(2:end))/2;
hist(lambda, Nbin)
plot(xm, plot_chebhist(cf,x), 'r.', 'markersize', 60);
xlim([-1 1]);
xlabel('\lambda')
ylabel('Count')
set(gca,'xtick',linspace(-1,1,11),'FontSize',40,'FontWeight','bold');
hold off;