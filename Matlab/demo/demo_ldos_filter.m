% demo_ldos_filter(dname, Nprobe, Nmoment, Nbin)
%
% Motif filtering for LDOS using Chebyshev approximation
%
% Input:
%   dname:  Data set name, or an adjacency matrix A.
%   Nprobe: Number of probe vectors for moment estimation. Default: 20.
%   Nmoment:  Number of Chebyshev moments. Default: 1000.
%   Nbin:   Number of histogram bins. Default: 50.
%
function demo_ldos_filter(dname, Nprobe, Nmoment, Nbin)

if nargin < 1, dname = 'Erdos02-cc'; end
if nargin < 2, Nprobe = 20;  end
if nargin < 3, Nmoment = 200; end
if nargin < 4, Nbin  = 50;   end

% Load data and precompute.
tic;
if ischar(dname)
    [A, lambda] = load_graph(dname);
else
    A = dname;
    if length(A) > 1e4
        warning('A might be too large for exact computation.');
    end
end
N = matrix_normalize(A, 's');
fprintf('Time to load and convert: %g\n', toc);

filter_value = [0, -1/2, -1/3, -1/4];

tic;
filter = [];
for i = 1:length(filter_value)
    filter(i).value = filter_value(i);
    if filter(i).value == 0
        filter(i).Q = zero_filter(A);
    elseif filter(i).value == -1/2
        filter(i).Q = clique_filter(A, 2);
    elseif filter(i).value == -1/3
        filter(i).Q = clique_filter(A, 3);
    elseif filter(i).value == -1/4
        filter(i).Q = clique_filter(A, 4);
    end
end
fprintf('Time to compute filter: %g\n', toc);

% Compute spectral density as Chebyshev moments.
tic;
cl = moments_filter_ldos(filter, N, Nprobe, Nmoment);
% Apply Jackson filter to remove Gibbs phenomenon.
cf = filter_jackson(cl);
fprintf('Time to compute Chebyshev moments: %g\n', toc);

plot_cheb_ldos(cf,filter);

function [yy,idx] = plot_cheb_ldos(c, filter)

% Parse arguments
xx0 = linspace(-1+1e-8,1-1e-8,50);
xx = xx0;


% Run the recurrence to compute CDF
[nmoment, nnodes] = size(c);
txx = acos(xx);
yy = c(1,:)'*(txx-pi)/2;
for np = 2:nmoment
    n = np-1;
    yy = yy + c(np,:)' * sin(n*txx)/n;
end
yy = -2/pi * yy;
yy = yy(:,2:end)-yy(:,1:end-1);

for i = 1:length(filter)
    ind = find(xx>filter(i).value,1,'first') - 1;
    yy(:, ind) = yy(:, ind) + sum(filter(i).Q.^2, 2);
end

% Compute sorted histogram
if nargout ~= 1
    [U,S,V] = svds(yy);
    [~,idx] = sort(U(:,1));
end

% Plot if appropriate
if nargout < 1
    yr = [1, nnodes];
    xr = [xx0(1)+xx0(2), xx0(end-1)+xx0(end)]/2;
    figure('outerposition',[0 0 1050 900]);
    bot = min(min(yy)); top = max(max(yy));

    imagesc(xr, yr, yy(idx,:));
    colormap('jet');
    caxis manual
    caxis([bot top])
    xlabel('\lambda')
    ylabel('Node Index')
    xlim([-1 1])
    set(gca,'xtick',linspace(-1,1,11),'FontSize',40,'FontWeight','bold');
end