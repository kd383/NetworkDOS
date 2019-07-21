% demo_HepTh(method, dname, Nprobe, Nmoment, Nbin)
%
% Apply motif filtering on the Arxiv High-Energy-Physics Theory Collbaration
% Netowrk. We use Kernel Polynomial Method with 200 moments and 20 probe
% vectors.
%
function demo_HepTh

Nprobe = 20;
Nmoment = 200;
Nbin = 50;

tic;
path = fileparts(mfilename('fullpath'));
datapath=strrep(path,'demo','data/HepTh.mat');
load(datapath);
N = matrix_normalize(A, 's');
fprintf('Time to load and convert: %g\n', toc);

% The eigenvalues we will be applying filters at.
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

c = moments_filter_dos('cheb', filter, N, Nprobe, Nmoment);
% Apply Jackson filter to remove Gibbs phenomenon.
cf = filter_jackson(c);
fprintf('Time to compute filtered Chebyshev moments: %g\n', toc);

% Plot
figure('outerposition',[0 0 1050 900]);
hold on
box on
lmin = max(min(lambda),-1);
lmax = min(max(lambda), 1);
x = linspace(lmin,lmax,Nbin+1);
xm = (x(1:end-1)+x(2:end))/2;
hist(lambda, Nbin)
yy = plot_chebhist(cf,x);
for i = 1:length(filter_value)
    ind = find(x>filter_value(i),1,'first') - 1;
    yy(ind) = yy(ind) + size(filter(i).Q, 2);
end
plot(xm, yy, 'r.', 'markersize', 60);
xlim([-1 1]);
xlabel('\lambda')
ylabel('Count')
set(gca,'xtick',linspace(-1,1,11),'FontSize',40,'FontWeight','bold');
hold off;