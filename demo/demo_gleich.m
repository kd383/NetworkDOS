% demo_gleich(dname, Nprobe, Ncheb, Nbin)
%
% Loads the graph in the 'gleich' collection (with reference eigenvalues)
% and plots histograms computed with true eigenvalue distribution adn with
% the filtered Chebyshev estimate.
%
% Input:
%   dname: Data set name
%   Nprobe: Number of probe vectors for moment estimation
%   Ncheb: Number of Chebyshev moments
%   Nbin:  Number of histogram bins
%
function demo_gleich(dname, Nprobe, Ncheb, Nbin)

  if nargin < 2, Nprobe = 20;  end
  if nargin < 3, Ncheb = 100; end
  if nargin < 4, Nbin  = 50;   end

  tic;
  if ischar(dname)
    [A,lambda] = load_graph(dname);
    lambda = 1-lambda;
  else
      A = dname;
  end
  
  N = matrix_normalize(A, 's');
  if ~exist('lambda','var')
      lambda = eig(full(N));
  end
  fprintf('Time to load and convert: %g\n', toc);

  tic;
  c = moments_cheb_dos(N, Nprobe, Ncheb);
  cf = filter_jackson(c);
  fprintf('Time to compute filtered moments: %g\n', toc);

  figure('outerposition',[0 0 900 900]);
  hold on
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

end
