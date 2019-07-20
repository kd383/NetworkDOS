% demof(dname, Ncheb, Nbin,filter)
% 
% Loads the graph named by dname, computes Ncheb chebyshev moments
% filter by the kernel
% (default: 1000) and plots the true and (filtered) Chebyshev estimates 
% of a histogram with Nbin bins (default:50).
%
% last update: 12-23-2014 (use a new compare_chebhistf method to output error) 

function demo_ldosf(dname,filter, randInit, Ncheb, Nbin )

  if nargin < 4, Ncheb = 500; end
  if nargin < 5, Nbin  = 50;   end

  tic; 
%   [A,lambda] = load_graph(dname);
  p = load('../data/as199911152.mat');
  A = p.A;
  lambda = 1-p.lambda;
  N = matrix_normalize(A);
  fprintf('Time to load and convert: %g\n', toc);
  
  if nargin < 2
      filter = zero_filter(A);
  end
 
  tic;
  if nargin<3
      c = moments_cheb_ldosf(N, Ncheb, 20, 1, filter);
  else
      c = moments_cheb_ldosf(N, Ncheb, 20, 1, filter, randInit);
  end
  c = filter_jackson(c);
  fprintf('Time to compute moments: %g\n', toc);

  plot_cheb_ldos(c,filter);
  
end
  % c = moments_cheb(A, N, num_samples, kind, filter)
%
% Estimate Chebyshev moments c_0, ... c_{N-1} where
%   c_n = sum_i T_n(lambda_i) = tr(T_n(A))
% using an estimator based on num_samples Gaussian vectors.
% filtered by the kernel
%
% last updated: 12-23-2014 (give the option to use specific random vectors)

function [c] = moments_cheb_ldosf(A, N, num_samples, kind, filter, randInit)

  if nargin < 5, flag = 0; else flag = 1;       end
  if nargin < 4, kind = 1;          end
  if nargin < 3, num_samples = 100; end
  if nargin < 2, N = 10;            end

  m = length(A);
  % use given random vectors or generate them
  if nargin<6
      Z = sign(randn(m,num_samples));
  else
      Z = randInit;
  end
  % taking orthogonal component of filter
  if 1
    Z = Z - filter*(filter'*Z);
  end
  
  c = zeros(N,m);

  TZp = Z;
  X = Z.*TZp;
  c(1,:) = mean(X,2).';

  TZk = kind*A*Z;
  X = Z.*TZk;
  c(2,:) = mean(X,2).';

  for k = 3:N
    TZ = 2*A*TZk - TZp;
    TZp = TZk;
    TZk = TZ;
    X = Z.*TZk;
    c(k,:) = mean(X,2).';
  end
end

% yy = plot_cheb_ldos(c,xx,ab)
%
% Given a set of first-kind Chebyshev moments, compute the associated
% density.   If no output argument is assigned, make a plot.
%
% Inputs:
%   c:  Array of Chebyshev moments (on [-1,1])
%   xx: evaluation points (defaults to mesh of 1001 pts)
%   ab: mapping parameters (defaults to identity)
%
% Output:
%   yy: Density evaluated at xx mesh (size nnodes-by-nmesh)
%   idx: Index for spectral re-ordering

function [yy,idx] = plot_cheb_ldos(c,filter)

  % Parse arguments
  xx0 = linspace(-1+1e-8,1-1e-8,52);
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
  temp = find(xx>0,1,'first');
  % Difference the CDF to compute histogram
  yy = yy(:,2:end)-yy(:,1:end-1);
  %yy(:,temp-1) = yy(:,temp-1)+sum(filter.^2,2);

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

end