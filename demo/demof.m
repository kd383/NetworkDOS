% demof(dname, Ncheb, Nbin,filter)
% 
% Loads the graph named by dname, computes Ncheb chebyshev moments
% filter by the kernel
% (default: 1000) and plots the true and (filtered) Chebyshev estimates 
% of a histogram with Nbin bins (default:50).
%
% last update: 12-23-2014 (use a new compare_chebhistf method to output error) 

function error=demof(dname,filter, randInit, Ncheb, Nbin )

  if nargin < 4, Ncheb = 1000; end
  if nargin < 5, Nbin  = 50;   end

  tic; 
%   [A,lambda] = load_graph(dname);
  p = load('../data/hepth.mat');
  A = p.A;
  lambda = 1-p.lambda;
  filter2 = p.filter2;
  filter3 = p.filter3;
  filter4 = p.filter4;
  N = matrix_normalize(A);
  fprintf('Time to load and convert: %g\n', toc);
  
%   if nargin < 2
%       filter = zero_filter(A);
%   end
    filter = zero_filter(A);
%   t = sum((N*filter).^2);
%   filter = filter(:,t<1e-10);
%   nf = [size(filter,2)];
%   filter = [filter];
  nf = [size(filter,2), size(filter2,2),size(filter4,2)];
  filter = [filter,filter2,filter4];
% filter=[];
  tic;
  if nargin<3
      c = moments_chebf(N, Ncheb, 20, 1, filter);
  else
      c = moments_chebf(N, Ncheb, 20, 1, filter, randInit);
  end
  fprintf('Time to compute moments: %g\n', toc);
  error = compare_chebhistf(1-lambda, filter_jackson(c), nf, Nbin);
  
end
  % c = moments_cheb(A, N, num_samples, kind, filter)
%
% Estimate Chebyshev moments c_0, ... c_{N-1} where
%   c_n = sum_i T_n(lambda_i) = tr(T_n(A))
% using an estimator based on num_samples Gaussian vectors.
% filtered by the kernel
%
% last updated: 12-23-2014 (give the option to use specific random vectors)

function [c] = moments_chebf(A, N, num_samples, kind, filter, randInit)

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
  if ~isempty(filter)
    Z = Z - filter*(filter'*Z);
  end
  %%%%%%%%%%
  c = zeros(N,1);
  
  c(1) = m-size(filter,2);
  c(2) = trace(A);
  P0 = Z;
  P1 = kind*(A*Z);

  for np = 3:N
    Pn = 2*(A*P1) - P0;
    for j = 1:num_samples
      c(np) = c(np) + Z(:,j)'*Pn(:,j);
    end
    c(np) = c(np) / num_samples;
    P0 = P1;
    P1 = Pn;
  end
end
  
  % error = compare_chebhistf(lambda, c, Nbin)
%
% Compare a histogram of the eigenvalue list lambda with an
% estimated histogram from the integrated Chebyshev density estimate.
%
% An alternate version that output the error as well
% Need to manually change the index of zero eigenvalues
%
% last updated: 12-23-2014 (created as a version of the original to output error)

function error=compare_chebhistf(lambda, c, nf, Nbin)

  if nargin < 4, Nbin = 50; end
  lmin = max(min(lambda),-1);
  lmax = min(max(lambda), 1);

  x = linspace(lmin,lmax,Nbin+1);
  y = plot_chebint(c,x);
  u = (x(2:end)+x(1:end-1))/2;
  v = y(2:end)-y(1:end-1);

  figure('outerposition',[0 0 1050 900]);
  hold on
  xlim([-1 1])
  hist(lambda, Nbin)
  % get the count at each histogram bar
  counts=hist(lambda, Nbin);
  temp = find(x>0,1,'first');
  v(temp-1) = v(temp-1) + nf(1);
  temp = find(x>-1/3,1,'first');
  v(temp-1) = v(temp-1) + nf(2);
%   temp = find(x>-1/2,1,'first');
%   v(temp-1) = v(temp-1) + nf(3);
  temp = find(x>-1/4,1,'first');
  v(temp-1) = v(temp-1) + nf(3);
  plot(u, v, 'r.', 'markersize', 60);
  xlabel('\lambda')
  ylabel('Count')
  set(gca,'xtick',linspace(-1,1,11),'FontSize',40,'FontWeight','bold');
  % 25 is usually the index for zero eigenvalues
%   counts(25)=[];
%   v(25)=[];
  error=sum(abs(counts-v));
  xlim([-1,1])
  hold off
%   close all
end