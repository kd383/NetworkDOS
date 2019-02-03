% [c, cs] = moments_cheb_ldos(Afun, n, nZ, N, kind)
% [c, cs] = moments_cheb_ldos(A, nZ, N, kind)
%
% Compute a column vector (or vectors) of Chebyshev moments of
% the form c(k,j) = [T_k(A)]_jj for k = 0 to N-1.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
% The diagonal entries are computed by a stochastic estimator
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    kind: 1 or 2 for first or second kind Chebyshev functions
%          (default is 1)
%
% Output:
%    c: an N-by-n matrix of moments
%    cs: standard deviation of the moments
%
function [c,cs,X] = moments_cheb_ldos(varargin)

  defaults = {'Afun', NaN, 'n', NaN, ...
              'nZ', 100, 'N', 10, 'kind', 1};
  [Afun, n, nZ, N, kind] = mfuncify(defaults, varargin{:});
  if N < 2, N = 2; end

  % Set up random probe vectors (we allow them to be passed in, too)
  if length(nZ) > 1
    Z = nZ;
    nZ = size(Z,2);
  else
    Z = sign(randn(n,nZ));
  end

  % Run three-term recurrence to estimate moments.
  % Uses the stochastic diagonal estimator of Bekas and Saad
  %  http://www-users.cs.umn.edu/~saad/PDF/umsi-2005-082.pdf

  c = zeros(N,n);
  if nargout > 1, cs = zeros(N,n); end

  TZp = Z;
  X = Z.*TZp;
  c(1,:) = mean(X,2).';
  if nargout > 1, cs(1,:) = std(X,0,2).'; end

  TZk = kind*Afun(Z);
  X = Z.*TZk;
  c(2,:) = mean(X,2).';
  if nargout > 1, cs(2,:) = std(X,0,2).'; end

  for k = 3:N
      tic;
    TZ = 2*Afun(TZk) - TZp;
    TZp = TZk;
    TZk = TZ;
    X = Z.*TZk;
    c(k,:) = mean(X,2).';
    toc
    if nargout > 1, cs(k,:) = std(X,0,2).'; end
  end
  if nargout > 1, cs = cs/sqrt(nZ); end

end
