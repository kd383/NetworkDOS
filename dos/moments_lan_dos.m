% [theta, wts] = moments_lan_dos(Afun, n, nZ, N, kmax, btol)
% [theta, wts] = moments_lan_dos(A, nZ, N, kmax, btol)
%
% Compute a column vector (or vectors) of Chebyshev moments of
% the form c(k) = tr(T_k(A)) for k = 0 to N-1.  This routine
% does no scaling; the spectrum of A should already lie in [-1,1].
% The traces are computed via a stochastic estimator with nZ probes.
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    kmax: maximum Lanczos steps (default is 100)
%    btol: Tolerance on Lanczos residual (default is 1e-6)
%
% Output:
%    c:  Chebyshev moments
%    cs: Standard deviations
%
function [c, cs] = moments_lan_dos(varargin)

  defaults = {'Afun', NaN, 'n', NaN, ...
              'nZ', 100, 'N', 10, 'kmax', 100, 'btol', 1e-6};
  [Afun, n, nZ, N, kmax, btol] = mfuncify(defaults, varargin{:});
  if N < 2, N = 2; end

  % Set up random probe vectors (we allow them to be passed in, too)
  if length(nZ) > 1
    Z = nZ;
    nZ = size(Z,2);
  else
    Z = randn(n,nZ);
  end

  % Estimate moments for each probe vector
  cZ = zeros(N,nZ);
  for k = 1:nZ
    [theta,wts] = moments_lanczos(Afun, n, Z(:,k), kmax, btol);
    cZ(:,k) = moments_quad2cheb(theta, wts, N);
  end

  % Get summary statistics
  c = mean(cZ,2);
  if nargout > 1
    cs = std(cZ,0,2)/sqrt(nZ);
  end

end
