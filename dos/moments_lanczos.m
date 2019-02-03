% [theta, wts] = moments_lanczos(Afun, n, q1, kmax, btol)
% [theta, wts] = moments_lanczos(A, q1, kmax, btol)
%
% Run a basic Lanczos iteration until either beta_k < btol
% or k == kmax.  Return weights wts and nodes theta for
% a quadrature formula associated with the DoS starting from q1.
%
% Inputs:
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    q1: Starting vector
%    kmax: Maximum number of steps to take
%    btol: Tolerance on off-diagonal entry beta
%
% Output:
%    theta: quadrature nodes
%    wts: quadrature weights
%
function [theta, wts] = moments_lanczos(varargin)

  defaults = {'Afun', NaN, 'n', NaN, ...
              'q1', NaN, 'kmax', 100, 'btol', 1e-6};
  [Afun, n, q1, kmax, btol] = mfuncify(defaults, varargin{:});

  k  = 0;
  qk = 0;
  n1 = norm(q1);
  r  = q1/n1;
  b  = 1;
  while (b > btol) & (k < kmax)
    k        = k+1;
    qkm1     = qk;
    qk       = r/b;
    Aqk      = Afun(qk);
    alpha(k) = qk'*Aqk;
    r        = Aqk-qk*alpha(k)-qkm1*b;
    b        = norm(r);
    beta(k)  = b;
  end

  % NB: Would be great to call through to LAPACK and do the tridiag solve
  T = diag(alpha) + diag(beta(1:end-1),1) + diag(beta(1:end-1),-1);
  [V,D] = eig(T);
  wts = (V(1,:).').^2 * n1^2;
  theta = diag(D);

end
