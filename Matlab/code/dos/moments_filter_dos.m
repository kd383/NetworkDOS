% [c, cs] = moments_filter_dos(method, filter, Afun, n, nZ, N, kind)
% [c, cs] = moments_filter_dos(method, filter, A, nZ, N, kind)
%
% Compute Chebyshev moments with the help of motif filtering.
%
% Inputs:
%    method: 'cheb' or 'lan'.
%    filter: An array of struct with two fields, .value and .Q.
%    A: Matrix or function to apply matrix (to multiple RHS)
%    n: Dimension of the space (if A is a function)
%    nZ: Number of probe vectors with which we want to compute moments
%    N: Number of moments to compute
%    kind: 1 or 2 for first or second kind Chebyshev functions
%          (default is 1)
%
% Output:
%    c: an column vector of N moment estimates
%    cs: standard deviation of the moment estimator (std/sqrt(nZ))
%
function [c, cs] = moments_filter_dos(method, filter, varargin)

if nargin < 1, method = 'cheb'; end
if nargin < 2, filter = []; end
defaults = {'Afun', NaN, 'n', NaN, 'nZ', 100, 'N', 10, 'kind', 1};
[Afun, n, nZ, N, kind] = mfuncify(defaults, varargin{:});
if N < 2, N = 2; end

% Set up random probe vectors (we allow them to be passed in, too)
if length(nZ) > 1
    Z = nZ;
    nZ = size(Z, 2);
else
    Z = sign(randn(n, nZ));
end

% Filter the probe vectors
if ~isempty(filter) > 0
    Q = [filter(:).Q];
end

Z = Z - Q*(Q'*Z);

switch method
    case 'cheb'
        [c, cs] = moments_cheb_dos(Afun, n, Z, N, kind);
    case 'lan'
        [c, cs] = moments_lan_dos(Afun, n, Z, N, kind);
    otherwise
        error('Method unrecognized.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expanding the filtered eigenvalues as delta functions into Chebshev
% moments will essentially give back the coefficients without filtering.
% We highly recommand manually add those back when plotting spectral
% histogram.
% 
% theta = [filter(:).value];
% wts = theta;
% for i = 1:length(wts)
%     wts(i) = size(filter(i).Q, 2);
% end
% cf = moments_quad2cheb(theta, wts, N);
% c = c + cf;