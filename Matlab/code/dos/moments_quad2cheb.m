% [c] = moments_quad2cheb(thetas, wts, N, ab)
%
% Compute the first N moments in a dual Chebyshev expansion for a
% density defined by a set of weighted delta functions.
%
% Inputs:
%   thetas: Delta centers
%   wts:    Weights
%   N:      Number of Chebyshev moments
%   ab:     Linear scaling parameters
%
function [c] = moments_quad2cheb(thetas, wts, N, ab)

if nargin > 3, thetas = (thetas-ab(2))/ab(1); end
c = zeros(N,1);
for k = 0:N-1
    c(k+1) = wts'*cos(k*acos(thetas));
end
