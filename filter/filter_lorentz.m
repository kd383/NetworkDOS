% [c] = filter_lorentz(c, lambda)
%
% Apply the Lorentz filter to a sequence of Chebyshev moments.
% The moments should be arranged column by column.  Weisse et al
% suggest lambda between 3 and 5.
%
function [c] = filter_lorentz(c, lambda)

  if nargin < 2, lambda = 4; end
  N = size(c,1);
  n = 0:N-1;
  g = sinh(lambda*(1-n/N))/sinh(lambda);
  c = bsxfun(@times, g', c);
