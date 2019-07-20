% [c] = filter_jackson(c)
%
% Apply the Jackson filter to a sequence of Chebyshev moments.
% The moments should be arranged column by column.
%
function [c] = filter_jackson(c)

  N = size(c,1);
  n = 0:N-1;
  tau = pi/(N+1);
  g = ( (N-n+1).*cos(tau*n) + sin(tau*n)*cot(tau) )/(N+1);
  c = bsxfun(@times, g', c);
