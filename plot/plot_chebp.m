% yy = plot_chebp(c,xx,ab)
%
% Given a set of first-kind Chebyshev moments, compute the associated
% polynomial (*not* a density).   If no output argument is assigned,
% make a plot.
%
% Inputs:
%   c:  Chebyshev moments (on [-1,1])
%   xx: evaluation points (defaults to mesh of 1001 pts)
%   ab: mapping parameters (defaults to identity)
%
% Output:
%   yy: Density evaluated at xx mesh

function [xx0,yy] = plot_chebp(varargin)

  % Parse arguments
  [c,xx,xx0,ab] = plot_cheb_argparse(1001, varargin{:});

  % Run the recurrence
  kind = 1;
  N = length(c);
  P0 = xx*0+1;
  P1 = kind*xx;
  yy = c(1) + c(2)*xx;
  for np = 3:N
    Pn = 2*(xx.*P1) - P0;
    yy = yy + c(np)*Pn;
    P0 = P1;
    P1 = Pn;
  end

  % Plot if appropriate
  if nargout < 1
    plot(xx0, yy);
  end

end
