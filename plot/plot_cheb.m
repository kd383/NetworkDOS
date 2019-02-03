% yy = plot_cheb(c,xx,ab)
%
% Given a set of first-kind Chebyshev moments, compute the associated
% density.   If no output argument is assigned, make a plot.
%
% Inputs:
%   c:  Chebyshev moments (on [-1,1])
%   xx: evaluation points (defaults to mesh of 1001 pts)
%   ab: mapping parameters (defaults to identity)
%
% Output:
%   yy: Density evaluated at xx mesh

function yy = plot_cheb(varargin)

  % Parse arguments
  [c,xx,xx0,ab] = plot_cheb_argparse(1001, varargin{:});

  % Run the recurrence
  kind = 1;
  N = length(c);
  P0 = xx*0+1;
  P1 = kind*xx;
  yy = c(1)/(3-kind) + c(2)*xx;
  for np = 3:N
    Pn = 2*(xx.*P1) - P0;
    yy = yy + c(np)*Pn;
    P0 = P1;
    P1 = Pn;
  end

  % Do the normalization
  if kind == 1
    yy = (2/pi/ab(1))*( yy./(1e-12+sqrt(1-xx.^2)) );
  else
    yy = (2/pi/ab(1))*( yy.*(sqrt(1-xx.^2)) );
  end

  % Plot if appropriate
  if nargout < 1
    plot(xx0, yy);
  end

end
