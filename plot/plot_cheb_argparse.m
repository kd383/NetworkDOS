% [c,xx,xx0,ab] = plot_cheb_argparse(npts,c,xx0,ab)
%
% Handle argument parsing for plotting routines.  Should not be
% called directly by users.
%
% Inputs:
%    npts: Number of points in a default mesh
%    c:    Vector of moments
%    xx0:  Input sampling mesh (original coordinates)
%    ab:   Scaling map parameters
%
% Outputs:
%    c:    Vector of moments
%    xx:   Input sampling mesh ([-1,1] coordinates)
%    xx0:  Input sampling mesh (original coordinates)
%    ab:   Scaling map parameters
%
function [c,xx,xx0,ab] = plot_cheb_argparse(npts,c,xx0,ab)

  if nargin < 3

    % Only c is specified
    ab = [1, 0];
    xx0 = linspace(-1+1e-8,1-1e-8,npts);
    xx = xx0;

  elseif nargin < 4
    if length(xx0) == 2

      % Parameters are c, ab
      ab = xx0;
      xx = linspace(-1+1e-8,1-1e-8,npts);
      xx0 = ab(1)*xx + ab(2);

    else

      % Parameters are c, xx0
      ab = [1, 0];
      xx = xx0;

    end
  else

    % All parameters specified
    xx = (xx0-ab(2))/ab(1);

  end
end
