% yy = plot_chebhist(c,xx,ab)
%
% Given a (possibly filtered) set of first-kind Chebyshev moments,
% compute the integral of the density
%
%   int_{0}^s (2/pi)*sqrt(1-x^2)) * ( c(0)/2 + sum_{n=1}^{N-1} c_n T_n(x) )
%
% If no output argument is assigned, make a plot.
%
% Inputs:
%   c:  Chebyshev moments (on [-1,1])
%   xx: evaluation points (defaults to mesh of 11 pts)
%   ab: mapping parameters (defaults to identity)
%
% Output:
%   yy: Estimated counts on buckets between xx points

function yy = plot_chebhist(varargin)

% Parse arguments
[c,xx,xx0,ab] = plot_cheb_argparse(21, varargin{:});

% Plot
yy = plot_chebint(c,xx0,ab);
yy = yy(2:end)-yy(1:end-1);
xm = (xx0(2:end)+xx0(1:end-1))/2;
if nargout < 1
    bar(xm, yy);
end