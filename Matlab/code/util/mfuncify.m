% [varargout] = mfuncify(defaults, varargin)
%
% Standardize an argument list of the form
%   (A, rest)
% to
%   (Afun, n, rest)
% where Afun is a function that multiplies by a vector and n
% is the dimension.
%
% Inputs:
%    defaults: Cell array of the form {'name', default, 'name', default}
%              A default of NaN indicates a required argument.
%    varargin: Array of arguments starting with a matrix spec
%
function [varargout] = mfuncify(defaults, varargin)
if isa(varargin{1}, 'function_handle')
    varargout = varargin;
else
    A = varargin{1};
    Afun = @(X) A * X;
    n = size(A, 1);
    varargout{1} = Afun;
    varargout{2} = n;
    varargout(3:nargin) = varargin(2:nargin-1);
end

nout = length(varargout);
for k=2*nout+1:2:length(defaults)
    if isnan(defaults{k+1})
        error(sprintf('Missing required argument %s', defaults{k}));
    else
        varargout{(k+1)/2} = defaults{k+1};
    end
end

if 2*nout > length(defaults)
    error('Too many arguments');
end