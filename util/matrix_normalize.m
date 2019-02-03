% [N] = matrix_normalize(W, mode)
%
% Normalize a weighted adjacency matrix.
%
% Input:
%   W: weighted adjacency matrix
%   mode: string indicating the style of normalization; this may be
%     's': Symmetric scaling by the degrees (default)
%     'r': Normalize to row-stochastic
%     'c': Normalize to col-stochastic
%
% Output:
%   N: a normalized adjacency matrix or stochastic matrix (in sparse form)
%
function [N] = matrix_normalize(W, mode)

  dc = full(sum(W,1)).';
  dr = full(sum(W,2));
  [i,j,wij] = find(W);

  if nargin < 2, mode = 's'; end
  if any(strcmp(mode, {'s', 'l'}))
    wij = wij./sqrt(dr(i).*dc(j));
  elseif strcmp(mode, 'r')
    wij = wij./dr(i);
  elseif strcmp(mode, 'c')
    wij = wij./dc(j);
  else
    error('Unknown mode!');
  end

  N = sparse(i,j,wij,size(W,1),size(W,2));

end
