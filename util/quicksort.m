function index = quicksort(input, inl_fcn)

% QUICKSORT
% syntax:
%  Index = quicksort(Input, Cmp_func)
%
% Input may be any Matlab data type multidimensional matrix or cell array.
%
% Input may be a multidimensional matrix or cell array.
%  The sort is performed along the first non-singleton
%  dimension of Input. The trivial case of a 1x1 input
%  will not cause a crash.
%
% Cmp_func may be an inline function or function handle.
%  Function handles are much faster than inlines.
%  Cmp_func must take as arguments the input array and
%  two indices. It then returns 1 if the element(s)
%  associated with the first index rank higher than that
%  (those) of the second index. It returns -1 if the same
%  ranks lower and 0 if they are considered equal.
%
% Index is a 1xN index vector which can be used to sort
%  the Input argument.
%
% Example:
%   Sort a two-column matrix of integers by pairs with the
%    first column in ascending order and the second column
%    in descending order when the first column elements are
%    equal.
%
%  % matrix to be sorted
%  x = floor(rand(10000,2)*10);
%
%  % inline comparison function
%  i_cmp = inline( ...
%    'sign(sign(m(x,1)-m(y,1))*10-sign(m(x,2)-m(y,2)))', ...
%    'm','x','y');
%
%  % calculate the sorted index
%  order = quicksort(x, i_cmp);
%
%  % obtain the sorted matrix
%  x_sorted = x(order,:);
%

sz = size(input);
while (sz(1) == 1) && (length(sz) > 1)
    sz = sz(2:end);
end

index = cumsum(ones(1,sz(1)));

if isa(inl_fcn, 'function_handle')
    index = quicksort_handle(input, inl_fcn, index, 1, sz(1));
else
    index = quicksort_inline(input, inl_fcn, index, 1, sz(1));
end

%
% quicksort with an inline 
%
function indx = quicksort_inline(inpt, cmp, indx, l, r)
% l and r remain unchanged. they are the left and right
%  bounds of this sorting level
i = l + 1; % leftmost unknown
j = r; % rightmost unknown
p = l; % rightmost equal

if l >= r
    return
end

% this while loop only runs while p == i - 1
while i <= j
    switch cmp(inpt, indx(i), indx(l))
        case 1
            tmp = indx(j);
            indx(j) = indx(i);
            indx(i) = tmp;
            j = j - 1;
        case -1
            i = i + 1;
            break
        otherwise
            p = p + 1;
            i = i + 1;
    end
end

% this is actually the main while loop
%  in this loop i > p + 1
while i <= j
    switch cmp(inpt, indx(i), indx(l))
        case 1
            tmp = indx(j);
            indx(j) = indx(i);
            indx(i) = tmp;
            j = j - 1;
        case -1
            i = i + 1;
        otherwise
            p = p + 1;
            tmp = indx(p);
            indx(p) = indx(i);
            indx(i) = tmp;
            i = i + 1;
    end
end

% swap "less thans" with "equals"
indx(l:j) = [ indx((p + 1):j) indx(l:p) ];

indx = quicksort_inline(inpt, cmp, indx, l, l + j - p);
indx = quicksort_inline(inpt, cmp, indx, i, r);

return;

%
% quicksort with a function handle 
%
function indx = quicksort_handle(inpt, cmp, indx, l, r)
% l and r remain unchanged. they are the left and right
%  bounds of this sorting level
i = l + 1; % leftmost unknown
j = r; % rightmost unknown
p = l; % rightmost equal

if l >= r
    return
end

% this while loop only runs while p == i - 1
while i <= j
    switch feval(cmp, inpt, indx(i), indx(l))
        case 1
            tmp = indx(j);
            indx(j) = indx(i);
            indx(i) = tmp;
            j = j - 1;
        case -1
            i = i + 1;
            break
        otherwise
            p = p + 1;
            i = i + 1;
    end
end

% this is actually the main while loop
%  in this loop i > p + 1
while i <= j
    switch feval(cmp, inpt, indx(i), indx(l))
        case 1
            tmp = indx(j);
            indx(j) = indx(i);
            indx(i) = tmp;
            j = j - 1;
        case -1
            i = i + 1;
        otherwise
            p = p + 1;
            tmp = indx(p);
            indx(p) = indx(i);
            indx(i) = tmp;
            i = i + 1;
    end
end

% swap "less thans" with "equals"
indx(l:j) = [ indx((p + 1):j) indx(l:p) ];

indx = quicksort_handle(inpt, cmp, indx, l, l + j - p);
indx = quicksort_handle(inpt, cmp, indx, i, r);

return;
