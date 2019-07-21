% Q = filter_construct(vinfo, sqd, n)
% construct a projection matrix Q representing the a filter.
%
% Input:
%   vinfo: Sets of nodes that form the filter.
%   sqd:   Square root of degrees.
%   n:     Number of vertices.
% 
% Output:
%   Q:     A matrix with eigenvectors of the target eigenvalue as columns.
%
function Q = filter_construct(vinfo, sqd, n)

nf = size(vinfo,1);
vfl = zeros(length(vinfo), 1);
for i = 1:length(vinfo)
    vfl(i) = length(vinfo{i});
end
nnzvf = sum(vfl);
npair = 2*(nnzvf - nf);
ind = zeros(npair,3);
count1 = 1;
count2 = 1;
% allocate based on vinfo
for i = 1:nf
    for j = 2:length(vinfo{i})
        ind(count2, :)=[count1, vinfo{i}(1), sqd(vinfo{i}(1))];
        ind(count2+1,:)=[count1, vinfo{i}(j), -sqd(vinfo{i}(j))];
        count1 = count1 + 1;
        count2 = count2 + 2;
    end
end
% construct sparse kernel
Q = sparse(ind(:, 2), ind(:, 1), ind(:, 3), n, nnzvf - nf);
% find the block sizes. blocks are naturally orthornormal to each other.
ind2 = [0; cumsum(vfl - ones(nf, 1))];
% orthonormalization by blocks
for j = 2:nf+1
    [Q(:, ind2(j-1)+1:ind2(j)), ~] = qr(Q(:, ind2(j-1)+1:ind2(j)),0);
end