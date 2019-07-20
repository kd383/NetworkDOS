% vdf = filter_construct(vd,sqd,nsize)
% construct kernel from doubling vertices 
%
% vd - list of indices of double vertices
% sqd - square root of degrees
% nsize - # of vertices
% 
% last updated: 12-22-2014 (renaming variable and update comments)

function vdf = filter_construct(vd,sqd,nsize)
    n=size(vd,1);
    vdl = zeros(length(vd),1);
    for i = 1:length(vd)
        vdl(i) = length(vd{i});
    end
    nnzvd = sum(vdl);
    npair=2*(nnzvd-n);
    ind=zeros(npair,3);
    count_1=1;
    count_2=1;
    % allocate based on v_double
    for i=1:n
        for j=2:length(vd{i})
            ind(count_2,:)=[count_1,vd{i}(1),sqd(vd{i}(1))];
            ind(count_2+1,:)=[count_1,vd{i}(j),-sqd(vd{i}(j))];
            count_1=count_1+1;
            count_2=count_2+2;
        end
    end
    % construct sparse kernel
    vdf=sparse(ind(:,2),ind(:,1),ind(:,3),nsize,nnzvd-n);
    % find the block sizes. blocks are naturally orthornormal to each other.
    ind2=[0;cumsum(vdl-ones(n,1))];
    % orthonormalization by blocks
    for j=2:n+1
        [vdf(:,ind2(j-1)+1:ind2(j)) ~]=qr(vdf(:,ind2(j-1)+1:ind2(j)),0);
    end
end