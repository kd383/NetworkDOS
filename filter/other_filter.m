function Q = other_filter(A)
        vd=vt_hash(A);
        sqd = sqrt(sum(A).^(1/2));
        nsize = size(A,1);
        Q = filter_construct(vd,sqd,nsize);
        % apply the filter
        %demof(dname,Q);
end