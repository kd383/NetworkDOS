function vd = vd_hash(A)
    n = size(A,1);
    h = A*randn(n,1);
    [~,~,IC] = uniquetol(h);
    k = max(IC);
    vd = cell(k,1);
    kt = 0;
    for i = 1:k
        idx = find(IC==i);
        if length(idx) > 1
            kt = kt+1;
            vd{kt} = idx;
        end
    end
    vd(kt+1:end) = [];
end