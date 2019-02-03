function vt = vt_hash(A)
    n = size(A,1);
    s = randn(n,1);
    h = A*s+s;
    [~,~,IC] = uniquetol(h);
    k = max(IC);
    vt = cell(k,1);
    kt = 0;
    for i = 1:k
        idx = find(IC==i);
        if length(idx) > 1
            kt = kt+1;
            vt{kt} = idx;
        end
    end
    vt(kt+1:end) = [];
end