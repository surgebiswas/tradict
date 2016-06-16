function [ ca ] = relabel_base_clustering( c )
% Assumes cluster labels in c are from 1,2,...,k


ca = zeros(size(c));

[~,rcind] = max(max(c));
rc = c(:,rcind);
ca(:,rcind) = rc;

for i = setdiff(1:size(c,2),rcind)
    A = munkres(-partition_similarity(rc,c(:,i)));
    ca(:,i) = relabel(c(:,i), A);
end

    function y = relabel(x,A)
        y = zeros(size(x));
        for j = 1 : max(x)
            y(x == j) = find(A(:,j));
        end
    end


end

