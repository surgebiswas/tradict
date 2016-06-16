function [ c ] = remove_zero_cluster_labels( c )
% Removes cluster indices with 0 representation in the cluster index
% vector.
% 
% Assumes cluster indices run from 1,2,...k.

t = tabulate(c);

while any(t(:,2) == 0)
    for i = 1 : size(t,1)
        if t(i,2) == 0
            c(c > i) = c(c > i) - 1;
            break
        end
    end
    
    t = tabulate(c);
end


end

