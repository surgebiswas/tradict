function [cidx, cidx_soft] = conclust( c )
% Consensus clustering
%
% Multiple cluster index assignments. n-observations x k-clusterings matrix
% of cluster indices. See apcluster_boot

ca = relabel_base_clustering(c);

cidx_soft = zeros(size(ca,1), max(max(ca)));
cidx = zeros(size(ca,1),1);
for i = 1 : size(cidx_soft,1)
    for j = 1 : size(cidx_soft,2)
        cidx_soft(i,j) = sum(ca(i,:) == j)/size(ca,2);
    end
    [~,cidx(i)] = max(cidx_soft(i,:));
end


end

