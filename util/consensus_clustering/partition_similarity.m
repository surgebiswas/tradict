function S = partition_similarity( rp, qp )
% Calculates the similarity matrix between a reference partition rp and a
% query partition qp
%
% S = max(rp) x max(qp) similarity matrix where S(i,j) denotes the
% intersection percentage (|intersection|/(|set1| + |set2|)) between
% cluster indices i and j. 
%
% Assumes (for now) that rp and qp cluster indices are 1,2,...,k.

S = zeros(max(rp),max(qp));
for i = 1 : max(rp)
    for j = 1 : max(qp)
        rpidx = rp == i;
        qpidx = qp == j;
        
        S(i,j) = sum(rpidx & qpidx) / ( sum(rpidx) + sum(qpidx) - sum(rpidx & qpidx) );
    end
end

 

end

