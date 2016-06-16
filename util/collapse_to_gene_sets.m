function [ sy_sets, num_genes, effnum_genes, pexp, coef ] = collapse_to_gene_sets( sy, tids, sets, varargin )
% sy = standardizedd [samples x genes]
% tids = genes long vector of transcript ids.
% sets = cell array, where sets{i} = is a cell array of gene ids for set i.

method = setParam(varargin, 'method', 'eigengene');

sy_sets = zeros(size(sy,1), length(sets));
num_genes = zeros(1, length(sets));
effnum_genes = num_genes;
pexp = num_genes;
coef = zeros(size(sy,2), length(sets));
for i = 1 : length(sets)
    mask = steq(tids, sets{i});
    [c, s, l, ~, pex] = pca(sy(:,mask), 'NumComponents', 1);
    m = mean(sy(:,mask),2);
    pexp(i) = pex(1);
    coef(mask,i) = c(:,1);
    
    S = sqrt(l);
    p = S/sum(abs(S));
    effnum_genes(i) = exp( -sum( p.*log(p) ) );
    num_genes(i) = sum(mask);
    
    if strcmpi(method, 'eigengene')
        rr = sign(corr(s(:,1),sy(:,mask))); 
        sm = 1;
        if mean(rr) < 0
            sm = -1;
        end
        coef(mask,i) = sm*coef(mask,i);
        sy_sets(:,i) = sy(:,mask)*coef(mask,i); %s(:,1);
    elseif strcmpi(method, 'average')
        sy_sets(:,i) = m;
    end
    
    
end





end

