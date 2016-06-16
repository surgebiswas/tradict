function [ stats ] = geneset_cluster( sY, tids, sets, varargin )

stats = setParam(varargin, 'stats', []);

% Collapse to gene sets. 
[ stats.geneset.sy_sets, stats.geneset.num_genes, stats.geneset.effnum_genes,  ... 
    stats.geneset.pexp, stats.geneset.coef ] = collapse_to_gene_sets( sY, tids, sets );
stats.tids = tids;
stats.geneset.sets = sets;

% Cluster gene sets using consensus clustering.
[stats.conclust.c, stats.conclust.csoft, stats.conclust.cmat] = ...
    kmeans_conclust(standardize(stats.geneset.sy_sets)');

end

