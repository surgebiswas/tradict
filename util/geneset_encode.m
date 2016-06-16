function [ stats ] = geneset_encode( sY, nmarkers, stats, varargin )
% stats = struct ouptut from geneset_cluster.m

% Expression optimization
% expdelta is allowed neighborhood size in # markers. 
% if expdelta is non-zero then mean_expression needs to be suplied. This is
% a gene long vector of average expressions for each gene.
expdelta = setParam(varargin, 'expression_delta', 0);
meanexp = setParam(varargin, 'mean_expression', []);
topN = setParam(varargin, 'topN', false);

memmat = stats.geneset.coef ~= 0;

S = zeros(1, length(nmarkers));
R = standardize(stats.geneset.sy_sets);
Y = R; 
tvy = tv(Y); tvr = zeros(1,length(nmarkers));
selected_clusters = zeros(1,length(nmarkers));

if ~topN
    for i = 1 : nmarkers
        cpunexp = calculate_cluster_unexplained_variance(R,stats);

        [maxpex, mind] = max(cpunexp);
        selected_clusters(i) = mind;

        [S, R] = add_marker_from_geneset_cluster(sY, Y, S,  R, mind, memmat, stats, expdelta, meanexp);

        tvr(i) = tv(R);

        fprintf('%0.0f\tPercent unexplained variance: %0.2f\n', i, 100*tvr(i)/tvy);
    end
else
    assert(~isempty(meanexp))
    [~,sind] = sort(meanexp, 'descend');
    S = sind(1:nmarkers);
end

stats.S = S;
stats.selected_clusters = selected_clusters;
stats.punexp = tvr./tvy;

    function cp = calculate_cluster_unexplained_variance(r,stats)
        cpall = zeros(1, size(r,2));
        for j = 1 : length(cpall)
            cpall(j) = var( r(:,j) );
        end
        
        cp = zeros(1, max(stats.conclust.c));
        for j = 1 : length(cp)
            cp(j) = mean( cpall(stats.conclust.c == j) );
        end
        
    end

    function [s, r] = add_marker_from_geneset_cluster(sY, y, s, r, idx, memmat, stats, expdelta, meanexp)
        
        % rsub is a subset of the residual belonging to genesets 
        rsub = r(:, stats.conclust.c == idx);
        
        % Pick a gene from those that are contained in the cluster of genesets
        % defined by idx.
        allowedgenes = find(sum(memmat(:, stats.conclust.c == idx),2) > 0);
        
        pc = corr(sY(:, allowedgenes), rsub);
        
        pcavgabs = mean(abs(pc),2); 
        
        if expdelta == 0
            [~,bestind] = max(pcavgabs);
        else
            [~,sidx] = sort(pcavgabs, 'descend');
            
            topind = sidx(1:min(length(pcavgabs),expdelta));
            mexps = meanexp(allowedgenes(topind));
            mtarg = prctile(meanexp, 75);
            
            % pick the gene within the neighborhood
            % with expression closest to the average expression.
            [~,bi] = min(abs(mexps - mtarg)); 
            bestind = topind(bi);
        end
        s(i) = allowedgenes(bestind); % 'i' in main code body.
        
        
        r = update_residual(y, sY, s);
        
    end

    function r = update_residual(y, sY, s)
        
        % Provide a linear update for now. 
        x = sY(:, s);
        
        b = x'*x \ x'*y;
        r = y - x*b;
    end

    function v = tv(x)
        v = sum(var(x));
    end




    
end

