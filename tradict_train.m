function model = tradict_train( T, o, tids, sets, varargin )
% model = tradict_train( T, o, tids, sets, varargin )
% 
% Description: Trains a tradict model.
% 
% INPUT
% T =   [n-samples x p-genes] matrix of training transcriptomes. T_{ij} is 
%       number of transcripts measured for gene j in sample i. 
% o =   [n-samples x 1] vector of sequencing depths for each sample.
% tids =    p-genes long vector of transcript IDs. These should be standard
%           gene names.
% sets =    # transcriptional programs long cell array. sets{i} is another
%           cell array that contains the transcript IDs for the transcripts
%           contained in transcriptional program i. 
% 
% OPTIONAL INPUTS
% 'nmarkers' =  number of markers to respresent the transcriptome with 
%               (Default = 100).
% 'topN' =      Use the most abundant 'nmarkers' as the selected marker
%               panel? [true | false] (Defaults = false);
%
% OUTPUT
% model =   A struct that contains all the relevant information from 
%           training that is needed for prediction. It has the following 
%           fields:
%   model.lag_priors =  #-genes long cell array of structs containing the 
%                       prior means and variances learned from the 
%                       lag-transformation.
%   model.train_mu  =   #-genes long vector of the mean of the 
%                       lag-transformed abundances for each gene.
%   model.train_sig =   #-genes long vector of the standard deviations of the 
%                       lag-transformed abundances for each gene.
%   model.geneset   =   struct output from geneset_cluster.m. Structure
%                       containing information regarding the collapsing of
%                       genes into transcriptional programs. Includes
%                       standardized expression of each transcriptional
%                       program across the training set, the number of
%                       genes contained in each program, the effective
%                       number of genes based on the size of the
%                       eigenvalues of each gene, the percent variance
%                       explained of all genes in that transcriptional
%                       program, the [genes x programs] matrix of PC1
%                       coefficients, and a copy of the input 'sets'. See
%                       geneset cluster for more information. 
%   model.tids  =       Copy of input 'tids'.
%   model.conclust =    Struct containing cluster assignments for each
%                       geneset. Clustering is required to prioritize which
%                       genes are selected in each iteration by the SOMP
%                       algorithm.
%   model.S =           Indices of selected markers. tids(model.S) gives
%                       the gene IDs of each selected marker.
%   model.selected_clusters = The SOMP algorithm iteratively selects genes
%                             from program clusters with still high residual
%                             variance. model.selected_clusters gives the
%                             cluster indices of the selected clusters.
%   model.punexp   =    Cumulative unexplained variance of the entire
%                       transcriptional program matrix as a function of
%                       iteration.
%   model.fit      =    struct containing information regarding fit of the
%                       Multivariate Normal Continuous-Poisson model.
%                       Includes information regardin the marker mean and
%                       covariance, and the mean and cross-covariance of
%                       the non-markers and transcriptional programs to the
%                       the markers.

    t = T; clear T;

    nmarkers = setParam(varargin, 'nmarkers', 100);
    expdelta = setParam(varargin, 'expression_delta', 0);
    topN = setParam(varargin, 'topN', false);

    % Perform lag
    [zlag, model.lag_priors] = lag_dataset(t, o);
    
    % perform the endcoding
    meanexp = mean(zlag);
    [sY, model.train_mu, model.train_sig] = standardize(zlag);
    model = geneset_cluster( sY, tids, sets, 'stats', model );
    model = geneset_encode(sY, nmarkers, model, 'expression_delta', expdelta, 'mean_expression', meanexp, 'topN', topN);
    markers = model.S;
    
    % Learn z_m, \mu^{(m)}, and \Sigma^{(m)}
    [model.fit.markers.z, model.fit.markers.mu, model.fit.markers.Sigma] = ...
            learn_pmvn(t(:,markers), o, zlag(:,markers));
    zlag(:,markers) = model.fit.markers.z; % update lag estimates for markers.

    % Learn mean and cross covariances for gene sets (processes)
    model.fit.geneset.mu = mean(model.geneset.sy_sets);
    model.fit.geneset.cov = cov(model.geneset.sy_sets);
    model.fit.geneset.cross_cov = cross_cov(model.fit.markers.z, model.geneset.sy_sets);
    
    % Learn mean and cross covariances for all genes
    model.fit.genes.mu = mean(zlag);
    model.fit.genes.cov = cov(zlag);
    model.fit.genes.cross_cov = cross_cov(model.fit.markers.z, zlag);
    
    function c = cross_cov(x,y)
        c = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);
    end

end

