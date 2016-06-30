function [ s_hat, t_hat, z_hat ] = tradict_predict( T_m, o, model, varargin )
% [ s_hat, t_hat, z_hat ] = tradict_predict( T_m, o, model, varargin )
% 
% Description: Uses the expression measurements of the selected marker genes, 
% to formulate predictions of the expression of transcriptional programs and
% the remaining non-marker genes.
%
% INPUT
% T_m =     [n-samples x m-marker-genes] matrix of measured expression values 
%           for the selected markers.
% o =       [n-samples x 1] vector of sequencing depths (in millions of reads)
%           associated with each sample.
% model =   model object output from tradict_train.m.
%
% OUTPUT
% s_hat =   [n-samples x s-transcriptional-programs] matrix of predicted 
%           expression values of the transcriptional programs defined during
%           training.
% t_hat =   [n-samples x #-genes] matrix of predicted expression values in TPM
%           of all genes in the transcriptome.
% z_hat =   [n-samples x #-genes] matrix of predicted expression values in log-
%           TPM of all genes in the transcriptome. 


    t_m = T_m; clear T_m;
    
    ldiag = setParam(varargin, 'learn_latent_diag', false);

    % First learn the marker latent abundances.
    mu_m = model.fit.markers.mu;
    Sigma_m = model.fit.markers.Sigma;
    if ldiag
        z_m = learn_pmvn_z(t_m, mu_m, diag(diag(Sigma_m)), o);
    else
        z_m = learn_pmvn_z(t_m, mu_m, Sigma_m, o);
    end
    
    % Predict gene-set/pathway scores
    mu_s = model.fit.geneset.mu;
    sig_ms = model.fit.geneset.cross_cov;
    s_hat = repmat(mu_s, size(z_m,1), 1) + ...
        bsxfun(@minus, z_m, mu_m)*(Sigma_m \ sig_ms);
    
    
    % Predict TPM for all remaining genes.
    mu_j = model.fit.genes.mu;
    sig_j = model.fit.genes.var;
    sig_mj = model.fit.genes.cross_cov;
    
    H = (Sigma_m \ sig_mj);
    mu_j_given_m = repmat(mu_j, size(z_m,1), 1) + ...
        bsxfun(@minus, z_m, mu_m)*H;
    sig_j_given_m = sig_j - sum(sig_mj .* H);
    
    z_hat = mu_j_given_m;
    t_hat = exp( mu_j_given_m + 0.5*repmat(sig_j_given_m, size(t_m,1), 1) );
    

end


