function pred = tradict_predict( T_m, o, model, varargin )
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
    calc_credible_intervals = setParam(varargin, 'calc_credible_intervals', true);
    cred_interval_size = 0.95;

    % First learn the marker latent abundances.
    mu_m = model.fit.markers.mu;
    Sigma_m = model.fit.markers.Sigma;
    if ldiag
        z_m = learn_pmvn_z(t_m, mu_m, diag(diag(Sigma_m)), o);
    else
        z_m = learn_pmvn_z(t_m, mu_m, Sigma_m, o);
    end
    
    % Predict transcriptional program/gene-set/pathway scores
    SIL_gs = model.fit.markers.Sigma \ model.fit.geneset.cross_cov;
    [s_hat, s_cond_var] = predict_program_expression(z_m, model, SIL_gs);
    pred.programs.s_hat = s_hat;
    
    % Predict gene expression
    SIL_gene = model.fit.markers.Sigma \ model.fit.genes.cross_cov;
    [z_hat, z_cond_var] = predict_gene_expression(z_m, model, SIL_gene);
    pred.genes.z_hat = z_hat;
    
    
    
    if calc_credible_intervals
        left_pctile = (1 - cred_interval_size)/2;
        right_pctile = 1 - left_pctile;
        
        % First draw from the posterior of z_m
        zm_smpl = sample_marker_latent_abundances(t_m, z_m, o, ...
            model.fit.markers.mu, model.fit.markers.Sigma);
        
        % Given these draws from the posterior of z_m, we can sample for
        % the transcriptional program and gene expression values.
        s_left = zeros(size(s_hat));
        s_right = s_left;
        z_left = zeros(size(z_hat));
        z_right = z_left;
        for i = 1 : length(zm_smpl)
            
            mu_gs = repmat(model.fit.geneset.mu, size(zm_smpl{i},1), 1) + ...
                bsxfun(@minus, zm_smpl{i}, model.fit.markers.mu)*SIL_gs;
            
            mu_gene = repmat(model.fit.genes.mu, size(zm_smpl{i},1), 1) + ...
                bsxfun(@minus, zm_smpl{i}, model.fit.markers.mu)*SIL_gene;
            
            
            pw_smpl = mvnrnd(mu_gs, s_cond_var);
            gene_smpl = mu_gene + repmat(sqrt(z_cond_var), size(zm_smpl{i},1), 1).*randn(size(mu_gene));
              
            % tr programs
            p = prctile(pw_smpl, 100*[left_pctile, right_pctile]);
            s_left(i,:) = p(1,:);
            s_right(i,:) = p(2,:);
%             for j = 1 : size(pw_smpl,2)
%                 f = ksdensity(pw_smpl(:,j), [left_pctile, right_pctile], ...
%                     'function', 'icdf');
%                 
%                 s_left(i,j) = f(1);
%                 s_right(i,j) = f(2);
%             end
            
            % genes
            p = prctile(gene_smpl, 100*[left_pctile, right_pctile]);
            z_left(i,:) = p(1,:);
            z_right(i,:) = p(2,:);
            
%             for j = 1 : size(gene_smpl,2)
%                 f = ksdensity(gene_smpl(:,j), [left_pctile, right_pctile], ...
%                     'function', 'icdf');
%                 
%                 z_left(i,j) = f(1);
%                 z_right(i,j) = f(2);
%             end
        end
        
        pred.programs.cred_left = s_left;
        pred.programs.cred_right = s_right;
        pred.genes.cred_left = z_left;
        pred.genes.cred_right = z_right;
        
    end
    
    
   
    
    % Accessory functions
    function [gh, gcv] = predict_gene_expression(z_m, model, sil)
        mum = model.fit.markers.mu;
        mu_g = model.fit.genes.mu;
        sig_mg = model.fit.genes.cross_cov;
        
        
        % conditional mean
        gh = repmat(mu_g, size(z_m,1), 1) + ...
            bsxfun(@minus, z_m, mum)*sil;
        
        
        % conditional covariance
        % Note that computing and sampling using the full covariance matrix
        % is too demanding. Importantly, we don't need to because we've
        % assumed that the cross covariance relationships between markers
        % and genes is adequately captured within the markers x genes cross
        % covariance matrix. Thus we only need to store the conditional
        % variance of each gene.
        gcv = model.fit.genes.var - sum(sig_mg.*sil);
        
        
    end
    
    function [sh, scv] = predict_program_expression(z_m, model, sil)
        mum = model.fit.markers.mu;
        mu_s = model.fit.geneset.mu;
        sig_ms = model.fit.geneset.cross_cov;
        
        % conditional mean
        sh = repmat(mu_s, size(z_m,1), 1) + ...
            bsxfun(@minus, z_m, mum)*sil;
        
        % conditional covariance
        % Note that this computation will induce minor computational errors
        % within machine precision causing the matrix to be slightly
        % asymmetric. We therefore symmetrize it. 
        scv = model.fit.geneset.cov - sig_ms'*sil;
        scv = (scv + scv')/2;
        
    end
    

end


