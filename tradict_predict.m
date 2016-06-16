function [ s_hat, t_hat, z_hat ] = tradict_predict( t_m, o, model, varargin )

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


