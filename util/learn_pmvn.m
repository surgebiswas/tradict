function [z, mu, Sigma] = learn_pmvn( t, o, z_init)

    z = z_init;
    regularizer = 0.1*diag(diag(cov(z)));
    old_Sig = regularizer;

    delta = Inf;
    old_obj = -Inf;
    while delta > 0.01
        Sigma = cov(z) + regularizer;
        mu = mean(z);
        [z,obj] = learn_pmvn_z( t, mu, Sigma, o, 'z', z );
        
        fprintf('Objective: %0.8f\n', obj);
        delta = max(max(abs( Sigma - old_Sig) )) / (norm(old_Sig,'fro')/size(old_Sig,2)); disp(delta)
        old_obj = obj;
        old_Sig = Sigma;
        
    end
   
    
    
end

