function o = pmvn_objective(t, z, o, mu, S)
    lambda = exp(z).*repmat(o, 1, size(z,2));

    poiss_comp = sum(sum( t.*log(lambda) - lambda  ));
    gauss_comp = sum(logmvnpdf(z, mu, S));
    o = poiss_comp + gauss_comp;
end