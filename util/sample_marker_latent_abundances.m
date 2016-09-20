function zm_smpl = sample_marker_latent_abundances( t_m, z_m, o, marker_mu, marker_Sigma )

verbose = true;
if verbose
    fprintf('Tuning acceptance rate of proposal distribution using a randomly chosen observation ... \n');
end

TARGET_ACCEPTANCE_RATE = 0.234;
eidx = randsample(size(z_m,1),1);
pwmin = -3.5; pwmax = -1;
propwidth = logspace(pwmin,pwmax,50);
acc = zeros(length(propwidth),1);
logpdf = @(x) pmvn_objective(t_m(eidx,:), x, o(eidx,:), marker_mu, marker_Sigma);
for i = 1 : length(propwidth)
    proprnd = @(x) propwidth(i)*randn(size(x)) + x;
    [~,acc(i)] = mhsample(z_m(eidx,:), 100, 'logpdf', logpdf, ...
        'symmetric', true, 'proprnd', proprnd, 'burnin', 100);
    fprintf('Proposal sampling width: %0.6f\tAcceptance Rate: %0.6f\n', ...
        propwidth(i), acc(i));
end

% Determine the best proposal width
fprintf('Predicting best proposal width ...\n');
pp = spline(log(propwidth), acc);
pwfine = logspace(pwmin, pwmax, 1000);
eacc = ppval(pp, log(pwfine));

[~,mind] = min(abs(eacc - TARGET_ACCEPTANCE_RATE));
bestwidth = pwfine(mind);
proprnd = @(x) bestwidth*randn(size(x)) + x;
[~,acc_final] = mhsample(z_m(eidx,:), 500, 'logpdf', logpdf, ...
        'symmetric', true, 'proprnd', proprnd, 'burnin', 100);
fprintf('Achieved acceptance rate = %0.6f at width = %0.6f\n', acc_final, bestwidth);


fprintf('Performing sampling at this width for each observation.\n');
zm_boot = cell(size(z_m,1),1);
for i = 1 : length(zm_boot)
    fprintf('Sampling latent marker abundances for observation %0.0f ... ', i);
    logpdf = @(x) pmvn_objective(t_m(i,:), x, o(i,:), marker_mu, marker_Sigma);
    [zm_boot{i}, acc] = mhsample(z_m(i,:), 200, 'logpdf', logpdf, ...
        'symmetric', true, 'proprnd', proprnd, 'burnin', 100, ...
        'thin', 100);
    fprintf('Acceptance rate: %0.6f\n', acc);
end
zm_smpl = zm_boot;


end

