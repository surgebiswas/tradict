function [ xs, mu, sig] = standardize( x, varargin)
% z-score transformation of x colum wise

centeronly = setParam(varargin, 'centeronly', false);
mu = setParam(varargin, 'mu', mean(x));
sig = setParam(varargin, 'std', std(x));

xs = bsxfun(@minus, x, mu);
if ~centeronly
    xs = bsxfun(@rdivide, xs, sig);
end

q = any(isnan(xs));
xs(:,q) = 0; 

end

