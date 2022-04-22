function [X1, X2, y1, y2] = simulate(n1, n2, p, k, rho, sigma)
%SIMULATE Summary of this function goes here
%   Detailed explanation goes here

X1 = randn(n1, p);
X2 = randn(n2, p);

beta1 = randn(p,1);
theta = [ones(k,1) .* ((rand(k,1)<.5)*2 - 1) ;zeros(p-k,1)];
theta = theta / norm(theta) * rho;
beta2 = beta1 + theta;

y1 = X1 * beta1 + randn(n1,1) * sigma;
y2 = X2 * beta2 + randn(n2,1) * sigma;

end

