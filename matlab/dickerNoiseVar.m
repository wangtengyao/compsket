function [sigma] = dickerNoiseVar(W,z)
%DICKERNOISEVAR To estimate the noise variance in a standard
%regression model.
%   W -- regression design   z -- response
%   The method is from Dicker 2014.
%   (https://academic.oup.com/biomet/article-abstract/101/2/269/194884)

m = size(W,1);
p = size(W,2);

gram_W_norm = W.' * W / m;
m1_hat = trace(gram_W_norm)/p;
m2_hat = sum(gram_W_norm.^2,'all')/p - (p/m) * m1_hat^2;
sigma_tilde_square = (1+p*m1_hat^2/(m+1)/m2_hat) *sum(z.^2)/m - m1_hat * sum((W.' * z).^2)/m/(m+1)/m2_hat;

sigma = sqrt(sigma_tilde_square);

end

