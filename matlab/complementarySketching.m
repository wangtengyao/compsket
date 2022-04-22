function [test_stat, test_result] = complementarySketching(X1, X2, y1, y2, sigma, sparse)
%COMPLEMENTARYSKETCHING Testing for equality of high-dimensional regression coefficients via complementary sketching
%   X1, y1 --- Design 1, response 1
%   X2, y2 --- Design 2, response 2
%   Model: y1 = X1 beta1 + eps1; y2 = X2 beta2 + eps2
%   Test: H0: beta1 = beta2; H1: beta1 - beta2 is sparse (if sparse=True) or beta1 != beta2 (if sparse=False)
%   Provide a sigma if you know it, otherwise it will be calculated by
%   function dickerNoiseVar.
%   Give a logical value to sparse, true or false, which corresponds to the
%   sparse or dense tests as in Gao and Wang 2020
%   (https://arxiv.org/abs/2011.13624)Testing for equality of high-dimensional regression coefficients via complementary sketching




n1 = size(X1, 1); n2 = size(X2,1); p = size(X1,2); n = n1+n2; m = n-p;
X = [X1; X2]; y = [y1; y2];



[Q,~] = qr(X);
A = Q(:,p+1:n).';
B = A; B(:,n1+1:n) = -B(:,n1+1:n);
W = B * X; z = A * y;

if isnan(sigma)
    sigma = dickerNoiseVar(W,z);
    disp(sigma)
end
W = W / sigma; z = z / sigma;
W_tilde = W ./ repmat(sqrt(sum(W.^2)),size(W,1),1);

if sparse
    lambda = sqrt(4*log(p));
    tau = 3 * log(p);
    Q = W_tilde.' * z;
    % Q_thresh = wthresh(Q,'h',lambda);
    Q_thresh = Q .* double(Q>lambda);
    test_stat = sum(Q_thresh.^2);
    test_result = test_stat > tau;
else 
    tau = m + sqrt(8*m*log(p)) + 4 * log(p);
    test_stat = sum(z.^2);
    test_result = test_stat > tau;
end

end

