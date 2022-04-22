function [nodes,test_stat,test_result, interacting_partners] = differentialNetworkAnalysisIntercerpt(X1,X2,num_partners,nodes,sparse,trace)
%DIFFERENTIALNETWORKANALYSIS Differential network testing for two covariate matrices
% Different from the version without intercept, an extra constant column is
% added to the design as the intercept term and all the origianl coolumns
% are centered. 
% Differential network testing for two covariate matrices X1 (shape n1 x p) and X2 (shape n2 x p)
% This has essentially the same effect as applying complementarySketching() p times comparing
% nodewise regressions X1[:, j] ~ X1[:, -j] and X2[:, j] ~ X2[:, -j], though code is optimised
% to compute the sketching matrix only once.
% It also can output top num_partners interacting partners of significant nodes found. The nodes
% parameter allows for faster computation if only a subset of nodes is interested. Set trace=True
% to output computational progress.


addpath('./imm3897/') % for lars implemented by Karl Sjöstrand
% addpath('./lars/') % for lars implemented by Sung Soo Kim

%% centre the columns of the design and add an extra constant column as intercept 

X = [X1;X2]; 
n = size(X,1);
X = X - repmat(sum(X)/n,n,1);
X = [X ones(n,1)];

p = size(X,2); % note that p-1 is the number of total genes
n1 = size(X1,1);



%% precompute a sketching matrix; the sketching matrix for individual nodes can be obtained
% by augmenting A0 with a single column
if trace
    disp("'Computing complementary sketches...'")
end

%tic;
[Q, R] = qr(X);
A0 = Q(:,p+1:n);
B0 = Q(:,1:p);
R0 = R(1:p,:);
%disp("First time QR decomposition")
%toc;

if isnan(nodes)
    nodes = (1:(p-1)).';
end
if isrow(nodes)
    nodes = nodes.';
end

test_stat = ones(size(nodes,1),1) * (-1);
test_result = ones(size(nodes,1),1)*(-1);
interacting_partners = ones(size(nodes,1),num_partners)*(-1);

no_nodes = size(nodes,1);

for i = 1:no_nodes
% parfor (i= 1:no_nodes)
    %tic;
    j = nodes(i);
    X1mj = X(1:n1, [1:j-1 j+1:end]); X2mj = X(n1+1:n, [1:j-1 j+1:end]); %
    X1j = X(1:n1,j); X2j = X(n1+1:n,j);
    Xj = [X1j; X2j];
    if trace
        fprintf("Computing %d / %d node...\n",i,no_nodes)
    end
    %disp('obtaining node-wide design')
    %toc;

    %% update the sketching matrix for each node
    %tic;
    [Q1, ~] = qrdelete(B0,R0,j,'col');
    %disp('QR deletion')
    %toc;
    %tic;
    u = Xj - Q1 * (Q1.' * Xj);
    w = u / norm(u);
    A = [A0 w];
    %disp('projection')
    %toc;


    %% construct complementarily sketched design W and response z
    % tic;
    A1 = A(1:n1,:).';
    A2 = A(n1+1:n,:).';
    W = A1 * X1mj - A2 * X2mj;
    z = A1 * X1j + A2 * X2j;
    % disp('sketching')
    % toc;


    %% compute sigma, the noise variance, using W and z
    sigma = dickerNoiseVar(W,z);
    if ~isreal(sigma)
        % df(i,:) = {j, 0, 0, {[]}}; 
        test_stat(i) = 0;
        test_result(i) = 0;
        continue
    end

    
    %% standardize W and z
    W = W / sigma; z = z / sigma;
    W_tilde = W ./ repmat(sqrt(sum(W.^2)),size(W,1),1);
    % W_tilde = normc(W);


    %% compute test statistics based on sparsity knoweldge
    if sparse
        lam = sqrt(4*log(p-1));
        tau = 3*log(p-1);
        Q = W_tilde.' * z;
        % Q_thresh = wthresh(Q,'h',lam);
        Q_thresh = Q .* double(abs(Q) > lam);
        test_stat(i) = sum(Q_thresh.^2);
    else
        tau = m + sqrt(8*log(p-1)) + 4*log(p);
        test_stat(i) = sum(z.^2);
    end

    test_result(i) = double(test_stat(i) > tau);



    %% compute interacting partners using Least Angle Regression of z 
    % against W
    %tic;
    if test_stat(i) && num_partners>0
        other_nodes = [1:j-1 j+1:p];
        b_lars = lars(W,z,'lars',-num_partners,0,[],0);
        active_filter = ~(b_lars(end,:)==0);
        interacting_partners(i,:) = other_nodes(active_filter);
    end
    %disp('lars')
    %toc;
end

    
if trace
    fprintf(" ")
    fprintf("Finished: %d significant nodes found.\n", sum(test_result))
end

end




