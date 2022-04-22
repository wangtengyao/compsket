% main script to analyse the gene data given in CD4_goodTREG_in_thymus.csv

% the data is stored in 'CD4_goodTREG_in_thymus.mat' and 'CTL4A.mat'

n1 = 400;
n2 = 420;
p = 600;
k = 100;
rho = 2;
sigma = 1;


[X1, X2, y1, y2] = simulate(n1,n2,p,k,rho,sigma);
tic;
[stat, result] = complementarySketching(X1,X2,y1,y2,1,true);
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%% analyzing the single-cell data %%%%%%%%%%%%%%%%
% CD4_goodTREG_in_thymus = readtable('CD4_goodTREG_in_thymus.csv');
tic;
load("CD4_goodTREG_in_thymus.mat")
gene_names = dat.Properties.VariableNames(3:end);
X = table2array(dat(:, 3:end));
CD4_filter = ismember(dat.anno_lvl_2_final_clean, 'CD4+T');
X1 = X(CD4_filter, :);
X2 = X(~CD4_filter, :);
toc;


%%%%%%%%%%%%%%%%%%%%%%%%% trial with CTL4A gene %%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
j = find(strcmp(gene_names, 'CTLA4'));
X1mj = X1(:, [1:j-1 j+1:end]); X2mj = X2(:, [1:j-1 j+1:end]); %
X1j = X1(:,j); X2j = X2(:,j);
toc;
tic;
complementarySketching(X1mj, X2mj, X1j, X2j,nan,true)
toc;





%%%%%%%%%%%%%%%%%%%%%%% comprehensive analysis to the whole dataset (or its
%%%%%%%%%%%%%%%%%%%%%%% subset)
tic;
nodes = [0, 133, 180] + 1;
% nodes = 2822;
% nodes = [133] +1;
% nodes = 100:200;
% nodes = 1: size(X2,2);
% You may turn on paralle computing instrument by edit the function
% differentialNetworkAnalysis.m by changing parfor to for. 
[nodes,test_stat, test_result,interacting_partners] = differentialNetworkAnalysis(X1,X2,8,nodes,true,true);
% [nodes,test_stat, test_result,interacting_partners] = differentialNetworkAnalysisIntercerpt(X1,X2,8,nodes,true,true);
toc;
test_result_filter = ~test_result==0;
nodes_identified = nodes(test_result_filter);
interacting_partners_identified = interacting_partners(test_result_filter,:);
genes_identified = gene_names(nodes_identified);
genes_interacting_identified = gene_names(interacting_partners_identified);

pvals = ones(size(nodes_identified,2),8) * (-1);

for i = 1:size(nodes_identified,2)
    j = nodes_identified(i); 
    pvals(i) = ranksum(X1(:,j), X2(:,j));
end

% nodes_already_identified = nodes_identified;

        
