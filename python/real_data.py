import numpy as np
import scipy as sp
import pandas as pd
from sklearn import linear_model
import sklearn
import time
from compsket import *


"""
read in data
"""
print('read in data...')
dat = pd.read_csv('CD4_TREG_in_thymus.csv.gz', index_col=0)
X = np.array(dat.iloc[:,1:])
gene_names = np.array(dat.columns[1:])
CD4_filter = dat['anno_lvl_2_final_clean'] == 'CD4+T'

"""
Perform differential network analysis
"""
X1 = np.array(X[CD4_filter, :])
X2 = np.array(X[~CD4_filter, :])
#X1 = X1 - np.average(X1, axis=0)
#X2 = X2 - np.average(X2, axis=0)
t = time.time()
#nodes = [1, 133, 180]
nodes = None
result = differentialNetworkAnalysis(X1, X2, nodes=nodes, num_partners=8, trace=True)
print('{} seconds elapsed.'.format(time.time() - t))

"""
save results
"""
print(result)
result = result.loc[result['test_result'] == 1, :]
result.index = gene_names[result.index]
result['interacting_genes'] = [gene_names[v] for v in result['interacting_partners']]

def DGE(x1, x2, alpha):
    pvals = [sp.stats.mannwhitneyu(x1[:,j], x2[:,j])[1] for j in range(x1.shape[1])]
    return np.where(np.array(pvals) < alpha, 1, 0)
        
result['DGE'] = [DGE(X1[:,v], X2[:,v], 0.05 / X1.shape[1]) for v in result['interacting_partners']]
result.pop('interacting_partners');
result.to_csv('result.csv', index=True, header=True, sep=',')

