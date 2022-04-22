import numpy as np
import scipy as sp
import pandas as pd
from sklearn import linear_model
import sklearn
import time
from datetime import datetime
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
nodes = [1, 133, 180]
# nodes = range(100,200)
# nodes = None
result = differentialNetworkAnalysisIntercept(X1, X2, nodes=nodes, num_partners=8, intercept=False, trace=True)
#result = differentialNetworkAnalysis(X1, X2, nodes=nodes, num_partners=8, trace=True)
print('{} seconds elapsed.'.format(time.time() - t))

"""
save results
"""
# print(result)
result_identified = result.loc[result['test_result'] == 1, :]
result_identified.index = np.array(gene_names)[result_identified.index]
result_identified.insert(len(result_identified.columns), 'interacting_genes', [gene_names[v] for v in result_identified['interacting_partners']])

def DGE(x1, x2, alpha):
    pvals = [sp.stats.mannwhitneyu(x1[:,j], x2[:,j])[1] for j in range(x1.shape[1])]
    return np.where(np.array(pvals) < alpha, 1, 0)

result_identified.insert(len(result_identified.columns),'DGE', [DGE(X1[:,v], X2[:,v], 0.05 / X1.shape[1]) for v in result_identified['interacting_partners']])
result_identified.pop('interacting_partners');
date = datetime.now().strftime("%Y%m%d-%H%M%S")
csv_filename = f"result_{date}.csv"
result_identified.to_csv(csv_filename, index=True, header=True, sep=',')

