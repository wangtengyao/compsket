import numpy as np
import scipy as sp
import pandas as pd
from sklearn import linear_model

"""
Compute the noise standard deviation in a sparse linear regression model with response z and covariates W
"""
def dickerNoiseSD(W, z):
    m, p = W.shape
    W_Gram = np.matmul(W.T, W) / m
    m1hat = np.sum(np.diagonal(W_Gram)) / p
    m2hat = np.sum(W_Gram**2)/p - p * (m1hat ** 2) / m
    V = (1 + p * m1hat**2 / (m+1) / m2hat) * np.sum(np.power(z, 2)) / m  - m1hat * np.sum(np.power(np.matmul(W.T, z), 2)) / m / (m+1) / m2hat
    return np.sqrt(V)

"""
Testing for equality of high-dimensional regression coefficients via complementary sketching
Model: y1 = X1 beta1 + eps1; y2 = X2 beta2 + eps2
Test: H0: beta1 = beta2; H1: beta1 - beta2 is sparse (if sparse=True) or beta1 != beta2 (if sparse=False)
If noise variance is known, supply it in the sigma parameter, otherwise, noise variance will be computed
via dickerNoiseSD()
"""
def complementarySketching(X1, X2, y1, y2, sparse=True, sigma=None):
    # stack covariates and responses into X and y and get shape parameters
    X = np.concatenate((X1, X2), axis=0)
    y = np.concatenate((y1, y2))
    n, p = X.shape
    m = n - p
    n1, _ = X1.shape
    
    # compute the sketching matrix A
    Q, _ = np.linalg.qr(X, mode='complete')
    A = Q[:, range(p, n)].T
    A1 = A[:, range(n1)]
    A2 = A[:, range(n1, n)]
    
    # construct complementarily sketched design W and response z
    W = np.matmul(A1, X1) - np.matmul(A2, X2)
    z = np.matmul(A1, y1) + np.matmul(A2, y2)
    
    # compute sigma, the noise variance, using W and z if necessary
    if sigma is None:
        sigma = dickerNoiseSD(W, z)

    if sigma <= 0:
        return (np.nan, np.nan)
    
    # standardise W and z
    W = W / sigma; z = z / sigma
    Wtilde = W / np.sqrt(np.sum(W ** 2, axis=0))
    
    # compute test statistics based on sparsity knowledge
    # sparse is true if the number of nonzero entries in the difference of two regression coefficients is smaller than square root of p, with p as the dimension of coefficients
    if sparse:
        lam = np.sqrt(4 * np.log(p))
        tau = 3 * np.log(p)
        Q = np.matmul(Wtilde.T, z)
        Q_thresh = Q * (np.abs(Q) > lam)
        test_stat = np.sum(Q_thresh ** 2)
    else:
        tau = m + np.sqrt(8 * np.log(p)) + 4 * np.log(p)
        test_stat = np.sum(z ** 2)
    
    test_result = 1 if test_stat > tau else 0
    return (test_stat, test_result)


"""
Differential network testing for two covariate matrices X1 (shape n1 x p) and X2 (shape n2 x p)
This has essentially the same effect as applying complementarySketching() p times comparing 
nodewise regressions X1[:, j] ~ X1[:, -j] and X2[:, j] ~ X2[:, -j], though code is optimised
to compute the sketching matrix only once. 

It also can output top num_partners interacting partners of significant nodes found. The nodes
parameter allows for faster computation if only a subset of nodes is interested. Set trace=True
to output computational progress.
"""

def differentialNetworkAnalysis(X1, X2, num_partners=0, nodes=None, sparse=True, trace=False):
    # stack covariates and responses into X and y and get shape parameters
    X = np.concatenate((X1, X2), axis=0)
    n, p = X.shape
    n1, _ = X1.shape
    
    # precompute a sketching matrix; the sketching matrix for individual nodes can be obtained
    # by augmenting A0 with a single column
    if trace: 
        print('Computing complementary sketches...')
    Q, R = np.linalg.qr(X, mode='complete')
    A0 = Q[:, range(p, n)]
    B0 = Q[:, range(p)]
    R0 = R[range(p),:]
    
    # results will be stored in df
    if nodes is None: 
        nodes = range(p)
    nodes = list(nodes)
    df = pd.DataFrame(index=nodes, columns=['test_stat', 'test_result', 'interacting_partners'])
    
    for j in nodes:
        i = nodes.index(j)
        if trace: 
            print('Computing {current}/{total} node...'.format(current=i+1, total=len(nodes)), end='\r')
            
        # update the sketching matrix for each node
        Q1, R1 = sp.linalg.qr_delete(B0, R0, k=j, which='col')
        u = X[:, j] - np.matmul(Q1, np.matmul(Q1.T, X[:, j]))
        w = u / np.linalg.norm(u)
        A = np.concatenate((A0, w.reshape((-1,1))), axis=1)
        
        # construct complementarily sketched design W and response z
        A1 = A[range(n1), :].T
        A2 = A[range(n1, n), :].T
        W = np.matmul(A1, np.delete(X1, j, axis=1)) - np.matmul(A2, np.delete(X2, j, axis=1))
        z = np.matmul(A1, X1[:, j]) + np.matmul(A2, X2[:, j])
        
        # compute sigma, the noise variance, using W and z
        sigma = dickerNoiseSD(W, z)
        if sigma <= 0:
            df.iloc[i, :] = [0, 0, '']
            continue
        
        # standardise W and z
        W = W / sigma; z = z / sigma
        Wtilde = W / np.sqrt(np.sum(W ** 2, axis=0))
        
        # compute test statistics based on sparsity knoweldge
        if sparse:
            lam = np.sqrt(4 * np.log(p-1))
            tau = 3 * np.log(p-1)
            Q = np.matmul(Wtilde.T, z)
            Q_thresh = Q * (np.abs(Q) > lam)
            test_stat = np.sum(Q_thresh ** 2)
        else:
            tau = m + np.sqrt(8 * np.log(p)) + 4 * np.log(p)
            test_stat = np.sum(z ** 2)

        test_result = 1 if test_stat > tau else 0
        df.iloc[i, :] = [test_stat, test_result, '']
        
        # compute interacting partners using Least Angle Regression of z against W
        if test_result==1 and num_partners>0:
            other_nodes = np.delete(range(p), j)
            lars = linear_model.Lars(n_nonzero_coefs=num_partners, fit_intercept=True, normalize=False)
            lars.fit(W, z)
            df.iloc[i, 2] = other_nodes[lars.active_]
        
    if trace:
        print()
        print('Finished: {} significant node found.'.format(np.sum(df['test_result'])))
        
    return df


"""
Differential network testing for two covariate matrices X1 (shape n1 x p) and X2 (shape n2 x p)
This has essentially the same effect as applying complementarySketching() p times comparing 
nodewise regressions X1[:, j] ~ X1[:, -j] and X2[:, j] ~ X2[:, -j], though code is optimised
to compute the sketching matrix only once. 

It also can output top num_partners interacting partners of significant nodes found. The nodes
parameter allows for faster computation if only a subset of nodes is interested. Set trace=True
to output computational progress.

The difference between this one and the one without Intercept is that it allows an extra intercept term to be added to the design such that centerisation of columns in the design could make sense.
"""


def differentialNetworkAnalysisIntercept(X1, X2, num_partners=0, nodes=None, sparse=True, intercept= True, trace=False):
    # stack covariates and responses into X and y and get shape parameters
    X = np.concatenate((X1, X2), axis=0)
    n, _ = X.shape
    n1, _ = X1.shape
    if intercept: # if intercept term is necessary, centre X's columns first, and then add an extra columns of ones
        X = X - X.sum(axis= 0)/n
        X = np.concatenate((X, np.ones((n,1))), axis=1)
    _, p = X.shape

    
    # precompute a sketching matrix; the sketching matrix for individual nodes can be obtained
    # by augmenting A0 with a single column
    if trace: 
        print('Computing complementary sketches...')
    Q, R = np.linalg.qr(X, mode='complete')
    A0 = Q[:, range(p, n)]
    B0 = Q[:, range(p)]
    R0 = R[range(p),:]
    
    # results will be stored in df
    if nodes is None:
        if intercept:
            nodes = range(p-1)
        else:
            nodes = range(p)
    nodes = list(nodes)
    df = pd.DataFrame(index=nodes, columns=['test_stat', 'test_result', 'interacting_partners'])
    
    for j in nodes:
        i = nodes.index(j)
        if trace: 
            print('Computing {current}/{total} node...'.format(current=i+1, total=len(nodes)), end='\r')
            
        # update the sketching matrix for each node
        Q1, R1 = sp.linalg.qr_delete(B0, R0, k=j, which='col')
        u = X[:, j] - np.matmul(Q1, np.matmul(Q1.T, X[:, j]))
        w = u / np.linalg.norm(u)
        A = np.concatenate((A0, w.reshape((-1,1))), axis=1)
        
        # construct complementarily sketched design W and response z
        A1 = A[range(n1), :].T
        A2 = A[range(n1, n), :].T
        W = np.matmul(A1, np.delete(X1, j, axis=1)) - np.matmul(A2, np.delete(X2, j, axis=1))
        z = np.matmul(A1, X1[:, j]) + np.matmul(A2, X2[:, j])
        
        # compute sigma, the noise variance, using W and z
        sigma = dickerNoiseSD(W, z)
        if sigma <= 0:
            df.iloc[i, :] = [0, 0, '']
            continue
        
        # standardise W and z
        W = W / sigma; z = z / sigma
        Wtilde = W / np.sqrt(np.sum(W ** 2, axis=0))
        
        # compute test statistics based on sparsity knoweldge
        if sparse:
            lam = np.sqrt(4 * np.log(p-1))
            tau = 3 * np.log(p-1)
            Q = np.matmul(Wtilde.T, z)
            Q_thresh = Q * (np.abs(Q) > lam)
            test_stat = np.sum(Q_thresh ** 2)
        else:
            tau = m + sqrt(8 * np.log(p)) + 4 * np.log(p)
            test_stat = np.sum(z ** 2)

        test_result = 1 if test_stat > tau else 0
        df.iloc[i, :] = [test_stat, test_result, '']
        
        # compute interacting partners using Least Angle Regression of z against W
        if test_result==1 and num_partners>0:
            other_nodes = np.delete(range(p), j)
            lars = linear_model.Lars(n_nonzero_coefs=num_partners, fit_intercept=True, normalize=False)
            lars.fit(W, z)
            df.iloc[i, 2] = other_nodes[lars.active_]
        
    if trace:
        print()
        print('Finished: {} significant node found.'.format(np.sum(df['test_result'])))
        
    return df
