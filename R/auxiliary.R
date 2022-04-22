ql <- function(A){
    B <- A[,ncol(A):1]
    tmp <- qr(B)
    Q <- qr.Q(tmp)
    R <- qr.R(tmp)
    Q <- Q[,ncol(Q):1]
    L <- R[nrow(R):1,ncol(R):1]
    return(list(Q=Q, L=L))
}

sqrtm <- function(A){
    tmp <- svd(A)
    tmp$u%*%diag(sqrt(tmp$d), nrow=nrow(A))%*%t(tmp$v)
}

standardise <- function(X){
    X <- sweep(X, 2, colMeans(X))
    X <- apply(X, 2, vector.normalise) * sqrt(nrow(X))
    return(X)
}

