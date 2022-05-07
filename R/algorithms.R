################ Implementation for Complementary Sketching ############
library(glmnet)

#' Computing the orthonormal matrix spanning the orthogonal complement of the column span of X
#' @param X an n x p matrix with n > p
#' @return matrix A of size n x (n-p) with orthonormal columns
#' @description it seems that QR decomposition in R is slow when the dimension is large e.g. n and p > 5000; when working with large matrices, we recommend use either the python or matlab version of this package.
#' @export
orthogonalProjection <- function(X){
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
    tmp <- matrix(rnorm(n*m), n, m)
    tmp <- tmp - X%*%(solve(t(X)%*%X)%*%(t(X)%*%tmp))
    return(qr.Q(qr(tmp)))
}

#' Computing the noise variance in a high-dimensional linear model
#' @param W design matrix for high-dimensional linear model
#' @param z response in a high-dimensional linear model
#' @return a nonnegative scalar of estimated noise standard deviation
#' @description Assume data are generated from the model z = W beta + e, where e has independent components N(0, sigma^2), we estimate sigma^2 by fitting a cross-validated Lasso estimator of beta and then estimate sigma from the residual sum of squares. 
#' @export
noise_sd <- function(W, z){
    fit.cv <- glmnet::cv.glmnet(W,z,intercept=F,standardize=F)
    betahat <- coef(fit.cv, s='lambda.min')
    Xbetahat <- predict(fit.cv, newx=W, s='lambda.min')
    sqrt(sum((z - Xbetahat)^2) / (nrow(W) - sum(betahat!=0)))
}

#' Computing the noise variance in a high-dimensional linear model
#' @param W design matrix for high-dimensional linear model
#' @param z response in a high-dimensional linear model
#' @return a nonnegative scalar of estimated noise standard deviation
#' @description Assume data are generated from the model z = W beta + e, where e has independent components N(0, sigma^2), we estimate sigma^2 by using the method in Dicker (2014).
#' @references Dicker, L. H. (2014). Variance estimation in high-dimensional linear models. Biometrika, 101(2), 269-284. 
#' @export
dicker_noise_sd <- function(W, z)
{
    m <- nrow(W); p <- ncol(W)
    gram.W.norm <- t(W) %*% W/m
    m1.hat <- sum(diag(gram.W.norm))/p
    m2.hat <- sum(gram.W.norm^2)/p - p * m1.hat^2 / m
    sigma.tilde.square <- (1 + p*m1.hat^2/(m+1)/m2.hat)*sum(z^2)/m - m1.hat*sum((t(W) %*% z)^2)/m/(m+1)/m2.hat
     if (sigma.tilde.square > 0){
         return(sqrt(sigma.tilde.square))
     }
     else return(noise_sd(W,z))
    return(sqrt(sigma.tilde.square))
}

#' Main function implementing the complementary sketching algorithm in Gao and Wang (2020)
#' @param X1 design matrix of the first linear regression model
#' @param X2 design matrix of the second linear regression model
#' @param y1 response vector of the first linear regression model
#' @param y2 response vector of the second linear regression model
#' @param sparse whether to test a sparse difference in the regression coefficients of the two models, default is TRUE
#' @param sigma noise standard deviation of the regression models; if unknown, set to NULL and will be estimated.
#' @return a list containing the test statistics and test results (0 for nonrejection and 1 for rejection).
#' @description From the model y1 = X1 beta1 + eps1 and y2 = X2 beta2 + eps2, where X1 and X2 are n1 x p and n2 x p design matrices respectively, we test the null H0: beta1 = beta2 against the alternative that beta1 - beta2 is nonzero (if sparse=FALSE) or beta1 - beta2 is nonzero and has at most sqrt(p) (if sparse=TRUE). The test is performed using the complementary sketching algorithm from Gao and Wang (2020).
#' @references Gao, F. and Wang, T. (2020) Two-sample testing of high-dimensional linear regression coefficients via complementary sketching. Preprint, arxiv:2011.13624.
#' @export
#' @examples 
#' # problem parameters
#' n1 <- 600 # sample size of first sample
#' n2 <- 600 # sample size of second sample
#' p <- 800 # dimension of covariates
#' k <- 10 # sparsity of difference of the two regression coefficients
#' rho <- 1 # difference in l_2 norm of the two regression coefficients
#' # generate design matrices
#' X1 <- matrix(rnorm(n1 * p), n1, p)
#' X2 <- matrix(rnorm(n2 * p), n2, p)
#' # generate regression coefficients
#' beta1 <- rnorm(p)
#' theta <- c(rnorm(k), rep(0, p-k)) 
#' theta <- theta / sqrt(sum(theta^2)) * rho
#' beta2 <- beta1 + theta
#' # generate response vectors
#' y1 <- X1 %*% beta1 + rnorm(n1)
#' y2 <- X2 %*% beta2 + rnorm(n2)
#' # test for difference in beta1 and beta2
#' complementarySketching(X1,X2,y1,y2)
complementarySketching <- function(X1, X2, y1, y2, sparse=TRUE, sigma=NULL){
    X <- rbind(X1, X2); y <- c(y1, y2)
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p; n1 <- dim(X1)[1];

    A <- orthogonalProjection(X)
    B <- A; B[(n1+1):n, ] <- -A[(n1+1):n, ]

    W <- t(B)%*%X; z <- t(A)%*%y
    if (is.null(sigma)) sigma <- dicker_noise_sd(W,z)
    if (sigma==0) return(list(stat=NA, result=NA))
    W <- W/sigma; z <- z/sigma
    Wtilde <- sweep(W, 2, sqrt(colSums(W^2)), '/')

    if (sparse) {
        lambda <- sqrt(4*log(p))
        tau <- 3*log(p)
        Q <- t(Wtilde) %*% z
        Q_thresh <- Q * (abs(Q) > lambda)
        test.stat <- sum(Q_thresh^2)
        test.result <- ifelse(test.stat > tau, 1, 0)
    } else {
        tau <- m + sqrt(8*m*log(p)) + 4*log(p)
        test.stat <- sum(z^2)
        test.result <- ifelse(test.stat > tau, 1, 0)
    }
    return(list(stat=test.stat, result=test.result))
}

