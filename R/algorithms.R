library(putils)
library(fastclime)
library(MASS)
library(lars)
library(glmnet)

##### likelihood ratio test #####
LRT <- function(X1, X2, y1, y2){
    X <- rbind(X1, X2); y <- c(y1, y2)
    n <- nrow(X); p <- ncol(X)

    lm1 <- lm(y1 ~ X1 - 1)
    lm2 <- lm(y2 ~ X2 - 1)
    lm <- lm(y ~ X - 1)

    RSS_null <- deviance(lm); RSS_alt <- deviance(lm1) + deviance(lm2)
    stat <- ((RSS_null - RSS_alt) / p) / (RSS_alt / (n - p * 2))
    pval <- 1 - pf(stat, p, n - p * 2)
    return(list(pval=pval, result=ifelse(pval<0.05, 1, 0)))
}

##### our test #####
orthogonalProjection <- function(X){
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
    # tmp <- cbind(X, matrix(rnorm(n*m), n, m))
    # tmp <- qr.Q(qr(tmp))
    # return(t(tmp[, (p+1):n]))

    tmp <- matrix(rnorm(n*m), n, m)
    tmp <- tmp - X%*%(solve(t(X)%*%X)%*%(t(X)%*%tmp))
    return(t(qr.Q(qr(tmp))))
}


orthogonalProjection.qr <- function(X){
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
    # tmp <- cbind(X, matrix(rnorm(n*m), n, m))
    # tmp <- qr.Q(qr(tmp))
    # return(t(tmp[, (p+1):n]))
    Xc <- matrix(rnorm(m*n), nrow = n)
    return(t(qr.Q(qr(cbind(X,Xc)))[,(p+1):n]))
}

noise_var <- function(W, z){
    #mad(z)
    fit.lasso <- glmnet(W,z,intercept=F,standardize=F)
    fit.cv <- cv.glmnet(W,z,intercept=F,standardize=F)
    betahat <- coef(fit.cv, s='lambda.min')
    Xbetahat <- predict(fit.cv, newx=W, s='lambda.min')
    sqrt(sum((z - Xbetahat)^2) / (nrow(W) - sum(betahat!=0)))
}

dicker.noise.var <- function(W, z)
{
    m <- nrow(W); p <- ncol(W)
    gram.W.norm <- t(W) %*% W/m
    m1.hat <- sum(diag(gram.W.norm))/p
    m2.hat <- sum(diag( gram.W.norm %*% gram.W.norm))/p - p * m1.hat^2 / m
    sigma.tilde.square <- (1 + p*m1.hat^2/(m+1)/m2.hat)*sum(z^2)/m - m1.hat*sum((t(W) %*% z)^2)/m/(m+1)/m2.hat
     if (sigma.tilde.square > 0){
         return(sqrt(sigma.tilde.square))
     }
     else return(noise_var(W,z))
    return(sqrt(sigma.tilde.square))
}

# gaussian.second.truncated.moment <- function(u)
# {
#      # for alpha_s in Carpentier et al
#      integrand <- function(x) {x^2 * exp(-x^2/2)/sqrt(2*pi)}
#      return(exp(log(integrate(integrand, lower = u , upper = Inf)$value) -  log( pnorm(u, lower.tail = F))))
# }

gaussian.second.truncated.moment <- function(u){
    1 + u * dnorm(u) / pnorm(u, lower.tail = F)
}

complementarySketching <- function(X1, X2, y1, y2, sparse=T, sigma=NULL){
    X <- rbind(X1, X2); y <- c(y1, y2)
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p; n1 <- dim(X1)[1];

    A <- orthogonalProjection(X)
    B <- A; B[,(n1+1):n] <- -A[,(n1+1):n]

    W <- B%*%X; z <- A%*%y
    if (is.null(sigma)) sigma <- dicker.noise.var(W,z)
    if (sigma==0) return(list(stat=NA, result=NA))
    W <- W/sigma; z <- z/sigma
    Wtilde <- sweep(W, 2, sqrt(colSums(W^2)), '/')

    if (sparse) {
        lambda <- sqrt(4*log(p))
        tau <- 3*log(p)
        Q <- t(Wtilde) %*% z
        Q_thresh <- vector.hard.thresh(Q, lambda)
        test.stat <- sum(Q_thresh^2)
        test.result <- ifelse(test.stat > tau, 1, 0)
    } else {
        tau <- m + sqrt(8*m*log(p)) + 4*log(p)
        test.stat <- sum(z^2)
        test.result <- ifelse(test.stat > tau, 1, 0)
    }
    return(list(stat=test.stat, result=test.result))
}

##### Method of Chabonnier et al #####
# inverting x+1/x-2=a for x > 1
ginv <- function(a){
    x <- max(1/a, a)
    repeat{
        g <- x+1/x-2-a
        x <- x - g / (1-1/x^2)
        if (g < 1e-10) break
    }
    return(x)
}

# Eq (15)
Qtilde <- function(u, a, n1, s){
    a1 <- sum(a); aoo <- max(a); a2 <- sum(a^2)
    b <- a1*u/(aoo*(n1-s)) + u + a2 / aoo - a1
    Delta <- b^2 - 4*u*(u-a1)/(n1-s)/aoo*(a1-a2/aoo)
    lambda_star <- (b-sqrt(Delta)) / (4*u/(n1-s)*(a1-a2/aoo))
    exp(-sum(log(1-lambda_star*a)) / 2 - (n1-s)/2*log(1+2*lambda_star*u/(n1-s)))
}

Chabonnier <- function(X1, X2, y1, y2){
    W <- rbind(cbind(X1, X1), cbind(X2, -X2))
    y <- c(y1, y2)
    ret <- lars::lars(W, y, normalize=F, intercept=F)
    action <- unlist(ret$actions)
    #plot(cumsum(abs(action <= k)) / cumsum(sign(action)), xlab='actions',
    #     ylab='prop signal vars', type='l')
    S <- numeric(0)
    n1 <- nrow(X1); n2 <- nrow(X2)
    D <- min(n1, n2) %/% 2
    p <- ncol(X1)
    test.result <- 0
    L <- match(D, cumsum(sign(action)), nomatch = length(action))
    support <- rep(0, p)

    for (i in seq_along(action)){
        v <- action[i]
        if (v > p) v <- v - p
        if (v < -p) v <- v + p
        if (v>0) {
            support[v] <- support[v] + 1
        } else {
            support[-v] <- support[-v] - 1
        }
        S <- which(support > 0)
        s <- length(S)
        if (s > D) break
        X1S <- X1[,S,drop=F]
        X2S <- X2[,S,drop=F]

        lm1 <- lm(y1~X1S-1)
        lm2 <- lm(y2~X2S-1)
        F_SV <- -2 + (deviance(lm1)/n1)/(deviance(lm2)/n2) +
            (deviance(lm2)/n2)/(deviance(lm1)/n1) # Eq (5)
        F_S1 <- sum((X2S%*%(lm1$coefficients - lm2$coefficients))^2) / n2 /
            (deviance(lm1)/n1) # Eq (6)
        F_S2 <- sum((X1S%*%(lm1$coefficients - lm2$coefficients))^2) / n1 /
            (deviance(lm2)/n2) # Eq (6)

        q_SV <- pf(ginv(F_SV)*n1*(n2-s)/n2/(n1-s), n1-s, n2-s) +
            pf(ginv(F_SV)*n2*(n1-s)/n1/(n2-s), n2-s, n1-s) # Eq (14)

        tmp <- n1/n2/(n1-s)*X2S%*%(solve(t(X1S)%*%X1S) + solve(t(X2S)%*%X2S))%*%t(X2S) # Defn 3.2
        a <- eigen(tmp, symmetric=T, only.values=T)$values[1:s] # Defn 3.2
        q_S1 <- Qtilde(F_S1, a, n1, s) # Eq (18)
        q_S2 <- Qtilde(F_S2, a, n2, s) # Eq (18)
        if (is.nan(q_S1)) q_S1 <- 1
        if (is.nan(q_S2)) q_S2 <- 1

        alpha <- 0.025 / D / choose(p, s) # Eq (20)
        if ((q_SV < alpha) || (q_S1 < alpha/2) || (q_S2 < alpha/2)) {
            test.result <- 1
            break
        }
        # printPercentage(i, length(action))
    }
    return(list(result=test.result))
}

##### Method of Zhu and Bradic; we give oracle sigma and Pi #####
# Dantzig_Selector <- function(X, y, lambda){
#     # CVXR optimisation
#     b <- CVXR::Variable(p)
#     objective <- CVXR::Minimize(CVXR::norm1(b))
#     n <- nrow(X)
#     constraint <- CVXR::norm_inf(t(X)%*%(y - X%*%b)) <= lambda * n
#     prob <- CVXR::Problem(objective, constraints=list(constraint))
#     result <- solve(prob)
#     thetahat <- result$getValue(b)
#     return(thetahat)
# }

Dantzig_Selector <- function(X, y, lambda){
    capture.output(ret <- dantzig(X, y, lambda, 100000), file='/dev/null')
    thetahat <- dantzig.selector(ret$lambdalist, ret$BETA0, lambda)
    return(thetahat)
}

ZhuBradic <- function(X1, X2, y1, y2, sigma=1, Pi=NULL){
    n <- min(nrow(X1), nrow(X2)); p <- ncol(X1)
    X1 <- X1[1:n, ]; X2 <- X2[1:n, ];
    y1 <- y1[1:n]; y2 <- y2[1:n];
    W <- X1 + X2; Z <- X1 - X2; y <- y1 + y2
    eta <- sqrt(2*log(p)) * sqrt(max(colSums(W^2))) / n # eta as in p.23

    thetahat <- Dantzig_Selector(W, y, eta * sigma * n) # (2.4), using oracle sigma

    if (is.null(Pi)){
        Pi <- matrix(0, p, p)
        for (j in 1:p){
            Pi[,j] <- Dantzig_Selector(W, W[,j], eta * sigma * n)
            # printPercentage(j,p)
        }
    }

    stat <- max(abs(t(Z-W%*%Pi)%*%(y-W%*%thetahat)))/norm(y-W%*%thetahat, 'F') # use oracle Pi, (2.10)
    Q <- t(Z - W%*%Pi) %*% (Z - W%*%Pi) / n # (2.12)
    tmp <- mvrnorm(1000, rep(0, p), Q)
    pval <- sum(stat < apply(abs(tmp), 1, max)) / 1000

    return(list(stat=stat, pval=pval))
}


Carpentier <- function(X1, X2, y1, y2, k, Cs=1, sigma=NULL){
    # must supply the sparsity parameter k, which is their s
    # Cs, C_* in their context
    X <- rbind(X1, X2); y <- c(y1, y2)
    n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p; n1 <- dim(X1)[1];

    A <- orthogonalProjection(X)
    B <- A; B[,(n1+1):n] <- -A[,(n1+1):n]

    W <- B%*%X; z <- A%*%y
    if (is.null(sigma)) sigma <- dicker.noise.var(W,z)
    if (sigma==0) return(list(stat=NA, result=NA))
    Wtilde <- sweep(W, 2, sqrt(colSums(W^2)), '/')

    Ginv <- solve(t(W) %*% W) # (X^T X)^(-1) as in Carpentier et al on page 4
    yi <- Ginv %*% t(W) %*% z # y as in Carpentier et al on page 4,
    # to differentiate from original y, so yi

    sparse <- ifelse(k<sqrt(p), T, F)
    test.size <- 0.05 # epsilon as in A_epsilon on p. 1821

    # Cs <- 1 # what the hell is C_*?
    if (sparse) { # as in s < sqrt(p) in Carpentier et al on page 4
        psi <- k * log(1+p/k/k) / m # psi(s,p) in Carpentier
        yi.sigma <- ifelse(yi^2 > 2 * sigma^2 * diag(Ginv) * log(1+p/k/k), 1,0)
        alphas <- gaussian.second.truncated.moment(sqrt(2*log(1 + p/k/k)))
        Qhat <- sum((yi^2 - sigma^2 * diag(Ginv) * alphas)*yi.sigma)

    } else { # as in s >= sqrt(p) in Carpentier et al on page 1820
        psi <- sqrt(p) / m
        Qhat <- sum(yi^2) - sigma^2 * sum(diag(Ginv))
    }
    test.stat <- sqrt(max(Qhat,0))
    lambda <- sigma * sqrt(psi)
    tau <- sqrt(Cs / test.size) *lambda / 2 # testing boundary defined on top of p. 1821
    test.result <- ifelse(test.stat > tau, 1, 0)
    return(list(stat=test.stat, tau=tau, result=test.result))
}



##### debugging #####
if (sys.nframe()==0L){
#    for (i in 1:10){
    n1 <- 1200; n2 <- 1200; p <- 1000; k <- 100;
    X1 <- random.GaussianMatrix(n1, p)
    X2 <- random.GaussianMatrix(n2, p)
    beta1 <- rnorm(p)
    rho <- 0
    theta <- vector.normalise(c(rnorm(k), rep(0, p-k))) * rho
    beta2 <- beta1 + theta
    y1 <- X1 %*% beta1 + rnorm(n1)
    y2 <- X2 %*% beta2 + rnorm(n2)

#    us <- complementarySketching(X1,X2,y1,y2)$result
#    us_dense <- complementarySketching(X1,X2,y1,y2,sparse=F)$result
#    println(us,' ',us_dense)
#    }
#    print(us <- complementarySketching(X1,X2,y1,y2,sparse=F)$result)
#    LRT <- ifelse(LRT(X1,X2,y1,y2)$pval < 0.05, 1, 0)
#    ZhuBradic(X1,X2,y1,y2,Pi=diag(p))$pval
    ZB <- ifelse(ZhuBradic(X1,X2,y1,y2,Pi=diag(p))$pval < 0.05, 1, 0)
    Chab <- Chabonnier(X1, X2, y1, y2)$result
    Carp <- Carpentier(X1, X2, y1, y2, k, 1)
    println(ZB, ' ', Chab)
}
