Ftest <- function(lm.r, lm.ur) {
    ## conduct standard F-test of linear restrictions
    ## Arguments:
    ##  lm.r - linear model with restrictions
    ##  lm.ur - linear model without restrictions
    SSR.r <- sum(lm.r$res^2)
    SSR.ur <- sum(lm.ur$res^2)
    F <- (SSR.r-SSR.ur)/SSR.ur*(lm.ur$df/(lm.r$df-lm.ur$df))
    pval <- pf(F, lm.r$df-lm.ur$df, lm.ur$df, lower.tail=FALSE)
    return(list(F = F, pval = pval))
}

simVAR1 <- function(T, mu, Phi, Sigma) {
    T0 <- 500 # warm-up sample (alternatively draw from uncond. distr.)
    N <- length(mu)
    Xsim <- matrix(0, T0+T, N)
    for (t in 2:(T0+T)) {
        Xsim[t,] <- mu + Phi %*% Xsim[t-1,] + Sigma %*% rnorm(N)
    }
    return(Xsim[(T0+1):(T0+T),])
}

getReturns <- function(Y, n, h, mats) {
    T <- nrow(Y)
    return(-(n-h)*Y[(1+h):T, mats==n-h] + n*Y[1:(T-h), mats==n] - h*Y[1:(T-h), mats==h])
}

createTbl <- function(rows, cols) {
    tbl <- matrix(NA, length(rows), length(cols))
    rownames(tbl) <- rows
    colnames(tbl) <- cols
    return(tbl)
}

CrossProduct3D <- function(x, y, i=1:3) {
    ## Project inputs into 3D, since the cross product only makes sense in 3D.
    To3D <- function(x) head(c(x, rep(0, 3)), 3)
    x <- To3D(x)
    y <- To3D(y)

    ## Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
    ## so on).  Index3D() lets us do that using R's convention of 1-based (rather
    ## than 0-based) arrays.
    Index3D <- function(i) (i - 1) %% 3 + 1

    ## The i'th component of the cross product is:
    ## (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
    ## as long as we treat the indices cyclically.
    return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
                x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

vcovNW <- function(lm., k = 18)
    NeweyWest(lm., k, prewhite=FALSE)

getTest <- function(lm.ur, lm.r, vcov.fn) {
    require(lmtest)
    require(sandwich)
    if (missing(vcov.fn))
        vcov.fn <- vcovNW
    waldtest(lm.ur, lm.r, vcov=vcov.fn, test="Chisq")
}

getPval <- function(lm2, lm1)
    return(getTest(lm2, lm1)$Pr[2])

getTestStat <- function(lm2, lm1)
    return(getTest(lm2, lm1)$Chisq[2])

## PCA
makePC <- function(Y) {
    V <- cov(Y)
    print(cumsum(eigen(V)$values/sum(eigen(V)$values)))
    (w <- eigen(V)$vectors[,1])
    w <- w/sum(w)
    return(as.matrix(Y) %*% w)
}

checkKKT <- function(theta, obj, ...) {
    require(numDeriv)
    y <- obj(theta, ...)
    kkttol <- 10*.Machine$double.eps^(1/4)
    kkt2tol <- 100* (.Machine$double.eps^(1/4))
    ngatend <- grad(obj, theta, method="Richardson", side=NULL, method.args=list(), ...)
    cat("Gradient:")
    print(ngatend)
    kkt1 <- max(abs(ngatend)) <= kkttol*(1.0+abs(y))
    cat("kkt1 = ", kkt1, "\n")
    nhatend <- hessian(obj, theta,  method="Richardson", method.args=list(), ...)
    hev <- eigen(nhatend)$values
    cat("Eigenvalues:", hev, "\n")
    negeig <- (hev <= -kkttol*(1+abs(y)))
    cat("negeig = ", negeig, "\n")
    evratio <- tail(hev, 1)/hev[1]
    cat("evratio =", evratio, "\n")
    cat("evratio requirement >", kkt2tol,"\n")
    kkt2 <- (evratio > kkt2tol) && (!negeig)
    cat("kkt2 =", kkt2, "\n")
}

getOptim <- function(theta, obj, ..., trace=0) {
    obj <- match.fun(obj)
    if (trace>0) {
        cat('Starting optimization...\n')
        cat('Function value at starting point = ',
            sprintf("%10.4f", obj(theta, ...)), '\n')
    }
    i <- 1; improvement <- Inf; prev.llk <- 0
    myparscale <- 10^round(log10(abs(theta)))
    while (improvement>.1) {
        res <- optim(theta, obj, gr=NULL, ..., control=list(trace=trace, parscale=myparscale) )
        improvement <- abs(res$value-prev.llk)
        prev.llk <- res$value
        theta <- res$par
        if (trace > 0)
            cat('iteration ', i,', likelihood = ', sprintf("%10.4f", res$value),'\n')
        i <- i + 1
    }
    if (trace > 0)
        cat('improvement = ', improvement, ' -- proceed to final step\n')
    res <- optim(theta, obj, gr=NULL, ..., control=list(trace=trace, maxit=50000, parscale=myparscale) )
    if (trace > 0) {
        cat('final Nelder-Mead step, likelihood = ',
            sprintf("%10.4f", res$value),'\n')
        cat('Convergence:', res$convergence, "\n")
        print(res$message)
    }
    return(res$par)
}

plot.recessions <- function(yrange) {
    rec.90.from <- 1990+7/12
    rec.90.to <- 1991+3/12
    rec.01.from <- 2001+3/12
    rec.01.to <- 2001+11/12
    rec.07.from <- 2007+12/12
    rec.07.to <- 2009+6/12
    polygon(x=c( rec.90.from, rec.90.from, rec.90.to, rec.90.to),
            y=c(yrange, rev(yrange)),
            density=NA, col="gray", border=NA)
    polygon(x=c( rec.01.from, rec.01.from, rec.01.to, rec.01.to),
            y=c(yrange, rev(yrange)),
            density=NA, col="gray", border=NA)
    polygon(x=c( rec.07.from, rec.07.from, rec.07.to, rec.07.to),
            y=c(yrange, rev(yrange)),
            density=NA, col="gray", border=NA)
}

matrix.power <- function(x, y) {
    ## calculate matrix power, allowing for non-integer exponents
    ## x^y where x is a matrix, y a scalar
    if (length(x)>1) {
        eig <- eigen(x)
        x.y <- eig$vectors%*%diag(eig$values^y)%*%solve(eig$vectors)
        dimnames(x.y) <- dimnames(x)
    } else {
        x.y <- x^y
    }
    return(Re(x.y))
}
