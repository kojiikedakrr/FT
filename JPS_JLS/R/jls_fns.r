lam.LB <- -0.5
lam.UB <- 0

## y is constrained to be between yLB and yUB, while thetea is unrestricted
## y2theta <- function(y, a = lam.LB, b = lam.UB)
##     log(y - a) - log(b - y)
## theta2y <- function(theta, a = lam.LB, b = lam.UB)
##     a + (b-a)*exp(theta)/(1+exp(theta))
## alternatives: use rescaling of y = tanh(theta) (see Nash p. 134) or y = theta/sqrt(1+theta^2)
stretch <- 1
y2theta <- function(y, a = lam.LB, b = lam.UB) {
    theta <- atanh(2*(y - a)/(b-a)-1)/stretch
    if (any(is.infinite(theta))) {
        ind <- which(is.infinite(theta))
        theta[ind] <- sign(theta[ind])*10  # this will make it close enough to boundary
    }
    theta
}
theta2y <- function(theta, a = lam.LB, b = lam.UB)
    a + .5*(b-a)*(1+tanh(stretch*theta))

scale.gam <- 1
scale.Sigma <- 1

theta2pars.jls <- function(theta) {
    cM <- ncol(M.o)
    cN <- cM+cL
    pars <- list(lamQ = sort(theta2y(theta[1:cN]), decr=TRUE))
    pars$gam0 <- matrix(theta[(cN+1):(cN+cM)], cM, 1)/scale.gam
    pars$gam1 <- matrix(theta[(cN+1+cM):(cN+cM+cM*cN)], cM, cN)/scale.gam
    pars$Sigma <- matrix(0, cN, cN);
    pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)] <- tail(theta, cN*(cN+1)/2)/scale.Sigma
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)
    pars
}

pars2theta.jls <- function(pars) {
    cM <- ncol(M.o)
    cN <- cM+cL
    if (length(pars$lamQ)!=cN) stop("lamQ has wrong length")
    Sigma.vec <- pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)]
    c(y2theta(pars$lamQ), scale.gam*pars$gam0, scale.gam*as.numeric(pars$gam1), Sigma.vec*scale.Sigma)
}

pars2theta.jls.onlyLam <- function(pars) {
    cM <- ncol(M.o)
    cN <- cM+cL
    cL <- cN-cM
    if (length(pars$lamQ)!=cN) stop("lamQ has wrong length")
    dlamQ <- c(pars$lamQ[1], diff(pars$lamQ))
    y2theta(dlamQ)
    ## y2theta(pars$lamQ)
}

theta2pars.jls.onlyLam <- function(theta) {
    ## global: Sigma.fix, gam0.fix, gam1.fix
    cM <- ncol(M.o)
    cN <- cM+cL
    pars <- list(gam0 = gam0.fix,
                 gam1 = gam1.fix,
                 Sigma = Sigma.fix)
    ## pars$lamQ <- sort(theta2y(theta), decr=TRUE)
    pars$dlamQ <- theta2y(theta)
    pars$lamQ <- cumsum(pars$dlamQ)
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)
    pars
}

theta2pars.jls.onlyGam <- function(theta) {
    ## global: Sigma.fix, lamQ.fix
    cM <- ncol(M.o)
    cN <- cM+cL
    pars <- list(gam0 = matrix(theta[1:cM], cM, 1)/scale.gam,
                 gam1 = matrix(theta[(1+cM):(cM+cM*cN)], cM, cN)/scale.gam,
                 Sigma = Sigma.fix,
                 lamQ = lamQ.fix)
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)
    pars
}

pars2theta.jls.onlyGam <- function(pars)
    c(pars$gam0, as.numeric(pars$gam1))*scale.gam

theta2pars.jls.onlySig <- function(theta) {
    cM <- ncol(M.o)
    cN <- cM+cL
    pars <- list(lamQ = lamQ.fix,
                 gam0 = gam0.fix,
                 gam1 = gam1.fix)
    pars$Sigma <- matrix(0, cN, cN);
    pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)] <- theta/scale.Sigma
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)
    pars
}

pars2theta.jls.onlySig <- function(pars)
    as.numeric(pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)])*scale.Sigma

checkPars <- function(pars) {
    ## check parameter restrictions common to all models
    valid <- TRUE
    if (min(diff(sort(pars$lamQ)))<1e-6) valid <- FALSE
    if (any(pars$lamQ < lam.LB)) valid <- FALSE
    if (any(pars$lamQ > lam.UB)) valid <- FALSE
    return(valid)
}

obj.jls <- function(theta) {
    ## objective function for ML estimation of JLS model -- observed risk factors
    ## Arguments:
    ##   theta - vector with model parameters
    ## Value: negative of log-likelihood
    ## Globals: W, Y, M.o, cL, mats
    cN <- ncol(M.o) + cL
    pars <- theta2pars.jls(theta)
    if (checkPars(pars)) {
        res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
        obj <- sum(res.llk$llk)
    } else {
        return(1e6)
    }
    ## } else {
    ##     obj <- 1000
    ## }
}

obj.jls.fixedSigma <- function(theta) {
    cN <- ncol(M.o) + cL
    pars <- theta2pars.jls.fixedSigma(theta)
    res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
    sum(res.llk$llk)
}

obj.jls.onlyLam <- function(theta) {
    cN <- ncol(M.o) + cL
    pars <- theta2pars.jls.onlyLam(theta)
    if (checkPars(pars)) {
        res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
        sum(res.llk$llk)
    } else {
        return(1e6)
    }
}

obj.jls.onlySig <- function(theta) {
    cN <- ncol(M.o) + cL
    pars <- theta2pars.jls.onlySig(theta)
    res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
    sum(res.llk$llk)
}

obj.jls.onlyGam <- function(theta) {
    cN <- ncol(M.o) + cL
    pars <- theta2pars.jls.onlyGam(theta)
    res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
    sum(res.llk$llk)
}

estJLS <- function(Y, W, M.o, mats, cL, pars.start) {
    getStartingValuesJLS <- function(Sigma) {
        ## get good starting values for JLS model estimation using random seeds
        ## Arguments: Sigma from VAR
        ## Value: list with starting values
        ## Globals: W, Y, M.o, n.per
        cM <- ncol(M.o)
        cN <- ncol(Sigma)
        J <- ncol(Y)
        WN <- matrix(W[1:cN,], cN, J)

        ## starting values for gam0 and gam1: regress macro on cN yield factors
        gam0 <- matrix(NA, cM, 1)
        gam1 <- matrix(NA, cM, cN)
        PN.o <- Y %*% t(WN)
        xdat <- PN.o
        for (i in 1:cM) {
            ydat <- M.o[,i]
            res <- lm(ydat ~ xdat)
            gam0[i] <- res$coef[1]
            gam1[i,] <- res$coef[2:(cN+1)]
        }

        ## random starting values for lamQ and kinfQ
        n.seeds <- 500
        best.llk <- Inf
        pars <- list(gam0 = gam0, gam1 = gam1, Sigma = Sigma)
        for (i in 1:n.seeds) {
            pars$lamQ <- -sort(abs(.1*rnorm(cN)))
            theta <- pars2theta.jls(pars)
            llk <- obj.jls(theta)
            if (llk < best.llk) {
                cat('Improved seed llk to ', llk, '\n')
                best.llk <- llk
                best.pars <- pars
            }
        }
        return(best.pars)
    }
    theta2pars.jls <- function(theta) {
        cM <- ncol(M.o)
        cN <- cM+cL
        pars <- list(dlamQ = theta[1:cN])
        pars$lamQ <- cumsum(pars$dlamQ)
        pars$gam0 <- matrix(theta[(cN+1):(cN+cM)], cM, 1)
        pars$gam1 <- matrix(theta[(cN+cM+1):(cN+cM+cM*cN)], cM, cN)
        pars$Sigma <- matrix(0, cN, cN);
        pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)] <- tail(theta, cN*(cN+1)/2)
        pars$Omega <- pars$Sigma %*% t(pars$Sigma)
        return(pars)
    }

    pars2theta.jls <- function(pars) {
        cM <- nrow(pars$gam0)
        cN <- ncol(pars$gam1)
        cL <- cN-cM
        dlamQ <- c(pars$lamQ[1],diff(pars$lamQ))
        if (length(pars$lamQ)!=cN) stop("lamQ has wrong length")
        Sigma.vec <- pars$Sigma[lower.tri(pars$Sigma,diag=TRUE)]
        theta <- c(dlamQ, pars$gam0, as.numeric(pars$gam1), Sigma.vec)
        return(theta)
    }
    checkPars <- function(pars) {
        ## check parameter restrictions common to all models
        valid <- TRUE
        ## diagonal elements of Sigma positive and bounded away from zero
        if (any(diag(pars$Sigma)<1e-7)) valid <- FALSE
        if (any(diag(pars$Sigma)>1)) valid <- FALSE
        ## Q-eigenvalues not explosive
        if (any(pars$lamQ>0)) valid <- FALSE
        if (any(pars$lamQ< -.4)) valid <- FALSE
        ## Q-eigenvalues sorted
        if (any(pars$dlamQ>0)) valid <- FALSE
        return(valid)
    }

    obj.jls <- function(theta) {
        ## objective function for ML estimation of JLS model -- observed risk factors
        ## Arguments:
        ##   theta - vector with model parameters
        ## Value: negative of log-likelihood
        ## Globals: W, Y, M.o, cL, mats
        cN <- ncol(M.o) + cL
        pars <- theta2pars.jls(theta)

        ## check restrictions on parameter space
        valid <- checkPars(pars)
        if (valid) {
            res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
            obj <- sum(res.llk$llk)
        } else {
            obj <- 1e6
        }
    }

    cM <- ncol(M.o)
    cN <- cM + cL
    J <- ncol(Y)

    WL <- matrix(W[1:cL,], cL, J)

    PL.o <- Y %*% t(WL)
    Z <- cbind(M.o, PL.o)

    ## starting values
    if (missing(pars.start)) {
        ## VAR for Z
        ar1 <- ar.ols(Z, aic=FALSE, order.max=1, intercept=TRUE, demean=FALSE)
        K1P <- ar1$ar[,,] - diag(cN)
        K0P <- ar1$x.intercept
        Omega <- ar1$var.pred
        Sigma <- t(chol(Omega))
        Sigma.ols <- Sigma
        ## various tries for other params
        pars.start <- getStartingValuesJLS(Sigma)
    }

    ## cat("starting values for gam0/gam1:\n")
    ## print(cbind(pars.start$gam0, pars.start$gam1))

    theta.start <- pars2theta.jls(pars.start)
    stopifnot(isTRUE(all.equal(theta.start, pars2theta.jls(theta2pars.jls(theta.start)), check.attributes=FALSE)))

    ## optimization over (lamQ, gam0, gam1, Sigma)
    theta <- getOptim(theta.start, obj.jls)
    pars <- theta2pars.jls(theta)

    ## cat("optimal values for gam0/gam1:\n")
    ## print(cbind(pars$gam0, pars$gam1))

    pars$K0P <- K0P
    pars$K1P <- K1P
    pars$Omega <- pars$Sigma %*% t(pars$Sigma)

    ## value of likelihood
    res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
    pars$sigma.e <- res.llk$sigma.e
    pars$kinfQ <- res.llk$kinfQ
    pars$llkP <- -sum(res.llk$llkP)
    pars$llkQ <- -sum(res.llk$llkQ)
    pars$llk <- pars$llkP + pars$llkQ
    pars$mu <- pars$K0P
    pars$Phi <- pars$K1P + diag(cN)
    pars$rho0 <- as.numeric(res.llk$rho0)
    pars$rho1 <- as.numeric(res.llk$rho1)
    pars$muQ <- as.numeric(res.llk$K0Q)
    pars$PhiQ <- res.llk$K1Q + diag(cN)
    pars$cL <- cL
    pars$A <- res.llk$A
    pars$B <- res.llk$B
    pars$Gam1 <- rbind(pars$gam1, cbind(diag(cL), matrix(0, cL, cM)))
    pars$Gam0 <- rbind(pars$gam0, matrix(0, cL,1))
    pars$Z <- Z
    cat("LLK          =", -sum(res.llk$llk), " \n")
    cat("check        =", -obj.jls(theta),"\n")

    pars2 <- pars
    pars2$Sigma <- Sigma.ols
    theta2 <- pars2theta.jls(pars2)
    cat("w/ Sigma.ols =", -obj.jls(theta2), "\n")

    return(pars)
}

estJLSnew <- function(Y, W, M.o, mats, cL, pars.start) {
    cM <- ncol(M.o)
    cN <- cM + cL
    J <- ncol(Y)
    WL <- matrix(W[1:cL,], cL, J)
    WN <- matrix(W[1:cN,], cN, J)
    PL.o <- Y %*% t(WL)
    PN.o <- Y %*% t(WN)
    Z <- cbind(M.o, PL.o)

    ## starting values
    ## VAR for Z
    ar1 <- ar.ols(Z, aic=FALSE, order.max=1, intercept=TRUE, demean=FALSE)
    K1P <- ar1$ar[,,] - diag(cN)
    K0P <- ar1$x.intercept
    Omega <- ar1$var.pred
    Sigma.ols <- t(chol(Omega))
    ## gam0 and gam1
    xdat <- cbind(1, PN.o)
    ydat <- M.o
    A <- t(ydat)%*%xdat%*%solve(t(xdat) %*% xdat)
    gam0.ols <- matrix(A[,1], cN-cL, 1)
    gam1.ols <- A[,-1]

    gam0.fix <<- gam0.ols
    gam1.fix <<- gam1.ols
    Sigma.fix <<- Sigma.ols

    ## Sobol quasi-random sequence
    N0 <- 100000
    N <- 20

    cat("## find best starting values for lamQ\n")
    lamQ <- lam.LB + (lam.UB-lam.LB)*sobol(N0, cN)
    lamQ <- lamQ[!apply(lamQ, 1, function(x) length(unique(x))<5),] ## no repeated eigenvalues
    lamQ <- t(apply(lamQ, 1, sort, decr=TRUE))
    obj1 <- function(lamQ) {
        res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=lamQ, gam0=gam0.fix, gam1=gam1.fix, Sigma=Sigma.fix, mats=mats, dt=1)
        sum(res.llk$llk)
    }
    llks1 <- apply(lamQ, 1, obj1)
    I1 <- head(order(llks1), N)
    ## lamQ.fix <- lamQ[I1[1],]
    ## I1 <- sample(I1)  # reorder - don't pair best lamQ with best gam starting values
    print(apply(lamQ[I1,], 1, obj1))

    ## optimization

    tbl <- matrix(NA, N, 7)
    colnames(tbl) <- c("llk start", "llk opt", paste("lamQ", 1:5))
    best.llk <- Inf
    for (i in 1:N) {
        cat("## starting values", i, "\n")
        pars <- list(lamQ = lamQ[I1[i],],
                     gam0 = gam0.ols,
                     gam1 = gam1.ols,
                     Sigma = Sigma.ols)
        theta.start <- pars2theta.jls(pars)
        last.llk <- obj.jls(theta.start)
        tbl[i, 1] <- last.llk
        cat("LLK at starting values:", last.llk, "\n")
        gam0.fix <<- pars$gam0
        gam1.fix <<- pars$gam1
        Sigma.fix <<- pars$Sigma
        cat("iterative optimization over parameter blocks\n")
        llk.decr <- -1
        theta1 <- pars2theta.jls.onlyLam(pars)
        theta2 <- pars2theta.jls.onlyGam(pars)
        theta3 <- pars2theta.jls.onlySig(pars)
        while (llk.decr < -.001) {
            ## optimize lamQ for given gam0, gam1, Sigma
            llk <- obj.jls.onlyLam(theta1)
            cat("before lamQ opt: LLK =", llk, "\n")
            if ( llk > last.llk) {
                warning(paste("neg-LLK jumped up: from", last.llk, "to", llk))
            }
            theta1 <- getOptim(theta1, obj.jls.onlyLam)
            llk <- obj.jls.onlyLam(theta1)
            pars <- theta2pars.jls.onlyLam(theta1)
            lamQ.fix <<- pars$lamQ
            cat("after lamQ opt: LLK =", llk, "\n")
            ## optimize gam0,gam1 for given lamQ, Sigma
            theta2 <- getOptim(theta2, obj.jls.onlyGam)
            llk <- obj.jls.onlyGam(theta2)
            pars <- theta2pars.jls.onlyGam(theta2)
            gam0.fix <<- pars$gam0
            gam1.fix <<- pars$gam1
            cat("after gam0/gam1 opt: LLK =", llk, "\n")
            ## optimize Sigma for given lamQ, gam0, gam1
            theta3 <- getOptim(theta3, obj.jls.onlySig)
            llk <- obj.jls.onlySig(theta3)
            pars <- theta2pars.jls.onlySig(theta3)
            Sigma.fix <<- pars$Sigma
            cat("after Sigma opt: LLK =", llk, "\n")
            llk.decr <- llk - last.llk
            cat("decrease in this iteration:", llk.decr, "\n")
            last.llk <- llk
        }
        theta <- pars2theta.jls(pars)
        cat("after again jointly optimizing over all parameters\n")
        theta <- getOptim(theta, obj.jls)
        llk <- obj.jls(theta)
        print(llk)
        tbl[i, 2] <- llk
        pars <- theta2pars.jls(theta)
        tbl[i, 3:7] <- pars$lamQ
        if (llk < best.llk) {
            best.llk <- llk
            best.pars <- pars
        }
    }
    print(tbl)
    pars <- best.pars

    ## value of likelihood
    res.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
    cat("Check:", sum(res.llk$llk), "\n")
    pars$K0P <- K0P
    pars$K1P <- K1P
    pars$sigma.e <- res.llk$sigma.e
    pars$kinfQ <- res.llk$kinfQ
    pars$llkP <- -sum(res.llk$llkP)
    pars$llkQ <- -sum(res.llk$llkQ)
    pars$llk <- pars$llkP + pars$llkQ
    pars$mu <- pars$K0P
    pars$Phi <- pars$K1P + diag(cN)
    pars$rho0 <- as.numeric(res.llk$rho0)
    pars$rho1 <- as.numeric(res.llk$rho1)
    pars$muQ <- as.numeric(res.llk$K0Q)
    pars$PhiQ <- res.llk$K1Q + diag(cN)
    pars$cL <- cL
    pars$A <- res.llk$A
    pars$B <- res.llk$B
    pars$Gam1 <- rbind(pars$gam1, cbind(diag(cL), matrix(0, cL, cM)))
    pars$Gam0 <- rbind(pars$gam0, matrix(0, cL,1))
    pars$Z <- Z
    pars
}


loadJLSmodel <- function(name, filename) {
    cat("*** loading JLS model", name, "\n")
    cat("*** results file:", filename, "\n")
    load(filename)  # pars, est.llk, Y, W, M.o, cL, mats, dates, n.per
    ## combine saved objects to model
    m <- append(pars, est.llk)
    m$name <- name
    m$Y <- Y
    m$W <- W
    m$M.o <- M.o
    m$cL <- cL
    m$mats <- mats
    m$dates <- dates
    ## additional objects
    m$cM <- ncol(M.o)
    cN <- m$cM+cL
    m$cN <- cN
    m$N <- cN
    J <- length(m$mats)
    m$WL <- matrix(W[1:cL,], cL, J)
    m$WN <- matrix(W[1:cN,], cN, J)
    ## short rate
    m$r <- m$rho0 + m$Z %*% m$rho1
    ## fitted yields
    T <- nrow(Y)
    m$Yhat <- rep(1,T) %*% m$A + m$Z %*% m$B
    ## term premium
    loads.rn <- gaussian.loadings(m$mats, m$mu, m$Phi - diag(cN),
                                  m$Omega, m$rho0, m$rho1)
    m$Yrn <- rep(1, T) %*% loads.rn$A + m$Z %*% loads.rn$B
    m$Ytp <- m$Yhat - m$Yrn
    ## rotated risk factors
    m$PN <- (m$Z - matrix(rep(m$Gam0, each=T), T, cN)) %*% t(solve(m$Gam1))
    ## alternativ: PN = m$Yhat %*% t(m$WN)
    ## population moments
    m$Cov.Z <- matrix(solve(diag(cN^2) - kronecker(m$Phi, m$Phi)) %*% as.numeric(m$Omega), cN, cN)
    m$E.Z <- as.vector(solve(diag(cN) - m$Phi) %*% m$mu)
    m <- append(m, getForwardRates(m))
    ##    m <- append(m, getForwardRates2(m))
    return(m)
}

getSummary <- function(model) {
    cN <- length(model$rho1)
    EQ.r <- 100*n.per*(model$rho0 + crossprod(model$rho1, solve(diag(cN) - model$PhiQ) %*% model$muQ))
    EP.r <- 100*n.per*(model$rho0 + crossprod(model$rho1, solve(diag(cN) - model$Phi) %*% model$mu))
    r.bar <- 100*n.per*(model$rho0 + crossprod(model$rho1, colMeans(model$Z)))

    ## fit
    T <- nrow(Y)
    Y.hat <- rep(1,T) %*% model$A + model$Z %*% model$B
    rmse <- 10000*n.per*sqrt( mean((Y-Y.hat)^2) )

    x <- c(model$llk, model$kinfQ, EQ.r, EP.r, r.bar,
           1+model$lamQ[1:3], diag(model$Phi)[1:3],
           abs(eigen(model$Phi)$values[1:3]),
           irf.var1(model$Phi, max.lag=120, g=3, h=3)[120],
           model$sigma.e*120000, rmse)

    names(x) <- c("LLK", "kinfQ", "EQ.r", "EP.r", "r.bar",
                  "lamQ.1", "lamQ.2", "lamQ.3",
                  "Phi.11", "Phi.22", "Phi.33",
                  "ev-P (1)", "ev-P (2)", "ev-P (3)", "IRF-P(10y)",
                  "sigma.e", "RMSE")
    return(x)
}

checkSpanning <- function(m) {
    ## (Q1) compare model-implied risk factors to PCs of yields
    T <- nrow(Y)
    J <- ncol(Y)
    cM <- ncol(M.o)
    cN <- cM+cL
    WN <- matrix(W[1:cN,], cN, J)

    PN <- (m$Z - matrix(rep(m$Gam0, each=T), T, cN)) %*% t(solve(m$Gam1))
    T <- nrow(Y)
    Yhat <- rep(1,T) %*% m$A + m$Z %*% m$B
    PN.hat <- Yhat %*% t(WN)
    PN.o <- Y %*% t(WN)

    tbl <- cbind(PN[1,], PN.hat[1,], PN.o[1,])
    colnames(tbl) <- c("P^N", "PCs of Yhat", "PCs of Y")
    print(tbl)
    ## first cL elements are the same:

    ## (Q3) Check spanning condition
    cat("Spanning condition:\n")
    print(cbind(m$gam0 + m$gam1 %*% PN[1,], M.o[1,]))

    ## check parameters
    print(m$gam0 + m$gam1 %*% WN %*% t(m$A))
    print(m$gam1 %*% WN %*% t(m$B))
}

jls.llk.kinfQ <- function (yields.o, M.o, W, cL, kinfQ=NA, lamQ, gam0, gam1, K0P=NA, K1P=NA, Sigma, mats, dt, sigma.e=NA) {
    ## Setup
    T <- nrow(yields.o)-1
    J <- ncol(yields.o)
    cM <- ncol(M.o)
    cN <- cL + cM
    WN <- matrix(W[1:cN,], cN, J)
    WL <- matrix(W[1:cL,], cL, J)
    PL.o <- ts( yields.o %*% t(WL))
    Z <- cbind(M.o, PL.o)
    Omega <- Sigma %*% t(Sigma)

    if (is.na(kinfQ)) {
        ## concentrate out kinfQ
        ## AZ = alpha0_Z*kinf + alpha1_Z
        rho0 <- 0
        ## AZ0, AX0 will be the loadings with rho0_Z = 0, which won't be correct
        loads <- jls.loadings.rho0(W, rho0, lamQ, gam0, gam1, Omega, mats, dt)
        B <- loads$B
        alpha0.Z <- loads$alpha0.Z
        alpha1.Z <- loads$alpha1.Z

        ## back out kinfQ that fits average yields
        require(MASS)
        V <- t(Null(t(WL)))
        kinfQ <- as.numeric(t(colMeans(yields.o[2:(T+1),]) - t(alpha1.Z) - t(B)%*%colMeans(Z[2:(T+1),]))%*%(t(V)%*%V%*%t(alpha0.Z)) / (alpha0.Z%*%t(V)%*%V%*%t(alpha0.Z)))

        ## get correct loadings
        A <- alpha0.Z*kinfQ + alpha1.Z;

        ## get these to return to caller
        AX <- loads$alpha0.X*kinfQ + loads$alpha1.X;
        BX <- loads$BX
        rho1 <- loads$rho1
        K1Q <- loads$K1Q
        U0 <- loads$Gam0 + loads$Gam1 %*% WN %*% t(AX)
        rho0 <- -crossprod(rho1, U0)
        K0Q.X <- matrix(0, cN, 1)
        K0Q.X[loads$m1] <- kinfQ
        K0Q <- loads$U1inv %*% K0Q.X - K1Q %*% U0
    } else {
        loads <- jls.loadings.kinfQ(W, kinfQ, lamQ, gam0, gam1, Omega, mats, dt)
        B <- loads$B; A <- loads$A
        AX <- loads$AX; BX <- loads$BX
        K0Q <- loads$K0Q; K1Q <- loads$K1Q
        rho0 <- loads$rho0; rho1 <- loads$rho1
    }

    yields.m <- rep(1,T+1)%*%A + Z %*% B # (T+1)*J, model-implied yields
    yield.errors <- yields.o[2:(T+1),] - yields.m[2:(T+1),]; # T*J
    squared.errors <- yield.errors^2; # T*J, but cL-dimensional projection onto W is always 0, so effectively (J-cL) dimensional

    ## cat("RMSE =", sqrt(mean(squared.errors))*120000, "\n")

    ## Compute optimal sigma.e if it is not supplied
    if (is.na(sigma.e))
        sigma.e <- sqrt( sum(squared.errors)/(T*(J-cL)) )

    llkQ <- .5*rowSums(squared.errors)/sigma.e^2 + (J-cL)*.5*log(2*pi) + .5*(J-cL)*log(sigma.e^2) # 1*T

    if (missing(K0P)|missing(K1P)) {
        ## Run OLS to obtain maximum likelihood estimates of K0P, K1P
        var1 <- ar.ols(Z, order.max=1, aic=FALSE, demean=FALSE, intercept=TRUE)
        K1P <- var1$ar[,,] - diag(cN)
        K0P <- var1$x.intercept
    }

    innovations = t(Z[2:(T+1),]) - (K0P%*%matrix(1,1,T) + (K1P+diag(cN))%*%t(Z[1:T,])) # N*T

    llkP = .5*cN*log(2*pi) + .5*log(det(Omega)) + .5*colSums(innovations*solve(Omega, innovations)) # 1*T

########################################################################

    jsz.llk <- list(llk=t(llkQ + llkP), A=A, B=B, AX=AX, BX=BX, K0P=K0P, K1P=K1P, sigma.e=sigma.e, llkQ=llkQ, llkP=llkP, kinfQ = kinfQ, K0Q=K0Q, K1Q=K1Q, rho0=rho0, rho1=rho1)

}

jls.loadings.kinfQ <- function(W, kinfQ, lamQ, gam0, gam1, Omega, mats, dt) {
    ## Inputs:
    ##   W          : J*J,      vector of portfolio weights
    ##   kinfQ      : scalar,   long run mean
    ##   lamQ       : N*1,      Q-eigenvalues
    ##   gam0       : cM*1      spanning parameters -- intercept
    ##   gam1       : cM*cN     spanning parameters -- coefficients
    ##   Omega      : cN*cN,    cov of innovations to Z
    ##   mats       : 1*J,      maturities in years
    ##   dt         : scalar,   length of period in years
    ##
    ## Returns:
    ##   A    : 1*J
    ##   B    : N*J
    ##   K0Q  : N*1
    ##   K1Q  : N*N
    ##   rho0 : scalar
    ##   rho1 : N*1

    J <- length(mats)
    cN <- length(lamQ)
    cM <- 2 ## default
    cL <- cN - cM
    WN <- matrix(W[1:cN,], cN, J)
    mats.periods <- round(mats/dt)
    M <- max(mats.periods)
    rho0.X <- 0;   rho1.X <- rep(1, cN)

    ## 1. based on primitive parameters, find (rho0, rho1, K0Q, K1Q)
    Gam1 <- rbind(gam1, cbind(diag(cL), matrix(0, cL, cM)))
    Gam0 <- rbind(gam0, matrix(0, cL,1))

    ## 1.1. find loadings of yields on Xt (Jordan-normalized factors)
    K1Q.X <- diag(lamQ)
    ##  adjK1QX <- jszAdjustK1QX(K1Q.X)
    ##  K1Q.X <- adjK1QX$K1Q.X
    ##  m1 <- adjK1QX$m1
    m1 <- 1
    K0Q.X <- matrix(0, cN, 1);
    K0Q.X[m1] <- kinfQ

    ## 1.1.1. first compute the loadings ignoring convexity term
    loads.X.prelim <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, matrix(0, cN, cN), rho0.X*dt, rho1.X*dt, dt)
    BX <- loads.X.prelim$B  ## cN * J
    U1inv <- Gam1 %*% WN %*% t(BX)
    U1 <- solve(U1inv)
    ## 1.1.2. calc. Omega.X and calculate correct loadings
    Omega.X <- U1 %*% Omega %*% t(U1)
    loads.X <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, Omega.X, rho0.X*dt, rho1.X*dt, dt)
    AX <- loads.X$A ## 1 * J

    ## 1.2. calculate remaining parameters
    U0 <- Gam0 + Gam1 %*% WN %*% t(AX)
    rho1 <- t(U1) %*% rep(1, cN)
    rho0 <- - crossprod(rho1, U0)
    K1Q <- U1inv %*% K1Q.X %*% U1
    K0Q <- U1inv %*% K0Q.X - K1Q %*% U0

    ## 2. compute affine loadings for Z
    loads.Z <- gaussian.loadings(mats.periods, K0Q, K1Q, Omega, rho0*dt, rho1*dt, dt)

############################################################
    jls.loadings <- list(A=loads.Z$A, B=loads.Z$B, AX=AX, BX=BX, K0Q=K0Q, K1Q=K1Q, rho0=rho0, rho1=rho1, U0=U0, U1=U1, Gam0=Gam0, Gam1=Gam1, AX=AX, BX=BX)

}

jls.loadings.rho0 <- function(W, rho0, lamQ, gam0, gam1, Omega, mats, dt) {
    ## like jls.loadings.kinfQ but parameterized in terms of rho0 instead of kinfQ
    ## Inputs:
    ##   W          : J*J,      vector of portfolio weights
    ##   kinfQ      : scalar,   long run mean
    ##   lamQ       : N*1,      Q-eigenvalues
    ##   gam0       : cM*1      spanning parameters -- intercept
    ##   gam1       : cM*cN     spanning parameters -- coefficients
    ##   Omega      : cN*cN,    cov of innovations to Z
    ##   mats       : 1*J,      maturities in years
    ##   dt         : scalar,   length of period in years
    ##
    ## Returns:
    ##   A    : 1*J
    ##   B    : N*J
    ##   K0Q  : N*1
    ##   K1Q  : N*N
    ##   rho0 : scalar
    ##   rho1 : N*1

    J <- length(mats)
    cN <- length(lamQ)
    cM <- 2 ## default
    cL <- cN - cM
    WN <- matrix(W[1:cN,], cN, J)
    mats.periods <- round(mats/dt)
    M <- max(mats.periods)
    rho0.X <- 0;   rho1.X <- rep(1, cN)

    ## 1. based on primitive parameters, find (rho0, rho1, K0Q, K1Q)
    Gam1 <- rbind(gam1, cbind(diag(cL), matrix(0, cL, cM)))
    Gam0 <- rbind(gam0, matrix(0, cL,1))

    ## 1.1. find loadings of yields on Xt (Jordan-normalized factors)
    K1Q.X <- diag(lamQ)
    m1 <- 1
    K0Q.X <- matrix(0, cN, 1);
    K0Q.X[m1] <- 1        ## changed -- 1 instead of kinfQ

    ## 1.1.1. first compute the loadings ignoring convexity term
    loads.X.prelim <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, matrix(0, cN, cN), rho0.X*dt, rho1.X*dt, dt)
    BX <- loads.X.prelim$B  ## cN * J
    alpha0.X <- loads.X.prelim$A    ### added
    U1inv <- Gam1 %*% WN %*% t(BX)
    U1 <- solve(U1inv)

    ## 1.1.2. calc. Omega.X and calculate correct loadings
    Omega.X <- U1 %*% Omega %*% t(U1)
    loads.X <- gaussian.loadings(mats.periods, K0Q.X, K1Q.X, Omega.X, rho0.X*dt, rho1.X*dt, dt)
    AX1 <- loads.X$A ## 1 * J     ### changed -- AX1 instead of AX
    alpha1.X <- AX1 - alpha0.X    ### added

    ## Need to find what kinf should be to get the desired rho0:
### all the following added
    a0 <- -matrix(1,1,cN) %*% U1 %*% Gam1 %*% WN %*% t(alpha0.X)
    a1 <- -matrix(1,1,cN) %*% U1 %*% (Gam0 + Gam1 %*% WN %*% t(alpha1.X))
    kinfQ <- as.numeric((rho0 - a1)/a0);  ## a0*kinfQ + a1 = rho0
    K0Q.X[m1] <- kinfQ;
    AX <- alpha0.X*kinfQ + alpha1.X
    ## same as gaussian.loadings(mats.periods, K0Q.X, K1Q.X, Omega.X, rho0.X*dt, rho1.X*dt, dt)$A

    ## this seems to be wrong
    C <- diag(J) - t(BX) %*% U1 %*% Gam1 %*% WN
    alpha0.Z <- t(C %*% t(alpha0.X))
    alpha1.Z <- t(C %*% t(alpha1.X) - t(BX) %*% U1 %*% Gam0)

    ## 1.2. calculate remaining parameters
    U0 <- Gam0 + Gam1 %*% WN %*% t(AX)
    rho1 <- t(U1) %*% rep(1, cN)
    ## check
    ## -crossprod(rho1,U0) == rho0
    K1Q <- U1inv %*% K1Q.X %*% U1
    K0Q <- U1inv %*% K0Q.X - K1Q %*% U0

    ## 2. compute affine loadings for Z
    loads.Z <- gaussian.loadings(mats.periods, K0Q, K1Q, Omega, rho0*dt, rho1*dt, dt)

    ## identical:
    ## print(AZ <- AX - t(t(BX) %*% U1 %*% U0))
    ## print(alpha0.Z*kinfQ + alpha1.Z)
    ## print(loads.Z$A)
    ## print(jls.loadings.kinfQ(W, kinfQ, lamQ, gam0, gam1, Omega, mats, dt)$A)

    ## identical:
    ## print(BZ <- t(U1) %*% BX)
    ## print(loads.Z$B)
    ## print(all.equal(BZ, loads.Z$B))

############################################################
    jls.loadings <- list(A=loads.Z$A, B=loads.Z$B, AX=AX, BX=BX, K0Q=K0Q, K1Q=K1Q, rho0=rho0, rho1=rho1, U0=U0, U1=U1, U1inv=U1inv, Gam0=Gam0, Gam1=Gam1, AX=AX, BX=BX,  alpha0.X=alpha0.X, alpha1.X=alpha1.X, alpha0.Z=alpha0.Z, alpha1.Z=alpha1.Z, m1=m1)

}
