## functions to compare models and analyze their spanning implications

yieldFit <- function(models) {
    m <- models[[1]]
    J <- length(m$mats)
    tbl <- matrix(NA, length(models), J+2)
    colnames(tbl) <- c("sigma.e", "All", m$mats)
    rownames(tbl) <- names(models)
    for (i in 1:length(models)) {
        m <- models[[i]]
        sig.e <- 120000*m$sigma.e
        rmse <- 10000*12*sqrt( mean((m$Y-m$Yhat)^2) )
        rmse.mats <- 10000*12*sqrt( colMeans((m$Y-m$Yhat)^2) )
        tbl[i, ] <- c(sig.e, rmse, rmse.mats)
    }
    tbl
}

errorSerialCorr <- function(models) {
    m <- models[[1]]
    J <- length(m$mats)
    tbl <- matrix(NA, length(models), J)
    colnames(tbl) <- m$mats
    rownames(tbl) <- names(models)
    for (i in 1:length(models)) {
        m <- models[[i]]
        e <- m$Y-m$Yhat
        for (j in 1:length(m$mats))
            tbl[i, j] <- acf(e[,j], plot=FALSE)$acf[2]
    }
    tbl
}


calcSpanning <- function(m, yieldFactorsOnly = FALSE, sigmas=c(6, 0)/120000, h=12) {
    ## calculate spanning metrics analytically
    tblVar <- matrix(NA, length(sigmas), ncol(m$M.o))
    tblRisk <- matrix(NA, length(sigmas), 3)       ## R^2 without and with macro, rel. RMSE
    tblForecasts <- matrix(NA, length(sigmas), ncol(m$M.o))  ## rel. RMSE
    rownames(tblVar) <- paste(paste("sigma =", round(sigmas*120000, digi=2)))
    rownames(tblRisk) <- rownames(tblVar)
    rownames(tblVar) <- rownames(tblVar)
    spanned.model <- "gam0" %in% names(m)
    N <- m$N
    B <- t(m$B)  # J x N
    mats <- m$mats
    if (spanned.model) {
        L <- m$cL ## number of yield factors
        Omega <- m$Omega
        iM <- cbind(diag(N-L), matrix(0, N-L, L)) # Mt = iM * Zt
        ## in spanned model, macro variables come first
    } else {
        L <- m$R
        yieldFactorsOnly <- TRUE
        Omega <- m$Omega.Z
        iM <- cbind(matrix(0, N-L, L), diag(N-L))
        ## in unspanned model, yield factors come first
    }
    if (yieldFactorsOnly) {
        W <- m$W[1:L,]
        R <- L
    } else {
        W <- m$W[1:N,]
        R <- N
    }
    Cov.Y <- B %*% m$Cov.Z %*% t(B)
    Cov.X <- W %*% Cov.Y %*% t(W)
    Cov.ZX <- m$Cov.Z %*% t(W %*% B)
    Cov.MX <- iM %*% Cov.ZX
    Cov.M <- iM %*% m$Cov.Z %*% t(iM)

    ## return regression
    Phi1.h <- matrix.power(m$Phi, h)
    cBn <- function(n) -n * B[mats==n, ]
    cB <- matrix(rowMeans(sapply(mats[mats>h]-h, cBn)), N, 1)   # avg. of cB(n-h) across mats
    beta <- function(n)
        cBn(n-h) %*% Phi1.h - cBn(n) + cBn(h)
    beta.bar <- matrix(rowMeans(sapply(mats[mats>h], beta)), N, 1)  # true coefficients in reg on Z
    Cov.nu <- matrix(0, N, N) # Cov(Z(t+h) - E_t Z(t+h)) = Cov(nu(t+h))
    for (i in 1:h) {
        Phi.hmi <- matrix.power(m$Phi, h-i)
        Cov.nu <- Cov.nu + Phi.hmi %*% Omega %*% t(Phi.hmi)
    }
    Var.rx <- t(beta.bar) %*% m$Cov.Z %*% beta.bar + t(cB) %*% Cov.nu %*% cB  # Var(Erx) + Var(FE)
    Cov.rxX <- t(beta.bar) %*% Cov.ZX
    Cov.rxM <- t(beta.bar) %*% m$Cov.Z %*% t(iM)
    ## V = (Xtilde', M)'
    Cov.rxV <- cbind(Cov.rxX, Cov.rxM)

    for (sigma in sigmas) {
        Cov.Xtilde <- Cov.X + sigma^2 * W %*% t(W)
        ## Unspanned macro risk
        ## restricted forecast: rx(t+h) on X(t)
        R2ret.r <- Cov.rxX %*% solve(Cov.Xtilde) %*% t(Cov.rxX) / Var.rx
        ## unrestricted forecast: rx(t+h) on M(t) and X(t)
        Cov.V <- rbind(cbind(Cov.Xtilde, t(Cov.MX)),
                       cbind(Cov.MX, Cov.M))
        R2ret.ur <- R2ret.r
        try(R2ret.ur <- Cov.rxV %*% solve(Cov.V) %*% t(Cov.rxV) / Var.rx, silent=TRUE)
        ## Cov.V is non-invertible if M is spanned by X
        tblRisk[sigmas==sigma, ] <- c(R2ret.r, R2ret.ur, sqrt((1-R2ret.ur)/(1-R2ret.r)))

        ## role of measurement error
        if ((spanned.model) && !yieldFactorsOnly && sigma==sigmas[1]) {
            cat("# Table D1: variance of model-implied yield PCs\n")
            cat("The role of measurement error\n")
            cat(m$name, " -- sigma =", sigma*120000, "\n")
            tbl <- cbind(diag(Cov.Xtilde), diag(Cov.X), diag(sigma^2*W%*%t(W)))
            tbl <- cbind(tbl, tbl[,3]/tbl[,1])
            rownames(tbl) <- paste("PC", 1:nrow(Cov.X))
            colnames(tbl) <- c("Cov(Xtilde)", "Cov(X)", "sigma^2 WW'", "% meas. error.")
            print(round(tbl*100, digi=2))
            ## xtbl <- xtable(tbl*100, digits=2)
            ## print(xtbl, include.rownames=TRUE, include.colname=TRUE, only.contents=TRUE, sanitize.text.function=function(x){x})
        }

        for (i in 1:ncol(m$M.o)) {
            im <- matrix(0, N, 1); if (spanned.model) { im[i] <- 1 } else { im[L+i] <- 1} # sel. vector
            Cov.mX <- t(im) %*% Cov.ZX
            Var.m <- t(im) %*% m$Cov.Z %*% im

            ## Unspanned macro variation -- R^2 for reg. of m on X^tilde
            tblVar[sigmas==sigma, i] <- Cov.mX %*% solve(Cov.Xtilde) %*% t(Cov.mX)/Var.m


            ## Unspanned macro forecasts
            Cov.m1X <- t(im) %*% m$Phi %*% Cov.ZX # Cov(m(t+1), X(t))
            ## restricted forecast: m(t+1) on X(t)
            R2m.r <- Cov.m1X %*% solve(Cov.Xtilde) %*% t(Cov.m1X)/Var.m
            ## unrestricted forecast: m(t+1) on m(t) and X(t)
            Cov.m1m <- t(im) %*% m$Phi %*% m$Cov.Z %*% im
            ## V(t) = (m(t), X(t)')'
            Cov.m1V <- cbind(Cov.m1m, Cov.m1X)
            Cov.V <- rbind(cbind(Var.m, Cov.mX),
                           cbind(t(Cov.mX), Cov.Xtilde))
            R2m.ur <- R2m.r
            try(R2m.ur <- Cov.m1V %*% solve(Cov.V) %*% t(Cov.m1V)/Var.m, silent=TRUE)
            ## Cov.V is non-invertible if m is spanned by X

            ## relative RMSE
            tblForecasts[sigmas==sigma, i] <- sqrt((1-R2m.ur)/(1-R2m.r))

        }
    }
    cbind(tblVar, tblRisk, tblForecasts)
}

analyzeSpanning <- function(m, T=nrow(m$Y), M=1000, R=length(m$lamQ), sigmas=c(6, 0)/120000, h=12) {
    ## analyze spanning implications of estimated MTSMs using simulations
    ## Inputs:
    ##  m - model, list with estimated model parameters and data
    ##  T - sample size
    ##  M - number of replications
    ##  R - number of PCs included in regressions (maximum N for spanned, N-M for unspanned)
    ##  sigmas - SDs of measurement error -- carry out simulation for each value
    ## Output: table with
    ##  rows for (1) data, (2) median statistics, and (3) p-values
    ##   - (2) and (3) for each sigma
    ##  columns for each spanning-statistic -- see getSpanningStats
    ##  notes: - p-values are fraction of sim. samples with statistic below data
    ##         - except for 3rd (Wald) and 7th (dR^2) column: fraction of samples with statistics above data

    set.seed(616)
    N <- m$N
    if (R>length(m$lamQ))
        stop("Can't regress on more PCs than yield factors (near collinearity)\n")

    cat("## Model", m$name, "- regression on", R, "PCs of yields\n")
    modelSpanningPCs(m, R)

    ## Data
    PCs <- m$Y %*% t(m$W[1:R,])  # first R PCs
    data.stats <- getSpanningStats(m$M.o, PCs, m$Y, m$mats, h)

    ## Table for results
    tbl <- matrix(NA, 1+2*length(sigmas), length(data.stats))
    colnames(tbl) <- names(data.stats)
    tbl[1, ] <- data.stats
    row <- 2
    tmp <- cbind(m$name, paste0("R =", R),
                 c("Data", rep(c("mean", "p-val"), length(sigmas))),
                 c("", rep(round(sigmas*120000, digi=1), each=2)))
    rownames(tbl) <- apply(tmp, 1, function(row) paste(row, collapse="-"))
    ## Simulations
    for (sigma in sigmas) {
        ## get spanning results in simulated data
        sim.stats <- simSpanning(m, sigma, T, M, h, R)
        ## tbl[row, ] <- colMeans(sim.stats)
        tbl[row, ] <- apply(sim.stats, 2, median)
        tbl[row+1, ] <- colMeans(sim.stats < matrix(rep(data.stats, each=M), M, length(data.stats)))
        tbl[row+1, 3] <- mean(sim.stats[,3] > data.stats[3])
        tbl[row+1, 7] <- mean(sim.stats[,7] > data.stats[7])
        row <- row+2
    }
    ##print(tbl)
    return(tbl)
}

getSpanningStats <- function(M, PCs, Y, mats, h) {
    ## given macro/yield data set, run spanning regressions and obtain summary statistics
    statnames <- c(paste("UMV", colnames(M)),
                   "Wald", "pval", "R2.1", "R2.2", "R2.2-R2.1",  "RMSE",
                   paste("UMF", colnames(M)))
    stats <- setNames(numeric(length(statnames)), statnames)
    T <- nrow(Y)
    ## UMV - Unspanned Macro Variation - 2 statistics (R^2)
    for (i in 1:2) {
        lm.data <- lm(M[,i] ~ PCs)
        stats[i] <- summary(lm.data)$r.squared
    }

    ## UMR - Unspanned Macro Risks - 6 statistics on predictability of returns
    ydat <- rowMeans(sapply(mats[mats>h], function(n) getReturns(Y, n, h, mats)))
    xdat1 <- PCs[1:(T-h), ] ## PCs
    xdat2 <- M[1:(T-h),] ## macro variables
    lm.y <- lm(ydat ~ xdat1)
    lm.my <- lm(ydat ~ xdat1 + xdat2)
    ## if (all(is.na(tail(lm.my$coef, 2)))) {
    if (is.na(tail(lm.my$coef, 1))) {
        ## this is a shortcut
        ## - in some cases for SM, R=5, sigma=0, 5 PCs are collinear -- must be numerical issue
        stats[3] <- 0    # Wald test statistic
        stats[4] <- 1  # p-value
    } else {
        rval <- getTest(lm.my, lm.y)
        stats[3] <- rval$Chisq[2]
        stats[4] <- rval$Pr[2]
    }
    stats[5] <- summary(lm.y)$r.squared
    stats[6] <- summary(lm.my)$r.squared
    stats[7] <- summary(lm.my)$r.squared - summary(lm.y)$r.squared
    stats[8] <- sqrt(mean(lm.my$residuals^2))/sqrt(mean(lm.y$residuals^2))

    ## UMF - Unspanned Macro Forecasts - 2 statistics (rel. RMSEs)
    for (i in 1:2) {
        ydat <- M[2:T, i]
        xdat.m <- M[1:(T-1), i]
        xdat.y <- PCs[1:(T-1), ]
        ## lm.m <- lm(ydat ~ xdat.m)
        lm.y <- lm(ydat ~ xdat.y)
        lm.my <- lm(ydat ~ xdat.y + xdat.m)
        ## if (is.na(tail(lm.my$coef, 1))) {
        ##     stats[ind+i] <- NA
        ## } else {
        ##     SEs <- sqrt(diag(vcov.fn(lm.my)))
        ##     stats[ind+i] <- tail(abs(lm.my$coef)/SEs,1)
        ## }
        stats[8 + i] <- sqrt(mean(lm.my$residuals^2))/sqrt(mean(lm.y$residuals^2))
    }
    stats
}

simSpanning <- function(m, sigma, T, M, h, R) {
    ## simulate yield and macro data from model and run spanning regressions
    J <- length(m$mats)
    N <- m$N
    cat("* Simulating", M, "data sets of length", T, " from model", m$name, "\n")
    cat("* N =", N, " --  R =", R, " --  sigma =", sigma*120000, "bps\n")

    for (b in 1:M) {
        ## simulate yields
        Zsim <- simVAR1(T, m$mu, m$Phi, m$Sigma)
        Yhat <- rep(1, T) %*% m$A + Zsim %*% m$B
        ## i.i.d errors
        errors <- matrix(rnorm(J*T, mean=0, sd=sigma), T, J)
        ## robustness check: serially correlated errors
        ## rho <- 0.8
        ## sig2 <- (1-rho^2)*sigma^2 # innovation variance
        ## errors <- filter(matrix(rnorm(T*J, mean=0, sd=sqrt(sig2)), T, J),
        ##                  rho, method = "recursive",
        ##                  init = matrix(rnorm(J, mean=0, sd=sigma), 1, J))
        Ysim <- Yhat + errors

        if ("gam0" %in% names(m)) {
            ## spanned model
            Msim <- Zsim[, 1:2]
        } else {
            ## unspanned model
            Msim <- Zsim[, (N-1):(N)]
        }

        ## construct PCs from simulated yields
        ## estimated loadings:
        Wsim <- t(eigen(cov(Ysim))$vectors)[1:R,]
        ## fixed loadings:
        ## Wsim <- m$W[1:R,]
        PCsim <- Ysim %*% t(Wsim)

        stats <- getSpanningStats(Msim, PCsim, Ysim, m$mats, h)
        if (b==1) {
            sim.stats <- matrix(NA, M, length(stats))
            sim.stats[1, ] <- stats
        } else {
            sim.stats[b, ] <- stats
        }
    }
    ## print("check: SD of errors and fitted yields")
    ## tbl <- rbind(apply(errors, 2, sd),
    ##             apply(Yhat, 2, sd))
    ## print(round(120000*tbl, digi=2))

    return(sim.stats)
}

modelSpanningPCs <- function(m, R) {
    ## summarize cross-sectional spanning in the data
    cat("Regressing macro variable on", R, "PCs...\n")
    tbl <- matrix(NA, 2, 2)
    rownames(tbl) <- c("of observed yields", "of fitted yields")
    colnames(tbl) <- colnames(m$M.o)
    ## observed PCs
    PC <- m$Y %*% t(m$W[1:R,])
    for (i in 1:2)
        tbl[1, i] <- summary(lm(m$M.o[,i] ~ PC))$r.squared
    ## ## model-implied PN -- this only works for spanned models
    ## for (i in 1:2)
    ##     tbl[2, i] <- summary(lm(m$M.o[,i] ~ m$PN))$r.squared
    ## PCs of Yhat
    PC <- m$Yhat %*% t(eigen(cov(m$Yhat))$vectors[1:R,]) # t(m$W[1:N,])
    for (i in 1:2)
        tbl[2, i] <- summary(lm(m$M.o[,i] ~ PC))$r.squared
    print(round(tbl, 4))
}

testKnifeEdge <- function(SM, USM, flag.details=FALSE) {
    chisq <- -2*(USM$llk - SM$llk)
    cat("# Testing knife-edge restrictions of model", USM$name, "vs. spanned model", SM$name,"\n")
    cat("Chi-squared statistic:", chisq, "\n")
    ##df <- SM$cM*(1+SM$cL)
    J <- length(SM$mats)
    df <- (J - SM$cL)*(1+SM$cN) - 1 - SM$cL
    cat("Degrees of freedom:", df, "\n")
    cat("Critical value:", qchisq(.05, df, lower.tail=F), "\n")
    cat("p-value:", pchisq(chisq, df, lower.tail=F), "\n")

    rval <- list(SM.llk=SM$llk, SM.llkP=SM$llkP, SM.llkQ=SM$llkQ,
                 USM.llk=USM$llk, USM.llkP=USM$llkP, USM.llkQ=USM$llkQ,
                 chisq=chisq, df=df, cv=qchisq(.05, df), pval=pchisq(chisq, df, lower.tail=F))
    if (flag.details) {
        cN <- SM$cN
        cL <- SM$cL
        ## compare fit across maturities
        fit.tbl <- 120000*rbind(sqrt(colMeans((SM$Y-SM$Yhat)^2)),
                                sqrt(colMeans((USM$Y-USM$Yhat)^2)))
        rownames(fit.tbl) <- c("SM", "USM")
        colnames(fit.tbl) <- SM$mats
        print(round(fit.tbl, digi=2))

        ## comparison of models
        row.names <- c("llk", "P-llk", "Q-llk", "RMSE", "sigma.e",
                       sapply(1:cN, function(i) paste("P-ev", i)),
                       sapply(1:cN, function(i) paste("Q-ev", i)))
        tbl <- matrix(NA, length(row.names), 2)
        rownames(tbl) <- row.names
        colnames(tbl) <- c("SM", "USM")

        ## compare likelihoods
        tbl[1, ] <- c(SM$llk, USM$llk)
        tbl[2, ] <- c(SM$llkP, USM$llkP)
        tbl[3, ] <- c(SM$llkQ, USM$llkQ)
        row <- 4

        ## pricing errors
        tbl[row, ] <- 120000*c(sqrt(mean((SM$Y-SM$Yhat)^2)), sqrt(mean((SM$Y-USM$Yhat)^2)))
        tbl[row+1, ] <- 120000*c(SM$sigma.e, USM$sigma.e)
        row <- row+2
        ## compare parameter estimates
        ## P-eigenvalues
        tbl[row:(row+cN-1), ] <- cbind(sort(abs(eigen(SM$K1P)$values+1), dec=TRUE),
                                       sort(abs(eigen(USM$KP.ZZ)$values), dec=TRUE))
        row <- row+cN
        ## Q-eigenvalues
        tbl[row:(row+cN-1), 1] <- SM$lamQ+1
        tbl[row:(row+cL-1), 2] <- USM$lamQ
        print(round(tbl, digi=4))
    }
    return(rval)
}

getForwardRates <- function(m) {
    ## forward rates of arbitrary maturities
    ## Note: cannot get ACTUAL forward rates, only MODEL-IMPLIED
    T <- nrow(m$Y)
    mats4f <- 12*(1:10)
    h1 <- mats4f[c(2,5,9)]
    h2 <- mats4f[c(3,6,10)]
    h1.mat <- matrix(rep(h1, each=T), T, 3)
    h2.mat <- matrix(rep(h2, each=T), T, 3)
    makeForwards <- function(Y) {
        F <- (h2.mat * Y[,match(h2, mats4f)] - h1.mat * Y[,match(h1, mats4f)]) / (h2.mat-h1.mat)
        colnames(F) <- c("2-to-3y", "5-to-6y", "9-to-10y")
        return(F)
    }
    ## fitted
    if ("KQ.PP" %in% names(m)) {
        ## unspanned/JPS and yields-only/JSZ models
        ## -> fitted yields only function of yield factors
        N <- length(m$KQ.0P)
        loads <- gaussian.loadings(mats4f, m$KQ.0P, m$KQ.PP - diag(N), m$Omega.cP,
                                   m$rho0.cP, m$rho1.cP)
        Ytmp <- rep(1, T) %*% loads$A + m$cP %*% loads$B
    } else {
        ## spanned/JLS model
        loads <- gaussian.loadings(mats4f, m$muQ, m$PhiQ - diag(m$cN), m$Omega,
                                   m$rho0, m$rho1)
        Ytmp <- rep(1, T) %*% loads$A + m$Z %*% loads$B
    }
    Fhat <- makeForwards(Ytmp)
    ## risk-neutral
    if ("K1P.cP" %in% names(m)) {
        ## yields-only/JSZ
        loads <- gaussian.loadings(mats4f, m$K0P.cP, m$K1P.cP, m$Omega.cP,
                                   m$rho0.cP, m$rho1.cP)
        Ytmp <- rep(1, T) %*% loads$A + m$cP %*% loads$B
    } else if ("KP.ZZ" %in% names(m)) {
        ## unspanned/JPS
        loads <- gaussian.loadings(mats4f, m$KP.0Z, m$KP.ZZ - diag(m$N), m$Omega.Z,
                                   m$rho0.Z, m$rho1.Z)
        Ytmp <- rep(1, T) %*% loads$A + m$Z %*% loads$B
    } else {
        ## spanned/JLS
        loads <- gaussian.loadings(mats4f, m$mu, m$Phi - diag(m$cN),
                                   m$Omega, m$rho0, m$rho1)
        Ytmp <- rep(1, T) %*% loads$A + m$Z %*% loads$B
    }
    Frn <- makeForwards(Ytmp)
    ## forward term premia
    Ftp <- Fhat - Frn
    return(list(Fhat=Fhat, Frn=Frn, Ftp=Ftp))
}

## getForwardRates2 <- function(m) {
##     ## forward rates
##     ## -> required maturities need to be available in m$mats
##     T <- nrow(m$Y)
##     h1 <- c(24, 48, 108); h2 <- h1+12
##     if (any(is.na(match(c(h1, h2), m$mats))))
##         stop("getForwardRates: Maturity not available")
##     h1.mat <- matrix(rep(h1, each=T), T, length(h1))
##     h2.mat <- matrix(rep(h2, each=T), T, length(h2))
##     makeForwards <- function(Y) {
##         F <- (h2.mat * Y[,match(h2, m$mats)] - h1.mat * Y[,match(h1, m$mats)]) / (h2.mat-h1.mat)
##         colnames(F) <- c("2-to-3y", "4-to-5y", "9-to-10y")
##         return(F)
##     }
##     F <- makeForwards(m$Y)
##     Fhat <- makeForwards(m$Yhat)
##     Frn <- makeForwards(m$Yrn)
##     Ftp <- Fhat - Frn
##     return(list(F2=F, Fhat2=Fhat, Frn2=Frn, Ftp2=Ftp))
## }

analyzeFTP <- function(sm, usm, ym, flag.tofile=FALSE) {
    cat("Forward term premia -- model", sm$name, "\n")

    tbl <- matrix(NA, 3, 8)
    colnames(tbl) <- c("R^2 for PCs", "R^2 for all", "const", "b(level)", "b(slope)", "b(curve)", "b(cycle)", "b(infl)")
    rownames(tbl) <- colnames(usm$F)
    T <- nrow(usm$Y)
    Ftp.proj <- matrix(NA, T, 3)
    for (i in 1:3) {
        ## project FTP on yield-curve factors
        lm.r <- lm(usm$Ftp[,i] ~ usm$cP)
        tbl[i, 1] <- summary(lm.r)$r.squared
        Ftp.proj[,i] <- lm.r$fitted.values
        ## loadings on ALL factors
        lm.ur <- lm(usm$Ftp[,i] ~ usm$Z)
        tbl[i, 2] <- summary(lm.ur)$r.squared
        tbl[i, 3:8] <- lm.ur$coef
    }
    print(round(tbl, digi=4))

    ## focus on 2-to-3y-TP
    i <- 1
    ## correlations
    cat("Correlations for", colnames(usm$Ftp)[i], "FTP\n")
    M <- cbind(sm$Ftp[,i], usm$Ftp[,i], ym$Ftp[,i], Ftp.proj[,i])
    colnames(M) <- c("SM", "USM", "JSZ", "Proj")
    print(cor(M))

    ## plot FTP
    start.date <- c(1985, 1)
    ftp.sm <- ts(1200*sm$Ftp[,i], start=start.date, freq=12)
    ftp.usm <- ts(1200*usm$Ftp[,i], start=start.date, freq=12)
    ftp.jsz <- ts(1200*ym$Ftp[,i], start=start.date, freq=12)
    ##    ftp.proj <- ts(1200*Ftp.proj[,i], start=start.date, freq=12)
    cols <- c("black", "red", "blue", "green")
    lwds <- c(1,2,2,2)
    ltys <- c(1,2,1,1)
    if (flag.tofile) {
        filename <- paste0("figures/tp_", colnames(sm$M.o)[2], ".pdf")
        print(filename)
        ## postscript(filename, width=7, height=5,
        ##            horizontal=FALSE, onefile=FALSE, paper="special", pointsize=10)
        pdf(filename, width=7, height=5, pointsize=12)
    } else {
        dev.new()
    }
    par(mar=c(4,4,2,1))
    yrange <- c(-.5,6)
    plot(ftp.sm, type="l", ylim=yrange, col=cols[1], lwd=lwds[1], lty=ltys[1], ylab="Percent", xlab="Year",
         yaxp=c(0,6,6), yaxs="i")
    plot.recessions(yrange)
    lines(ftp.sm, col=cols[1], lwd=lwds[1], lty=ltys[1])
    lines(ftp.usm, col=cols[2], lwd=lwds[2], lty=ltys[2])
    lines(ftp.jsz, col=cols[3], lwd=lwds[3], lty=ltys[3])
    ##    lines(ftp.proj, col=cols[4], lwd=lwds[4], lty=ltys[4])
    abline(h=0, lty=2)
    ##    legend("top", c("SM", "USM", "JSZ", "USM on PCs"), col=cols, lty=ltys, lwd=lwds)
    legend("topright", c("SM(3,2)", "USM(3,2)", "yields-only"), col=cols, lty=ltys, lwd=lwds, bg="white")
    if (flag.tofile) {
        dev.off()
    } else {
        title(paste("Forward term premia - ", colnames(sm$M.o)[2], "-", colnames(sm$Ftp)[i]))
    }

}

testModelBased <- function(models) {
    makeTable <- function(rval) {
        tbl <- matrix(NA, 3, 3)
        rownames(tbl) <- c("\\quad LLK $SM(3,2)$", "\\quad LLK $USM(3,2)$", "\\quad Likelihood-ratio statistic")
        tbl[1,] <- unlist(rval[c("SM.llkQ", "SM.llkP", "SM.llk")])
        tbl[2,] <- unlist(rval[c("USM.llkQ", "USM.llkP", "USM.llk")])
        tbl[3,3] <- rval$chisq
        tbl
    }
    tbl1 <- makeTable(testKnifeEdge(models$SM.GRO, models$USM.GRO))
    tbl2 <- makeTable(testKnifeEdge(models$SM.UGAP, models$USM.UGAP))
    print(round(tbl1,digi=0))
    print(round(tbl2,digi=0))
}

makePopTbl <- function(tbls, i=1) {
    tbl.gro <- rbind(tbls$gro.usm, tbls$gro.sm3, tbls$gro.sm5)
    tbl.ugap <- rbind(tbls$ugap.usm, tbls$ugap.sm3, tbls$ugap.sm5)
    sigma.labels <- c("sigma=6bp", "sigma=0")
    rownames(tbl.gro) <- c(paste("USM(3,2), 3 PCs,", sigma.labels),
                           paste("SM(3,2), 3 PCs,", sigma.labels),
                           paste("SM(3,2), 5 PCs,", sigma.labels))
    rownames(tbl.ugap) <- rownames(tbl.gro)
    tbl1 <- cbind(tbl.gro[,1:2], tbl.ugap[,1:2])
    tbl2 <- cbind(tbl.gro[,3:5], tbl.ugap[,3:5])
    tbl3 <- cbind(tbl.gro[,6:7], tbl.ugap[,6:7])
    colnames(tbl1) <- c("GRO", "INF", "CPI", "UGAP")
    colnames(tbl2) <- rep(c("R2 w/o", "R2 w/", "RMSE"), 2)
    colnames(tbl3) <- c("GRO", "INF", "CPI", "UGAP")
    if (i==1) {
        tbl <- tbl1
    } else if (i==2) {
        tbl <- tbl2
    } else if (i==3) {
        tbl <- tbl3
    }
    print(round(tbl, digi=2))
    ## xtbl <- xtable(tbl, digits=3)
    ## print(xtbl, include.rownames=TRUE, include.colname=TRUE, only.contents=TRUE, sanitize.text.function=function(x){x})
}
