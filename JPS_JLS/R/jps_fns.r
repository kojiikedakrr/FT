jps.llk <- function (Y, M.o, W, R, KQ.XX, kinfQ=NA, KP.0Z=NA, KP.ZZ=NA, Omega.Z, mats, dt, sigma.e=NA, restr=0, gamma, Lam0, Lam1) {
    ## Compute the likelihood for a Gaussian term structure.
    ##
    ## INPUTS:
    ## Y   : (T+1)*J,  matrix of observed yields (first row are t=0 observations, which likelihood conditions on)
    ## mats       : 1*J,      maturities in years
    ## dt         : scalar,   length of period in years
    ## W          : N*J,      vector of portfolio weights to fit without error.
    ## KQ.XX      : N*N,      normalized latent-model matrix (does not have to be diagonal, see form below)
    ## kinfQ      : scalar,   when the model is stationary, the long run mean of the annualized short rate under Q is -kinfQ/K1(m1,m1)
    ## Omega.cP   : N*N,      positive definite matrix that is the covariance of innovations to cP
    ##
    ## OPTIONAL INPUTS -- concentrated out if not supplied:
    ## KP.0Z      : N*1       intercept in VAR for Z
    ## KP.ZZ      : N*N       mean reversion matrix in VAR for Z
    ## sigma.e    : scalar    standard error of yield observation errors
    ## restr      : scalar    what type of restrictions on market prices of risk or VAR params
    ##                        0 -- no restrictions
    ##                        1 -- zero restrictions indicated by vector gamma
    ##                        2 -- reduced rank restriction (JSZ)
    ##                        3 -- eigenvalue restr. lam_1^Q = lam_1^P (JPS)
    ##                        planned -- zero restrictions and eigenvalue restriction
    ## gamma  : N^2*1     for case 1 -- vector with indicators for zero restrictions on Lam1
    ## ev.max     : scalar    for cases 3 and 4 -- largest eigenvalue of Phi -- NA means Q-ev
    ##
    ## OUTPUT:
    ## llk        : T*1       time series of -log likelihoods (includes 2-pi constants)
    ## AcP        : 1*J       yt = AcP' + BcP'*Xt  (yt is J*1 vector)
    ## BcP        : N*J       AcP, BcP satisfy internal consistency condition that AcP*W' = 0, BcP*W' = I_N
    ## AX         : 1*J       yt = AX' + BX'*Xt
    ## BX         : N*J       Xt is the 'jordan-normalized' latent state
    ## ...
    ##
    ## The model takes the form:
    ##   r(t) = rho0.cP + rho1.cP'*cPt
    ##        = 1'*Xt  (Xt is the 'jordan-normalized' state
    ##        = 1 period discount rate (annualized)
    ##
    ## Under Q:
    ##   X(t+1)  = KQ.0X + KQ.XX*X(t)  + eps_X(t+1),   cov(eps_X(t+1)) = Omega.X
    ##   cP(t+1) = KQ.0P + KQ.PP*cP(t) + eps_cP(t+1),  cov(eps_cP(t+1)) = Omega.cP
    ##   where Omega.X is chosen to match Omega.cP
    ## and K0Q_X(m1) = kinfQ where m1 is the multiplicity of the highest eigenvalue (typically 1)

    ## Under P:
    ##   Z(t+1) = KP.0Z + KP.ZZ*Z(t) + eps_Z(t+1),  cov(eps_cP(t+1)) = Omega.Z

    ## Setup
    T <- nrow(Y)-1
    J <- ncol(Y)
    if (R != nrow(KQ.XX))
        stop("R not equal to nrow(KQ.XX)")
    if ((ncol(W)!=J))
        stop("W should have J columns")
    WR <- matrix(W[1:R, ], R, J)
    N <- R + ncol(M.o)
    cP <- Y %*% t(WR) # (T+1)*R
    Z <- cbind(cP, M.o)
########################################################################

########################################################################
    ## COMPUTE THE Q-LIKELIHOOD:
    ## First find the loadings for the model:
    ## yt = AcP' + BcP'*cPt, AcP is 1*J, BcP is N*J

    Omega.cP <- Omega.Z[1:R, 1:R]
    K1Q.X <- KQ.XX - diag(R)   ## JSZ parameterization

    if (is.na(kinfQ)) {
        ## concentrate out kinfQ
        ## AcP = alpha0_cP*kinf + alpha1_cP
        rho0.cP <- 0
        ## AcP0, AX0 will be the loadings with rho0_cP = 0, which won't be correct
        loads <- jsz.loadings.rho0cP(WR, K1Q.X, rho0.cP, Omega.cP, mats, dt)
        ## [BcP, AcP0, K0Q_cPx, K1Q_cP, rho0_cP, rho1_cP, K0Q_X, K1Q_X, AX0, BX, Sigma_X, alpha0_cP, alpha1_cP, alpha0_X, alpha1_X, m1]
        BcP <- loads$BcP; BX <- loads$BX
        alpha0.X <- loads$alpha0.X; alpha1.X <- loads$alpha1.X
        alpha0.cP <- loads$alpha0.cP; alpha1.cP <- loads$alpha1.cP
        ##K1Q.X <- loads$K1Q.X;
        m1 <- loads$m1;
        ## back out kinfQ that fits average yields
        require(MASS)
        V <- t(Null(t(WR)))
        kinfQ <- t(colMeans(Y[2:(T+1),]) - t(alpha1.cP) - t(BcP)%*%colMeans(cP[2:(T+1),]))%*%(t(V)%*%V%*%t(alpha0.cP)) / (alpha0.cP%*%t(V)%*%V%*%t(alpha0.cP))
        kinfQ <- as.numeric(kinfQ)

        ## get correct loadings
        AX <- alpha0.X*kinfQ + alpha1.X;
        AcP <- alpha0.cP*kinfQ + alpha1.cP;

        ## get these to return to caller (not used in jsz.llk)
        K0Q.X <- matrix(0,R,1);
        K0Q.X[m1] <- kinfQ;
        params <- jsz.rotation(WR, K1Q.X, K0Q.X, dt, BX, AX);
        KQ.0X <- K0Q.X
        KQ.0P <- params$K0Q.cP; KQ.PP <- params$K1Q.cP + diag(R);
        rho0.cP <- params$rho0.cP; rho1.cP <- params$rho1.cP
    } else {
        loads <- jsz.loadings(WR, K1Q.X, kinfQ, Omega.cP, mats, dt)
        BcP <- loads$BcP; AcP <- loads$AcP; KQ.0P <- loads$K0Q.cP; KQ.PP <- loads$K1Q.cP + diag(R);
        rho0.cP <- loads$rho0.cP; rho1.cP <- loads$rho1.cP; KQ.0X <- loads$K0Q.X; ##KQ.XX <- loads$K1Q.X;
        AX <- loads$AX; BX <- loads$BX;
    }
    yields.m <- rep(1,T+1)%*%AcP + cP %*% BcP # (T+1)*J, model-implied yields
    yield.errors <- Y[2:(T+1),] - yields.m[2:(T+1),]; # T*J
    square_orthogonal_yield.errors <- yield.errors^2; # T*J

    ## Compute optimal sigma.e if it is not supplied
    if (is.na(sigma.e))
        sigma.e <- sqrt( sum(square_orthogonal_yield.errors)/(T*(J-R)) )

    term1 <- .5*rowSums(square_orthogonal_yield.errors)/sigma.e^2
    term2 <- (J-R)*.5*log(2*pi)
    term3 <- .5*(J-R)*log(sigma.e^2)
    llkQ <- term1 + term2 + term3 # 1*T

    ## MATLAB: NP=>N=5, NQ=>R=3, M=>J
    ## - 1/2*T*(M-NQ)*log(2*pi) => term2
    ## - 1/2*T*sum(log(model.derived.sPC)) => term3
    ## - 1/2*T*(M-NQ); => term1

    ## note: sum(term1)=.5*T*(J-R)

########################################################################

########################################################################
    ## COMPUTE THE P-LIKELIHOOD:

    if (any(is.na(KP.0Z))||any(is.na(KP.ZZ))) {
        ## VAR parameters not supplied
        if (restr==0) {
            ## unrestricted MPR -- run OLS for unconstrained VAR
            var1 <- ar.ols(Z, order=1, aic=FALSE, demean=FALSE, intercept=TRUE)
            KP.0Z <- var1$x.intercept
            KP.ZZ <- var1$ar[,,]
###########################
        } else if (restr==1) {
            ## zero restrictions on MPR
            if (missing(Lam0)|missing(Lam1)) {
                ## Lam0/Lam1 not provided -- concentrated out of LLK
                if (missing(gamma))
                    stop("jps.llk: when restr=1 need to provide either (Lam0,Lam1) or gamma")
                if (length(gamma)!=R*(1+N))
                    stop("jps.llk: number of indicators not equal to elements of [Lam0, Lam1]")
                if (sum(gamma)>0) {
                    mats <- getLambdaJPS(gamma, cP, M.o, KQ.0P, KQ.PP, solve(Omega.cP))
                    Lam0 <- mats$Lam0
                    Lam1 <- mats$Lam1
                } else {
                    Lam0 <- matrix(0, R, 1)
                    Lam1 <- matrix(0, R, N)
                }
            } else {
                if (missing(gamma))
                    warning("jsz.llk: restr=1, Lam0/Lam1 and gamma provided -> ignoring gamma")
            }

            ## construct P-dynamics
            KP.0P <- KQ.0P + Lam0
            KP.PZ <- cbind(KQ.PP, matrix(0, R, N-R)) + Lam1
            ## obtain KP.0M and KP.MZ from regression of M(t) to Z(t-1)
            var1 <- ar.ols(Z, order=1, aic=FALSE, demean=FALSE, intercept=TRUE)
            KP.0M <- tail(var1$x.intercept, N-R)
            KP.MZ <- var1$ar[,,][(R+1):N,]
            KP.0Z <- matrix(c(KP.0P, KP.0M), N, 1)
            KP.ZZ <- rbind(KP.PZ, KP.MZ)
            ## debug
            ## print(Lam0)
            ## print(Lam1)

###########################
        }  else {
            stop("not yet implemented")
        }

    } else {
        if (!restr==0) stop("can't provide VAR params AND have restrictions on MPR!")
    }

    innovations = t(Z[2:(T+1),]) - (KP.0Z %*% matrix(1,1,T) + KP.ZZ %*% t(Z[1:T,])) # N*T

    term1 <- .5*N*log(2*pi)
    term2 <- .5*log(det(Omega.Z))
    term3 <- .5*colSums(innovations*solve(Omega.Z, innovations))
    llkP = term1 + term2 + term3  # 1*T

    ## MATLAB
    ## - 1/2*(T-1)*NP*log(2*pi) => term1
    ## - 1/2*(T-1)*2*log(det(model.canonical.L)) +. term2
    ## - 1/2*trace(eX.'*eX) => term3 (?)

########################################################################

    return(list(llk=t(llkQ + llkP), AcP=AcP, BcP=BcP, AX=AX, BX=BX,
                KP.0Z=KP.0Z, KP.ZZ=KP.ZZ, sigma.e=sigma.e,
                KQ.0P=KQ.0P, KQ.PP=KQ.PP, rho0.cP=rho0.cP, rho1.cP=rho1.cP,
                KQ.0X=KQ.0X, KQ.XX=KQ.XX, K1Q.X=K1Q.X, cP=cP, llkP=llkP, llkQ=llkQ))

}

objJPS <- function(theta, Y, M.o, W, R, mats, gamma) {
    ## objective function is sum of negative log likelihoods
    N <- R + ncol(M.o)
    pars <- theta2pars.jps(theta, R, N)

    if (missing(gamma))
        stop("objJPS: need to provide gamma")

    ## check restrictions on parameter space
    valid <- TRUE
    ## diagonal elements of Sigma positive and bounded away from zero
    if (any(diag(pars$L.Z)<1e-7)) valid <- FALSE
    if (any(diag(pars$L.Z)>1)) valid <- FALSE
    ## eigenvalues of Phi.Q not explosive
    if (any(pars$lamQ>1)) valid <- FALSE
    ## eigenvalues sorted in decreasing order
    if (any(pars$dlamQ>0)) valid <- FALSE

    ## if parameters satisfy restriction on param space
    if (valid) {
        ## evaluate likelihood function and return sum of negative logliks
        res.llk <- jps.llk(Y, M.o, W, R, KQ.XX=diag(pars$lamQ),
                           Omega.Z=pars$Omega.Z, mats=mats, dt=1,
                           restr=1, gamma=gamma)
        return(sum(res.llk$llk))
    } else {
        ## else return penalty value
        return(1e6)
    }
}

theta2pars.jps <- function(theta, R, N) {
    ## convert theta vector to list of individual parameters
    ## Q parameters
    pars <- list(dlamQ=theta[1:R])
    pars$lamQ=cumsum(pars$dlamQ+c(1, rep(0, R-1)))
    ## P-innovation covariance matrix
    pars$L.Z <- matrix(0,N,N)
    pars$L.Z[lower.tri(pars$L.Z, diag=TRUE)] <- tail(theta, -R)
    pars$Omega.Z <- pars$L.Z %*% t(pars$L.Z)
    return(pars)
}

pars2theta.jps <- function(pars) {
    ## convert list of individual parameters to theta vector
    dlamQ <- c(pars$lamQ[1]-1,diff(pars$lamQ));
    return(c(dlamQ, pars$L.Z[lower.tri(pars$L.Z, diag=TRUE)]))
}

gamma2pars.jps <- function(gamma, R, N) {
    ## convert gamma vector to list of individual parameters
    ## same as theta2pars but with all parameters, and without reparameterization
    pars <- list(kinfQ=gamma[1], lamQ=gamma[2:(R+1)]); gamma <- tail(gamma, -R-1)
    pars$KP.0Z <- matrix(gamma[1:N], N, 1); gamma <- tail(gamma, -N)
    pars$KP.ZZ <- matrix(gamma[1:N^2], N, N); gamma <- tail(gamma, -N^2)
    pars$L.Z <- matrix(0, N, N)
    pars$L.Z[lower.tri(pars$L.Z, diag=TRUE)] <- gamma[1:(N*(N+1)/2)]
    pars$Omega.Z <- pars$L.Z %*% t(pars$L.Z)
    pars$sigma.e <- tail(gamma, 1)
    return(pars)
}

pars2gamma.jps <- function(pars)
    ## convert list of individual parameters to gamma vector
    c(pars$kinfQ, pars$lamQ, pars$KP.0Z, c(pars$KP.ZZ), pars$L.Z[lower.tri(pars$L.Z, diag=TRUE)], pars$sigma.e)
## length: 1+R+N+N^2+N*(N+1)/2+1 = 1+3+30+15+1 = 50

estJPS <- function(Y, M.o, W, R, mats, gamma) {
    ## Obtain MLE
    ## Arguments:
    ##  gamma -- risk price specification: vector with length N+N^2 indicating which elements of lambda are unrestricted (1) or restricted to zero (0)
    ## Value:
    ##  pars -- list with parameter estimates

    getStartingValuesForMLE <- function(L.Z) {
        ## obtain starting values for lamQ for MLE
        ##  -> random seeds for Q-eigenvalues
        ##
        ## Arguments: L.Z
        ## Value: list with starting values

        if (missing(L.Z))
            error("L.Z needs to be provided")
        Omega.Z <- L.Z %*% t(L.Z)
        ## nSeeds <- 500;  # how many random seeds to try
        ## best.llk <- Inf
        ## for (i in 1:nSeeds) {
        ##     (lamQ <- 1-abs(sort(runif(R, 0, .1))))
        ##     res.llk <- jps.llk(Y, M.o, W, R, KQ.XX=diag(lamQ), Omega.Z = Omega.Z,
        ##                        mats=mats, dt=1)
        ##     llk <- sum(res.llk$llk)
        ##     if (llk<best.llk) {
        ##         cat('Improved seed llk to ', llk, '\n')
        ##         best.llk <- llk
        ##         best.lamQ <- lamQ
        ##     }
        ## }
        ##best.lamQ <- c(0.9971, 0.9650, 0.8868)
        best.lamQ <- c(0.997, 0.95, 0.9)[1:R]
        return(list(lamQ=best.lamQ, L.Z=L.Z))
    }

    N <- R + ncol(M.o)
    J <- ncol(Y)
    cP <- Y %*% t(W[1:R,])
    Z <- cbind(cP, M.o)

    if (missing(mats)||length(mats)!=J)
        stop("estML: mats needs to be provided and have length consistent with yield data")
    if (missing(W)||ncol(W)!=J)
        stop("estML: W needs to be provided and have dimensions consistent with yield data")

    if (missing(gamma))
        gamma <- rep(1, R*(N+1))

    cat("*** MLE ***\n")

    ## (1) estimate VAR -- just to get Sigma.hat
    lm <- ar.ols(Z, aic=FALSE, order.max=1, intercept=TRUE, demean=FALSE)

    Omega.Z <- lm$var.pred
    L.Z <- t(chol(Omega.Z))

    ## numerical optimization
    pars.start <- getStartingValuesForMLE(L.Z)
    theta.start <- pars2theta.jps(pars.start)
    cat("at starting values:", objJPS(theta.start, Y, M.o, W, R, mats, gamma), "\n")
    theta <- getOptim(theta.start, objJPS, Y, M.o, W, R, mats, gamma)
    cat("at optimal values: ", objJPS(theta, Y, M.o, W, R, mats, gamma), "\n")
    pars <- theta2pars.jps(theta, R, N)  ## lamQ, Omega

    pars$gamma <- gamma
    pars$Omega.cP <- pars$Omega.Z[1:R, 1:R]

    ## value of likelihood
    res.llk <- jps.llk(Y, M.o, W, R, KQ.XX=diag(pars$lamQ), Omega.Z=pars$Omega.Z, mats=mats, dt=1, restr=1, gamma=gamma)
    pars$kinfQ <- res.llk$KQ.0X[which(res.llk$KQ.0X!=0)]
    pars$KP.0Z <- res.llk$KP.0Z
    pars$KP.ZZ <- res.llk$KP.ZZ
    pars$sige2 <- res.llk$sigma.e^2
    pars$llkP <- -sum(res.llk$llkP)
    pars$llkQ <- -sum(res.llk$llkQ)
    pars$llk <- pars$llkP + pars$llkQ

    return(pars)
}


getLambdaJPS <- function(gamma, cP, M, KQ.0P, KQ.PP, OmegaInv) {
    ## restricted VAR estimation to obtain risk prices
    ## KQ.0P and KQ.PP are givenn
    ## restriction: vec( [KP.0P, KP.PZ] = R * lambda + vec( [KQ.0P, KQ.PP, 0(3x2)] )
    ## estimate lambda -- Luetkepohl (5.2.6)
    R <- ncol(cP)
    N <- R + ncol(M)
    T <- nrow(cP)

    ## restriction:  beta = R*lambda + r
    S <- matrix(0, R*(N+1), sum(gamma==1))
    S[ cbind(which(gamma==1), seq(1, sum(gamma))) ] <- 1
    r <- c( KQ.0P, KQ.PP, matrix(0, R, N-R) )    ## vec( [KQ.0P, KQ.PP, 0(3x2)] )

    Z <- t(cbind(1, cbind(cP, M)[1:(T-1),]));  ## (1+N) x (T-1)
    Y <- t(cP[2:T,])                          ## R x (T-1)
    y <- as.numeric(Y)                         ## N*(T-1) x 1

    ## just for fun -- unrestricted VAR identical to c(mu.hat, as.numeric(Phi.hat))
    ## beta.hat <- kronecker(solve(tcrossprod(Z))%*%Z, diag(N)) %*% y

    z <- y - kronecker(t(Z), diag(R)) %*% r
 ##                    T-1 x (N+1)   RxR
    ##                  R(T-1) x R(N+1)
    inv.cov.mat <- t(S) %*% kronecker(tcrossprod(Z), OmegaInv) %*% S
    cov.mat <- solve(inv.cov.mat)
    lambda.hat <- as.numeric( cov.mat %*% t(S) %*% kronecker(Z, OmegaInv) %*% z )
    lambda.full <- S %*% lambda.hat
    return(list(lambda.hat=lambda.hat, Lam0=lambda.full[1:R], Lam1=matrix(tail(lambda.full, -R), R, N), cov.mat=cov.mat, inv.cov.mat=inv.cov.mat, R=R))
}

loadJPSmodel <- function(name, filename) {
    # require(jsz)
    cat("*** loading JPS model", name, "\n")
    cat("*** results file:", filename, "\n")
    load(filename) # pars, est.llk, Y, W, M.o, R, N, mats,
    m <- append(pars, est.llk)
    m$name <- name
    m$Y <- Y
    m$M.o <- M.o
    m$W <- W
    m$N <- N
    m$R <- R
    m$mats <- mats
    ## additional objects
    m$Z <- cbind(m$cP, m$M.o)
    m$mu <- m$KP.0Z
    m$Phi <- m$KP.ZZ
    m$Sigma <- m$L.Z
    ## loadings on Z
    m$A <- m$AcP
    m$B <- rbind(m$BcP, 0, 0)  ## zero loadings on macro variables
    ## fitted yields
    T <- nrow(Y)
    m$Yhat <- rep(1, T) %*% m$AcP + m$cP %*% m$BcP
    ## term premium
    m$rho0.Z <- m$rho0.cP
    m$rho1.Z <- c(m$rho1.cP, 0, 0)
    loads.rn <- gaussian.loadings(m$mats, m$KP.0Z, m$KP.ZZ - diag(N),
                                  m$Omega.Z, m$rho0.Z, m$rho1.Z)
    m$Yrn <- rep(1, T) %*% loads.rn$A + m$Z %*% loads.rn$B
    m$Ytp <- m$Yhat - m$Yrn
    ## population moments
    m$Cov.Z <- matrix(solve(diag(N^2) - kronecker(m$Phi, m$Phi)) %*% as.numeric(m$Omega.Z), N, N)
    m$E.Z <- as.vector(solve(diag(N) - m$Phi) %*% m$mu)
    m <- append(m, getForwardRates(m))
##    m <- append(m, getForwardRates2(m))
    return(m)
}

analyzeTaylorRule <- function(m, B=500, T=ncol(m$Y)) {
    ## analyze Taylor-rule implications of JPS models
    ## -> simulate from models and estimate Taylor rules
    require(sandwich)
    R2 <- numeric(B)
    coef <- matrix(NA, B, 2)
    tstats <- matrix(NA, B, 2)
    pvals <- matrix(NA, B, 2)
    J <- ncol(m$Y)

    ## a) regression of model-implied short rate on M.t
##    r.data <- 12*(m$Y[,1]) # 6m rate
    r.data <- 12*(m$rho0.Z + m$Z %*% m$rho1.Z)
    lm.data <- lm(r.data ~ m$M)

    getBeta <- function(Cov.Z) {
        Cov.M <- Cov.Z[4:5, 4:5]
        Cov.MP <- Cov.Z[4:5, 1:3]
        return(solve(Cov.M) %*% Cov.MP %*% m$rho1.Z[1:3])
    }

    ## b) calculate beta using sample moments
    beta.smpl <- getBeta(cov(m$Z))
    ## c) calculate beta using population moments
    cN <- 5
    beta.pop <- getBeta(matrix(solve(diag(cN^2) - kronecker(m$KP.ZZ, m$KP.ZZ)) %*% as.numeric(m$Omega.Z), cN, cN))

    ## d) regression on simulated yields
    for (b in 1:B) {
        ## simulate
        Zsim <- simVAR1(T, m$KP.0Z, m$KP.ZZ, m$L.Z)
        ##Yhat <- rep(1, T) %*% m$AcP + Zsim[,1:3] %*% m$BcP
        ##errors <- matrix(rnorm(J*T, mean=0, sd=sqrt(m$sige2)), T, J)
        ##Ysim <- Yhat + errors
        ## short rate
        rsim <- 12*(m$rho0.Z + Zsim %*% m$rho1.Z)
        ## alternative: fitted or actual short yield
        ## macro vars
        Msim <- Zsim[, 4:5]
        ## run Taylor-rule regression
        lm1 <- lm(rsim ~ Msim)
        R2[b] <- summary(lm1)$r.squared
        coef[b,] <- lm1$coef[2:3]
        tstats[b,] <- (lm1$coef/sqrt(diag(vcovHAC(lm1))))[2:3]
        pvals[b,] <- pt(abs(tstats[b,]), df=nrow(m$Y)-3, lower.tail=FALSE)*2
    }
    tbl <- createTbl(c("data", "data-tstat", "analytical smpl.", "analytical pop.",
                       "sim. mean", "sim. mean-tstat", "5%", "95%", "% sig"),
                     c("b(cycle)", "b(infl)", "R^2"))
    tbl[1,] <- c(lm.data$coef[2:3], summary(lm.data)$r.squared)
    tbl[2, 1:2] <- lm.data$coef[2:3]/sqrt(diag(vcovHAC(lm.data)))[2:3]
    tbl[3, 1:2] <- beta.smpl*12
    tbl[4, 1:2] <- beta.pop*12
    tbl[5,] <- colMeans(cbind(coef, R2))
    tbl[6,1:2] <- colMeans(tstats)
    tbl[7:8,] <- apply(cbind(coef, R2), 2, quantile, c(0.05, 0.95))
    tbl[9,1:2] <- colMeans(pvals<.05)
    print(round(tbl, digi=4))
    ##print(summary(pvals))

    ## plot
    ## (1) T-Bill rate
    ## (2) fitted short rate
    ## (3) desired rate
    start.date <- c(1985, 1)
    r.TBill <- ts( m$Y[,1]*1200, start=start.date, frequency=12)
    r.fitted <- ts( r.data*100, start=start.date, frequency=12)
    r.taylor <- ts( lm.data$fitted*100, start=start.date, frequency=12)
    dev.new()
    ##postscript(paste("figures/taylor_rule.eps",sep=""), horizontal=FALSE, onefile=FALSE, paper="special", width=7, height=6)
    ##pdf("figures/taylor_rule.pdf", width=7, height=5, pointsize=12)
    par(mar = c(4,4,2,1)+.1)
    yrange <- range(0, r.TBill, r.taylor, r.fitted)
    lwds <- c(1,2,2,2)
    ltys <- c(1,1,2,4)
    colors <- c("black", "blue", "red")
    plot(r.TBill, type='l', ylim=yrange, xlab="Year", ylab="Percent", col=colors[1], lwd=lwds[1])
    plot.recessions(yrange)
    lines(r.TBill, col=colors[1], lwd=lwds[1], lty=ltys[1])
    lines(r.fitted, col=colors[2], lwd=lwds[2], lty=ltys[2])
    lines(r.taylor, col=colors[3], lwd=lwds[3], lty=ltys[3])
    legend("bottomleft", c("3m T-Bill rate", "fitted short rate", "Taylor rule"), lwd=lwds, col=colors, lty=ltys, bg="white")
    ##dev.off()
}

analyzeSpannedTaylorRule <- function(m, B=500, T=ncol(m$Y)) {
    ## analyze implications of JPS models for SPANNED Taylor rule
    require(sandwich)
    R2 <- numeric(B)
    coef <- matrix(NA, B, 2)
    tstats <- matrix(NA, B, 2)
    pvals <- matrix(NA, B, 2)
    J <- ncol(m$Y)

    ## a) regression of model-implied short rate on M.span
    r.data <- 12*(m$rho0.Z + m$Z %*% m$rho1.Z)
    M1.span <- lm(m$M[,1] ~ m$Z[,1:3])$fitted.values
    M2.span <- lm(m$M[,2] ~ m$Z[,1:3])$fitted.values
    lm.data <- lm(r.data ~ M1.span + M2.span)

    ## analytical: rotation of risk factors into M.span and U (orthogonal)
    getBeta <- function(Cov.Z) {
        Cov.M <- Cov.Z[4:5, 4:5]
        Cov.MP <- Cov.Z[4:5, 1:3]
        Cov.P <- Cov.Z[1:3, 1:3]
        gam1 <- Cov.MP %*% solve(Cov.P)  ## 2x3 * 3x3 = 2x3
        ## orthogonality: d' Cov(PN.t) gam1' = 0
        F <- Cov.P %*% t(gam1)
        d <- CrossProduct3D(F[,1], F[,2])
        d <- d/d[1]
        ##print(Cov.P)
        ##print(gam1)
        ##print(F)
        D <- rbind(gam1, d)
        print(D)
        rho1.star <- m$rho1.Z[1:3] %*% solve(D)
        return(rho1.star[1:2])
    }

    ## b) calculate beta using sample moments
    beta.smpl <- getBeta(cov(m$Z))
    ## c) calculate beta using population moments
    cN <- 5
    beta.pop <- getBeta(matrix(solve(diag(cN^2) - kronecker(m$KP.ZZ, m$KP.ZZ)) %*% as.numeric(m$Omega.Z), cN, cN))

    ## d) regression on simulated yields
    for (b in 1:B) {
        ## simulate
        Zsim <- simVAR1(T, m$KP.0Z, m$KP.ZZ, m$L.Z)
        ##Yhat <- rep(1, T) %*% m$AcP + Zsim[,1:3] %*% m$BcP
        ##errors <- matrix(rnorm(J*T, mean=0, sd=sqrt(m$sige2)), T, J)
        ##Ysim <- Yhat + errors
        ## short rate
        rsim <- 12*(m$rho0.Z + Zsim %*% m$rho1.Z)
        ## alternative: fitted or actual short yield
        ## macro vars
        Msim <- Zsim[, 4:5]
        ## run Taylor-rule regression
        M1.span <- lm(Msim[,1] ~ Zsim[,1:3])$fitted.values
        M2.span <- lm(Msim[,2] ~ Zsim[,1:3])$fitted.values
        lm1 <- lm(rsim ~ M1.span + M2.span)
        R2[b] <- summary(lm1)$r.squared
        coef[b,] <- lm1$coef[2:3]
        tstats[b,] <- (lm1$coef/sqrt(diag(vcovHAC(lm1))))[2:3]
        pvals[b,] <- pt(abs(tstats[b,]), df=nrow(m$Y)-3, lower.tail=FALSE)*2
    }
    tbl <- createTbl(c("data", "data-tstat", "analytical smpl.", "analytical pop.",
                       "sim. mean", "sim. mean-tstat", "5%", "95%", "% sig"),
                     c("b(cycle)", "b(infl)", "R^2"))
    tbl[1,] <- c(lm.data$coef[2:3], summary(lm.data)$r.squared)
    tbl[2, 1:2] <- lm.data$coef[2:3]/sqrt(diag(vcovHAC(lm.data)))[2:3]
    tbl[3, 1:2] <- beta.smpl*12
    tbl[4, 1:2] <- beta.pop*12
    tbl[5,] <- colMeans(cbind(coef, R2))
    tbl[6,1:2] <- colMeans(tstats)
    tbl[7:8,] <- apply(cbind(coef, R2), 2, quantile, c(0.05, 0.95))
    tbl[9,1:2] <- colMeans(pvals<.05)
    print(round(tbl, digi=4))
    ##print(summary(pvals))
}

getTblRules <- function(m, macro.name, span=FALSE) {
    ## analyze Taylor-rule implications of JPS models
    require(sandwich)
    J <- ncol(m$Y)

    ## a) regression of model-implied short rate on M.t
##    r.data <- 12*(m$rho0.Z + m$Z %*% m$rho1.Z)
    r.data <- 12*(m$Y[,1]) # 6m rate
    if (span) {
        ## regression on spanned M.t
        M1.span <- lm(m$M[,1] ~ m$Z[,1:3])$fitted.values
        M2.span <- lm(m$M[,2] ~ m$Z[,1:3])$fitted.values
        lm.data <- lm(r.data ~ M1.span + M2.span)
    } else {
        lm.data <- lm(r.data ~ m$M)
    }

    getPars <- function(Cov.Z, E.Z) {
        Cov.M <- Cov.Z[4:5, 4:5]
        Cov.MP <- Cov.Z[4:5, 1:3]
        beta <- solve(Cov.M) %*% Cov.MP %*% m$rho1.Z[1:3]
        E.r <- m$rho0.Z + m$rho1.Z %*% E.Z
        alpha <- E.r - crossprod(beta, E.Z[4:5])
        explVar <- t(beta) %*% Cov.Z[4:5, 4:5] %*% beta
        totVar <- m$rho1.Z %*% Cov.Z %*% m$rho1.Z
        return(list(alpha=alpha, beta=beta, R2=explVar/totVar))
    }

    getParsSpanned <- function(Cov.Z, E.Z) {
        ## rotation of risk factors into M.span and U (orthogonal): L = (M.span, U)
        Cov.M <- Cov.Z[4:5, 4:5]
        Cov.MP <- Cov.Z[4:5, 1:3]
        Cov.P <- Cov.Z[1:3, 1:3]
        gam1 <- Cov.MP %*% solve(Cov.P)  ## 2x3 * 3x3 = 2x3
        gam0 <- E.Z[4:5] - gam1 %*% E.Z[1:3]
        ## orthogonality: d' Cov(PN.t) gam1' = 0
        F <- Cov.P %*% t(gam1)
        d <- CrossProduct3D(F[,1], F[,2])
        d <- d/d[1]
        D <- rbind(gam1, d)
        rho1.star <- m$rho1.Z[1:3] %*% solve(D)  ## r = rho0.star + rho1.star'L
        beta <- rho1.star[1:2]                   ## r = const + beta'M.span + resid
        E.r <- m$rho0.Z + m$rho1.Z %*% E.Z
        alpha <- E.r - t(beta) %*% (gam0 + gam1 %*% E.Z[1:3]) ## beta' gam1 P(t)
        explVar <- t(beta) %*% gam1 %*% Cov.Z[1:3, 1:3] %*% t(gam1) %*% beta
        totVar <- m$rho1.Z %*% Cov.Z %*% m$rho1.Z
        return(list(beta=beta, alpha=alpha, R2=explVar/totVar))
    }

    ## b) calculate beta using sample moments
    if (span) {
        pars.smpl <- getParsSpanned(cov(m$Z), colMeans(m$Z))
    } else {
        pars.smpl <- getPars(cov(m$Z), colMeans(m$Z))
    }

    ## c) calculate beta using population moments
    cN <- 5
    E.Z.pop <- as.vector(solve(diag(cN) - m$KP.ZZ) %*% m$KP.0Z)
    Cov.Z.pop <- matrix(solve(diag(cN^2) - kronecker(m$KP.ZZ, m$KP.ZZ)) %*% as.numeric(m$Omega.Z), cN, cN)
    if (span) {
        pars.pop <- getParsSpanned(Cov.Z.pop, E.Z.pop)
    } else {
        pars.pop <- getPars(Cov.Z.pop, E.Z.pop)
    }

    tbl <- createTbl(c(paste("OLS", macro.name, sep="-"),
                       paste("OLS-", macro.name, ", t-stat", sep=""),
                       paste("JPS-", macro.name, ", smpl.", sep=""),
                       paste("JPS-", macro.name, ", pop.", sep="")),
                     c("Intercept", "b(cycle)", "b(infl)", "R^2"))
    tbl[1,] <- c(lm.data$coef[1]*100, lm.data$coef[2:3], summary(lm.data)$r.squared)
    tbl[2, 1:3] <- abs(lm.data$coef/sqrt(diag(vcovHAC(lm.data))))
    tbl[3, 1:3] <- c(pars.smpl$alpha*100, pars.smpl$beta)*12
    tbl[3, 4] <- pars.smpl$R2
    tbl[4, 1:3] <- c(pars.pop$alpha*100, pars.pop$beta)*12
    tbl[4, 4] <- pars.pop$R2
    return(tbl)
}

gaussian.loadings <- function(maturities, K0d, K1d, H0d, rho0d, rho1d, timestep=1) {
  ## maturities: M*1
  ## K0d      : N*1
  ## K1d      : N*1
  ## H0d      : N*N
  ## rho0d    : scalar
  ## rho1d    : N*1
  ## timestep : optional argument.
  ##
  ## By : N*M
  ## Ay : 1*M  (faster to not compute with only one output argument)
  ##
  ## r(t)   = rho0d + rho1d'Xt
  ##        = 1 period discount rate
  ## P(t)   =  price of  t-period zero coupon bond
  ##        = EQ0[exp(-r0 - r1 - ... - r(t-1)]
  ##        = exp(A+B'X0)
  ## yields = Ay + By'*X0
  ##   yield is express on a per period basis unless timestep is provided.
  ##   --For example, if the price of a two-year zero is exp(-2*.06)=exp(-24*.005),
  ##   --and we have a monthly model, the function will return Ay+By*X0=.005
  ##   --unless timestep=1/12 is provided in which case it returns Ay+By*X0=.06
  ##
  ## Where under Q:
  ##   X(t+1) - X(t) = K0d + K1d*X(t) + eps(t+1),  cov(eps(t+1)) = H0d
  ##
  ## A1 = -rho0d
  ## B1 = -rho1d
  ## At = A(t-1) + K0d'*B(t-1) .5*B(t-1)'*H0d*B(t-1) - rho0d
  ## Bt = B(t-1) + K1d'*B(t-1) - rho1d
  ##
  ## maturities: 1*M # of periods
  
  M = length(maturities)
  N = length(K0d)
  Atemp = 0
  Btemp = matrix(0,N,1)
  Ay = matrix(NA,1,M)
  By = matrix(NA,N,M)
  
  curr_mat = 1
  K0dp <- t(K0d)
  K1dp <- t(K1d)
  for (i in 1:maturities[M]) {
    Atemp <- Atemp + K0dp%*%Btemp +.5%*%t(Btemp)%*%H0d%*%Btemp - rho0d
    Btemp <- Btemp + K1dp%*%Btemp - rho1d
    
    if (i==maturities[curr_mat]) {
      Ay[1,curr_mat] <- -Atemp/maturities[curr_mat]
      By[,curr_mat] <- -Btemp/maturities[curr_mat]
      curr_mat <- curr_mat + 1
    }
  }
  
  gaussian.loadings <- list(A = Ay/timestep, B = By/timestep)
}
