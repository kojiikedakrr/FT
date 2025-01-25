## Table 3: Tests of knife-edge unspanned macro restrictions in reduced-form regressions

rm(list=ls())
library(xtable)
library(lmtest)
library(sandwich)
library(car)
source("R/data_fns.r")

estSystem <- function(Y, X, q=NA) {
    ## arguments
    ##  Y is TxM
    ##  X is TxN
    N <- ncol(X)
    X <- cbind(1, X)  ## T x N+1
    Chat <- crossprod(Y, X) %*% solve(crossprod(X)) ## M x N+1
    ## print(round(Chat, digi=2))
    ehat <- Y - X %*% t(Chat)
    Omega.hat <- 1/nrow(ehat) * crossprod(ehat)
    ##Omega.hat <- diag(diag(Omega.hat))
    ## print(det(Omega.hat))
    if (det(Omega.hat)<0)
        print(eigen(Omega.hat)$values)
    beta.hat <- as.numeric(t(Chat))  ## OLS estimates stacked - first N+1 elements correspond to first equation
    vcov.std <- kronecker(Omega.hat, solve(crossprod(X)))
    rval <- list(Omega.hat=Omega.hat, beta.hat=beta.hat, vcov=vcov.std, resids=ehat)
    ## Newey-West standard errors
    if (!is.na(q)) {
        M <- ncol(Y)
        T <- nrow(Y)
        bread <- kronecker(diag(M), solve(crossprod(X)))
        Gamma.j <- function(j) {
            tmp <- matrix(0, (N+1)*M, (N+1)*M)
            for (t in (j+1):T)
                tmp <- tmp + kronecker(ehat[t,] %o% ehat[t-j,], X[t,] %o% X[t-j,])
            tmp
        }
        weights <- seq(1, 0, by = -(1/(q + 1)))
        meat <- 0.5 * Gamma.j(0)
        for (j in 1:q)
            meat <- meat + (1 - j/(q+1)) * Gamma.j(j)
        meat <- meat + t(meat)
        rval <- c(rval, list(vcov.nw=bread %*% meat %*% bread))
    }
    rval
}

summaryNW <- function(model, q=12) UseMethod("summaryNW")
summaryNW.lm <- function(model, q=12){
    if (!require(sandwich)) stop("Required sandwich package is missing.")
    V <- NeweyWest(model, lag=q, prewhite=FALSE, adjust=FALSE)
    sumry <- summary(model)
    table <- coef(sumry)
    table[,2] <- sqrt(diag(V))
    table[,3] <- table[,1]/table[,2]
    table[,4] <- 2*pt(abs(table[,3]), df.residual(model), lower.tail=FALSE)
    sumry$coefficients <- table
    sumry
}

cat("## likelihood-ration tests\n")
tbl1 <- matrix(NA, 5, 4)
tbl2 <- matrix(NA, 5, 6)
colnames(tbl1) <- c("GRO-INF", "UGAP-CPI", 'c.v. 5%', 'c.v. 0.1%')
colnames(tbl2) <- c("GRO", "INF", "UGAP", "CPI", 'c.v. 5%', 'c.v. 0.1%')
rownames(tbl1) <- 1:5
rownames(tbl2) <- 1:5
col1 <- 1
col2 <- 1
for (flag.jpsmacro in c(TRUE, FALSE)) {
    loadDTSMdata(flag.jpsmacro)
    J <- ncol(Y)
    eig <- eigen(cov(Y))
    W <- t(eig$vectors)  ## rows contain loadings, first row is first PC
    ## W <- diag(J)
    ## W <- W[c(1,5,12, 2:4,6:11),]
    for (N in 1:5) {
        P <- Y %*% t(W[1:N,,drop=FALSE])
        Y2 <- Y %*% t(W[(N+1):J,])  ## can't use Y because otherwise linear combinations of errors are zero
        Om1 <- estSystem(Y2, cbind(P, M.o))$Omega.hat   ## unrestricted
        Om0 <- estSystem(Y2, P)$Omega.hat               ## knife-edge restrictions
        Om0a <- estSystem(Y2, cbind(P, M.o[,2]))$Omega.hat  ## only first macro variable excluded
        Om0b <- estSystem(Y2, cbind(P, M.o[,1]))$Omega.hat  ## only second macro variable excluded
        LR <- (nrow(Y)-1-N-2)*(log(det(Om0)) - log(det(Om1)))    # with Sims small-sample correction
        LRa <- (nrow(Y)-1-N-2)*(log(det(Om0a)) - log(det(Om1)))
        LRb <- (nrow(Y)-1-N-2)*(log(det(Om0b)) - log(det(Om1)))
        tbl1[N, col1] <- LR
        tbl1[N, 3] <- qchisq(0.05, (J-N)*2, lower.tail=FALSE)
        tbl1[N, 4] <- qchisq(0.001, (J-N)*2, lower.tail=FALSE)
        tbl2[N, col2] <- LRa
        tbl2[N, col2+1] <- LRb
        tbl2[N, 5] <- qchisq(.05, J-N, lower.tail=FALSE)
        tbl2[N, 6] <- qchisq(.001, J-N, lower.tail=FALSE)
    }
    col1 <- col1+1
    col2 <- col2+2
}
tbl <- cbind(tbl1, tbl2)
print(round(tbl, digi=1))

