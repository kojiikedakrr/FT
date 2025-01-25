analyzeSpanning <- function(data) {
    ## spanning of macro information by yield curve
    cat("## R-squared in spanning regressions: \n")
    tbl <- matrix(NA, length(macro.varnames), 8)
    rownames(tbl) <- macro.varnames
    colnames(tbl) <- c("3 PCs", "5 PCs", "FFR", "level", "slope", "curve", "PC4", "PC5")
    lms <- list()
    lms.5 <- list()
    lms.ffr <- list()
    for (depvar in macro.varnames) {
        lms[[depvar]] <- lm(get(depvar) ~ PC1+PC2+PC3, data=data)
        lms.5[[depvar]] <- lm(get(depvar) ~ PC1+PC2+PC3+PC4+PC5, data=data)
        lms.ffr[[depvar]] <- lm(get(depvar) ~ ff_avg, data=data)
        for (i in 1:5)
            tbl[depvar, 3+i] <- summary(lm(get(depvar) ~ get(paste("PC",i,sep="")), data=data))$r.squared
    }
    tbl[,1] <- sapply(lms, function(lm) summary(lm)$r.squared)
    tbl[,2] <- sapply(lms.5, function(lm) summary(lm)$r.squared)
    tbl[,3] <- sapply(lms.ffr, function(lm) summary(lm)$r.squared)

    return(tbl)
}

analyzeTaylorRules <- function(data) {
    ### R^2 in Taylor rules
    pol.tbl <- matrix(NA, length(macro.varnames), 6)
    rownames(pol.tbl) <- macro.varnames
    colnames(pol.tbl) <- c("R2", "tstat", "marg. R2", "partial R^2", "rel.RMSE", "ratio of SER")
    lm.corecpi <- lm(ff_avg ~ corecpi, data=data)
    lm.ugap <- lm(ff_avg ~ ugap, data=data)
    for (i in macro.varnames) {
        if (i %in% cycle.varnames) {
            j <- "corecpi"
        } else {
            j <- "ugap"
        }
        taylor <- lm(ff_avg ~ get(i) + get(j), data=data)
        lm.XZ <- lm(get(i) ~ get(j), data=data)
        lm.YZ <- lm(ff_avg ~ get(j), data=data)
        pol.tbl[i,1] <- summary(taylor)$r.squared
        pol.tbl[i,2] <- abs(taylor$coef[2]/sqrt(diag(vcovHAC(taylor)))[2])
        pol.tbl[i,3] <- summary(taylor)$r.squared - summary(lm.corecpi)$r.squared
        pol.tbl[i,4] <- cor(lm.XZ$residuals, lm.YZ$residuals)^2
        pol.tbl[i,5] <- sqrt(mean(taylor$resid^2))/sqrt(mean(lm.YZ$resid^2))
        pol.tbl[i,6] <- summary(taylor)$sigma/summary(lm.YZ)$sigma
    }
    cat("*** UNIVARIATE R^2 of CORECPI for FFR:", summary(lm(ff_avg ~ corecpi, data=data))$r.squared, '\n')
    cat("*** UNIVARIATE R^2 of UGAP for FFR:", summary(lm(ff_avg ~ ugap, data=data))$r.squared, '\n')
    return(pol.tbl)
}


analyzeReturns <- function(data) {
    cat("## Predictability of returns\n")
    tbl <- matrix(NA, length(macro.varnames), 8)
    rownames(tbl) <- macro.varnames
    colnames(tbl) <- c("Rets on PCs+M", "t-stat", "pval", "rel.RMSE",
                       "Rets on 5PCs+M", "t-stat", "pval", "rel.RMSE")
    lm.3 <- lm(xr.avg ~ PC1 + PC2 + PC3, data=data)
    lm.5 <- lm(xr.avg ~ PC1 + PC2 + PC3 + PC4 + PC5, data=data)
    cat("3 PCs only: R^2 =", summary(lm.3)$r.squared, "\n")
    cat("5 PCs only: R^2 =", summary(lm.5)$r.squared, "\n")
    for (varname in macro.varnames) {
        lm.3m <- lm(xr.avg ~ PC1 + PC2 + PC3 + get(varname), data=data)
        lm.5m <- lm(xr.avg ~ PC1 + PC2 + PC3 + PC4 + PC5 + get(varname), data=data)
        tbl[varname, 1] <- summary(lm.3m)$r.squared - summary(lm.3)$r.squared
        tstats3 <- abs(lm.3m$coef)/sqrt(diag(vcovNW(lm.3m)))
        tstats5 <- abs(lm.5m$coef)/sqrt(diag(vcovNW(lm.5m)))
        tbl[varname, 2] <- tail(tstats3,1)
        tbl[varname, 3] <- pt(tail(tstats3,1), df=lm.3m$df, lower.tail=FALSE)*2
        tbl[varname, 4] <- sqrt(mean(lm.3m$resid^2))/sqrt(mean(lm.3$resid^2))
        tbl[varname, 5] <- summary(lm.5m)$r.squared - summary(lm.5)$r.squared
        tbl[varname, 6] <- tail(abs(lm.5m$coef)/sqrt(diag(vcovNW(lm.5m))),1)
        tbl[varname, 7] <- pt(tail(tstats5,1), df=lm.5m$df, lower.tail=FALSE)*2
        tbl[varname, 8] <- sqrt(mean(lm.5m$resid^2))/sqrt(mean(lm.5$resid^2))
    }
    tbl
}

analyzeMacro <- function(data, h=1) {
    ## if variable is spanned, then including macro lag does not improve RMSE
    ##  relative RMSE = 1
    ## if variable is not spanned, then including macro lag lowers RMSE a lot
    ##  relative RMSE = rmse(lm.my)/rmse(lm.y)
    tbl <- matrix(NA, length(macro.varnames), 5)
    rownames(tbl) <- macro.varnames
    colnames(tbl) <- c("Autocorr.", "R2(y)", "R2(my)", "t(m==0)", "rel.RMSE")
    lags=12
    for (varname in macro.varnames) {
        lead.varname <- paste('lead', varname, sep='.')
        data[lead.varname] <- c(data[[varname]][(1+h):nrow(data)], rep(NA, h))
        lm.m <- lm(get(lead.varname) ~ get(varname), data=data)
        lm.y <- lm(get(lead.varname) ~ PC1 + PC2 + PC3 + PC4 + PC5, data=data)
        lm.my <- lm(get(lead.varname) ~ PC1 + PC2 + PC3 + PC4 + PC5 + get(varname), data=data)
        ## lm.y <- lm(get(lead.varname) ~ PC1 + PC2 + PC3, data=data)
        ## lm.my <- lm(get(lead.varname) ~ PC1 + PC2 + PC3 + get(varname), data=data)
        tbl[varname, 1] <- acf(data[varname], plot=FALSE)$acf[2] ## summary(lm.m)$r.squared
        tbl[varname, 2] <- summary(lm.y)$r.squared
        tbl[varname, 3] <- summary(lm.my)$r.squared
        ## lags <- floor(bwNeweyWest(lm.my))
        ## cat("#", varname, " - Newey-West lags:", lags, "\n")
        Vhat <- NeweyWest(lm.my, lag=lags, prewhite=FALSE)
        tbl[varname, 4] <- tail(abs(lm.my$coef)/sqrt(diag(Vhat)),1)
        tbl[varname, 5] <- sqrt(mean(lm.my$residuals^2))/sqrt(mean(lm.y$residuals^2))
    }
    return(tbl)
}

plotGroUgap <- function(data) {
    ## plot PC2, GRO, UGAP
    ## (i) check that PC2 here is the same as PC2 in my JLS model
    ## (ii) check that UGAP here is the same as in my JLS model
    loadDTSMdata()

    W <- eigen(cov(Y))$vectors
    PC2.JLS <- Y %*% W[,2]
    cat("# checking correlations of macro and yield data between two different data sources\n")
    print(cor(PC2.JLS, data$PC2)) # perfect
    print(cor(M.o[,2], data$ugap)) # negative, almost perfect
    print(cor(M.o[,1], data$corecpi)) # perfect

    ##print(cbind(M.o[,2], data$ugap))

    plot.pc2 <- ts((data$PC2-mean(data$PC2))/sd(data$PC2), start=c(1985, 1), freq=12)
    plot.ugap <- ts((data$ugap-mean(data$ugap))/sd(data$ugap), start=c(1985, 1), freq=12)
    plot.GRO <- ts((data$GRO-mean(data$GRO))/sd(data$GRO), start=c(1985, 1), freq=12)
    yrange <- range(plot.pc2, plot.ugap, plot.GRO)
    cols <- c("black", "blue", "red")
    ltys <- c(1,2,1)
    lwds <- c(2,2,1)

    dev.new()
    par(mar = c(4,4,2,1)+.1)
    plot(plot.pc2, type="l", lwd=lwds[1], col=cols[1], lty=ltys[1], ylim=yrange, ylab="Percent", xlab="Year", yaxs="i")
    plot.recessions(yrange)
    lines(plot.pc2, col=cols[1], lty=ltys[1], lwd=lwds[1])
    lines(plot.ugap, col=cols[2], lty=ltys[2], lwd=lwds[2])
    lines(plot.GRO, col=cols[3], lty=ltys[3], lwd=lwds[3])
    legend("bottom", c("Slope", "UGAP", "GRO"), lwd=lwds, lty=ltys, col=cols)
}
