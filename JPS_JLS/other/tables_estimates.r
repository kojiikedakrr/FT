rm(list=ls())

printTblEstJLS <- function(file.estimates, file.out) {
    ## create Latex table with parameter estimates for JLS model
    load(file.estimates)
    cN <- length(pars$mu)
    cM <- nrow(pars$gam1)
    cL <- cN - cM
    ## for JLS estimation, M.o are first two variables in VAR
    ## - need to reorder to make comparable to JPS/consistent with notation in paper - first two rows become last two rows
    r <- c(3:5, 1:2)
    mu <- pars$mu[r]
    Phi <- pars$Phi[r,r]
    Omega <- pars$Omega[r,r]
    Sigma <- t(chol(Omega))
    cat("\\begin{tabular}{rrrrrr} \\toprule \n")
    cat("\\multicolumn{6}{l}{Spanned model: $SM(3,2)$} \\\\ \\midrule \n")
    cat("\\multicolumn{1}{c}{$k_\\infty^\\mathds{Q}$} & \\multicolumn{1}{c}{$\\lambda^\\mathds{Q}$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){2-6}\n")
    cat(sprintf("%2.3f", c(1200*pars$kinfQ, 1+sort(pars$lamQ, dec=TRUE))), sep="&")
    cat("\\\\ \n")
    cat("\\multicolumn{1}{c}{$\\gamma_0$} & \\multicolumn{1}{c}{$\\gamma_1$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){2-6}\n")
    for (i in 1:cM) {
        cat(sprintf("%2.3f", c(pars$gam0[i], pars$gam1[i,])), sep="&")
        cat("\\\\ \n")
    }
    cat("\\multicolumn{1}{c}{$\\mu$} & \\multicolumn{1}{c}{$\\Phi$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){2-6}\n")
    for (i in 1:cN) {
        cat(sprintf("%2.3f", c(mu[i], Phi[i,])), sep="&")
        cat("\\\\ \n")
    }
    cat("\\multicolumn{1}{c}{$\\Sigma$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-5}\n")
    for (i in 1:cN) {
        for (j in 1:cN) {
            if (j>i) cat("0") else cat(sprintf("%2.3f", Sigma[i,j]))
            if (j<cN) cat("&")
        }
        cat("\\\\ \n")
    }
    cat("\\multicolumn{1}{c}{$\\sigma_e$} && \\multicolumn{2}{c}{LLK} && \\multicolumn{1}{c}{$\\Phi$ ev.}\\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){3-4} \\cmidrule(lr{.25em}){6-6}\n")
    cat(sprintf("%2.3f", pars$sigma.e*1200))
    cat("&& \\multicolumn{2}{c}{", sprintf("%5.1f", pars$llk), "}")
    cat("&&\\multicolumn{1}{c}{", sprintf("%2.3f", max(abs(eigen(pars$Phi)$values))), "}\\\\ \n")
    cat("\\midrule\n")
}

printTblEstJPS <- function(file.estimates) {
    ## create Latex table with parameter estimates for JPS model
    load(file.estimates)
    cat("\\multicolumn{6}{l}{Unspanned model: $USM(3,2)$} \\\\ \\midrule \n")
    cat("\\multicolumn{1}{c}{$k_\\infty^\\mathds{Q}$} & \\multicolumn{1}{c}{$\\lambda^\\mathds{Q}$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){2-4}\n")
    cat(sprintf("%2.3f", c(1200*pars$kinfQ, sort(pars$lamQ, dec=TRUE))), sep="&")
    cat("\\\\ \n")
    cat("\\multicolumn{1}{c}{$\\mu$} & \\multicolumn{1}{c}{$\\Phi$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){2-6}\n")
    for (i in 1:N) {
        cat(sprintf("%2.3f", c(pars$KP.0Z[i], pars$KP.ZZ[i,])), sep="&")
        cat("\\\\ \n")
    }
    cat("\\multicolumn{1}{c}{$\\Sigma$} \\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-5}\n")
    for (i in 1:N) {
        for (j in 1:N) {
            if (j>i) cat("0") else cat(sprintf("%2.3f", pars$L.Z[i,j]))
            if (j<N) cat("&")
        }
        cat("\\\\ \n")
    }
    cat("\\multicolumn{1}{c}{$\\sigma_e$} && \\multicolumn{2}{c}{LLK} && \\multicolumn{1}{c}{$\\Phi$ ev.}\\\\ \n")
    cat("\\cmidrule(lr{.25em}){1-1} \\cmidrule(lr{.25em}){3-4} \\cmidrule(lr{.25em}){6-6}\n")
    cat(sprintf("%2.3f", sqrt(pars$sige2)*1200))
    cat("&& \\multicolumn{2}{c}{", sprintf("%5.1f", -sum(est.llk$llk)), "}")
    cat("&& \\multicolumn{1}{c}{", sprintf("%2.3f", max(abs(eigen(pars$KP.ZZ)$values))), "} \\\\ \n")
    cat("\\bottomrule \n\\end{tabular}\n")
}

cat("# Table BI - Parameter estimates for macro-finance models using data set with GRO/INF\n")
printTblEstJLS("estimates/jls_gro_L3_published.RData")
printTblEstJPS("estimates/jps_gro_R3_published.RData")

cat("# Table BII -  Parameter estimates for macro-finance models using data set with UGAP/CPI\n")
printTblEstJLS("estimates/jls_ugap_L3_published.RData")
printTblEstJPS("estimates/jps_ugap_R3_published.RData")
