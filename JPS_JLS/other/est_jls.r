## estimate JLS model
## - risk factors are first L PCs of yields and M macrofactors
## - PCs are observed
## - macro factors are observed
## - normalization as in JLS

rm(list=ls())
set.seed(616)
source("R/jsz_fns.r")
source("R/jls_fns.r")
source("R/data_fns.r")
source("R/util_fns.r")

library(randtoolbox) # for sobol random sequence

cL <- 3

## estimate JLS model with GRO/INF
cat("## GRO/INF\n")
loadDTSMdata(flag.jpsmacro=TRUE)
pars <- estJLS(Y, W, M.o, mats, cL)
est.llk <- jls.llk.kinfQ(Y, M.o, W, cL, kinfQ=pars$kinfQ, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
cat("LLK, explicit kinfQ =", -sum(est.llk$llk), "\n")
est.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
save(pars, est.llk, Y, W, M.o, cL, mats, dates, n.per,
     file = paste("estimates/jls_gro_L", cL, ".RData", sep=""))

## estimate JLS model with UGAP/CPI
cat("## UGAP/CPI\n")
loadDTSMdata(flag.jpsmacro=FALSE)
pars <- estJLS(Y, W, M.o, mats, cL)
est.llk <- jls.llk.kinfQ(Y, M.o, W, cL, kinfQ=pars$kinfQ, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
cat("LLK, explicit kinfQ =", -sum(est.llk$llk), "\n")
est.llk <- jls.llk.kinfQ(Y, M.o, W, cL, lamQ=pars$lamQ, gam0=pars$gam0, gam1=pars$gam1, Sigma=pars$Sigma, mats=mats, dt=1)
save(pars, est.llk, Y, W, M.o, cL, mats, dates, n.per,
     file = paste("estimates/jls_ugap_L", cL, ".RData", sep=""))

