## estimate JSZ model
rm(list=ls())
set.seed(616)
setwd("C:/Users/kojii/OneDrive/office/term_premium/JPS_JLS")
source("R/jsz_fns.r")
source("R/data_fns.r")
source("R/util_fns.r")

loadDTSMdata()

N <- 3

pars <- estJSZ(Y, W, N, mats)
est.llk <- jsz.llk(Y, W[1:N,], K1Q.X=diag(pars$lamQ-1), Sigma.cP = pars$Omega.cP,
                   mats=mats, dt=1)
save(data, pars, est.llk, Y, W, mats,
     file=paste("estimates/jsz_N", N, ".RData", sep=""))


