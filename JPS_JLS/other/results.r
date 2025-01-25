rm(list=ls())
source("R/jls_fns.r")
source("R/jps_fns.r")
source("R/jsz_fns.r")
source("R/data_fns.r")
source("R/util_fns.r")
source("R/spanning_fns.r")

library(xtable)

models <- list(
    SM.GRO = loadJLSmodel("SM.GRO", "estimates/jls_gro_L3_published.RData"),
    SM.UGAP = loadJLSmodel("SM.UGAP", "estimates/jls_ugap_L3_published.RData"),
    USM.GRO = loadJPSmodel("USM.GRO", "estimates/jps_gro_R3_published.RData"),
    USM.UGAP = loadJPSmodel("USM.UGAP", "estimates/jps_ugap_R3_published.RData"),
    JSZ = loadJSZmodel("JSZ3", "estimates/jsz_N3_published.RData")
)

cat("# Table 1 - cross-sectional fit:\n")
tbl <- yieldFit(models)[,c(1, 2,3,4,5,6,8,10,13,14)]
print(round(tbl, d=1))
## xtbl <- xtable(tbl, digits=1)
## print(xtbl, sanitize.text.function=function(x){x}, only.contents=T, include.colnames=T, include.rownames=T)

cat("# error serial correlation:\n")
tbl <- errorSerialCorr(models)
print(round(tbl, d=4))

cat("# Table 2 - test knife-edge restrictions\n")
testModelBased(models)

## Tables 5, 7, 8 - small-sample simulations
T <- nrow(models$USM.GRO$Y)
M <- 10000
tbl.gro.usm <- analyzeSpanning(models$USM.GRO, T, M)
tbl.gro.sm3 <- analyzeSpanning(models$SM.GRO, T, M, R=3)
tbl.gro.sm5 <- analyzeSpanning(models$SM.GRO, T, M)
tbl.ugap.usm <- analyzeSpanning(models$USM.UGAP, T, M)
tbl.ugap.sm3 <- analyzeSpanning(models$SM.UGAP, T, M, R=3)
tbl.ugap.sm5 <- analyzeSpanning(models$SM.UGAP, T, M)
## tables with all results: data, simulations SM, simulations USM
tbl.gro <- rbind(tbl.gro.sm3[1,], tbl.gro.sm5[1,], tbl.gro.usm[-1,], tbl.gro.sm3[-1,], tbl.gro.sm5[-1,])
rownames(tbl.gro)[1:2] <- c("Data, 3 PCs", "Data, 5 PCs")
tbl.ugap <- rbind(tbl.ugap.sm3[1,], tbl.ugap.sm5[1,], tbl.ugap.usm[-1,], tbl.ugap.sm3[-1,], tbl.ugap.sm5[-1,])
rownames(tbl.ugap)[1:2] <- c("Data, 3 PCs", "Data, 5 PCs")

cat("# Table 5: unspanned macro variation\n")
tbl <- cbind(tbl.gro[,1:2], tbl.ugap[,1:2])
print(round(tbl, digi=2))

cat("# Table 7: unspanned macro risk\n")
tbl <- cbind(tbl.gro[,c(3,7)], tbl.ugap[,c(3,7)])
print(round(tbl, digi=2))

cat("# Table 8: unspanned macro forecasts\n")
tbl <- cbind(tbl.gro[,9:10], tbl.ugap[,9:10])
print(round(tbl, digi=2))

## large-sample results -- for Appendix
cat("## Using model-implied population moments\n")
tbls.an <- list(gro.usm = calcSpanning(models$USM.GRO),
                ugap.usm = calcSpanning(models$USM.UGAP),
                gro.sm5 = calcSpanning(models$SM.GRO),
                ugap.sm5 = calcSpanning(models$SM.UGAP),
                gro.sm3 = calcSpanning(models$SM.GRO, yieldFactorsOnly=TRUE),
                ugap.sm3 = calcSpanning(models$SM.UGAP, yieldFactorsOnly=TRUE))
cat("# Table C1: Unspanned macro variation in MTSMs - population moments results\n")
makePopTbl(tbls.an, i=1)
cat("# Table C2. Unspanned macro risk in MTSMs—population moments\n")
makePopTbl(tbls.an, i=2)
cat("# Table C3. Unspanned macro forecasts in MTSMs—population moments\n")
makePopTbl(tbls.an, i=3)

## ## CHECK
## cat("## Simulation of long sample\n")
## ## note: need to use fixed W - then results are consistent
## T <- 100000
## M <- 2
## tbls.sim <- list(gro.usm = analyzeSpanning(models$USM.GRO, T, M)[c(2,4,6),-c(3,4,7)],
##                  gro.sm5 = analyzeSpanning(models$SM.GRO, T, M)[c(2,4,6),-c(3,4,7)],
##                  gro.sm3 = analyzeSpanning(models$SM.GRO, T, M, R=3)[c(2,4,6),-c(3,4,7)],
##                  ugap.usm = analyzeSpanning(models$USM.UGAP, T, M)[c(2,4,6),-c(3,4,7)],
##                  ugap.sm5 = analyzeSpanning(models$SM.UGAP, T, M)[c(2,4,6),-c(3,4,7)],
##                  ugap.sm3 = analyzeSpanning(models$SM.UGAP, T, M, R=3)[c(2,4,6),-c(3,4,7)])
## cat("UMV - simulation\n")
## makePopTbl(tbls.sim, i=1)
## cat("UMR - simulation\n")
## makePopTbl(tbls.sim, i=2)
## cat("UMF - simulation\n")
## makePopTbl(tbls.sim, i=3)

## term premia
## Figure 2
analyzeFTP(models$SM.GRO, models$USM.GRO, models$JSZ, flag.tofile=FALSE)
## Figure 3
analyzeFTP(models$SM.UGAP, models$USM.UGAP, models$JSZ, flag.tofile=FALSE)


