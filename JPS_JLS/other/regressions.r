### macro-spanning
## How much information in macro variables is spanned?
## Does unspanned macro risk matter in Taylor rule?
## Does unspanned macro risk predict excess returns?

rm(list=ls())
graphics.off()

source("R/data_fns.r")
source("R/reg_fns.r")
source("R/util_fns.r")

library(xtable)
library(lmtest)
library(sandwich)

##data <- constructDataSet(source.yields = "GSW")
##data <- constructDataSet(source.yields = "JPS") # three obs missing
data <- constructDataSet(source.yields = "Le")

## alternative measures of real activity
gap.varnames <- c('ugap', 'outgap', 'capu_mfg', 'confb_coinc', 'confb_lead')
growth.varnames <- c('GRO', 'CFNAI_ma3', 'gdp_logch_ma3', 'ip_logch_ma3', 'empl_logch_ma3', 'ur_ch_ma3', 'confb_coinc_ch_ma3', 'confb_lead_ch_ma3', 'CFNAI_ma12', 'gdp_yoy', 'ip_yoy', 'empl_yoy', 'ur_yoy', 'confb_coinc_yoy', 'confb_lead_yoy')
cycle.varnames <- c(gap.varnames, growth.varnames)
## cycle.varnames <- c('ugap', 'outgap', 'capu_mfg', 'GRO', 'ip_logch', 'empl_logch', 'ur_ch', 'confb_coinc_ch', 'confb_lead_ch')

## alternative measures of inflation
core.infl.varnames <- c('corecpi', 'corepce', 'trim.CPI', 'trim.PCE')
data$inf.PC.core <- makePC(data[core.infl.varnames])
infl.varnames <- c('INF', core.infl.varnames, 'inf.PC.core', 'cpi', 'pce')
print(round(cor(data[infl.varnames]), digi=3))

macro.varnames <- c(cycle.varnames, infl.varnames)
macro.rowlabels <- c("$UGAP$", "$OUTGAP$", "$CAPU$", "$GRO$", "$IP$", "$EMP$", "$\\DeltaUR$", "$\\DeltaCBCI$", "$\\DeltaCBLI$", "$INF$", "$CORECPI$", "$COREPCE$", "$TRCPI$", "$TRPCE$", "$COREINFPC$", "$CPI$", "$PCE$")

policy.varnames <- c('ugap', 'outgap', 'INF', 'corecpi', 'corepce')
policy.rownames <- c("Unemp.~gap", "Output gap", "INF (JPS)", "Core CPI (yoy)", "Core PCE (yoy)")
nonpolicy.varnames <- c('GRO', 'gdp_logch_ma3', 'ip_logch_ma3', 'empl_logch_ma3', 'gdp_yoy')
nonpolicy.rownames <- c("GRO (JPS)", "Real GDP (ma3)", "Real GDP (yoy)", "IP (ma3)", "Payroll (ma3)")

## plotMacro(data, cycl.meas="GRO", infl.meas="INF")
plotGroUgap(data)

print(cor(data[c('ugap', 'GRO', 'PC2')]))
print(cor(data[c('corecpi', 'INF', 'PC1')]))

##################################################
### analysis

## spanning of macro by the yield curve
span.tbl <- analyzeSpanning(data)

## policy rules
pol.tbl <- analyzeTaylorRules(data)

## excess returns on PCs and macro
pred.tbl <- analyzeReturns(data)

## future macro on PCs and macro
macro.tbl <- analyzeMacro(data, h=1)

## ## future yields
## ypred.tbl <- analyzeYieldPred(data)
## ## current yields
## test.tbl <- testKnifeEdge(data)

cat("# Table 4: Monetary policy rules and unspanned macro variation\n")
tbl <- cbind(pol.tbl[,c(1,4)], span.tbl[,c(1,4,5,6)])[c(policy.varnames, nonpolicy.varnames),]
rownames(tbl) <- c(policy.rownames, nonpolicy.rownames)
print(round(tbl, digi=2))
## xtbl <- xtable(tbl, digits=2)
## print(xtbl, include.rownames=TRUE, include.colname=TRUE, only.contents=TRUE,
##       sanitize.text.function=function(x){x},
##       hline.after=c(-1,0, length(policy.varnames), nrow(xtbl), nrow(xtbl)))

cat("# Table 6: Unspanned macro risk and unspanned macro forecasts\n")
tbl <- cbind(pred.tbl[,1:4], macro.tbl[,c(1,4,5)])[c(policy.varnames, nonpolicy.varnames),]
rownames(tbl) <- c(policy.rownames, nonpolicy.rownames)
print(round(tbl, digi=2))
## xtbl <- xtable(tbl, digits=2)
## print(xtbl, include.rownames=TRUE, include.colname=TRUE, only.contents=TRUE,
##       sanitize.text.function=function(x){x},
##       hline.after=c(-1, 0, length(policy.varnames), nrow(xtbl), nrow(xtbl)))
