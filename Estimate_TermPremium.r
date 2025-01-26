####パッケージインストール####\
setwd("C:/Users/kojii/OneDrive/office/term_premium")
library(tidyverse)
library(lubridate)
source("./JPS_JLS/R/jsz_fns.r")
source("./JPS_JLS/R/jps_fns.r")
source("./JPS_JLS/R/data_fns.r")
source("./JPS_JLS/R/util_fns.r")
source("./JPS_JLS/R/spanning_fns.r")

####共通で使用する設定####
set.seed(19921002)
FREQ = 12
MATURITY = c(3/12, 6/12, 1, 2, 4, 7, 10)
START_DATE = ymd("1990/7/1")
END_DATE = ymd("2024/4/30")
SVENSSON_PARAMS = c("BETA0","BETA1","BETA2","BETA3","TAU1","TAU2")

#Svenssonモデルからゼロクーポンレートを算出する関数
svensson_ZCR = function(svensson_parameter,t){
  # パラメータ設定
  beta0 = svensson_parameter[1]
  beta1 = svensson_parameter[2]
  beta2 = svensson_parameter[3]
  beta3 = svensson_parameter[4]
  tau1 = 1/svensson_parameter[5]
  tau2 = 1/svensson_parameter[6]
  
  # Svenssonのゼロクーポンレート計算
  SV = c()
  for(t_iter in t){
    SV1 = beta0
    SV2 = beta1 * (1-exp(-t_iter*tau1))/(t_iter*tau1)
    SV3 = beta2 * ((1-exp(-t_iter*tau1))/(t_iter*tau1) - exp(-t_iter*tau1))
    SV4 = beta3 * ((1-exp(-t_iter*tau2))/(t_iter*tau2) - exp(-t_iter*tau2))
    SV = c(SV, (SV1+SV2+SV3+SV4))
  }
  
  # 返値
  names(SV) = paste0("SV",t)
  return(unlist(SV))
}


# データの準備
dat = read_csv("./data/feds200628.csv")
dat = dat %>% mutate(Date = ymd(Date)) %>% filter(Date>=START_DATE,Date<=END_DATE) #日付のフォーマット変更と抽出
dat = fill(dat,all_of(SVENSSON_PARAMS))
dat = dat %>% mutate(month_index = floor_date(Date, "month")) %>% group_by(month_index) %>% filter(Date == max(Date)) %>% ungroup()

date_index = dat[,"Date"]
svensson_parameters = dat[,SVENSSON_PARAMS]
svensson_yield_dat = t(apply(svensson_parameters,1,svensson_ZCR,t = MATURITY)) / 100

PCA_Yield = prcomp(svensson_yield_dat,center = TRUE, scale. = TRUE)
PCA_Rotation = t(-PCA_Yield$rotation)


econ_dat = read_csv("./data/econ_data.csv")
econ_dat = econ_dat %>% mutate(Date = ymd(Date) - 1) %>% drop_na() %>% 
  mutate(US_CPI = log(US_CPI), Industrial_Production = log(Industrial_Production))

econ_dat = econ_dat %>% filter(Date>=START_DATE,Date<=END_DATE) #日付のフォーマット変更と抽出
econ_dat = as.matrix(econ_dat[,-1])

## estimate JSZ model
N <- 3
mats = MATURITY * FREQ

Y = svensson_yield_dat
W = PCA_Rotation

pars <- estJSZ(Y, W, N, mats)
est.llk <- jsz.llk(Y, W[1:N,], K1Q.X=diag(pars$lamQ-1), Sigma.cP = pars$Omega.cP,
                   mats=mats, dt=1)
save(data, pars, est.llk, Y, W, mats,
     file=paste("./JPS_JLS/estimates/jsz_N", N, ".RData", sep=""))

tmp = loadJSZmodel(name = "TEST",filename = paste("./JPS_JLS/estimates/jsz_N", N, ".RData", sep=""))
JSZ_TP = tmp$Ytp
colnames(JSZ_TP) = paste0("JSZ_TP_",MATURITY)

## estimate unspanned-macro-risk (UMR) models
R <- 3
N <- R+2
M.o = econ_dat

### estimate JPS model with GRO/INF
# cat("## GRO/INF\n")
# loadDTSMdata(flag.jpsmacro=TRUE)
pars <- estJPS(Y, M.o, W, R, mats, gamma=rep(1, R*(N+1)))
cat("optimal LLK:", pars$llk, "\n")
est.llk <- jps.llk(Y, M.o, W, R, KQ.XX=diag(pars$lamQ), Omega.Z = pars$Omega.Z, mats=mats, dt=1)
cat("check       ", -sum(est.llk$llk), "\n")
save(pars, est.llk, Y, W, M.o, R, N, mats,
     file=paste("./JPS_JLS/estimates/jps_gro_R", R, ".RData", sep=""))

tmp = loadJPSmodel("TEST",paste("./JPS_JLS/estimates/jps_gro_R", R, ".RData", sep=""))
JPS_TP = tmp$Ytp
colnames(JPS_TP) = paste0("JPS_TP_",MATURITY)


## plot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

compare_dat = tibble(cbind(date_index,Term_Premium,JSZ_TP,JPS_TP))

TP_plot2 = ggplot(compare_dat) + 
  geom_line(aes(x=Date,y=SV2_Premium * 100, col = "KW"),linewidth = 2) + 
  geom_line(aes(x=Date,y=JSZ_TP_2 * 100, col = "JSZ"),linewidth = 2) +
  geom_line(aes(x=Date,y=JPS_TP_2 * 100, col = "JPS"),linewidth = 2) +
  theme_bw() +
  scale_y_continuous(breaks=seq(-2,10,2),limits=c(-2,4),labels = str_c(seq(-2,10,2),'%')) +
  labs(x='Date', y='Term Premium 2Y', col='Model')

TP_plot10 = ggplot(compare_dat) + 
  geom_line(aes(x=Date,y=SV10_Premium * 100, col = "KW"),linewidth = 2) + 
  geom_line(aes(x=Date,y=JSZ_TP_10 * 100, col = "JSZ"),linewidth = 2) +
  geom_line(aes(x=Date,y=JPS_TP_10 * 100, col = "JPS"),linewidth = 2) +
  theme_bw() +
  scale_y_continuous(breaks=seq(-2,10,2),limits=c(-2,8),labels = str_c(seq(-2,10,2),'%')) +
  labs(x='Date', y='Term Premium 10Y', col='Model')

multiplot(TP_plot2,TP_plot10)
