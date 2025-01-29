####パッケージインストール####
library(tidyverse)
library(lubridate)
library(dlm)
library(matrixcalc)
library(MASS)
library(expm)
library(DEoptim)

####共通で使用する設定####
LATENT_DIM = 3
FREQ = 12
MATURITY = c(3/12, 6/12, 1, 2, 4, 7, 10)
START_DATE = ymd("1990/7/1")
END_DATE = ymd("2024/12/31")
SVENSSON_PARAMS = c("BETA0","BETA1","BETA2","BETA3","TAU1","TAU2")


####補助関数を定義####
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

#exp(X)を計算する関数（Xは行列)
matexp = function(X){
  mat = expm(X)
  return(mat)
}

# modify_K = function(K){
#   values = eigen(K)$values
#   if(is.complex(values)){
#     Re_term = Re(values)
#     Im_term = Im(values)
#     Re_term = ifelse(Re_term>0,Re_term,0.001)
#     values = complex(Re_term,Im_term)
#   }else{
#     values = ifelse(values>0,values,0.001)
#   }
#   vecs = eigen(K)$vectors
#   mod_K = vecs %*% diag(values) %*% ginv(vecs)
#   return(mod_K)
# }

modify_K = function(K){
  values = eigen(K)$values
  values = ifelse(values>0,values,0.001)
  vecs = eigen(K)$vectors
  mod_K = vecs %*% diag(values) %*% ginv(vecs)
  return(mod_K)
}

#vec作用素でベクトル化されたデータを行列に戻す作用素
#参考：https://math.stackexchange.com/questions/805565/what-is-the-inverse-of-the-mboxvec-operator
vec_inverse = function(vectorized_matrix,row_dim,col_dim){
  I_m = diag(1, row_dim)
  I_n = diag(1, col_dim)
  matrix_from_vec = kronecker(t(vec(I_n)),I_m) %*% kronecker(I_n,vectorized_matrix)
  return(matrix_from_vec)
}

#ATSMの係数更新用関数
##kim and Orphanides(2004) (50)式 P.37 
Update_M1 = function(K,tau){
  I = diag(1,nrow(K))
  M1_list = list()
  for(tau_iter in tau){
    M1 = -t(ginv(K)) %*% (matexp(-t(K)*tau_iter) - I)
    M1_list = append(M1_list,list(M1))
  }
  return(M1_list)
}

#ATSMの係数更新用関数
##kim and Orphanides(2004) (51)式 P.37 
Update_M2 = function(K,Sigma,tau){
  I = diag(1,nrow(K))
  M2_list = list()
  for(tau_iter in tau){
    M2 = -vec_inverse(ginv(kronecker(K,I)+kronecker(I,K)) %*% vec(matexp(-K*tau_iter)%*%Sigma%*%t(Sigma)%*%matexp(-t(K)*tau_iter) - Sigma%*%t(Sigma)),LATENT_DIM,LATENT_DIM)
    M2_list = append(M2_list,list(M2))
  }
  return(M2_list)
}

#パラメータを整理して出力する関数
get_param_list = function(params){
  #事前準備（パラメータ変数から取得する範囲を決める）
  params_list = list()

  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  # K_P_params_range = 1:LATENT_DIM
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  # V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim-3)
  V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
 
  #状態方程式のパラメータ設計
  # K_P = diag(params[K_P_params_range])
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  K_P = modify_K(K_P)
  
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  W_sqrt = sqrt(W)
  
  lambda_a = params[MRP_a_params_range]
  # lambda_b = diag(params[MRP_b_params_range])
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  K_Q = K_P - W_sqrt %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  K_Q = modify_K(K_Q)
  mu_Q = -ginv(K_Q) %*% W_sqrt %*% lambda_a #kim and Orphanides(2004) (6)式 P.6
  
  #観測方程式のパラメータ設計
  # V_params = c(params[V_params_range[1]],params[V_params_range[2]],-30,params[V_params_range[3]],-30,params[V_params_range[4]],-30)
  # V = diag(exp(V_params))
  V = diag(exp(params[V_params_range]))
  M1 = Update_M1(K_Q,MATURITY*FREQ)
  M2 = Update_M2(K_Q,W_sqrt,MATURITY*FREQ)
  rho_zero = params[rho_zero_params_range]
  rho = rep(1,LATENT_DIM)
  
  ##kim and Orphanides(2004) (48),(49)式 P.37 
  A = c()
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * FREQ
    tmp_B = M1[[iter]] %*% rho / tau
    tmp_A = -(t(K_Q%*%mu_Q)%*%(M1[[iter]]-tau*diag(1,nrow(M1[[iter]])))%*%t(ginv(K_Q))%*%rho + 0.5 * t(rho) %*% ginv(K_Q) %*% (M2[[iter]] - W_sqrt%*%t(W_sqrt)%*%M1[[iter]] - M1[[iter]]%*%W_sqrt%*%t(W_sqrt) + tau*W_sqrt%*%t(W_sqrt)) %*% t(ginv(K_Q)) %*% rho - tau*rho_zero)/tau
    A = c(A,tmp_A)
    B[iter,] = t(tmp_B)
  }
  
  #init params(推定後のスムージング結果から算出)
  m0 = c(1.0103401,0.2787958,-1.3209884)
  C0 = diag(exp(c(-7.626178,-6.201878,-5.773171)), LATENT_DIM)
  
  #出力整理
  params_list["K_P"] = list(K_P)
  params_list["W"] = list(W)
  params_list["lambda_a"] = list(lambda_a)
  params_list["lambda_b"] = list(lambda_b)
  params_list["K_Q"] = list(K_Q)
  params_list["mu_Q"] = list(mu_Q)
  params_list["V"] = list(V)
  params_list["rho_zero"] = list(rho_zero)
  params_list["rho"] = list(rho)
  params_list["A"] = list(A)
  params_list["B"] = list(B)
  params_list["m0"] = list(m0)
  params_list["C0"] = list(C0)
  
  return(params_list)
}

#KWモデルのLog Likelihoodsを返す関数
build_dlm = function(params){
  param_list = get_param_list(params)
  ndim = length(MATURITY)
  
  #dlmパッケージ用パラメータ設定
  input_dat = as.matrix(svensson_yield_dat/100)
  dat_wo_A = input_dat - matrix(param_list$A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = param_list$B
  V = param_list$V
  GG = matexp(-param_list$K_P)
  W = param_list$W
  m0 = param_list$m0
  C0 = param_list$C0
  LL = dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0))
  LL = ifelse(is.nan(LL),999999,LL)
  return(LL)
}

# #KWモデルのLog Likelihoodsを返す関数
# build_dlm_forDE = function(params){
#   param_list = get_param_list(params)
#   ndim = length(MATURITY)
#   
#   #dlmパッケージ用パラメータ設定
#   input_dat = as.matrix(svensson_yield_dat/100)
#   dat_wo_A = input_dat - matrix(param_list$A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
#   FF = param_list$B
#   V = param_list$V
#   GG = matexp(-param_list$K_P)
#   W = param_list$W
#   m0 = param_list$m0
#   C0 = param_list$C0
#   LL = dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0))
#   LL = ifelse(is.nan(LL),999999,LL)
#   return(LL)
# }

#KWモデルのDLM Filtering
Filtered_dlm = function(params,input_dat){
  param_list = get_param_list(params)
  ndim = length(MATURITY)
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(param_list$A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = param_list$B
  V = param_list$V
  GG = matexp(-param_list$K_P)
  W = param_list$W
  m0 = param_list$m0
  C0 = param_list$C0
  return(dlmFilter(dat_wo_A,mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}

Estimate_TermPremium = function(opt_params, input_dat, state_variable){
  
  ndim = length(MATURITY)
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K = opt_params$K_P
  
  short_rate_mat = matrix(NA,nrow=nrow(state_variable)-1,ncol=ndim)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * FREQ
    expect_x = state_variable[-1,] %*% matexp(K*tau)
    short_rate = rho_zero + expect_x %*% rho
    short_rate_mat[,iter] = short_rate
  }
  
  term_premium = input_dat - short_rate_mat
  return(term_premium)
  
}

Estimate_TermPremium2 = function(opt_params, input_dat, state_variable){
  
  ndim = length(MATURITY)
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K = opt_params$K_Q
  mu_Q = opt_params$mu_Q
  W_sqrt = opt_params$W %>% sqrt()
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  
  M1 = Update_M1(K,MATURITY*FREQ)
  M2 = Update_M2(K,W_sqrt,MATURITY*FREQ)
  
  
  ##kim and Orphanides(2004) (48),(49)式 P.37 
  A = c()
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * FREQ
    tmp_B = M1[[iter]] %*% rho / tau
    tmp_A = -(t(K%*%mu_Q)%*%(M1[[iter]]-tau*diag(1,nrow(M1[[iter]])))%*%t(ginv(K))%*%rho + 0.5 * t(rho) %*% ginv(K) %*% (M2[[iter]] - W_sqrt%*%t(W_sqrt)%*%M1[[iter]] - M1[[iter]]%*%W_sqrt%*%t(W_sqrt) + tau*W_sqrt%*%t(W_sqrt)) %*% t(ginv(K)) %*% rho - tau*rho_zero)/tau
    A = c(A,tmp_A)
    B[iter,] = t(tmp_B)
  }
  
  term_premium = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T) -  state_variable[-1,] %*% t(B)
  return(term_premium)
  
}

Estimate_TermPremium3 = function(opt_params, input_dat, state_variable){
  
  ndim = length(MATURITY)
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K = opt_params$K_P
  
  short_rate_mat = matrix(NA,nrow=nrow(state_variable)-1,ncol=ndim)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * FREQ
    expect_x = state_variable[-1,] %*% matexp(K)
    short_rate = rho_zero + expect_x %*% rho
    for(tau_iter in seq(2,tau)){
      expect_x = state_variable[-1,] %*% matexp(K*tau_iter)
      short_rate = short_rate + rho_zero + expect_x %*% rho
    }
    short_rate_mat[,iter] = short_rate / tau
  }
  
  term_premium = input_dat - short_rate_mat
  return(term_premium)
  
}

#--------------------MAIN PROGRAMM---------------------

####データ取得・加工####
setwd("C:/Users/kojii/OneDrive/office/term_premium")
dat = read_csv("./data/feds200628.csv")
dat = dat %>% mutate(Date = ymd(Date)) %>% filter(Date>=START_DATE,Date<=END_DATE) #日付のフォーマット変更と抽出
dat = fill(dat,all_of(SVENSSON_PARAMS))
dat = dat %>% mutate(month_index = floor_date(Date, "month")) %>% group_by(month_index) %>% filter(Date == max(Date)) %>% ungroup()

date_index = dat[,"Date"]
svensson_parameters = dat[,SVENSSON_PARAMS]
svensson_yield_dat = as_tibble(t(apply(svensson_parameters,1,svensson_ZCR,t = MATURITY)))


####パラメータ推計####
set.seed(19921002)
# lower_bounds = c(rep(0,LATENT_DIM),rep(-2,LATENT_DIM),rep(-2,LATENT_DIM),rep(-20,LATENT_DIM),rep(-20,length(MATURITY)-3),0)
lower_bounds = c(rep(-3,LATENT_DIM*(LATENT_DIM+1)/2),rep(-3,LATENT_DIM),rep(-3,LATENT_DIM*(LATENT_DIM+1)/2),rep(-20,LATENT_DIM),c(-30,-30,-30,-30,-30,-30,-30),-0.3)

upper_bounds = c(rep(3,LATENT_DIM*(LATENT_DIM+1)/2),rep(3,LATENT_DIM),rep(3,LATENT_DIM*(LATENT_DIM+1)/2),rep(-5,LATENT_DIM),c(-30,-30,-30,-30,-30,-30,-30),0.3)
# upper_bounds = c(rep(0.5,LATENT_DIM),rep(2,LATENT_DIM),rep(2,LATENT_DIM),rep(-7,LATENT_DIM),rep(-7,length(MATURITY)-3),1)

optim_params = DEoptim(build_dlm, lower = lower_bounds, upper = upper_bounds,control=list(NP = 20*50, itermax = 300, F = 1, CR = 1))
opt_par = get_param_list(optim_params$optim$bestmem)
# 
tmp = c(-0.926475   , 1.108033    ,0.799890  , -2.826338 ,  -2.474364 ,  -1.076704   , 0.000607  , -0.965778   , 0.489243  , -2.285978   , 1.798633 ,  -2.609103 ,  -1.650448   , 0.179481,   -1.930279  , -8.828258   ,
        -8.703104   ,-7.551265 , -11.145819 , -11.669333 , -30.000000 , -14.719069  ,-30.000000 , -12.316741  ,-30.000000   , 0.021872)
opt_par = get_param_list(tmp)
filtered_dat = Filtered_dlm(tmp,as.matrix(svensson_yield_dat/100))
smoothed_dat = dlmSmooth(filtered_dat)

#DEで算出したパラメータ群
# init_params = c(0.010525, 0.046761, 0.009284, -1.070896, 1.376122, -0.637748, 1.542830, 0.163931,, 
#                 -0.717650, -8.409622, -6.454401, -8.426557, -6.844543, -9.668102, -5.383698, -8.267461, 0.649877)
# 
# optim_params = optim(init_params,build_dlm, lower = lower_bounds, upper = upper_bounds, method = "L-BFGS-B", control = list(trace = 3,maxit = 1000))
# opt_par = get_param_list(optim_params$par)

filtered_dat = Filtered_dlm(optim_params$optim$bestmem,as.matrix(svensson_yield_dat/100))
smoothed_dat = dlmSmooth(filtered_dat)

filter_dat = as_tibble(filtered_dat$f + matrix(opt_par$A, ncol = length(MATURITY), nrow = nrow(filtered_dat$f), byrow = T))
names(filter_dat) = paste0(names(filter_dat),"_filter")
Term_Premium = Estimate_TermPremium(opt_params = opt_par, state_variable = smoothed_dat$s, input_dat = as.matrix(svensson_yield_dat/100))
Term_Premium2 = Estimate_TermPremium2(opt_params = opt_par, state_variable = smoothed_dat$s, input_dat = as.matrix(svensson_yield_dat/100))
Term_Premium3 = Estimate_TermPremium3(opt_params = opt_par, state_variable = smoothed_dat$s, input_dat = as.matrix(svensson_yield_dat/100))
colnames(Term_Premium) = paste0(colnames(Term_Premium), '_Premium')

plot_dat = as_tibble(cbind(date_index, filter_dat*100, svensson_yield_dat, Term_Premium*100))

g_spot2 = ggplot(plot_dat) +
  geom_line(aes(x = Date, y = SV2), linewidth=1, colour = 'black') +
  geom_line(aes(x = Date, y = SV2_filter), linewidth=1, colour = 'red') +
  geom_line(aes(x = Date, y = SV2_Premium), linewidth=1, colour = 'blue') +
  geom_line(aes(x = Date, y = JPS_TP[,4]*100), linewidth=1, colour = 'green') +
  labs(x = 'Date', y = 'spot rate 2Y')+
  scale_y_continuous(breaks=seq(-2,10,2),limits=c(-2,8),labels = str_c(seq(-2,10,2),'%'))
plot(g_spot2)

g_spot10 = ggplot(plot_dat) +
  geom_line(aes(x = Date, y = SV10), linewidth=1, colour = 'black') +
  geom_line(aes(x = Date, y = SV10_filter), linewidth=1, colour = 'red') +
  geom_line(aes(x = Date, y = SV10_Premium), linewidth=1, colour = 'blue') +
  geom_line(aes(x = Date, y = JPS_TP[,7]*100), linewidth=1, colour = 'green') +
  labs(x = 'Date', y = 'spot rate 10Y')+
  scale_y_continuous(breaks=seq(-2,10,2),limits=c(-2,10),labels = str_c(seq(-2,10,2),'%'))
plot(g_spot10)

