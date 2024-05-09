####分析の前提####
# 単位時間を１週間に設定
# 課題：
# MRPの設定に制約をつけるか否か。下三角or対角。
# rhoの設定
# 初期値の設定
# Term Structureの数式チェック（特に符号など）
# 定常性の仮定
# Qを基準として考えるか、Pを基準として考えるか
# →MRPを対角化することでK_PとK_Qが下三角になるので最適化しやすくなるはず。
# 対角化 or Qを基準にすることを検討

# ある程度安定する推定設定
# 分散共分散W、Vはdiag, lambda_bは下三角、methodはNelder-Mead、QもPも定常を強制的に持たせる、データは/100して%ではなく実数

####パッケージインストール####
library(tidyverse)
library(lubridate)
library(dlm)
library(matrixcalc)

####共通で使用する設定####
LATENT_DIM = 3
MATURITY = c(3/12, 6/12, 1, 2, 4, 7, 10)
START_DATE = ymd("1990/7/1")
END_DATE = ymd("2005/7/1")
SVENSSON_PARAMS = c("BETA0","BETA1","BETA2","BETA3","TAU1","TAU2")

####補助関数を定義####
#Svenssonモデルからゼロクーポンレートを算出する関数
svensson_ZCR = function(t, svensson_parameter){
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
#固有値分解を使用した関数。solveを使っているためエラーがよく発生する
matexp = function(X){
  y <- eigen(X)
  mat = y$vectors %*% diag(exp(y$values)) %*% solve(y$vectors)
  return(Re(mat))
}

#exp(X)を計算する関数（Xは行列)の近似値
#solveのエラーが発生するためにテイラー展開した結果を出力する関数
mat_exp_appr = function(X, inverse = T){
  mat = diag(1,nrow = nrow(X)) + X + (1/2) * matrix.power(X,2) + (1/6) * matrix.power(X,3) + (1/24) * matrix.power(X,4) + 
    (1/120) * matrix.power(X,5) + (1/720) * matrix.power(X,6) + (1/5040) * matrix.power(X,7) + (1/factorial(8)) * matrix.power(X,8) + (1/factorial(9)) * matrix.power(X,9) + 
    (1/factorial(10)) * matrix.power(X,10) + (1/factorial(11)) * matrix.power(X,11) + (1/factorial(12)) * matrix.power(X,12) + (1/factorial(13)) * matrix.power(X,13) + (1/factorial(14)) * matrix.power(X,14) +
    (1/factorial(15)) * matrix.power(X,15) + (1/factorial(16)) * matrix.power(X,16) + (1/factorial(17)) * matrix.power(X,17) + (1/factorial(18)) * matrix.power(X,18) + (1/factorial(19)) * matrix.power(X,19) +
    (1/factorial(20)) * matrix.power(X,20) + (1/factorial(21)) * matrix.power(X,21) + (1/factorial(22)) * matrix.power(X,22) + (1/factorial(23)) * matrix.power(X,23) + (1/factorial(24)) * matrix.power(X,24) +
    (1/factorial(25)) * matrix.power(X,25) + (1/factorial(26)) * matrix.power(X,26) + (1/factorial(27)) * matrix.power(X,27) + (1/factorial(28)) * matrix.power(X,28) + (1/factorial(29)) * matrix.power(X,29)
  if(inverse){mat = solve(mat)}
  return(mat)
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
    M1 = -t(solve(K)) %*% (mat_exp_appr(t(K)*tau_iter,T) - I)
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
    M2 = -solve(vec_inverse(solve(kronecker(K,I)+kronecker(I,K)) %*% vec(mat_exp_appr(K*tau_iter,T)%*%Sigma%*%t(Sigma)%*%mat_exp_appr(t(K)*tau_iter,T) - Sigma%*%t(Sigma)),LATENT_DIM,LATENT_DIM))
    M2_list = append(M2_list,list(M2))
  }
  return(M2_list)
}

#最適化したパラメータを整理して出力する関数
get_optim_param = function(params){
  #事前準備（パラメータ変数から取得する範囲を決める）
  params_list = list()

  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM**2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  # rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1)
  # rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  
  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  # diag(K_P) = diag(K_P)^2 / (1 + diag(K_P)^2) #定常性の仮定
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  # lambda_b = diag(params[MRP_b_params_range],nrow=LATENT_DIM,ncol=LATENT_DIM)
  # lambda_b = matrix(params[MRP_b_params_range],nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  # diag(K_Q) = diag(K_Q)^2 / (1 + diag(K_Q)^2) #定常性の仮定
  diag(K_Q) = (exp(diag(K_Q))-exp(-diag(K_Q))) / (exp(diag(K_Q))+exp(-diag(K_Q))) #定常性の仮定
  mu_Q = -solve(K_Q) %*% W %*% lambda_a #kim and Orphanides(2004) (6)式 P.6
  
  #観測方程式のパラメータ設計
  V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  # V = diag(0.0001,nrow=ndim,ncol=ndim)
  M1 = Update_M1(K_Q,MATURITY*52)
  M2 = Update_M2(K_Q,W,MATURITY*52)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  # rho = rep(1,length(rho_params_range))
  
  ##kim and Orphanides(2004) (48),(49)式 P.37 
  A = c()
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * 52
    tmp_B = M1[[iter]] %*% rho / tau
    tmp_A = -(t(K_Q%*%mu_Q)%*%(M1[[iter]]-tau*diag(1,nrow(M1[[iter]])))%*%t(solve(K_Q))%*%rho + 0.5 * t(rho) %*% solve(K_Q) %*% (M2[[iter]] - W%*%t(W)%*%M1[[iter]] - M1[[iter]]%*%W%*%t(W) + tau*W%*%t(W)) %*% t(solve(K_Q)) %*% rho - tau*rho_zero)/tau
    A = c(A,tmp_A)
    B[iter,] = t(tmp_B)
  }
  
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
  
  return(params_list)
}

#KWモデルのLog Likelihoodsを返す関数
build_dlm = function(params,input_dat){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  
  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  # diag(K_Q) = diag(K_Q)^2 / (1 + diag(K_Q)^2) #定常性の仮定
  # diag(K_Q) = (exp(diag(K_Q))-exp(-diag(K_Q))) / (exp(diag(K_Q))+exp(-diag(K_Q))) #定常性の仮定
  mu_Q = -solve(K_Q) %*% W %*% lambda_a #kim and Orphanides(2004) (6)式 P.6

  #観測方程式のパラメータ設計
  V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  # V = diag(0.0001,nrow=ndim,ncol=ndim)
  M1 = Update_M1(K_Q,MATURITY*52)
  M2 = Update_M2(K_Q,W,MATURITY*52)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  # rho = rep(1,length(rho_params_range))
  
  ##kim and Orphanides(2004) (48),(49)式 P.37 
  A = c()
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * 52
    tmp_B = M1[[iter]] %*% rho / tau
    tmp_A = -(t(K_Q%*%mu_Q)%*%(M1[[iter]]-tau*diag(1,nrow(M1[[iter]])))%*%t(solve(K_Q))%*%rho + 0.5 * t(rho) %*% solve(K_Q) %*% (M2[[iter]] - W%*%t(W)%*%M1[[iter]] - M1[[iter]]%*%W%*%t(W) + tau*W%*%t(W)) %*% t(solve(K_Q)) %*% rho - tau*rho_zero)/tau
    A = c(A,tmp_A)
    B[iter,] = t(tmp_B)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = mat_exp_appr(K_P,T)
  # GG = diag(GG)^2 / (1 + diag(GG)^2)
  # GG = (exp(diag(GG))-exp(-diag(GG))) / (exp(diag(GG))+exp(-diag(GG))) #定常性の仮定
  W = -vec_inverse(solve(kronecker(K_P,diag(1,LATENT_DIM))+kronecker(diag(1,LATENT_DIM),K_P))%*%vec(mat_exp_appr(K_P,T)%*%W%*%t(W)%*%mat_exp_appr(t(K_P),T)-W%*%t(W)),LATENT_DIM,LATENT_DIM)
  m0 = c(0,0,0)
  C0 = diag(10.0, LATENT_DIM)
  return(dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}

#KWモデルのDLM Filtering
Filtered_dlm = function(params,input_dat){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)

  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  diag(K_P) = diag(K_P)^2 / (1 + diag(K_P)^2) #定常性の仮定
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  diag(K_Q) = diag(K_Q)^2 / (1 + diag(K_Q)^2) #定常性の仮定
  mu_Q = -solve(K_Q) %*% W %*% lambda_a #kim and Orphanides(2004) (6)式 P.6
  
  #観測方程式のパラメータ設計
  V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  M1 = Update_M1(K_Q,MATURITY*52)
  M2 = Update_M2(K_Q,W,MATURITY*52)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  
  ##kim and Orphanides(2004) (48),(49)式 P.37 
  A = c()
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  for(iter in 1:ndim){
    tau = MATURITY[iter] * 52
    tmp_B = M1[[iter]] %*% rho / tau
    tmp_A = -(t(K_Q%*%mu_Q)%*%(M1[[iter]]-tau*diag(1,nrow(M1[[iter]])))%*%t(solve(K_Q))%*%rho + 0.5 * t(rho) %*% solve(K_Q) %*% (M2[[iter]] - W%*%t(W)%*%M1[[iter]] - M1[[iter]]%*%W%*%t(W) + tau*W%*%t(W)) %*% t(solve(K_Q)) %*% rho - tau*rho_zero)/tau
    A = c(A,tmp_A)
    B[iter,] = t(tmp_B)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = mat_exp_appr(K_P,T)
  W = -vec_inverse(solve(kronecker(K_P,diag(1,LATENT_DIM))+kronecker(diag(1,LATENT_DIM),K_P))%*%vec(mat_exp_appr(K_P,T)%*%W%*%t(W)%*%mat_exp_appr(t(K_P),T)-W%*%t(W)),LATENT_DIM,LATENT_DIM)
  m0 = rep(0.0,LATENT_DIM)
  C0 = diag(10.0, LATENT_DIM)
  return(dlmFilter(dat_wo_A,mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}


#--------------------MAIN PROGRAMM---------------------

####データ取得・加工####
setwd("C:/Users/kojii/OneDrive/office/term_premium")
dat = read_csv("./data/feds200628.csv")
dat = dat %>% mutate(Date = ymd(Date)) %>% filter(Date>=START_DATE,Date<=END_DATE) #日付のフォーマット変更と抽出

svensson_parameter = dat[,SVENSSON_PARAMS]
svensson_yield_dat = as_tibble(t(apply(svensson_parameter,1,svensson_ZCR,t = MATURITY)))
svensson_yield_dat = na.omit(svensson_yield_dat) #NAを削除で処理（課題：削除でよいのか？KFはNA込みで推定できないか？）

####パラメータ推計####
set.seed(19921002)
tmp = optim(rep(1,29),build_dlm,input_dat = as.matrix(svensson_yield_dat/100), method = "Nelder-Mead")
for(i in 1:10){
  tmp = optim(tmp$par,build_dlm, input_dat = as.matrix(svensson_yield_dat/100), method = "Nelder-Mead")
  print(paste(i,"回目：",tmp$value))
}
tmp = optim(tmp$par,build_dlm, input_dat = as.matrix(svensson_yield_dat/100), method = "Nelder-Mead", control = list(maxit=10000))
opt_par = get_optim_param(tmp$par)

V_non_const_retult = tmp2
V_non_const_params = opt_par2
V_non_const_filtered_dat = Filtered_dlm(V_non_const_retult$par,as.matrix(svensson_yield_dat/100))

tmp_filter = dlmFilter(as.matrix(svensson_yield_dat[1:100,]),mod = build_dlm(tmp$par))
