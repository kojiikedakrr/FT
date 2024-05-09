####メモ####
# 単位時間を１週間（水曜日）に設定（論文まま）
# 定常性を持たせると定数項が期間全体の平均値と近しい値となり、平均的に誤差を抑えようとする動きが強くなる。
# データの差分をとるなどの加工せずにそのまま使用している→非定常であることは当たり前であるため菊池の非定常性の仮定は動作しない可能性が高い
# しかし、非定常とするとパラメータ測定がうまくいかない可能性が高まるか？→初期値設定次第で何とかなりそう。
# →パラ推計の結果を見るに、対角成分は1に近い数値になっている。定常性の仮定をおいても問題はなさそう？むしろ推定を安定化させる可能性あり。
# 観測方程式のボラは固定したほうがパラ推定が安定。状態方程式のボラが吸収？菊池の論文からはreparametrizeすることで観測のボラ固定は正当化されるようなことが記載されている。
# 状態変数のプロットからは金利推移の傾向はしっかり推計できている。K_Pなどの係数の大小により結果がうまく出ていないように見える。
# →やはりK_Pの初期設定を調整することで金利を説明することができるようになってきた。パラメータの正負・大きさの考察などを行う必要がありそう。
# また、最適化を行うとぼちぼち異なる結果が算出されるが金利自体は説明できている。いくつかの設定で計測してその結果の平均をとるなどの工夫は必要そう。
# 短期はよく説明できているが、長期が微妙。短期の積み重ねによる制約が影響しているか？平均値がキーになる可能性？つまりは、λの数値の精緻化？観測誤差は関係なし？
# →λのみの最適化も行ったが効果なし。最初の設定による影響が強いか→対角成分だけで推定した場合はタームプレミアムは安定するか要確認。
# 状態変数の初期値がかなり重要な役割を果たす。ばらけるように初期値をもっと離すなどの工夫は必要かも。
# Vとλを一緒に推定することは安定せず難しい。Vを0.0001に固定すると推定は安定化する（ここは菊池に倣っている）。そのあとにVのみの推定をトライしてみてもよいかも。
# 菊池のいうところのreparametrizeについては調査が必要そう。

# 残りの課題：
# KWモデル周りの数理的ロジックの整理
# そのほかの推定法SANNやBGFSなどによるパラメータ推計調査
# タームプレミアムの推計・確認・考察
# 報告書作成
# 将来の課題のチェック・時間があれば実装
# アロケーションのKWモデルとACMモデルの関係性について理数的に調べる

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

#菊池　式(8)　→　服部ペーパーより
update_B = function(rho,K,n){
  ndim = nrow(K)
  B = rep(0,ndim)
  for(i in 1:n){
    B = -rho + t(K) %*% B
  }
  return(B)
}

#菊池　式(9)　→　服部ペーパーより
update_A = function(rho0,rho1,lambda,K,W,n){
  if(n == 1){
    A = -rho0
    B = -rho1
  }else{
    A = -rho0
    B = -rho1
    for(i in 2:n){
      B = update_B(rho1,K,i-1)
      A = A - rho0 + 0.5 * t(B) %*% W %*% t(W) %*% B - t(B) %*% W %*% lambda
    } 
  }
  return(A)
}

#KWモデルのLog Likelihoodsを返す関数
build_dlm = function(params,input_dat,opt_init){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  # V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  # rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  # rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  
  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  # diag(K_P) = (exp(diag(K_P))-exp(-diag(K_P))) / (exp(diag(K_P))+exp(-diag(K_P))) 
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  # lambda_b = diag(params[MRP_b_params_range])
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6

  #観測方程式のパラメータ設計
  # V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  V = diag(0.00005,ndim,ndim)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B[tau_iter,] = update_B(rho,K_Q,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,lambda_a,K_Q,W,tau) /tau * -1
    A = c(A,tmp_A)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = K_P
  W = W
  m0 = opt_init[1:LATENT_DIM]
  C0 = diag(exp(opt_init[(LATENT_DIM+1):(2*LATENT_DIM)]))
  return(dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}

#最適化したパラメータを整理して出力する関数
get_optim_param = function(params){
  #事前準備（パラメータ変数から取得する範囲を決める）
  params_list = list()
  
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  # V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  # rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  # rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  
  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  # diag(K_P) = (exp(diag(K_P))-exp(-diag(K_P))) / (exp(diag(K_P))+exp(-diag(K_P))) 
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  # lambda_b = diag(params[MRP_b_params_range])
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  
  #観測方程式のパラメータ設計
  # V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  V = diag(0.00005,ndim,ndim)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B[tau_iter,] = update_B(rho,K_Q,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,lambda_a,K_Q,W,tau) /tau * -1
    A = c(A,tmp_A)
  }
  
  #出力整理
  params_list["K_P"] = list(K_P)
  params_list["W"] = list(W)
  params_list["lambda_a"] = list(lambda_a)
  params_list["lambda_b"] = list(lambda_b)
  params_list["K_Q"] = list(K_Q)
  params_list["V"] = list(V)
  params_list["rho_zero"] = list(rho_zero)
  params_list["rho"] = list(rho)
  params_list["A"] = list(A)
  params_list["B"] = list(B)
  
  return(params_list)
}

#KWモデルのDLM Filtering
Filtered_dlm = function(params,input_dat){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  K_P_params_range = 1:((LATENT_DIM+1)*LATENT_DIM/2)
  MRP_a_params_range = (length(K_P_params_range)+1):(length(K_P_params_range)+LATENT_DIM)
  MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM*(LATENT_DIM+1)/2)
  # MRP_b_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+LATENT_DIM)
  W_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+LATENT_DIM)
  # V_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+ndim)
  # rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+1)
  # rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(V_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  rho_zero_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+1)
  rho_params_range = (length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+1):(length(K_P_params_range)+length(MRP_a_params_range)+length(MRP_b_params_range)+length(W_params_range)+length(rho_zero_params_range)+LATENT_DIM)
  
  #状態方程式のパラメータ設計
  K_P = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  K_P[lower.tri(K_P,diag = T)] = params[K_P_params_range]
  # diag(K_P) = (exp(diag(K_P))-exp(-diag(K_P))) / (exp(diag(K_P))+exp(-diag(K_P))) 
  W = diag(exp(params[W_params_range]),nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_a = params[MRP_a_params_range]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[MRP_b_params_range]
  # lambda_b = diag(params[MRP_b_params_range])
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  
  #観測方程式のパラメータ設計
  # V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  V = diag(0.00005,ndim,ndim)
  rho_zero = params[rho_zero_params_range]
  rho = params[rho_params_range]
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B[tau_iter,] = update_B(rho,K_Q,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,lambda_a,K_Q,W,tau) /tau * -1
    A = c(A,tmp_A)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = K_P
  W = W
  # m0 = c(-10,0,10)
  m0 = tmp_smooth$s[1,]
  C0 = diag(apply(tmp_smooth$s,2,sd), LATENT_DIM)
  return(dlmFilter(dat_wo_A,mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}

#KWモデルのDLM Filtering
Estimate_TermPremium = function(opt_params, state_variable){

  #推計されたパラメータの整理
  ndim = length(MATURITY)
  lambda_a = opt_params$lambda_a
  lambda_b = opt_params$lambda_b
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K = opt_params$K_P
  W = opt_params$W
  A = opt_params$A
  B = opt_params$B
  
  B_RF = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A_RF = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B_RF[tau_iter,] = update_B(rho,K,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,rep(0,3),K,W,tau) /tau * -1
    A_RF = c(A_RF,tmp_A)
  }
  
  #イールド推計
  estimate_yield = matrix(A, ncol = 7, nrow = nrow(state_variable), byrow = T) + state_variable %*% t(B)
  estimate_RF_yield = matrix(A_RF, ncol = 7, nrow = nrow(state_variable), byrow = T) + state_variable %*% t(B_RF)
  term_premium = estimate_yield - estimate_RF_yield
  
  return(term_premium)
}


############################################################
#--------------------MAIN PROGRAMM---------------------
############################################################

####データ取得・加工####
setwd("C:/Users/kojii/OneDrive/office/term_premium")
dat = read_csv("./data/feds200628.csv")
dat = dat %>% mutate(Date = ymd(Date)) %>% filter(Date>=START_DATE,Date<=END_DATE) #日付のフォーマット変更と抽出

svensson_parameter = dat[,SVENSSON_PARAMS]
svensson_yield_dat = as_tibble(t(apply(svensson_parameter,1,svensson_ZCR,t = MATURITY)))
# svensson_yield_dat = na.omit(svensson_yield_dat) #NAを削除で処理（課題：削除でよいのか？KFはNA込みで推定できないか？）

# 水曜日の週次データ取得
date_index = dat[,"Date"]
svensson_yield_dat = bind_cols(date_index,svensson_yield_dat) %>% mutate(youbi = wday(Date)) %>% filter(youbi == 5)
svensson_yield_dat = svensson_yield_dat %>% select(-c("Date","youbi")) %>% drop_na()

####パラメータ推計####

# パラメータ推計準備
set.seed(19921002)
init_params = c(0.9,0,0,0.9,0,0.9,rep(0,3),rep(0,6),rep(1,3),0,rep(0,3))

# Nelder-Mead
tmp = optim(init_params,build_dlm,input_dat = as.matrix(svensson_yield_dat/100), opt_init = c(-10,0,10,4,4,4))
opt_par = get_optim_param(tmp$par)
tmp_filter = Filtered_dlm(tmp$par,input_dat = as.matrix(svensson_yield_dat/100))
tmp_smooth = dlmSmooth(tmp_filter)
for(i in 1:3){
  tmp = optim(tmp$par,build_dlm, input_dat = as.matrix(svensson_yield_dat/100), opt_init = c(tmp_smooth$s[1,],log(apply(tmp_smooth$s,2,sd))))
  opt_par = get_optim_param(tmp$par)
  tmp_filter = Filtered_dlm(tmp$par,input_dat = as.matrix(svensson_yield_dat/100))
  tmp_smooth = dlmSmooth(tmp_filter)
  print(paste(i,"回目：",tmp$value))
}
# tmp = optim(tmp$par,build_dlm, input_dat = as.matrix(svensson_yield_dat/100), opt_init = c(-10,0,10,4,4,4))
# tmp = optim(tmp$par,build_dlm, input_dat = as.matrix(svensson_yield_dat/100), opt_init = c(-10,0,10,4,4,4),control = list(maxit=5*10^3))
# opt_par = get_optim_param(tmp$par)
# tmp_filter = Filtered_dlm(tmp$par,input_dat = as.matrix(svensson_yield_dat/100))
# tmp_smooth = dlmSmooth(tmp_filter)


filter_dat = tmp_filter$f + matrix(opt_par$A, ncol = 7, nrow = nrow(tmp_filter$f), byrow = T)
filter_dat = filter_dat %>% as_tibble() %>% mutate(SV0.25_filter = SV0.25, SV0.5_filter = SV0.5,SV1_filter = SV1,SV2_filter = SV2,SV4_filter = SV4,SV7_filter = SV7,SV10_filter = SV10)
filter_dat = filter_dat %>% select(contains("filter"))

Term_Premium = Estimate_TermPremium(opt_params = opt_par, state_variable = tmp_smooth$s)

####グラフ描写####
plot_dat = as_tibble(cbind(filter_dat*100, svensson_yield_dat))
maturity = 10
g = ggplot(plot_dat) + geom_line(aes(x = 1:length(unlist(plot_dat[paste0("SV",maturity)])),y=unlist(plot_dat[paste0("SV",maturity)]))) + geom_line(aes(x = 1:length(unlist(plot_dat[paste0("SV",maturity,"_filter")])),y=unlist(plot_dat[paste0("SV",maturity,"_filter")])))
plot(g)

#状態変数のフィルタリング
plot(tmp_filter$m[,1], type="l")
plot(tmp_filter$m[,2], type="l")
plot(tmp_filter$m[,3], type="l")

#状態変数のスムージング
plot(tmp_smooth$s[,1], type="l")
plot(tmp_smooth$s[,2], type="l")
plot(tmp_smooth$s[,3], type="l")

#Term Premium
plot(Term_Premium[,1], type="l")
plot(Term_Premium[,2], type="l")
plot(Term_Premium[,3], type="l")
plot(Term_Premium[,4], type="l")
plot(Term_Premium[,5], type="l")
plot(Term_Premium[,6], type="l")
plot(Term_Premium[,7], type="l")




#############################
#-------没コード
#############################

#KWモデルのLog Likelihoodsを返す関数
build_dlm_lambda = function(params,input_dat,opt_params,opt_init){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K_P = opt_params$K_P
  W = opt_params$W
  
  #状態方程式のパラメータ設計
  lambda_a = params[1:LATENT_DIM]
  lambda_b = matrix(0,nrow=LATENT_DIM,ncol=LATENT_DIM)
  lambda_b[lower.tri(lambda_b,diag = T)] = params[(LATENT_DIM+1):(2*LATENT_DIM)]
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  
  #観測方程式のパラメータ設計
  # V = diag(exp(params[V_params_range]),nrow=ndim,ncol=ndim)
  V = diag(0.0001,ndim,ndim)
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B[tau_iter,] = update_B(rho,K_Q,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,lambda_a,K_Q,W,tau) /tau * -1
    A = c(A,tmp_A)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = K_P
  W = W
  m0 = opt_init[1:LATENT_DIM]
  C0 = diag(exp(opt_init[(LATENT_DIM+1):(2*LATENT_DIM)]))
  return(dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}

#KWモデルのLog Likelihoodsを返す関数
build_dlm_V = function(params,input_dat,opt_params,opt_init){
  #事前準備（パラメータ変数から取得する範囲を決める）
  ndim = length(MATURITY)
  rho_zero = opt_params$rho_zero
  rho = opt_params$rho
  K_P = opt_params$K_P
  W = opt_params$W
  
  #状態方程式のパラメータ設計
  lambda_a = opt_params$lambda_a
  lambda_b = opt_params$lambda_b
  K_Q = K_P - W %*% lambda_b #kim and Orphanides(2004) (5)式 P.6
  
  #観測方程式のパラメータ設計
  V = diag(exp(params),ndim,ndim)
  B = matrix(NA, nrow = ndim, ncol = LATENT_DIM)
  A = c()
  for(tau_iter in 1:ndim){
    tau = MATURITY[tau_iter] * 52
    B[tau_iter,] = update_B(rho,K_Q,tau) / tau * -1
    tmp_A = update_A(rho_zero,rho,lambda_a,K_Q,W,tau) /tau * -1
    A = c(A,tmp_A)
  }
  
  #dlmパッケージ用パラメータ設定
  dat_wo_A = input_dat - matrix(A, ncol = ndim, nrow = nrow(input_dat), byrow = T)
  FF = B
  V = V
  GG = K_P
  W = W
  m0 = opt_init[1:LATENT_DIM]
  C0 = diag(exp(opt_init[(LATENT_DIM+1):(2*LATENT_DIM)]))
  return(dlmLL(y = dat_wo_A, mod = dlm(FF=FF, V=V, GG=GG, W=W, m0=m0, C0=C0)))
}
