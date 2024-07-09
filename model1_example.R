######################
source('./model1_main.R')
library(parallel)
library(doParallel)
cl<-makeCluster(32)
registerDoParallel(cl)
ll <- Sys.time()
print("flag0")
model1_N100_n50_R5_T50_s0.1<-foreach(i=1:50,.combine='rbind') %dopar% { 
  library(huge)
  library(tensor)
  library(rTensor)
  library(rARPACK)
  library(MASS)
  library(gplots)
  library(splines)
  n<-50
  N<-100
  T<-50
  D<-bs(seq(0,1,length=T),degree=3,knots=seq(0.1,0.9,length=5))
  K<-dim(D)[2] 
  R<-5
  s<-0.1
  EPS = 1e-3
  x<-rnorm(N)
  xt <- matrix(rnorm(N*T),N,T)
  w1<-rep(1,R)
  b1<-matrix(rnorm(n*R),n,R)
  b3<-matrix(rnorm(K*R),K,R)
  tw1<-apply(b1,2,norm_vec)*apply(b1,2,norm_vec)*apply(b3,2,norm_vec)
  tb1<-apply(b1,2,scale_vec); tb3<-apply(b3,2,scale_vec)
  talpha<-getTensor(tw1,tb1,tb1,tb3)
  tbeta1<-array(1,c(n,n,K))
  slice_vec<-rep(0,n*(n-1)/2)
  slice_vec[sample(c(1:(n*(n-1)/2)),ceiling(n*n*s/2))]<-1
  slice<-matrix(0,n,n)
  slice[upper.tri(slice)]<-slice_vec
  slice<-slice+t(slice)
  tbeta1<-tbeta1*outer(slice,rep(1,K))
  
  tbeta2<-array(1,c(n,n,K))
  slice_vec<-rep(0,n*(n-1)/2)
  slice_vec[sample(c(1:(n*(n-1)/2)),ceiling(n*n*s/2))]<-1
  slice<-matrix(0,n,n)
  slice[upper.tri(slice)]<-slice_vec
  slice<-slice+t(slice)
  tbeta2<-tbeta2*outer(slice,rep(1,K))
  
  Z<-generateTensor_ex1(x,xt,D,talpha,tbeta1,tbeta2)
  
  
  ###elementwise_spline
  ENB<-EW0_ex1(Z,x,xt,D)
  
  Terror1<-sqrt(sum((ENB[[1]]-talpha)^2))
  Berror1<-sqrt(sum((ENB[[2]]-tbeta1)^2))
  Berror2<-sqrt(sum((ENB[[3]]-tbeta2)^2))
  mse1<-0
  
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(ENB[[3]], D, 3, 2)[,,t]  
    C[,,,t] <- outer(xt[,t],tensor_t)
  }
  
  eta <- tensor(outer(rep(1,N),ENB[[1]])+outer(x,ENB[[2]]),D,4,2)+C
  
  tC <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(tbeta2, D, 3, 2)[,,t]  
    tC[,,,t] <- outer(xt[,t],tensor_t)
  }
  
  
  teta <- tensor(outer(rep(1,N),talpha)+outer(x,tbeta1),D,4,2)+tC
  
  
  mse1 <- 0
  for(j in 1:N){
    mse1<-mse1+sqrt(sum((logistic(vec(as.tensor(teta[j,,,])))-logistic(vec(as.tensor(eta[j,,,]))))^2))
  }
  mse1<-mse1/N
  Perror1<-mse1
  ######
  
  
  ###proposed
  initial_sim<-findinitial(Z,x,D)
  result<-VCNR0(Z,x,D=D,R=R,initial=ENB[[1]],step=1/300/300/300/100)
  lambda.result<-lambda.select_ex1(Z,x,xt,D,w1=result$w,b1=result$b1,b3=result$b3,beta0=array(0,c(n,n,K)),R=R,lambda.vec=seq(20,120,length=10))
  
  lambda<-lambda.result$lambda
  output1<-VCNR2_ex1(Z,x,xt,D,R=R,w1=result$w,b1=result$b1,b3=result$b3,beta1=lambda.result$beta1,beta2=lambda.result$beta2,lambda=lambda)
  Terror<-sqrt(sum((output1$alpha-talpha)^2))
  Berror<-sqrt(sum((output1$beta1-tbeta1)^2))
  B2error<-sqrt(sum((output1$beta2-tbeta2)^2))
  
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(output1$beta2, D, 3, 2)[,,t]  
    C[,,,t] <- outer(xt[,t],tensor_t)
  }
  
  eta <- tensor(outer(rep(1,N),output1$alpha)+outer(x,output1$beta1),D,4,2)+C
  
  
  mse <- 0
  for(j in 1:N){
    mse<-mse+sqrt(sum((logistic(vec(as.tensor(teta[j,,,])))-logistic(vec(as.tensor(eta[j,,,]))))^2))
  }
  mse<-mse/N
  Perror<-mse
  
  tpr1<-sum((output1$beta1!=0)*(tbeta1!=0))/sum(tbeta1!=0)#
  precision1<-sum((output1$beta1!=0)*(tbeta1!=0))/sum(output1$beta1!=0)
  fpr1<-(sum(output1$beta1!=0)-sum((output1$beta1!=0)*(tbeta1!=0)))/sum(tbeta1==0)
  F1_1<-2*precision1*tpr1/(tpr1+precision1)
  
  
  tpr2<-sum((output1$beta2!=0)*(tbeta2!=0))/sum(tbeta2!=0)#
  precision2<-sum((output1$beta2!=0)*(tbeta2!=0))/sum(output1$beta2!=0)
  fpr2<-(sum(output1$beta2!=0)-sum((output1$beta2!=0)*(tbeta2!=0)))/sum(tbeta2==0)
  F1_2<-2*precision2*tpr2/(tpr2+precision2)
  ######
  entry_re<-EW_BASIC_ex1(Z,x,xt)
  
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    C[,,,t] <- outer(xt[,t],entry_re$B2[,,t])
  }
  eta <-outer(rep(1,N),entry_re$B0)+outer(x,entry_re$B1)+C
  
  
  mse2<-0
  for(j in 1:N){
    mse2<-mse2+sqrt(sum((logistic(vec(as.tensor(teta[j,,,])))-logistic(vec(as.tensor(eta[j,,,]))))^2))
  }
  mse2<-mse2/N
  Perror2<-mse2
  en_tpr1<-sum((entry_re$signal1!=0)*(tbeta1[,,1]!=0))/(n*n*s)
  en_precision1<-sum((entry_re$signal1!=0)*(tbeta1[,,1]!=0))/sum(entry_re$signal1!=0)
  en_fpr1<-(sum(entry_re$signal1!=0)-sum((entry_re$signal1!=0)*(tbeta1[,,1]!=0)))/(n*n*(1-s))
  en_F1_1<-2*en_precision1*en_tpr1/(en_tpr1+en_precision1)
  
  en_tpr2<-sum((entry_re$signal2!=0)*(tbeta2[,,1]!=0))/(n*n*s)
  en_precision2<-sum((entry_re$signal2!=0)*(tbeta2[,,1]!=0))/sum(entry_re$signal2!=0)
  en_fpr2<-(sum(entry_re$signal2!=0)-sum((entry_re$signal2!=0)*(tbeta2[,,1]!=0)))/(n*n*(1-s))
  en_F1_2<-2*en_precision2*en_tpr2/(en_tpr2+en_precision2)
  
  a<-rbind(Terror,Berror,B2error,Perror,tpr1,fpr1,F1_1,tpr2,fpr2,F1_2,Terror1,Berror1,Berror2,Perror1,Perror2,en_tpr1,en_fpr1,en_F1_1,en_tpr2,en_fpr2,en_F1_2)
}
e <- Sys.time()
print(e-ll)
save(model1_N100_n50_R5_T50_s0.1,file="./model1_N100_n50_R5_T50_s0.1.RData")
