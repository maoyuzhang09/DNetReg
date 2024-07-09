######################
source('./model2_main.R')
library(parallel)
library(doParallel)
cl<-makeCluster(32)
registerDoParallel(cl)
ll <- Sys.time()
print("flag0")
model2_N100_n50_T50_s005_R2<-foreach(i=1:50,.combine='rbind') %dopar% {
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
  R<-2
  s<-0.05
  EPS = 1e-3
  x<-rnorm(N)
  w1<-rep(1,R)
  b1<-matrix(rnorm(n*R),n,R)
  b3<-matrix(rnorm(K*R),K,R)
  tw1<-apply(b1,2,norm_vec)*apply(b1,2,norm_vec)*apply(b3,2,norm_vec)
  tb1<-apply(b1,2,scale_vec); tb3<-apply(b3,2,scale_vec)
  talpha<-getTensor(tw1,tb1,tb1,tb3)
  
  b_w1<-rep(1,R)
  b_b1<-matrix(rnorm(n*R),n,R)
  b_b1[sample(1:n,ceiling((1-sqrt(s))*n)),]=0
  b_b3<-matrix(rnorm(K*R),K,R)
  b_tw1<-apply(b_b1,2,norm_vec)*apply(b_b1,2,norm_vec)*apply(b_b3,2,norm_vec)
  b_tb1<-apply(b_b1,2,scale_vec); b_tb3<-apply(b_b3,2,scale_vec)
  tbeta<-getTensor(b_tw1,b_tb1,b_tb1,b_tb3)
  tbeta[tbeta!=0]<-1
  
  
  Z<-generateTensor(x,D,talpha,tbeta)
  
  
  ENB<-EW0(Z,x,D)
  Terror1<-sqrt(sum((ENB[[1]]-talpha)^2))
  Berror1<-sqrt(sum((ENB[[2]]-tbeta)^2))
  mse1<-0
  for(j in 1:N){
    mse1<-mse1+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-logistic(vec(as.tensor(tensor(ENB[[1]]+x[j]*ENB[[2]],D,3,2)))))^2))
  }
  mse1<-mse1/N
  Perror1<-mse1
  
  #####
  
  result<-VCNR0_ex2(Z,x,D=D,R=R,initial_alpha=ENB[[1]],initial_beta = ENB[[2]],eta=1/300/300/300/100)
  lambda.result<-lambda.select_ex2(Z,x,D,w1=result$w,b1=result$b1,b3=result$b3,b_w1=result$b_w1,b_b1=result$b_b1,b_b3=result$b_b3,R=R,lambda.vec=seq(2000,4000,length=50))
  lambda<-lambda.result$lambda
  
  output1<-VCNR1_ex2(Z,x,D,R=R,w1=result$w,b1=result$b1,b3=result$b3,b_w1=lambda.result$b_w1,b_b1=lambda.result$b_b1,b_b3=lambda.result$b_b3,lambda=lambda)
  Terror<-sqrt(sum((output1$alpha-talpha)^2))
  Berror<-sqrt(sum((output1$beta-tbeta)^2))
  mse<-0
  for(j in 1:N){
    mse<-mse+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-logistic(vec(as.tensor(tensor(output1$alpha+x[j]*output1$beta,D,3,2)))))^2))
  }
  mse<-mse/N
  Perror<-mse
  tpr<-sum((output1$beta!=0)*(tbeta!=0))/sum(tbeta!=0)
  precision<-sum((output1$beta!=0)*(tbeta!=0))/sum(output1$beta!=0)
  fpr<-(sum(output1$beta!=0)-sum((output1$beta!=0)*(tbeta!=0)))/sum(tbeta==0)
  F1<-2*precision*tpr/(tpr+precision)
  
  ########
  entry_re<-EW_BASIC(Z,x)
  mse2<-0
  for(j in 1:N){
    mse2<-mse2+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-logistic(vec(as.tensor(entry_re[[1]]+x[j]*entry_re[[2]]))))^2))
  }
  mse2<-mse2/N
  Perror2<-mse2
  tpr2<-sum((entry_re[[5]]!=0)*(tbeta[,,1]!=0))/(n*n*s)
  precision2<-sum((entry_re[[5]]!=0)*(tbeta[,,1]!=0))/sum(entry_re[[5]]!=0)
  fpr2<-(sum(entry_re[[5]]!=0)-sum((entry_re[[5]]!=0)*(tbeta[,,1]!=0)))/(n*n*(1-s))
  F1_2<-2*precision2*tpr2/(tpr2+precision2)
  a<-rbind(Terror,Berror,Perror,precision,tpr,fpr,F1,Terror1,Berror1,Perror1,Perror2,precision2,tpr2,fpr2,F1_2)
}
e <- Sys.time()
print(e-ll)
save(model2_N100_n50_T50_s005_R2,file="./model2_N100_n50_T50_s005_R2.RData")

