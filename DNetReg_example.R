##########Example###########
source("DNetReg.R")
library(huge)
library(tensor)
library(rTensor)
library(rARPACK)
library(MASS)
library(gplots)
library(splines)
n <- 50
N <- 50
T <- 100
D <- bs(seq(0,1,length=T),degree=3,knots=seq(0.1,0.9,length=5))
K <- dim(D)[2] 
R <- 5
s <- 0.1
EPS = 1e-3
x <- rnorm(N)
w1 <- rep(1,R)
b1 <- matrix(rnorm(n*R),n,R)
b3 <- matrix(rnorm(K*R),K,R)
tw1 <- apply(b1,2,norm_vec)*apply(b1,2,norm_vec)*apply(b3,2,norm_vec)
tb1 <- apply(b1,2,scale_vec); tb3 <- apply(b3,2,scale_vec)
talpha <- getTensor(tw1,tb1,tb1,tb3)
tbeta <- array(1,c(n,n,K))
slice_vec <- rep(0,n*(n-1)/2)
slice_vec[sample(c(1:(n*(n-1)/2)),ceiling(n*n*s/2))] <- 1
slice <- matrix(0,n,n)
slice[upper.tri(slice)] <- slice_vec
slice <- slice+t(slice)
tbeta <- tbeta*outer(slice,rep(1,K))
Z <- generateTensor(x,D,talpha,tbeta)
#initial_sim <- findinitial(Z,x,D)

######DEdgeReg########
ENB <- EW0(Z,x,D)
Terror1 <- sqrt(sum((ENB[[1]]-talpha)^2))
Berror1 <- sqrt(sum((ENB[[2]]-tbeta)^2))
mse1 <- 0
for(j in 1:N){
  mse1 <- mse1+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-
                           logistic(vec(as.tensor(tensor(ENB[[1]]+x[j]*ENB[[2]],D,3,2)))))^2))
}
mse1 <- mse1/N

#######DNetReg###########
###Rank tuning for the baseline coefficient tensor B_0
Rrange <- c(1:10)
eBIC <- c()
for(iter in 1:10){
  R <- Rrange[iter]
  output <- VCNR0(Z,x,D=D,R=R,initial=ENB[[1]],eta=1/300/300/300/100)
  eBIC[iter] <- output$BIC
  plot(eBIC,main="eBIC for rank selection")
}
R <- which.min(eBIC)
result <- VCNR0(Z,x,D=D,R=R,initial=ENB[[1]],eta=1/300/300/300/100)
######lambda tuning
lambda.result <- lambda.select(Z,x,D,w1=result$w,b1=result$b1,b3=result$b3,
                               beta0=array(0,c(n,n,K)),R=R,lambda.vec=seq(20,120,length=10))
lambda <- lambda.result[[1]]
output1 <- VCNR2(Z,x,D,R=R,w1=result$w,b1=result$b1,b3=result$b3,
                 beta=lambda.result[[2]],lambda=lambda)
Terror <- sqrt(sum((output1$alpha-talpha)^2))
Berror <- sqrt(sum((output1$beta_new_entry-tbeta)^2))
mse <- 0
for(j in 1:N){
  mse <- mse+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-
                         logistic(vec(as.tensor(tensor(output1$alpha+x[j]*output1$beta_new_entry,D,3,2)))))^2))
}
mse <- mse/N
tpr <- sum((output1$beta_new_entry!=0)*(tbeta!=0))/sum(tbeta!=0)
fpr <- (sum(output1$beta!=0)-sum((output1$beta!=0)*(tbeta!=0)))/sum(tbeta==0)


##########EdgeReg###########
entry_re <- EW_BASIC(Z,x)
mse2 <- 0
for(j in 1:N){
  mse2 <- mse2+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-
                           logistic(vec(as.tensor(entry_re[[1]]+x[j]*entry_re[[2]]))))^2))
}
mse2  <-  mse2/N
tpr2  <-  sum((entry_re[[3]]!=0)*(tbeta[,,1]!=0))/(n*n*s)
fpr2  <-  (sum(entry_re[[3]]!=0)-sum((entry_re[[3]]!=0)*(tbeta[,,1]!=0)))/(n*n*(1-s))




