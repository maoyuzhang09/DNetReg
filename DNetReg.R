######################################
######## binary networks #############
######################################
library(huge)
library(tensor)
library(rTensor)
library(rARPACK)
library(MASS)
library(gplots)
library(splines)
EPS = 1e-3




norm_vec <- function(x) sqrt(sum(x^2))
scale_vec <- function(x) x/sqrt(sum(x^2))

logit <- function(x){
  log(x/(1-x))
}

logistic <- function(x){
  1-1/(1+exp(x))
}

###get tensor from the CP decomposition
getTensor <- function(w,A,B,C){
  R <- length(w)
  n <- dim(A)[1]
  K <- dim(C)[1]
  T <- array(0,c(n,n,K))
  for (r in 1:R){
    T <- T+outer((A[,r]%o%B[,r]),C[,r])*w[r]
  }
  T
}

loglikelihood <- function(Z,x,D,alpha,beta){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))/N
}


####find the initial B0 (This is another method to find the initial B0, not the Initialization introdunced in the paper)
findinitial <- function(Z,x,D){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  K <- dim(D)[2]
  ZD <- logit(apply(Z,c(2,3,4),sum)/N)
  ZD[ZD>4] <- 4; ZD[ZD<(-4)] <- (-4)
  alpha <- array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      alpha[j,jj,] <- lm(ZD[j,jj,]~D-1)$coefficients
      alpha[jj,j,] <- alpha[j,jj,]
    }
  }
  alpha
}

#####the gradient of u_1#######
b1_gr <- function(b1,D,w1,b3,diff){
  R <- length(w1)
  n <- dim(b1)[1]
  gb1 <- matrix(0,n,length(w1))
  for (r in 1:length(w1)){
    gb1[,r] <- 2*apply(diff*(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r]))),1,sum)
  }
  -gb1
  
}

#####the gradient of u_1#######
b1_max <- function(b1,Z,x,D,w1,b3,beta){
  fn <- 0
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  R <- length(w1)
  K <- dim(b3)[1]
  
  b1 <- matrix(b1,ncol=R)
  alpha <- array(0,c(n,n,K))
  for (r in 1:R){
    alpha <- alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  
  fn <- -sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff <- apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)),c(2,3,4),sum)
  attr(fn,"gradient") <- b1_gr(b1,D,w1,b3,diff)
  fn
}

#####the gradient of u_3
b3_gr <- function(b3,K,T,D,w1,b1,diff){
  R <- length(w1)
  gb3 <- matrix(0,K,length(w1))
  for (r in 1:length(w1)){
    gb3[,r] <- apply(diff*outer(b1[,r]%o%(b1[,r]*w1[r]),rep(1,T)),3,sum)%*%D
  }
  -gb3
}

#####the gradient of u_3######
b3_max <- function(b3,K,T,Z,x,D,w1,b1,beta){
  fn <- 0
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  R <- length(w1)
  
  b3 <- matrix(b3,ncol=R)
  alpha <- array(0,c(n,n,K))
  for (r in 1:R){
    alpha <- alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  fn <- -sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff <- apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)),c(2,3,4),sum)
  attr(fn,"gradient") <- b3_gr(b3,K,T,D,w1,b1,diff)
  fn
  
}

#####the shrinkage operator for the penalty
glasso <- function(beta,lambda){
  K <- dim(beta)[3]
  m <- (1-lambda/apply(beta,c(1,2),norm_vec))
  ind <- m*(m>0)
  beta*outer(ind,rep(1,K))
}



#### Iteration of baseline tensor B_0#######
VCNR0 <- function(Z,x,D,R,initial,eta=1/300/300/10){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  K <- dim(D)[2]
  ## initialize
  Niter <- 20; outeriter <- 0;COND <- FALSE
  while(!COND){
    outeriter <- outeriter+1
    output <- cp(as.tensor(initial),R,max_iter = 50)
    w1 <- output[[1]];b1 <- output$U[[1]];b2 <- output$U[[2]];b3 <- output$U[[3]]
    b3 <- sweep(b3,2,sign(b1[1,])*sign(b2[1,]),"*")
    beta <- array(0,c(n,n,K))
    
    niter <- 20; iter <- 0; cond <- FALSE
    while(!cond){
      iter <- iter+1
      alpha_old <- getTensor(w1,b1,b1,b3)
      b1_temp <- matrix(0,n,R); b3_temp <- matrix(0,K,R)
      diff <- apply(Z-logistic(tensor(outer(rep(1,N),alpha_old)+outer(x,beta),D,4,2)),c(2,3,4),sum)
      for (r in 1:R){
        gb1 <- 2*apply(diff*(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r]))),1,sum)
        gb3 <- apply(diff*outer(b1[,r]%o%(b1[,r]*w1[r]),rep(1,T)),3,sum)%*%D
        b1_temp[,r] <- b1[,r]+gb1*eta
        b3_temp[,r] <- b3[,r]+gb3*eta
      }
      w1 <- w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
      b1 <- apply(b1_temp,2,scale_vec)
      b3 <- apply(b3_temp,2,scale_vec)
      alpha <- getTensor(w1,b1,b1,b3)
      cond <- (iter > niter) | sum((alpha-alpha_old)^2)> 10^10
      print(c(iter,sum((alpha-alpha_old)^2)))
    }
    COND <- (outeriter > Niter)| sum((alpha-alpha_old)^2)< 10^10
  }
  
  ## optimize
  niter <- 20; iter <- 0; cond <- FALSE
  while(!cond){
    iter <- iter+1
    alpha_old <- getTensor(w1,b1,b1,b3)
    optim_result <- nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp <- matrix(optim_result$estimate,ncol=R)
    optim_result <- nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp <- matrix(optim_result$estimate,ncol=R)
    w1 <- w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1 <- apply(b1_temp,2,scale_vec)
    b3 <- apply(b3_temp,2,scale_vec)
    alpha <- getTensor(w1,b1,b1,b3)
    ascent <- sum((alpha-alpha_old)^2)
    cond <- (iter > niter) | ascent<=EPS
    print(c(iter,ascent))
  }
  BIC <- -N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K))*(R*(n+K))
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha))  
}






####the DNetReg with FISTA method######
VCNR2 <- function(Z,x,D,R,w1,b1,b3,beta,lambda){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  K <- dim(D)[2]
  p <- 1
  alpha <- getTensor(w1,b1,b1,b3)
  
  ## optimize
  Niter <- 40; outeriter <- 0;COND <- FALSE
  mu_t <- 4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]
  
  while(!COND){
    outeriter <- outeriter+1
    alpha_prev <- alpha
    beta_prev <- beta
    optim_result <- nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp <- matrix(optim_result$estimate,ncol=R)
    optim_result <- nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp <- matrix(optim_result$estimate,ncol=R)
    w1 <- w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1 <- apply(b1_temp,2,scale_vec)
    b3 <- apply(b3_temp,2,scale_vec)
    alpha <- getTensor(w1,b1,b1,b3)
    
    niter <- 30; iter <- 0; cond <- FALSE;h=1
    while(!cond){
      h_old <- h
      iter <- iter+1
      beta_old <- beta
      diffx <- (Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta_old),D,4,2)))*outer(x,rep(1,n)%o%rep(1,n)%o%rep(1,T))
      gbeta <- apply(tensor(diffx,D,4,1),c(2,3,4),sum)
      beta <- glasso(beta_old+gbeta*mu_t,mu_t*lambda)   
      h=(1+sqrt(1+4*h_old^2))/2
      beta=beta-(h_old-1)/h*(beta-beta_old)
      cond <- (iter > niter) #| ascent<=EPS
    }
    
    ascent <- loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha_prev,beta_prev)
    COND <- (outeriter > Niter) | ascent<=EPS
    print(c(outeriter,ascent))
  }
  beta[which(is.na(beta))]=0
  beta_new_entry <- beta_old_entry <- array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      if(sum(beta[j,jj,])!=0){
        beta_new_entry[j,jj,] <- rep(0,K)
        beta_old_entry[j,jj,] <- beta_new_entry[j,jj,]
        beta_new_result <- nlm(b0_max_e,beta_old_entry[j,jj,],Z=Z[,j,jj,],x=x,D=D,alpha=alpha[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        beta_new_entry[j,jj,] <- beta_new_result$estimate
        beta_new_entry[jj,j,] <- beta_new_entry[j,jj,]
      }
    }
  }
  
  BIC <- -N*loglikelihood(Z,x,D,alpha,beta_new_entry)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta)>0)/2)
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha,beta=beta,beta_new_entry=beta_new_entry))  
}





#######Parameter tuning of lambda
lambda.select <- function(Z,x,D,w1,b1,b3,beta0,R,lambda.vec){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  K <- dim(D)[2]
  p=1
  alpha <- getTensor(w1,b1,b1,b3)  
  mu_t <- 4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]
  m <- length(lambda.vec)
  beta.list <- list()
  eBIC.vec <- c()
  for(i in 1:m){
    beta <- beta0
    lambda <- lambda.vec[i]
    niter <- 30; iter <- 0;cond <- FALSE
    while(!cond){
      iter <- iter+1
      beta_old <- beta
      diffx <- (Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)))*outer(x,rep(1,n)%o%rep(1,n)%o%rep(1,T))
      gbeta <- apply(tensor(diffx,D,4,1),c(2,3,4),sum)
      beta <- glasso(beta+gbeta*mu_t,mu_t*lambda) 
      ascent <- loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha,beta_old)
      cond <- (iter > niter) | (abs(ascent)<=EPS)
    }
    beta.list[[i]] <- beta
    eBIC.vec[i] <- -N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta)>0)/2)
  }
  plot(lambda.vec, eBIC.vec,xlab="sparsity", main="eBIC for sparsity selection")
  return(list(lambda.vec[which.min(eBIC.vec)],beta.list[[which.min(eBIC.vec)]],beta.list))
}

####generate the dynamic networks
generateTensor <- function(x,D,alpha,beta){
  N <- length(x)
  n <- dim(alpha)[1]
  T <- dim(D)[1]
  l <- n*n*T
  
  Z <- array(0,c(N,n,n,T))
  for(i in 1:N){
    Zi <- new("Tensor",3L,c(n,n,T),data=rbinom(n=l,size=1,prob=logistic(vec(as.tensor(tensor(alpha+x[i]*beta,D,3,2))))))
    for(t in 1:T){
      slice <- matrix(0,n,n)
      slice[upper.tri(slice)] <- Zi[,,t]@data[upper.tri(Zi[,,t]@data)]
      slice <- slice+t(slice) 
      diag(slice) <- diag(Zi[,,t]@data)
      Zi[,,t]@data <- slice
    }
    Z[i,,,] <- Zi@data
  }
  Z
}



###DEdgeReg####3
EW0 <- function(Z,x,D){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  K <- dim(D)[2]
  alpha_new_entry <- beta_new_entry <- alpha_old_entry <- beta_old_entry <- array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      alpha_new_entry[j,jj,] <- rep(0,K)
      beta_new_entry[j,jj,] <- rep(0,K)
      niter <- 40; iter <- 0; cond <- FALSE
      while (!cond) {
        iter <- iter+1
        alpha_old_entry[j,jj,] <- alpha_new_entry[j,jj,]
        beta_old_entry[j,jj,] <- beta_new_entry[j,jj,]
        alpha_new_result <- nlm(b0_max_e,alpha_old_entry[j,jj,],Z=Z[,j,jj,],x=x,D=D,beta=beta_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        alpha_new_entry[j,jj,] <- alpha_new_result$estimate
        alpha_new_entry[jj,j,] <- alpha_new_entry[j,jj,]
        beta_new_result <- nlm(b0_max_e,beta_old_entry[j,jj,],Z=Z[,j,jj,],x=x,D=D,alpha=alpha_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        beta_new_entry[j,jj,] <- beta_new_result$estimate
        beta_new_entry[jj,j,] <- beta_new_entry[j,jj,]
        ascent <- loglikelihood_entry(Z[,j,jj,],x,D,alpha_new_entry[j,jj,],beta_new_entry[j,jj,])-loglikelihood_entry(Z[,j,jj,],x,D,alpha_old_entry[j,jj,],beta_old_entry[j,jj,])
        cond <- (iter > niter) | ascent<=EPS
        #print(c(iter,ascent))
      }
    }
  }
  return(list(alpha=alpha_new_entry,beta=beta_new_entry))
}


#####EdgeReg#####
EW_BASIC <- function(Z,x){
  N <- dim(Z)[1]
  n <- dim(Z)[2]
  T <- dim(Z)[4]
  B0 <- B1 <- B1_pvalue <- array(0,c(n,n,T))
  signal <- matrix(0,n,n)
  B1_pvalue_adjust <- array(0,c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      for (t in 1:T){
        fit <- glm(Z[,j,jj,t]~x,family=binomial(link = logit))
        B0[j,jj,t] <- fit$coefficients[1]
        B0[jj,j,t] <- B0[j,jj,t]
        B1[j,jj,t] <- fit$coefficients[2]
        B1[jj,j,t] <- B1[j,jj,t]
        B1_pvalue[j,jj,t] <- summary(fit)$coefficient[2,4]
        B1_pvalue[jj,j,t] <-  B1_pvalue[j,jj,t]
      }
    }
  }
  B1_pvalue_adjust <- array(p.adjust(B1_pvalue,"bonferron"),c(n,n,T))
  #B1_pvalue_adjust <- array(p.adjust(B1_pvalue,"fdr"),c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      if(min(B1_pvalue_adjust[j,jj,])<=0.05){
        signal[j,jj] <- 1
        signal[jj,j] <- signal[j,jj]
      }
    }
  }
  return(list(B0=B0,B1=B1,signal=signal))
}




##########Example###########
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
  plot(eBIC,main="BIC for rank selection")
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
fpr <- (sum(output1$beta!=0)-sum((output1$beta!=0)*(tbeta!=0)))/sum(tbeta!=0)


##########EdgeReg###########
entry_re <- EW_BASIC(Z,x)
mse2 <- 0
for(j in 1:N){
  mse2 <- mse2+sqrt(sum((logistic(vec(as.tensor(tensor(talpha+x[j]*tbeta,D,3,2))))-
                         logistic(vec(as.tensor(entry_re[[1]]+x[j]*entry_re[[2]]))))^2))
}
mse2  <-  mse2/N
tpr2  <-  sum((entry_re[[3]]!=0)*(tbeta[,,1]!=0))/(n*n*s)
fpr2  <-  (sum(entry_re[[3]]!=0)-sum((entry_re[[3]]!=0)*(tbeta[,,1]!=0)))/(n*n*s)




