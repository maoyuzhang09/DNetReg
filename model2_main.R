
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

logit<-function(x){
  log(x/(1-x))
}

logistic<-function(x){
  1-1/(1+exp(x))
}


getTensor<-function(w,A,B,C){
  R<-length(w)
  n<-dim(A)[1]
  K<-dim(C)[1]
  T<-array(0,c(n,n,K))
  for (r in 1:R){
    T<-T+outer((A[,r]%o%B[,r]),C[,r])*w[r]
  }
  T
}


generateTensor<-function(x,D,alpha,beta){
  N<-length(x)
  n<-dim(alpha)[1]
  T<-dim(D)[1]
  l<-n*n*T
  
  Z<-array(0,c(N,n,n,T))
  for(i in 1:N){
    Zi<-new("Tensor",3L,c(n,n,T),data=rbinom(n=l,size=1,prob=logistic(vec(as.tensor(tensor(alpha+x[i]*beta,D,3,2))))))
    for(t in 1:T){
      slice<-matrix(0,n,n)
      slice[upper.tri(slice)]<-Zi[,,t]@data[upper.tri(Zi[,,t]@data)]
      slice<-slice+t(slice) 
      diag(slice)<-diag(Zi[,,t]@data)
      Zi[,,t]@data<-slice
    }
    Z[i,,,]<-Zi@data
  }
  Z
}




loglikelihood<-function(Z,x,D,alpha,beta){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))/N
}


###cp decomposition of alpha
b1_gr<-function(b1,D,w1,b3,diff){
  R<-length(w1)
  n<-dim(b1)[1]
  gb1<-matrix(0,n,length(w1))
  for (r in 1:length(w1)){
    gb1[,r]<-2*apply(diff*(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r]))),1,sum)
  }
  -gb1
  
}


b1_max<-function(b1,Z,x,D,w1,b3,beta){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(w1)
  K<-dim(b3)[1]
  
  b1<-matrix(b1,ncol=R)
  alpha<-array(0,c(n,n,K))
  for (r in 1:R){
    alpha<-alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  
  fn<--sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff<-apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)),c(2,3,4),sum)
  attr(fn,"gradient")<-b1_gr(b1,D,w1,b3,diff)
  fn
}


b3_gr<-function(b3,K,T,D,w1,b1,diff){
  R<-length(w1)
  gb3<-matrix(0,K,length(w1))
  for (r in 1:length(w1)){
    gb3[,r]<-apply(diff*outer(b1[,r]%o%(b1[,r]*w1[r]),rep(1,T)),3,sum)%*%D
  }
  -gb3
}


b3_max<-function(b3,K,T,Z,x,D,w1,b1,beta){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(w1)
  
  b3<-matrix(b3,ncol=R)
  alpha<-array(0,c(n,n,K))
  for (r in 1:R){
    alpha<-alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  fn<--sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff<-apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)),c(2,3,4),sum)
  attr(fn,"gradient")<-b3_gr(b3,K,T,D,w1,b1,diff)
  fn
  
}









######the CP decomposition of beta
b_b1_gr_ex2<-function(b_b1,x,D,b_w1,b_b3,diff){
  R<-length(b_w1)
  n<-dim(b_b1)[1]
  gb1<-matrix(0,n,length(b_w1))
  
  for (r in 1:length(b_w1)){
    gb1[,r]<-2*apply(diff*outer(x,(rep(1,n)%o%((b_b1[,r]*b_w1[r])%o%as.vector(D%*%b_b3[,r])))),2,sum)
  }
  -gb1
  
}

b_b1_max_ex2<-function(b_b1,Z,x,D,b_w1,b_b3,alpha){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(b_w1)
  K<-dim(b_b3)[1]
  
  b_b1<-matrix(b_b1,ncol=R)
  beta<-array(0,c(n,n,K))
  for (r in 1:R){
    beta<-beta+outer((b_b1[,r]%o%b_b1[,r]),b_b3[,r])*b_w1[r]
  }
  
  fn<--sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff<-(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)))
  attr(fn,"gradient")<-b_b1_gr_ex2(b_b1,x,D,b_w1,b_b3,diff)
  fn
}


b_b3_gr_ex2<-function(b_b3,K,T,x,D,b_w1,b_b1,diff){
  R<-length(b_w1)
  gb3<-matrix(0,K,length(b_w1))
  
  for (r in 1:length(b_w1)){
    gb3[,r]<-apply(diff*outer(x,outer(b_b1[,r]%o%(b_b1[,r]*b_w1[r]),rep(1,T))),4,sum)%*%D
  }
  -gb3
}


b_b3_max_ex2<-function(b_b3,K,T,Z,x,D,b_w1,b_b1,alpha){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(b_w1)
  
  b_b3<-matrix(b_b3,ncol=R)
  beta<-array(0,c(n,n,K))
  for (r in 1:R){
    beta<-beta+outer((b_b1[,r]%o%b_b1[,r]),b_b3[,r])*b_w1[r]
  }
  fn<--sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))
  diff<-(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)))
  attr(fn,"gradient")<-b_b3_gr_ex2(b_b3,K,T,x,D,b_w1,b_b1,diff)
  fn
  
}



######
lasso<-function(b_b1,lambda){
  K<-dim(beta)[3]
  m<-(1-lambda/abs(b_b1))
  ind<-m*(m>0)
  b_b1*ind
}


#######
VCNR0_ex2<-function(Z,x,D,R,initial_alpha,initial_beta,eta=1/300/300/10){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  ## initialize
  Niter<-20; outeriter<-0;COND<-FALSE
  while(!COND){
    outeriter<-outeriter+1
    output<-cp(as.tensor(initial_alpha),R,max_iter = 50)
    w1<-output[[1]];b1<-output$U[[1]];b2<-output$U[[2]];b3<-output$U[[3]]
    b3<-sweep(b3,2,sign(b1[1,])*sign(b2[1,]),"*")
    
    output_beta<-cp(as.tensor(initial_beta),R,max_iter = 50)
    b_w1<-output_beta[[1]];b_b1<-output_beta$U[[1]];b_b2<-output_beta$U[[2]];b_b3<-output_beta$U[[3]]
    b_b3<-sweep(b_b3,2,sign(b_b1[1,])*sign(b_b2[1,]),"*")
    
    #beta<-array(0,c(n,n,K))
    
    niter<-20; iter<-0; cond<-FALSE
    while(!cond){
      iter<-iter+1
      alpha_old<-getTensor(w1,b1,b1,b3)
      beta_old <-getTensor(b_w1,b_b1,b_b1,b_b3)
      b1_temp<-matrix(0,n,R); b3_temp<-matrix(0,K,R)
      b_b1_temp<-matrix(0,n,R); b_b3_temp<-matrix(0,K,R)
      #optim_result<-nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
      #b1_temp<-matrix(optim_result$estimate,ncol=R)
      #optim_result<-nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
      #b3_temp<-matrix(optim_result$estimate,ncol=R)
      diff<-(Z-logistic(tensor(outer(rep(1,N),alpha_old)+outer(x,beta_old),D,4,2)))
      for (r in 1:R){
        gb1<-2*apply(diff*outer(rep(1,N),(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r])))),2,sum)
        gb3<-apply(diff*outer(rep(1,N),outer(b1[,r]%o%(b1[,r]*w1[r]),rep(1,T))),4,sum)%*%D
        b1_temp[,r]<-b1[,r]+gb1*eta
        b3_temp[,r]<-b3[,r]+gb3*eta
        
        b_gb1<-2*apply(diff*outer(x,(rep(1,n)%o%((b_b1[,r]*b_w1[r])%o%as.vector(D%*%b_b3[,r])))),2,sum)
        b_gb3<-apply(diff*outer(x,outer(b_b1[,r]%o%(b_b1[,r]*b_w1[r]),rep(1,T))),4,sum)%*%D
        b_b1_temp[,r]<-b_b1[,r]+b_gb1*eta
        b_b3_temp[,r]<-b_b3[,r]+b_gb3*eta
      }
      w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
      b1<-apply(b1_temp,2,scale_vec)
      b3<-apply(b3_temp,2,scale_vec)
      alpha<-getTensor(w1,b1,b1,b3)
      
      b_w1<-b_w1*apply(b_b1_temp,2,norm_vec)*apply(b_b1_temp,2,norm_vec)*apply(b_b3_temp,2,norm_vec)
      b_b1<-apply(b_b1_temp,2,scale_vec)
      b_b3<-apply(b_b3_temp,2,scale_vec)
      beta<-getTensor(b_w1,b_b1,b_b1,b_b3)
      
      cond<-(iter > niter) | sum((alpha-alpha_old)^2)> 10^10
      print(c(iter,sum((alpha-alpha_old)^2)))
    }
    COND<-(outeriter > Niter)| sum((alpha-alpha_old)^2)< 10^10
  }
  
  ## optimize
  niter<-20; iter<-0; cond<-FALSE
  while(!cond){
    iter<-iter+1
    alpha_old<-getTensor(w1,b1,b1,b3)
    beta_old <-getTensor(b_w1,b_b1,b_b1,b_b3)
    optim_result<-nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta_old,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp<-matrix(optim_result$estimate,ncol=R)
    optim_result<-nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta_old,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp<-matrix(optim_result$estimate,ncol=R)
    w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1<-apply(b1_temp,2,scale_vec)
    b3<-apply(b3_temp,2,scale_vec)
    alpha<-getTensor(w1,b1,b1,b3)
    
    optim_result<-nlm(b_b1_max_ex2,b_b1,Z=Z,x=x,D=D,b_w1=b_w1,b_b3=b_b3,alpha=alpha,check.analyticals=FALSE,gradtol=1e-3)
    b_b1_temp<-matrix(optim_result$estimate,ncol=R)
    optim_result<-nlm(b_b3_max_ex2,b_b3,K=K,T=T,Z=Z,x=x,D=D,b_w1=b_w1,b_b1=b_b1,alpha=alpha,check.analyticals=FALSE,gradtol=1e-3)
    b_b3_temp<-matrix(optim_result$estimate,ncol=R)
    
    b_w1<-b_w1*apply(b_b1_temp,2,norm_vec)*apply(b_b1_temp,2,norm_vec)*apply(b_b3_temp,2,norm_vec)
    b_b1<-apply(b_b1_temp,2,scale_vec)
    b_b3<-apply(b_b3_temp,2,scale_vec)
    beta<-getTensor(b_w1,b_b1,b_b1,b_b3)
    
    
    #ascent<-loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha_old,beta_old)
    ascent <- sum((alpha-alpha_old)^2)#loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha_old,beta)#sum((alpha-alpha_old)^2)#
    cond<-(iter > niter) | ascent<=EPS
    print(c(iter,ascent))
  }
  
  BIC<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K))*(2*R*(n+K))
  # BIC1<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N)+log(n*n*K))*(R*(n+K))
  #loglike<--N*loglikelihood(Z,x,D,alpha,beta)
  #penalty<-(log(N*n*n*T)+log(n*n*K))*(R*(n+K))
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha,b_w1=b_w1,b_b1=b_b1,b_b3=b_b3,beta=beta))  
}



VCNR1_ex2<-function(Z,x,D,R,w1,b1,b3,b_w1,b_b1,b_b3,lambda){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  p<-1
  alpha<-getTensor(w1,b1,b1,b3)
  beta<-getTensor(b_w1,b_b1,b_b1,b_b3)
  ## optimize
  Niter<-40; outeriter<-0;COND<-FALSE
  mu_t<-4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]/10000
  
  while(!COND){
    outeriter<-outeriter+1
    alpha_prev<-alpha
    beta_prev<-beta
    
    optim_result<-nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp<-matrix(optim_result$estimate,ncol=R)
    optim_result<-nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp<-matrix(optim_result$estimate,ncol=R)
    w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1<-apply(b1_temp,2,scale_vec)
    b3<-apply(b3_temp,2,scale_vec)
    alpha<-getTensor(w1,b1,b1,b3)
    
    niter<-20; iter<-0; cond<-FALSE
    while(!cond){
      #b1_temp<-matrix(0,n,R); b3_temp<-matrix(0,K,R)
      b_b1_temp<-matrix(0,n,R); b_b3_temp<-matrix(0,K,R)
      iter<-iter+1
      beta_old<-beta
      diff<-(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta_old),D,4,2)))
      for (r in 1:R) {
        b_gb1<-2*apply(diff*outer(x,(rep(1,n)%o%((b_b1[,r]*b_w1[r])%o%as.vector(D%*%b_b3[,r])))),2,sum)
        b_gb3<-apply(diff*outer(x,outer(b_b1[,r]%o%(b_b1[,r]*b_w1[r]),rep(1,T))),4,sum)%*%D
        b_b1_temp[,r]<-lasso(b_b1[,r]+b_gb1*mu_t,mu_t*lambda)
        b_b3_temp[,r]<-b_b3[,r]+b_gb3*mu_t
      }
      b_w1<-b_w1*apply(b_b1_temp,2,norm_vec)*apply(b_b1_temp,2,norm_vec)*apply(b_b3_temp,2,norm_vec)
      b_b1<-apply(b_b1_temp,2,scale_vec)
      b_b3<-apply(b_b3_temp,2,scale_vec)
      beta<-getTensor(b_w1,b_b1,b_b1,b_b3)
      cond<-(iter > niter) #| ascent<=EPS
    }
    
    ascent<-loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha_prev,beta_prev)
    COND<-(outeriter > Niter) | ascent<=EPS
    print(c(outeriter,ascent))
  }
  
  
  BIC<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta)>0)/2)
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha,beta=beta))  
}


lambda.select_ex2<-function(Z,x,D,w1,b1,b3,b_w1,b_b1,b_b3,R,lambda.vec){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  p=1
  alpha<-getTensor(w1,b1,b1,b3)  
  beta0<-getTensor(b_w1,b_b1,b_b1,b_b3)
  
  mu_t<-4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]/10000
  m<-length(lambda.vec)
  beta.list<-list()
  b_w1.list <- list()
  b_b1.list <- list()
  b_b3.list <- list()
  eBIC.vec<-c()
  for(i in 1:m){
    beta<-beta0
    lambda<-lambda.vec[i]
    niter<-20; iter<-0;cond<-FALSE
    while(!cond){
      iter<-iter+1
      b_b1_temp<-matrix(0,n,R); b_b3_temp<-matrix(0,K,R)
      beta_old<-beta
      diff<-(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta_old),D,4,2)))
      for (r in 1:R) {
        b_gb1<-2*apply(diff*outer(x,(rep(1,n)%o%((b_b1[,r]*b_w1[r])%o%as.vector(D%*%b_b3[,r])))),2,sum)
        b_gb3<-apply(diff*outer(x,outer(b_b1[,r]%o%(b_b1[,r]*b_w1[r]),rep(1,T))),4,sum)%*%D
        b_b1_temp[,r]<-lasso(b_b1[,r]+b_gb1*mu_t,mu_t*lambda)
        b_b3_temp[,r]<-b_b3[,r]+b_gb3*mu_t
      }
      b_w1<-b_w1*apply(b_b1_temp,2,norm_vec)*apply(b_b1_temp,2,norm_vec)*apply(b_b3_temp,2,norm_vec)
      b_b1<-apply(b_b1_temp,2,scale_vec)
      b_b3<-apply(b_b3_temp,2,scale_vec)
      beta<-getTensor(b_w1,b_b1,b_b1,b_b3)
      ascent<-sum((beta-beta_old)^2)
      #loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha,beta_old)
      cond<-(iter > niter) | ((ascent)<=EPS)
      print(c(iter,ascent))
    }
    
    b_w1.list[[i]] <-b_w1 
    b_b1.list[[i]] <-b_b1
    b_b3.list[[i]] <-b_b3 
    beta.list[[i]]<-beta
    eBIC.vec[i]<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta)>0)/2)
  }
  plot(lambda.vec, eBIC.vec,xlab="sparsity", main="eBIC for sparsity selection")
  return(list(lambda=lambda.vec[which.min(eBIC.vec)], eBIC.vec=eBIC.vec,b_w1=b_w1.list[[which.min(eBIC.vec)]],b_b1=b_b1.list[[which.min(eBIC.vec)]],b_b3=b_b3.list[[which.min(eBIC.vec)]],beta=beta.list[[which.min(eBIC.vec)]],beta.list=beta.list,b_w1.list=b_w1.list,b_b1.list=b_b1.list,b_b3.list=b_b3.list))
}






b0_max_e<-function(alpha,Z,x,D,beta){
  fn<-0
  N<-dim(Z)[1]
  #K<-dim(D)[2]
  #alpha<-rep(0,K)
  fn<--sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2))))
  #diff<-apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2)),2,sum)
  #attr(fn,"gradient")<-b1_gr(b1,D,w1,b3,diff)
  fn
}

loglikelihood_entry<-function(Z,x,D,alpha,beta){
  N<-dim(Z)[1]
  T<-dim(Z)[4]
  sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2))))/N
}



EW0<-function(Z,x,D){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  alpha_new_entry<-beta_new_entry<-alpha_old_entry<-beta_old_entry<-array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      alpha_new_entry[j,jj,]<-rep(0,K)
      beta_new_entry[j,jj,]<-rep(0,K)
      niter<-40; iter<-0; cond<-FALSE
      while (!cond) {
        iter<-iter+1
        alpha_old_entry[j,jj,]<-alpha_new_entry[j,jj,]
        beta_old_entry[j,jj,]<-beta_new_entry[j,jj,]
        alpha_new_result<-nlm(b0_max_e,alpha_old_entry[j,jj,],Z=Z[,j,jj,],x=x,D=D,beta=beta_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        if(max(abs(alpha_new_result$estimate))>10){
          alpha_new_entry[j,jj,]<-rep(0,K)
        }else{
          alpha_new_entry[j,jj,]<-alpha_new_result$estimate}
        alpha_new_entry[jj,j,]<-alpha_new_entry[j,jj,]
        beta_new_result<-nlm(b0_max_e,beta_old_entry[j,jj,],Z=Z[,j,jj,],x=x,D=D,alpha=alpha_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        if(max(abs(beta_new_result$estimate)>10)){
          beta_new_entry[j,jj,]<-rep(0,K)
        }else{
          beta_new_entry[j,jj,]<-beta_new_result$estimate}
        beta_new_entry[jj,j,]<-beta_new_entry[j,jj,]
        ascent<-loglikelihood_entry(Z[,j,jj,],x,D,alpha_new_entry[j,jj,],beta_new_entry[j,jj,])-loglikelihood_entry(Z[,j,jj,],x,D,alpha_old_entry[j,jj,],beta_old_entry[j,jj,])
        cond<-(iter > niter) | ascent<=EPS
        print(c(iter,ascent))
      }
    }
  }
  return(list(alpha_new_entry,beta_new_entry))
}


#####
EW_BASIC<-function(Z,x){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  B0<-B1<-B1_pvalue<-array(0,c(n,n,T))
  signal<-matrix(0,n,n)
  B1_pvalue_min<-matrix(0,n,n)
  B1_pvalue_adjust<-array(0,c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      for (t in 1:T){
        fit<-glm(Z[,j,jj,t]~x,family=binomial(link = logit))
        B0[j,jj,t]<-fit$coefficients[1]
        B0[jj,j,t]<-B0[j,jj,t]
        B1[j,jj,t]<-fit$coefficients[2]
        B1[jj,j,t]<-B1[j,jj,t]
        B1_pvalue[j,jj,t]<-summary(fit)$coefficient[2,4]
        B1_pvalue[jj,j,t]<- B1_pvalue[j,jj,t]
      }
    }
  }
  B1_pvalue_adjust<-array(p.adjust(B1_pvalue,"bonferron"),c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      if(min(B1_pvalue_adjust[j,jj,])<=0.05){
        signal[j,jj]<-1
        signal[jj,j]<-signal[j,jj]
      }
    }
  }
  return(list(B0=B0,B1=B1,B1_pvalue_min=B1_pvalue_min,B1_pvalue_adjust=B1_pvalue_adjust,signal=signal))
}



