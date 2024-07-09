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





loglikelihood<-function(Z,x,D,alpha,beta){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  sum(Z*tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2)-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,4,2))))/N
}

loglikelihood_ex1<-function(Z,x,xt,D,alpha,beta1,beta2){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(beta2, D, 3, 2)[,,t]  
    C[,,,t] <- outer(xt[,t],tensor_t)
  }
  eta <- tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,4,2)+C
  sum(Z*eta-log(1+exp(eta)))/N
}


generateTensor_ex1<-function(x,xt,D,alpha,beta1,beta2){
  N<-length(x)
  n<-dim(alpha)[1]
  T<-dim(D)[1]
  l<-n*n*T
  
  Z<-array(0,c(N,n,n,T))
  for(i in 1:N){
    Ci <- array(0,c(n,n,T))
    for(t in 1:T){
      Ci[,,t] <- tensor(beta2,D,3,2)[,,t]*xt[i,t]
    }
    Zi<-new("Tensor",3L,c(n,n,T),data=rbinom(n=l,size=1,prob=logistic(vec(as.tensor(tensor(alpha+x[i]*beta1,D,3,2)+Ci)))))
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


####find the initial B0
findinitial<-function(Z,x,D){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  ZD<-logit(apply(Z,c(2,3,4),sum)/N)
  ZD[ZD>4]<-4; ZD[ZD<(-4)]<-(-4)
  alpha<-array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      alpha[j,jj,]<-lm(ZD[j,jj,]~D-1)$coefficients
      alpha[jj,j,]<-alpha[j,jj,]
    }
  }
  return(alpha)
}



b1_gr<-function(b1,D,w1,b3,diff){
  R<-length(w1)
  n<-dim(b1)[1]
  gb1<-matrix(0,n,length(w1))
  for (r in 1:length(w1)){
    gb1[,r]<-2*apply(diff*(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r]))),1,sum)
  }
  -gb1
  
}






b1_max_ex1<-function(b1,Z,x,xt,D,w1,b3,beta1,beta2){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(w1)
  K<-dim(b3)[1]
  T<-dim(Z)[4]
  
  b1<-matrix(b1,ncol=R)
  alpha<-array(0,c(n,n,K))
  for (r in 1:R){
    alpha<-alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(beta2, D, 3, 2)[,,t]  
    C[,,,t] <- outer(xt[,t],tensor_t)
  }
  eta <- tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,4,2)+C
  sum(Z*eta-log(1+exp(eta)))/N
  
  fn<--sum(Z*eta-log(1+exp(eta)))
  
  diff<-apply(Z-logistic(eta),c(2,3,4),sum)
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


b3_max_ex1<-function(b3,K,T,Z,x,xt,D,w1,b1,beta1,beta2){
  fn<-0
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  R<-length(w1)
  T<-dim(Z)[4]
  
  b3<-matrix(b3,ncol=R)
  alpha<-array(0,c(n,n,K))
  for (r in 1:R){
    alpha<-alpha+outer((b1[,r]%o%b1[,r]),b3[,r])*w1[r]
  }
  
  C <- array(0,c(N,n,n,T))
  for (t in 1:T) {
    tensor_t <- tensor(beta2, D, 3, 2)[,,t]  # 假设这是一个合适的函数或矩阵
    C[,,,t] <- outer(xt[,t],tensor_t)
  }
  eta <- tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,4,2)+C
  
  fn<--sum(Z*eta-log(1+exp(eta)))
  
  diff<-apply(Z-logistic(eta),c(2,3,4),sum)
  attr(fn,"gradient")<-b3_gr(b3,K,T,D,w1,b1,diff)
  fn
  
}


glasso<-function(beta,lambda){
  K<-dim(beta)[3]
  m<-(1-lambda/apply(beta,c(1,2),norm_vec))
  ind<-m*(m>0)
  beta*outer(ind,rep(1,K))
}

VCNR0<-function(Z,x,D,R,initial,step=1/300/300/10){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  ## initialize
  Niter<-30; outeriter<-0;COND<-FALSE
  while(!COND){
    outeriter<-outeriter+1
    output<-cp(as.tensor(initial),R,max_iter = 50)
    w1<-output[[1]];b1<-output$U[[1]];b2<-output$U[[2]];b3<-output$U[[3]]
    b3<-sweep(b3,2,sign(b1[1,])*sign(b2[1,]),"*")
    beta<-array(0,c(n,n,K))
    
    niter<-20; iter<-0; cond<-FALSE
    while(!cond){
      iter<-iter+1
      alpha_old<-getTensor(w1,b1,b1,b3)
      b1_temp<-matrix(0,n,R); b3_temp<-matrix(0,K,R)
      
      diff<-apply(Z-logistic(tensor(outer(rep(1,N),alpha_old)+outer(x,beta),D,4,2)),c(2,3,4),sum)
      for (r in 1:R){
        gb1<-2*apply(diff*(rep(1,n)%o%((b1[,r]*w1[r])%o%as.vector(D%*%b3[,r]))),1,sum)
        gb3<-apply(diff*outer(b1[,r]%o%(b1[,r]*w1[r]),rep(1,T)),3,sum)%*%D
        b1_temp[,r]<-b1[,r]+gb1*step
        b3_temp[,r]<-b3[,r]+gb3*step
      }
      w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
      b1<-apply(b1_temp,2,scale_vec)
      b3<-apply(b3_temp,2,scale_vec)
      alpha<-getTensor(w1,b1,b1,b3)
      cond<-(iter > niter) | sum((alpha-alpha_old)^2)> 10^10
      print(c(iter,sum((alpha-alpha_old)^2)))
    }
    COND<-(outeriter > Niter)| sum((alpha-alpha_old)^2)< 10^10
  }
  
  ## optimize
  niter<-30; iter<-0; cond<-FALSE
  while(!cond){
    iter<-iter+1
    alpha_old<-getTensor(w1,b1,b1,b3)
    optim_result<-nlm(b1_max,b1,Z=Z,x=x,D=D,w1=w1,b3=b3,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp<-matrix(optim_result$estimate,ncol=R)
    optim_result<-nlm(b3_max,b3,K=K,T=T,Z=Z,x=x,D=D,w1=w1,b1=b1,beta=beta,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp<-matrix(optim_result$estimate,ncol=R)
    w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1<-apply(b1_temp,2,scale_vec)
    b3<-apply(b3_temp,2,scale_vec)
    alpha<-getTensor(w1,b1,b1,b3)
    ascent<-sum((alpha-alpha_old)^2)#loglikelihood(Z,x,D,alpha,beta)-loglikelihood(Z,x,D,alpha_old,beta)#sum((alpha-alpha_old)^2)#
    cond<-(iter > niter) | ascent<=EPS
    print(c(iter,ascent))
  }
  
  BIC<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N*n*n*T)+log(n*n*K))*(R*(n+K))
  #BIC1<--N*loglikelihood(Z,x,D,alpha,beta)+(log(N)+log(n*n*K))*(R*(n+K))
  #loglike<--N*loglikelihood(Z,x,D,alpha,beta)
  #penalty<-(log(N*n*n*T)+log(n*n*K))*(R*(n+K))
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha))  
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

b0_max_ex1<-function(alpha,Z,x,xt,D,beta1,beta2){
  fn<-0
  N<-dim(Z)[1]
  #K<-dim(D)[2]
  #alpha<-rep(0,K)
  fn<--sum(Z*(tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,2,2)+(xt*outer(rep(1,N),as.vector(t(beta2)%*%t(D)))))-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,2,2)+(xt*outer(rep(1,N),as.vector(t(beta2)%*%t(D)))))))
  #diff<-apply(Z-logistic(tensor(outer(rep(1,N),alpha)+outer(x,beta),D,2,2)),2,sum)
  #attr(fn,"gradient")<-b1_gr(b1,D,w1,b3,diff)
  fn
}



####FISTA without refit
VCNR2_ex1<-function(Z,x,xt,D,R,w1,b1,b3,beta1,beta2,lambda){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  p<-2
  alpha<-getTensor(w1,b1,b1,b3)
  
  ## optimize
  Niter<-40; outeriter<-0;COND<-FALSE
  mu_t<-4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]
  
  while(!COND){
    outeriter<-outeriter+1
    alpha_prev<-alpha
    beta1_prev<-beta1
    beta2_prev<-beta2
    
    optim_result<-nlm(b1_max_ex1,b1,Z=Z,x=x,xt=xt,D=D,w1=w1,b3=b3,beta1=beta1,beta2=beta2,check.analyticals=FALSE,gradtol=1e-3)
    b1_temp<-matrix(optim_result$estimate,ncol=R)
    optim_result<-nlm(b3_max_ex1,b3,K=K,T=T,Z=Z,x=x,xt=xt,D=D,w1=w1,b1=b1,beta1=beta1,beta2=beta2,check.analyticals=FALSE,gradtol=1e-3)
    b3_temp<-matrix(optim_result$estimate,ncol=R)
    w1<-w1*apply(b1_temp,2,norm_vec)*apply(b1_temp,2,norm_vec)*apply(b3_temp,2,norm_vec)
    b1<-apply(b1_temp,2,scale_vec)
    b3<-apply(b3_temp,2,scale_vec)
    alpha<-getTensor(w1,b1,b1,b3)
    
    niter<-30; iter<-0; cond<-FALSE;h=1
    while(!cond){
      h_old<-h
      iter<-iter+1
      beta1_old<-beta1
      beta2_old<-beta2
      
      C <- array(0,c(N,n,n,T))
      for (t in 1:T) {
        tensor_t <- tensor(beta2, D, 3, 2)[,,t] 
        C[,,,t] <- outer(xt[,t],tensor_t)
      }
      
      eta <- tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,4,2)+C
      
      
      diffx1<-(Z-logistic(eta))*outer(x,rep(1,n)%o%rep(1,n)%o%rep(1,T))
      gbeta1<-apply(tensor(diffx1,D,4,1),c(2,3,4),sum)
      beta1<-glasso(beta1_old+gbeta1*mu_t,mu_t*lambda)   
      h=(1+sqrt(1+4*h_old^2))/2
      beta1=beta1-(h_old-1)/h*(beta1-beta1_old)
      
      
      diffx2<-(Z-logistic(eta))*aperm(outer(xt,rep(1,n)%o%rep(1,n)),c(1,3,4,2))
      gbeta2<-apply(tensor(diffx2,D,4,1),c(2,3,4),sum)
      beta2<-glasso(beta2_old+gbeta2*mu_t,mu_t*lambda)   
      h=(1+sqrt(1+4*h_old^2))/2
      beta2=beta2-(h_old-1)/h*(beta2-beta2_old)
      
      cond<-(iter > niter) #| ascent<=EPS
    }
    
    ascent<-loglikelihood_ex1(Z,x,xt,D,alpha,beta1,beta2)-loglikelihood_ex1(Z,x,xt,D,alpha_prev,beta1_prev,beta2_prev)
    COND<-(outeriter > Niter) | ascent<=EPS
    print(c(outeriter,ascent))
  }
  
  BIC<--N*loglikelihood_ex1(Z,x,xt,D,alpha,beta1,beta2)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta1)>0)/2+sum(abs(beta2)>0)/2)
  return(list(BIC=BIC,w=w1,b1=b1,b3=b3,alpha=alpha,beta1=beta1,beta2=beta2))  
}

############FISTA with refit






######lambda select
lambda.select_ex1<-function(Z,x,xt,D,w1,b1,b3,beta0,R,lambda.vec){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  p=2
  alpha<-getTensor(w1,b1,b1,b3)  
  # print(-N*loglikelihood(Z,x,D,alpha,beta0))
  mu_t1<-4/as.vector(t(x)%*%x)/eigen(t(D)%*%D)$values[1]
  mu_t2<-mu_t1
  #mu_t<-5e-2
  m<-length(lambda.vec)
  beta1.list<-beta2.list<-list()
  eBIC.vec<-c()
  for(i in 1:m){
    beta1<-beta0
    beta2<-beta0
    lambda<-lambda.vec[i]
    niter<-30; iter<-0;cond<-FALSE
    while(!cond){
      iter<-iter+1
      beta1_old<-beta1
      beta2_old<-beta2
      
      C <- array(0,c(N,n,n,T))
      for (t in 1:T) {
        tensor_t <- tensor(beta2, D, 3, 2)[,,t]  # 假设这是一个合适的函数或矩阵
        C[,,,t] <- outer(xt[,t],tensor_t)
      }
      eta <- tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,4,2)+C
      
      
      diffx1<-(Z-logistic(eta))*outer(x,rep(1,n)%o%rep(1,n)%o%rep(1,T))
      gbeta1<-apply(tensor(diffx1,D,4,1),c(2,3,4),sum)
      beta1<-glasso(beta1_old+gbeta1*mu_t1,mu_t1*lambda)   
      #h=(1+sqrt(1+4*h_old^2))/2
      #beta1=beta1-(h_old-1)/h*(beta1-beta1_old)
      
      
      diffx2<-(Z-logistic(eta))*aperm(outer(xt,rep(1,n)%o%rep(1,n)),c(1,3,4,2))
      gbeta2<-apply(tensor(diffx2,D,4,1),c(2,3,4),sum)
      beta2<-glasso(beta2_old+gbeta2*mu_t2,mu_t2*lambda)   
      #h=(1+sqrt(1+4*h_old^2))/2
      #beta2=beta2-(h_old-1)/h*(beta2-beta2_old)
      

      ascent<-loglikelihood_ex1(Z,x,xt,D,alpha,beta1,beta2)-loglikelihood_ex1(Z,x,xt,D,alpha,beta1_old,beta2_old)
      cond<-(iter > niter) | (abs(ascent)<=EPS)
      print(c(iter,ascent))
    }
    
    beta1.list[[i]]<-beta1
    beta2.list[[i]]<-beta2
    eBIC.vec[i]<--N*loglikelihood_ex1(Z,x,xt,D,alpha,beta1,beta2)+(log(N*n*n*T)+log(n*n*K*(p+1)))*(R*(n+K)+sum(abs(beta1)>0)/2+sum(abs(beta2)>0)/2)
  }
  #plot(lambda.vec, eBIC.vec,xlab="sparsity", main="eBIC for sparsity selection")
  return(list(lambda=lambda.vec[which.min(eBIC.vec)], eBIC.vec=eBIC.vec,beta1=beta1.list[[which.min(eBIC.vec)]],beta2=beta2.list[[which.min(eBIC.vec)]],beta1.list=beta1.list,beta2.list=beta2.list))
}




#####


loglikelihood_entry_ex1<-function(Z,x,xt,D,alpha,beta1,beta2){
  N<-dim(Z)[1]
  T<-dim(Z)[4]
  sum(Z*(tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,2,2)+(xt*outer(rep(1,N),as.vector(t(beta2)%*%t(D)))))-log(1+exp(tensor(outer(rep(1,N),alpha)+outer(x,beta1),D,2,2)+(xt*outer(rep(1,N),as.vector(t(beta2)%*%t(D)))))))/N
}


EW0_ex1<-function(Z,x,xt,D){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  K<-dim(D)[2]
  alpha_new_entry<-beta1_new_entry<-beta2_new_entry<-alpha_old_entry<-beta1_old_entry<-beta2_old_entry<-array(0,c(n,n,K))
  for(j in 1:n){
    for(jj in j:n){
      alpha_new_entry[j,jj,]<-rep(0,K)
      beta1_new_entry[j,jj,]<-rep(0,K)
      beta2_new_entry[j,jj,]<-rep(0,K)
      niter<-40; iter<-0; cond<-FALSE
      while (!cond) {
        iter<-iter+1
        alpha_old_entry[j,jj,]<-alpha_new_entry[j,jj,]
        beta1_old_entry[j,jj,]<-beta1_new_entry[j,jj,]
        beta2_old_entry[j,jj,]<-beta2_new_entry[j,jj,]
        alpha_new_result<-nlm(b0_max_ex1,alpha_old_entry[j,jj,],Z=Z[,j,jj,],x=x,xt=xt,D=D,beta1=beta1_new_entry[j,jj,],beta2=beta2_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-2)
        if(max(abs(alpha_new_result$estimate))>10){
          alpha_new_entry[j,jj,]<-rep(0,K)
        }else{
          alpha_new_entry[j,jj,]<-alpha_new_result$estimate}
        alpha_new_entry[jj,j,]<-alpha_new_entry[j,jj,]
        beta1_new_result<-nlm(b0_max_ex1,beta1_old_entry[j,jj,],Z=Z[,j,jj,],x=x,xt=xt,D=D,alpha=alpha_new_entry[j,jj,],beta2=beta2_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        #if(max(abs(beta_new_result$estimate)>10)){
         # beta_new_entry[j,jj,]<-rep(0,K)
        #}else{
        beta1_new_entry[j,jj,]<-beta1_new_result$estimate
        beta1_new_entry[jj,j,]<-beta1_new_entry[j,jj,]
        
        beta2_new_result<-nlm(b0_max_ex1,beta2_old_entry[j,jj,],Z=Z[,j,jj,],x=x,xt=xt,D=D,alpha=alpha_new_entry[j,jj,],beta1=beta1_new_entry[j,jj,],check.analyticals=FALSE,gradtol=1e-3)
        #if(max(abs(beta_new_result$estimate)>10)){
        # beta_new_entry[j,jj,]<-rep(0,K)
        #}else{
        beta2_new_entry[j,jj,]<-beta2_new_result$estimate
        beta2_new_entry[jj,j,]<-beta2_new_entry[j,jj,]
        ascent<-loglikelihood_entry_ex1(Z[,j,jj,],x,xt,D,alpha_new_entry[j,jj,],beta1_new_entry[j,jj,],beta2_new_entry[j,jj,])-loglikelihood_entry_ex1(Z[,j,jj,],x,xt,D,alpha_old_entry[j,jj,],beta1_old_entry[j,jj,],beta2_old_entry[j,jj,])
        cond<-(iter > niter) | ascent<=EPS
        print(c(iter,ascent))
      }
    }
  }
  return(list(alpha_new_entry,beta1_new_entry,beta2_new_entry))
}



############
EW_BASIC_ex1<-function(Z,x,xt){
  N<-dim(Z)[1]
  n<-dim(Z)[2]
  T<-dim(Z)[4]
  B0<-B1<-B2<-B1_pvalue<-B2_pvalue<-array(0,c(n,n,T))
  signal1<-signal2<-matrix(0,n,n)
  B1_pvalue_adjust<-B2_pvalue_adjust<-array(0,c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      for (t in 1:T){
        fit<-glm(Z[,j,jj,t]~x+xt[,t],family=binomial(link = logit))
        B0[j,jj,t]<-fit$coefficients[1]
        B0[jj,j,t]<-B0[j,jj,t]
        B1[j,jj,t]<-fit$coefficients[2]
        B1[jj,j,t]<-B1[j,jj,t]
        B2[j,jj,t]<-fit$coefficients[3]
        B2[jj,j,t]<-B2[j,jj,t]
        
        B1_pvalue[j,jj,t]<-summary(fit)$coefficient[2,4]
        B1_pvalue[jj,j,t]<- B1_pvalue[j,jj,t]
        B2_pvalue[j,jj,t]<-summary(fit)$coefficient[3,4]
        B2_pvalue[jj,j,t]<- B2_pvalue[j,jj,t]
      }
    }
  }
  B1_pvalue_adjust<-array(p.adjust(B1_pvalue,"bonferron"),c(n,n,T))
  B2_pvalue_adjust<-array(p.adjust(B2_pvalue,"bonferron"),c(n,n,T))
  for(j in 1:n){
    for(jj in j:n){
      if(min(B1_pvalue_adjust[j,jj,])<=0.05){
        signal1[j,jj]<-1
        signal1[jj,j]<-signal1[j,jj]
      }else if(min(B2_pvalue_adjust[j,jj,])<=0.05){
        signal2[j,jj]<-1
        signal2[jj,j]<-signal2[j,jj]
      }
    }
  }
  return(list(B0=B0,B1=B1,B2=B2,B1_pvalue_adjust=B1_pvalue_adjust,B2_pvalue_adjust=B2_pvalue_adjust,signal1=signal1,signal2=signal2))
}



