require(MASS)
require(psych)
require(matrixcalc)
require(Matrix)
library(plyr)
library(LaplacesDemon)
library(mnormt)
library(OpenMx)

BFA.EM.BE.P.2 <- function(x,v,b=NULL,q=2,eps=0.0001,it=50,seed=4,scaling=FALSE,init=NULL,hyper=list(theta=0.5,aS = 1,bS = 1,l0=0.1,l1=100,at=1,bt=1,epsM=0.001),varianceBE=FALSE,intercept=FALSE){
  
  set.seed(seed)
  
  if (scaling==TRUE) x<-scale(x)
  
  n<-nrow(x)
  p<-ncol(x)
  p_b<-1
  w_b<-1
  if (!is.null(b)) {p_b=ncol(b)
  n_b=apply(b,2,sum)
  w_b=n_b/n
  if (intercept==FALSE){v=cbind(v,b)}
  }
  if (intercept==TRUE){
    v_0=rep(1,n)
    if(!is.null(b)){
      v=cbind(v_0,v,b[,-p_b])
    }else{
      v=cbind(v_0,v)
    }
  }
  p_v<-ncol(v) 
  
  #Hyperparameters
  g.theta = hyper$theta
  aS = hyper$aS
  bS = hyper$bS
  l0 = hyper$l0
  l1 = hyper$l1
  at = hyper$at
  bt = hyper$bt
  epsM = hyper$epsM
  
  ## initialization of M Psi and theta
  if (is.null(init)) { 
    fm1 <- lm(x ~ v)
    theta = t(fm1$coefficients[-1,]) 
    #Replacing NA for zeros
    theta[is.na(theta)] <- 0
    stima<-try(factanal(fm1$residuals,r,rotation="none"),silent=TRUE) 
    psi<-diag(p)
    svd.X <- svd(fm1$residuals)
    svd.U <- svd.X$u
    svd.V <- svd.X$v
    Lambda <- diag(svd.X$d)
    sigma2 <- 1/(p-q) * sum(svd.U[(q+1):p])                                                          # variance; average variance associated with discarded dimensions
    aux_p <- min(n,p)
    if(is.positive.definite(Lambda-sigma2*diag(aux_p))==FALSE){
      aux.chol = nearPD(Lambda-sigma2*diag(aux_p))$mat
    }else{aux.chol = Lambda-sigma2*diag(aux_p)}
    M<-(svd.V %*% chol(aux.chol) %*% diag(aux_p))[,1:q]
    z<-mvrnorm(n = n, rep(0,q), diag(q))
    #g.theta <- rep(g.theta,p)
    g.theta <- rep(g.theta,q)
    gamma <- E.gamma.2(M,psi,l0,l1,g.theta)
    D <- gamma$d
    p.gamma <- gamma$p
    
  } else {
    theta<-init$theta
    psi<-init$psi
    M<-init$M
    z<-init$z
    D <- init$d
    p.gamma <- init$d
    g.theta <- initi$g.theta}
  
  ##traces
  trace_psi <- matrix(nrow =1, ncol = p*p_b,diag(psi))
  trace_theta <- matrix(nrow =1, ncol = p*p_v,c(theta))
  trace_M <- matrix(nrow = 1, ncol = p*q,c(M))
  trace_D <- matrix(nrow = 1, ncol = p*q,c(D))
  trace_p <- matrix(nrow = 1, ncol = p*q,c(p.gamma))
  trace_g.theta <- matrix(nrow = 1, ncol = q,c(g.theta))
  
  #####EM#####
  count<-1
  likelihood<-NULL 
  change<-1000
  lik<--100000000000
  x_o <- x
  x<-x_o-v%*%t(theta)
  I.q <- diag(q)
  
  trace_like <- lik
  
  if(varianceBE==FALSE){
    psi_aux=vector()
  }else{psi_aux=matrix(1,p,p_b)}
  
  while ((count < it) & (change > eps )) {
    ##E step
    Ez<-E.Z(x,as.matrix(M),psi,q)
    Ez.x<-Ez$z
    Ezz.x<-Ez$zz
    W<-Ez$W
    
    gamma <- E.gamma.2(M,psi,l0,l1,g.theta)
    D <- gamma$d
    p.gamma <- gamma$p
    #LOF
    D <- lof(D,p.gamma)
    p.gamma <- lof(p.gamma,p.gamma)
    
    ##M-step
    M<-M.new(x,Ez.x,Ezz.x,D)
    tM <- t(M)
    
    if (varianceBE==FALSE){
      psi<-Psi.new(x,Ez.x,Ezz.x,M,D,aS,bS)
      psi_aux<-diag(psi)
    }else{
      for (i in 1:p_b){
        x_b<-x[c(1:n)[b[,i]==1],]
        z_b<-Ez.x[c(1:n)[b[,i]==1],]
        psiInv <- solve(diag(psi_aux[,i]))
        W_b<-solve(I.q + tM%*%psiInv%*%M)
        zz_b<-t(z_b)%*%z_b+n_b[i]*W_b
        psi_aux[,i]<-diag(Psi.new(x_b,z_b,zz_b,M,D,aS,bS))
      }
      psi<-diag(c(psi_aux %*% w_b))
    }
    
    theta<-Theta.new(x_o,v,Ez.x,M)
    
    g.theta <- theta.gamma.new.2(p.gamma,at,bt)
    
    x<-x_o-v%*%theta
    
    #Saving values
    trace_psi <- rbind(trace_psi,c(psi_aux))
    trace_theta <- rbind(trace_theta,c(t(theta)))
    trace_M <- rbind(trace_M,c(M))
    trace_D <- rbind(trace_D,c(D))
    trace_p <- rbind(trace_p,c(p.gamma))
    trace_g.theta <- rbind(trace_g.theta,c(g.theta))
    
    ### Rotation
    #A<-chol(W+(t(Ez.x)%*%(Ez.x)*1/n))
    #M <- M%*%A
    
    ##Likelihood
    likelihood<-likelihoodFA(x,M,psi)
    change<-likelihood-lik
    lik<-likelihood
    trace_like<-c(trace_like,likelihood)
    count<-count+1
    #print(count)
    #print(change)
  }
  
  return(list(M=M,
              Psi=psi_aux,
              Ez=Ez.x,
              Theta=t(theta),
              gTheta = g.theta,
              gamma = p.gamma,
              D = D,
              tracePsi=trace_psi,
              traceTheta=trace_theta,
              traceM=trace_M,
              traceD=trace_D,
              traceGtheta=trace_g.theta,
              traceP=trace_p,
              iterations=count,
              like=trace_like[-1]))
}