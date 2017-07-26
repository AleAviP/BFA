require(MASS)
require(psych)
require(matrixcalc)
require(Matrix)
BFA.EM.BE.P <- function(x,v,b=NULL,q=2,eps=0.0001,it=50,seed=4,scaling=FALSE,init=NULL,hyper=list(theta=0.5,aS = 1,bS = 1,l0=0.1,l1=100,at=1,bt=1,epsM=0.001),varianceBE=FALSE,intercept=FALSE){
  
  set.seed(seed)
  
  if (scaling==TRUE) x<-scale(x)
  
  n<-nrow(x)
  p<-ncol(x)
  p_b<-1
  w_b<-1
  if (!is.null(b)) {p_b=ncol(b)
  w_b=apply(b,2,sum)/n
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
    stima=try(factanal(fm1$residuals,r,rotation="none"),silent=TRUE) 
    psi=diag(p)
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
    
  } else {
    theta<-init$theta
    psi<-init$psi
    M<-init$M
    z<-init$z}
  
  ##traces
  trace_psi <- matrix(nrow =1, ncol = p*p_b,diag(psi))
  trace_theta <- matrix(nrow =1, ncol = p*p_v,c(theta))
  trace_M <- matrix(nrow = 1, ncol = p*q,c(M))
  
  #####EM#####
  count<-1
  likelihood<-NULL 
  change<-1000
  lik<--100000000000
  x_o <- x
  x<-x_o-v%*%t(theta)
  
  trace_like <- lik
  
  if(varianceBE==FALSE){
    psi_aux=vector()
  }else{psi_aux=matrix(1,p,p_b)}
  
  while ((count < it) & (change > eps )) {
    ##E step
    Ez<-E.Z(x,as.matrix(M),psi,q)
    Ez.x<-Ez$z
    W<-Ez$W
    
    ##M-step
    M<-M.new(x,Ez.x,W,n)
    
    if (varianceBE==FALSE){
      psi<-Psi.new(x,Ez.x,M,aS,bS)
      psi_aux<-diag(psi)
    }else{
      for (i in 1:p_b){
        x_b<-x[c(1:n)[b[,i]==1],]
        z_b<-Ez.x[c(1:n)[b[,i]==1],]
        psi_aux[,i]<-diag(Psi.new(x_b,z_b,M,aS,bS))
      }
      psi<-diag(c(psi_aux %*% w_b))
    }
    
    theta<-Theta.new(x_o,v,Ez.x,M)
    
    x<-x_o-v%*%theta
    
    #Saving values
    trace_psi = rbind(trace_psi,c(psi_aux))
    trace_theta = rbind(trace_theta,c(t(theta)))
    trace_M = rbind(trace_M,c(M))
    
    ### Rotation
    #A<-chol(W+(t(Ez.x)%*%(Ez.x)*1/n))
    #M <- M%*%A
    
    ##Likelihood
    likelihood<-likelihoodFA(x,M,psi)
    change<-likelihood-lik
    lik<-likelihood
    trace_like<-c(trace_like,likelihood)
    count<-count+1
  }
  
  return(list(M=M,
              Psi=psi_aux,
              Ez=Ez.x,
              Theta=t(theta),
              tracePsi=trace_psi,
              traceTheta=trace_theta,
              traceM=trace_M,
              iterations=count,
              like=trace_like[-1]))
}