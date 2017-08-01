#####E-step###
E.Z <- function(x,M,psi,q,p){
  n<-nrow(x)
  I.q <- diag(q)
  tM <- t(M)
  #psiInv <- solve(psi)
  psiInv <- diag(1/diag(psi),p)
  #psiInv <- ginv(psi)
  W <- solve(I.q + tM%*%psiInv%*%M)
  Z <- W%*%tM%*%psiInv%*%t(x)
  Ezz.x <- Z%*%t(Z)+n*W
  return(list(z=t(Z),W=W,zz=Ezz.x))
}

E.gamma.i <- function(m.i,psi.i,l0,l1,theta.i){
  dens0 = dnorm(m.i,0,sqrt(psi.i*l0))
  dens1 = dnorm(m.i,0,sqrt(psi.i*l1))
  p.i =dens1*theta.i/(dens0*(1-theta.i)+dens1*theta.i)
  d.i = (1-p.i)/l0 +  p.i/l1
  return(list(p=p.i,d=d.i))
}

E.gamma.2 <- function(M,psi,l0,l1,theta){
  p = nrow(M)
  Gamma.aux = llply(.data = 1:p, .fun = function(y){
    E.gamma.i(M[y,],diag(psi)[y],l0,l1,theta)
  }, .parallel = FALSE)
  return(Reduce(function(x,y) Map(rbind, x, y),Gamma.aux))
}

#####M-step####
Theta.new<-function(x_o,v,Ez.x,M){
  tv <- t(v)
  solve(tv%*%v)%*%tv%*%(x_o-Ez.x%*%t(M))
}

M.pMOM<-function(x,Ez.x,Ezz.x,D,M.old,p.gamma,sigma){
  p=ncol(x)
  xz <- t(x)%*%Ez.x
  M.old.inv<-1/M.old
  W.M <- 2*M.old.inv*p.gamma
  #Replace inf for ones
  M.old.inv[is.infinite(M.old.inv)] <- 1
  Maux = llply(.data = 1:p, .fun = function(y){
    (xz[y,]+W.M[y,]*sigma[y])%*%solve(diag(D[y,])+Ezz.x)
  }, .parallel = FALSE)
  Mnew = do.call("rbind", Maux)
  #Mnew <- t(x)%*%Ez.x%*%solve(Ezz.x)
  return(Mnew)
}

Psi.pMOM<-function(x,Ez.x,Ezz.x,M,D,a,b,p.gamma){
  n = nrow(x)
  q = ncol(M)
  p = nrow(M)
  W =  apply(p.gamma,1,sum)
  #diag((1/(n+a+2+q))*diag(t(x)%*%x-2*t(x)%*%Ez.x%*%t(M)+M%*%Ezz.x%*%t(M)+a*b+diag(apply(D*M^2,1,sum))))
  diag((1/((n+a+2+q)+W))*diag(t(x)%*%x-2*t(x)%*%Ez.x%*%t(M)+M%*%Ezz.x%*%t(M)+a*b+diag(apply(D*M^2,1,sum))))
}

theta.gamma.new.2 <- function(P,a,b){
  p = nrow(P)
  q = ncol(P)
  j = c(1:q)
  theta = (apply(P,2,sum)+a/j-1)/(a/j+b+p-2)
  #Restrict values between 0 and 1
  theta[theta>1]<-1
  theta[theta<0]<-0
  return(theta)
}

lof<-function(matrix,reference){
  matrix[,order(apply(reference,2,sum),decreasing=TRUE)]
}

likelihoodFA <- function(x,M,psi){
  p = ncol(x)
  n = nrow(x)
  q = ncol(M)
  C = M %*% t(M) + psi
  #l = -n/2 *(p*log(2 * pi )+log(det(C)) + tr(solve(C)%*%cov(x)))
  l = sum(dmnorm(x, mean = rep(0, p), C, log = TRUE))
  return(l)
}

logpriorFApMOM <- function(M,psi,gtheta,gamma,D,hyper,varianceBE=FALSE,wb=1){
  q=ncol(M)
  p=nrow(M)
  l.psi=0
  l.psi=sum(log(dinvgamma(psi,hyper$aS/2,hyper$aS*hyper$bS/2)))
  sum.gamma=apply(gamma,2,sum)
  BetaBin.aux = llply(.data = 1:q, .fun = function(y){
    log(beta(sum.gamma[y]+hyper$at/y,p-sum.gamma[y]+hyper$bt))-log(beta(hyper$at/y,hyper$at))
  }, .parallel = FALSE)
  BetaBin.aux = do.call("rbind", BetaBin.aux)
  # Stirling approximation
  BetaBin.approx = llply(.data = 1:q, .fun = function(y){
    0.5*log(2*pi)+(sum.gamma[y]+hyper$at/y-0.5)*log(sum.gamma[y]+hyper$at/y)+
      +(p-sum.gamma[y]+hyper$bt-0.5)*log(p-sum.gamma[y]+hyper$bt)-
      (p+hyper$at/y+hyper$bt-0.5)*log(p+hyper$at/y+hyper$bt)
  }, .parallel = FALSE)
  BetaBin.approx = do.call("rbind", BetaBin.approx)
  BetaBin.aux[is.infinite(BetaBin.aux)]=BetaBin.approx[is.infinite(BetaBin.aux)]
  l.gamma.gtheta=sum(BetaBin.aux)
  if(varianceBE==TRUE){psi=psi%*%wb}
  l.M=do.call("rbind",llply(.data = 1:p, .fun = function(y){
    dmnorm(M[y,],mean=rep(0,q),diag(psi[y]*D[y,]),log=TRUE)
  }, .parallel = FALSE))
  l.M=sum(l.M)
  l.pMOM = sum(log(M^2)[,!gtheta==0]*gamma[,!gtheta==0])-sum(t(log(psi*hyper$l1))%*%gamma)
  l.total=l.psi+l.gamma.gtheta+l.M+l.pMOM
  return(list(lPsi=l.psi,
              lpMOM =l.pMOM,
              lGammaGtheta=l.gamma.gtheta,
              lM=l.M,
              lTotal=l.total))
}

Mreplace<-function(Mold,Mnew,i){
  M2<-Mold
  M2[i,]<-Mnew[i,]
  return(M2)
}