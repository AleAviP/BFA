#####E-step###
E.Z <- function(x,M,psi,q){
  I.q <- diag(q)
  tM <- t(M)
  psiInv <- solve(psi)
  #psiInv <- ginv(psi)
  W <- solve(I.q + tM%*%psiInv%*%M)
  Z <- W%*%tM%*%psiInv%*%t(x)
  return(list(z=t(Z),W=W))
}

#####M-step####
Theta.new<-function(x_o,v,Ez.x,M){
  tv <- t(v)
  solve(tv%*%v)%*%tv%*%(x_o-Ez.x%*%t(M))
}

M.new<-function(x,Ez.x,W,n){
  Ezz.x <- t(Ez.x)%*%Ez.x+n*W
  M <- t(x)%*%Ez.x%*%solve(Ezz.x)
  return(M)
}

Psi.new<-function(x,Ez.x,M,a,b){
  n = nrow(x)
  diag((1/(n+a+2))*diag(t(x)%*%x-t(x)%*%Ez.x%*%t(M)+a*b))
}

likelihoodFA <- function(x,M,psi){
  p = ncol(x)
  n = nrow(x)
  q = ncol(M)
  C = M %*% t(M) + psi
  l = -n/2 *(p*log(2 * pi )+log(det(C)) + tr(solve(C)%*%cov(x)))
  return(l)
}