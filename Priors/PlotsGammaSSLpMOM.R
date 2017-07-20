#####SSL####
##Including m_ij and sigma^2 as a different variable
dgammaSSL1 <- function(m,sigma2,l0,l1,g.theta){
  ratio <- (1-g.theta)/g.theta * sqrt(l1/l0) * exp(-1/2*(m^2/sigma2)*(1/l0-1/l1))
  density <- 1/(1+ratio)
  #density <- 1/((1+ratio)*sigma2)
  return(density)
}

inflex.gammaSSL <- function(m,sigma2,l0,l1,g.theta){
  a <- (1-g.theta)/g.theta*sqrt(l1/l0)
  b <- (1/sigma2)*(1/l0-1/l1)
  x <- (a*b*m*exp(.5*b*m^2))/(a+exp(.5*b*m^2))^2
  return(x)
}

second.gammaSSL <- function(m,sigma2,l0,l1,g.theta){
  a <- (1-g.theta)/g.theta*sqrt(l1/l0)
  b <- (1/sigma2)*(1/l0-1/l1)
  x <- (a*b*(a*exp(.5*b*m^2)*(b*m^2+1)-exp(b*m^2)*(b*m^2-1)))/(a+exp(.5*b*m^2))^3
  return(x)
}


##Including m_ij and sigma^2 as one variable
dgammaSSL <- function(m,l0,l1,g.theta){
  ratio <- (1-g.theta)/g.theta * sqrt(l1/l0) * exp(-1/2*(m^2)*(1/l0-1/l1))
  density <- 1/(1+ratio)
  return(density)
}

##Obtaining NOT relevant values of m_ij (P=0.05)
mij05 <-function(l0,l1,g.theta){
  ratio.p <- (1-g.theta)/(g.theta)
  ratio.l <- sqrt(l1/l0)
  m<-sqrt((2*l1*l0*log(19/(ratio.p*ratio.l)))/(l0-l1))
  return(m)
}

mij05(0.01,30,.5)
dgammaSSL(mij05(0.01,30,.5),.01,30,.5)

mij95 <-function(l0,l1,g.theta){
  ratio.p <- (1-g.theta)/(g.theta)
  ratio.l <- sqrt(l1/l0)
  m<-sqrt((2*l1*l0*log(1/(19*(ratio.p*ratio.l)))/(l0-l1)))
  return(m)
}

mij95(0.1,1,.5)
dgammaSSL(mij95(0.01,1,.5),.01,1,.5)

#####pMOM####
##Including m_ij and sigma^2 as a different variable
dgammapMOM1 <- function(m,sigma2,l0,l1,g.theta){
  #ratio <- (1-g.theta)/g.theta * sqrt(l1/l0) *sqrt(sigma2*l1)/(m^2)* exp(-1/2*(m^2/sigma2)*(1/l0-1/l1))
  ratio <- (1-g.theta)/g.theta * sqrt(l1/l0) *(sigma2*l1)/(m^2)* exp(-1/2*(m^2/sigma2)*(1/l0-1/l1))
  density <- 1/(1+ratio)
  return(density)
}

##Including m_ij and sigma^2 as one variable
dgammapMOM <- function(m,l0,l1,g.theta){
  ratio <- (1-g.theta)/g.theta * sqrt(l1/l0) *l1/(m^2)* exp(-1/2*(m^2)*(1/l0-1/l1))
  density <- 1/(1+ratio)
  return(density)
}
