library(LaplacesDemon)
#####SSL####
#Solving the equation system
l0=3429/(34295*log(19))
l1=3429/(5*log(19))
g.theta=0.5
seq <- seq(-2,2,length=1000)
gamma1<-dgammaSSL(seq,l0,l1,g.theta)
plot(seq, gamma1, type="l", xlab="m_ij", ylab="probability",
     main="Gamma density", col="turquoise4")
abline(h=.95,lty=3,col="indianred3")
abline(h=.05,lty=3,col="indianred3")
abline(v=sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.1),lty=2,col="darkblue")
abline(v=sqrt(.1),lty=2,col="darkblue")

#Other inclusion probability
l1_2=3429/(10*log(19))
l0_2=3429/(68590*log(19))
gamma2<-dgammaSSL(seq,l0_2,l1_2,g.theta)
lines(seq,gamma2,col="mediumorchid4")
abline(v=sqrt(.05),lty=2,col="darkblue")
abline(v=sqrt(.25),lty=2,col="darkblue")
abline(v=-sqrt(.05),lty=2,col="darkblue")
abline(v=-sqrt(.25),lty=2,col="darkblue")

#####pMOM####
#Matching non-important values
gammapMOM2<-dgammapMOM(seq,0.00490529,l1,.5)
plot(seq, gammapMOM2, type="l", xlab="m_ij", ylab="probability",
     main="Gamma",col="turquoise4")
lines(seq,gamma1,lty=1,col="mediumorchid4")
abline(h=.95,lty=2,col="indianred3")
abline(h=.05,lty=2,col="indianred3")
abline(v=sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.1),lty=2,col="darkblue")
abline(v=sqrt(.1),lty=2,col="darkblue")

#KLD distance
plot(seq,KLD(gammapMOM2,gamma2)$KLD.px.py, type="l", col="darkblue")

#Matching important values
gammapMOM3<-dgammapMOM(seq,0.0180885,l1,.5)
plot(seq, gammapMOM3, type="l", xlab="m_ij", ylab="probability",
     main="Gamma",col="turquoise4")
lines(seq,gamma1,lty=1,col="mediumorchid4")
abline(h=.95,lty=2,col="indianred3")
abline(h=.05,lty=2,col="indianred3")
abline(v=sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.5),lty=2,col="darkblue")
abline(v=-sqrt(.1),lty=2,col="darkblue")
abline(v=sqrt(.1),lty=2,col="darkblue")

plot(seq,KLD(gammapMOM3,gamma2)$KLD.px.py, type="l", col="darkblue")
