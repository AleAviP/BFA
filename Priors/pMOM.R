library(mombf)
tau <- 1
thseq <- seq(-3,3,length=1000)
plot(thseq,demom(thseq,tau=tau),type='l',ylab='Prior density')
lines(thseq,dimom(thseq,tau=tau),lty=2,col=2)
lines(thseq,dmom(thseq,tau=tau),lty=2,col="turquoise")
