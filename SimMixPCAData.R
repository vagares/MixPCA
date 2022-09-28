library(pracma)

### Paramètres

K = 3
q = 4
p = 10 # a voir si p et q doivent etre egaux

pik = c(0.2,0.35,0.45)
muk = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = p,ncol=K)
#matplot(muk,type="l")

geom = runif(4,-1,1)

beta = matrix(rnorm(q*4,mean=0,sd=2),ncol=4)

SNR1 = 100 
bkGi = as.vector(beta%*%geom)
sigma2 = (sd(bkGi)/SNR1)^2
alphaik = bkGi + rnorm(q,mean=0,sd=sqrt(sigma2))

#plot(alphaik,type="l")
#lines(bkGi,type="l",col=2)
library(fda)

bfourier = create.fourier.basis(rangeval=c(0,1),nbasis=q)
plot(bfourier)
fmat = getbasismatrix(evalarg=seq(0,1,length=p),bfourier)

#Qk = gramSchmidt(matrix(runif(p*q,-1,1),nrow=p,ncol=q))$Q
Qk = gramSchmidt(fmat[,-1])$Q
SNR2 = 3
theta2 = (mean(apply(muk,2,sd))/SNR2)^2 

matplot(muk+as.vector(Qk%*%alphaik)+matrix(rnorm(p*K,mean=0,sd=sqrt(theta2)),ncol=K,nrow=p),type="l",lty=3)
matlines(muk+as.vector(Qk%*%alphaik),lty=1)
matlines(muk,lty=1)

