library(mvtnorm) 
library(fda)
library(pracma)
#########################################################################

data = function(n=c(100,100,100),K = 3,q = 4,p = 10,nx=4,
                pi = c(0.2,0.35,0.45),
                mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3),
                beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4),
                SNR1 = 100,
                SNR2 = 3){
  
  # a voir si p et q doivent etre egaux
  
  theta2 = (mean(apply(mu,2,sd))/SNR2)^2
  y = NULL
  x = NULL
  for (k in (1:K)){
    nk=n[k]
    betak=beta
    xk = matrix(runif(nx*nk,-1,1),ncol=nx)
    bkGi = t(betak%*%t(xk))
    sigma2 = (sd(bkGi)/SNR1)^2
    pik = pi[k]
    muk = mu[,k]
    alphaik = bkGi + matrix(rnorm(nk*q,mean=0,sd=sqrt(sigma2)),ncol=q)
    yk = matrix(p*nk,ncol=p)
    bfourier = create.fourier.basis(rangeval=c(0,1),nbasis=q)
    #plot(bfourier)
    fmat = getbasismatrix(evalarg=seq(0,1,length=p),bfourier)
    
    Qk = gramSchmidt(fmat[,-1])$Q
    
    yk=rep(1,nk)%*%t(muk)+t(Qk%*%t(alphaik))+matrix(rnorm(p*nk,mean=0,sd=sqrt(theta2)),ncol=p,nrow=nk)
    y=rbind(y,yk)
    x=rbind(x,xk)
    #matplot(muk+as.vector(Qk%*%alphaik)+matrix(rnorm(p*K,mean=0,sd=sqrt(theta2)),ncol=K,nrow=p),type="l",lty=3)
    #matlines(muk+as.vector(Qk%*%alphaik),lty=1)
    #matlines(muk,lty=1)
  }
  return(data)
}

data1=data(n=c(100,100,100),K = 3,q = 4,p = 10,nx=4,
                pi = c(0.2,0.35,0.45),
                mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3),
                beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4),
                SNR1 = 100,
                SNR2 = 3)

matplot(t(data1[,1:p]),col=data1$g,type="l")
