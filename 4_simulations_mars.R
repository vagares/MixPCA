# library(simsalapar)
# library(tibble)
# library(nlme)
# library(dplyr)
# library(tidyr)
# library(tidyverse)
# setwd('C:\\Users\\vagares\\Documents\\MixPCA')

source("1_simdata.r")
source("2_estimators.r")


B   = 100 # nombre de repetitions
q   = 3   # nombre d'axes principaux 
nx  = 4   # nombre de covariables
p   = 10
K   = 3
pii = c(0.2,0.35,0.45)

mu        = matrix(c((0:(p-1))^2/(2*p),2*cos((0:(p-1))/2)+1, (1+exp(-(0:(p-1))))),nrow = p,ncol=K)
s         = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
beta      = list()
beta[[1]] = matrix(c(-1,-1,-1,1),q,nx)
beta[[2]] = matrix(c(1,1,0,0),q,nx)
beta[[3]] = matrix(c(-2,0,2,2),q,nx)

# parametres a faire varier 
nind = c(500,1000,10000)
SNR2 = c(3,10)
SNR1 = c(3,10)

tol = 1e-4
maxit = 100


for (b in (56:B)){

  for (n in nind){
    snr1 = 5
    snr2 = 5
    set.seed(b*n*snr1*snr2)

    data1=data_gen(n,K,q,p,nx,
                   s = s,
                   pii = pii,
                   mu = mu,
                   beta = beta,
                   SNR1 = snr1,
                   SNR2 = snr2)

    test.ind = sample(n,size=round(0.2*n))
    data1$data$app = rep(1,n)
    data1$data$app[test.ind] = 0
    #save(data1, file=paste("simulations0323/donneesSim/donnees-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    X=data1$data[data1$data$app==1,]

    est1  =   tryCatch(
      expr  = {est0 = estimates(X,K,par_init = NULL,maxits=100,
                                tol=tol,
                                q,
                                p,
                                nx,
                                cste = FALSE,
                                verbose=TRUE)
      est0}, error  =  function(cond) {
        mu=list()
        mu[[1]] = rep(NA,p)
        mu[[2]] = rep(NA,p)
        mu[[3]] = rep(NA,p)
        beta=list()
        beta[[1]] = matrix(rep(NA,nx*q),nx,q)
        beta[[2]] = matrix(rep(NA,nx*q),nx,q)
        beta[[3]] = matrix(rep(NA,nx*q),nx,q)
        Q=list()
        Q[[1]] = matrix(rep(NA,p*q),p,q)
        Q[[2]] = matrix(rep(NA,p*q),p,q)
        Q[[3]] = matrix(rep(NA,p*q),p,q)
        theta2 = NA
        sigma2i = list()
        sigma2i[[1]] = NA
        sigma2i[[2]] = NA
        sigma2i[[3]] =NA
        G=rep(NA,n)
        alphai = list()
        alphai[[1]] = matrix(rep(NA,n*q),q,n)
        alphai[[2]] = matrix(rep(NA,n*q),q,n)
        alphai[[3]] = matrix(rep(NA,n*q),q,n)
        tau = matrix(rep(NA,n*3),n,3)
        list(piik=rep(NA,3),mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai,tau=tau)
      })
    save(est1, file=paste("simulations0323/results/results-sanscste-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
  }
  
  for (snr1 in SNR1){
    snr2 = 5
    n    = 1000
    set.seed(b*n*snr1*snr2)
    
    data1=data_gen(n,K,q,p,nx,
                   s = s, 
                   pii = pii,
                   mu = mu,
                   beta = beta,
                   SNR1 = snr1,
                   SNR2 = snr1)
    
    test.ind = sample(n,size=round(0.2*n))
    data1$data$app = rep(1,n)
    data1$data$app[test.ind] = 0
    #save(data1, file=paste("simulations0323/donneesSim/donnees-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    X=data1$data[data1$data$app==1,]
    
    est1  =   tryCatch(
      expr  = {est0 = estimates(X,K,par_init = NULL,maxits=100,
                                tol=tol, 
                                q,
                                p,
                                nx,
                                cste=FALSE,
                                verbose=TRUE)
      est0}, error  =  function(cond) {
        mu=list()
        mu[[1]] = rep(NA,p)
        mu[[2]] = rep(NA,p)
        mu[[3]] = rep(NA,p)
        beta=list()
        beta[[1]] = matrix(rep(NA,nx*q),nx,q)
        beta[[2]] = matrix(rep(NA,nx*q),nx,q)
        beta[[3]] = matrix(rep(NA,nx*q),nx,q)
        Q=list()
        Q[[1]] = matrix(rep(NA,p*q),p,q)
        Q[[2]] = matrix(rep(NA,p*q),p,q)
        Q[[3]] = matrix(rep(NA,p*q),p,q)
        theta2 = NA
        sigma2i = list()
        sigma2i[[1]] = NA
        sigma2i[[2]] = NA
        sigma2i[[3]] =NA
        G=rep(NA,n)
        alphai = list()
        alphai[[1]] = matrix(rep(NA,n*q),q,n)
        alphai[[2]] = matrix(rep(NA,n*q),q,n)
        alphai[[3]] = matrix(rep(NA,n*q),q,n)
        tau = matrix(rep(NA,n*3),n,3)
        list(piik=rep(NA,3),mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai,tau=tau)
      })
    save(est1, file=paste("simulations0323/results/results-sanscste-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
  }
  
  for (snr2 in SNR2){
    snr1 = 5
    n    = 1000
    set.seed(b*n*snr1*snr2)
    
    data1=data_gen(n,K,q,p,nx,
                   s = s, 
                   pii = pii,
                   mu = mu,
                   beta = beta,
                   SNR1 = snr1,
                   SNR2 = snr1)
    
    test.ind = sample(n,size=round(0.2*n))
    data1$data$app = rep(1,n)
    data1$data$app[test.ind] = 0
    #save(data1, file=paste("simulations0323/donneesSim/donnees-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    X=data1$data[data1$data$app==1,]
    
    est1  =   tryCatch(
      expr  = {est0 = estimates(X,K,par_init = NULL,maxits=100,
                                tol=tol, 
                                q,
                                p,
                                nx,
                                cste=FALSE,
                                verbose=TRUE)
      est0}, error  =  function(cond) {
        mu=list()
        mu[[1]] = rep(NA,p)
        mu[[2]] = rep(NA,p)
        mu[[3]] = rep(NA,p)
        beta=list()
        beta[[1]] = matrix(rep(NA,nx*q),nx,q)
        beta[[2]] = matrix(rep(NA,nx*q),nx,q)
        beta[[3]] = matrix(rep(NA,nx*q),nx,q)
        Q=list()
        Q[[1]] = matrix(rep(NA,p*q),p,q)
        Q[[2]] = matrix(rep(NA,p*q),p,q)
        Q[[3]] = matrix(rep(NA,p*q),p,q)
        theta2 = NA
        sigma2i = list()
        sigma2i[[1]] = NA
        sigma2i[[2]] = NA
        sigma2i[[3]] =NA
        G=rep(NA,n)
        alphai = list()
        alphai[[1]] = matrix(rep(NA,n*q),q,n)
        alphai[[2]] = matrix(rep(NA,n*q),q,n)
        alphai[[3]] = matrix(rep(NA,n*q),q,n)
        tau = matrix(rep(NA,n*3),n,3)
        list(piik=rep(NA,3),mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai,tau=tau)
      })
    save(est1, file=paste("simulations0323/results/results-sanscste-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
  }
}

    
    
    
 

