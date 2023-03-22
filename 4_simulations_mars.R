library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
setwd('C:\\Users\\vagares\\Documents\\MixPCA')

source("1_simdata.r")
source("2_estimators.r")
#source("3_parameters.r")


#res  =  runSims1(vList = vlis1(nsim = 1000),  doOne  =  doOne, parallel = TRUE, seedList = NULL)






B=1
q = 4
nx=4
p=10
n = 10000
pii = c(0.2,0.35,0.45) # n'est pas utilisé
mu = matrix(c((0:(p-1))^2/(2*p),2*cos((0:(p-1))/2)+1, rep(1,p)),nrow = p,ncol=K)
s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
betat = matrix(rnorm(nx*q,mean=0,sd=2),ncol=nx)
beta = matrix(rnorm(nx*q,mean=0,sd=2),ncol=nx)
SNR2 = 3
SNR1 = 100
sig2 = .01
tol = 1e-4
maxit = 100
K=3
covXpca = matrix(0,p*p,B)
covXest = matrix(0,p*p,B)
covX = matrix(0,p*p,B)
err= numeric(B)
pike = matrix(0,K,B)
mu1e = matrix(0,p,B)
mu2e = matrix(0,p,B)
mu3e = matrix(0,p,B)
beta1e = matrix(0,nx*q,B)
beta2e = matrix(0,nx*q,B)
beta3e = matrix(0,nx*q,B)
Qe1 = matrix(0,p*q,B)
Qe2 = matrix(0,p*q,B)
Qe3 = matrix(0,p*q,B)
theta2e = numeric(B)
sigma2e = matrix(0,K,B)

for (b in (1:B)){
  data1=data_gen(n,K,q,p,nx,
                 s = s, 
                 pii = pii,
                 mu = mu,
                 beta = beta,
                 SNR1 = SNR1,
                 SNR2 = SNR2)
  X=data1$data
  
  est1  =   tryCatch(
    expr  = {est0 = estimates(X,K,maxits=100,
                              tol=tol, 
                              q,
                              p,
                              nx,
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
  sim = list(datasim = data1,est = est1)
  save(sim, file= paste("simulations0323/SimMars-b=",b,".Rdata",sep=""))
}
