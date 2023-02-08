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







q = 4
nx=4
p=10
n = 10000
pii = c(0.2,0.35,0.45) # n'est pas utilisé
mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3)
s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
betat = matrix(rnorm(4*4,mean=0,sd=2),ncol=4)
beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4)
SNR2 = 3
SNR1 = 100
sig2 = .01
tol = 1e-4
maxit = 100

covXpca = matrix(0,p*p,M)
covXest = matrix(0,p*p,M)
covX = matrix(0,p*p,M)
err= numeric(M)


for (i in (1:M)){
  data1=data_gen(n=n,K = 3,q = 4,p = 10,nx=4,
                 s = s, 
                 pii = pii,
                 mu = mu,
                 beta = beta,
                 SNR1 = SNR1,
                 SNR2 = SNR2)
  X=data1$data
  covX[,i] = as.vector(cov(X[1:p]))
  
  est1  =   tryCatch(
    expr  = {est0 = estimates(X,K=3,maxits=100,
                              tol=1e-4, 
                              q = 4,
                              p=10,
                              nx=4,
                              verbose=TRUE)
    est0}, error  =  function(cond) {
      mu=list()
      mu[[1]] = rep(NA,p)
      mu[[2]] = rep(NA,p)
      mu[[3]] = rep(NA,p)
      beta=list()
      beta[[1]] = rep(NA,nx*q)
      beta[[2]] = rep(NA,nx*q)
      beta[[3]] = rep(NA,nx*q)
      Q=list()
      Q[[1]] = rep(NA,p*q)
      Q[[2]] = rep(NA,p*q)
      Q[[3]] = rep(NA,p*q)
      theta2 = NA
      sigma2i = list()
      sigma2i[[1]] = NA
      sigma2i[[2]] = NA
      sigma2i[[3]] =NA
      G=rep(NA,nn)
      Q[[1]] = rep(NA,p*q)
      Q[[2]] = rep(NA,p*q)
      Q[[3]] = rep(NA,p*q)
      alphai[[1]] = rep(NA,p*q)
      alphai[[2]] = rep(NA,p*q)
      alphai[[3]] = rep(NA,p*q)
      list(pik=rep(NA,3),mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai)
    })
  est0 = estimates(X,K=3,maxits=100,
                   tol=1e-4, 
                   q = 4,
                   p=10,
                   nx=4,
                   verbose=TRUE)
  if (is.na(est0$piik[1]) ==FALSE){
    n1 = as.numeric(which.max(table(est0$G,data1$data$g)[1,]))
    n2 = as.numeric(which.max(table(est0$G,data1$data$g)[2,]))
    n3 = as.numeric(which.max(table(est0$G,data1$data$g)[3,]))
    err[i] = (max(table(est0$G,data1$data$g)[1,])+max(table(est0$G,data1$data$g)[2,])+max(table(est0$G,data1$data$g)[3,]))/n}else{
      n1=1;n2=2;n3=3;err[i]=NA
    }
  Xest0 = matrix(0,n,p)
  nn1=length(est0$G[est0$G == 1])
  nn2=length(est0$G[est0$G == 2])
  nn3=length(est0$G[est0$G == 3])
  x1 = data1$data[est0$G == 1,(p+1):(p+nx)]
  alphai1 = t(est0$beta[[1]]%*%t(x1)) + matrix(rnorm(nn1*q,mean=0,sd=sqrt(est0$sigma2)),ncol=q)
  Xest0[which(est0$G == 1),1:p]=rep(1,nn1)%*%t(est0$mu[[1]])+t(est0$Q[[1]]%*%t(alphai1))+matrix(rnorm(p*nn1,mean=0,sd=sqrt(est0$theta2)),ncol=p,nrow=nn1)
  x2 = data1$data[est0$G == 2,(p+1):(p+nx)]
  alphai2 = t(est0$beta[[2]]%*%t(x2)) + matrix(rnorm(nn2*q,mean=0,sd=sqrt(est0$sigma2)),ncol=q)
  Xest0[est0$G == 2,]=rep(1,nn2)%*%t(est0$mu[[2]])+t(est0$Q[[2]]%*%t(alphai2))+matrix(rnorm(p*nn2,mean=0,sd=sqrt(est0$theta2)),ncol=p,nrow=nn2)
  x3 = data1$data[est0$G == 3,(p+1):(p+nx)]
  alphai3 = t(est0$beta[[3]]%*%t(x3)) + matrix(rnorm(nn3*q,mean=0,sd=sqrt(est0$sigma2)),ncol=q)
  Xest0[est0$G == 3,]=rep(1,nn3)%*%t(est0$mu[[3]])+t(est0$Q[[3]]%*%t(alphai3))+matrix(rnorm(p*nn3,mean=0,sd=sqrt(est0$theta2)),ncol=p,nrow=nn3)
  
  covXest[,i] = as.vector(cov(Xest0))
  #table(est1$G,data1$data$g)
  
  
}
covXmean = matrix(apply(covX,1,mean),p,p)
covXpcamean = matrix(apply(covXpca,1,mean),p,p)
covXestmean = matrix(apply(covXest,1,mean),p,p)

image(covXmean[,rev(1:p)])
title("Cov empirique")
image(covXpcamean[,rev(1:p)])
title("Cov estimée PCA")
image(covXestmean[,rev(1:p)])
title("Cov estimée PPCA")