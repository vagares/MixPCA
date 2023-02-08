library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
setwd('C:\\Users\\vagares\\Documents\\MixPCA')

source("1_simdata.r")
source("2_estimators_mixACP.r")
#source("3_parameters.r")


#res  =  runSims1(vList = vlis1(nsim = 1000),  doOne  =  doOne, parallel = TRUE, seedList = NULL)




K=3
M=100
n = 10000
p = 10
q = 4
nx=4
s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
pii = c(0.2,0.35,0.45) # n'est pas utilisé
mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3)
sigma2 = 0.05
SNR2 = 3
sig2 = .01
tol = 1e-4
maxit = 100
covXpca = matrix(0,p*p,M)
covXest = matrix(0,p*p,M)
covX = matrix(0,p*p,M)
err= numeric(M)
#nn=sum(n)
for (i in (1:M)){
  data1=data_gen_mixAcp(n=n,K = 3,q = 4,p = 10,
                        s = s, 
                        pii = pii,
                        mu = mu,
                        SNR2 = SNR2,
                        sigma2 = sigma2)
  X=data1$data
  covX[,i] = as.vector(cov(X[1:p]))

  ### Si on utilise le code EM 
  
  est0 = estimates(X,K,100,1e-4, q,p)
  est1  =   tryCatch(
    expr  = {est0 = estimates(X,K,100,1e-4, q,p)
    est0}, error  =  function(cond) {
      mu=list()
      mu[[1]] = rep(NA,p)
      mu[[2]] = rep(NA,p)
      mu[[3]] = rep(NA,p)
      Q=list()
      Q[[1]] = rep(NA,p*q)
      Q[[2]] = rep(NA,p*q)
      Q[[3]] = rep(NA,p*q)
      theta2 = NA
      sigma2i = list()
      sigma2i[[1]] = NA
      sigma2i[[2]] = NA
      sigma2i[[3]] =NA
      G=rep(NA,n)
      alphai=list()
      alphai[[1]] = rep(NA,n*q)
      alphai[[2]] = rep(NA,n*q)
      alphai[[3]] = rep(NA,n*q)
      list(pik=rep(NA,3),mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai)
    })
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
  x1 = data1$data[est0$G == 1,1:p]
  Xest0[which(est0$G == 1),1:p]=t(est0$alphai[[1]][,which(est0$G == 1)])%*%t(est0$Q[[1]])+matrix(rep(est0$mu[[1]],nn1),nn1,p,byrow=TRUE)+sqrt(est0$sigma2i[1])*matrix(rnorm(nn1*p),nn1,p)
  x2 = data1$data[est0$G == 2,1:p]
  Xest0[which(est0$G == 2),1:p]=t(est0$alphai[[2]][,which(est0$G == 2)])%*%t(est0$Q[[2]])+matrix(rep(est0$mu[[2]],nn2),nn2,p,byrow=TRUE)+sqrt(est0$sigma2i[2])*matrix(rnorm(nn2*p),nn2,p)
  x3 = data1$data[est0$G == 3,1:p]
  Xest0[which(est0$G == 3),1:p]=t(est0$alphai[[3]][,which(est0$G == 3)])%*%t(est0$Q[[3]])+matrix(rep(est0$mu[[3]],nn3),nn3,p,byrow=TRUE)+sqrt(est0$sigma2i[3])*matrix(rnorm(nn3*p),nn3,p)
  
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