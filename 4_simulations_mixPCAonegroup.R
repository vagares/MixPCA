library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
library(FactoMineR)
setwd('C:\\Users\\vagares\\Documents\\MixPCA')

source("1_simdata.r")
source("2_estimators_mixACP_ungroup.r")
#source("3_parameters.r")


#res  =  runSims1(vList = vlis1(nsim = 1000),  doOne  =  doOne, parallel = TRUE, seedList = NULL)
dec= read.csv("C:\\Users\\vagares\\Documents\\MixPCA\\decathlon.csv",sep=";") 
print(dec)

tab = as.matrix(dec[2:10])

est = estimates(tab,maxits=100,tol=0.01, q = 4,p=9)






M=100
n = 10000
p = 10
q=4
s = c(0.7,-0.4)
mu = c(-3,-3,-3,-1,0,0,1,2,2,2)
sig2 = .01
tol = 1e-4
maxit = 100
covXpca = matrix(0,p*p,M)
covXest = matrix(0,p*p,M)
covX = matrix(0,p*p,M)
#nn=sum(n)
for (i in (1:M)){
  data1=data_gen_mixAcponegroup(n=n,q = 4,p=p,
                                s=s,
                                mu=mu,
                                sig2=sig2)
  X=data1$data
  covX[,i] = as.vector(cov(X))
  # est1  =   tryCatch(
  #   expr  = {est0 = estimates(data,
  #                             maxits=100,
  #                             tol=1e-4, 
  #                             q = 4,
  #                             p=10,
  #                             verbose=TRUE)
  #   est0}, error  =  function(cond) {
  #     mu=list()
  #     mu[[1]] = rep(NA,p)
  #     Q=list()
  #     Q[[1]] = rep(NA,p*q)
  #     theta2 = NA
  #     sigma2i = list()
  #     sigma2i[[1]] = NA
  #     G=rep(NA,n)
  #     list(mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i)
  #   })
  
  pca = PCA(X,nc=q)
  U = pca$svd$V
  K = diag(pca$eig[1:q])
  What = U%*%sqrt(K-sig2*diag(rep(1,q)))
  x = matrix(rnorm(n*q),n,q)
  Xhat = x%*%t(What)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)
  # les vraies valeurs de paramètres sont mises
  covXpca[,i] = as.vector(cov(Xhat))
  ### Si on utilise le code EM 
  #est0 = estimates(X,100,1e-4, q,p)
  est1  =   tryCatch(
    expr  = {est0 = estimates(X,100,1e-4, q,p)
    est0}, error  =  function(cond) {
      mu=list()
      mu[[1]] = rep(NA,p)
      Q=list()
      Q[[1]] = matrix(rep(NA,p*q),p,q)
      theta2 = NA
      sigma2i = NA
      alphai=list()
      alphai[[1]] = matrix(rep(NA,n*q),q,n)
      list(mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i,alphai=alphai)
    })
  Xest0 = t(est1$alphai[[1]])%*%t(est1$Q[[1]])+matrix(rep(est1$mu[[1]],n),n,p,byrow=TRUE)+sqrt(est1$sigma2i)*matrix(rnorm(n*p),n,p)
  covXest[,i] = as.vector(cov(Xest0))
  #table(est1$G,data1$data$g)
 
  
}
covXmean = matrix(apply(covX,1,mean),p,p)
covXpcamean = matrix(apply(covXpca,1,mean),p,p)
covXestmean = matrix(apply(covXest,1,mean),p,p)

image(covXmean[,rev(1:p)])
title("Empirical covariance")
image(covXpcamean[,rev(1:p)])
title("Covariance estimated by PCA")
image(covXestmean[,rev(1:p)])
title("Covariance estimated by PPCA")
