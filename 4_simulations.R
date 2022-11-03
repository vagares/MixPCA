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




K=3
maxits=100
tol=0.0001
q = 4
nx=4
p=10
n=c(1000,1000,1000)
pi = c(0.2,0.35,0.45) # n'est pas utilisé
mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3)
betat = matrix(rnorm(4*4,mean=0,sd=2),ncol=4)
beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4)
SNR2 = 3
SNR1 = 100

M=100
pike = matrix(0,K,M)
mu1e = matrix(0,p,M)
mu2e = matrix(0,p,M)
mu3e = matrix(0,p,M)
beta1e = matrix(0,nx*q,M)
beta2e = matrix(0,nx*q,M)
beta3e = matrix(0,nx*q,M)
Qe1 = matrix(0,p*q,M)
Qe2 = matrix(0,p*q,M)
Qe3 = matrix(0,p*q,M)
theta2e = numeric(M)
sigma2e = matrix(0,K,M)
err = numeric(M)
nn=sum(n)
for (i in (1:M)){
  data1=data_gen(n,K,q,p,nx,
             pi,
             mu,
             beta,
             SNR1 ,
             SNR2)
  data=data1$data
  est1  =   tryCatch(
    expr  = {est0 = estimates(data,K,maxits,tol, q,nx,p)
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
      list(pik=rep(NA,3),mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G)
    })
  #table(est1$G,data1$data$g)
  if (is.na(est1$pik[1]) ==FALSE){
  n1 = as.numeric(which.max(table(est1$G,data1$data$g)[1,]))
  n2 = as.numeric(which.max(table(est1$G,data1$data$g)[2,]))
  n3 = as.numeric(which.max(table(est1$G,data1$data$g)[3,]))
  err[i] = (max(table(est1$G,data1$data$g)[1,])+max(table(est1$G,data1$data$g)[2,])+max(table(est1$G,data1$data$g)[3,]))/nn}else{
    n1=1;n2=2;n3=3;err[i]=NA
  }
  pike[,i] = est1$pik
  mu1e[,i] = est1$mu[[n1]] - mu[,1]
  mu2e[,i] = est1$mu[[n2]] - mu[,2]
  mu3e[,i] = est1$mu[[n3]] - mu[,3]
  beta1e[,i] = as.vector(est1$beta[[n1]]) - as.vector(betat)
  beta2e[,i] = as.vector(est1$beta[[n2]]) - as.vector(betat)
  beta3e[,i] = as.vector(est1$beta[[n3]]) - as.vector(betat)
  Qe1[,i] = as.vector(est1$Q[[n1]]) -  as.vector(data1$Q[[1]])
  Qe2[,i] = as.vector(est1$Q[[n2]])-  as.vector(data1$Q[[2]])
  Qe3[,i] = as.vector(est1$Q[[n3]])-  as.vector(data1$Q[[3]])
  theta2e[i] = est1$theta2 - data1$theta2
  sigma2e[1,i] = est1$sigma2i[[n1]] - data1$sigma2[[1]]
  sigma2e[2,i] = est1$sigma2i[[n2]] - data1$sigma2[[2]]
  sigma2e[3,i] = est1$sigma2i[[n3]] - data1$sigma2[[3]]
}


dfpi1  =  data.frame(seq(1:M), pike[1,])
dfpi2  =  data.frame(seq(1:M), pike[2,])
dfpi3  =  data.frame(seq(1:M), pike[3,])
dfmu1  =  data.frame(seq(1:M), t(mu1e))
dfmu2  =  data.frame(seq(1:M), t(mu2e))
dfmu3  =  data.frame(seq(1:M), t(mu3e))
dfbeta1  =  data.frame(seq(1:M), t(beta1e))
dfbeta2  =  data.frame(seq(1:M), t(beta2e))
dfbeta3  =  data.frame(seq(1:M), t(beta3e))
dfQ1  =  data.frame(seq(1:M), t(Qe1))
dfQ2  =  data.frame(seq(1:M), t(Qe2))
dfQ3  =  data.frame(seq(1:M), t(Qe3))
dftheta  =  data.frame(seq(1:M), theta2e)
dfsigma1  =  data.frame(seq(1:M),sigma2e[1,])
dfsigma2  =  data.frame(seq(1:M),sigma2e[2,])
dfsigma3  =  data.frame(seq(1:M),sigma2e[3,])
dferr  =  data.frame(seq(1:M),err)
dfpi1$G  =  "G1"
colnames(dfpi1)=  c("n.sim","value","G")
dfpi2$G  =  "G2"
colnames(dfpi2)=  c("n.sim","value","G")
dfpi3$G  =  "G3"
colnames(dfpi3)=  c("n.sim","value","G")
dfpi  =  rbind(dfpi1,dfpi2,dfpi3)
dfpi$coefs="pi"
dfpi$coefs2="pi"
colnames(dfpi)  =  c("n.sim","value","G","coefs","coefs2")
dfmu1$G  =  "G1"
dfmu2$G  =  "G2"
dfmu3$G  =  "G3"
dfmu  =  rbind(dfmu1,dfmu2,dfmu3)
colnames(dfmu)  =  c("n.sim",paste("mu[",1:p,"]",sep=""),"G")
dfmu  =  dfmu %>%
  gather(coefs,  value,   - c("n.sim", "G"))
dfmu$coefs2="mu"
dfmu = dfmu[,c(1,4,2,3,5)]
dfbeta1$G  =  "G1"
dfbeta2$G  =  "G2"
dfbeta3$G  =  "G3"
dfbeta  =  rbind(dfbeta1,dfbeta2,dfbeta3)
colnames(dfbeta)  =  c("n.sim",paste("beta[",1:(nx*q),"]",sep=""),"G")
dfbeta  =  dfbeta %>%
  gather(coefs,  value,   - c("n.sim", "G"))
dfbeta$coefs2="beta"
dfbeta = dfbeta[,c(1,4,2,3,5)]
dfQ1$G  =  "G1"
dfQ2$G  =  "G2"
dfQ3$G  =  "G3"
dfQ  =  rbind(dfQ1,dfQ2,dfQ3)
colnames(dfQ)  =  c("n.sim",paste("Q[",1:(p*q),"]",sep=""),"G")
dfQ  =  dfQ %>%
  gather(coefs,  value,   - c("n.sim", "G"))
dfQ$coefs2="Q"
dfQ = dfQ[,c(1,4,2,3,5)]
dfsigma1$G  =  "G1"
colnames(dfsigma1)=  c("n.sim","value","G")
dfsigma2$G  =  "G2"
colnames(dfsigma2)=  c("n.sim","value","G")
dfsigma3$G  =  "G3"
colnames(dfsigma3)=  c("n.sim","value","G")
dfsigma  =  rbind(dfsigma1,dfsigma2,dfsigma3)
colnames(dfsigma)  =  c("n.sim","value","G")
dfsigma$coefs="sigma"
dfsigma$coefs2="sigma"
colnames(dfsigma)  =  c("n.sim","value","G","coefs","coefs2")
df = rbind(dfpi,dfmu,dfQ,dfbeta,dfsigma)
