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






M=1
q = 4
nx=4
p=10
K=3
n = 10000
pii = c(0.2,0.35,0.45) # n'est pas utilisé
mustar = matrix(c((0:(p-1))^2/(2*p),2*cos((0:(p-1))/2)+1,  1/(1+exp(-(0:(p-1))))),nrow = p,ncol=K)
s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
betat = matrix(rnorm(q*nx,mean=0,sd=2),ncol=nx)
betastar = matrix(rnorm(q*nx,mean=0,sd=2),ncol=nx)
SNR2 = 3
SNR1 = 100
sig2 = .01
tol = 1e-4
maxit = 100

covXpca = matrix(0,p*p,M)
covXest = matrix(0,p*p,M)
covX = matrix(0,p*p,M)
err= numeric(M)
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

for (i in (1:M)){
  data1=data_gen(n=n,K = 3,q = 4,p = 10,nx=4,
                 s = s, 
                 pii = pii,
                 mu = mustar,
                 beta = betastar,
                 SNR1 = SNR1,
                 SNR2 = SNR2)
  X=data1$data
  
  est1  =   tryCatch(
    expr  = {est0 = estimates(X,K=3,maxits=100,
                              tol=tol, 
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
  # est0 = estimates(X,K=3,maxits=100,
  #                  tol=1e-4, 
  #                  q = 4,
  #                  p=10,
  #                  nx=4,
  #                  verbose=TRUE)
  if (is.na(est1$piik[1]) ==FALSE){
    n1 = as.numeric(which.max(table(est1$G,data1$data$g)[1,]))
    n2 = as.numeric(which.max(table(est1$G,data1$data$g)[2,]))
    n3 = as.numeric(which.max(table(est1$G,data1$data$g)[3,]))
    err[i] = (max(table(est1$G,data1$data$g)[1,])+max(table(est1$G,data1$data$g)[2,])+max(table(est1$G,data1$data$g)[3,]))/n}else{
      n1=1;n2=2;n3=3;err[i]=NA
    }
  if (is.na(est1$piik[1]) ==FALSE){
  Xest1 = matrix(0,n,p)
  nn1=length(est1$G[est1$G == 1])
  nn2=length(est1$G[est1$G == 2])
  nn3=length(est1$G[est1$G == 3])
  x1 = data1$data[est1$G == 1,(p+1):(p+nx)]
  alphai1 = t(est1$beta[[1]]%*%t(x1)) + matrix(rnorm(nn1*q,mean=0,sd=sqrt(est1$sigma2)),ncol=q)
  Xest1[which(est1$G == 1),1:p]=rep(1,nn1)%*%t(est1$mu[[1]])+t(est1$Q[[1]]%*%t(alphai1))+matrix(rnorm(p*nn1,mean=0,sd=sqrt(est1$theta2)),ncol=p,nrow=nn1)
  x2 = data1$data[est1$G == 2,(p+1):(p+nx)]
  alphai2 = t(est1$beta[[2]]%*%t(x2)) + matrix(rnorm(nn2*q,mean=0,sd=sqrt(est1$sigma2)),ncol=q)
  Xest1[which(est1$G == 2),1:p]=rep(1,nn2)%*%t(est1$mu[[2]])+t(est1$Q[[2]]%*%t(alphai2))+matrix(rnorm(p*nn2,mean=0,sd=sqrt(est1$theta2)),ncol=p,nrow=nn2)
  x3 = data1$data[est1$G == 3,(p+1):(p+nx)]
  alphai3 = t(est1$beta[[3]]%*%t(x3)) + matrix(rnorm(nn3*q,mean=0,sd=sqrt(est1$sigma2)),ncol=q)
  Xest1[which(est1$G == 3),1:p]=rep(1,nn3)%*%t(est1$mu[[3]])+t(est1$Q[[3]]%*%t(alphai3))+matrix(rnorm(p*nn3,mean=0,sd=sqrt(est1$theta2)),ncol=p,nrow=nn3)
  
  covXest[,i] = as.vector(cov(Xest1))}else{covXest[,i] = rep(NA,p*p)}
  
  if (is.na(est1$piik[1]) ==FALSE){
    Xest1 = matrix(0,n,p)
    for (i in (1:n)){
    x1 = data1$data[i,(p+1):(p+nx)]
    alphai1 = t(est1$beta[[1]]%*%t(x1)) +rnorm(q,mean=0,sd=sqrt(est1$sigma2))
    X11=t(est1$mu[[1]])+t(est1$Q[[1]]%*%t(alphai1))+rnorm(p,mean=0,sd=sqrt(est1$theta2))
    alphai2 = t(est1$beta[[2]]%*%t(x1)) + rnorm(q,mean=0,sd=sqrt(est1$sigma2))
    X12=t(est1$mu[[2]])+t(est1$Q[[2]]%*%t(alphai2))+rnorm(p,mean=0,sd=sqrt(est1$theta2))
    alphai3 = t(est1$beta[[3]]%*%t(x1)) +rnorm(q,mean=0,sd=sqrt(est1$sigma2))
    X13=t(est1$mu[[3]])+t(est1$Q[[3]]%*%t(alphai3))+rnorm(p,mean=0,sd=sqrt(est1$theta2))
    
    Xest1[i,1:p] = est1$tau[i,1] * X11 + est1$tau[i,2] *X12 + est1$tau[i,3] *X13
    }
    covXest[,i] = as.vector(cov(Xest1))}else{covXest[,i] = rep(NA,p*p)}
  #table(est1$G,data1$data$g)
  pike[,i] = est1$piik
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
covXmean = matrix(apply(covX,1,mean),p,p)
covXestmean = matrix(apply(covXest,1,mean),p,p)

sqrt(sum((apply(mu1e,1,mean,na.rm=TRUE))^2))/p
sqrt(sum((apply(mu2e,1,mean,na.rm=TRUE))^2))/p
sqrt(sum((apply(mu3e,1,mean,na.rm=TRUE))^2))/p

sqrt(mean((1/p)*apply(((t(mu1e) - rep(1,100)%*%t(apply(mu1e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 
sqrt(mean((1/p)*apply(((t(mu2e) - rep(1,100)%*%t(apply(mu2e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 
sqrt(mean((1/p)*apply(((t(mu3e) - rep(1,100)%*%t(apply(mu3e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 


sqrt(sum((apply(beta1e,1,mean,na.rm=TRUE))^2))/16
sqrt(sum((apply(beta2e,1,mean,na.rm=TRUE))^2))/16
sqrt(sum((apply(beta3e,1,mean,na.rm=TRUE))^2))/16

sqrt(mean((1/16)*apply(((t(beta1e) - rep(1,100)%*%t(apply(beta1e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 
sqrt(mean((1/16)*apply(((t(beta2e) - rep(1,100)%*%t(apply(beta2e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 
sqrt(mean((1/16)*apply(((t(beta3e) - rep(1,100)%*%t(apply(beta3e,1,mean,na.rm=TRUE))))^2,1,sum,na.rm=TRUE))) 


mean(sigma2e[1,],na.rm=TRUE)
mean(sigma2e[2,],na.rm=TRUE)
mean(sigma2e[3,],na.rm=TRUE)

sqrt(mean((sigma2e[1,] - mean(sigma2e[1,],na.rm=TRUE))^2,na.rm=TRUE))
sqrt(mean((sigma2e[2,] - mean(sigma2e[2,],na.rm=TRUE))^2,na.rm=TRUE))
sqrt(mean((sigma2e[3,] - mean(sigma2e[3,],na.rm=TRUE))^2,na.rm=TRUE))


mean(theta2e,na.rm=TRUE)
sqrt(mean((theta2e - mean(theta2e,na.rm=TRUE))^2,na.rm=TRUE))

x11()
par(col=2)
image(covXmean[,rev(1:p)])
title("Cov empirique")
image(covXestmean[,rev(1:p)])
title("Cov estimée PPCA")



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
