library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
setwd('C:\\Users\\vagares\\Documents\\MixPCA')

source("1_simdata.r")
source("2_estimators_mixACP_ungroup.r")
#source("3_parameters.r")


#res  =  runSims1(vList = vlis1(nsim = 1000),  doOne  =  doOne, parallel = TRUE, seedList = NULL)




K=1
maxits=1000
tol=0.0001
q = 4

p=10
n=1000
mu = c((0:9)^2/20)
sigma2 = 3
SNR2 = 3

M=100
mu1e = matrix(0,p,M)
Qe1 = matrix(0,p*q,M)
theta2e = numeric(M)
sigma2e = matrix(0,K,M)
err = numeric(M)
#nn=sum(n)
for (i in (1:M)){
  data1=data_gen_mixAcponegroup(n,q,p,
             mu,
             SNR2,
             sigma2)
  data=data1$data
  est1  =   tryCatch(
    expr  = {est0 = estimates(data,maxits,tol, q,p)
    est0}, error  =  function(cond) {
      mu=list()
      mu[[1]] = rep(NA,p)
      Q=list()
      Q[[1]] = rep(NA,p*q)
      theta2 = NA
      sigma2i = list()
      sigma2i[[1]] = NA
      G=rep(NA,n)
      list(mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i)
    })
  #table(est1$G,data1$data$g)
 
  mu1e[,i] = est1$mu[[1]] - mu
  Qe1[,i] = as.vector(est1$Q[[1]]) -  as.vector(data1$Q[[1]])
  theta2e[i] = est1$theta2 - data1$theta2
  sigma2e[1,i] = est1$sigma2i[[1]] - data1$sigma2[[1]]
}


dfmu1  =  data.frame(seq(1:M), t(mu1e))
dfQ1  =  data.frame(seq(1:M), t(Qe1))
dftheta  =  data.frame(seq(1:M), theta2e)
dfsigma1  =  data.frame(seq(1:M),sigma2e[1,])


dfmu  =  dfmu1
colnames(dfmu)  =  c("n.sim",paste("mu[",1:p,"]",sep=""))
dfmu  =  dfmu %>%
  gather(coefs,  value,   - c("n.sim"))
dfmu$coefs2="mu"

dfQ  =  dfQ1
colnames(dfQ)  =  c("n.sim",paste("Q[",1:(p*q),"]",sep=""))
dfQ  =  dfQ %>%
  gather(coefs,  value,   - c("n.sim"))
dfQ$coefs2="Q"

colnames(dfsigma1)=  c("n.sim","value")

dfsigma  =  dfsigma1
colnames(dfsigma)  =  c("n.sim","value")
dfsigma$coefs="sigma"
dfsigma$coefs2="sigma"

colnames(dfsigma)  =  c("n.sim","value","coefs","coefs2")
dfsigma=dfsigma[,c(1,3,2,4)]
df = rbind(dfmu,dfQ,dfsigma)
