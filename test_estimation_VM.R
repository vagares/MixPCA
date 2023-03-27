library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
#setwd('C:\\Users\\vagares\\Documents\\MixPCA')
setwd('~/OneDrive - Universit?? de Rennes 1/RECHERCHE/MixPCA-main/')

source("1_simdata.r")
source("2_estimators.r")
#source("3_parameters.r")


M=1 # M=100
K = 3
q = 2
nx = 2
p=10
n = 10000
pii = c(0.2,0.35,0.45) # n'est pas utilis??
mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3)
s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
beta = list()
beta[[1]] = matrix(c(-1,-1,-1,1),q,nx)
beta[[2]] = matrix(c(1,1,0,0),q,nx)
beta[[3]] = matrix(c(-2,0,2,2),q,nx)

SNR2 = 1
SNR1 = 1
sig2 = .01


tol = 1e-4
maxit = 100

data1=data_gen(n=n,K = K,q = q,p = 10,nx=nx,
               s = s, 
               pii = pii,
               mu = mu,
               beta = beta,
               SNR1 = SNR1,
               SNR2 = SNR2)
X=data1$data

### Initialisation
par_init = list()
par_init$piik = pii
par_init$theta2 =  10*data1$theta2
par_init$sigma2 =  2#mean(unlist(data1$sigma2))#unlist(data1$sigma2)
par_init$Q = data1$Q
beta_init = list()
for (k in 1:K){
  beta_init[[k]] = matrix(0,q,nx+1) 
  beta_init[[k]][,2:(nx+1)] = beta[[k]]
}
par_init$beta = beta_init
par_init$mu = list()
for (k in 1:K){par_init$mu[[k]] = mu[,k]}

est0 = estimates(X,
                 K=3,
                 par_init=par_init,
                 maxits=5,
                 tol=tol, 
                 q = q,
                 p=10,
                 nx=nx,
                 verbose=TRUE)

