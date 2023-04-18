library(simsalapar)
library(tibble)
library(nlme)
library(dplyr)
library(tidyr)
library(tidyverse)
#setwd('C:\\Users\\vagares\\Documents\\MixPCA')
setwd('~/OneDrive - Université de Rennes 1/RECHERCHE/MixPCA-main/')

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
beta[[3]] = 0*matrix(c(-2,0,2,2),q,nx)
beta[[3]] = matrix(c(-2,0,2,2),q,nx)

SNR2 = 2
SNR1 = 20#.5


tol = 1e-4
maxit = 100

set.seed("123")
data1=data_gen(n=n,
               K = K,
               q = q,
               p = 10,
               nx=nx,
               s = s, 
               pii = pii,
               mu = mu,
               beta = beta,
               SNR1 = SNR1,
               SNR2 = SNR2)
X=data1$data
matplot(t(data1$data[,1:p]),col=data1$data$g,type="l",lty=1)

est0 = estimates(X,
                 K=3,
                 par_init=NULL,
                 maxits=30,
                 tol=tol, 
                 q = q,
                 p=10,
                 nx=nx,
                 verbose=TRUE)

### Validation par simulation =========================
library(fields)

B = 1
Xest0 = array(0,c(n,p,B))

x = X[,(p+1):(p+nx)]
x1 = matrix(1,dim(x)[1],nx+1)
x1[,2:(nx+1)] = as.matrix(x)

# Dans les simulations, la loi de x ne dépend pas des classes

for (b in 1:B){
  g = numeric(n)
  u=runif(n)
  piic=cumsum(est0$piik)
  for (i in (1:n)){
    if (u[i]<piic[1]){g[i]=1}else
      for (k in (1:(K-1))){
        if (piic[k]<u[i] & u[i]<piic[k+1]){g[i]=k+1}}
  }
  for (k in (1:K)){
    nk=sum(as.numeric(g==k)) 
    betak=est0$beta[[k]] 
    bkGi = t(betak%*%t(x1[g==k,])) # comment choisir xk???
    #theta2 = (mean(apply(bkGi,1,sd))/SNR1)^2
    #if (min(theta2)>0) {theta22[[k]] = theta2
    #} else {theta22[[k]] = sigma2}
    piik = pii[k]
    muk = mu[,k]
    alpha_sim = bkGi + matrix(rnorm(nk*q,0,sd=sqrt(est0$theta2)),nk,q)
    Xest0[which(g == k),1:p,b]=alpha_sim%*%t(est0$Q[[k]])+matrix(rep(est0$mu[[k]],nk),nk,p,byrow=TRUE)+sqrt(est0$sigma2i[k]/p)*matrix(rnorm(nk*p),nk,p)
  }
  
}

# comparaison des covariances (pour un tirage)
par(mfrow=c(1,2))
CX = cov(X[,1:p])
CXsim = cov(Xest0[,,1])
image.plot(CX,zlim=range(CX,CXsim))
title("Observed")
image.plot(CXsim,zlim=range(CX,CXsim))
title("Simulated")



par(mfrow=c(3,2),mar=c(3,3,3,3))
for (j in c(1,4,10)){
  hist(X[,j],main=paste("Obs., var",j))
  hist(Xest0[,j,b],main=paste("Sim., var",j))
}

# paramètres estimés : Q%$%beta
# attention il faudrait mettre les classes dans l'ordre
par(mfrow=c(3,2))
jmi = which.min(est0$piik)
jma = which.max(est0$piik)
j = setdiff(1:3,c(jmi,jma))
P = data1$Q[[1]]%*%beta[[1]]
plot(P,pch=20)
plot(est0$Q[[jmi]]%*%est0$beta[[jmi]][,2:3],pch=20,
     xlim = range(P[,1]),ylim=range(P[,2]))

P = data1$Q[[2]]%*%beta[[2]]
plot(P,pch=20)
plot(est0$Q[[j]]%*%est0$beta[[j]][,2:3],pch=20,
     xlim = range(P[,1]),ylim=range(P[,2]))


P = data1$Q[[3]]%*%beta[[3]]
plot(P,pch=20)
plot(est0$Q[[jma]]%*%est0$beta[[jma]][,2:3],pch=20,
     xlim = range(P[,1]),ylim=range(P[,2]))




N = dim(X)[1]
y = as.matrix(X[,1:p],N,p)
G = X$g


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




