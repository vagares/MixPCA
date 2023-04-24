library(fields)

B   = 100   # nombre de repetitions
q   = 3   # nombre d'axes principaux 
nx  = 4   # nombre de covariables
p   = 10
K   = 3
pii = c(0.2,0.35,0.45)

mu        = matrix(c((0:(p-1))^2/(2*p),2*cos((0:(p-1))/2)+1, rep(1,p)),nrow = p,ncol=K)
s         = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
beta      = list()
beta[[1]] = matrix(c(-1,-1,-1,1),q,nx)
beta[[2]] = matrix(c(1,1,0,0),q,nx)
beta[[3]] = matrix(c(-2,0,2,2),q,nx)

nind = c(500,1000,10000)
SNR2 = c(3,5,10)
SNR1 = c(3,5,10)

############## lecture resultats


res.simu = rep(0,24)

for (b in (1:100)){
  for (n in nind){
    snr1 = 5
    snr2 = 5
    load(paste("simulations0323/results/results-n=",n,"-snr1=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    load(paste("simulations0323/donneesSim/donnees-n=",n,"-snr=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    
    llik.dec = any(diff(est1$llik)>0)
    piik     = sort(est1$piik)
    ordregpe = order(est1$piik)
    
    theta2.hat = est1$theta2
    theta2.sim = data1$theta2
    
    sigma2.hat = sum(est1$sigma2i)/p
    sigma2.sim = data1$sigma2
    
    ## MSE sur mu et sur Q, beta, et leur produit
    mu.hat = est1$mu[ordregpe]
    MSE.mu = sapply(1:K,FUN=function(k){mean((mu.hat[[k]] - mu[,k])^2)} )
    
    Q.hat = est1$Q[ordregpe]
    MSE.Q = sapply(1:K,FUN=function(k){mean((Q.hat[[k]] - data1$Q[[k]])^2)} )
    
    beta.hat = est1$beta[ordregpe]
    MSE.beta = sapply(1:K,FUN=function(k){mean((beta.hat[[k]][,-1] - beta[[k]])^2)} )
    
    MSE.betaQ = sapply(1:K,FUN=function(k){mean((Q.hat[[k]]%*%beta.hat[[k]][,-1] - data1$Q[[k]]%*%beta[[k]])^2)} )
    
    res.simu = rbind(res.simu,c(n,snr1,snr2,b,llik.dec,piik,theta2.sim,theta2.hat,sigma2.sim,sigma2.hat,MSE.mu,MSE.Q,MSE.beta,MSE.betaQ))
  }
  
  for (snr1 in SNR1){
    snr2 = 5
    n    = 1000
    load(paste("simulations0323/donneesSim/donnees-n=",n,"-snr=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    load(paste("simulations0323/results/results-n=",n,"-snr=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    
    llik.dec = any(diff(est1$llik)>0)
    piik     = sort(est1$piik)
    ordregpe = order(est1$piik)
    
    theta2.hat = est1$theta2
    theta2.sim = data1$theta2
    
    sigma2.hat = sum(est1$sigma2i)/p
    sigma2.sim = data1$sigma2
    
    ## MSE sur mu et sur Q, beta, et leur produit
    mu.hat = est1$mu[ordregpe]
    MSE.mu = sapply(1:K,FUN=function(k){mean((mu.hat[[k]] - mu[,k])^2)} )
    
    Q.hat = est1$Q[ordregpe]
    MSE.Q = sapply(1:K,FUN=function(k){mean((Q.hat[[k]] - data1$Q[[k]])^2)} )
    
    beta.hat = est1$beta[ordregpe]
    MSE.beta = sapply(1:K,FUN=function(k){mean((beta.hat[[k]][,-1] - beta[[k]])^2)} )
    
    MSE.betaQ = sapply(1:K,FUN=function(k){mean((Q.hat[[k]]%*%beta.hat[[k]][,-1] - data1$Q[[k]]%*%beta[[k]])^2)} )
    
    res.simu = rbind(res.simu,c(n,snr1,snr2,b,llik.dec,piik,theta2.sim,theta2.hat,sigma2.sim,sigma2.hat,MSE.mu,MSE.Q,MSE.beta,MSE.betaQ))
  }
  
  for (snr2 in SNR2){
    snr1 = 5
    n    = 1000
    load(paste("simulations0323/donneesSim/donnees-n=",n,"-snr=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))
    load(paste("simulations0323/results/results-n=",n,"-snr=",snr1,"-snr2=",snr2,"-b=",b,".Rdata",sep=""))

    llik.dec = any(diff(est1$llik)>0)
    piik     = sort(est1$piik)
    ordregpe = order(est1$piik)
    
    theta2.hat = est1$theta2
    theta2.sim = data1$theta2
    
    sigma2.hat = sum(est1$sigma2i)/p
    sigma2.sim = data1$sigma2
    
    ## MSE sur mu et sur Q, beta, et leur produit
    mu.hat = est1$mu[ordregpe]
    MSE.mu = sapply(1:K,FUN=function(k){mean((mu.hat[[k]] - mu[,k])^2)} )
    
    Q.hat = est1$Q[ordregpe]
    MSE.Q = sapply(1:K,FUN=function(k){mean((Q.hat[[k]] - data1$Q[[k]])^2)} )
    
    beta.hat = est1$beta[ordregpe]
    MSE.beta = sapply(1:K,FUN=function(k){mean((beta.hat[[k]][,-1] - beta[[k]])^2)} )
    
    MSE.betaQ = sapply(1:K,FUN=function(k){mean((Q.hat[[k]]%*%beta.hat[[k]][,-1] - data1$Q[[k]]%*%beta[[k]])^2)} )
    
    res.simu = rbind(res.simu,c(n,snr1,snr2,b,llik.dec,piik,theta2.sim,theta2.hat,sigma2.sim,sigma2.hat,MSE.mu,MSE.Q,MSE.beta,MSE.betaQ))
  }
}

res.simu = res.simu[-1,]
colnames(res.simu) = c("n","SNR1","SNR2","b","llik.dec",paste("piik",1:3,sep=""),"theta2.sim","theta2.hat","sigma2.sim","sigma2.hat",paste("MSE.mu",1:3,sep=""),paste("MSE.Q",1:3,sep=""),paste("MSE.beta",1:3,sep=""),paste("MSE.betaQ",1:3,sep=""))

write.table(res.simu,file="simulations0323/analyse/ResultsSummary.txt",row.names=FALSE)
  


    
    

B=50
Xest0 = array(0,c(n,p,B))

for (b in 1:B){
  load(paste("simulations0323/SimMars-i=",b,".Rdata",sep=""))
  est0=sim$est
  X = sim$datasim$data
  g = numeric(n)
  u=runif(n)
  piic=cumsum(est0$piik)
  for (i in (1:n)){
    if (u[i]<piic[1]){g[i]=1}else
      for (k in (1:(K-1))){if (piic[k]<u[i] & u[i]<piic[k+1]){g[i]=k+1}}
  }
  for (k in (1:K)){
    nk=sum(as.numeric(g==k)) 
    alpha_sim = matrix(rnorm(nk*4,0,sd=sqrt(est0$theta2)),nk,4)
    Xest0[which(g == k),1:p,b]=alpha_sim%*%t(est0$Q[[k]])+matrix(rep(est0$mu[[k]],nk),nk,p,byrow=TRUE)+sqrt(est0$sigma2i[k]/p)*matrix(rnorm(nk*p),nk,p)
  }
  
}

# comparaison des covariances (pour un tirage)
par(mfrow=c(1,2))

CX = cov(X[,-11])
CXsim = cov(Xest0[,,1])
image.plot(CX,zlim=range(CX,CXsim))
title("Observed")
image.plot(CXsim,zlim=range(CX,CXsim))
title("OSimulated")



par(mfrow=c(3,2),mar=c(3,3,3,3))
for (j in c(1,4,10)){
  hist(X[,j],main=paste("Obs., var",j))
  hist(Xest0[,j,b],main=paste("Sim., var",j))
}


# lois marginales
par(mfrow=c(1,3),mar=c(3,3,3,3))
for (j in c(1,4,10)){
qqp = qqplot(X[,j],Xest0[,j,1],pch=20,col="white",xlab='Observations',ylab='Simulations',cex=.6)
abline(a=0,b=1)
q = matrix(0,B,length(X[,j]))
s.data = sort(X[,j])
lens = length(s.data)/10
for (b in 1:B) {
    q[b,] = sort(Xest0[,j,b])
}
for (b in 1:B){lines(s.data[seq(1,10000,10)],approx(1:length(X[,1]),q[b,], n = lens)$y,col="gray")}
}

  