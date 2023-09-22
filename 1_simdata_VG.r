library(mvtnorm) 
library(fda)
library(pracma)
library(FactoMineR)
#########################################################################

data_gen = function(n=1000,K = 3,q = 3,p = 10,nx=4,
                    s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2),
                    pii = c(0.2,0.35,0.45),
                    mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, 1/(1+exp(-(0:9)))),nrow = 10,ncol=3),
                    beta,
                    SNR1 = 5,
                    SNR2 = 3){
  
  sigma2 = (mean(apply(mu,2,sd))/SNR2)^2
  y = NULL
  x = NULL
  g = NULL
  Q = list()
  theta22 = list()
  
  g = sample(1:K,size=n,replace=TRUE,prob=pii)
  
  for (k in (1:K)){
    nk=sum(as.numeric(g==k))
    #betak=beta
    xk = matrix(runif(nx*nk,-1,1),ncol=nx)
    betak=beta[[k]] 
    bkGi = t(betak%*%t(xk))
    theta2 = (mean(apply(bkGi,1,sd))/SNR1)^2
    if (min(theta2)>0) {theta22[[k]] = theta2
    } else {theta22[[k]] = sigma2}
  }
  theta2 = mean(sqrt(unlist(theta22)))^2
  
  
  gg=NULL
  for (k in (1:K)){
    nk=sum(as.numeric(g==k))
    xk = matrix(runif(nx*nk,-1,1),ncol=nx)
    betak=beta[[k]] 
    bkGi = t(betak%*%t(xk))
    piik = pii[k]
    muk = mu[,k]
    alphaik = bkGi + matrix(rnorm(nk*q,mean=0,sd=sqrt(theta2)),ncol=q)
    Sigma = diag(rep(1,p))
    for (j in 1:(p-1)){
      Sigma[j,j+1] = s[1,k]
      Sigma[j+1,j] = s[1,k]
    }
    for (j in 1:(p-2)){
      Sigma[j,j+2] = s[2,k]
      Sigma[j+2,j] = s[2,k]
    }
    Qk = eigen(Sigma)$vectors[,1:q]  # Est-ce que ça suffit?
    xxk=rep(1,nk)%*%t(muk)+t(Qk%*%t(alphaik))+matrix(rnorm(p*nk,mean=0,sd=sqrt(sigma2)),ncol=p,nrow=nk)
    pca = PCA(xxk-rep(1,nk)%*%t(muk),scale.unit=FALSE,nc=q,graph=FALSE)
    U = pca$svd$V
    K = diag(pca$eig[1:q])
    Qk = (U%*%sqrt(K-sigma2*diag(rep(1,q))))#/sqrt(theta2)
    
    Q[[k]]=Qk
    yk=rep(1,nk)%*%t(muk)+t(Qk%*%t(alphaik))+matrix(rnorm(p*nk,mean=0,sd=sqrt(sigma2)),ncol=p,nrow=nk)
    y=rbind(y,yk)
    x=rbind(x,xk)
    gg=c(gg,rep(k,nk))
    #matplot(muk+as.vector(Qk%*%alphaik)+matrix(rnorm(p*K,mean=0,sd=sqrt(theta2)),ncol=K,nrow=p),type="l",lty=3)
    #matlines(muk+as.vector(Qk%*%alphaik),lty=1)
    #matlines(muk,lty=1)
  }
  data = data.frame(y = y,x=x,g=gg)
  colnames(data) = c(paste("y",1:p,sep=""), paste("x",1:nx,sep=""),"g")
  return(list(data = data,Q=Q,sigma2=sigma2, theta2=theta2))
}

# data1=data_gen(n=1000,K = 3,q = 4,p = 10,nx=4,
                # s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2), 
                # pii = c(0.2,0.35,0.45),
                # mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, 1/(1+exp(-(0:9)))),nrow = 10,ncol=3),
                # beta = matrix(rnorm(4*4,mean=0,sd=2),ncol=4),
                # SNR1 = 100,
                # SNR2 = 3)

#p=10
#matplot(t(data1$data[,1:p]),col=data1$data$g,type="l")


data_gen_mixAcp = function(n=1000,K = 3,q = 4,p = 10,
                    s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2), 
                    pii = c(0.2,0.35,0.45),
                    mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3),
                    SNR2 = 3,
                    sigma2 = 0.01){
  
  # a voir si p et q doivent etre egaux
  
  theta2 = (mean(apply(mu,2,sd))/SNR2)^2
  y = NULL
  x = NULL
  g = NULL
  Q = list()
  sigma22 = list()
  g = numeric(n)
  u=runif(n)
  piic=cumsum(pii)
  for (i in (1:n)){
    if (u[i]<piic[1]){g[i]=1}else
      for (k in (1:(K-1))){if (piic[k]<u[i] & u[i]<piic[k+1]){g[i]=k+1}}
  }
  gg=NULL
  for (k in (1:K)){
    nk=sum(as.numeric(g==k)) 
    sigma22[[k]] = sigma2 
    piik = pii[k]
    muk = mu[,k]
    alphaik =  matrix(rnorm(nk*q,mean=0,sd=sqrt(sigma2)),ncol=q)
    yk = matrix(p*nk,ncol=p)
    Sigma = diag(rep(1,p))
    for (j in 1:(p-1)){
      Sigma[j,j+1] = s[1,k]
      Sigma[j+1,j] = s[1,k]
    }
    for (j in 1:(p-2)){
      Sigma[j,j+2] = s[2,k]
      Sigma[j+2,j] = s[2,k]
    }
    Qk = eigen(Sigma)$vectors[,1:q]  # Est-ce que ça suffit? 
    xxk=rep(1,nk)%*%t(muk)+t(Qk%*%t(alphaik))+matrix(rnorm(p*nk,mean=0,sd=sqrt(theta2)),ncol=p,nrow=nk)
    
    
    pca = PCA(xxk,nc=q,graph=FALSE)
    U = pca$svd$V
    K = diag(pca$eig[1:q])
    Qk = U%*%sqrt(K-sigma2*diag(rep(1,q)))
    
    Q[[k]]=Qk
    yk=rep(1,nk)%*%t(muk)+t(Qk%*%t(alphaik))+matrix(rnorm(p*nk,mean=0,sd=sqrt(theta2)),ncol=p,nrow=nk)
    y=rbind(y,yk)
    gg=c(gg,rep(k,nk))
    #matplot(muk+as.vector(Qk%*%alphaik)+matrix(rnorm(p*K,mean=0,sd=sqrt(theta2)),ncol=K,nrow=p),type="l",lty=3)
    #matlines(muk+as.vector(Qk%*%alphaik),lty=1)
    #matlines(muk,lty=1)
  }
  data = data.frame(y = y,g=gg)
  colnames(data) = c(paste("y",1:p,sep=""),"g")
  return(list(data = data,Q=Q,sigma2=sigma22, theta2=theta2))
}
# matplot(t(data1$data[,1:p]),col=data1$data$g,type="l")
# data1=data_gen_mixAcp(n=1000,K = 3,q = 4,p = 10,
#                s = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2), 
#                pii = c(0.2,0.35,0.45),
#                mu = matrix(c((0:9)^2/20,2*cos((0:9)/2)+1, rep(1,10)),nrow = 10,ncol=3),
#                SNR2 = 3,
#                sigma2 = 0.01)

data_gen_mixAcponegroup = function(n=1000,q = 4,p = 10,
                           s = c(0.7,-0.4),
                           mu = c(-3,-3,-3,-1,0,0,1,2,2,2),
                           sig2 = .01){
  
  Sigma = diag(rep(1,p))
  for (j in 1:(p-1)){
    Sigma[j,j+1] = s[1]
    Sigma[j+1,j] = s[1]
  }
  for (j in 1:(p-2)){
    Sigma[j,j+2] = s[2]
    Sigma[j+2,j] = s[2]
  }
  
  
  W = eigen(Sigma)$vectors[,1:q]  # Est-ce que ça suffit? 
  x = matrix(rnorm(n*q),n,q)
  X = x%*%t(W)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)
  
  pca = PCA(X,nc=q,graph=F)
  U = pca$svd$V
  K = diag(pca$eig[1:q])
  W = U%*%sqrt(K-sig2*diag(rep(1,q)))
  
  y = x%*%t(W)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)
  
    #matplot(muk+as.vector(Qk%*%alphaik)+matrix(rnorm(p*K,mean=0,sd=sqrt(theta2)),ncol=K,nrow=p),type="l",lty=3)
    #matlines(muk+as.vector(Qk%*%alphaik),lty=1)

  data = data.frame(y = y)
  colnames(data) = c(paste("y",1:p,sep=""))
  return(list(data = data,W=W,x=x))
}

# data1=data_gen_mixAcponegroup(n=1000,q = 4,p = 10,
#                               s = c(0.7,-0.4),
#                               mu = c(-3,-3,-3,-1,0,0,1,2,2,2),
#                               sig2 = .01)
# matplot(t(data1$data[,1:p]),type="l")
