
#Packages


#########################################################################

Estep = function(n,K = 3,q = 4,p = 10,nx=4,
                 pik,mu,beta, theta2=1,
                 sigma2=1,
                 x=x,y=y,
                 Q=Q,
                 C=C){
    alphai =list()
    alphai2 =list()
    alphai22 =list()
    alphai2Q =list()
    nn=sum(n)
    tau = matrix(0,nn,K)
   
    for (k in 1:K){
      betak=beta[[k]]
      Qk=Q[[k]]
      muk=mu[[k]]
      ez= muk%*%t(rep(1,nn)) + Qk%*%betak%*%t(x)
      vz = theta2 * (Qk%*%t(Qk))+sigma2*diag(p)
      for (i in (1:nn)){
      tau[i,k] = dmvnorm(y[i,], mean = ez[,i], sigma = vz)*pik[k]}
      alphai[[k]] = betak%*%t(x) + theta2*t(Qk) %*%solve(theta2*Qk%*%(t(Qk))+sigma2*diag(p)) %*% (t(y)-muk%*%t(rep(1,nn))-Qk%*%betak%*%t(x))
      alphai22[[k]] = matrix(0,ncol=nn,nrow=q*q)
      alphai2[[k]] = numeric(nn)
      alphai2Q[[k]] = numeric(nn)
      for (i in (1:nn)){
        alphai2[[k]][i] = sum(diag(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)))+sum(alphai[[k]][,i]^2)
        alphai2Q[[k]][i] = sum(diag(theta2*Qk%*%solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)%*%t(Qk)))+sum((Qk%*%alphai[[k]][,i])^2)
        alphai22[[k]][,i] = as.vector(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)+alphai[[k]][,i]%*%t(alphai[[k]][,i]))
      }
    }
    sommetau = apply(tau,1,sum)
    tau = tau/sommetau
    return(list(alphai=alphai,alphai2=alphai2,alphai22=alphai22,alphai2Q=alphai2Q,tau=tau))
}

Mstep = function(n=c(100,100,100),K = 3,q = 4,p = 10,nx=4,
                 x,y,z,C,muold,alphai,alphai2,alphai22, alphai2Q,tau){
  pik = apply(tau,2,mean)
  mu = list()
  Q = list()
  sigma2i = numeric(K)
  nn=sum(n)
  theta2i = matrix(0,nn,K)
  for (k in 1:K){
    tauik=tau[,k]
    alphaik = alphai[[k]]
    alphai2k = alphai2[[k]]
    alphai22k = alphai22[[k]]
    alphai2Qk = alphai2Q[[k]]
    muoldk=muold[[k]]
    #tau  = as.numeric(z==k)
    x=as.matrix(x)
    y=as.matrix(y)
    S1=0
    S2=0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (x[i,] %*%t(x[i,]))
    S2 =S2 + tauik[i] * x[i,] %*%t(alphaik[,i])
    }
    beta[[k]]= solve(S1)%*%S2
    #betak = matrix(as.numeric(beta[[k]]),q,q)
    for (i in (1:nn)){
    theta2i[i,k] =  tauik[i] * (alphai2k[i]  -2* t(x[i,])%*%t(beta[[k]])%*%alphaik[,i]+t(x[i,])%*%t(beta[[k]])%*%beta[[k]]%*%x[i,])
    }
    
    S1=0
    S2 = 0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (y[i,] - muoldk) %*%t(alphaik[,i])
    S2 = S2 + tauik[i] * matrix(alphai22k[,i],nrow=q,ncol=q)}
    Q[[k]]= S1 %*% solve(S2)
    
    mu[[k]] = apply(rep(1,p)%*%t(tauik) * (t(y)-Q[[k]]%*%alphaik),1,sum)/sum(tauik)
    #veci = numeric(nn)
    #for (i in 1:nn){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum(tauik * (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  theta2 = (1/(nn*q))*sum(theta2i)
  sigma2 = (1/(nn*p))*sum(sigma2i)
  sigma2i = sigma2i / nn
  return(list(pik=pik,beta = beta, mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,K=3,maxits=100,tol=0.01, q = 4,nx=4,p=10){
  y = data[,1:p]
  x = data[,(p+1):(p+nx)]
  #initialisation
  groupe=1:K
  C= sample(groupe,dim(y)[1],replace=TRUE)
  pik=rep(1/K,K)
  n = as.numeric(summary(as.factor(C)))

  mu = list()
  for (k in (1:K)){mu[[k]]=apply(y[which(C==k),],2,mean)}
  beta = list()
  for (k in (1:K)){beta[[k]]= matrix(rnorm(q*nx,mean=0,sd=2),ncol=nx)}
  theta2=1
  sigma2=1
  Q=list()
  #for (k in (1:K)){Q[[k]]=diag(p)[,1:q]}
  for (k in (1:K)){Q[[k]]=matrix(rep(1,p*q),ncol=q)}
  diff = 10
  iter=0
  while (diff > tol && iter < maxits){
    old.pik <- pik
    old.mu <- mu
    old.Q <- Q
    old.beta <- beta
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,nx,pik,mu,beta,theta2,sigma2,x,y,Q,C)
    alphai=Estepresults$alphai;alphai2=Estepresults$alphai2;alphai22=Estepresults$alphai22;alphai2Q=Estepresults$alphai2Q;tau=Estepresults$tau
    
    Mstepresults = Mstep(n,K,q,p,nx=4,x,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,tau)
    pik=Mstepresults$pik;beta=Mstepresults$beta;mu=Mstepresults$mu;Q=Mstepresults$Q;theta2=Mstepresults$theta2;sigma2=Mstepresults$sigma2;sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((pik[k]-old.pik[k])^2)+sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) + sum((beta[[k]] - old.beta[[k]])^2)}
    diff = sum(diff1) + sum((theta2-old.theta2)^2) + sum((sigma2-old.sigma2)^2)
    print(diff)
    iter = iter +1
  }
  return(list(pik=pik,mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i))
  }
  
  