
#Packages


#########################################################################

Estep = function(n,K = 1,q = 4,p = 10,
                 mu, theta2=1,
                 sigma2=1,
                 y=y,
                 Q=Q,
                 C=C){
    alphai =list()
    alphai2 =list()
    alphai22 =list()
    alphai2Q =list()

    for (k in 1:K){
      Qk=Q[[k]]
      muk=mu[[k]]
      ez= muk%*%t(rep(1,n)) 
      vz = theta2 * (Qk%*%t(Qk))+sigma2*diag(p)
      #for (i in (1:n)){
      #tau[i,k] = dmvnorm(y[i,], mean = ez[,i], sigma = vz)*pik[k]}
      alphai[[k]] =  theta2*t(Qk) %*%solve(theta2*Qk%*%(t(Qk))+sigma2*diag(p)) %*% (t(y)-muk%*%t(rep(1,n)))
      alphai22[[k]] = matrix(0,ncol=n,nrow=q*q)
      alphai2[[k]] = numeric(n)
      alphai2Q[[k]] = numeric(n)
      for (i in (1:n)){
        alphai2[[k]][i] = sum(diag(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)))+sum(alphai[[k]][,i]^2)
        alphai2Q[[k]][i] = sum(diag(theta2*Qk%*%solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)%*%t(Qk)))+sum((Qk%*%alphai[[k]][,i])^2)
        alphai22[[k]][,i] = as.vector(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)+alphai[[k]][,i]%*%t(alphai[[k]][,i]))
      }
    }
    #sommetau = apply(tau,1,sum)
    #tau = tau/sommetau
    return(list(alphai=alphai,alphai2=alphai2,alphai22=alphai22,alphai2Q=alphai2Q))
}

Mstep = function(n=c(100,100,100),K = 1,q = 4,p = 10,
                 y,z,C,muold,alphai,alphai2,alphai22, alphai2Q){
  #pik = apply(tau,2,mean)
  mu = list()
  Q = list()
  sigma2i = numeric(K)
  theta2i = matrix(0,n,K)
  for (k in 1:K){
    #tauik=tau[,k]
    alphaik = alphai[[k]]
    alphai2k = alphai2[[k]]
    alphai22k = alphai22[[k]]
    alphai2Qk = alphai2Q[[k]]
    muoldk=muold[[k]]
    #tau  = as.numeric(z==k)
    y=as.matrix(y)
    S1=0
    S2=0
   
    for (i in (1:n)){
    theta2i[i,k] =   (alphai2k[i])
    }
    
    S1= 0
    S2 = 0
    for (i in (1:n))
    {S1 = S1 +  (y[i,] - muoldk) %*%t(alphaik[,i])
    S2 = S2 +  matrix(alphai22k[,i],nrow=q,ncol=q)}
    Q[[k]]= S1 %*% solve(S2)
    
    mu[[k]] = apply( (t(y)-Q[[k]]%*%alphaik),1,sum)
    #veci = numeric(n)
    #for (i in 1:n){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum( (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  theta2 = (1/(n*q))*sum(theta2i)
  sigma2 = (1/(n*p))*sum(sigma2i)
  sigma2i = sigma2i / n
  return(list(mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,maxits=100,tol=0.01, q = 4,p=10){
  K=1
  y = data[,1:p]
  #initialisation
  n = dim(data)[1]

  mu = list()
  mu[[1]]=apply(y,2,mean)
  theta2=1
  sigma2=1
  Q=list()
  #for (k in (1:K)){Q[[k]]=diag(p)[,1:q]}
  Q[[1]]=matrix(rep(1,p*q),ncol=q)
  diff = 10
  iter=0
  while (diff > tol && iter < maxits){
    old.mu <- mu
    old.Q <- Q
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,mu,theta2,sigma2,y,Q,C)
    alphai=Estepresults$alphai;alphai2=Estepresults$alphai2;alphai22=Estepresults$alphai22;alphai2Q=Estepresults$alphai2Q;#tau=Estepresults$tau
    
    Mstepresults = Mstep(n,K,q,p,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q)
    mu=Mstepresults$mu;Q=Mstepresults$Q;theta2=Mstepresults$theta2;sigma2=Mstepresults$sigma2;sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) }
    diff = sum(diff1) + sum((theta2-old.theta2)^2)
    print(diff)
    iter = iter +1
  }
  #G=numeric(n)
  #G=apply(tau,1,which.max)
  return(list(mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i))
  }
  
  