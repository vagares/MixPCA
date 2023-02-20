
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
   
    theta2i[,k] =   alphai2k
    
    S1= 0
    S2 = 0
    for (i in (1:n))
    {S1 = S1 +  (y[i,] - muoldk) %*%t(alphaik[,i])
    S2 = S2 +  matrix(alphai22k[,i],nrow=q,ncol=q)
	}
    Q[[k]]= S1 %*% solve(S2)
    
    mu[[k]] = apply( (t(y)-Q[[k]]%*%alphaik),1,sum)/n
    #veci = numeric(n)
    #for (i in 1:n){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum( (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  theta2 = (1/(n*q))*sum(theta2i)
  sigma2 = (1/(n*p))*sum(sigma2i)
  sigma2i = sigma2i / n
  return(list(mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,
                     maxits=100,
                     tol=1e-4, 
                     q = 4,
                     p=10,
                     verbose=TRUE){
  K=1
  
  #initialisation
  n = dim(data)[1]
  y = as.matrix(data[,1:p],n,p)
  mu = list()
  mu[[1]]=apply(y,2,mean)
  theta2=1
  sigma2=1
  Q=list()
  #for (k in (1:K)){Q[[k]]=diag(p)[,1:q]}
  Q[[1]]=matrix(rep(1,p*q),ncol=q)
  yc = y-matrix(mu[[1]],n,p,byrow=TRUE)
  Q[[1]] = svd(yc,nu=q,nv=q)$v
  diff = 10
  iter=0
  loglik=-Inf
  cvce=FALSE  
  while (!cvce && iter < maxits){
    old.mu <- mu
    old.Q <- Q
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,mu,theta2,sigma2,y,Q,C)
    alphai=Estepresults$alphai;
	alphai2=Estepresults$alphai2;
	alphai22=Estepresults$alphai22;
	alphai2Q=Estepresults$alphai2Q;
    Mstepresults = Mstep(n,K,q,p,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q)
    mu=Mstepresults$mu;
	Q=Mstepresults$Q;
	theta2=Mstepresults$theta2;
	sigma2=Mstepresults$sigma2;
	sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) }
    diff = sum(diff1) + sum((theta2-old.theta2)^2)
    loglik_old = loglik
    k = 1
    if (sigma2i<0){loglik = -Inf}else{
    loglik = ll_mixPCA_onegroup(y,
                                mu=mu,
                                Q=Q,
                                alphai=alphai,
                                sigma2i=sigma2i)}
    if (verbose) {print(paste("iteration",iter,", LL = ",round(loglik,2)))}
    
    cvce = EM_converged(loglik,loglik_old)$converged
    iter = iter +1
  }
  #G=numeric(n)
  #G=apply(tau,1,which.max)
  return(list(mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i,alphai=alphai))
  }

ll_mixPCA_onegroup = function(X,mu=mu,Q=Q,alphai=alphai,sigma2i=sigma2i){
  N = nrow(X)
  p = ncol(X)
  #tmp = (X-matrix(rep(mu,N),N,p,byrow=TRUE)-t(alphai[[1]])%*%t(Q))
  #LL = -N*p*(log(2*pi)+log(sigma2i))-sum(diag(tmp%*%t(tmp)))
  tmp = matrix(rep(mu[[1]],N),N,p,byrow=TRUE)
  vv =  Q[[1]]%*%t(Q[[1]])+ diag(sigma2i,p)
  LL =  dmvnorm(X-tmp,rep(0,p),vv,log=FALSE)
  L=-sum(log((LL)))
  return(L)
}
  
  
EM_converged <-
  function(loglik, previous_loglik, threshold = 1e-4) {
    
    converged = 0;
    decrease = 0;
    if(!(previous_loglik==-Inf)){
      if (loglik - previous_loglik < -1e-2) # allow for a little imprecision 
      {
        print(paste("******likelihood decreased from ",previous_loglik," to ", loglik,sep=""),quote = FALSE)
        decrease = 1;
      }
      
      delta_loglik = abs(loglik - previous_loglik);
      avg_loglik = (abs(loglik) + abs(previous_loglik) + threshold)/2;
      bb = ((delta_loglik/avg_loglik) < threshold)
      if (bb) {converged = 1}
    }
    
    res <- NULL
    res$converged <- converged
    res$decrease <- decrease
    return(res)
    
  }

  























  
