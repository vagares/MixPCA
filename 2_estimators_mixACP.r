
#Packages


#########################################################################

Estep = function(n,K = 3,q = 4,p = 10,
                 piik,mu, theta2=1,
                 sigma2=1,
                 y=y,
                 Q=Q,
                 C=C){
    alphai =list()
    alphai2 =list()
    alphai22 =list()
    alphai2Q =list()
    nn=sum(n)
    
    tau = matrix(0,nn,K)
    
    for (k in 1:K){
      Qk=Q[[k]]
      muk=mu[[k]]
      ez= muk%*%t(rep(1,nn)) 
      vz = theta2 * (Qk%*%t(Qk))+sigma2*diag(p)
      for (i in (1:nn)){
      tau[i,k] = dmvnorm(y[i,], mean = ez[,i], sigma = vz)*piik[k]}
      alphai[[k]] =  theta2*t(Qk) %*%solve(theta2*Qk%*%(t(Qk))+sigma2*diag(p)) %*% (t(y)-muk%*%t(rep(1,nn)))
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

Mstep = function(n=c(100,100,100),K = 3,q = 4,p = 10,
                 y,z,C,muold,alphai,alphai2,alphai22, alphai2Q,tau){
  piik = apply(tau,2,mean)
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
    y=as.matrix(y)
    S1=0
    S2=0
   
    for (i in (1:nn)){
    theta2i[i,k] =  tauik[i] * (alphai2k[i])
    }
    
    S1= 0
    S2 = 0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (y[i,] - muoldk) %*%t(alphaik[,i])
    S2 = S2 + tauik[i] * matrix(alphai22k[,i],nrow=q,ncol=q)}
    Q[[k]]= S1 %*% solve(S2)
    
    mu[[k]] = apply(rep(1,p)%*%t(tauik) * (t(y)-Q[[k]]%*%alphaik),1,sum)/sum(tauik)
    #veci = numeric(n)
    #for (i in 1:n){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum(tauik * (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  theta2 = (1/(nn*q))*sum(theta2i)
  sigma2 = (1/(nn*p))*sum(sigma2i)
  sigma2i = sigma2i / nn
  return(list(piik=piik, mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,
                     K=3,
                     maxits=100,
                     tol=1e-4, 
                     q = 4,
                     p=10,
                     verbose=TRUE){
  N = dim(data)[1]
  y = as.matrix(data[,1:p],N,p)
  #initialisation
  groupe=1:K
  km=kmeans(y,3,nstart = 10)
  C= km$cluster
  piik=km$size/N
  n = as.numeric(summary(as.factor(C)))

  mu = list()
  for (k in (1:K)){
    nk=n[k]
	  mu[[k]]=apply(y[which(C==k),],2,mean)}
  yc = matrix(N*p,N,p)
  Q=list()
  for (k in (1:K)){
    nk=n[k]
  	yc[which(C==k),]=y[which(C==k),]-matrix(rep(mu[[k]],nk),nk,p,byrow=TRUE)
	  Q[[k]] = svd(yc[which(C==k),],nu=q,nv=q)$v
	  }
  theta2=1
  sigma2=1
  #Q=list()
  #for (k in (1:K)){Q[[k]]=diag(p)[,1:q]}
  #for (k in (1:K)){Q[[k]]=matrix(rep(1,p*q),ncol=q)}
  diff = 10
  iter=0
  loglik=-Inf
  cvce=FALSE
  while (!cvce && iter < maxits){
    old.piik <- piik
    old.mu <- mu
    old.Q <- Q
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,piik,mu,theta2,sigma2,y,Q,C)
    alphai=Estepresults$alphai;
	alphai2=Estepresults$alphai2;
	alphai22=Estepresults$alphai22;
	alphai2Q=Estepresults$alphai2Q;
	tau=Estepresults$tau
    
    Mstepresults = Mstep(n,K,q,p,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,tau)
    piik=Mstepresults$piik;
	mu=Mstepresults$mu;
	Q=Mstepresults$Q;
	theta2=Mstepresults$theta2;
	sigma2=Mstepresults$sigma2;
	sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((piik[k]-old.piik[k])^2)+sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) }
    diff = sum(diff1) + sum((theta2-old.theta2)^2)
    print(diff)
	
	loglik_old = loglik
    k = 1
    if ((sigma2i[1]<0) || (sigma2i[2]<0) || (sigma2i[3]<0) || (is.na(sigma2i[1])==TRUE)){loglik = -Inf}else{loglik = ll_mixPCA(y,
						mu=mu,
                        Q=Q,
                        alphai=alphai,
                        sigma2i=sigma2i,K,tau,piik=piik)}
    if (verbose) {print(paste("iteration",iter,", LL = ",round(loglik,2)))}
    
    cvce = EM_converged(loglik,loglik_old)$converged
	
    iter = iter +1
  }
  G=numeric(N)
  G=apply(tau,1,which.max)
  return(list(piik=piik,mu=mu,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai))
  }


ll_mixPCA = function(X,mu=mu,Q=Q,alphai=alphai,sigma2i=sigma2i,K=K,tau=tau,piik=piik){
  N = nrow(y)
  p = ncol(y)
  tmp=numeric(K)
  LL=matrix(0,N,K)
  for (k in (1:K)){
    tmp = matrix(rep(mu[[k]],N),N,p,byrow=TRUE)
    vv = theta2 * Q[[k]]%*%t(Q[[k]])+ diag(sigma2i[k],p)
    LL[,k] = piik[k]* dmvnorm(y-tmp,rep(0,p),vv,log=FALSE)
  }
  L=sum(log(apply(LL,2,sum)))
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
  
  