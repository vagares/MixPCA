
#Packages

library(mvtnorm)
#########################################################################

Estep = function(n,K = 3,q = 4,p = 10,nx=4,
                 piik,mu,beta, theta2=1,
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
      tau[i,k] = dmvnorm(y[i,], mean = ez[,i], sigma = vz)*piik[k]}
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
                 x,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,tau){
  piik = apply(tau,2,mean)
  mu = list()
  Q = list()
  beta=list()
  sigma2i = numeric(K)
  nn=sum(n)
  theta2i = matrix(0,nn,K)
  for (k in 1:K){
    tauik=tau[,k]
    alphaik = alphai[[k]]
    alphai2k = alphai2[[k]]
    alphai22k = alphai22[[k]]
    alphai2Qk = alphai2Q[[k]]
    old.muk=old.mu[[k]]
    #tau  = as.numeric(z==k)
    x=as.matrix(x)
    y=as.matrix(y)
    S1=0
    S2=0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (x[i,] %*%t(x[i,]))
    S2 =S2 + tauik[i] * x[i,] %*%t(alphaik[,i])
    }
    if (svd(S1)$d[nx]>10^(-8)) {beta[[k]]= t(solve(S1)%*%S2)}else{beta[[k]]= t(solve(S1+diag(10^(-8),nx))%*%S2)}
    #betak = matrix(as.numeric(beta[[k]]),q,q)
    for (i in (1:nn)){
    theta2i[i,k] =  tauik[i] * (alphai2k[i]  -2* t(x[i,])%*%t(beta[[k]])%*%alphaik[,i]+t(x[i,])%*%t(beta[[k]])%*%beta[[k]]%*%x[i,])
    }
    
    S1= 0
    S2 = 0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (y[i,] - old.muk) %*%t(alphaik[,i])
    S2 = S2 + tauik[i] * matrix(alphai22k[,i],nrow=q,ncol=q)}
    Q[[k]]= S1 %*% solve(S2)
    
    mu[[k]] = apply(rep(1,p)%*%t(tauik) * (t(y)-Q[[k]]%*%alphaik),1,sum)/sum(tauik)
    #veci = numeric(nn)
    #for (i in 1:nn){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum(tauik * (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  theta2 = (1/(nn*q))*sum(theta2i)
  sigma2 = (1/(nn*p))*sum(sigma2i)
  sigma2i = sigma2i / (nn*p)
  return(list(piik=piik,beta = beta, mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,K=3,maxits=100,
                     tol=1e-4, 
                     q = 4,
                     p=10,
                     nx=4,
                     verbose=TRUE){
  N = dim(data)[1]
  y = as.matrix(data[,1:p],N,p)
  x = data[,(p+1):(p+nx)]
  #initialisation
  groupe=1:K
  #C= sample(groupe,dim(y)[1],replace=TRUE)
  km=kmeans(y,3,nstart = 10)
  C= km$cluster
  piik=km$size/N
  n = as.numeric(summary(as.factor(C)))

  mu = list()
  for (k in (1:K)){mu[[k]]=apply(y[which(C==k),],2,mean)}
  
  Q=list()
  sigma2_t=numeric(K)
  theta2_t=numeric(K)
  beta = list()
  #for (k in (1:K)){Q[[k]]=diag(p)[,1:q]}
  for (k in (1:K)){
      nk = length(C[which(C==k)])
      Z = t(t(y[which(C==k),]) - mu[[k]])
      pca = PCA(Z,ncp = q, graph = FALSE)
      V = pca$svd$V
      U = pca$svd$U
      KK = diag(pca$eig[1:q])
      Q[[k]] = (1/sqrt(nk))*(V%*%sqrt(KK))
      sigma2_t[k] = (sd(Z -t(V%*%sqrt(KK)%*%t(U))))^2
      xk = as.matrix(x[which(C==k),])
      xx = U
      lm_xk = lm(xx ~ xk)
      beta[[k]] = t(lm_xk$coefficients[-1,])
      theta2_t[k] = (sd(lm_xk$residuals))^2
  }
  theta2 = mean(theta2_t)
  sigma2 = mean(sigma2_t)
  #beta = list()
  #for (k in (1:K)){beta[[k]]= matrix(rnorm(q*nx,mean=0,sd=2),ncol=nx)}
  #theta2=1
  #sigma2=1
  diff = 10
  iter=0
  loglik=-Inf
  while (diff > tol && iter < maxits){
    old.piik <- piik
    old.mu <- mu
    old.Q <- Q
    old.beta <- beta
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,nx,piik,mu,beta,theta2,sigma2,x,y,Q,C)
    alphai=Estepresults$alphai;alphai2=Estepresults$alphai2;alphai22=Estepresults$alphai22;alphai2Q=Estepresults$alphai2Q;tau=Estepresults$tau
    
    Mstepresults = Mstep(n,K,q,p,nx,x,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,tau)
    piik=Mstepresults$piik;beta=Mstepresults$beta;mu=Mstepresults$mu;Q=Mstepresults$Q;theta2=Mstepresults$theta2;sigma2=Mstepresults$sigma2;sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((piik[k]-old.piik[k])^2)+sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) + sum((beta[[k]] - old.beta[[k]])^2)}
    diff = sum(diff1) + sum((theta2-old.theta2)^2) + sum((sigma2-old.sigma2)^2)
    print(diff)
	  loglik_old = loglik
    k = 1
    if ((sigma2i[1]<0) || (sigma2i[2]<0) || (sigma2i[3]<0)) {loglik = -Inf}else{loglik = ll_mixPCA(y,x,
						            mu=mu,
						            beta=beta,
                        Q=Q,
                        alphai=alphai,
                        sigma2i=sigma2i,
						            theta2 = theta2,
						            K,tau,piik=piik)}
    
    if (verbose) {print(paste("iteration",iter,", LL = ",round(loglik,2)))}
    
    cvce = EM_converged(loglik,loglik_old)$converged				
    iter = iter +1
  }
  G=numeric(N)
  G=apply(tau,1,which.max)
  return(list(piik=piik,mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai,tau=tau))
  }

ll_mixPCA = function(y,x,mu=mu,beta=beta,Q=Q,alphai=alphai,sigma2i=sigma2i,theta2=theta2,K=K,tau=tau,piik=piik){
  N = nrow(y)
  p = ncol(y)
  tmp=numeric(K)
  LL=matrix(0,N,K)
  for (k in (1:K)){
  	tmp = matrix(rep(mu[[k]],N),N,p,byrow=TRUE)+ t(Q[[k]]%*%beta[[k]]%*%t(x))
  	vv = theta2 * Q[[k]]%*%t(Q[[k]])+ diag(sigma2i[k],p)
	  LL[,k] = piik[k]* dmvnorm(y-tmp,rep(0,p),vv,log=FALSE)
	}
	L=-sum(log(apply(LL,2,sum)))
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
      if (is.na(bb)==FALSE){if (bb) {converged = 1}}
    }
    
    res <- NULL
    res$converged <- converged
    res$decrease <- decrease
    return(res)
    
  }																		 
  
  
  