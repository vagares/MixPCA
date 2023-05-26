
#Packages

library(mvtnorm)
#########################################################################

Estep = function(n,K = 3,q = 4,p = 10,nx=4,
                 piik,mu,beta, theta2=1,
                 sigma2=1,
                 x=x,y=y,
                 Q=Q,
                 C=C,cste=cste){
  if (cste){X = cbind(rep(1,nrow(x)),as.matrix(x))
  }else{ X = as.matrix(x)}
  
 # X = matrix(1,dim(x)[1],nx+1)
  #X[,2:(nx+1)] = x
    alphai =list()
    alphai2 =list()
    alphai22 =list()
    alphai2Q =list()
    nn = nrow(x)
    tau = matrix(0,nn,K)
   
    for (k in 1:K){
      betak=beta[[k]]
      Qk=Q[[k]]
      muk=mu[[k]]
      ez= muk%*%t(rep(1,nn)) + Qk%*%betak%*%t(X)
      vz = theta2 * (Qk%*%t(Qk))+sigma2*diag(p)
      for (i in (1:nn)){ tau[i,k] = dmvnorm(y[i,], mean = ez[,i], sigma = vz)}
      tau[,k] = tau[,k]*piik[k]
      alphai[[k]] = betak%*%t(X) + theta2*t(Qk) %*%solve(theta2*Qk%*%(t(Qk))+sigma2*diag(p)) %*% (t(y)-muk%*%t(rep(1,nn))-Qk%*%betak%*%t(X))
      alphai22[[k]] = matrix(0,ncol=nn,nrow=q*q)
      alphai2[[k]] = numeric(nn)
      alphai2Q[[k]] = numeric(nn)
      S1k = solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)
      T1k = sum(diag(theta2*S1k)) # trace
      T2k = sum(diag(theta2*Qk%*%S1k%*%t(Qk)))

      alphai2[[k]] = apply(alphai[[k]]^2,2,sum)+T1k
      alphai2Q[[k]] = apply((Qk%*%alphai[[k]])^2,2,sum)+T2k
      for (i in (1:nn)){
        alphai22[[k]][,i] = as.vector(theta2*S1k+alphai[[k]][,i]%*%t(alphai[[k]][,i]))
      }
   }
    sommetau = apply(tau,1,sum)
    tau = tau/sommetau
    return(list(alphai=alphai,alphai2=alphai2,alphai22=alphai22,alphai2Q=alphai2Q,tau=tau))
}

Mstep = function(n=c(100,100,100),K = 3,q = 4,p = 10,nx=4,
                 x,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,tau,cste=cste){

  if (cste){X = cbind(rep(1,nrow(x)),as.matrix(x))
  }else{ X = as.matrix(x)}
  y=as.matrix(y)
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

    S1=0
    S2=0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (X[i,] %*%t(X[i,]))
    S2 =S2 + tauik[i] * X[i,] %*%t(alphaik[,i])
    }
    if (svd(S1)$d[nx]>10^(-8)) {beta[[k]]= t(solve(S1)%*%S2)}else{beta[[k]]= t(solve(S1+diag(10^(-8),nx))%*%S2);print(paste("plus petite vp =",svd(S1)$d[nx]))}
    #betak = matrix(as.numeric(beta[[k]]),q,q)
    for (i in (1:nn)){
      theta2i[i,k] =  tauik[i]*(alphai2k[i]-2*t(X[i,])%*%t(beta[[k]])%*%alphaik[,i]+t(X[i,])%*%t(beta[[k]])%*%beta[[k]]%*%X[i,])
    }
    theta2 = (1/(nn*q))*sum(theta2i)
    
    S1 = 0
    S2 = 0
    for (i in (1:nn)){
      S1 = S1 + tauik[i] * (y[i,] - old.muk) %*%t(alphaik[,i])
      S2 = S2 + tauik[i] * matrix(alphai22k[,i],nrow=q,ncol=q)
    }
    Q[[k]] = S1 %*% solve(S2)
    
    mu[[k]] = apply(rep(1,p)%*%t(tauik) * (t(y)-Q[[k]]%*%alphaik),1,sum)/sum(tauik)
    #veci = numeric(nn)
    #for (i in 1:nn){ veci[i] = t(t(y[i,]-mu[[k]]))%*%Q[[k]]%*%(alphaik[,i])}
    sigma2i[k] = sum(tauik * (diag(t(t(y)-mu[[k]])%*%(t(y)-mu[[k]]))-2*diag(t(t(y)-mu[[k]])%*%Q[[k]]%*%alphaik)+alphai2Q[[k]]))
  }
  sigma2 = (1/(nn*p))*sum(sigma2i)
  sigma2i = sigma2i / nn
  return(list(piik=piik,beta = beta, mu=mu,Q=Q,theta2=theta2,sigma2=sigma2,sigma2i=sigma2i))
}

estimates = function(data,K=3, par_init = NULL, maxits=100,
                     tol=1e-4, 
                     q = 4,
                     p=10,
                     nx=4,
                     cste = TRUE,
                     verbose=TRUE){
  N = dim(data)[1]
  y = as.matrix(data[,1:p],N,p)
  x = data[,(p+1):(p+nx)]
  #x = matrix(1,dim(X)[1],nx+1)
  #x[,2:(nx+1)] = as.matrix(X)
  #initialisation
  #groupe=1:K
  #C= sample(groupe,dim(y)[1],replace=TRUE)
  if (is.null(par_init)){
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
       # ACP pour déterminer axes principaux (Q), scores indiv (alphaik) et erreur moyenne sigma
      pca = PCA(Z,scale.unit=F,ncp = q, graph = FALSE)
      Q[[k]] = pca$svd$V
      alphaik.init = pca$ind$coord
      sigma2_t[k] = mean(apply(Z-t(Q[[k]]%*%t(alphaik.init)),2,var))
      # regression des x sur les alphaik pour trouver betak (coefficients) et theta2 (residus)
      xk = as.matrix(x[which(C==k),])
      if (cste){  lm_xk = lm(alphaik.init ~ xk)
      }else{
        lm_xk = lm(alphaik.init ~ xk-1)
      }
      beta[[k]] = t(lm_xk$coefficients)
      theta2_t[k] = mean(sapply(1:ncol(alphaik.init),FUN=function(i){tmp=summary(lm_xk)[[i]];return(tmp$sigma)}))
  }
  theta2 = mean(theta2_t)
  sigma2 = mean(sigma2_t)
  } else {
    piik = par_init$piik
    mu = par_init$mu
    theta2 = par_init$theta2
    sigma2 = par_init$sigma2
    Q = par_init$Q
    beta = par_init$beta
  }  
  diff = 10
  iter=0
  loglik=-Inf
  llik = c()
  while (diff > tol && iter < maxits){
    old.piik <- piik
    old.mu <- mu
    old.Q <- Q
    old.beta <- beta
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,nx,piik,mu,beta,theta2,sigma2,x,y,Q,C,cste)
    alphai=Estepresults$alphai;alphai2=Estepresults$alphai2;alphai22=Estepresults$alphai22;alphai2Q=Estepresults$alphai2Q;tau=Estepresults$tau
    
    Mstepresults = Mstep(n,K,q,p,nx,x,y,z,C,old.mu,alphai,alphai2,alphai22,alphai2Q,tau,cste)
    piik=Mstepresults$piik;beta=Mstepresults$beta;mu=Mstepresults$mu;Q=Mstepresults$Q;theta2=Mstepresults$theta2;sigma2=Mstepresults$sigma2;sigma2i=Mstepresults$sigma2i
    diff1=numeric(K)
    for (k in (1:K)){diff1[k] = sum((piik[k]-old.piik[k])^2)+sum((mu[[k]] - old.mu[[k]])^2)+sum((Q[[k]] - old.Q[[k]])^2) + sum((beta[[k]] - old.beta[[k]])^2)}
    diff = sum(diff1) + sum((theta2-old.theta2)^2) + sum((sigma2-old.sigma2)^2)
	  loglik_old = loglik
    k = 1
    if ((sigma2i[1]<0) || (sigma2i[2]<0) || (sigma2i[3]<0)) {loglik = -Inf}else{loglik = ll_mixPCA(y,x,
						            mu=mu,
						            beta=beta,
                        Q=Q,
                        alphai=alphai,
                        sigma2=sigma2,
						            theta2 = theta2,
						            K,tau,piik=piik,cste=cste)}
    
    if (verbose) {print(paste("iteration",iter,", LL = ",round(loglik,2)))}
    
    llik = c(llik,loglik)
    cvce = EM_converged(loglik,loglik_old)$converged				
    iter = iter +1
  }
  G=numeric(N)
  G=apply(tau,1,which.max)
  return(list(piik=piik,mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2i=sigma2i,G=G,alphai=alphai,tau=tau,llik=llik))
  }

ll_mixPCA = function(y,x,mu=mu,beta=beta,Q=Q,alphai=alphai,sigma2=sigma2,theta2=theta2,K=K,tau=tau,piik=piik,cste = cste){
  if (cste){X = cbind(rep(1,nrow(x)),as.matrix(x))
  }else{ X = as.matrix(x)}
  
  N = nrow(y)
  p = ncol(y)
  tmp=numeric(K)
  LL=matrix(0,N,K)
  for (k in (1:K)){
  	tmp = matrix(rep(mu[[k]],N),N,p,byrow=TRUE)+ t(Q[[k]]%*%beta[[k]]%*%t(X))
  	vv = theta2 * Q[[k]]%*%t(Q[[k]])+ sigma2*diag(1,p)# + diag(sigma2i[k],p)
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
      if (loglik - previous_loglik < -1e-4) # allow for a little imprecision 
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
  
  
  