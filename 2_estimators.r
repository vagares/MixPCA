
#Packages


#########################################################################

Estep = function(n,K = 3,q = 4,p = 10,nx=4,
                 pi,mu,beta, theta2=1,
                 sigma2=1,
                 x=x,
                 Q=Q,
                 C=C){
    alphai =list()
    alphai2 =list()
    alphai22 =list()
    alphai2Q =list()
    taui =list()
    nn=sum(n)
    for (k in 1:k){
      betak=beta[[k]]
      Qk=Q[[k]]
      muk=mu[[k]]
      ez= muk%*%t(rep(1,nn)) + Qk%*%betak%*%t(x)
      vz = theta2 * (Qk%*%t(Qk))+sigma2*diag(p)
      taui[[k]] = numeric(nn)
      for (i in (1:nn)){
      taui[[k]][i] = (dmvnorm(y[i,], mean = ez[,i], sigma = vz))*pi[k]}
      alphai[[k]] = betak%*%t(x) + theta2*t(Qk) %*%solve(theta2*Qk%*%(t(Qk))+sigma2*diag(p)) %*% (t(y)-muk%*%t(rep(1,nn))-Qk%*%betak%*%t(x))
      alphai22[[k]] = matrix(0,ncol=nn,nrow=q*q)
      alphai2[[k]] = numeric(nn)
      alphai2Q[[k]] = numeric(nn)
      for (i in (1:nn)){
        alphai2[[k]][i] = sum(diag(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)))+t(alphai[[k]][,i])%*%alphai[[k]][,i]
        alphai2Q[[k]][i] = sum(diag(theta2*Qk%*%solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)%*%t(Qk)))+t(alphai[[k]][,i])%*%t(Qk)%*%Qk%*%alphai[[k]][,i]
        alphai22[[k]][,i] = as.vector(theta2*solve(diag(q)+(theta2/sigma2)*t(Qk)%*%Qk)+alphai[[k]][,i]%*%t(alphai[[k]][,i]))
      }
    }
    for (i in (1:nn)){
      Staui = 0
      for (k in 1:k){
        Staui = Staui +taui[[k]][i] 
      }
      for (k in 1:k){taui[[k]][i]=taui[[k]][i]/Staui}
    }
    return(list(alphai=alphai,alphai2=alphai2,alphai22=alphai22,alphai2Q=alphai2Q,taui=taui))
}
Mstep = function(n=c(100,100,100),K = 3,q = 4,p = 10,nx=4,
                 x,y,z,C,muold,alphai,alphai2,alphai22, alphai2Q,taui){
  pi = unlist(lapply(taui,mean))
  mu = list()
  Q = list()
  theta2i = numeric(K)
  sigma2i = numeric(K)
  nn=sum(n)
  for (k in 1:k){
    tauik=taui[[k]]
    alphaik = alphai[[k]]
    alphai2k = alphai2[[k]]
    alphai22k = alphai22[[k]]
    alphai2Qk = alphai2Q[[k]]
    muoldk=muold[[k]]
    #tau  = as.numeric(z==k)
    x=as.matrix(x)
    y=as.matrix(y)
    theta2i[k] =  sum(tauik * (alphai2k  -2* x%*%t(betak)%*%alphaik+x%*%t(betak)%*%betak%*%t(x)))
    S1=0
    S2=0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (x[i,] %*%t(x[i,]))
    S2 =S2 + tauik[i] * x[i,] %*%t(alphaik[,i])}
    beta[[k]]= solve(S2)%*%S1

    S1=0
    S2 = 0
    for (i in (1:nn))
    {S1 = S1 + tauik[i] * (y[i,] - muoldk) %*%t(alphaik[,i])
    S2 = S2 + tauik[i] * matrix(alphai22k[,i],nrow=q,ncol=q)}
    Q[[k]]= t(solve(S2) %*% t(S1))
    mu[[k]] = apply(rep(1,p)%*%t(tauik) * (t(y)-Q[[k]]%*%alphaik),1,sum)/sum(tauik)
    sigma2i[k] = sum(tauik * (diag((y-t(mu[[k]]%*%t(rep(1,nn))))%*%t(y-t(mu[[k]]%*%t(rep(1,nn)))))-2*diag((y-mu[[k]])%*%Qk%*%alphaik)))
  }
  theta2 = (1/(nn*q))*sum(theta2i)
  sigma2 = (1/(nn*p))*sum(sigma2i)
  return(list(pi=pi,beta = beta, mu=mu,Q=Q,theta2=theta2,sigma2=sigma2))
}

estimates = function(data,K,maxits,tol, q = 4,p = 10,nx=4){
  y = data[,1:p]
  x = data[,(p+1):(p+nx)]
  #initialisation
  groupe=1:K
  C= sample(groupe,dim(y)[1],replace=TRUE)
  pi=rep(1/K,K)
  n = NULL
  for (k in (1:K)){n = c(n,sum(as.numeric(C==k)))}
  nx=dim(x)[2]
  #q=dim(y)[2]
  p=dim(y)[2]
  pi = c(0.3,0.3,0.3)
  mu1 = matrix(c(rep(1,10),rep(1,10), rep(1,10)),nrow = p,ncol=K)
  mu = list()
  for (k in (1:K)){mu[[k]]=mu1[,k]}
  beta1 = matrix(rnorm(4*4,mean=0,sd=2),ncol=nx)
  beta = list()
  for (k in (1:K)){beta[[k]]=beta1}
  theta2=1
  sigma2=1
  bfourier = create.fourier.basis(rangeval=c(0,1),nbasis=q)
  fmat = getbasismatrix(evalarg=seq(0,1,length=p),bfourier)
  Q1 = gramSchmidt(fmat[,-1])$Q
  Q=list()
  for (k in (1:K)){Q[[k]]=Q1}
  diff = 10
  iter=0
  while (diff > tol && iter < maxits){
    old.pi <- pi
    old.mu <- mu
    old.Q <- Q
    old.beta <- beta
    old.theta2 <- theta2
    old.sigma2 <- sigma2
    Estepresults = Estep(n,K,q,p,nx,pi,mu,beta,theta2,sigma2,x,Q,C)
    alphai=Estepresults$alphai;alphai2=Estepresults$alphai2;alphai22=Estepresults$alphai22;alphai2Q=Estepresults$alphai2Q;taui=Estepresults$taui
    
    Mstepresults = Mstep(n,K,q,p,nx=4,x,y,z,C,old.mu,alphai,alphai2,alphai22, alphai2Q,taui)
    pi=Mstepresults$pi;beta=Mstepresults$beta;mu=Mstepresults$mu;Q=Mstepresults$Q;theta2=Mstepresults$theta2;sigma2=Mstepresults$sigma2
    for (k in (1:K)){diff1 = sum(pi[k]-old.pi[k])+sum(mu[[k]] - old.mu[[k]])+sum(Q[[k]] - old.Q[[k]]) + sum(beta[[k]] - old.beta[[k]])}
    diff = sum(unlist(diff1)) + sum(theta2-old.theta2) + sum(sigma2-old.sigma2)
    iter = iter +1
  }
  return(list(pi=pi,mu=mu,beta=beta,Q=Q,theta2=theta2,sigma2=sigma2))
  
  