library(fields)




B=50
Xest0 = array(0,c(n,p,B))

for (b in 1:B){
  load(paste("simulations0323/SimMars-i=",b,".Rdata",sep=""))
  est0=sim$est
  X = sim$datasim
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

  