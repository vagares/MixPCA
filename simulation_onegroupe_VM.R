library(FactoMineR)
source("~/Library/CloudStorage/OneDrive-UniversitédeRennes1/RECHERCHE/COOPER/GEOMETRIE(1)/MixPCA-main/2_estimators_mixACP_ungroup_VM.R")


n = 10000
p = 10

# je construit une matrice bande diagonale (par simplicité)
Sigma = diag(rep(1,p))
for (j in 1:(p-1)){
  Sigma[j,j+1] = .7
  Sigma[j+1,j] = .7
}
for (j in 1:(p-2)){
  Sigma[j,j+2] = -.4
  Sigma[j+2,j] = -.4
}
image(Sigma[,rev(1:p)])

W = eigen(Sigma)$vectors[,1:q]  # Est-ce que ça suffit? 
x = matrix(rnorm(n*q),n,q)

sig2 = .01
mu = c(-3,-3,-3,-1,0,0,1,2,2,2)
X = x%*%t(W)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)

pca = PCA(X,nc=q)
U = pca$svd$V
K = diag(pca$eig[1:q])
Wtilde = U%*%sqrt(K-sig2*diag(rep(1,q)))
eigen(t(W)%*%W)$values # on trouve (1,1,1,1), c'est logique puisqu'on a construit une matrice orthonormée? 
eigen(t(Wtilde)%*%Wtilde)$values
pca$eig[1:q]
#eigen(t(pca$svd$V)%*%pca$svd$V)$values


Xtilde = x%*%t(Wtilde)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)
pca_tilde = PCA(Xtilde,nc=q)
pca_tilde$eig[1:q]
U = pca_tilde$svd$V
K = diag(pca_tilde$eig[1:q])
Wtilde = U%*%sqrt(K-sig2*diag(rep(1,q)))
What = U%*%sqrt(K-sig2*diag(rep(1,q)))

Xhat = x%*%t(What)+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(sig2)*matrix(rnorm(n*p),n,p)
par(mfrow=c(1,3))
image(Sigma[,rev(1:p)])
image(cov(X)[,rev(1:p)])
image(cov(Xhat)[,rev(1:p)])

### Si on utilise le code EM 
tol = 1e-4
maxit = 100
est0 = estimates(Xtilde,100,1e-4, q,p)
str(est0)
eigen(t(est0$Q[[1]])%*%est0$Q[[1]])$values
Xest0 = t(est0$alphai[[1]])%*%t(est0$Q[[1]])+matrix(rep(mu,n),n,p,byrow=TRUE)+sqrt(est0$sigma2i)*matrix(rnorm(n*p),n,p)
x11()
par(mfrow=c(2,2))
image(Sigma[,rev(1:p)])
title("Cov théorique")
image(cov(X)[,rev(1:p)])
title("Cov empirique")
image(cov(Xhat)[,rev(1:p)])
title("Cov estimée PCA")
image(cov(Xest0)[,rev(1:p)])
title("Cov estimée PPCA")

