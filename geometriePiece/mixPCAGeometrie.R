## packages
library(readxl)
library(fda)
library(tidyfun)
library(viridis)
library(tidyverse)

setwd("~/GitHub/MixPCA")
source("2_estimators.r")

### chargement jeu de donnees Geometrie

geom = readxl::read_xlsx("geometriePiece/ResultsAddedForExcel.xlsx",col_names=TRUE,skip=1)
colnames(geom)[1] = "type"
geom$type = as.factor(geom$type)
#summary(geom)

### on prend 1 cahoutchouc (le premier)
geom = split(geom,geom$type)$P225M909
Xgeom = geom[,2:4]
Y = geom[,5:ncol(geom)]
Gap = seq(30,7,by=-0.5)
# matplot(Gap,t(Y),type="l")

### on remplace les donnees manquantes par un plateau (justification physique)
NAcurve = sort(unique(which(is.na(Y),arr.ind = TRUE)[,1]))
for (i in NAcurve){Y[i,which(is.na(Y[i,]))] = Y[i,min(which(is.na(Y[i,])))-1]}
# summary(Y) # OK, plus de donnees manquantes

################ graphes tidyfun ##############################################
geom = tibble(
  height = X$Height, 
  shape = X$Shape, 
  thickness = X$Thickness)
geom$effort = tfd(as.matrix(Y), arg = Gap)

geom_curve =
  geom %>% 
  ggplot(aes(tf = effort)) + geom_spaghetti(alpha=0.15,color = "steelblue") 
geom_curve+xlab("Gap (mm)")+ylab("Efforts (N)")
###############################################################################


### projection dans une base de splines pour obtenir des coefficients splines (smoothing penalise)
splbasis = create.bspline.basis(rangeval=range(Gap),norder=4,breaks=seq(min(Gap),max(Gap),length=12))
gcv = 1:21
for (i in 1:21){
  lambda = exp(i-10)
  fdparTemp = fdPar(splbasis,Lfdobj = 2,lambda=lambda)
  smoothdata = smooth.basis(argvals=Gap,t(Y),fdParobj = fdparTemp)
  gcv[i] = mean(smoothdata$gcv)
}

lambda = exp(which.min(gcv)-10)
fdparTemp = fdPar(splbasis,Lfdobj = 2,lambda=lambda)
smoothdata = smooth.basis(argvals=Gap,t(Y),fdParobj = fdparTemp)
# plotfit.fd(t(Y),Gap,smoothdata$fd,pch=20,cex=0.5) #on est tres proche de l'interpolation

Ycoef = t(smoothdata$fd$coefs)
dim(Ycoef) # 14 coefficients pour n individus

### on estime avec mixPCA
data = cbind(Ycoef,Xgeom)
K=3;maxits=100;tol=1e-4; q = 3;p=14;nx=3

est0 = estimates(cbind(Ycoef,Xgeom),K=3,par_init=NULL,maxits=100,
                 tol=1e-4, 
                 q = 3,
                 p=14,
                 nx=3,
                 verbose=TRUE)

save(est0,file="geometriePiece/resultsGeom-31032023.Rdata")
#load("geometriePiece/resultsGeom-08032023.Rdata")

est0$piik  # taille des groupes
est0$theta2 # on obtient une valeur bizarre (très grande)
est0$sigma2i

apply(est0$tau,1,max) # les probas d'appartenance sont plutot tranchées, 17 indiv sur 498 avec proba max <0.95

est0$G #groupe d'appartenance, MAP

#### CLUSTERING
# graphique des courbes de départ par groupe
matplot(Gap,t(Y),type="l",col=est0$G)
par(mfrow=c(1,3))
for (i in 1:K){
matplot(Gap,t(Y[which(est0$G==i),]),type="l",col=i,ylab="Efforts",ylim=c(0,42))
}

#### MOYENNES MU
est0$mu #moyenne par groupe (en coefficients splines, sous forme de liste)

#reconstruction des courbes moyennes fonctionnelles 
splbasis = create.bspline.basis(rangeval=range(Gap),norder=4,breaks=seq(min(Gap),max(Gap),length=12))
mufct = sapply(est0$mu,FUN=function(x) {ft =fd(x,splbasis) ; return(eval.fd(Gap,ft))})
# representation
layout(matrix(c(1,1,1,2,3,4),ncol=2,nrow=3))
matplot(Gap,mufct,type="l",ylab="moyennes fonctionnelles",lty=1,col=c("coral1","palegreen3","steelblue"))
abline(h=0,lwd=0.5,lty=3)
couleur = c("coral1","palegreen3","steelblue")
for (i in 1:K){
  par(mar=c(2,2,1,1))
  matplot(Gap,t(Y[which(est0$G==i),]),type="l",col=couleur[i],ylab="Efforts",ylim=c(0,42),lwd=0.5)
  lines(Gap,mufct[,i],lwd=2)
}

###### ACP FONCTIONNELLE - COMPOSANTES PRINCIPALES
est0$Q

## reconstruction fonctionnelle
Qfct = lapply(est0$Q,FUN=function(x) {ft =fd(x,splbasis) ; return(eval.fd(Gap,ft))})
par(mfrow=c(1,3))
for (i in 1:K){
  matplot(Gap,Qfct[[i]],type="l",col=2:5)
}
### dans Q[[i]], est ce qu'on a les 3 composantes du groupe i ou la composante i pour les 3 groupes (on pourrait croire qu'on a la deuxieme) 

## representation comme variabilite autour de la moyenne (par contre, on a pas la valeur propre pour ajuster le coeff)

#groupe 1
par(mfrow=c(1,3))
k=1
for(i in 1:3){
  matplot(Gap,t(Y[which(est0$G==k),]),type="l",col="gray",ylab="Efforts",ylim=c(0,42),lwd=0.5)
  lines(Gap,mufct[,k],col=1)
  lines(Gap,mufct[,k]+4*Qfct[[k]][,i],type="p",pch="+",col=i)
  lines(Gap,mufct[,k]-4*Qfct[[k]][,i],type="p",pch="-",col=i,lwd=2)
}

#groupe 2
par(mfrow=c(1,3))
k=2
for(i in 1:3){
  matplot(Gap,t(Y[which(est0$G==k),]),type="l",col="gray",ylab="Efforts",ylim=c(0,42),lwd=0.5)
  lines(Gap,mufct[,k],col=1)
  lines(Gap,mufct[,k]+4*Qfct[[k]][,i],type="p",pch="+",col=i)
  lines(Gap,mufct[,k]-4*Qfct[[k]][,i],type="p",pch="-",col=i,lwd=2)
}

#groupe 3
par(mfrow=c(1,3))
k=3
for(i in 1:3){
  matplot(Gap,t(Y[which(est0$G==k),]),type="l",col="gray",ylab="Efforts",ylim=c(0,42),lwd=0.5)
  lines(Gap,mufct[,k],col=1)
  lines(Gap,mufct[,k]+4*Qfct[[k]][,i],type="p",pch="+",col=i)
  lines(Gap,mufct[,k]-4*Qfct[[k]][,i],type="p",pch="-",col=i,lwd=2)
}

## les composantes semblent plus interpretables maintenant, on constate
## que dans certains groupes, elles ne semblent pas être dans le bon ordre...

##### coefficients de regression beta

est0$beta



