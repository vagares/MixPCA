## packages
library(readxl)
library(fda)

setwd("~/GitHub/MixPCA")
source("2_estimators.r")

### chargement jeu de donnees Geometrie

geom = readxl::read_xlsx("C:\\Users\\vagares\\Documents\\MixPCA\\geometriePiece\\ResultsAddedForExcel.xlsx",col_names=TRUE,skip=1)
colnames(geom)[1] = "type"
geom$type = as.factor(geom$type)
#summary(geom)

### on prend 1 cahoutchouc (le premier)
geom = split(geom,geom$type)$P225M909
X = geom[,2:4]
Y = geom[,5:ncol(geom)]
Gap = seq(30,7,by=-0.5)
# matplot(Gap,t(Y),type="l")

### on remplace les donnees manquantes par un plateau (justification physique)
NAcurve = sort(unique(which(is.na(Y),arr.ind = TRUE)[,1]))
for (i in NAcurve){Y[i,which(is.na(Y[i,]))] = Y[i,min(which(is.na(Y[i,])))-1]}
# summary(Y) # OK, plus de donnees manquantes

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
est0 = estimates(cbind(Ycoef,X),K=3,maxits=100,
                 tol=1e-4, 
                 q = 3,
                 p=14,
                 nx=3,
                 verbose=TRUE)
