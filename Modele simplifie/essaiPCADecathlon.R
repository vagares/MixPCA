library(FactoMineR)

data("decathlon")

dim(decathlon)
Tn = scale(decathlon[,1:10],center=TRUE,scale=FALSE)
res.pca = PCA(Tn, scale.unit=F, ncp=5, graph=F)
V = res.pca$svd$V # matrice W
t(V)%*%V
x1 = res.pca$ind$coord[1,] #coordonnees indiv 1
x1
t(V)%*%Tn[1,]
# on obtient bien le meme resultat 

U = res.pca$svd$U
Tn%*%U
KK = diag(res.pca$eig[1:5])
