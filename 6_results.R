
###########################################################
library(simsalapar)
library(tibble)
library(nlme)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(remotes)
###########################################################

B   = 100   # nombre de repetitions
q   = 3   # nombre d'axes principaux 
nx  = 4   # nombre de covariables
p   = 10
K   = 3
pii = c(0.2,0.35,0.45)

mu        = matrix(c((0:(p-1))^2/(2*p),2*cos((0:(p-1))/2)+1, rep(1,p)),nrow = p,ncol=K)
s         = matrix(c(0.7,-0.4,0.7,0.4,0.8,0.2),ncol=3,nrow=2)
beta      = list()
beta[[1]] = matrix(c(-1,-1,-1,1),q,nx)
beta[[2]] = matrix(c(1,1,0,0),q,nx)
beta[[3]] = matrix(c(-2,0,2,2),q,nx)

res.simu = read.table("simulations0323/analyse/ResultsSummary.txt",header=TRUE)

## decroissance vraisemblance
sum(res.simu$llik.dec)
table(as.factor(res.simu$llik.dec),as.factor(res.simu$n))
table(as.factor(res.simu$llik.dec),as.factor(res.simu$SNR1))
table(as.factor(res.simu$llik.dec),as.factor(res.simu$SNR2))

## estimation piik
par(mfrow=c(1,3))
for (k in 1:3){
  summary(res.simu[,5+k])
  boxplot(res.simu[,5+k])
  abline(h=pii[k],lwd=1.5,col=2)
}


## estimation theta2
MSE.theta2 = (res.simu$theta2.hat-res.simu$theta2.sim)^2
boxplot(MSE.theta2)
abline(h=0,lwd=1.5,col=2)



##########################################
df <- df %>% mutate(biais = value,
                      biais2 = (value)^2,
)


df$G = factor(df$G)
df$coefs = factor(df$coefs)
df$coefs2 = factor(df$coefs2)


results <- df %>% group_by(G,coefs2,coefs) %>% dplyr::summarize(
  B = n(),
  coefs.mean = mean(value,na.rm = TRUE),
  coefs.sd = sd(value,na.rm =  TRUE),
  coefs.inf = coefs.mean - qnorm(0.975)*coefs.sd,
  coefs.sup = coefs.mean + qnorm(0.975)*coefs.sd,

  biais.mean = mean(biais,na.rm =  TRUE),
  biais.median = median(biais,na.rm =  TRUE)
  )



save(results,file="results.Rdata")
save(df,file="df.Rdata")
