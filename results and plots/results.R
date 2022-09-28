source("simuls.r")
source("results.r")
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
setwd('C:\\Users\\vagares\\Documents\\MixPCA\\results and plots')
load("res.Rdata")
############################################################
# use with package simsalar par
val <- getArray(res)

dimnames(val)

names(dimnames(val))=c("var","cont","n.sim")


dimnames(val)[[1]]=list("Scopt_beta0","Scopt_beta1","Scopt_sigma0","Scopt_sigma1")                  )
df <- array2df(val)
df <- df %>% separate(var, c("method","coefs"), sep = "_")




df$coefsT <- ifelse(df$coefs == "beta0", 1, NA)
df$coefsT <- ifelse(df$coefs == "beta1", 1, 1)
df$coefsT <- ifelse(df$coefs == "sigma0", 1, 1)
df$coefsT <- ifelse(df$coefs == "sigma1", 1, 1)


df <- df %>% mutate(biais = value - coefsT,
                      biais2 = (value - coefsT)^2,
                      biais.relative = 100*(value- coefsT)/coefsT)


df$method = factor(df$method)
df$coefs = factor(df$coefs)
df$cont = factor(df$cont)


results <- df %>% group_by(method,cont,coefs) %>% dplyr::summarize(
  B = n(),
  coefs.mean = mean(value,na.rm = TRUE),
  coefs.sd = sd(value,na.rm =  TRUE),
  coefs.inf = coefs.mean - qnorm(0.975)*coefs.sd,
  coefs.sup = coefs.mean + qnorm(0.975)*coefs.sd,

  biais.mean = mean(biais,na.rm =  TRUE),
  biais.median = median(biais,na.rm =  TRUE)
  )




results$prettycoefs = factor(results$coefs, c("beta0","beta1","sigma0","sigma1"), c("hat(beta)[0]","hat(beta)[1]","hat(sigma)[0]","hat(sigma)[1]"))
results$prettycont = factor(results$cont, c(0,1,2,3),  c("(0/0)","(0.1/0)","(0/0.1)","(0.1/0.1)"))
results$prettymethod = factor(results$method, c("Scopt"),  c("Scopt"))

df$prettycoefs = factor(df$coefs, c("beta0","beta1","sigma0","sigma1"), c("hat(beta)[0]","hat(beta)[1]","hat(sigma)[0]^2","hat(sigma)[1]^2"))
df$prettycont = factor(df$cont, c(0,1,2,3),  c("(0/0)","(0.1/0)","(0/0.1)","(0.1/0.1)"))
df$prettymethod = factor(df$method,   c("Scopt"),  c("Scopt"))


setwd('C:\\Users\\vagares\\Documents\\MixPCA\\results and plots')
save(results,file="results.Rdata")
save(df,file="df.Rdata")
