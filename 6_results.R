
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
#setwd('C:\\Users\\vagares\\Dropbox\\Dropbox\\Analyse robuste modele mixte\\code\\code_VG\\Rik\\Results')
#load("res.Rdata")
############################################################
# use with package simsalar par
#val <- getArray(res)

#dimnames(val)

# names(dimnames(val))=c("var","cont","n.sim")
# 
# 
# dimnames(val)[[1]]=list("Scopt_beta0","Scopt_beta1","Scopt_sigma0","Scopt_sigma1")                  )
# df <- array2df(val)
# df <- df %>% separate(var, c("method","coefs"), sep = "_")




# df$coefsT <- ifelse(df$coefs == "beta0", 1, NA)
# df$coefsT <- ifelse(df$coefs == "beta1", 1, 1)
# df$coefsT <- ifelse(df$coefs == "sigma0", 1, 1)
# df$coefsT <- ifelse(df$coefs == "sigma1", 1, 1)


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
