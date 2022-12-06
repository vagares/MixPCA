library(ggplot2)
library(ggpubr)


#df$prettycoefs = factor(df$coefs,  c("beta0", "beta1", "sigma0", "sigma1"),  c("hat(beta)[0]", "hat(beta)[1]", "hat(sigma)[0]^2", "hat(sigma)[1]^2"))
#df$prettycont = df$cont
#df$prettymethod = df$method


#cols2  =  c("Scopt" = "#009E73")
#lt2  =  c("Scopt" = "solid")
#shaps2  =  c("Scopt" = 2)



X11()
mu  =  ggplot(data = as.data.frame(df[df$coefs2 == "mu", ]),  aes(x = coefs2,  y = value)) +
  facet_grid( ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  #geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(mu)) +
  ylab(quote(hat(mu))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")

mu

X11()
Q  =  ggplot(data = as.data.frame(df[df$coefs2 == "Q", ]),  aes(x = coefs2,  y = value)) +
  facet_grid( ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab(quote(Q)) +
  ylab(quote(hat(Q))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
Q
X11()
sigma  =  ggplot(data = as.data.frame(df[df$coefs2 == "sigma", ]),   aes(x = coefs2,  y = value)) +
  facet_grid( ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_hline(yintercept = 0) +
  xlab(quote(sigma)) +
  ylab(quote(hat(sigma))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
sigma


