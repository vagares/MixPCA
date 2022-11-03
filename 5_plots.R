library(ggplot2)
library(ggpubr)


#df$prettycoefs = factor(df$coefs,  c("beta0", "beta1", "sigma0", "sigma1"),  c("hat(beta)[0]", "hat(beta)[1]", "hat(sigma)[0]^2", "hat(sigma)[1]^2"))
#df$prettycont = df$cont
#df$prettymethod = df$method


#cols2  =  c("Scopt" = "#009E73")
#lt2  =  c("Scopt" = "solid")
#shaps2  =  c("Scopt" = 2)



x11()
pi  =  ggplot(data = as.data.frame(df[df$coefs2 == "pi", ]),  aes(x = coefs2,  y = value,  group = G)) +
  facet_grid(. ~ G,
             labeller = label_parsed,
             scales = "free_y") +
  geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(pi)) +
  ylab(quote(hat(pi))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
pi
X11()
mu  =  ggplot(data = as.data.frame(df[df$coefs2 == "mu", ]),  aes(x = coefs2,  y = value,  group = G)) +
  facet_grid(G ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(mu)) +
  ylab(quote(hat(mu))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")

mu
X11()
beta  =  ggplot(data = as.data.frame(df[df$coefs2 == "beta", ]),  aes(x = coefs2,  y = value,  group = G)) +
  facet_grid(G ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(beta)) +
  ylab(quote(hat(beta))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
beta
X11()
Q  =  ggplot(data = as.data.frame(df[df$coefs2 == "Q", ]),  aes(x = coefs2,  y = value,  group = G)) +
  facet_grid(G ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(Q)) +
  ylab(quote(hat(Q))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
Q
X11()
sigma  =  ggplot(data = as.data.frame(df[df$coefs2 == "sigma", ]),   aes(x = coefs2,  y = value,  group = G)) +
  facet_grid(G ~ coefs,
             labeller = label_parsed,
             scales = "free_y") +
  geom_boxplot(aes(fill = G)) + labs(fill = "Groupe") +
  geom_hline(yintercept = 0) +
  xlab(quote(sigma)) +
  ylab(quote(hat(sigma))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
sigma
X11()
theta  =  ggplot(data = dftheta,  aes(  y = theta2e)) +
  geom_boxplot(fill='#00008B')  +
  geom_hline(yintercept = 0) +
  xlab(quote(theta)) +
  ylab(quote(hat(theta))) +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
theta
X11()
err  =  ggplot(data = dferr,  aes(  y = err)) +
  geom_boxplot(fill='#00008B')  +
  xlab("") +
  ylab("Accuracy") +
  theme_bw() + theme(axis.text.x = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position = "top")
err
ggarrange(pi,mu, beta,Q,sigma,theta, 
          ncol = 5,  nrow = 1,  common.legend = TRUE)
ggsave("figure1.pdf",  width = 8,  height = 10)