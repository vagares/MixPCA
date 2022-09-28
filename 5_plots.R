library(ggplot2)
library(ggpubr)
library(Ipaper)

df$prettycoefs = factor(df$coefs,  c("beta0", "beta1", "sigma0", "sigma1"),  c("hat(beta)[0]", "hat(beta)[1]", "hat(sigma)[0]^2", "hat(sigma)[1]^2"))
df$prettycont = df$cont
df$prettymethod = df$method


cols2  =  c("Scopt" = "#009E73")
lt2  =  c("Scopt" = "solid")
shaps2  =  c("Scopt" = 2)



x11()
beta_0  =  ggplot(data = as.data.frame(df[df$coefs == "beta0", ]),  aes(x = prettycont,  y = value,  group = prettycont)) +
  geom_boxplot(aes(fill = prettycont)) + labs(fill = "Contamination") +
  geom_hline(yintercept = 1) +
  xlab("Contamination") +
  ylab(quote(hat(beta)[0])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) +
  theme(legend.position = "top")
beta_1  =  ggplot(data = as.data.frame(df[df$coefs == "beta1", ]),  aes(x = prettycont,  y = value,  group = prettycont)) +
  geom_boxplot(aes(fill = prettycont)) + labs(fill = "Contamination") +
  geom_hline(yintercept = 1) +
  xlab("Contamination") +
  ylab(quote(hat(beta)[1])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) +
  theme(legend.position = "top")
sigma_0  =  ggplot(data = as.data.frame(df[df$coefs == "sigma0", ]),  aes(x = prettycont,  y = value,  group = prettycont)) +
  geom_boxplot(aes(fill = prettycont)) + labs(fill = "Contamination") +
  geom_hline(yintercept = 1) +
  xlab("Contamination") +
  ylab(quote(hat(sigma)[0])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) +
  theme(legend.position = "top")
sigma_1  =  ggplot(data = as.data.frame(df[df$coefs == "sigma1", ]),  aes(x = prettycont,  y = value,  group = prettycont)) +
  geom_boxplot(aes(fill = prettycont)) + labs(fill = "Contamination") +
  geom_hline(yintercept = 1) +
  xlab("Contamination") +
  ylab(quote(hat(sigma)[0])) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1)) +
  theme(legend.position = "top")

ggarrange(beta_0,  beta_1, sigma_0, sigma_1, 
          ncol = 4,  nrow = 1,  common.legend = TRUE)
ggsave("figure1.pdf",  width = 8,  height = 10)