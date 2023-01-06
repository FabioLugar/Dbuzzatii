setwd("/Users/fmachado/Google Drive/Alas Drosophila/FULL DIM WD/Skewers/")
library(plyr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(evolqg)
library(doParallel)
library(cowplot)
library(foreach)
library(psych)
library(abind)

load("../Gcalc/Gs_11_weakPrior.Rdata")
load("../GPA/GPA.Rdata")


dims<-dim(Females_Fixed)[2]
betas<-apply(matrix(rnorm((dims*2)*1000),(dims*2),1000),2,Normalize)

apply(full_G,1:2,median) %>% write.csv(.,file = "G_median.csv")
apply(full_G,1:2,median) %>% cov2cor %>% write.csv(.,file = "G_median_cor.csv")

cors<-adply(full_G,3,function(x) diag(cov2cor(x)[1:dims,1:dims+dims]))
round(quantile(unlist(cors[,-1]), prob=c(0.025,0.975),na.rm = TRUE),3)

cors<-adply(full_G,3,function(x) {
  C<-cov2cor(x)[1:dims,1:dims+dims]
  diag(C)<-NA
  C
})
round(quantile(unlist(cors[,-1]), prob=c(0.025,0.975),na.rm = TRUE),3)


R_rand<-foreach(i=1:10000,.combine = "c") %do% {
  I<-sample(1:500,1)
  G<-full_G[,,I]
  ev<-Evolvability(G,betas[,I])
  G[1:dims,1:dims+dims]<-0
  G[1:dims+dims,1:dims]<-0
  ev_0<-Evolvability(G,betas[,I])
  ev/ev_0
}

R_obs<-foreach(i=1:1000,.combine = "c") %do% {
  I<-sample(1:500,1)
  G<-full_G[,,I]
  G<-ExtendMatrix(G,ret.dim = 10)$ExtMat
  b<-solve(G) %*% t(cbind(Females_Fixed[,,I],Males_Fixed[,,I]))
  b<-apply(b,2,Normalize)
  ev<-Evolvability(G,b)
  G[1:dims,1:dims+dims]<-0
  G[1:dims+dims,1:dims]<-0
  ev_0<-Evolvability(G,b)
  ev/ev_0
}

R.plot<-
  rbind(data.frame(R=R_rand,type="rand"),
        data.frame(R=R_obs,type="obs")) %>%
  ggplot(.,aes(R,fill=type))+
  geom_density(adjust = 2,alpha=0.5, show.legend = F,)+
  scale_fill_manual(values=c("darkorchid4", "black"))+
  geom_vline(xintercept=1, linetype=2)+
  theme_bw()
  
ggsave("R.pdf" ,R.plot,width = 8,height = 5, dpi=600, units="cm")  
 # Vsm_f %*% t(shapePCA$rotation[,1:5])