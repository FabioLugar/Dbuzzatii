setwd("/Users/fmachado/Google Drive/Alas Drosophila/FULL DIM WD/Multivar_Qst/")
library(geomorph)
library(psych)
library(reshape2)
library(plyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(evolqg)
library(mvtnorm)
library(doParallel)
library(ggbeeswarm)

confint95.stat <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.975))
  names(r) <- c("ymin","ymax")
  r
}

load("../Gcalc/Gs_11_weakPrior.Rdata"); load("../RAFM/afm.Rdata"); load("../GPA/GPA.Rdata")
afm$theta<-afm$theta[,,501:1000]
afm$fst<-afm$fst[501:1000]
afm$kappa<-afm$kappa[,,501:1000]
afm$alpha<-afm$alpha[,501:1000]

links <- matrix(
  c(1,  2,
    2,  3,
    3,  6,
    6,  7,
    7, 13,
    13,15,
    15,14,
    14,12,
    12, 4,
    4,  1,
    2,  4,
    2,  5,
    5, 12,
    5,  8,
    8,  9,
    9, 10,
    10,11,
    11,13,
    8, 14,
    10,15,
    11, 6,
    3,  9,
    1,  3,
    3,  7),
  24,2, byrow = T)

dims<-dim(Males_Fixed)[2]
full_Fixed<-array(NA,dim=c(12,dims*2,500),
                  dimnames = list(dimnames(Males_Fixed)[[1]],
                                  paste0(rep(c("f","m"),each=dims),
                                         dimnames(Males_Fixed)[[2]]),NULL))
full_Fixed[,1:dims,]<-Females_Fixed
full_Fixed[,dims+1:dims,]<-Males_Fixed



# full multi-Qst ----------------------------------------------------------
iF<-0.250/2
scale<- 1/(2*iF)

full_GHi<-full_G*scale
# full_GHi<-full_G
# for(i in 1:dim(full_GHi)[3]) full_GHi[,,i]<-full_G[,,i]* runif(1,2,5)

B_median<-var(apply(full_Fixed,1:2,median))
G_median<-apply(full_GHi,1:2,median)

Im<-grep(":M",colnames(G_median))
If<-grep(":F",colnames(full_GHi))

var.plot<-
  rbind(
    data.frame(ExtendMatrix(G_median[Im,Im])$var, sex="Males"),
    data.frame(ExtendMatrix(G_median[If,If])$var, sex="Females"),
    data.frame(ExtendMatrix(G_median)$var, sex="Both")) %>%
  ggplot(., aes(PC,Gradient.Variance))+
  geom_line(aes(color=sex))+
  geom_point(aes(color=sex))+
  scale_color_brewer(name="", type = "qual",palette = 2)+
  theme_bw()+
  annotate("segment", x = 10.5, xend = 10, 
           y = 2e-04, yend = 2e-05,
           arrow = arrow(angle = 30, length = unit(.2,"cm")))


ggsave("gradVar_Gs.pdf", var.plot, height=5,width = 6)
#full multi Qst
FSTq_dist <-
  foreach(i=1:500) %do% {
    B<-var(full_Fixed[,,i])
    G<-full_GHi[,,i]
    G<-ExtendMatrix(G,ret.dim = 10)$ExtMat
    # G<-ExtendMatrix(G,var.cut.off = 2.207600e-06)$ExtMat
    solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  }
#median FSTq
FSTq_median<-
  abind::abind(FSTq_dist,along=3) %>% apply(., 1:2, median)


#getting eigenvectors of the median to stablish as base line
RWA<-eigen(FSTq_median)
Vs<-RWA$vectors[,1:dims] 
colnames(Vs)<-paste0("Axis",1:(dims))
Vsf<-t(t(Vs) %*% 
         expm::sqrtm(B_median+2*G_median))

#calculating divergence based on FSTq as mean eigenvalue
FSTq_distHi<-FSTq_dist

#obtaining eigenvalues for FSTq
FSTq_eigen<-ldply(FSTq_dist, function(x) {
  x <- t(Vs) %*% x %*% (Vs) # rotating sample Qst to the median FSTq
  ev<-eigen(x)$values
  data.frame(M=mean(ev),V=t(ev))
})

#obtaining eigenvalues for FSTq
FSTq_NULL <- foreach(i=1:500,.combine = "rbind") %do% {
  G<-full_GHi[,,i]
  G<-t(Vs) %*% G %*% (Vs) # rotating G to the median FSTq
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  G<-ExtendMatrix(G,ret.dim = 10)$ExtMat
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-eigen(Qstm)$value[1:dims]
  data.frame(M=mean(ev),V=t(ev))
}

confints<-ldply(FSTq_NULL,quantile, probs=c(0.025,0.975)) %>% melt
colnames(confints)[1]<-"ID"
axis<-1:11
FSTq_eigenHi<-FSTq_eigen
confintsHi<-confints

##Plots-------------------------------
diver.plot<-
  ldply(FSTq_distHi, function(C) {
    data.frame(female=mean(eigen(C[1:dims,1:dims])$values),
               male=mean(eigen(C[1:dims+dims,1:dims+dims])$values))
  }) %>% melt %>%
  ggplot(., aes(variable,value))+
  geom_boxplot()+
  xlab("")+
  ylab("divergence")+
  theme_bw()
diver.plot
ggsave("diver.pdf", diver.plot,width = 5,height = 4)

multiQst_FULL.plot<-
  FSTq_eigenHi %>% 
  melt%>% 
  subset(., variable %in% paste0("V.",axis)) %>%
  ggplot(., aes(variable, value))+
  stat_summary(fun.data=confint95.stat, geom = "errorbar", width=.1,position = position_dodge(0.5))+
  stat_summary(fun="median",size = 2, geom = "point",position = position_dodge(0.5))+
  geom_line(aes(1:11, value),
            confintsHi  %>%
              subset(., ID %in% paste0("V.",axis) & variable=="2.5%"),
            linetype=2)+
  geom_line(aes(1:11, value),
            confintsHi %>% 
              subset(., ID %in% paste0("V.",axis) & variable=="97.5%"),
            linetype=2)+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14))+
  scale_x_discrete(labels=c('V.1'=expression(lambda[1]),
                            'V.2'=expression(lambda[2]),
                            'V.3'=expression(lambda[3]),
                            'V.4'=expression(lambda[4]),
                            'V.5'=expression(lambda[5]),
                            'V.6'=expression(lambda[6]),
                            'V.7'=expression(lambda[7]),
                            'V.8'=expression(lambda[8]),
                            'V.9'=expression(lambda[9]),
                            'V.10'=expression(lambda[10]),
                            'V.11'=expression(lambda[11])))
multiQst_FULL.plot
ggsave("multiQst_FULL.pdf", multiQst_FULL.plot,width = 5,height = 4)


#Sexos separados------------------------
##GHi----
I<-grep(":M",colnames(full_GHi))
FSTq_m <- foreach(i=1:500) %do% {
  B<-var(full_Fixed[,I,i])
  G<-full_GHi[I,I,i]
  solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
}


FSTq_eigen_m<-ldply(FSTq_m, function(x) {
  x <- t(Vs[I,]) %*% x %*% (Vs[I,]) # rotating sample Qst to the median FSTq
  ev<-eigen(x)$values
  data.frame(M=mean(ev),V=t(ev))
})

FSTq_m_NULL <- foreach(i=1:500,.combine = "rbind") %do% {
  G<-full_GHi[I,I,i]
  G<-t(Vs[I,]) %*% G %*% (Vs[I,])
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-eigen(Qstm)$values
  data.frame(M=mean(ev),V=t(ev))
}

I<-grep(":F",colnames(full_GHi))
FSTq_f <- foreach(i=1:500) %do% {
  B<-var(full_Fixed[,I,i])
  G<-full_GHi[I,I,i]
  solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
}

FSTq_eigen_f<-ldply(FSTq_m, function(x) {
  x <- t(Vs[I,]) %*% x %*% (Vs[I,]) # rotating sample Qst to the median FSTq
  ev<-eigen(x)$values
  data.frame(M=mean(ev),V=t(ev))
})

FSTq_f_NULL <- foreach(i=1:500,.combine = "rbind") %do% {
  G<-full_GHi[I,I,i]
  G<-t(Vs[I,]) %*% G %*% (Vs[I,])
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-eigen(Qstm)$values
  data.frame(M=mean(ev),V=t(ev))
}

confintsF<-ldply(FSTq_f_NULL,quantile, probs=c(0.025,0.975)) %>% melt %>% data.frame(sex="female")
names(confintsF)[1]<-"ID"

confintsM<-ldply(FSTq_m_NULL,quantile, probs=c(0.025,0.975)) %>% melt %>% data.frame(sex="male")
names(confintsM)[1]<-"ID"

confints<-rbind(confintsF,confintsM)

PCQst_m<-ldply(FSTq_m, function(x) eigen(x)$values)
PCQst_f<-ldply(FSTq_f, function(x) eigen(x)$values)

FSTq_eigen_mHi<-FSTq_eigen_m
FSTq_eigen_fHi<-FSTq_eigen_f
confintsHi<-confints

##Plots--------
GQst.plot<-
  rbind(data.frame(sex="male",FSTq_eigen_mHi),
        data.frame(sex="female",FSTq_eigen_fHi)) %>%
  select(., -M, -paste0("V.", 7:dims)) %>%
  melt %>% 
  ggplot(., aes(variable, value))+
  facet_grid(~sex)+
  stat_summary(fun.data=confint95.stat, geom = "errorbar", width=.1,position = position_dodge(0.5))+
  stat_summary(fun="mean",size = 2, geom = "point",position = position_dodge(0.5))+
  # geom_line(aes(c(1:11,1:11), value),
  #           subset(confintsHi, variable=="2.5%" & ID %in% paste0("V.", 1:11)),
  #           linetype=2)+
  # geom_line(aes(c(1:11,1:11), value),
  #           subset(confintsHi, variable=="97.5%" & ID %in% paste0("V.", 1:11)),
  #           linetype=2)+
  geom_line(aes(c(1:6,1:6), value),
            subset(confintsHi, variable=="2.5%" & ID %in% paste0("V.", 1:6)),
            linetype=2)+
  geom_line(aes(c(1:6,1:6), value),
            subset(confintsHi, variable=="97.5%" & ID %in% paste0("V.", 1:6)),
            linetype=2)+
  ylab("FSTq")+
  xlab("")+
  theme_bw()+
  scale_x_discrete(labels=c('V.1'=expression(lambda[1]),
                            'V.2'=expression(lambda[2]),
                            'V.3'=expression(lambda[3]),
                            'V.4'=expression(lambda[4]),
                            'V.5'=expression(lambda[5]),
                            'V.6'=expression(lambda[6]),
                            'V.7'=expression(lambda[7]),
                            'V.8'=expression(lambda[8]),
                            'V.9'=expression(lambda[9]),
                            'V.10'=expression(lambda[10]),
                            'V.11'=expression(lambda[11])))
GQst.plot
ggsave("FSTq_sep.pdf", GQst.plot, width = 16,height = 8, dpi=600, units="cm")

ggsave("multiQst_both.pdf",  
       plot_grid(multiQst_FULL.plot,
                 GQst.plot,
                 ncol = 1,vjust = T,hjust = T), 
       width = 5, height = 6)

save(FSTq_eigen, FSTq_eigen, FSTq_dist, file = "multiQst.Rdata")

