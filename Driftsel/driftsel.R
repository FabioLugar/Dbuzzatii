setwd("/Users/fmachado/Google Drive/Alas Drosophila/FULL DIM WD/Driftsel/")
library(plyr)
library(ggplot2)
library(dplyr)
library(foreach)
library(evolqg)
source("~/Dropbox/__R/__Config/driftsel/driftsel.r")
load("../Gcalc/Gs_11_weakPrior.Rdata")
load("../RAFM/afm.Rdata")
load("../Multivar_Qst/multiQst.Rdata")
load("../GPA/GPA.Rdata")
# for(i in 1:500) full_G[,,i]<-full_G[,,i]*runif(1,2,5)
sample_rg<-501:1000
afm$theta<-afm$theta[,,sample_rg]
afm$kappa<-afm$kappa[,,sample_rg]
afm$fst<-afm$fst[sample_rg]
afm$alpha<-afm$alpha[,sample_rg]
# full_G<-full_G[,,sample_rg]
# Females_Fixed<-Females_Fixed[,,sample_rg]
# Males_Fixed<-Males_Fixed[,,sample_rg]
dims<-dim(Males_Fixed)[2]
full_Fixed<-array(NA,dim=c(12,dims*2,500),
                  dimnames = list(dimnames(Males_Fixed)[[1]],
                                  paste0(rep(c("f","m"),each=dims),
                                         dimnames(Males_Fixed)[[2]]),NULL))
full_Fixed[,1:dims,]<-Females_Fixed
full_Fixed[,1:dims+dims,]<-Males_Fixed
full_Gext<-full_G
iF<-0.25/2
scale<-1/(2*iF)
for(i in 1:500) full_Gext[,,i]<-ExtendMatrix(full_G[,,i]*scale, ret.dim = 10)$ExtMat
shape_Stest_Full<-S.test(full_Fixed,
                         full_Gext, 
                         afm$theta, silent = F)

shape_Stest_Female<-S.test(Females_Fixed, 
                           full_G[1:dims+dims,1:dims+dims,]*scale,
                           afm$theta, silent = F)

shape_Stest_Male  <-S.test(Males_Fixed, 
                           full_G[1:dims,1:dims,]*scale,
                           afm$theta, silent = F)

svalue.plot<-
  data.frame(full=shape_Stest_Full, 
             male=shape_Stest_Male, 
             female=shape_Stest_Female) %>%
  reshape2::melt() %>%
  ggplot(., aes(variable, value))+
  # geom_violin()+
  geom_hline(yintercept=c(0.025,0.975),linetype=2)+
  geom_jitter()+
  stat_summary(aes(variable, value),fun="mean",color="gray")+
  theme_bw()+
  xlab("")+
  ylab("S value")
svalue.plot
ggsave("svalue.pdf",svalue.plot,width = 5,height = 5)

## individual axis

FSTq_median<-
  abind::abind(FSTq_distHi,along=3) %>% apply(., 1:2, median)
B_median<-var(apply(full_Fixed,1:2,median))
G_median<-apply(full_G*scale,1:2,median)

#getting eigenvectors of the median to stablish as base line
RWA<-eigen(FSTq_median)
Vs<-RWA$vectors[,1:dims] 
colnames(Vs)<-paste0("Axis",1:(dims))
Vsf<-t(t(Vs) %*% 
         expm::sqrtm(B_median+2*G_median))


pc_Stest<-foreach(i=1:11) %do% {
  S.test(
    array(aaply(full_Fixed, 3, {
    function(x) x %*% Vsf[,i]
    }),
    dim = c(12,1,500)),
    array(aaply(full_Gext*scale, 3, {
      function(x) t(Vsf[,i]) %*% x %*% Vsf[,i]
    }),dim = c(1,1,500)),
    afm$theta, silent = F)
}

laply(pc_Stest, mean)

pc_Stest<-foreach(i=1:11) %do% {
  S.test(full_Fixed[,c(i,i+dims),],
         full_Gext[c(i,i+dims),c(i,i+dims),]*scale,
         afm$theta, silent = F)
}



pc_Stestf<-foreach(i=1:11) %do% {
  S.test(array(aaply(full_Fixed[,1:dims,], 3, {
    function(x) x %*% Vsf[1:dims,i]
  }) %>% t,
  dim = c(12,1,500)),
  array(aaply(full_G[1:dims,1:dims,]*scale, 3, {
    function(x) t(Vsf[1:dims,i]) %*% x %*% Vsf[1:dims,i]
  }),dim = c(1,1,500)),
  afm$theta, silent = F)
}

pc_Stestm<-foreach(i=1:11) %do% {
  S.test(array(aaply(full_Fixed[,1:dims+dims,], 3, {
    function(x) x %*% Vsf[1:dims+dims,i]
  }) %>% t,
  dim = c(12,1,500)),
  array(aaply(full_G[1:dims+dims,1:dims+dims,]*scale, 3, {
    function(x) t(Vsf[1:dims+dims,i]) %*% x %*% Vsf[1:dims+dims,i]
  }),dim = c(1,1,500)),
  afm$theta, silent = F)
}

laply(pc_Stestf, mean)
laply(pc_Stestm, mean)
