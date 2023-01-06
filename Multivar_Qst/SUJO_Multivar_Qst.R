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
full_Fixed<-array(NA,dim=c(12,dims*2,1000),
                  dimnames = list(dimnames(Males_Fixed)[[1]],
                                  paste0(rep(c("f","m"),each=dims),
                                         dimnames(Males_Fixed)[[2]]),NULL))
full_Fixed[,1:dims,]<-Females_Fixed
full_Fixed[,dims+1:dims,]<-Males_Fixed


##########full multi-Qst---------
full_G<-full_G*5
multiQst_obs <-
  foreach(i=501:1000) %do% {
  B<-var(full_Fixed[,,i])
  G<-full_G[,,i-500]
  solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
}

median_MultiQst<-
  abind::abind(multiQst_obs,along=3) %>% apply(., 1:2, median)
median_B<-var(apply(full_Fixed,1:2,median))
median_G<-apply(full_G,1:2,median)

RWA<-eigen(median_MultiQst)
Vs<-RWA$vectors[,1:dims] 
colnames(Vs)<-paste0("Axis",1:(dims))
Vsf<-t(t(Vs) %*% 
        expm::sqrtm(median_B+2*median_G))

# MQst_full <-
#   foreach(i=1:1000) %do% {
#   B<-var(full_Fixed[,,i])
#   G<-full_G[,,i]
#   Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
#   Qstm
# }

diver.plot<-
  ldply(multiQst_obs, function(C) {
  data.frame(female=tr(C[1:dims,1:dims])/dims, 
             male=tr(C[(dims+1):(dims*2),(dims+1):(dims*2)])/dims)
    }) %>% melt %>% 
  ggplot(., aes(variable,value))+
  geom_boxplot()+
  xlab("")+
  ylab("divergence")+
  theme_bw()
diver.plot

ggsave("diver.pdf", diver.plot,width = 5,height = 4)

eigen_multiQst<-ldply(multiQst_obs, function(x) {
  x <- t(Vs) %*% x %*% (Vs)
  ev<-eigen(x)$values
  data.frame(M=mean(ev),V=t(ev))
  })

MQst_NULL <- foreach(i=501:1000,.combine = "rbind") %do% {
  G<-full_G[,,i-500]
  G<-t(Vs) %*% G %*% (Vs)
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-Re(eigen(Qstm)$value[1:dims])
  # Qstm_r<- t(RWA$vectors) %*% Qstm %*% RWA$vectors
  # data.frame(M=mean(diag(Qstm_r)),t(diag(Qstm_r)))
  data.frame(M=mean(ev),V=t(ev))
}

quantile(eigen_multiQst$M, c(0.025,0.975))
quantile(MQst_NULL$M, c(0.025,0.975))
quantile(afm$fst, c(0.025,0.975))

# PCQst_full<-ldply(MQst_full, function(Qstm) {
#   Qstm_r<- t(RWA$vectors) %*% Qstm %*% RWA$vectors
#   diag(Qstm_r)
#   })

confints<-ldply(MQst_NULL,quantile, probs=c(0.025,0.975)) %>% melt
colnames(confints)[1]<-"ID"
axis<-1:11
multiQst_FULL.plot<-
  eigen_multiQst %>%
  melt %>% 
  subset(., variable %in% c("M", paste0("V.",axis))) %>%
  ggplot(., aes(variable, value))+
  # geom_violin(draw_quantiles = c(0.025,0.975))+
  stat_summary(fun.data=confint95.stat, geom = "errorbar", width=.1)+
  stat_summary(fun="median",size = 2, geom = "point")+
  geom_line(aes(1:12, value),
            confints %>% subset(., ID %in% c("M", paste0("V.",axis)) & variable=="2.5%"),
             linetype=2)+
  geom_line(aes(1:12, value),
            confints %>% subset(., ID %in% c("M", paste0("V.",axis)) & variable=="97.5%"),
            linetype=2)+
  ylab("FSTq")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14))+
  scale_x_discrete(labels=c('M'=expression(bar(lambda)),
                            'V.1'=expression(lambda[1]),
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
###################

# PCalign_sample<-
#   foreach(i=1:1000,.combine = "rbind") %do% {
#     mqst<-MQst_full[[i]]
#     b<-var(full_Fixed[,,i])
#     g<-full_G[,,i]
#     RWA<-eigen(mqst)
#     Vsm<-RWA$vectors
#     colnames(Vsm)<-paste0("Axis",1:10)
#     Vsm<-t(Vsm) %*% expm::sqrtm(b+2*g)
#     t(apply(Vsm, 1, Normalize))
#   }
# # 
# PCalign_sample <- PCalign_sample %*% t(PCalign_sample)
# diag(PCalign_sample)
# 
# PCalign_dist<-sqrt(2*(1-abs(PCalign_sample)))
# 
# i<-!apply(PCalign_dist,1,function(x) any(is.na(x)))
# PCalign_dist<-PCalign_dist[i,i]
# PCoA<-cmdscale(na.omit(PCalign_dist))
# save(PCoA,file = "PCoA.Rdata")
# load("PCoA.Rdata")
# 
# PCoA.plot<-ggplot(data.frame(axis=rownames(PCoA),PCoA))+
#   geom_point(aes(X1,X2),alpha=0.1)+
#   geom_point(aes(X1,X2, color=axis),subset(data.frame(axis=rownames(PCoA),PCoA), 
#                                axis %in% c("Axis1","Axis2","Axis3")))+
#   coord_fixed()+
#   xlab("PCoA1")+
#   ylab("PCoA2")+
#   scale_color_manual(values=c("darkred","darkgreen","goldenrod"))
# ggsave("PCoA.pdf", PCoA.plot)

###############Alignment----

diffalign_sample_full<-
  foreach(i=1:500,.combine = "rbind") %do% {
  mqst<-multiQst_obs[[i]]
  b<-var(full_Fixed[,,i])
  g<-full_G[,,i]
  RWA<-eigen(mqst)
  Vsm<-Re(RWA$vectors[,1:dims])
  colnames(Vsm)<-paste0("Axis",1:(dims))
  Vsm<-t(t(Vsm)*RWA$values)
  Vsm_f<-(t(Vsm) %*% 
            expm::sqrtm(b+2*g))[,1:dims] 
  Vsm_f<-Vsm_f %*% t(shapePCA$rotation[,1:dims])
  Vsm_m<-(t(Vsm) %*% 
            expm::sqrtm(b+2*g))[,(dims+1):(dims*2)] 
  Vsm_m<-Vsm_m %*% t(shapePCA$rotation[,1:dims])
  
  data.frame(axis=1:(dims*2),
             size=Re(RWA$values),
             diff=Re(apply(Vsm_m-Vsm_f,1, Norm)),
             align=Re(diag(t(apply(Vsm_m,1,Normalize)) %*% apply(Vsm_f,1,Normalize))))
}

diffalign_sample<-diffalign_sample_full %>% subset(., size>0.25)
ctest<-cor.test(diffalign_sample$diff,fisherz(diffalign_sample$align))
diffalign_sample.plot<-
  ggplot(diffalign_sample,aes(diff, fisherz(align)))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_y_continuous(breaks=c(-2:3),
                     labels=round(fisherz2r(c(-2:3)),3))+
  xlab("Procrustes distance")+
  ylab("Vector correlation")+
  annotate("text",y = 2.5,x = 0.004, 
           label=paste0("cor=",round(ctest$estimate,3),", p< 2.2e-16"))
diffalign_sample.plot
ggsave("diffalign_sample.pdf", diffalign_sample.plot)

###################
# Vsm<-RWA$vectors
# colnames(Vsm)<-paste0("Axis",1:10)
# Vsm<-t(t(Vsm)*RWA$values)
# Vsm_f<-(t(Vsm) %*% 
#   expm::sqrtm(Full_m_B+2*Full_m_G))[,1:5] 
# Vsm_f<-Vsm_f %*% t(shapePCA$rotation[,1:5])
# Vsm_m<-(t(Vsm) %*% 
#   expm::sqrtm(Full_m_B+2*Full_m_G))[,6:10] 
# Vsm_m<-Vsm_m %*% t(shapePCA$rotation[,1:5])

# diffalign_mean.df<-data.frame(diff=apply(Vsm_m-Vsm_f,1, Norm),
#            align=diag(t(apply(Vsm_m,1,Normalize)) %*% apply(Vsm_f,1,Normalize)))
# 
# diffalign_mean.plot<-
#   ggplot(diffalign_mean.df,aes(diff, fisherz(align)))+
#   geom_smooth(se=F,method="lm",color="black", size=0.5, linetype=2)+
#   geom_point()+
#   geom_text_repel(aes(label=1:10))+
#   theme_classic()+
#   scale_y_continuous(breaks=fisherz(c(0.45,0.75,0.9,0.96)),
#                      labels=c(0.45,0.75,0.90,0.96))+
#   theme(legend.position = "none")+
#   xlab("Procrustes distance")+
#   ylab("Vector correlation")+
#   annotate("text",y = 1.75,x = 0.0015, label=paste0("cor=",
#                                                     round(cor(diffalign_mean.df)[1,2],3)))
# 
# ggsave("diffalign_mean.pdf", diffalign_mean.plot)



# Vsf<-Vsm_f*200
# Vsm<-Vsm_m*200
# Vsf<-Vsf*5
# Vsm<-Vsm*5

# round(diag(t(apply(Vsf, 1, Normalize)) %*% apply(Vsm, 1, Normalize)),3)
# 
# Qstvectors.plot<-
#   rbind(data.frame(sex="f",
#                    adply(Vsf,1,function(v) {
#                      data.frame(GPA$consensus,t(t(GPA$consensus)+v))})),
#         data.frame(sex="m",
#                    adply(Vsm,1,function(v) {
#                      data.frame(GPA$consensus,t(t(GPA$consensus)+v))}
#                    ))) %>%
#   # subset(., X1 %in% paste0("Axis",1:6)) %>%
#     ggplot(., aes(X,Y))+
#     geom_point()+
#     geom_segment(aes(x=X,y=Y,xend=X.1,yend=Y.1),
#                  arrow = arrow(length = unit(0.05, "npc")))+
#     geom_segment(aes(x=X,y=Y,xend=X.1,yend=Y.1),
#                  data.frame(GPA$consensus[links[,1],],
#                             GPA$consensus[links[,2],]), 
#                  alpha=0.5)+
#     coord_fixed()+
#     theme_bw()+
#     facet_grid(X1~sex)
# 
# ggsave("Qstvectors.pdf",Qstvectors.plot,width = 6,height = 15)


source('Piras et al Frontiers Supplementary material SECOND REVISION/Piras et al Frontiers ancillary functions.R')
Vsm<-(t(Vs) %*% 
        expm::sqrtm(median_B+2*median_G))[,(dims+1):(dims*2)] 
Vsm<-Vsm %*% t(shapePCA$rotation[,1:dims])

links<-alply(links, 1, identity)
pal<-"Blue-Red 2"
pdf("FSTq_axes.pdf",width = 16,height = 40)
layout(matrix(1:(dims*2),dims,2,byrow = T))
par(mar=c(0,0,0,0))
# zlim<-c(-1.5,1.5)
scale<-10
for(h in 1:dims){
  targ<-t(t(GPA$consensus)+(Vsf[h,]*scale))
  ref<-GPA$consensus
  psl2d(ref,targ,links,
        mag=1,
        ar=100,
        alpha=0.1,
        log=T,
        st=F,
        plotlegend=F,
        # zlim=zlim,
        interpfin=T,
        elli=F,
        elliax=F,plotgrid1 = F, plotgrid2 = F,plotboth=F,
        colors = hcl.colors(4,palette=pal), colorsleig=hcl.colors(4,palette=pal))
  lineplot(ref,links,col="gray",lwd = 1)
  lineplot(targ,links,col=1,lwd = 1.5)
  # gradientLegend2(zlim,depth=0.02,pos=0.6,notvals = 0,side=4,color=hcl.colors(4,palette=pal))
  
  targ<-t(t(GPA$consensus)+(Vsm[h,]*scale))
  ref<-GPA$consensus
  psl2d(ref,targ,links,
        mag=1,
        ar=100,
        alpha=0.1,
        log=T,
        st=F,
        plotlegend=F,
        # zlim=zlim,
        interpfin=T,
        elli=F,
        elliax=F,plotgrid1 = F, plotgrid2 = F,plotboth=F,
        colors = hcl.colors(4,palette=pal), colorsleig=hcl.colors(4,palette=pal))
  lineplot(ref,links,col="gray",lwd = 1)
  lineplot(targ,links,col=1,lwd = 1.5)
  # gradientLegend2(zlim,depth=0.02,pos=0.6,notvals = 0,side=4,color=hcl.colors(4,palette=pal))
}
dev.off()

#Sexos separados----
I<-grep(":M",colnames(full_G))
MQst_m <- foreach(i=1:500) %do% {
  B<-var(full_Fixed[,I,i])
  G<-full_G[I,I,i]
  solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
}

GQst_m<-laply(MQst_m, function(x) eigen(x)$values %>% mean)

MQst_N <- foreach(i=1:500,.combine = "rbind") %do% {
  G<-full_G[I,I,i]
  G<-t(Vs[I,]) %*% G %*% (Vs[I,])
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-eigen(Qstm)$values
  data.frame(M=mean(ev),V=t(ev))
}

I<-grep(":F",colnames(full_G))
MQst_f <- foreach(i=1:500) %do% {
  B<-var(full_Fixed[,I,i])
  G<-full_G[I,I,i]
  solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
}

GQst_f<-laply(MQst_f, function(x) eigen(x)$values %>% mean)

FQst_N <- foreach(i=501:1000,.combine = "rbind") %do% {
  G<-full_G[I,I,i-500]
  G<-t(Vs[I,]) %*% G %*% (Vs[I,])
  fst<-afm$fst[i]
  B<-(2*fst/(1-fst))*G
  x<-rmvnorm(12,sigma = B)
  B<-var(x)
  Qstm<-solve(expm::sqrtm(B+2*G)) %*% B %*% solve(expm::sqrtm(B+2*G))
  ev<-eigen(Qstm)$values
  data.frame(M=mean(ev),V=t(ev))
}

confintsF<-ldply(FQst_N,quantile, probs=c(0.025,0.975)) %>% melt %>% data.frame(sex="females")
names(confintsF)[1]<-"ID"

confintsM<-ldply(MQst_N,quantile, probs=c(0.025,0.975)) %>% melt %>% data.frame(sex="males")
names(confintsM)[1]<-"ID"

confints<-rbind(confintsF,confintsM)

# quantile(MQst_N,c(0.025,0.975))
# quantile(FQst_N,c(0.025,0.975))
quantile(afm$fst,c(0.025,0.975))

PCQst_m<-ldply(MQst_m, function(x) eigen(x)$values)
PCQst_f<-ldply(MQst_f, function(x) eigen(x)$values)

GQst.plot<-
  rbind(data.frame(M= GQst_m, PCQst_m, sex="males"),
        data.frame(M= GQst_f, PCQst_f, sex="females")) %>%
  melt %>% 
  subset(.,  variable %in% c("M", paste0("V", 1:dims))) %>%
  ggplot(., aes(variable, value))+
  facet_grid(~sex)+
  stat_summary(fun.data=confint95.stat, geom = "errorbar", width=.1)+
  stat_summary(fun="median",size = 2, geom = "point")+
  geom_line(aes(c(1:12,1:12), value),
            subset(confints, variable=="2.5%" & ID %in% c("M", paste0("V.", 1:11))),
            linetype=2)+
  geom_line(aes(c(1:12,1:12), value),
            subset(confints, variable=="97.5%" & ID %in% c("M", paste0("V.", 1:11))),
            linetype=2)+
  ylab("FSTq")+
  xlab("")+
  theme_bw()+
  scale_x_discrete(labels=c('M'=expression(bar(lambda)),
                            'V1'=expression(lambda[1]),
                            'V2'=expression(lambda[2]),
                            'V3'=expression(lambda[3]),
                            'V4'=expression(lambda[4]),
                            'V5'=expression(lambda[5]),
                            'V6'=expression(lambda[6]),
                            'V7'=expression(lambda[7]),
                            'V8'=expression(lambda[8]),
                            'V9'=expression(lambda[9]),
                            'V10'=expression(lambda[10]),
                            'V11'=expression(lambda[11])))
GQst.plot
ggsave("GQst_sep.pdf", GQst.plot, width = 16,height = 8, dpi=600, units="cm")

ggsave("multiQst_both.pdf",  
       plot_grid(multiQst_FULL.plot,diver.plot,ncol = 1,vjust = T,hjust = T),
       width = 5,height = 6)

save(file = "multiQst.Rdata", multiQst_obs)

#ploting males

pdf("FSTq_1axe_Males.pdf",width = 8,height = 5)
{targ<-t(t(GPA$consensus)+(Vsm[1,]*scale))
ref<-GPA$consensus
psl2d(ref,targ,links,
      mag=1,
      ar=100,
      alpha=0.1,
      log=T,
      st=F,
      plotlegend=F,
      # zlim=zlim,
      interpfin=T,
      elli=F,
      elliax=F,plotgrid1 = F, plotgrid2 = F,plotboth=F,
      colors = hcl.colors(4,palette=pal), colorsleig=hcl.colors(4,palette=pal))
# lineplot(ref,links,col="gray",lwd = 1)
lineplot(targ,links,col=1,lwd = 1.5)}
dev.off()