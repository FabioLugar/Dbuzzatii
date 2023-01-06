setwd("~/Google Drive/Alas Drosophila/FULL DIM WD/DIMORPH/")
library(geomorph)
library(doParallel)
library(plyr)
library(dplyr)
library(evolqg)
library(ggplot2)
library(psych)
library(cowplot)
load("../GPA/GPA.Rdata")
registerDoParallel(cores = 8)
gmdf<-geomorph.data.frame(shape=GPA$coords, classifier)

# popManova<-procD.lm(shape~Poblacion+Poblacion:Sexo, SS.type = "II", data = gmdf)
# anova(popManova)
popManova<-procD.lm(shape~Poblacion*Sexo,iter = 10000, SS.type = "II", data = gmdf,turbo = T,Parallel = T)
anova(popManova)
write.csv(anova(popManova)$table,file = "manova.csv")

dimorph.lm<-
  foreach(i=unique(classifier$Poblacion)) %do% {
    I<-classifier$Poblacion==i
    gmdf<-geomorph.data.frame(shape=GPA$coords[,,I], classifier[I,])
    procD.lm(shape~Sexo, iter = 10000, data = gmdf,turbo = T,Parallel = T)
  }

laply(dimorph.lm, function(x) anova(x)$table[1,"Z"]) %>% mean
laply(dimorph.lm, function(x) anova(x)$table[1,"Pr(>F)"])
#Alignment----
load("../Gcalc/Gs_11_weakPrior.Rdata")

diffalign_fixed<-
  foreach(i=1:500) %do% {
    data.frame(diff=apply(Females_Fixed[,,i]-Males_Fixed[,,i],1,Norm),
               align=diag(t(apply(Females_Fixed[,,i],1,Normalize)) %*% apply(Males_Fixed[,,i],1,Normalize)))
  }

diffalign<-
  foreach(i=1:500,.combine = "rbind") %do% {
    fShape<-Females_Fixed[,,i] %*% t(shapePCA$rotation[,1:11])
    mShape<-Males_Fixed[,,i] %*% t(shapePCA$rotation[,1:11])
    
    data.frame(axis=rownames(fShape),
               diff=apply(fShape-mShape,1,Norm),
               align=diag(t(apply(mShape,1,Normalize)) %*% apply(fShape,1,Normalize)))
  }

diffalign_cor<-
  diffalign_fixed %>% laply(., function(x){
    cor(x[,1],fisherz(x[,2]))
  })


diffalign_cor.plot<-
  ggplot(data.frame(x=1,cor=diffalign_cor), aes(x,cor))+
  geom_violin(fill="black", alpha=0.2,draw_quantiles = c(0.025,0.975))+
  theme_bw()+geom_hline(yintercept=0, linetype=2)+
  ylim(-1,1)+
  theme(axis.title.x=element_text(color="white"),
        axis.text.x=element_text(color="white"),
        axis.ticks.x=element_line(color="white"))


diffalign.plot<-
  ggplot(diffalign, aes(diff, fisherz(align)))+
  geom_point(color="gray66")+
  theme_bw()+
  scale_y_continuous(breaks=fisherz(c(0.99,0.95,0.7,0,-0.7,-0.95, -0.99)),
                     labels=c(0.99,0.95,0.7,0,-0.7,-0.95, -0.99))+
  xlab("Procrustes distance")+
  ylab("Vector correlation")#+
# geom_point(aes(diff, fisherz(align)),diffalign_FSTq)+
# ggrepel::geom_label_repel(aes(diff, fisherz(align), label=axis),diffalign_FSTq)

ggsave("diffalign.pdf", plot_grid(diffalign.plot,diffalign_cor.plot,rel_widths = c(1,0.5),labels = c("A.","B."), label_fontfamily = "Times"),width = 13,height = 7, dpi=600, units="cm")  
