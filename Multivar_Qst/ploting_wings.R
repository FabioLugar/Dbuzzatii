#ploting wings----
source('~/Google Drive/Alas Drosophila/FULL DIM WD/Multivar_Qst/Piras et al Frontiers Supplementary material SECOND REVISION/Piras et al Frontiers ancillary functions.R')


links<-alply(links, 1, identity)

Vsm<-(t(Vs) %*% 
        expm::sqrtm(B_median+2*G_median))[,1:dims+dims] 
Vsm<-Vsm[,1] %*% t(shapePCA$rotation[,1:dims])

Vsf<-(t(Vs) %*% 
        expm::sqrtm(B_median+2*G_median))[,1:dims] 
Vsf<-Vsf[,1] %*% t(shapePCA$rotation[,1:dims])

pal<-"Blue-Red"
# pal<-"Dark 3"

Vsm<-(t(Vs) %*% 
        expm::sqrtm(B_median+2*G_median))[,1:dims+dims]
Vsf<-(t(Vs) %*% 
        expm::sqrtm(B_median+2*G_median))[,1:dims] 

colMeans(apply(Males_Fixed, 1:2, median) %*% t(Vsm))*1000000

Vsm<-Vsm%*% t(shapePCA$rotation[,1:dims])
Vsf<-Vsf%*% t(shapePCA$rotation[,1:dims])

data.frame(cor=diag(t(apply(Vsm, 1, Normalize)) %*% apply(Vsf, 1, Normalize))[1:6],
           ldply(1:6,function(i){
             ang<-angleTest(Vsf[i,], Vsm[i,])
             data.frame(angle=ang$angle*(180/pi))
           }))


Vsf[c(2,4,5),]<-Vsf[c(2,4,5),]*-1
Vsm[c(2,4,5),]<-Vsm[c(2,4,5),]*-1
naxis<-6

source('~/Dropbox/__R/Rfun/_Morf_Geo/tps2d.Rfun.r')
source('~/Dropbox/__R/Rfun/_Morf_Geo/tps.Rfun.r')
library(sp)
{
  pal<-"Blue-Red 2"
  pal<-"Blue-Red"
  pdf("FSTq_axesTPS.pdf",width = 16,height = 30)
  
  layout(matrix(1:(naxis*2),naxis,2,byrow = T))
  par(mar=c(0,0,0,0))
  # zlim<-c(-1.5,1.5)
  sc<-20
  ref<-GPA$consensus
  for(h in 1:naxis){
    targ<-GPA$consensus
    ref<-GPA$consensus + matrix(Vsf[h,]*sc,ncol = 2, byrow = T)
    sFE1<-spsample(Polygon(ref[c(1,3,7,13,15,14,12,4),]),10000,type="regular")
    sR1<-sFE1@coords
    sT1<-tps2d(sR1,ref,targ)
    def<-sqrt(apply((sT1-sR1)^2,1,sum))
    xl<-length(unique(sR1[,1]))
    yl<-length(unique(sR1[,2]))
    im<-matrix(NA,xl,yl)
    xind<-(1:xl)[as.factor(rank(sR1[,1]))]
    yind<-(1:yl)[as.factor(rank(sR1[,2]))]
    n<-length(xind)
    for (i in 1:n){im[xind[i], yind[i]]<-def[i]}
    par(mar=c(0,0,0,0))
    image(sort(unique(sR1[,1])),sort(unique(sR1[,2])),im,
          col=hcl.colors(100,palette=pal),asp=T,axes=F,frame=F, xlab="", ylab="", 
          xlim=range(targ[,1]*1.2),
          ylim=range(targ[,2]*1.2))
    lineplot(ref,links,col=1,lwd = 1.5)
    
    ref<-GPA$consensus + matrix(Vsm[h,]*sc,ncol = 2, byrow = T)
    sFE1<-spsample(Polygon(ref[c(1,3,7,13,15,14,12,4),]),10000,type="regular")
    sR1<-sFE1@coords
    sT1<-tps2d(sR1,ref,targ)
    def<-sqrt(apply((sT1-sR1)^2,1,sum))
    xl<-length(unique(sR1[,1]))
    yl<-length(unique(sR1[,2]))
    im<-matrix(NA,xl,yl)
    xind<-(1:xl)[as.factor(rank(sR1[,1]))]
    yind<-(1:yl)[as.factor(rank(sR1[,2]))]
    n<-length(xind)
    for (i in 1:n){im[xind[i], yind[i]]<-def[i]}
    par(mar=c(0,0,0,0))
    image(sort(unique(sR1[,1])),sort(unique(sR1[,2])),im,
          col=hcl.colors(100,palette=pal),asp=T,axes=F,frame=F, xlab="", ylab="", 
          xlim=range(targ[,1]*1.2),
          ylim=range(targ[,2]*1.2))
    lineplot(ref,links,col=1,lwd = 1.5)
    
  }
  dev.off()
}

{
  pal<-"Blue-Red 2"
  pdf("FSTq_axes.pdf",width = 16,height = 30)
  naxis<-6
  layout(matrix(1:(naxis*2),naxis,2,byrow = T))
  par(mar=c(0,0,0,0))
  sc<-20
  ref<-GPA$consensus
  ar<-175
  for(h in 1:naxis){
    targ<-(GPA$consensus + matrix(Vsf[h,]*sc,ncol = 2, byrow = T))
    psl2d(ref,targ,links,
          mag=1,
          ar=ar,
          alpha=1,
          log=T,
          st=F,
          plotlegend=F,
          # zlim=zlim,
          interpfin=T,
          elli=F,
          elliax=F,plotgrid1 = F, plotgrid2 = F,plotboth=F,
          colors = hcl.colors(50,palette=pal), colorsleig=hcl.colors(50,palette=pal))
    # lineplot(ref,links,col="gray",lwd = 1)
    lineplot(targ,links,col=1,lwd = 1.5)
    
    targ<-(GPA$consensus + matrix(Vsm[h,]*sc,ncol = 2, byrow = T))
    psl2d(ref,targ,links,
          mag=1,
          ar=ar,
          alpha=1,
          log=T,
          st=F,
          plotlegend=F,
          # zlim=zlim,
          interpfin=T,
          elli=F,
          elliax=F,plotgrid1 = F, plotgrid2 = F,plotboth=F,
          colors = hcl.colors(50,palette=pal), colorsleig=hcl.colors(50,palette=pal))
    # lineplot(ref,links,col="gray",lwd = 1)
    lineplot(targ,links,col=1,lwd = 1.5)
  }
  dev.off()
}

{
  pdf("FSTq_axes_LAND.pdf",width = 16,height = 30)
  naxis<-6
  layout(matrix(1:(naxis*2),naxis,2,byrow = T))
  par(mar=c(0,0,0,0))
  sc<-20
  ref<-GPA$consensus
  for(h in 1:naxis){
    targ<-(GPA$consensus + matrix(Vsf[h,]*sc,ncol = 2, byrow = T))
    plotRefToTarget(ref,targ, method="vector", 
                    gridPars=gridPar(tar.pt.bg = "black", link.col="black", link.lwd=0,pt.bg="black", pt.size = 0.5))
    lineplot(ref,links,lwd = 1)
    
    targ<-(GPA$consensus + matrix(Vsm[h,]*sc,ncol = 2, byrow = T))
    plotRefToTarget(ref,targ, method="vector", 
                    gridPars=gridPar(tar.pt.bg = "black", link.col="black", link.lwd=0,pt.bg="black", pt.size = 0.5))
    lineplot(ref,links,lwd = 1)
  }
  dev.off()
}



