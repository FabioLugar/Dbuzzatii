source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

#### Create (trivially) the trapezoid that is the bilinear non affine deformation of the regular square
quad2d<-cbind(c(-0.5,0.5,0.5,-0.5),c(-0.5,-0.5,0.5,0.5)) ###create a centered regular square
quad3<-cbind(c(-1,1,0.3,-0.3),c(-0.25,-0.25,0.25,0.25))
##end
der1<-tpsjacpsl2d(quad2d,quad3,jitter(quad2d,factor=0.00000000000001),doopa=F)### evaluate tensors

signs<-NULL
codes<-NULL
princ<-NULL
princeigen<-NULL
for(i in 1:nrow(der1$cjacps)){
  signi<-NULL
  codi<-NULL
  for(k in 1:ncol((der1$cjacps))){
    signik<-sign(der1$cjacps[i,k]-1)
    codik<-ifelse(signik>0,2,1)
    signi<-c(signi,signik)
    codi<-c(codi,codik)
  }
  princi<-which.max(abs(1-der1$cjacps[i,]))
  signs<-rbind(signs,signi)
  codes<-rbind(codes,codi)
  princeigeni<-der1$cjacps[i,princi]
  princeigen<-c(princeigen,princeigeni)
  princ<-c(princ,princi)
}
matc<-NULL
matf<-NULL
for(i in 1:length(der1$jac1)){
  matfi<-cbind(rbind(der1$jac1[[1]],0),0)
  matci<-cbind(rbind(expm::sqrtm(der1$ccjac[[i]]),0),0)
  matc<-c(matc,list(matci))
  matf<-c(matf,list(matfi))
}
dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
par(mfrow=c(1,2))
###
thismag=0.2
tpsgridpaolo(quad2d,quad3,linksTT = conslinks(4,open=F),ext=0.3,lwdispl=4,linksYY=conslinks(4,open=F),ngrid=20,colgrid="grey",xbegin=-1, ybegin=-1,xwidth = 2)
for(i in 1:4){
  centro<-quad2d[i,]
  defelliadd(dat,centro,matc[[i]],col=makeTransparent("white",alpha=0),mag=thismag,border=1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]+(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]-(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]+(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]-(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
}


tpsgridpaolo(quad2d,quad3,linksTT = conslinks(4,open=F),ext=0.3,lwdispl=4,linksYY=conslinks(4,open=F),ngrid=20,colgrid="grey",xbegin=-1, ybegin=-1,xwidth = 2)
for(i in 1:4){
  centro<-quad3[i,]
  defelliadd(dat,centro,matf[[i]],col=makeTransparent("white",alpha=0),mag=thismag,border=1)
  finpsl1<-pslinfin1(der1$cjac1psl[i,1:2],der1$jac1[[i]])#### this function rotate PSD to plot them on the target configuration
  finpsl2<-pslinfin1(der1$cjac1psl[i,3:4],der1$jac1[[i]])
  arrows(centro[1],centro[2],(centro[1]+(finpsl1[1]*thismag)),(centro[2]+(finpsl1[2]*thismag)),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(finpsl1[1]*thismag)),(centro[2]-(finpsl1[2]*thismag)),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]+(finpsl2[1]*thismag)),(centro[2]+(finpsl2[2]*thismag)),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(finpsl2[1]*thismag)),(centro[2]-(finpsl2[2]*thismag)),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
}








