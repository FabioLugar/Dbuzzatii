source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

#Create a regular square
quad2d<-cbind(c(-0.5,0.5,0.5,-0.5),c(-0.5,-0.5,0.5,0.5))
##end


### Deform it into a parallelogram 
vecaff<-c(0.12132,-0.307107,2,0.12132)
thismag=0.2
defo<-trapgen(quad2d,a11=vecaff[1],a12=vecaff[3],a21=vecaff[2],a22=vecaff[4])
tensor<-defo$A+diag(2)########  this is the tensor
fin<-defo$newc
der1<-tpsjacpsl2d(quad2d,fin,centroids(quad2d),doopa=F) ### evaluate tensors at square' centroid
###end

## extract  F, C, PSD and C' eigenvalues and eigenvectors for plotting
sign1<-sign(der1$cjacps[1]-1)
code1<-ifelse(sign1>0,2,1)
sign2<-sign(der1$cjacps[2]-1)
code2<-ifelse(sign2>0,2,1)
col1<-ifelse(sign1>0,1,6)
col2<-ifelse(sign2>0,1,6)
matt<-der1[[1]][[1]]
mats<-expm::sqrtm(der1$ccjac[[1]])
mats2<-matrix(0,ncol=3,nrow=3)
mats2[1:2,1:2]<-mats
matt2<-matrix(0,ncol=3,nrow=3)
matt2[1:2,1:2]<-matt
### end

##Plot things
par(mfrow=c(1,3))
plotmyarrays(quad2d,links=c(conslinks(4,open=F)),col=1,xlim=c(-1.8,1.8))
plotmyarrays(fin,links=c(conslinks(4,open=F)),col=2,xlim=c(-1.8,1.8))
plotmyarrays(shapes::abind(quad2d,fin),links=c(conslinks(4,open=F)),xlim=c(-1.8,1.8))
dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
centro<-c(0,0)
defelliadd(dat,centro,mats2,col=makeTransparent("white",alpha=0),mag=thismag,border=1)
defelliadd(dat,centro,matt2,col=makeTransparent("white",alpha=0),mag=thismag,border=2)
arrows(0,0,(0+der1$cjac1psl[1]*sqrt(der1$cjacps[1]))*thismag,(0+der1$cjac1psl[2]*sqrt( der1$cjacps[1]))*thismag,col=col1,length=0.1,code=code1,lwd=0.1)
arrows(0,0,(0-der1$cjac1psl[1]*sqrt(der1$cjacps[1]))*thismag,(0-der1$cjac1psl[2]*sqrt( der1$cjacps[1]))*thismag,col=col1,length=0.1,code=code1,lwd=0.1)
arrows(0,0,(0+der1$cjac1psl[3]*sqrt(der1$cjacps[2]))*thismag,(0+der1$cjac1psl[4]*sqrt( der1$cjacps[2]))*thismag,col=col2,length=0.1,code=code2,lwd=0.1)
arrows(0,0,(0-der1$cjac1psl[3]*sqrt(der1$cjacps[2]))*thismag,(0-der1$cjac1psl[4]*sqrt( der1$cjacps[2]))*thismag,col=col2,length=0.1,code=code2,lwd=0.1)
finpsl1<-pslinfin1(der1$cjac1psl[1:2],der1$jac1[[1]])#### this function rotate PSD to plot them on the target configuration
finpsl2<-pslinfin1(der1$cjac1psl[3:4],der1$jac1[[1]])#### this function rotate PSD to plot them on the target configuration
arrows(0,0,(0+finpsl1[1])*thismag,(0+finpsl1[2])*thismag,col=col1,length=0.1,code=code1,lwd=0.1)
arrows(0,0,(0-finpsl1[1])*thismag,(0-finpsl1[2])*thismag,col=col1,length=0.1,code=code1,lwd=0.1)
arrows(0,0,(0+finpsl2[1])*thismag,(0+finpsl2[2])*thismag,col=col2,length=0.1,code=code2,lwd=0.1)
arrows(0,0,(0-finpsl2[1])*thismag,(0-finpsl2[2])*thismag,col=col2,length=0.1,code=code2,lwd=0.1)

