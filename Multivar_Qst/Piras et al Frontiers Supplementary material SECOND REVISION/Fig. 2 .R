source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

quad2d<-cbind(c(-0.5,0.5,0.5,-0.5),c(-0.5,-0.5,0.5,0.5)) ###create a centered regular square

#First row-left panel
randquad2d<-randomrot(quad2d,maxz=20)[,,1]
tpsgridpaolo(quad2d,randquad2d,linksTT=conslinks(4,open=F),linksYY=conslinks(4,open=F),xlim=c(-0.8,0.8),lwdispl=4) ###this function hacks tpsgrid() from "shapes" package
#### end

### Shear + rotation
vecaff<-c(0.12132,-0.307107,2,0.12132)
thismag=0.2
defo<-trapgen(quad2d,a11=vecaff[1],a12=vecaff[3],a21=vecaff[2],a22=vecaff[4]) #### function that creates affine and optionally non affine transformations
tensor<-defo$A+diag(2)########  this is the tensor
fin<-defo$newc
##### end 

#### evaluate gradient tensors
der1<-tpsjacpsl2d(quad2d,fin,jitter(quad2d,factor=0.00000000000001),doopa=F) ### this function evaluates F, C and second order gradients at specific points; when using the same landmarks source, an infinitiesimal (virtually 0) jitter is necessary for numerical reasons
###end

##### This piece of code extracts necessary objects and compute PSD for plotting
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
######### end

#Fig.2 first row-left panel
tpsgridpaolo(quad2d,fin,linksTT = conslinks(4,open=F),ext=0.3,lwdispl=4,linksYY=conslinks(4,open=F),ngrid=20,colgrid="grey",xbegin=-1, ybegin=-1,xwidth = 2) 
dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
for(i in 1:4){
  centro<-fin[i,]
  defelliadd(dat,centro,matf[[i]],col=makeTransparent("white",alpha=0),mag=thismag,border=1)
  finpsl1<-pslinfin1(der1$cjac1psl[i,1:2],der1$jac1[[i]])
  finpsl2<-pslinfin1(der1$cjac1psl[i,3:4],der1$jac1[[i]])
  arrows(centro[1],centro[2],(centro[1]+(finpsl1[1]*thismag)),(centro[2]+(finpsl1[2]*thismag)),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(finpsl1[1]*thismag)),(centro[2]-(finpsl1[2]*thismag)),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]+(finpsl2[1]*thismag)),(centro[2]+(finpsl2[2]*thismag)),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(finpsl2[1]*thismag)),(centro[2]-(finpsl2[2]*thismag)),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
}

#Fig.2  second row-left panel
tpsgridpaolo(quad2d,fin,linksTT = conslinks(4,open=F),ext=0.3,lwdispl=4,linksYY=conslinks(4,open=F),ngrid=20,colgrid="grey",xbegin=-1, ybegin=-1,xwidth = 2)
dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
for(i in 1:4){
  centro<-quad2d[i,]
  defelliadd(dat,centro,matc[[i]],col=makeTransparent("white",alpha=0),mag=thismag,border=1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]+(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]-(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]+(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]-(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
}


### The same as above but on shapes aligned via OPA 
der1<-tpsjacpsl2d(quad2d,rotonto(quad2d,fin)$yrot,jitter(quad2d,factor=0.00000000000001),doopa=F) ### this function evaluate F, C and second order gradient at specific points; when using the same landmarks source an infinitiesimal (virtually 0) jitter is necessary for computational reasons
###end
##### This piece of code extracts necessary objects and compute PSD for plotting
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
######### end

#Fig.2  second row-right panel
tpsgridpaolo(quad2d,rotonto(quad2d,fin)$yrot,linksTT = conslinks(4,open=F),ext=0.3,lwdispl=4,linksYY=conslinks(4,open=F),ngrid=20,colgrid="grey",xbegin=-1, ybegin=-1,xwidth = 2)
dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
for(i in 1:4){
  centro<-quad2d[i,]
  defelliadd(dat,centro,matc[[i]],col=makeTransparent("white",alpha=0),mag=thismag,border=1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]+(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,1]*sqrt(der1$cjacps[i,1]))*thismag),(centro[2]-(der1$cjac1psl[i,2]*sqrt( der1$cjacps[i,1]))*thismag),col=ifelse(signs[i,1]>0,1,"violet"),length=0.1,code=codes[i,1],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]+(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]+(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
  arrows(centro[1],centro[2],(centro[1]-(der1$cjac1psl[i,3]*sqrt(der1$cjacps[i,2]))*thismag),(centro[2]-(der1$cjac1psl[i,4]*sqrt( der1$cjacps[i,2]))*thismag),col=ifelse(signs[i,2]>0,1,"violet"),length=0.1,code=codes[i,2],lwd=0.1)
}
##########  end

###### Create two generic irregular polygons
irrindef<-centershapes(matrix(c(-49.2941176470588, -74.2941176470588, -97.2941176470588, -65.2941176470588, -100.294117647059, -123.294117647059, -68.2941176470588, 74.7058823529412, 110.705882352941, 46.7058823529412, 81.7058823529412, 119.705882352941, 89.7058823529412, 45.7058823529412, 20.7058823529412, 3.70588235294122, -15.2941176470588, 125.882352941176, 65.8823529411765, -4.11764705882354, -38.1176470588235, -75.1176470588235, -110.117647058824, -148.117647058824, -129.117647058824, -61.1176470588235, -56.1176470588235, -20.1176470588235, 23.8823529411765, 56.8823529411765, 84.8823529411765, 53.8823529411765, 97.8823529411765, 132.882352941176),ncol=2))[,,1]
defindef<-centershapes(matrix(c(-124.764705882353, -108.764705882353, -136.764705882353, -103.764705882353, -22.7647058823529, 84.2352941176471, 166.235294117647, 188.235294117647, 134.235294117647, 103.235294117647, 61.2352941176471, 51.2352941176471, 1.23529411764707, -33.7647058823529, -46.7647058823529, -88.7647058823529, -123.764705882353, 77.3529411764706, 30.3529411764706, -49.6470588235294, -120.647058823529, -158.647058823529, -129.647058823529, -69.6470588235294, -13.6470588235294, 13.3529411764706, -27.6470588235294, -29.6470588235294, 53.3529411764706, 36.3529411764706, 52.3529411764706, 113.352941176471, 115.352941176471, 107.352941176471),ncol=2))[,,1]
##end

#Fig. 2 third row-left panel
plotmyarrays(shapes::abind(irrindef,defindef),links=conslinks(17,open=F),txt=T,xlim=c(-200,200),ylim=c(-200,200),pch=c(19,19),cextext = 0.7)
#end
#Fig. 2 third row-right panel; the function mopa() align and center a target on a source using OPA or MOPA
plotmyarrays(shapes::abind(irrindef,mopa(irrindef,defindef,rot="opa")$opizzata[,,1]),links=conslinks(17,open=F),txt=T,xlim=c(-200,200),ylim=c(-200,200),pch=c(19,19),cextext = 0.7)
par(new=T)
plotmyarrays(mopa(irrindef,defindef,rot="mopa")$opizzata[,,1],links=conslinks(17,open=F),txt=T,xlim=c(-200,200),ylim=c(-200,200),pch=c(19),cextext = 0.7,col=4)
##end


