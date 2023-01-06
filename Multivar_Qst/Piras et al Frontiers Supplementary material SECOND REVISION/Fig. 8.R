source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

####Create a cube in 3d
config<-centershapes(t(cube3d()$vb[-4,])/2)[,,1]
###end

#### Create a tensor
tensor3d<-matrix(c(1.12132,-0.307107,0,2,1.12132,0,0,0,1),ncol=3)
###end

###Create a deformed state of the regular cube
fin3d<-config%*%t(tensor3d)
##end

der13d<-tpsjacpsl3d(config,fin3d,centroids(fin3d),doopa=F)###evaluate gradient tensors at target centroid
linkscube<-list(c(1,2),c(2,4),c(4,3),c(3,1),c(5,6),c(6,8),c(8,7),c(7,5),c(6,2),c(8,4),c(7,3),c(5,1))### links for a cube
plotmyarrays(shapes::abind(config,fin3d),links=linkscube)### plot shapes


####This function evaluates gradient tensors using discrete way by exploiting a homologous triangulation structure
cubtria<-matrix(c(1, 1, 2, 2, 4, 4, 3, 3, 1, 1, 5, 5, 2, 2, 6, 4, 8, 3, 7, 1, 2, 2, 7, 6, 5, 6, 8, 8, 7, 7, 5, 5, 3, 4, 6, 8),ncol=3)
obcub3d<-psltris(config,fin3d,cubtria,doopa=F)

open3d()
####Build Fig. 8 first row
mfrow3d(1,3,sharedMouse = T)
#plot3d(rbind(config,fin3d),axes=F,aspect=F,cex=0,col="white",xlab="",ylab="",zlab="")
crux(rbind(config,fin3d))
lineplot(config,linkscube)
text3d(config,texts=c(1:8))
next3d()
#plot3d(rbind(config,fin3d),axes=F,aspect=F,cex=0,col="white",xlab="",ylab="",zlab="")
crux(rbind(config,fin3d))
lineplot(fin3d,linkscube,col=2)
text3d(fin3d,texts=c(1:8),col=2)
next3d()
crux(rbind(config,fin3d))
#### this function plot a psltris object and plots both ellipses and ellipsoids with a bulk of options
plotpsltris(obcub3d,mag=0.2,usec=T,tpspoints=centroids(obcub3d$init),drawpoints=centroids(obcub3d$init),drawpelli=centrotri(obcub3d$ini,cubtria),plotelli=F,elliax=F,alphasph = 0.2,colsph=1,onlyprinc3d = F)
plotpsltris(obcub3d,mag=0.2,usec=F,tpspoints=centroids(obcub3d$init),drawpoints=centroids(obcub3d$init),drawpelli=centrotri(obcub3d$ini,cubtria),plotelli = F,elliax = F,colsph=2,alphasph = 0.2,onlyprinc3d = F)
lineplot(config,linkscube)
text3d(config,texts=c(1:8))
lineplot(fin3d,linkscube,col=2)
text3d(fin3d,texts=c(1:8),col=2)

open3d()
### Build Fig. 8 second row
quali<-which(duplicated((unlist(lapply(obcub3d$ff,sum))))==F) #### tensori unici
quads<-rbind(c(1,2,5,6),c(2,6,4,8),c(4,8,7,3),c(3,7,5,1),c(1,3,4,2),c(5,7,8,6))
centrotri(obcub3d$init,t(obcub3d$triang))
centrofacce<-centrotri(obcub3d$init,quads)
centrofaccefin<-centrotri(obcub3d$fin,quads)
indi<-c(1,3,5,7,9,11)
ffacce<-lapply(obcub3d$ff[indi],function(x) round(x,2))######## these are the projections of F on target faces
ccfacce<-lapply(lapply(lapply(obcub3d$ff[indi],function(x) round(x,2)),function(x) t(x)%*%x),eigen)
mfrow3d(1,2,sharedMouse = T)
plot3d(rbind(config,fin3d),size=0,type="s",axes=F,box=F,aspect=F,xlab="",ylab="",zlab="")
lineplot(config,linkscube)
text3d(config,texts=c(1:8))
mag=0.2
### this loop plots ellipses and PSD on source faces
for(i in 1: nrow(centrofacce)){
  defmat<-expm::sqrtm(obcub3d$cc[indi][[i]])
  mytrinit<-obcub3d$init[quads[i,],]
  mytri<-obcub3d$fin[quads[i,],]
  circ3d<-centershapes(cbind(circle2(radius=1,plot=F),0))[,,1]*mag
  circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytrinit,156))),mytrinit[1:4,],mytrinit[1,]),circ3d,reflection=T)
  forma<-centershapes(circ3d2$yrot)[,,1]
  tcirc3d<-t(defmat%*%t(forma))
  lines3d(tcirc3d+rep.row(centrofacce[i,],nrow(circ3d)),col=1)
  if(ccfacce[[i]]$values[1]>1){col1=1;code1=2}else{col1="violet";code1=1}
  if(ccfacce[[i]]$values[2]>1){col2=1;code2=2}else{col2="violet";code2=1}
  mag1<-sqrt(ccfacce[[i]]$values[1])
  mag2<-sqrt(ccfacce[[i]]$values[2])
  psl1i<-ccfacce[[i]]$vectors[,1]*mag1*mag
  psl2i<-ccfacce[[i]]$vectors[,2]*mag2*mag
  compositions::arrows3D(centrofacce[i,],centrofacce[i,]+psl1i,add=T,code=code1,col=col1)
  compositions::arrows3D(centrofacce[i,],centrofacce[i,]-psl1i,add=T,code=code1,col=col1)
  compositions::arrows3D(centrofacce[i,],centrofacce[i,]+psl2i,add=T,code=code2,col=col2)
  compositions::arrows3D(centrofacce[i,],centrofacce[i,]-psl2i,add=T,code=code2,col=col2)
}
plotpsltris(obcub3d,tpspoints=centroids(obcub3d$fin),drawpoints=centroids(obcub3d$fin),plotelli=F,elliax=F,alphasph=0.2,plotsph=T,tpspsl = T,mag=0.2,usec=T,colsph=1,onlyprinc3d = F)

next3d()
plot3d(rbind(config,fin3d),size=0,type="s",axes=F,box=F,aspect=F,xlab="",ylab="",zlab="")
lineplot(fin3d,linkscube,col=2)
text3d(fin3d,texts=c(1:8))
mag=0.2
indi<-c(1,3,5,7,9,11)
### this loop plots ellipses and PSD on target faces
for(i in 1: nrow(centrofaccefin)){
  defmat<-obcub3d$ff[indi][[i]]
  mytri<-obcub3d$init[quads[i,],]
  mytrinit<-obcub3d$fin[quads[i,],]
  circ3d<-centershapes(cbind(circle2(radius=1,plot=F),0))[,,1]*mag
  circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytri,156))),mytri[1:4,],mytri[1,]),circ3d,reflection=T)
  forma<-centershapes(circ3d2$yrot)[,,1]
  tcirc3d<-t(defmat%*%t(forma))
  lines3d(tcirc3d+rep.row(centrofaccefin[i,],nrow(circ3d)),col=2)
  if(ccfacce[[i]]$values[1]>1){col1=1;code1=2}else{col1="violet";code1=1}
  if(ccfacce[[i]]$values[2]>1){col2=1;code2=2}else{col2="violet";code2=1}
  psl1i<-c(pslinfin1(ccfacce[[i]]$vectors[,1],ffacce[[i]])*mag)
  psl2i<-c(pslinfin1(ccfacce[[i]]$vectors[,2],ffacce[[i]])*mag)
  compositions::arrows3D(centrofaccefin[i,],centrofaccefin[i,]+psl1i,add=T,code=code1,col=col1)
  compositions::arrows3D(centrofaccefin[i,],centrofaccefin[i,]-psl1i,add=T,code=code1,col=col1)
  compositions::arrows3D(centrofaccefin[i,],centrofaccefin[i,]+psl2i,add=T,code=code2,col=col2)
  compositions::arrows3D(centrofaccefin[i,],centrofaccefin[i,]-psl2i,add=T,code=code2,col=col2)
}
plotpsltris(obcub3d,tpspoints=centroids(obcub3d$fin),drawpoints=centroids(obcub3d$fin),plotelli=F,elliax=F,alphasph=0.2,plotsph=T,tpspsl = T,mag=0.2,usec=F,colsph=2,onlyprinc3d = F)
