source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

#### Create a triangle in 3d
sou<-centershapes(matrix(c(-0.456852841433138, -0.419287296338007, 0.495183987533674, 0.57672172951512, -0.391928912168369, 0.524729123292491, -2.42947177005, -2.40142457499169, -1.49298409295268),ncol=3))[,,1]
##end

##Create a deformed state using a randomly chosen tensor F using Y=XF'
mytensor<-matrix(c(0.560417981,-0.161174118,0.499914806,1.994749570,1.122686033,0.004679548,0.569923566,-0.148280066,0.492044566),ncol=3)
tar<-(sou%*%t(mytensor))
##end

prova<-tpsjacpsl3d(sou,tar,doopa=F,centroids(tar))#### evaluate tensor using TPS 
provadry<-tpsdry2(sou,tar,doopa=F)#### the same using an alternative function
psltri1ob<-psltri1(sou,tar) #### evaluate tensor using the discrete way

fdiscr<-matrix(psltri1(sou,tar)[27:35],ncol=3)### this is F found in the discrete way
fdiscr/mytensor###### the discrete way correctly recovers initial tensor

ftps<-prova$jac1[[1]]### this is F found using TPS

cdiscr<-matrix(psltri1(sou,tar)[12:20],ncol=3)### this is C found in the discrete way
ctps<-prova$ccjac[[1]]### this is C found using TPS
#####  as explained in the main text fdiscr and ftps do not coincide

########## moreover...the eigen-decomposition performed on cdiscr returns only two non-zero eigenvalues while we have full rank for ctps:
eigen(cdiscr)$values
eigen(ctps)$values

provadry$at### this is equal to ftps TRANSPOSED

##########This piece of code finds the projectors for C, projectors for F, normals and contravariant basis

a1s<-sou[2,]-centroids(sou)
a2s<-sou[3,]-centroids(sou)
a3s<-CrossProduct3D(a1s,a2s)
a3sn<-a3s/tensorA::norm(a3s,type="F")
Ps<-diag(3)-a3sn%*%t(a3sn)
tprojs<-t(Ps)%*%prova$ccjac[[1]]%*%Ps
matrix(psltri1(sou,tar)[12:20],ncol=3)-tprojs
a1t<-tar[2,]-centroids(tar)
a2t<-tar[3,]-centroids(tar)
a3t<-CrossProduct3D(a1t,a2t)
a3tn<-a3t/tensorA::norm(a3t,type="F")
Pt<-diag(3)-a3tn%*%t(a3tn)
tprojt<-t(Pt)%*%prova$jac1[[1]]%*%Ps
tprojt<-Pt%*%prova$jac1[[1]]%*%t(Ps)
matrix(psltri1(sou,tar)[27:35],ncol=3)-tprojt





##### build Figure 7
linkscube<-list(c(1,2),c(2,4),c(4,3),c(3,1),c(5,6),c(6,8),c(8,7),c(7,5),c(6,2),c(8,4),c(7,3),c(5,1))
mfrow3d(3,2,sharedMouse = T)
lineplot(sou,conslinks(3,open=F),col=1)
text3d(sou,texts=c(1:3),col=1)
lineplot(tar,conslinks(3,open=F),col=2)
text3d(tar,texts=c(1:3),col=2)
decorate3d()
next3d()
lineplot(sou,conslinks(3,open=F),col=1)
text3d(sou,texts=c(1:3),col=1)
lineplot(tar,conslinks(3,open=F),col=2)
text3d(tar,texts=c(1:3),col=2)
decorate3d()
plotpsltri1(sou,tar,psltri1ob,mag=0.2)
plotpsltri1(sou,tar,psltri1ob,usec=F,mag=0.2,colelli=2)
next3d()
lineplot(sou,conslinks(3,open=F),col=1)
text3d(sou,texts=c(1:3),col=1)
lineplot(tar,conslinks(3,open=F),col=2)
text3d(tar,texts=c(1:3),col=2)
decorate3d()
shade3d(traslamesh(scalemesh(defosph(creasph(),prova$jac1[[1]]),0.2,center="none"),centroids(tar)),alpha=0.2,col=2)
shade3d(traslamesh(scalemesh(defosph(creasph(),expm::sqrtm(prova$ccjac[[1]])),0.2,center="none"),centroids(sou)),alpha=0.2,col=1)
arrows3D(centroids(sou),centroids(sou)+a3sn)## add normal
arrows3D(centroids(tar),centroids(tar)+a3tn,col=2)## add normal
next3d()
lineplot(sou,conslinks(3,open=F),col=1)
text3d(sou,texts=c(1:3),col=1)
lineplot(tar,conslinks(3,open=F),col=2)
text3d(tar,texts=c(1:3),col=2)
decorate3d()
shade3d(traslamesh(scalemesh(defosph(creasph(),prova$jac1[[1]]),0.2,center="none"),centroids(tar)),alpha=0.2,col=2)
shade3d(traslamesh(scalemesh(defosph(creasph(),expm::sqrtm(prova$ccjac[[1]])),0.2,center="none"),centroids(sou)),alpha=0.2,col=1)
arrows3D(centroids(sou),centroids(sou)+a3sn)## add normal
arrows3D(centroids(tar),centroids(tar)+a3tn,col=2)## add normal
plotpsltri1(sou,tar,psltri1ob,mag=0.2)
plotpsltri1(sou,tar,psltri1ob,usec=F,mag=0.2,colelli=2)
next3d()
twoshapes1<-shapes::abind(t(cube3d()$vb[1:3,]),t(cube3d()$vb[1:3,])%*%provadry$at)
twoshapes2<-shapes::abind(t(cube3d()$vb[1:3,]),t(cube3d()$vb[1:3,])%*%t(matrix(psltri1(sou,tar)[27:35],ncol=3)))
ranx=c(-6,6)
rany=ranx
ranz=ranx
plotmyarrays(abind::abind(twoshapes1*5,twoshapes2*5),col="white",txt=F,cex=0)
plotmyarrays(abind::abind(twoshapes1,twoshapes2)[,,1:2],links=linkscube,cex=0,xlim=ranx,ylim=rany,zlim=ranz)
decorate3d()
next3d()
plotmyarrays(abind::abind(twoshapes1*5,twoshapes2*5),col="white",txt=F,cex=0)
plotmyarrays(abind::abind(twoshapes1,twoshapes2)[,,3:4],links=linkscube,cex=0,xlim=ranx,ylim=rany,zlim=ranz)
decorate3d()

