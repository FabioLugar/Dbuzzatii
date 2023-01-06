source("YOUR PATH/Piras et al Frontiers ancillary functions.R")  #### type in your path to the source file

### Create a regular square grid
ret2d<-coordinates(GridTopology(c(0,0), c(1/3,1/3), c(7,7)))
links=list(c(1,43),c(43,49),c(49,7),c(7,6),c(6,5),c(5,4),c(4,3),c(3,2),c(2,1))
###end

##Create manually a (non affinely) deformed state
ret2dna<-ret2d
ret2dna[c(2,9,16,23,30,37,44),2]<-ret2dna[c(2,9,16,23,30,37,44),2]*0.65
ret2dna[c(5,12,19,26,33,40,40,47),2]<-ret2dna[c(5,12,19,26,33,40,40,47),2]*1.5
####end

##### plot grid
tpsgridpaolo(ret2d,rotonto(ret2d,ret2dna)$yrot,ngrid=20,displ=F,cex=0,ylim=c(-0.5,3),linksTT=links,linksYY=links)
dev.off()

###Build Fig. 9 first row
par(mfrow=c(1,2))
plotmyarrays(shapes::abind(ret2d,ret2dna),ind=1,links=links,col=c(1,makeTransparent("white",alpha=0)),pch=19,,xlim=c(0,2),cex=c(0.6),axes=F,cextext=0.8)
axis(2)
axis(1)
plotmyarrays(shapes::abind(ret2d,ret2dna),ind=2,links=links,col=c(1),pch=19,cex=(0.6),axes=F,cextext=0.8,xlim=c(0,2))
axis(2)
axis(1)
########  end

#### Create Fig. 9 second row
#### pls2d() evaluate tensors, homologous grids, PSD and corresponding deformation ellipses in bot source and target shapes
#### here only grids are shown
obi0<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,bebs=T,interpinit=F,interpfin=F,elli=F,st=T,alpha=1,plotelliax1=F,plotelliax2=F,plotri2=T,plotri1=T,links=links)
par(mfg=c(1,1))
points(obi0$centrosinit,pch=19,cex=0.3)##add triangles centroids
par(mfg=c(1,2))
points(obi0$centrosfin,pch=19,cex=0.3,col=2)##add triangles centroids
