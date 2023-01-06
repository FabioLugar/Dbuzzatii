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

par(mfrow=c(2,2),mar=c(1,1,1,1))

obi1<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,plotelliax1 = F,plotelliax2 = F,elli=F,interpinit=T,interpfin=T,nadd=5000,zlim=c(-1,1))
lineplot(obi1$init,links,lwd=3)
lineplot(obi1$fin,col=2,links,lwd=3)
gradientLegend(valRange=c(-1,1), side=2, inside=FALSE,color=c("blue","cyan","yellow","red"),pos=0.5,depth = 0.05,length=0.5)

obi1<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,plotelliax1 = F,plotelliax2 = F,elli=F,interpinit=T,interpfin=T,nadd=5000,se=T,zlim=c(-15,2))
lineplot(obi1$init,links,lwd=3)
lineplot(obi1$fin,col=2,links,lwd=3)
gradientLegend(valRange=c(-15,2), side=2, inside=FALSE,color=c("blue","cyan","yellow","red"),pos=0.5,depth = 0.05,length=0.5)

obi1<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,plotelliax1 = F,plotelliax2 = F,elli=F,interpinit=T,interpfin=T,nadd=5000,bebs=T,zlim=c(-3,11))
lineplot(obi1$init,links,lwd=3)
lineplot(obi1$fin,col=2,links,lwd=3)
gradientLegend(valRange=c(-3,11), side=2, inside=FALSE,color=c("blue","cyan","yellow","red"),pos=0.5,depth = 0.05,length=0.5)

