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
#### Build Fig. 10 first row
par(mfrow=c(1,2))
obi0<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,zlim=c(-1,1),elli=T,st=T)
par(mfg=c(1,1))
lineplot(obi0$init,links,lwd=3)
lineplot(obi0$fin,col=2,links,lwd=3)
gradientLegend(valRange=c(-1,1), side=2, inside=FALSE,color=c("blue","cyan","yellow","red"),pos=0.3,fit.margin =F,length=0.5)

par(mfg=c(1,2))
lineplot(obi0$init,links,lwd=3)
lineplot(obi0$fin,col=2,links,lwd=3)
gradientLegend(valRange=c(-1,1), side=2, inside=FALSE,color=c("blue","cyan","yellow","red"),pos=0.3,fit.margin =F,length=0.5)

#### Build Fig. 10 second row
obi0<-psl2d(ret2d,ret2dna,mag=0.07,plotgrid1=T,plotgrid2=T,ar=180,plotlegend = F,log=T,zlim=c(-1,1),elli=F,st=T,elliax=T)
par(mfg=c(1,1))
lineplot(obi0$init,links,lwd=3)
lineplot(obi0$fin,col=2,links,lwd=3)
par(mfg=c(1,2))
lineplot(obi0$init,links,lwd=3)
lineplot(obi0$fin,col=2,links,lwd=3)
