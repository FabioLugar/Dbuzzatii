library(geomorph)

land<-readland.tps("COMBINADAS_LineasIguales.TPS", specID ="imageID")
options(contrasts=c("contr.sum","contr.poly"))


#En "classifier" est?n todos los individuos medidos de todas las lineas geneticas
classifier <- read.csv("Clas_TANDA_TODASv2_LineasIguales.CSV", header=T, sep=";")

#Procrustes analysis------------------------------------------------------------
GPA<-gpagen(land)

shapePCA<-gm.prcomp(GPA$coords)
summary(shapePCA)

rm(land)
save.image("GPA.Rdata")
