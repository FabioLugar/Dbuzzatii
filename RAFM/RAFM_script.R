library(readxl)
source("~/Dropbox/__R/__Config/rafm/RAFM.r")
microsat<-read_xlsx("Microsat.xlsx")
microsat<-as.matrix(microsat)
# microsat[microsat=="-1"]<-NA

afm<-do.all(microsat, 1010000, 10000, 1000)

save(afm,file = "afm.Rdata")


fst_mic<-read_xlsx("Matriz_MicFST.xlsx",col_names = TRUE)
fst_mic<-as.matrix(fst_mic)
rownames(fst_mic)<-fst_mic[,1]
fst_mic<-fst_mic[,-1]
fst_mic<-cbind(rbind(NA,fst_mic),NA)
rownames(fst_mic)[1]<-colnames(fst_mic)[1]
colnames(fst_mic)[dim(fst_mic)[1]]<-rownames(fst_mic)[dim(fst_mic)[1]]
fst_mic<-as.dist(fst_mic)
fst_mic<-as.matrix(fst_mic)
fst_mic[fst_mic<0]<-0


pdf("traces1000000.pdf")
for(i in 1:12) for (j in 1:12) if (i>j){
  plot(afm$theta[i,j,],type="l")
  abline(a = fst_mic[i,j], b=0, lty=2)
  title(paste("i=",i,"j=", j))
}
dev.off()

