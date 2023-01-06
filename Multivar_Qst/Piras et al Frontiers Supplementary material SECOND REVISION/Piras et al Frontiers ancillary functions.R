   #### install/load required packages ####

required<-c("Morpho","rgl","geomorph","calibrate","MASS","compositions","bayesm","tensorA",
            "robustbase","shapes","manipulate","Arothron","gdata","Rvcg","vegan","phytools",
            "autoimage","Biobase","nat","tensr","circular","expm","plotfunctions","sp","fields",
            "matrixcalc","ape","geometry","pracma","tripack","geoR","gstat","RTriangle", "ggm")

if(any(!required%in%installed.packages()[,1])){
  if("Biobase"%in%required[which(!required%in%installed.packages()[,1])]) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      BiocManager::install("Biobase")
  }
  install.packages(required[which(!required%in%installed.packages()[,1])])
} 
sapply(required,require,character.only=TRUE)

#### Functions ####

randomrot<-function(array,minx=0,miny=0,minz=0,maxx=360,maxy=360,maxz=360){
  if(class(array)=="matrix"){array<-array(array,dim=c(dim(array)[1],dim(array)[2],1))}
  anglessz<-runif(dim(array)[3],min=minz,max=maxz)
  anglessy<-runif(dim(array)[3],min=miny,max=maxy)
  anglessx<-runif(dim(array)[3],min=minx,max=maxx)
  arraymod<-NULL
  if(dim(array)[2]<3){
    for(i in 1:dim(array)[3]){
      arraymodi<-rotate2(array[,,i],0,0,anglessz[i]) #### rotate2() is a wrapper of rotate()
      arraymod<-abind::abind(arraymod,arraymodi)
    }
  }else{
    for(i in 1:dim(array)[3]){
      arraymodi<-rotate2(array[,,i],anglessx[i],anglessy[i],anglessz[i]) #### rotate2() is a wrapper of rotate()
      arraymod<-abind::abind(arraymod,arraymodi)
    }
    
  }
  return(arraymod)
}

rotate2<-function (x, thetax, thetay, thetaz){
  
  if(class(x)=="matrix"){x<-array(x,dim=c(nrow(x),ncol(x),1))}else{x<-x}
  
  thetax <- thetax/180 * pi
  thetay <- thetay/180 * pi
  thetaz <- thetaz/180 * pi
  Rx <- matrix(c(1, 0, 0, 0, cos(thetax), sin(thetax), 0, -sin(thetax), 
                 cos(thetax)), 3, 3)
  Ry <- matrix(c(cos(thetay), 0, sin(thetay), 0, 1, 0, -sin(thetay), 
                 0, cos(thetay)), 3, 3)
  Rz <- matrix(c(cos(thetaz), sin(thetaz), 0, -sin(thetaz), 
                 cos(thetaz), 0, 0, 0, 1), 3, 3)
  y <- NULL
  n <- dim(x)[3]
  for (i in 1:n) {
    
    if(dim(x)[2]<3){yi <- cbind(x[,  ,i],rep(0,dim(x)[1])) %*% Rx %*% Ry %*% Rz}else{yi <- x[,  ,i] %*% Rx %*% Ry %*% Rz}
    
    
    y2<-array(yi,dim=c(dim(x)[1],dim(x)[2],1))
    y<-abind::abind(y,y2)
  }
  y
}

tpsgridpaolo<-function (TT, YY, xbegin = -999, ybegin = -999, xlim=NULL,ylim=NULL,xwidth = -999, opt = 1, ext = 0.0, asp=1,ngrid = 22, cex = 1, pch = 20,colshift=1,zslice = 0, graphics=T,mag = 1, axes3 = FALSE,linksTT=NULL,linksYY=NULL,collinksTT=1,collinksYY=2,lwdtt=2,lwdyy=2,colgrid=1,axes2d=T,collandsTT=collinksTT,collandsYY=collinksYY,displ=T,lwdispl=1,coldispl=4){
  ######### some SMALL changes from Ian Dryden's function from package "shapes"
  k <- dim(TT)[1]
  m <- dim(TT)[2]
  YY <- TT + (YY - TT) * mag
  bb <- array(TT, c(dim(TT), 1))
  aa <- defplotsize2(bb)
  if (xwidth == -999) {
    xwidth <- aa$width
  }
  if (xbegin == -999) {
    xbegin <- aa$xl
  }
  if (ybegin == -999) {
    ybegin <- aa$yl
  }
  if (m == 3) {
    zup <- max(TT[, 3])
    zlo <- min(TT[, 3])
    zpos <- zslice
    for (ii in 1:length(zslice)) {
      zpos[ii] <- (zup + zlo)/2 + (zup - zlo)/2 * zslice[ii]
    }
  }
  xstart <- xbegin
  ystart <- ybegin
  ngrid <- trunc(ngrid/2) * 2
  kx <- ngrid
  ky <- ngrid - 1
  l <- kx * ky
  step <- xwidth/(kx - 1)
  r <- 0
  X <- rep(0, times = kx)
  Y2 <- rep(0, times = ky)
  for (p in 1:kx) {
    ystart <- ybegin
    xstart <- xstart + step
    for (q in 1:ky) {
      ystart <- ystart + step
      r <- r + 1
      X[r] <- xstart
      Y2[r] <- ystart
    }
  }
  TPS <- bendingenergy(TT)
  gamma11 <- TPS$gamma11
  gamma21 <- TPS$gamma21
  gamma31 <- TPS$gamma31
  W <- gamma11 %*% YY
  ta <- t(gamma21 %*% YY)
  B <- gamma31 %*% YY
  WtY <- t(W) %*% YY
  trace <- c(0)
  for (i in 1:m) {
    trace <- trace + WtY[i, i]
  }
  if(m==2){benergy <- (16 * pi) * trace}else{benergy <- (8 * pi) * trace}
  l <- kx * ky
  phi <- matrix(0, l, m)
  s <- matrix(0, k, 1)
  
  for (islice in 1:length(zslice)) {
    if (m == 3) {
      refc <- matrix(c(X, Y2, rep(zpos[islice], times = kx * 
                                    ky)), kx * ky, m)
    }
    if (m == 2) {
      refc <- matrix(c(X, Y2), kx * ky, m)
    }
    for (i in 1:l) {
      s <- matrix(0, k, 1)
      for (im in 1:k) {
        s[im, ] <- shapes::sigmacov(refc[i, ] - TT[im, ])
      }
      phi[i, ] <- ta + t(B) %*% refc[i, ] + t(W) %*% s
    }
    if(graphics==T){
      if (m == 3) {
        if (opt == 2) {
          shapes3d(TT, color = 2, axes3 = axes3, rglopen = FALSE)
          shapes3d(YY, color = 4, rglopen = FALSE)
          for (i in 1:k) {
            lines3d(rbind(TT[i, ], YY[i, ]), col = 1)
          }
          for (j in 1:kx) {
            lines3d(refc[((j - 1) * ky + 1):(ky * j), ], 
                    color = 6)
          }
          for (j in 1:ky) {
            lines3d(refc[(0:(kx - 1) * ky) + j, ], color = 6)
          }
        }
        shapes3d(TT, color = collandsTT, axes3 = axes3, rglopen = FALSE)
        shapes3d(YY, color = collandsYY, rglopen = FALSE)
        for (i in 1:k) {
          lines3d(rbind(TT[i, ], YY[i, ]), col = colshift)
        }
        for (j in 1:kx) {
          lines3d(phi[((j - 1) * ky + 1):(ky * j), ], color = colgrid)
        }
        for (j in 1:ky) {
          lines3d(phi[(0:(kx - 1) * ky) + j, ], color = colgrid)
        }
      }
      
    }
    
    
  }
  
  if (m == 2) {
    par(pty = "s")
    if (opt == 2) {
      par(mfrow = c(1, 2))
      order <- linegrid(refc, kx, ky)
      
      if(is.null(xlim)==T){xlim = c(xbegin -xwidth * ext, xbegin + xwidth * (1 + ext))}else{xlim=xlim}
      if(is.null(ylim)==T){ylim = c(ybegin - (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx * (1 + ext))}else{ylim=ylim}
      if(graphics==T){
        plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = xlim, ylim = ylim, xlab = " ", ylab = " ",asp=asp,,axes=axes2d)
        lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], type = "l",col=colgrid,xlim = xlim, ylim = ylim,asp=asp)
        points(TT, cex = cex, pch = pch, col = collandsTT)
      }}
    if(is.null(xlim)==T){xlim = c(xbegin -xwidth * ext, xbegin + xwidth * (1 + ext))}else{xlim=xlim}
    if(is.null(ylim)==T){ylim = c(ybegin - (xwidth * ky)/kx * ext, ybegin + (xwidth * ky)/kx * (1 + ext))}else{ylim=ylim}
    order <- linegrid(phi, kx, ky)
    if(graphics==T){plot(order[1:l, 1], order[1:l, 2], type = "l", xlim = xlim, ylim = ylim, xlab = " ", ylab = " ",col=colgrid,asp=asp,axes=axes2d)
      lines(order[(l + 1):(2 * l), 1], order[(l + 1):(2 * l), 2], type = "l",col=colgrid,xlim = xlim, ylim = ylim,asp=asp)}
    if(graphics==T){points(YY, cex = cex, pch = pch, col = collandsYY)}
    if(graphics==T){points(TT, cex = cex, pch = pch, col = collandsTT)}
    if(graphics==T){
      if(displ==T){
        for (i in 1:(k)) {
          arrows(TT[i, 1], TT[i, 2], YY[i, 1], YY[i, 2], col = coldispl, length = 0.1, lwd=lwdispl,angle = 20)
        }
      }}
    firstcol<-order[1:l,][1:kx,][-kx,]
    firstrow<-order[(l + 1):(2 * l),][((kx*ky-1)-ky):(kx*ky-1),][-c(1:2),]
    lastcol<-order[1:l,][((kx*ky)-ky):(kx*ky),][-1,]
    lastrow<-order[(l + 1):(2 * l),][2:ky,][order((ky:2)),]
    bound = rbind(firstcol,firstrow,lastcol,lastrow)
    them<-nrow(order[1:l,])/ngrid
  }
  if(m==2){grid<-list(ngrid=order[1:l,],ngrid2=order,m=them,n=ngrid,bound=bound,l=l,xlim=xlim,ylim=ylim,TT=TT,YY=YY)}else{
    grid<-list(n=ngrid,l=l,TT=TT,YY=YY)
  }
  out<-list(YY=YY,gamma11=gamma11,gamma21=gamma21,gamma31=gamma31,W=W,ta=ta,B=B,WtY=WtY,trace=trace,benergy=benergy,grid=grid,kx=kx,ky=ky)
  if(graphics==T){if(!is.null(linksTT)){lineplot(TT,linksTT,lwd=lwdtt,col=collinksTT)}}
  if(graphics==T){if(!is.null(linksYY)){lineplot(YY,linksYY,lwd=lwdyy,col=collinksYY)}}
  return(out)
}

conslinks<-function(number,open=T,plus=0){
  k=seq(1:(number-1))
  k=k+plus
  aw=NULL
  for(i in k){
    b=list(c(i,i+1))
    aw<-c(aw,b)
  }
  if(open==T){return(aw)}else{
    aw<-c(aw,list(c(min(unlist(aw)),max(unlist(aw)))))
  }
  aw
}

trapgen<-function(config,a11=0,a12=0,a21=0,a22=0,chi1=0,chi2=0){
  A=matrix(c(a11,a21,a12,a22),ncol=2)
  w<-matrix(c(0,-1,1,0),ncol=2)
  chi<-matrix(c(chi1,0,0,chi2),ncol=2)
  newc<-config+config%*%t(A)+(config%*%t(w))*(config%*%chi)
  #se<-1/24*h*l*(12*a11^2+12*a12^2+12*a21^2+12*a22^2+(h^2+h^2)*(chi1^2 + chi2^2))
  out=list(newc=newc,A=A,chi=chi)
  out
}

tpsjacpsl2d<-function(init,fin,m,areas=NULL,doopa=T){
  
  if(doopa==T){
    theopa<-rotonto(init,fin,scale=F,reflection=F)
    init<-init
    fin<-theopa$yrot
  }
  
  matr<-init
  M<-m
  
  
  
  
  tpsgrid<-tpsgridpaolo(init,fin,graphics=F)
  
  veclist<-NULL
  for(j in 1:nrow(M)){
    vecj<-NULL
    for(i in 1:nrow(matr)){
      
      vec1ji<-2*(M[j,1]-matr[i,1])+2*(M[j,1]-matr[i,1])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vec2ji<-2*(M[j,2]-matr[i,2])+2*(M[j,2]-matr[i,2])*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      vecji<-c(vec1ji,vec2ji)
      vecj<-rbind(vecj,vecji)
    }
    veclist<-c(veclist,list(vecj))
  }
  jac<-NULL
  for(i in 1: length(veclist)){
    jaci<-t(tpsgrid$B)+t(tpsgrid$W)%*%veclist[[i]]
    jac<-c(jac,list(jaci))
  }
  
  deltus<-NULL
  dpl<-NULL
  sens<-NULL
  for(i in 1:length(jac)){
    deltusi<-jac[[i]]-diag(nrow(jac[[i]]))
    dpli<-0.5*(deltusi+t(deltusi)) 
    seni<-0.5*matrix.trace(t(dpli)%*%dpli)
    deltus<-c(deltus,list(deltusi))
    dpl<-c(dpl,list(dpli))
    sens<-c(sens,seni)
  }
  
  if(is.null(areas)==F){sens2<-sum(sens*areas)}else{sens2<-c("no areas")}
  
  myj<-unlist(lapply(jac,det))
  
  myl<-tpsdry2(init,fin,doopa=F)
  veclist2<-NULL
  for(j in 1:nrow(M)){
    vecj<-NULL
    for(i in 1:nrow(matr)){
      
      vec1ji<-2+(4*(M[j,1]-matr[i,1])^2/((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2))+(2*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2))
      vec2ji<-2+(4*(M[j,2]-matr[i,2])^2/((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2))+(2*log((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2))
      vec12ji<-4*(M[j,1]-matr[i,1])*(M[j,2]-matr[i,2])/((M[j,1]-matr[i,1])^2+(M[j,2]-matr[i,2])^2)
      
      arrji<-array(matrix(c(vec1ji,vec12ji,vec12ji,vec2ji),ncol=2),dim=c(2,2,1))
      vecj<-abind::abind(vecj,arrji)
    }
    veclist2<-c(veclist2,list(vecj))
  }
  
  veclist2reshaped<-NULL
  for(z in 1:length(veclist2)){
    veclist2reshapedz<-NULL
    for(k in 1:dim(veclist2[[z]])[[3]]){
      veclist2reshapedzk<-c(veclist2[[z]][,,k])
      veclist2reshapedz<-rbind(veclist2reshapedz,veclist2reshapedzk)
    }
    veclist2reshaped<-c(veclist2reshaped,list(veclist2reshapedz))
  }
  
  jac2<-NULL
  for(i in 1: length(veclist2reshaped)){
    jac2i<-t(myl$W)%*%veclist2reshaped[[i]]
    jac2<-c(jac2,list(as.matrix(jac2i)))
  }
  
  beb<-NULL
  norms<-NULL
  for(i in 1:length(jac2)){
    normi<-norm(jac2[[i]],type="F")
    if(is.null(areas)==F){bebi<-normi^2*areas[[i]]}else{bebi<-c("no areas")}
    beb<-c(beb,bebi)
    norms<-c(norms,normi^2)
  }
  if(is.null(areas)==F){beb<-sum(beb)}else{beb<-c("no areas")}
  cjacps<-NULL
  cjac1psl<-NULL
  ccjac<-NULL
  for(i in 1:length(jac)){
    
    cjaci<-t(jac[[i]])%*%jac[[i]]
    jaceigen<-eigen(cjaci)
    
    
    psl1<-Re(jaceigen$vectors[,1])
    psl2<-Re(jaceigen$vectors[,2])
    psl1<-psl1/(norm(psl1,type="F"))
    psl2<-psl2/(norm(psl2,type="F"))
    psl<-c(psl1,psl2)
    cjac1psl<-rbind(cjac1psl,psl)
    
    cjacpsi<-Re(jaceigen$values)
    cjacps<-rbind(cjacps,cjacpsi)
    
    ccjac<-c(ccjac,list(cjaci))
    
  }
  out=list(jac1=jac,jac2=jac2,sens=sens,sentot=sens2,beb=beb,jac2norms=norms,cjac1psl=cjac1psl,cjacps=cjacps,ccjac=ccjac)
  out
}

tpsdry2<-function(init,fin,doopa=F,meth=c("mor","dm"),g11=T,g31=T){
  library(Morpho)
  if(doopa==T){
    theopa<-rotonto(init,fin,scale=F,reflection=F)
    init<-init
    fin<-theopa$yrot}
  k<-nrow(init)
  m<-ncol(init)
  if(meth=="dm"){
    myl<-CreateL(init)
    if(g11==T){gamma11<-myl$Lsubk}else{gamma11=c("You did not want g11")}
    appo<-myl$Linv[(k+1):(k+1+m),1:k]%*%fin
    W<-gamma11%*%fin
    ct<-appo[1,]
    at<-appo[-1,]
  }else{
    ######## copied from Stefan Schlager code
    trafo<-computeTransform(fin,init,type = "tps")
    if(g11==T){
      myl<-CreateL(init,output="Lsubk")
      gamma11<-myl[[1]]}else{gamma11=c("You did not want g11")}
    
    if(g31==T){
      myl<-CreateL(init,output="Linv")
      gamma31<-myl[[1]][(nrow(myl[[1]])-(m-1)):nrow(myl[[1]]),1:k]}else{gamma31=c("You did not want g31")}
    W<-trafo$coeff[1:k,]
    ct<-trafo$coeff[k+1,]
    at<-trafo$coeff[-c((1:(k+1))),]
  }
  out<-list(gamma11=as.matrix(gamma11),W=as.matrix(W),ct=ct,at=as.matrix(at),gamma31=as.matrix(gamma31))
  out
}

centershapes<-function(array){
  if(is.matrix(array)==T){array<-array(array,dim=c(nrow(array),ncol(array),1))}
  centros<-centroids(array)
  k<-dim(array)[1]
  m<-dim(array)[2]
  n<-dim(array)[3]
  prov<-array(rep(t(centros),each=k),dim=c(k,m,n))
  array2<-array-prov
  return(array2)
}

circle2<-function (radius = 1, origin = c(0, 0),plot=T,by=0.01) {
  t <- seq(-pi, pi, by = by)
  a <- origin[1]
  b <- origin[2]
  r <- radius
  x <- a + r * cos(t)
  y <- b + r * sin(t)
  if(plot==T){points(x, y, type = "l")}
  return(cbind(x,y))
}

centroids<-function(array){
  meanmat<-function(mat){
    mean<-apply(mat,2,mean)
    mean
  }
  if(class(array)=="matrix"){array<-array(array,dim=c(dim(array)[1],dim(array)[2],1))}
  centroidi<-t(apply(array,3,meanmat))
  return(centroidi)
}

defelliadd<-function(circledat,centro=c(0,0,0),mat,col=makeTransparent(1,0.3),mag=1,after=T,border=NULL,lwd=1){
  if(after==T){newdat<-centershapes(circledat%*%t(mat))[,1:2,1]}else{newdat<-centershapes(t(mat%*%t(circledat)))[,1:2,1]}
  points(newdat,asp=1,cex=0)
  newdat2<-newdat*mag+rep.row(centro,nrow(circledat))
  polygon(newdat2[,1],newdat2[,2],col=col,border=border,lwd=lwd)
  newdat2
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

makeTransparent<-function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

pslinfin1<-function(vec,FF,renorm=F){
  newvec<-((FF))%*%vec
  newvec2<-newvec/norm(newvec,type="F")
  if(renorm==F){newvec}else{newvec2}
}

plotmyarrays<-function(x,l=c(1:dim(x)[1]),v=c(1:dim(x)[2]),ind=c(1:dim(x)[3]),group=NULL,links=NULL,xlim=range(x),ylim=range(x),zlim=range(x),alpha=1,col=NULL,txt=T,lwd=1,pch=NULL,cex=NULL,asp=1,cextext=1,xlab="",ylab="",zlab="",xaxt="s",yaxt="s",axes=T){
  require(calibrate)
  require(Morpho)
  require(compositions)
  require(rgl)
  #### x is  a two- or three-dimensional array.
  #### l is a vector indicating the number of variables of 1th array dimension to plot
  #### v is a vector indicating the number of variables of 2th array dimension to plot
  #### ind is a vector indicating which elements of 3th array dimension should be plotted
  #### length(group) must equal length(dim(x)[3]) even if you plot just some individuals (the argument "ind")
  
  
  if(class(x)=="matrix"){x<-array(x,dim=c(dim(x)[1],dim(x)[2],1))}
  
  if(length(v)>2){
    
    if(!is.null(group)){
      
      
      open3d()
      for (i in ind){       
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=c(as.numeric(group)[i]),add=T,bbox=F,xlim=xlim,ylim=ylim,zlim=zlim,alpha=alpha)
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=c(as.numeric(group)[i]),add=T)
          if(txt==T){text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=c(as.numeric(group)[i]),add=T,cex=cextext)}
        }else{NULL}}
      
    }else{
      
      if(is.null(col)){col=c(1:length(ind))}else{col=col}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      for (i in ind){  
        plot3D(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),col=col[i],cex=cex[i],add=T,bbox=F,xlim=xlim,ylim=ylim,zlim=zlim)
        if(txt==T){ text3d(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),text=l,col=col[which(ind==i)],add=T,cex=cextext)}
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i],x[l,v[3],i]),links,lwd=lwd,col=col[which(ind==i)],add=T)}
      }
    }
  }
  
  
  else{
    
    if(!is.null(group)){
      if(is.null(pch)){pch==rep(1:length(ind))}else{pch=pch}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      
      for (i in ind){
        par(new=T)
        plot(x[l,v[1],i],x[l,v[2],i],col=c(as.numeric(group)[i]),xlim=xlim,ylim=ylim,pch=pch[i],asp=asp,cex=cex[i],xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,axes=axes)
        if(txt==T){textxy(x[l,v[1],i],x[l,v[2],i],labs=l,col=c(as.numeric(group)[i]),cex=cextext)}
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i]),links,lwd=lwd,col=c(as.numeric(group)[i]))}else{NULL}
        
      }
      
      
    }else{
      
      if(is.null(col)){col=c(1:length(ind))}else{col=col}
      if(is.null(pch)){pch=rep(1,length(ind))}else{pch=pch}
      if(is.null(cex)){cex=rep(1,length(ind))}else{cex=cex}
      
      
      for (i in ind){
        plot(x[l,v[1],i],x[l,v[2],i],col=col[which(ind==i)],xlim=xlim,ylim=ylim,pch=pch[which(ind==i)],asp=asp,cex=cex[which(ind==i)],xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,axes=axes)
        if(txt==T){textxy(x[l,v[1],i],x[l,v[2],i],labs=l,col=col[which(ind==i)],cex=cextext)}
        if(!is.null(links)){
          lineplot(cbind(x[l,v[1],i],x[l,v[2],i]),links,lwd=lwd,col=col[which(ind==i)])}else{NULL}
        par(new=T)
      }
    }
  }
  par(new=F)
}

mopa<-function(tt,yy,rot=c("mopa","opa"),CSinit = FALSE,volinit=FALSE,center=TRUE){
  warning("the SECOND input matrix will be rotated on the FIRST one")
  if(is.matrix(tt)==T){
    k<-nrow(tt)
    m<-ncol(tt)
    ttar<-array(tt,dim=c(k,m,1))
    yyar<-array(yy,dim=c(k,m,1))}else {
      ttar=tt
      yyar=yy
      k<-dim(tt)[1]
      m<-dim(tt)[2]
    }
  
  if(center==T){ttcs<-centershapes(ttar)[,,1]
  yycs<-centershapes(yyar)[,,1]} 
  
  if(CSinit==T){
    yycs<-scaleshapes(yycs)
    ttcs<-scaleshapes(ttcs)
  }else{yycs<-yycs;ttcs<-ttcs}
  
  if(volinit==T){
    at<-tpsdry2(ttcs,yycs,meth="mor",g11=F)$at
    detat<-det(at)
    if(ncol(yycs)>2){yycs<-yycs/(detat^(1/3))}else{yycs<-yycs/(detat^(1/2))}
  }
  
  if(rot=="opa"){rotprocstep1<-svd(t(ttcs)%*%yycs)}
  stef<-CreateL(ttcs)
  W<-(stef$Lsubk%*%yycs)
  appo<-stef$Linv[(k+1):(k+(m+1)),1:k]%*%yycs
  cT<-appo[1,]
  At<-appo[-1,]
  if(rot=="mopa"){
    rotprocstep1<-svd(At)}
  rotprocstep2<-rotprocstep1$v%*%t(rotprocstep1$u)
  opizzata1<-yycs%*%rotprocstep2
  opizzata2<-array(opizzata1,dim=c(k,m,1))[,,1,drop=F]####centershapes(array(opizzata1,dim=c(k,m,1)))[,,1,drop=F]
  Wafter<-stef$Lsubk%*%opizzata2[,,1]
  appoafter<-stef$Linv[(k+1):(k+(m+1)),1:k]%*%opizzata2[,,1]
  cTafter<-appoafter[1,]
  Atafter<-appoafter[-1,]
  out=list(opizzata=opizzata2,W=W,cT=cT,At=At,Wafter=Wafter,cTafter=cTafter,Atafter=Atafter)
  return(out)
}

tpsjacpsl3d<-function(source,target,m,vols=NULL,doopa=T){
  library(Morpho)
  library(geometry)
  library(pracma)
  mate<-source
  mate2<-target
  mate<-mate
  mate2<-mate+(mate2-mate)*1
  if(doopa==T){
    theopa<-rotonto(mate,mate2,scale=F,reflection=F)
    mate<-mate
    mate2<-theopa$yrot
  }
  der1<-function(a,b){
    res<- -(a-b)/sqrt(sum((a-b)^2))
    res
  }
  vec.new<-function(mat1,mat2=NULL){
    tt<-mat1
    yy<-mat2
    SGMr <- apply(tt,1,function(t)apply(yy,1,der1,t))
    #replace(SGMr,which(SGMr=="NaN",arr.ind=T),0)
    SGMr
  }
  M<-m
  source<-mate
  target<-mate2
  matr<-source
  veclist<-NULL
  resue<-vec.new(M,matr)
  dim(resue)
  for(i in 1:ncol(resue)){
    veclisti<-t(matrix(resue[,i],nrow=3))
    veclist<-c(veclist,list(veclisti))
  }
  myl<-tpsdry2(mate,mate2,doopa=F,meth="mor",g11=F)
  myl2<-myl
  jac<-NULL
  for(i in 1: length(veclist)){
    jaci<-myl2$at+t(myl2$W)%*%veclist[[i]]
    jac<-c(jac,list(t(jaci)))
  }
  jac<-lapply(jac,as.matrix)
  
  deltus<-NULL
  dpl<-NULL
  sens<-NULL
  for(i in 1:length(jac)){
    deltusi<-jac[[i]]-diag(nrow(jac[[i]]))
    dpli<-0.5*(deltusi+t(deltusi)) 
    seni<-0.5*matrix.trace(t(dpli)%*%dpli)
    deltus<-c(deltus,list(deltusi))
    dpl<-c(dpl,list(dpli))
    sens<-c(sens,seni)
  }
  
  if(is.null(vols)==F){sentot<-sum(sens*vols)}else{sentot<-c("no vols")}
  
  der2_1<-function(a,b){
    vec1ji<-(-(a[2] - b[2])^2 - (a[3] - b[3])^2)/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2)
    vec12ji<-((a[1] - b[1])*(a[2] - b[2]))/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2)
    vec13ji<-((a[1] - b[1])*(a[3] - b[3]))/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2)
    vec2ji<-(-(a[1] - b[1])^2 - (a[3] - b[3])^2)/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2) 
    vec23ji<-((a[2] - b[2])*(a[3] - b[3]))/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2)
    vec3ji<-(-(a[1] - b[1])^2 - (a[2] - b[2])^2)/((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)^(3/2)
    arrji<-array(matrix(c(vec1ji,vec12ji,vec13ji,vec12ji,vec2ji,vec23ji,vec13ji,vec23ji,vec3ji),ncol=3),dim=c(3,3,1))
    arrji
  }
  vec.new2<-function(mat1,mat2=NULL){
    tt<-mat1
    yy<-mat2
    SGMr <- apply(tt,1,function(t)apply(yy,1,der2_1,t))
    #replace(SGMr,which(SGMr=="NaN",arr.ind=T),0)
    SGMr
  }
  prova<-vec.new2(M,matr)
  prova2<-t(prova)
  prova3<-split(prova2, 1:nrow(prova2))
  prova4<-lapply(prova3,function(x) t(Reshape(c(x),9,nrow(matr))))
  veclist2reshaped<-prova4
  
  jac2<-NULL
  for(i in 1: length(veclist2reshaped)){
    jac2i<-t(myl$W)%*%veclist2reshaped[[i]]
    jac2<-c(jac2,list(as.matrix(jac2i)))
  }
  
  beb<-NULL
  norms<-NULL
  for(i in 1:length(jac2)){
    normi<-norm(as.matrix(jac2[[i]]),type="F")
    if(is.null(vols)==F){bebi<-normi^2*vols[[i]]}else{bebi<-c("no vols")}
    
    beb<-c(beb,bebi)
    norms<-c(norms,normi^2)
  }
  if(is.null(vols)==F){beb<-sum(beb)}else{beb<-c("no vols")}
  
  
  cjacps<-NULL
  cjac1psl<-NULL
  ccjac<-NULL
  for(i in 1:length(jac)){
    
    cjaci<-t(jac[[i]])%*%jac[[i]]
    jaceigen<-eigen(cjaci)
    
    
    psl1<-Re(jaceigen$vectors[,1])
    psl2<-Re(jaceigen$vectors[,2])
    psl3<-Re(jaceigen$vectors[,3])
    
    psl1<-psl1/(norm(as.matrix(psl1),type="F"))
    psl2<-psl2/(norm(as.matrix(psl2),type="F"))
    psl3<-psl3/(norm(as.matrix(psl3),type="F"))
    
    psl<-c(psl1,psl2,psl3)
    cjac1psl<-rbind(cjac1psl,psl)
    
    cjacpsi<-Re(jaceigen$values)
    cjacps<-rbind(cjacps,cjacpsi)
    
    ccjac<-c(ccjac,list(cjaci))
    
  }
  
  out=list(jac1=jac,jac2=jac2,sens=sens,sentot=sentot,beb=beb,jac2norms=norms,cjac1psl=cjac1psl,cjacps=cjacps,ccjac=ccjac,myl=myl)
  out
}

psltri1<-function(matri1,matri2){
  centro1<-centroids(matri1)
  centro2<-centroids(matri2)
  a1<-matri1[1,]-centro1
  a2<-matri1[2,]-centro1
  a3<-CrossProduct3D(a1,a2)
  a_1<-CrossProduct3D(a2,a3)/(mydot(a1,CrossProduct3D(a2,a3)))
  a_2<-CrossProduct3D(a3,a1)/(mydot(a1,CrossProduct3D(a2,a3)))
  a1f<-matri2[1,]-centro2
  a2f<-matri2[2,]-centro2
  FF1<-t(matrix(kronecker(a1f,a_1)+kronecker(a2f,a_2),nrow=3))
  CC1<-t(FF1)%*%FF1
  
  
  a1<-matri1[2,]-centro1
  a2<-matri1[3,]-centro1
  a3<-CrossProduct3D(a1,a2)
  a_1<-CrossProduct3D(a2,a3)/(mydot(a1,CrossProduct3D(a2,a3)))
  a_2<-CrossProduct3D(a3,a1)/(mydot(a1,CrossProduct3D(a2,a3)))
  a1f<-matri2[2,]-centro2
  a2f<-matri2[3,]-centro2
  FF2<-t(matrix(kronecker(a1f,a_1)+kronecker(a2f,a_2),nrow=3))
  CC2<-t(FF2)%*%FF2
  
  a1<-matri1[3,]-centro1
  a2<-matri1[1,]-centro1
  a3<-CrossProduct3D(a1,a2)
  a_1<-CrossProduct3D(a2,a3)/(mydot(a1,CrossProduct3D(a2,a3)))
  a_2<-CrossProduct3D(a3,a1)/(mydot(a1,CrossProduct3D(a2,a3)))
  a1f<-matri2[3,]-centro2
  a2f<-matri2[1,]-centro2
  FF3<-t(matrix(kronecker(a1f,a_1)+kronecker(a2f,a_2),nrow=3))
  CC3<-t(FF3)%*%FF3
  
  
  
  CC<-(CC1+CC2+CC3)/3
  FF<-(FF1+FF2+FF3)/3
  
  e1<-(matri1[2,]-matri1[1,])/norm(matri1[2,drop=F]-matri1[1,,drop=F],type="F")
  e3<-CrossProduct3D(matri1[2,]-matri1[1,],matri1[3,]-matri1[2,])/norm(t(as.matrix(CrossProduct3D(matri1[2,]-matri1[1,],matri1[3,]-matri1[2,]),type="F")))
  e2<-CrossProduct3D(e3,e1)
  
  CCred11<-mydot(CC%*%(e1),c(e1))
  CCred21<-mydot(CC%*%(e2),c(e1))
  CCred12<-mydot(CC%*%(e1),c(e2))
  CCred22<-mydot(CC%*%(e2),c(e2))
  CCred<-matrix(c(CCred11,CCred21,CCred12,CCred22),ncol=2)
  
  
  FFred11<-mydot(FF%*%(e1),c(e1))
  FFred21<-mydot(FF%*%(e2),c(e1))
  FFred12<-mydot(FF%*%(e1),c(e2))
  FFred22<-mydot(FF%*%(e2),c(e2))
  FFred<-matrix(c(FFred11,FFred21,FFred12,FFred22),ncol=2)
  
  EEred<-0.5*(CCred-diag(2))
  snl<-0.5*matrix.trace(t(EEred)%*%EEred)+24*matrix.trace(EEred)^2### il numero ? arbitrario...il valore varierebbe tra 0 e infinito a seconda della compressibilit? del materiale (qui platonica e costante tra le varie forme)
  
  
  
  
  eigs<-eigen(CC)
  ps1<-eigs$values[1]
  ps2<-eigs$values[2]
  #psl1<-FF%*%eigs$vectors[,1]/(norm(FF%*%eigs$vectors[,1],type="F"))
  #psl2<-FF%*%eigs$vectors[,2]/(norm(FF%*%eigs$vectors[,2],type="F"))
  
  psl1<-eigs$vectors[,1]/(norm(as.matrix(eigs$vectors[,1]),type="F"))
  psl2<-eigs$vectors[,2]/(norm(as.matrix(eigs$vectors[,2]),type="F"))
  
  
  
  if(ps1<1){alpha1=1}else{alpha1=0}
  if(ps2<1){alpha2=1}else{alpha2=0}
  
  
  
  
  ##### procedura angolo
  newpsl2<-CrossProduct3D(a3/norm(t(as.matrix(a3,type="F"))),eigs$vectors[,1])
  
  eigccred<-eigen(CCred)
  c1<-eigccred$vectors[1,1]
  c2<-eigccred$vectors[1,2]
  
  if(sign(c1)>0){
    if(sign(c2)>0){
      angolo<--(acos(c1)*180/pi)
    }else{angolo<-(acos(c1)*180/pi)}
  }else{
    
    if(sign(c2)>0){
      angolo<-(pi-(acos(c1)))*180/pi
    }else{angolo<-(-pi+(acos(c1)))*180/pi}
    
  }
  
  
  out<-list(ps1=ps1,ps2=ps2,psl1=c(psl1),psl2=c(psl2),alpha1=alpha1,alpha2=alpha2,angolo=angolo,CC=CC,CCred=CCred,e3norm=norm(t(as.matrix(e3)),type="F"),lasta3norm=norm(t(as.matrix(a3)),type="F"),FF=FF,FFred=FFred,snl=snl)
  unlist(out)
}

CrossProduct3D <- function(x, y, i=1:3) {#from:https://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)
  
  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1
  
  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

mydot<-function(a,b,G=diag(nrow(a))){
  library(matrixcalc)
  res<-matrix.trace(t(a)%*%G%*%b)
  res
}

plotpsltri1<-function(source,target,ob,usec=T,mag=1,colelli=1){
  if(usec==T){where=centroids(source)}else{where=centroids(target)}
  if(usec==T){defmat<-Re(expm::sqrtm(matrix(ob[12:20],ncol=3)))}else{defmat<-matrix(ob[27:35],ncol=3)}
  mytrinit<-source
  circ3d<-centershapes(cbind(circle2(radius=1,plot=F),0))[,,1]*mag
  circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytrinit,209))),mytrinit[1:2,]),circ3d,reflection=F)
  forma<-centershapes(circ3d2$yrot)[,,1]
  tcirc3d<-t(defmat%*%t(forma))
  lines3d(tcirc3d+rep.row(where,nrow(circ3d)),col=colelli)
  if(ob[1]>1){col1=1;code1=2}else{col1="violet";code1=1}
  if(ob[2]>1){col2=1;code2=2}else{col2="violet";code2=1}
  if(usec==T){mag1<-sqrt(ob[1])}else{mag1=1}
  if(usec==T){mag2<-sqrt(ob[2])}else{mag2=1}
  if(usec==T){psl1<-ob[3:5]*mag1*mag}else{psl1<-pslinfin1(ob[3:5],matrix(ob[27:35],ncol=3))*mag1*mag}
  if(usec==T){psl2<-ob[6:8]*mag2*mag}else{psl2<-pslinfin1(ob[6:8],matrix(ob[27:35],ncol=3))*mag2*mag}
  compositions::arrows3D(where,where+c(psl1),add=T,code=code1,col=col1)
  compositions::arrows3D(where,where-c(psl1),add=T,code=code1,col=col1)
  compositions::arrows3D(where,where+c(psl2),add=T,code=code2,col=col2)
  compositions::arrows3D(where,where-c(psl2),add=T,code=code2,col=col2)
}

list2matrix<-function(mylist){
  final<-NULL
  for(i in 1:length(mylist)){
    temp<-mylist[[i]]
    final<-rbind(final,temp)
  }
  #if(is.null(names(mylist))==F){rownames(final)<-names(mylist)}
  return(final)
}

array2list<-function(array){
  
  
  thelist<-NULL
  
  for(i in 1:dim(array)[3]){
    
    eli<-array[,,i]
    
    thelist<-c(thelist,list(eli))
  }
  if(is.null(dimnames(array)[[3]])==F){names(thelist)<-dimnames(array)[[3]]}
  
  return(thelist)
}

reparray<-function(array,n,type=c("sequ","block")){
  
  
  if(is.matrix(array)==T){array<-array(array,dim=c(nrow(array),ncol(array),1))}
  if(is.null(dimnames(array)[3])){dimnames(array)[3]<-list(c(1:dim(array)[3]))}else{dimnames(array)[3]<-dimnames(array)[3]}
  
  require(abind)
  myarray2<-NULL
  steps<-n
  
  if(type=="sequ"){
    
    
    for(i in 1:dim(array)[3]){
      temp1<-rep(list(array(array[,,i],dim=c(dim(array)[1],dim(array)[2],1))),steps)
      temp2<-list2array(temp1)
      myarray2<-abind::abind(myarray2,temp2)
    }
    return(myarray2)}else{NULL}
  
  
  if(type=="block"){
    temp1<-list2matrix(rep(array2list(array),n))
    temp2<-matrix2arrayasis(temp1,dim(array)[1])
    return(temp2)
  }else{NULL}
}

traslamesh<-function(mesh,c){
  newvb<-t(mesh$vb[-4,])+rep.row(c,nrow(t(mesh$vb[-4,])))
  newmesh<-list(vb=t(cbind(newvb,1)),it=mesh$it)
  class(newmesh)<-"mesh3d"
  newmesh
}

defosph<-function(sphere,mat,after=T){
  if(after==T){newvb<-t(sphere$vb[-4,])%*%t(mat)}else{newvb<-t(mat%*%sphere$vb[-4,])}
  newmesh<-list(vb=t(cbind(newvb,1)),it=sphere$it)
  class(newmesh)<-"mesh3d"
  newmesh
}

creasph<-function(radius=1,centroid=c(0,0,0),subdivision=3){
  temp_sphere<-vcgSphere(subdivision = subdivision)
  temp_sphere$vb[1,]<-temp_sphere$vb[1,]+centroid[1]
  temp_sphere$vb[2,]<-temp_sphere$vb[2,]+centroid[2]
  temp_sphere$vb[3,]<-temp_sphere$vb[3,]+centroid[3]
  final_sphere<-scalemesh(temp_sphere, radius, center = "none")
  return(final_sphere)
}

list2array<-function(mylist){
  require(abind)
  final<-NULL
  for(i in 1:length(mylist)){
    temp<-array(mylist[[i]],dim=c(nrow(mylist[[i]]),ncol(mylist[[i]]),1))
    final<-abind::abind(final,temp)
  }
  
  if(is.null(names(mylist))==F){dimnames(final)[[3]]<-names(mylist)}
  
  
  return(final)
}

psltris<-function(mat1,mat2,triang,doopa=F){
  if(nrow(triang)>3){triang<-t(triang)}
  if(doopa==T){
    theopa<-rotonto(mat1,mat2,scale=F,reflection=F)
    mat1<-mat1
    mat2<-theopa$yrot}
  listA <- mesh2listri(mat1,triang)
  listB <- mesh2listri(mat2,triang)
  part<-mcmapply(function(x,y) psltri1(x,y),x=listA,y=listB,SIMPLIFY = F)
  allres<-do.call(rbind,part)
  colnames(allres)<-c("ps1","ps2","psl1eigx","psl1eigy","psl1eigz","psl2eigx","psl2eigy","psl2eigz","alpha1","alpha2","angolo","CC11","CC12","CC13","CC21","CC22","CC23","CC31","CC32","CC33","CCred11","CCred12","CCred21","CCred22","e3norm","lasta3norm","FF11","FF12","FF13","FF21","FF22","FF23","FF31","FF32","FF33","FFred11","FFred12","FFred21","FFred22","snl")                                                                                      
  ff<-NULL
  cc<-NULL
  for(i in 1:nrow(allres)){
    ffi<-matrix(allres[i,27:35],nrow=3)
    cci<-matrix(allres[i,12:20],nrow=3)
    ff<-c(ff,list(ffi))
    cc<-c(cc,list(cci))
  }
  allres2<-list(init=mat1,fin=mat2,res=allres,triang=triang,ff=ff,cc=cc)
  allres2
}

plotpsltris<-function(psltrisob,col=1,coldist=1,colcomp="violet",wh=c("both","first","second","princ2d","none"),mag=1,lah=0,lwdsl=2,newscene=F,plotelli=T,elliax=T,colelli=2,lwd=2,plotsph=T,col2dd=F,colprinc2d=NULL,zlimcolprinc2d=NULL,col3dd=F,colprinc=NULL,zlimcolprinc=NULL,onlyprinc3d=T,alphasph=0.7,colsph=2,nex=F,tpspsl=T,tpspoints=NULL,drawpoints=NULL,drawpelli=NULL,usec=F){
  
  if(nrow(psltrisob$triang)<2){psltrisob$triang<-t(psltrisob$triang)}
  
  centros2<-NULL
  for(i in 1:ncol(psltrisob$triang)){
    centros2i<-centroids(psltrisob$fin[psltrisob$triang[,i],])
    centros2<-rbind(centros2,centros2i)
  }
  if(length(colelli<2)){colelli<-rep(colelli,nrow(centros2))}else{colelli<-colelli}
  if(length(colsph<2)){colsph<-rep(colsph,nrow(centros2))}else{colsph<-colsph}
  
  if(is.null(tpspoints)==T){tpspoints<-centrotri(psltrisob$init,t(psltrisob$triang))}else{tpspoints<-tpspoints}
  if(is.null(drawpoints)==T){drawpoints<-centrotri(psltrisob$fin,t(psltrisob$triang))}else{drawpoints<-drawpoints}
  if(is.null(drawpelli)==T){drawpelli<-centros2}else{drawpelli<-drawpelli}
  
  centros2<-drawpelli
  
  
  if(usec==F){pslfin<-Re(pslinfin(psltrisob,renorm=F))
  }else{
    pslfin<-cbind(psltrisob$res[,3:5]*sqrt(psltrisob$res[,1]),psltrisob$res[,6:8]*sqrt(psltrisob$res[,2]))}
  
  mag2=rep(mag,nrow(centros2))
  mag1=rep(mag,nrow(centros2))
  sign2d1<-sign(psltrisob$res[,1]-1)
  code2d1<-ifelse(sign2d1>0,2,1)
  colarr2d1=ifelse(sign2d1>0,coldist,colcomp)
  sign2d2<-sign(psltrisob$res[,2]-1)
  code2d2<-ifelse(sign2d2>0,2,1)
  colarr2d2=ifelse(sign2d2>0,coldist,colcomp)
  inimat<-centros2
  finmat2p<-centros2+pslfin[,4:6]*mag1
  finmat2m<-centros2-pslfin[,4:6]*mag1
  finmat1p<-centros2+pslfin[,1:3]*mag1
  finmat1m<-centros2-pslfin[,1:3]*mag1
  inimat1p<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  inimat1m<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  inimat2p<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  inimat2m<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  finmat1p2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  finmat1m2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  finmat2p2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  finmat2m2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
  for(i in 1:nrow(inimat)){
    if(code2d1[i]==2){inimat1p[i,]<-inimat[i,]}else{inimat1p[i,]<-finmat1p[i,]} 
    if(code2d1[i]==2){inimat1m[i,]<-inimat[i,]}else{inimat1m[i,]<-finmat1m[i,]} 
    if(code2d1[i]==2){finmat1p2[i,]<-finmat1p[i,]}else{finmat1p2[i,]<-inimat[i,]} 
    if(code2d1[i]==2){finmat1m2[i,]<-finmat1m[i,]}else{finmat1m2[i,]<-inimat[i,]} 
    if(code2d2[i]==2){inimat2p[i,]<-inimat[i,]}else{inimat2p[i,]<-finmat2p[i,]} 
    if(code2d2[i]==2){inimat2m[i,]<-inimat[i,]}else{inimat2m[i,]<-finmat2m[i,]} 
    if(code2d2[i]==2){finmat2p2[i,]<-finmat2p[i,]}else{finmat2p2[i,]<-inimat[i,]} 
    if(code2d2[i]==2){finmat2m2[i,]<-finmat2m[i,]}else{finmat2m2[i,]<-inimat[i,]} 
  }
  
  princ2d<-apply(abs(1-Re(psltrisob$res[,1:2])),1,which.max)
  codeprinc2d<-ifelse(princ2d<2,code2d1,code2d2)
  eigprinc2d<-NULL
  for(i in 1:length(princ2d)){
    eigprinc2di<-psltrisob$res[i,1:2][princ2d[i]]
    eigprinc2d<-c(eigprinc2d,eigprinc2di)
  }
  
  der13d<-tpsjacpsl3d(psltrisob$init,psltrisob$fin,tpspoints,doopa=F)
  if(plotsph==T){
    if(usec==F){defmats<-der13d$jac1}else{defmats<-lapply(der13d$ccjac,expm::sqrtm)}
    finpoints<-drawpoints
    for(i in 1:length(der13d$jac1)){
      sphi<-creasph(radius=mag,subdivision=2)
      shade3d(traslamesh(scalemesh(defosph(sphi,defmats[[i]],after=F),1,center="none"),finpoints[i,]),col=colsph[i],alpha=alphasph)
    }
  }
  dist1<-abs(1-der13d$cjacps[,1])
  sign1<-sign(der13d$cjacps[,1]-1)
  code1<-ifelse(sign1>0,2,1)
  colarr1=ifelse(sign1>0,coldist,colcomp)
  
  dist2<-abs(1-der13d$cjacps[,2])
  sign2<-sign(der13d$cjacps[,2]-1)
  code2<-ifelse(sign2>0,2,1)
  colarr2=ifelse(sign2>0,coldist,colcomp)
  
  dist3<-abs(1-der13d$cjacps[,3])
  sign3<-sign(der13d$cjacps[,3]-1)
  code3<-ifelse(sign3>0,2,1)
  colarr3=ifelse(sign3>0,coldist,colcomp)
  
  indprinc<-apply(cbind(dist1,dist2,dist3),1,which.max)
  
  eigprinc<-NULL
  for(i in 1:nrow(der13d$cjacps)){
    eigprinci<-der13d$cjacps[i,indprinc[i]]
    eigprinc<-c(eigprinc,eigprinci)
  }
  
  pslprinc<-NULL
  for(i in 1:length(indprinc)){
    if(indprinc[i]==1){pslprinci<-der13d$cjac1psl[i,1:3]}
    if(indprinc[i]==2){pslprinci<-der13d$cjac1psl[i,4:6]}
    if(indprinc[i]==3){pslprinci<-der13d$cjac1psl[i,7:9]}
    pslprinc<-rbind(pslprinc,pslprinci)
  }
  
  
  jacprinc<-NULL
  for(i in 1:nrow(der13d$cjacps)){
    jacprinci<-der13d$jac1[[i]]
    jacprinc<-c(jacprinc,list(jacprinci))
  }
  
  codeprinc<-NULL
  codemat<-cbind(code1,code2,code3)
  for(i in 1:nrow(codemat)){
    codeprinci<-codemat[i,indprinc[i]]
    codeprinc<-c(codeprinc,codeprinci)
  }
  
  coldprinc<-NULL
  codedd<-cbind(colarr1,colarr2,colarr3)
  for(i in 1:nrow(codedd)){
    coldprinci<-codedd[i,indprinc[i]]
    coldprinc<-c(coldprinc,coldprinci)
  }
  if(is.null(zlimcolprinc)==T){zlimcolprinc<-range(sqrt(eigprinc),sqrt(eigprinc2d))}else{zlimcolprinc=zlimcolprinc}
  if(is.null(colprinc)==T&col3dd==F){colprinc<-numcol(sqrt(eigprinc),zlim=zlimcolprinc,colors=c("black","cyan","violet","blue","red"))}
  if(is.null(colprinc)==T&col3dd==T){colprinc<-coldprinc}
  if(is.null(colprinc)==F&col3dd==F){colprinc<-colprinc}
  if(is.null(colprinc)==F&col3dd==T){colprinc<-coldprinc}
  
  coldprinc2d<-ifelse(princ2d<2,colarr2d1,colarr2d2)
  
  if(is.null(zlimcolprinc2d)==T){zlimcolprinc2d<-range(sqrt(eigprinc),sqrt(eigprinc2d))}else{zlimcolprinc2d=zlimcolprinc2d}
  if(is.null(colprinc2d)==T&col2dd==F){colprinc2d<-numcol(sqrt(eigprinc2d),zlim=zlimcolprinc2d,colors=c("black","cyan","violet","blue","red"))}
  if(is.null(colprinc2d)==T&col2dd==T){colprinc2d<-coldprinc2d}
  if(is.null(colprinc2d)==F&col2dd==F){colprinc2d<-colprinc2d}
  if(is.null(colprinc2d)==F&col2dd==T){colprinc2d<-coldprinc2d}
  
  par(mfrow=c(2,2))
  h<-heat3d(psltrisob$init,psltrisob$fin,t(psltrisob$triang),graphics=F)
  theseq2d<-seq(zlimcolprinc2d[1],zlimcolprinc2d[2],length.out=1000)
  theseq2d<-theseq2d[order(theseq2d)]
  
  plot(cbind(seq(range(h$obs)[1],range(h$obs)[2],length.out=1000),0),col=numcol(seq(range(h$obs)[1],range(h$obs)[2],length.out=1000)),color=c("blue","cyan","yellow","red"),pch=15,cex=4,xaxt="n",yaxt="n",xlab="",ylab="",axes=F)
  axis(1)
  
  theseq3d<-seq(zlimcolprinc[1],zlimcolprinc[2],length.out=1000)
  theseq3d<-theseq3d[order(theseq3d)]
  plot(cbind(theseq3d,0),col=numcol(theseq3d, zlim=zlimcolprinc,colors=c("black","cyan","violet","blue","red")),pch=15,cex=4,xaxt="n",yaxt="n",xlab="",ylab="",axes=F)
  axis(1)
  plot(density(h$obs),main="",xlab="",ylab="")
  densgm(cbind(sqrt(eigprinc2d),sqrt(eigprinc)))
  
  
  if(elliax==T){
    if(wh=="both"){
      compositions::arrows3D(inimat1p,finmat1p2,col=rep(colarr2d1,each=2),code=code2d1,length=lah,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat1m,finmat1m2,col=rep(colarr2d1,each=2),code=code2d1,length=lah,add=T,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat2p,finmat2p2,col=rep(colarr2d2,each=2),code=code2d2,length=lah,add=T,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat2m,finmat2m2,col=rep(colarr2d2,each=2),code=code2d2,length=lah,add=T,add=T,lwd=lwdsl)
    }
    
    if(wh=="first"){
      compositions::arrows3D(inimat1p,finmat1p2,col=rep(colarr2d1,each=2),code=code2d1,length=lah,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat1m,finmat1m2,col=rep(colarr2d1,each=2),code=code2d1,length=lah,add=T,add=T,lwd=lwdsl)
    }
    if(wh=="second"){
      compositions::arrows3D(inimat2p,finmat2p2,col=rep(colarr2d2,each=2),code=code2d2,length=lah,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat2m,finmat2m2,col=rep(colarr2d2,each=2),code=code2d2,length=lah,add=T,add=T,lwd=lwdsl)
    }
    
    if(wh=="princ2d"){
      princ2dcond1<-ifelse(as.numeric(princ2d)==1,TRUE,FALSE)
      princ2dcond2<-ifelse(princ2d==2,TRUE,FALSE)
      
      if(any(princ2dcond1)){
        compositions::arrows3D(inimat1p[princ2dcond1,],finmat1p2[princ2dcond1,],col=rep(colprinc2d[princ2dcond1],each=2),code=codeprinc2d[princ2dcond1],length=lah,add=T,lwd=lwdsl)
        compositions::arrows3D(inimat1m[princ2dcond1,],finmat1m2[princ2dcond1,],col=rep(colprinc2d[princ2dcond1],each=2),code=codeprinc2d[princ2dcond1],length=lah,add=T,lwd=lwdsl)
      }
      if(any(princ2dcond2)){
        compositions::arrows3D(inimat2p[princ2dcond2,],finmat2p2[princ2dcond2,],col=rep(colprinc2d[princ2dcond2],each=2),code=codeprinc2d[princ2dcond2],length=lah,add=T,lwd=lwdsl)
        compositions::arrows3D(inimat2m[princ2dcond2,],finmat2m2[princ2dcond2,],col=rep(colprinc2d[princ2dcond2],each=2),code=codeprinc2d[princ2dcond2],length=lah,add=T,lwd=lwdsl)
      }
      
    }
    
  }
  
  if(plotelli==T){
    for(i in 1:nrow(centros2)){
      if(usec==F){defmat<-(Re(matrix(psltrisob$res[i,27:35],ncol=3)))}else{defmat<-expm::sqrtm(Re(matrix(psltrisob$res[i,12:20],ncol=3)))}
      
      
      if(usec==T){
        mytrinit<-psltrisob$init[psltrisob$triang[,i],]
        mytri<-psltrisob$fin[psltrisob$triang[,i],]
      }else{
        mytrinit<-psltrisob$fin[psltrisob$triang[,i],]
        mytri<-psltrisob$init[psltrisob$triang[,i],]
        
      }
      
      
      circ3d<-centershapes(cbind(circle2(radius=1,plot=F),0))[,,1]*mag
      
      if(usec==T){circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytrinit,209))),mytrinit[1:2,]),circ3d,reflection=T)$yrot}else{
        circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytri,209))),mytri[1:2,]),circ3d,reflection=T)$yrot
        
      }
      
      
      forma<-centershapes(circ3d2)[,,1]
      
      tcirc3d<-t(defmat%*%t(forma))
      
      lines3d(tcirc3d+rep.row(drawpelli[i,],nrow(circ3d)),col=colelli[i],lwd=lwd)
      
      
      
      # 
      # circ3d<-centershapes(cbind(circle2(radius=1,plot=F,by=0.05),0))[,,1]*mag
      # defcirc<-t(defmat%*%t(circ3d))
      # 
      # circ3d2<-rotonto(rbind(list2matrix(array2list(reparray(mytrinit,40))),mytrinit[1:2,]),defcirc,reflection=T)$yrot
      # 
      # 
      # tcirc3d<-centershapes(circ3d2)[,,1]
      # 
      # lines3d(tcirc3d+rep.row(drawpelli[i,],nrow(circ3d)),col=colelli[i],lwd=lwd)
      # 
      # 
      # 
      
      
      
    }
  }
  
  if(tpspsl==T){
    pslfin21<-NULL
    pslfin22<-NULL
    pslfin23<-NULL
    pslfinprinc<-NULL
    for(i in 1:nrow(der13d$cjac1psl)){
      
      if(usec==F){pslfin21i<-c(pslinfin1(der13d$cjac1psl[i,1:3],der13d$jac1[[i]]))}else{pslfin21i<-c(der13d$cjac1psl[i,1:3])*sqrt(der13d$cjacps[i,])[1]}
      pslfin21<-rbind(pslfin21,pslfin21i)
      
      if(usec==F){pslfin22i<-c(pslinfin1(der13d$cjac1psl[i,4:6],der13d$jac1[[i]]))}else{pslfin22i<-c(der13d$cjac1psl[i,4:6])*sqrt(der13d$cjacps[i,])[2]}
      pslfin22<-rbind(pslfin22,pslfin22i)
      
      if(usec==F){pslfin23i<-c(pslinfin1(der13d$cjac1psl[i,7:9],der13d$jac1[[i]]))}else{pslfin23i<-c(der13d$cjac1psl[i,7:9])*sqrt(der13d$cjacps[i,])[3]}
      pslfin23<-rbind(pslfin23,pslfin23i)
      
      if(usec==F){pslfinprinci<-c(pslinfin1(pslprinc[i,],jacprinc[[i]]))}else{pslfinprinci<-c(pslprinc[i,])*sqrt(eigprinc[i])}
      pslfinprinc<-rbind(pslfinprinc,pslfinprinci)
      
      
    }
    
    mag2=rep(mag,nrow(drawpoints))
    mag1=rep(mag,nrow(drawpoints))
    mag3=rep(mag,nrow(drawpoints))
    
    
    inimat<-drawpoints
    finmat1p<-drawpoints+pslfin21*mag1
    finmat1m<-drawpoints-pslfin21*mag1
    
    finmat2p<-drawpoints+pslfin22*mag2
    finmat2m<-drawpoints-pslfin22*mag2
    
    finmat3p<-drawpoints+pslfin23*mag3
    finmat3m<-drawpoints-pslfin23*mag3
    
    finmatprincp<-drawpoints+pslfinprinc*mag3
    finmatprincm<-drawpoints-pslfinprinc*mag3
    
    inimat1p<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    inimat1m<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    inimat2p<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    inimat2m<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    
    inimat3p<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    inimat3m<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    inimatprincp<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    inimatprincm<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    
    finmat1p2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    finmat1m2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    finmat2p2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    finmat2m2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    finmat3p2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    finmat3m2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    finmatprincp2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    finmatprincm2<-matrix(0,nrow=nrow(inimat),ncol=ncol(inimat))
    
    
    
    for(i in 1:nrow(inimat)){
      
      if(code1[i]==2){inimat1p[i,]<-inimat[i,]}else{inimat1p[i,]<-finmat1p[i,]} 
      if(code1[i]==2){inimat1m[i,]<-inimat[i,]}else{inimat1m[i,]<-finmat1m[i,]} 
      
      if(code1[i]==2){finmat1p2[i,]<-finmat1p[i,]}else{finmat1p2[i,]<-inimat[i,]} 
      if(code1[i]==2){finmat1m2[i,]<-finmat1m[i,]}else{finmat1m2[i,]<-inimat[i,]} 
      
      if(code2[i]==2){inimat2p[i,]<-inimat[i,]}else{inimat2p[i,]<-finmat2p[i,]} 
      if(code2[i]==2){inimat2m[i,]<-inimat[i,]}else{inimat2m[i,]<-finmat2m[i,]} 
      
      if(code2[i]==2){finmat2p2[i,]<-finmat2p[i,]}else{finmat2p2[i,]<-inimat[i,]} 
      if(code2[i]==2){finmat2m2[i,]<-finmat2m[i,]}else{finmat2m2[i,]<-inimat[i,]} 
      
      
      if(code3[i]==2){inimat3p[i,]<-inimat[i,]}else{inimat3p[i,]<-finmat3p[i,]} 
      if(code3[i]==2){inimat3m[i,]<-inimat[i,]}else{inimat3m[i,]<-finmat3m[i,]} 
      
      if(code3[i]==2){finmat3p2[i,]<-finmat3p[i,]}else{finmat3p2[i,]<-inimat[i,]} 
      if(code3[i]==2){finmat3m2[i,]<-finmat3m[i,]}else{finmat3m2[i,]<-inimat[i,]} 
      
      
      
      if(codeprinc[i]==2){inimatprincp[i,]<-inimat[i,]}else{inimatprincp[i,]<-finmatprincp[i,]} 
      if(codeprinc[i]==2){inimatprincm[i,]<-inimat[i,]}else{inimatprincm[i,]<-finmatprincm[i,]} 
      
      if(codeprinc[i]==2){finmatprincp2[i,]<-finmatprincp[i,]}else{finmatprincp2[i,]<-inimat[i,]} 
      if(codeprinc[i]==2){finmatprincm2[i,]<-finmatprincm[i,]}else{finmatprincm2[i,]<-inimat[i,]} 
      
      
    }
    
    if(onlyprinc3d==T){
      
      compositions::arrows3D(inimatprincp,finmatprincp2,col=rep(colprinc,each=2),code=2,length=lah,add=T,lwd=lwdsl)
      compositions::arrows3D(inimatprincm,finmatprincm2,col=rep(colprinc,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      
      
    }else{
      
      compositions::arrows3D(inimat1p,finmat1p2,col=rep(colarr1,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat1m,finmat1m2,col=rep(colarr1,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      
      
      compositions::arrows3D(inimat2p,finmat2p2,col=rep(colarr2,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat2m,finmat2m2,col=rep(colarr2,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      
      compositions::arrows3D(inimat3p,finmat3p2,col=rep(colarr3,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
      compositions::arrows3D(inimat3m,finmat3m2,col=rep(colarr3,each=2),code=2,length=lah,add=T,add=T,lwd=lwdsl)
    }
  }
  if(newscene==F){
    if(nex==T){next3d()}
  }
  out<-list(princ2d=princ2d,eigprinc2d=eigprinc2d,princ3d=indprinc,eigprinc3d=eigprinc)
  out
}

mesh2listri<- function(mat, tri) {
  if (ncol(tri) > 3) tri <- t(tri)
  n <- nrow(tri)
  l <- vector(mode = "list", length = n)
  for (i in 1:n) {
    resi <- mat[tri[i, ], ]
    l[[i]] <- resi
  }
  l
}

crux<-function(x,axlab=NULL,origin=NULL,textsize=0.7,axis.len=0.1,bilat=T){
  library(compositions)
  library(rgl)
  xs2<-x
  xs<-x
  if(!is.null(axlab)){
    vlabs=c(axlab[1],axlab[2],axlab[3])
  }else{vlabs = c("x", "y", "z")}
  
  if(!is.null(origin)){origin=origin}else{origin<-c(mean(xs2[,1]),mean(xs2[,2]),mean(xs2[,3]))}
  
  
  poslengths<-c(abs(max(xs[,1])-origin[1]),abs((max(xs[,2])-origin[2])),abs((max(xs[,3])-origin[3])))
  neglengths<--c(abs(min(xs[,1])-origin[1]),abs((min(xs[,2])-origin[2])),abs((min(xs[,3])-origin[3])))
  
  
  
  axis3D(axis.origin=origin,axis.scale=poslengths,axis.col="black",vlabs=vlabs,axis.len=axis.len)
  
  if(bilat==T){axis3D(axis.origin=origin,axis.scale=neglengths,axis.col="black",axis.len=axis.len,vlabs=vlabs)}
  
  
  
  if(sum(origin-apply(x,2,mean))<1e-8){
    mytextx<-as.character(round(quantile(seq(from=min(xs[,1]),to=max(xs[,1]),length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
    mytexty<-as.character(round(quantile(seq(from=min(xs[,2]),to=max(xs[,2]),length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
    mytextz<-as.character(round(quantile(seq(from=min(xs[,3]),to=max(xs[,3]),length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
    text3d(as.numeric(mytextx),rep(origin[2],5),rep(origin[3],5),texts=mytextx,cex=textsize)
    text3d(rep(origin[1],5),as.numeric(mytexty),rep(origin[3],5),texts=mytexty,cex=textsize)
    text3d(rep(origin[1],5),rep(origin[2],5),as.numeric(mytextz),texts=mytextz,cex=textsize)}else{
      
      
      mytextx<-as.character(round(quantile(seq(from=poslengths[1],to=neglengths[1],length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
      mytexty<-as.character(round(quantile(seq(from=poslengths[2],to=neglengths[2],length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
      mytextz<-as.character(round(quantile(seq(from=poslengths[3],to=neglengths[3],length.out=5),probs = seq(0, 1,length.out=100)),2)[c(10,30,50,70,90)])
      text3d(as.numeric(mytextx),rep(origin[2],5),rep(origin[3],5),texts=mytextx,cex=textsize)
      text3d(rep(origin[1],5),as.numeric(mytexty),rep(origin[3],5),texts=mytexty,cex=textsize)
      text3d(rep(origin[1],5),rep(origin[2],5),as.numeric(mytextz),texts=mytextz,cex=textsize)}
  
  
  
}

centrotri<-function(mat,tria){
  centros<-NULL
  for(i in 1:nrow(tria)){
    centrosi<-centroids(mat[tria[i,],])
    centros<-rbind(centros,centrosi)
  }
  centros
}

numcol<-function(value,colors=c("blue","cyan","yellow","red"), tol=1e-10,num = 1000,zlim=NULL,NAcol = "white",below=1,beyond="violet") {
  
  if(diff(range(na.omit(value)))<tol){num=1}
  
  if(is.null(zlim)==T){zlim<-c(min(na.omit(value)),max(na.omit(value)))}else{zlim=zlim}
  from=zlim[1]
  to=zlim[2]
  cols = colorRampPalette(colors)(num)
  cols = cols[findInterval(value, vec = seq(from =from,to = to, length.out = num),all.inside =T)]
  cols[which(is.na(cols))]<-NAcol
  cols[value>to] <-beyond
  cols[value<from] <-below
  cols
}

heat3d<-function(source,target,triang,iter=0,dointerp=F,linkss=NULL,linkst=NULL,plotlands=F,legend=T,cols=1,colt=2,plotsource=T,plottarget=T,collinkss=1,lwds=2,collinkst=2,lwdt=2,cexs=0.5,cext=0.5,colors=c("blue","cyan","yellow","red"),alpha=1,ngrid=0,mag=1,graphics=T,to=NULL,from=NULL,scaleramp=F,lines=T){
  mate<-source
  mate2<-target
  mate2<-mate+(mate2-mate)*mag
  library(fields)
  mmate<-plotsurf(mate,t(triang),plot=F)
  mmate2<-plotsurf(mate2,t(triang),plot=F)
  
  
  mesh2listri<-function(mat,tri){
    if(ncol(tri)>3){tri<-t(tri)}
    res <- apply(tri, 1, function(x) mat[x,])
    res <- lapply(as.data.frame(res), function(x) matrix(x, nrow = 3, ncol = ncol(tri)))
    res
  }
  
  
  
  ar01<-vcgArea(mmate,perface=T)$pertriangle
  ar02<-vcgArea(mmate2,perface=T)$pertriangle
  
  
  obs0<-log(ar02/ar01)
  
  M0<-vcgBary(mmate)
  M02<-vcgBary(mmate2)
  
  tes<-tessell3d(triang,mate,iter)
  class(tes) <- "mesh3d"
  matr<-mate
  M<-vcgBary(tes)
  tes2<-tessell3d(triang,mate2,iter)
  class(tes2) <- "mesh3d"
  M2<-vcgBary(tes2)
  
  triangtes<-t(tes$it)
  if(iter>0){
    ar1<-vcgArea(tes,perface=T)$pertriangle
    ar2<-vcgArea(tes2,perface=T)$pertriangle
  }
  
  
  if(iter>0){obsref<-log(ar2/ar1)}else{obsref<-obs0}
  if(dointerp==T){fit<- fastTps(M02,obs0,theta=0.009)}else{fit<-NULL} 
  if(dointerp==T){obs2<-predict(fit,rbind(mate2,M2))}else{obs2<-NULL}
  if(dointerp==T){obs3<-c(obs2)}else{obs3<-NULL}
  if(dointerp==T){obs3[is.na(obs3)]<-mean(obs3,na.rm=T)}
  
  
  
  if(dointerp==T){obm<-meshDist(plotsurf(rbind(mate2,M2),mmate$it,plot=F),distv=obs3,add=T,rampcolors =colors,to=to,from=from,scaleramp=scaleramp,plot=F,alpha=alpha,meshColor="legacy")}else{obm<-NULL}
  
  if(graphics==T){
    if(plotlands==T){deformGrid3d(mate,mate2,ngrid=ngrid,lines=lines,col1=cols,col2=colt)}
    if(dointerp==T){shade3d(obm$colMesh,alpha=alpha)}
  }
  if(dointerp==F){
    if(is.null(from)&is.null(to)==T){
      colo<-numcol(obs0,colors=colors)
      if(graphics==T){shade3d(mmate2,col=colo,meshColor="faces",alpha=alpha)}
    }else{
      #colo1<-numcol(c(seq(from,to,length.out=100),obs0),colors=colors)
      #colo<-colo1[-c(1:100)]
      
      colo<-numcol(obs0,colors=colors,zlim=c(from,to))
      
      
      
      if(graphics==T){shade3d(mmate2,col=colo,meshColor="faces",alpha=alpha)}
    }
    
    obm<-list(colMesh=list(vb=mmate2$vb,it=mmate2$it,material=list(color=colo)))
    class(obm$colMesh)<-"mesh3d"
  }
  
  if(is.null(from)&is.null(to)==T){
    colo1<-numcol(seq(min(obs0),max(obs0),length.out=100),colors=colors)
    if(graphics==T){
      plot(cbind(seq(min(obs0),max(obs0),length.out=100),0),col=colo1,cex=10,pch=15,xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
      axis(1)
    }
  }else{
    
    colo1<-numcol(c(seq(from,to,length.out=100)),colors=colors,zlim=c(from,to))
    if(graphics==T){
      plot(cbind(c(seq(from,to,length.out=100)),0),col=colo1,cex=10,pch=15,xaxt="n",yaxt="n",axes=F,xlab="",ylab="")
      axis(1)
    }
  }
  
  if(graphics==T){
    if(!is.null(linkst)){lineplot(mate2,linkst,col=collinkst,lwd=lwds)}
    if(plotsource==T){
      shade3d(plotsurf(source,t(triang),plot=F),alpha=0.5)
      if(!is.null(linkss)){lineplot(mate,linkss,col=collinkss,lwd=lwdt)}}
  }
  if(iter>0){ar1<-ar1}else{ar1<-c("no iter")}
  
  if(iter>0){ar2<-ar2}else{ar2<-c("no iter")}
  out<-list(mate=mate,mate2=mate2,mmate=mmate,mmate2=mmate2,ar01=ar01,ar02=ar02,ar1=ar1,ar2=ar2,M0=M0,M02=M02,M=M,M2=M2,obsref=obsref,obs2=obs2,obs=obs0,obs3=obs3,fit=fit,cols=makeTransparent(colorRampPalette(colors)(n = length(obs3)),alpha=alpha),tes=tes,tes2=tes2,obm=obm)
  out
}

pslinfin<-function(ob,renorm=F){
  res<-NULL
  for(i in 1:nrow(ob$res)){
    vec1i<-ob$res[i,3:5]
    vec2i<-ob$res[i,6:8]
    FFi<-matrix(ob$res[i,27:35],ncol=3)
    res2i<-pslinfin1(vec2i,FFi,renorm=renorm)
    res1i<-pslinfin1(vec1i,FFi,renorm=renorm)
    
    res<-rbind(res,c(res1i,res2i))
  }
  res
}

plotsurf<-function(lmatrix,triang,col=1,alpha=0.5,plot=T){
  thesurf_1=t(cbind(lmatrix,1))
  triangolazioni=triang
  thesurf=list(vb=thesurf_1,it=triangolazioni)
  class(thesurf) <- "mesh3d"
  if(plot==T){shade3d(thesurf,col=col,alpha=alpha)}
  return(thesurf)
}

tessell3d<-function(tri,ver,iter){
  tri_t=tri
  ver_t=ver
  if(iter>0){
    for(j in 1:iter){
      nrows=nrow(tri_t)
      numb=range(tri_t)[2]
      for(i in 1:nrows){
        tri_i=ver_t[tri_t[i,],]
        cen_i=colMeans(tri_i) 
        tri_1=c(tri_t[i,c(1,2)],numb+i)
        tri_2=c(tri_t[i,1],numb+i,tri_t[i,3])
        tri_3=c(tri_t[i,c(2,3)],numb+i)
        tri_t=rbind(tri_t,tri_1,tri_2,tri_3)
        ver_t=rbind(ver_t,cen_i)
      }}
    vertici_out=cbind(ver_t,1)
    rownames(vertici_out)=NULL
    triangoli_out=tri_t[(dim(tri)[1]+1):dim(tri_t)[1],]
    rownames(triangoli_out)=NULL
  }else{
    vertici_out=cbind(ver,1)
    triangoli_out=tri}
  mesh.out=list(vb=t(vertici_out),it=t(triangoli_out))
}

densgm<-function(mat,cols=1:ncol(mat),ty=1:ncol(mat),xlim=NULL,ylim=NULL,lwd=rep(1,ncol(mat))){
  
  densis<-NULL
  for(i in 1:ncol(mat)){
    densisi<-density(mat[,i],na.rm=T)
    densis<-c(densis,list(densisi))
  }
  xx<-subListExtract(densis,"x")
  yy<-subListExtract(densis,"y")
  if(is.null(xlim)==T){xlim<-range(xx)}else{xlim<-xlim}
  if(is.null(ylim)==T){ylim<-range(yy)}else{ylim<-ylim}
  for(i in 1:ncol(mat)){
    plot(densis[[i]],main="",xlab="",ylab="",xlim=xlim,ylim=ylim,col=cols[i],lty=ty[i],xaxt="n",yaxt="n",lwd=lwd[i])
    par(new=T)
  }
  axis(1)
  axis(2)
}

psl2d<-function(init,fin,links=NULL,doopa=T,tri=NULL,ar=100,nadd=5000,tit1=NULL,tit2=NULL,interpfin=F,interpinit=F,cold=1,colp="violet",elliax=T,plotelliax1=T,plotelliax2=T,onlyprinc=F,discretepsl=F,collinks1=1,collinks2=2,lwd=2,mag=1,send=F,snl=F,log=F,bebs=F,plotri1=F,plotri2=F,plotgrid1=T,coltri1=1,plotgrid2=T,coltri2=2,elli=T,elliborder=NULL,vals=NULL,colelli=NULL,zlim=NULL,colors=c("blue","cyan","yellow","red"),colorsleig=c("blue","cyan","yellow","red"),alpha=0.4,plotboth=T,la=0,lwda=1,st=F,plotdens=F,plotlegend=T,slcoleig=F){
  if(doopa==T){
    theopa<-rotonto(init,fin,scale=F,reflection=F)
    init<-init
    fin<-theopa$yrot}
  tris<-triang2dto(init,fin,links=links,ar=ar,Y=F)
  
  if(is.null(tri)==T){thetria<-tris$tria}else{thetria<-tri}
  
  psltrisob<-psltris(cbind(tris$newmat1,0),cbind(tris$newmat2,0),thetria,doopa=doopa)
  pslenes<-enespsl(psltrisob)
  
  if(is.null(tri)==T){
    centrosinit<-tris$centrosnewmat1
    centrosfin<-tris$centrosnewmat2
  }else{
    centrosinit<-centrotri(init,thetria)
    centrosfin<-centrotri(fin,thetria)
    
  }
  newinit<-tris$newmat1
  newfin<-tris$newmat2
  
  
  detsf<-NULL
  areasinit<-NULL
  areasfin<-NULL
  for(i in 1:nrow(psltrisob$res)){
    
    areasiniti<-ar3dtri(cbind(newinit[thetria[i,],],0))
    
    areasinit<-c(areasinit,areasiniti)
    
    areasfini<-ar3dtri(cbind(newfin[thetria[i,],],0))
    
    areasfin<-c(areasfin,areasfini)
    
    detsfi<-det(matrix(psltrisob$res[i,27:35],ncol=3))
    detsf<-c(detsf,detsfi)
  }
  
  areasratio<-areasfin/areasinit
  
  objac2d<-tpsjacpsl2d(init,fin,centrosinit,areas=areasinit,doopa=F)######### where I evaluate jacobian infinitesimally
  beb<-objac2d$beb
  
  if(is.null(vals)==T){vals<-areasratio}else{vals<-vals}
  totsn<-objac2d$sentot
  if(send==T){
    if(discretepsl==F){vals<-objac2d$sens}else{vals<-pslenes[,3]}}
  if(bebs==T){vals<-objac2d$jac2norms}
  if(snl==T){vals<-psltrisob$res[,40]}
  if(log==T){vals<-log2(vals)}
  if(is.null(colelli)==T){colelli<-makeTransparent(numcol(vals,colors=colors,zlim=zlim),alpha=alpha)}else{colelli<-colelli}
  
  
  posimpnoh<-areasip(init,links=links,extpol=T,graph=F)
  
  maufin<-fin[unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)]),]
  mauinit<-init[unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)]),]
  matepro<-fin
  matepro[1:nrow(fin)%in%unique(posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1)])==F,]<-NA
  faclinks<-as.factor(is.na(matepro[,1]))
  newindlinks<-posfac(faclinks)
  maulinks<-posimpnoh$triangob$E[posimpnoh$triangob$EB%in%c(1),]
  #plotmyarrays(mate,links=array2list(matrix2arrayasis(maulinks,1)))
  newmaulinks<-maulinks
  newmaulinks[,1]<-c(1:nrow(maulinks))
  newmaulinks[,2]<-match(maulinks[,2],maulinks[,1])
  #plotmyarrays(mau,links=array2list(matrix2arrayasis(newmaulinks,1)),xlim=c(25,70))
  cycles2 <- function(links) {
    require(ggm)
    if (any(is.na(links))) stop("missing value/s in links")
    mat <- matrix(0L, max(links), max(links))
    mat[links] <- 1
    lapply(ggm::fundCycles(mat), function(xa) rev(xa[ , 1L]))
  }
  #plotmyarrays(mau[cycles2(newmaulinks)[[1]],],links=conslinks(nrow(mau),open=F))
  sfe1init<-spsample(Polygon(rbind(mauinit[cycles2(newmaulinks)[[1]],],mauinit[cycles2(newmaulinks)[[1]],][1,])),nadd,type="regular")
  sr1init<-sfe1init@coords
  
  sfe1fin<-spsample(Polygon(rbind(maufin[cycles2(newmaulinks)[[1]],],maufin[cycles2(newmaulinks)[[1]],][1,])),nadd,type="regular")
  sr1fin<-sfe1fin@coords
  
  if(interpfin==T){
    Mfin<-centrosfin
    fitfin<- Tps(Mfin, vals) 
    predfin<-predict(fitfin,sr1fin)
  }else{predfin<-c("no predfin")}
  if(interpinit==T){
    Minit<-centrosinit
    fitinit<- Tps(Minit, vals) 
    predinit<-predict(fitinit,sr1init)
  }else{predinit<-c("no predinit")}
  
  
  if(is.null(zlim)==T){zlim<-c(min(na.omit(vals)),max(na.omit(vals)))}else{zlim=zlim}
  
  
  
  if(st==T){
    par(mfrow=c(1,2))
    plotmyarrays(shapes::abind(newinit,newfin)[,,1],links=links,lwd=lwd,col=collinks1,txt=F,xlim=range(rbind(newinit,newfin)[,1]),ylim=range(rbind(newinit,newfin)[,2]),cex=0,xaxt="n",yaxt="n",axes=F)
    if(interpinit==T){
      par(new=T)
      imainit<-image(xyz2img(cbind(sr1init,predinit),tolerance=0.1),asp=1,xlim=range(rbind(newinit,newfin)[,1]),ylim=range(rbind(newinit,newfin)[,2]),col=colorRampPalette(colors)(n = length(vals)),xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)
      
    }
    title(tit1)
    if(plotboth==T){points(fin,cex=0);lineplot(fin,links,col=collinks2,lwd=lwd)}
    if(plotgrid1==T){deformGrid2d(init,fin,ngrid=40,cex1=0,cex2=0,lwd=1,col2=4,add=T,lines=F)}
    if(plotri1==T){plotri2d(newinit,thetria,col=coltri1)}
    for(i in 1:nrow(thetria)){
      if(discretepsl==T){ccpsl<-matrix(psltrisob$res[i,12:20],ncol=3)}else{
        ccpsl<-cbind(rbind(t(objac2d$jac1[[i]])%*%objac2d$jac1[[i]],0),0)
      }
      if(elli==T){
        
        
        mat<-expm::sqrtm(matrix(ccpsl,ncol=3))
        dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]
        centro<-centrosinit[i,]
        defelliadd(dat,centro,mat,col=colelli[i],mag=mag)
      }
      
      if(discretepsl==T){
        
        psl1x<-psltrisob$res[i,3]
        psl1y<-psltrisob$res[i,4]
        psl2x<-psltrisob$res[i,6]
        psl2y<-psltrisob$res[i,7]
        
        mag1=sqrt(psltrisob$res[i,1])*mag
        mag2=sqrt(psltrisob$res[i,2])*mag
        
        
        sign1<-sign(psltrisob$res[i,1]-1)
        code1<-ifelse(sign1>0,2,1)
        
        
        sign2<-sign(psltrisob$res[i,2]-1)
        code2<-ifelse(sign2>0,2,1)
        
        if(slcoleig==F){col1<-ifelse(sign1>0,cold,colp)}else{col1<-numcol(sqrt(psltrisob$res[,1]),colors=colorsleig,zlim=range(sqrt(psltrisob$res[,1:2])))[i]}
        if(slcoleig==F){col2<-ifelse(sign2>0,cold,colp)}else{col2<-numcol(sqrt(psltrisob$res[,2]),colors=colorsleig,zlim=range(sqrt(psltrisob$res[,1:2])))[i]}
        
      }else{
        
        eigccpsl<-eigen(ccpsl)
        psl1eig<-Re(eigccpsl$vectors)[,1]
        psl1x<-(psl1eig/(norm(psl1eig,type="F")))[1]
        psl1y<-(psl1eig/(norm(psl1eig,type="F")))[2]
        
        psl2eig<-Re(eigccpsl$vectors)[,2]
        psl2x<-(psl2eig/(norm(psl2eig,type="F")))[1]
        psl2y<-(psl2eig/(norm(psl2eig,type="F")))[2]
        
        mag1<-sqrt(Re(eigccpsl$values))[1]*mag
        mag2<-sqrt(Re(eigccpsl$values))[2]*mag
        
        
        sign1<-sign(Re(eigccpsl$values)[1]-1)
        code1<-ifelse(sign1>0,2,1)
        
        sign2<-sign(Re(eigccpsl$values)[2]-1)
        code2<-ifelse(sign2>0,2,1)
        
        if(slcoleig==F){col1<-ifelse(sign1>0,cold,colp)}else{col1<-numcol(sqrt(objac2d$cjacps[,1]),colors=colorsleig,zlim=range(sqrt(objac2d$cjacps[,1:2])))[i]}
        if(slcoleig==F){col2<-ifelse(sign2>0,cold,colp)}else{col2<-numcol(sqrt(objac2d$cjacps[,2]),colors=colorsleig,zlim=range(sqrt(objac2d$cjacps[,1:2])))[i]}
        
      }
      
      if(elliax==T){
        
        if(plotelliax2==T){ 
          arrows(centrosinit[i,1],centrosinit[i,2],(centrosinit[i,1]-psl2x*mag2),(centrosinit[i,2]-psl2y*mag2),col=col2,length=la,code=code2,lwd=lwda)
          arrows(centrosinit[i,1],centrosinit[i,2],(centrosinit[i,1]+psl2x*mag2),(centrosinit[i,2]+psl2y*mag2),col=col2,length=la,code=code2,lwd=lwda)
        }
        if(plotelliax1==T){ 
          arrows(centrosinit[i,1],centrosinit[i,2],(centrosinit[i,1]-psl1x*mag1),(centrosinit[i,2]-psl1y*mag1),col=col1,length=la,code=code1,lwd=lwda)
          arrows(centrosinit[i,1],centrosinit[i,2],(centrosinit[i,1]+psl1x*mag1),(centrosinit[i,2]+psl1y*mag1),col=col1,length=la,code=code1,lwd=lwda)
        }
        
      }
      
    } 
  }
  
  pslfin<-Re(pslinfin(psltrisob))
  mag22=rep(mag,nrow(centrosfin))
  mag11=rep(mag,nrow(centrosfin))
  
  plotmyarrays(shapes::abind(newinit,newfin)[,,2],col=collinks2,links=links,lwd=lwd,txt=F,xlim=range(rbind(newinit,newfin)[,1]),ylim=range(rbind(newinit,newfin)[,2]),cex=0,xaxt="n",yaxt="n",axes=F)
  if(interpfin==T){
    par(new=T)
    imafin<-image(xyz2img(cbind(sr1fin,predfin),tolerance=0.1),asp=1,xlim=range(rbind(newinit,newfin)[,1]),ylim=range(rbind(newinit,newfin)[,2]),col=colorRampPalette(colors)(n = length(predfin)),xaxs="r",yaxs="r",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",zlim=zlim)
  }
  
  title(tit2)
  if(plotboth==T){points(init,cex=0);lineplot(init,links,col=collinks1,lwd=lwd)}
  
  if(plotgrid2==T){deformGrid2d(init,fin,ngrid=40,cex1=0,cex2=0,lwd=1,col2=4,add=T,lines=F)}
  if(plotri2==T){plotri2d(newfin,thetria,col=2)}
  
  princ=NULL
  for(i in 1:nrow(thetria)){
    if(elli==T){
      
      
      if(discretepsl==T){
        mat<-Re((matrix(psltrisob$res[i,27:35],ncol=3)))}else{mat<-cbind(rbind(objac2d$jac1[[i]],0),0)}
      dat<-centershapes(cbind(circle2(plot=F,radius=1),0))[,,1]*mag
      centro<-centrosfin[i,]
      defelliadd(dat,centro,mat,col=colelli[i],after=F,border=elliborder)
    }
    
    if(discretepsl==T){
      psl1x<-pslfin[i,1]
      psl1y<-pslfin[i,2]
      psl2x<-pslfin[i,4]
      psl2y<-pslfin[i,5]
      
      
      sign1<-sign(psltrisob$res[i,1]-1)
      code1<-ifelse(sign1>0,2,1)
      sign2<-sign(psltrisob$res[i,2]-1)
      code2<-ifelse(sign2>0,2,1)
      
      if(slcoleig==F){col1<-ifelse(sign1>0,cold,colp)}else{col1<-numcol(sqrt(psltrisob$res[,1]),colors=colorsleig,zlim=range(sqrt(psltrisob$res[,1:2])))[i]}
      if(slcoleig==F){col2<-ifelse(sign2>0,cold,colp)}else{col2<-numcol(sqrt(psltrisob$res[,2]),colors=colorsleig,zlim=range(sqrt(psltrisob$res[,1:2])))[i]}
      
      princi<-which.max(abs(1-psltrisob$res[i,1:2]))
      
    }else{
      ccpsl<-cbind(rbind(t(objac2d$jac1[[i]])%*%objac2d$jac1[[i]],0),0)
      
      eigccpsl<-eigen(ccpsl)
      psl1eig<-Re(eigccpsl$vectors)[,1]
      psl1x<-(psl1eig/(norm(psl1eig,type="F")))[1]
      psl1y<-(psl1eig/(norm(psl1eig,type="F")))[2]
      psl1<-pslinfin1(c(psl1x,psl1y),objac2d$jac1[[i]])
      psl1x<-psl1[1]
      psl1y<-psl1[2]
      
      psl2eig<-Re(eigccpsl$vectors)[,2]
      psl2x<-(psl2eig/(norm(psl2eig,type="F")))[1]
      psl2y<-(psl2eig/(norm(psl2eig,type="F")))[2]
      psl2<-pslinfin1(c(psl2x,psl2y),objac2d$jac1[[i]])
      psl2x<-psl2[1]
      psl2y<-psl2[2]
      
      
      sign1<-sign(Re(eigccpsl$values)[1]-1)
      code1<-ifelse(sign1>0,2,1)
      sign2<-sign(Re(eigccpsl$values)[2]-1)
      code2<-ifelse(sign2>0,2,1)
      
      if(slcoleig==F){col1<-ifelse(sign1>0,cold,colp)}else{col1<-numcol(sqrt(objac2d$cjacps[,1]),colors=colorsleig,zlim=range(sqrt(objac2d$cjacps[,1:2])))[i]}
      if(slcoleig==F){col2<-ifelse(sign2>0,cold,colp)}else{col2<-numcol(sqrt(objac2d$cjacps[,2]),colors=colorsleig,zlim=range(sqrt(objac2d$cjacps[,1:2])))[i]}
      
      princi<-which.max(abs(1-Re(eigccpsl$values[1:2])))
      if(princi<2){
        pslprincx<-psl1x
        pslprincy<-psl1y
        colprinc<-col1
        codeprinc<-code1
      }else{
        pslprincx<-psl2x
        pslprincy<-psl2y
        colprinc<-col2
        codeprinc<-code2
      }
      
      
    }
    if(elliax==T){
      
      if(plotelliax2==T){ 
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]-psl2x*mag22[i]),(centrosfin[i,2]-psl2y*mag22[i]),col=col2,length=la,code=code2,lwd=lwda)
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]+psl2x*mag22[i]),(centrosfin[i,2]+psl2y*mag22[i]),col=col2,length=la,code=code2,lwd=lwda)
      }
      
      if(plotelliax1==T){ 
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]-psl1x*mag11[i]),(centrosfin[i,2]-psl1y*mag11[i]),col=col1,length=la,code=code1,lwd=lwda)
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]+psl1x*mag11[i]),(centrosfin[i,2]+psl1y*mag11[i]),col=col1,length=la,code=code1,lwd=lwda)
      }
      
      if(onlyprinc==T){
        
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]-pslprincx*mag22[i]),(centrosfin[i,2]-pslprincy*mag22[i]),col=colprinc,length=la,code=codeprinc,lwd=lwda)
        arrows(centrosfin[i,1],centrosfin[i,2],(centrosfin[i,1]+pslprincx*mag22[i]),(centrosfin[i,2]+pslprincy*mag22[i]),col=colprinc,length=la,code=codeprinc,lwd=lwda)
        
        
      }
      
      
    }
    
    princ<-c(princ,princi)
  }
  
  if(is.null(zlim)==T){zliml<-c(min(na.omit(vals)),max(na.omit(vals)))}else{zliml=zlim}
  
  theseq<-seq(zliml[1],zliml[2],length.out=1000)
  theseq<-theseq[order(theseq)]
  if(plotlegend==T){plot(cbind(theseq,0),col=numcol(theseq),pch=15,cex=4,xaxt="n",yaxt="n",xlab="",ylab="",axes=F)
    axis(1)}
  if(plotdens==T){plot(density(vals))}
  out=list(vals=vals,psltrisob=psltrisob,pslfin=pslfin,init=init,fin=fin,newinit=newinit,newfin=newfin,trisob=tris,gridob=tpsgrid,sign1=sign1,code1=code1,sign2=sign2,code2=code2,centrosinit=centrosinit,centrosfin=centrosfin,tria=thetria,jac2norms=objac2d$jac2norms,areasinit=areasinit,areasfin=areasfin,areasratio=areasratio,sens=objac2d$sens,totsen=objac2d$sentot,snl=psltrisob$res[,40],beb=beb,tpsjacob=objac2d,pslenes=pslenes,predinit=predinit,predfin=predfin,sr1init=sr1init,sr1fin=sr1fin,princ=princ,colors=colors,colorsleig=colorsleig)
  out
}

triang2dto<-function(mat1,mat2,links=NULL,ar=70,PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,extpol=F){
  ar1<-areasip(mat1,links=links,graph=F)
  ar2<-areasip(mat1,links=links,a=ar1$area/ar,graph=F,PB=PB,PA=PA,S=S,SB=SB,H=H,V=V,q=q,Y=Y,j=j,D=D,St=St,Q=Q,extpol=extpol)
  nonorig<-ar2$triangob$P[(nrow(mat1)+1):nrow(ar2$triangob$P),]
  tria<-ar2$triangob$T
  newp<-applytps(mat1,mat2,nonorig)
  newmat2=rbind(mat2,newp)
  centrosnewmat1<-centrotri(ar2$triangob$P,tria)
  centrosnewmat2<-centrotri(newmat2,tria)
  out=list(mat1=mat1,newmat1=ar2$triangob$P,newmat2=newmat2,tria=tria,areasipob=ar2,centrosnewmat1=centrosnewmat1,centrosnewmat2=centrosnewmat2)
  out
}

areasip<-function(matrix,links=NULL,PB=NA,PA=NA,S=NA,SB=NA,H=NA,V=3,a=NULL,q=NULL,Y=FALSE,j=FALSE,D=FALSE,St=Inf,Q=TRUE,graph=T,extpol=F){
  if(!is.null(links)){warning("Links **must** identify, among other structures, a closed contour")}
  library(geometry)
  library(tripack)
  library(geoR)
  library(gstat)
  library(sp)
  library(fields)
  library(RTriangle)
  mate<-matrix
  if(!is.null(links)&is.na(S)==T){S<-list2matrix(links)}
  
  if(extpol==T){
    step1<-pslg(mate,S=S,PB=PA,PA=PA,SB=SB,H=NA)
    posimpnoh<-RTriangle::triangulate(step1,V=V,a=NULL,S=St,q=q,Y=Y,j=j,D=D,Q=Q)
    step2<-matrix[unique(posimpnoh$E[ posimpnoh$EB%in%c(1)]),]
    maulinks<-posimpnoh$E[posimpnoh$EB%in%c(1),]
    newmaulinks<-maulinks
    newmaulinks[,1]<-c(1:nrow(maulinks))
    newmaulinks[,2]<-match(maulinks[,2],maulinks[,1])
    cycles2 <- function(links) {
      require(ggm)
      if (any(is.na(links))) stop("missing value/s in links")
      mat <- matrix(0L, max(links), max(links))
      mat[links] <- 1
      lapply(ggm::fundCycles(mat), function(xa) rev(xa[ , 1L]))
    }
    
    mypol<-step2[cycles2(newmaulinks)[[1]],]}else{mypol<-c("you do not want external polygon")}
  
  
  p<-pslg(mate,PB=PB,PA=PA,S=S,SB=SB,H=H)
  au<-RTriangle::triangulate(p,V=1,a=a,q=q,Y=Y,j=j,D=D,S=St,Q=Q)
  if(graph==T){
    plot(au,asp=1)
    lineplot(mate,links,col=3)
    points(au$P,pch=19)
    points(mate,pch=21,cex=2)
  }
  centros<-NULL
  areas<-NULL
  for(i in 1:nrow(au$T)){
    centrosi<-apply(au$P[au$T[i,],],2,mean)
    centros<-rbind(centros,centrosi)
    areasi<-convhulln(au$P[au$T[i,],],option="FA")$vol
    areas<-c(areas,areasi)
  }
  if(graph==T){points(centros,col=2)}
  M<-centros
  
  
  
  thecentr<-centroids(mate)
  dasomm<-NULL
  for(i in 1:nrow(M)){
    dasommi<-dist(rbind(thecentr,M[i,]))^2*areas[i]
    dasomm<-c(dasomm,dasommi)
  }
  ip<-sum(dasomm)
  are<-sum(areas)
  
  deltri<-au$T
  ptri<-au$P
  origcentros<-rbind(mate,M)
  
  out<-list(area=are,areas=areas,ip=ip,centros=centros,deltri=deltri,ptri=ptri,origcentros=origcentros,triangob=au,ext=mypol)
  out
}

applytps<-function(init,fin,init2,doopa=F,meth=c("mor","dm")){
  library(Morpho)
  if(doopa==T){
    theopa<-rotonto(init,fin,scale=F,reflection=F)
    init<-init
    fin<-theopa$yrot}
  if(meth=="dm"){
    esempio<-tpsdry2(init,fin,doopa=F)
    out<-as.matrix(rep(1,nrow(init2)))%*%esempio$ct+init2%*%esempio$at+sigm.new(init,init2)%*%esempio$W
  }else{
    #########  copied from Stefan Schlager code
    trafo<-computeTransform(fin,init,type = "tps")
    .fx <- function(refmat,M,coefs,time=TRUE,threads=1) {
      M <- cbind(1,M)
      splM <- .Call("tpsfx",refmat,  M, t(coefs),threads) 
      return(splM)
    }
    out <- .fx(trafo$refmat, init2, trafo$coeff, threads = 1)
    
  }
  as.matrix(out)
}

enespsl<-function(psltrisob){
  n=nrow(psltrisob$res)
  k<-ncol(psltrisob$init)
  if(mean(psltrisob$init[,3])==0){k<-2}
  trc<-NULL
  cdets<-NULL
  sl<-NULL
  for(i in 1:n){
    ci<-matrix(psltrisob$res[i,12:20],ncol=3)
    ci<-ci[1:k,1:k]
    trci<-matrix.trace(ci)
    detci<-det(ci)
    trc<-c(trc,trci)
    cdets<-c(cdets,detci)
    alphainit=1
    betainit=1
    E=0.5*(ci-diag(k))#### green lagrangian strain
    
    E2<-E%*%E
    trE<-matrix.trace(E)
    trE2<-matrix.trace(E2)
    lucioi<-(alphainit/2)*(trE^2+betainit*trE2)
    
    
    sl<-c(sl,lucioi)
  }
  enes<-cbind(trc,cdets,sl)
  enes
}

ar3dtri<-function(matrix){
  a<-matrix[1,]
  b<-matrix[2,]
  c<-matrix[3,]
  ab<-dist(rbind(a,b))[1]
  bc<-dist(rbind(b,c))[1]
  ca<-dist(rbind(c,a))[1]
  s<-(ab+bc+ca)/2
  ar<-sqrt((s*(s-ab)*(s-bc)*(s-ca)))
  ar
}

posfac<-function(factor){
  positions<-NULL
  for(k in 1:max(as.numeric(factor))){
    factor<-factor(factor,levels=unique(factor))
    positioni<-which(as.numeric(factor)==k)
    positions<-c(positions,list(positioni))
  }
  names(positions)<-levels(factor)
  sequ<-NULL
  for(i in 1:max(length(positions))){
    seqi<-c(1:length(positions[[i]]))
    sequ<-c(sequ,list(seqi))
  }
  
  pure<-rep(NA,length(unlist(positions)))
  for(j in 1:max(length(positions))){
    pure<-replace(pure,positions[[j]],c(sequ[[j]]))
  }
  
  pure}

plotri2d<-function(mat,tria,col=1){
  for(i in 1:nrow(tria)){
    plotri12d(mat[tria[i,],],col=col)
  }
}

plotri12d<-function(mat32d,col=1){
  lineplot(mat32d,conslinks(3,open=F),col=col)
}

gradientLegend2<-function (valRange, color = "topo", nCol = 30, pos = 0.5, 
                           side = 4, length = 0.25, depth = 0.05, inside = TRUE, coords = FALSE, 
                           pos.num = NULL, n.seg = 2, border.col = "black", dec = NULL, 
                           fit.margin = TRUE,notvals=NULL) {
  loc <- c(0, 0, 0, 0)
  if (is.null(pos.num)) {
    if (side %in% c(1, 3)) {
      pos.num = 3
    }
    else {
      pos.num = side
    }
  }
  if (length(pos) == 1) {
    pos.other <- ifelse(side > 2, 1, 0)
    if (side %in% c(1, 3)) {
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos - 0.5 * length, pos.other - 
                           switch * depth, pos + 0.5 * length, pos.other + 
                           (1 - switch) * depth), side = c(side, 2, side, 
                                                           2))
    }
    else if (side %in% c(2, 4)) {
      switch <- ifelse(inside, 0, 1)
      switch <- ifelse(side > 2, 1 - switch, switch)
      loc <- getCoords(c(pos.other - switch * depth, pos - 
                           0.5 * length, pos.other + (1 - switch) * depth, 
                         pos + 0.5 * length), side = c(1, side, 1, side))
    }
  }
  else if (length(pos) == 4) {
    if (coords) {
      loc <- pos
    }
    else {
      loc <- getCoords(pos, side = c(1, 2, 1, 2))
    }
  }
  mycolors <- c()
  if (length(color) > 1) {
    mycolors <- color
  }
  else if (!is.null(nCol)) {
    if (color == "topo") {
      mycolors <- topo.colors(nCol)
    }
    else if (color == "heat") {
      mycolors <- heat.colors(nCol)
    }
    else if (color == "terrain") {
      mycolors <- terrain.colors(nCol)
    }
    else if (color == "rainbow") {
      mycolors <- rainbow(nCol)
    }
    else {
      warning("Color %s not recognized. A palette of topo.colors is used instead.")
      mycolors <- topo.colors(nCol)
    }
  }
  else {
    stop("No color palette provided.")
  }
  vals <- seq(min(valRange), max(valRange), length = length(mycolors))
  if (!is.null(dec)) {
    vals <- round(vals, dec[1])
  }
  im <- as.raster(mycolors[matrix(1:length(mycolors), ncol = 1)])
  ticks <- c()
  if (side%%2 == 1) {
    rasterImage(t(im), loc[1], loc[2], loc[3], loc[4], col = mycolors, 
                xpd = T)
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
    ticks <- seq(loc[1], loc[3], length = n.seg)
    segments(x0 = ticks, x1 = ticks, y0 = rep(loc[2], n.seg), 
             y1 = rep(loc[4], n.seg), col = border.col, xpd = TRUE)
  }
  else {
    rasterImage(rev(im), loc[1], loc[2], loc[3], loc[4], 
                col = mycolors, xpd = T)
    rect(loc[1], loc[2], loc[3], loc[4], border = border.col, 
         xpd = T)
    ticks <- seq(loc[2], loc[4], length = n.seg)
    segments(x0 = rep(loc[1], n.seg), x1 = rep(loc[3], n.seg), 
             y0 = ticks, y1 = ticks, col = border.col, xpd = TRUE)
  }
  determineDec <- function(x) {
    out = max(unlist(lapply(strsplit(x, split = "\\."), 
                            function(y) {
                              return(ifelse(length(y) > 1, nchar(gsub("^([^0]*)([0]+)$", 
                                                                      "\\1", as.character(y[2]))), 0))
                            })))
    return(out)
  }
  labels = sprintf("%f", seq(min(valRange), max(valRange), 
                             length = n.seg))
  if (is.null(dec)) {
    dec <- min(c(6, determineDec(labels)))
  }
  eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                            paste("%.", dec, "f", sep = ""))))
  if (pos.num == 1) {
    if (fit.margin) {
      lab.height = max(strheight(labels)) * 0.8
      max.pos = getFigCoords()[3]
      if ((max.pos - loc[2]) < lab.height) {
        warning("Increase bottom margin, because labels for legend do not fit.")
      }
    }
    text(y = loc[2], x = ticks, labels = seq(min(valRange), 
                                             max(valRange), length = n.seg), col = border.col, 
         pos = 1, cex = 0.8, xpd = T)
  }
  else if (pos.num == 2) {
    if (fit.margin) {
      checkagain = TRUE
      while (checkagain == TRUE) {
        lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) * 
          0.8
        min.pos = getFigCoords()[1]
        if ((loc[1] - min.pos) < lab.width) {
          if (!is.null(dec)) {
            dec = max(c(0, dec - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          else {
            tmp = max(unlist(lapply(strsplit(labels, 
                                             split = "\\."), function(x) {
                                               return(ifelse(length(x) > 1, nchar(x[2]), 
                                                             0))
                                             })))
            dec = max(c(0, tmp - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                                    paste("%.", dec, "f", sep = ""))))
        }
        else {
          checkagain = FALSE
        }
      }
    }
    text(y = ticks, x = loc[1], labels = labels, pos = 2, 
         cex = 0.8, col = border.col, xpd = T)
  }
  else if (pos.num == 3) {
    if (fit.margin) {
      lab.height = max(strheight(labels)) * 0.8
      max.pos = getFigCoords()[4]
      if ((max.pos - loc[4]) < lab.height) {
        warning("Increase top margin, because labels for legend do not fit.")
      }
    }
    text(y = loc[4], x = ticks, labels = seq(min(valRange), 
                                             max(valRange), length = n.seg), col = border.col, 
         pos = 3, cex = 0.8, xpd = T)
  }
  else if (pos.num == 4) {
    if (fit.margin) {
      checkagain = TRUE
      while (checkagain == TRUE) {
        lab.width = (max(strwidth(labels)) + 0.5 * par()$cxy[1]) * 
          0.8
        max.pos = getFigCoords()[2]
        if ((max.pos - loc[3]) < lab.width) {
          if (!is.null(dec)) {
            dec = max(c(0, dec - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          else {
            tmp = max(unlist(lapply(strsplit(labels, 
                                             split = "\\."), function(x) {
                                               return(ifelse(length(x) > 1, nchar(x[2]), 
                                                             0))
                                             })))
            dec = max(c(0, tmp - 1))
            if (dec == 0) {
              warning("Decimal set to 0 (dec=0), but the labels still don't fit in the margin. You may want to add the color legend to another side, or increase the margin of the plot.")
              checkagain = FALSE
            }
          }
          eval(parse(text = sprintf("labels = sprintf('%s', round(seq(min(valRange), max(valRange), length = n.seg), dec) )", 
                                    paste("%.", dec, "f", sep = ""))))
        }
        else {
          checkagain = FALSE
        }
      }
    }
    text(y = ticks, x = loc[3], labels = labels, pos = 4, 
         col = border.col, cex = 0.8, xpd = T)
  }
  if(is.null(notvals)==F){
    value    <- notvals
    pos      <- pos    
    length   <- length
    depth    <- depth
    # calculations:
    minval <- min(valRange)
    maxval <- max(valRange)
    value.prop.scale <- (value - minval) / (maxval - minval) # proportion on gradient
    value.prop <- value.prop.scale * length + (pos - 0.5*length) # proportion on y-axis
    value.y.coord <- getCoords(value.prop, side=side, input='p')
    value.x0 <- getCoords(1, side=1, input='p')
    value.x1 <- getCoords(1-depth, side=1, input='p')
    # draw line in red:
    segments(x0=value.x0, x1=value.x1, y0=value.y.coord, y1=value.y.coord, col=1)
    text(loc[3],value.y.coord,labels=as.character(value),pos = 4, col = border.col, cex = 0.8, xpd = T)
  }
}

centermesh<-function(object){
  oldver<-object$vb[1:3,]
  thetri<-object$it
  newver<-t(centershapes(t(oldver))[,,1])
  thesurf=list(vb=rbind(newver,1),it=thetri)
  class(thesurf)<-"mesh3d"
  thesurf
}
