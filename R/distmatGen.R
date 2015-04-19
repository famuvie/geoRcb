distmatGen<-function(pts,costsurf){
  #generate cost based distances amoung spatialpoints
  #   library(geoRcb)
  #   data(noise)
  #   library(gdistance)
  #   library(rgdal)

  #set resolution to min distance between points
  coords<-coordinates(obs)
  coords.shift<-rbind(coords[-1,],NA)
  coords.diff<-coords-coords.shift
  minres<-min(c(abs(coords.diff[,1]),abs(coords.diff[,2])),na.rm=T)
  diag<-sqrt(minres^2+minres^2)

  r<-raster(nrows=36,ncols=36)
  r<-setExtent(r,extent(obs))
  res(r)<-minres
  r<-setValues(r,rep(1,ncell(r)))
  projection(r)<-proj4string(obs)
  tr<-transition(r,function(x) 1/mean(x),8)

  distmat<-list()
  for(i in 1:nrow(obs)){
    distmat[[i]]<-extract(accCost(tr,coordinates(obs)[i,]),obs)*diag
  }
  do.call("rbind",distmat)

}

#
#   res<-cbind(t(matrix(res[1,],nrow=1)),t(matrix(dd.distmat[1,],nrow=1)))
#   plot(res[,1],res[,2])
#   abline(0,1)
#
