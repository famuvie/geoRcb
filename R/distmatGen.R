#'@name distmatGen
#'@title Generate a cost-based distance matrix among SpatialPoints
#'@export
#'@import gdistance
#'@details
#'@param pts A SpatialPointsDataFrame with a defined projection
#'@param costsurf optional cost surface Raster
#'@details If a cost surface raster is not supplied a uniform (cost 1) grid will be generated to match to bounding box of pts
#'@example #res<-distmatGen(obs)
#'plt<-cbind(res[1,],dd.distmat[1,])
#'plot(plt)
#'abline(0,1,col="red")

distmatGen<-function(pts,...){

  coords<-coordinates(pts)

  if(!exists("costsurf")){
  ocoords<-cbind(coords[order(coords[,1]),1],coords[order(coords[,2]),2])
  ocoords.shift<-rbind(ocoords[-1,],NA)
  ocoords.diff<-ocoords-ocoords.shift
  ocoords.diff<-ocoords.diff[apply(ocoords.diff,1,function(x) sum(x!=0))==2,]
  #hist(abs(ocoords.diff[,2]))
  minres<-mean(c(abs(ocoords.diff[,1]),abs(ocoords.diff[,2])),na.rm=T)
  diag<-sqrt(minres^2+minres^2)

  r<-raster(nrows=36,ncols=36)
  r<-setExtent(r,extent(pts))
  res(r)<-minres
  r<-setExtent(r,extent(pts))
  r<-setValues(r,rep(1,ncell(r)))
  projection(r)<-proj4string(pts)
  }else{
    r<-costsurf
    diag<-sqrt(res(r)^2+res(r)^2)
  }

  tr<-transition(r,function(x) 1/mean(x),8)

  distmat<-list()
  for(i in 1:nrow(pts)){
    distmat[[i]]<-extract(accCost(tr,coordinates(pts)[i,]),pts)*diag
  }
  do.call("rbind",distmat)
}

