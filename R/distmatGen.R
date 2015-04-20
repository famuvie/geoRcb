#'@name distmatGen
#'@title Generate a cost-based distance matrix among SpatialPoints
#'@export
#'@import gdistance
#'@details
#'@param pts A SpatialPointsDataFrame with a defined projection
#'@param costsurf cost surface Raster
#'@param ret specify whether to return distances between obs-obs ("o"), obs-loc ("l"), or both ("b")
#'@details If a cost surface raster is not supplied a uniform (cost 1) grid will be generated to match to bounding box of pts
#'@example
#'data(noise)
#'r<-raster(nrows=40,ncols=40,resolution=1)
#'r<-setExtent(r,extent(obs),keepres=T)
#'aggfac<-5
#'costsurf<-rasterize(malilla,r,background=aggfac,field=rep(10000,length(polygons(malilla))))
#'costsurf<-reclassify(costsurf,c(aggfac+1,Inf,NA))
#'costsurf<-aggregate(costsurf,aggfac,expand=T,fun=max,na.rm=T)
#'
#'res<-distmatGen(obs,costsurf,ret="l")
#'plot(dd.distmat[,51],res[,51],ylim=c(0,600),xlim=c(0,600))
#'abline(0,1,col="red")
#'

distmatGen<-function(pts,costsurf,ret){

  #pts<-obs
  #ret="o"
  r<-costsurf

  if(!extent(r)>=extent(pts)){
    xn<-floor(extent(pts)[1]-1)
    xx<-ceiling(extent(pts)[2]+1)
    yn<-floor(extent(pts)[3]-1)
    yx<-ceiling(extent(pts)[4]+1)
    ext<-unlist(lapply(list(xn,xx,yn,yx),function(x) round(x,0)))
    r<-extend(r,extent(ext))
  }

  tr<-transition(r,function(x) 1/mean(x),8)

  if(ret=="b"|ret=="o"){
  oret<-diag(0,nrow=nrow(pts))
  oret.c<-matrix(costDistance(tr,coordinates(pts)))
  oret[lower.tri(oret)]<-oret.c
  oret<-oret+t(oret)
  }

  if(ret=="b"|ret=="l"){
    tocoords<-coordinates(xyFromCell(tr,which(values(r)==res(r)),spatial=TRUE))
    lret<-diag(0,nrow=nrow(tocoords))
    lret.c<-matrix(costDistance(tr,tocoords))#CANNOT ALLOCATE ENOUGH MEMORY
    lret[lower.tri(lret)]<-lret.c
    lret<-lret+t(lret)
    }

  if(ret=="b"){
    list(oret,lret)
  }
  if(ret=="o"){
    oret
  }else{
    lret
  }

}


##reimplementation of gdistance costDistance function (sparse matrices solve memory allocation error?)
#path distances are always less than or equal to straight-line distances

# cd2<-function(x,fromCoords){
#   #x<-tr
#   #fromCoords<-tocoords
#
#   fromCoords <- gdistance:::.coordsToMatrix(fromCoords)
#   fromCells <- cellFromXY(x, fromCoords)
#
#   if(!all(!is.na(fromCells))){
#     warning("some coordinates not found and omitted")
#     fromCells <- fromCells[!is.na(fromCells)]
#   }
#
#   costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(fromCoords[,1]))
#   rownames(costDist) <- rownames(fromCoords)
#   colnames(costDist) <- rownames(fromCoords)
#
#   if(isSymmetric(transitionMatrix(x))) {m <- "undirected"} else{m <- "directed"}
#
# }
