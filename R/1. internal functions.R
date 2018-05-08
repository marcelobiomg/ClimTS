# internal functions

################# Find the closest values to fill NAs
# from file 2.1 Annual Variables Migratory
# http://stackoverflow.com/questions/27562076/if-raster-value-na-search-and-extract-the-nearest-non-na-pixel

subs.na <- function(r, xy){
  sampled <- apply(X = xy, MARGIN = 1, FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
  if(any(is.na(sampled))==T){
    sampled <- apply(X = xy, MARGIN = 1, FUN = function(xy) getValues(r)[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
  }
  unlist(lapply(sampled, function(x) mean(x)))
}

# Find the closest values to fill NAs across layers of a stack
ClosestAV <- function(spdf, r, max.dist.km = max.dist.km){
  xy <- matrix(coordinates(spdf)[is.na(spdf@data[,1]),], ncol=2) #}
  for(i in 1:nlayers(r)){
    # print(paste("layer ", i))
    extrct <- subs.na(r[[i]], xy, max.dist.km)
    if(is.null(extrct)) extrct <- NA
    if(all(is.na(extrct))) extrct <- NA
    spdf@data[is.na(spdf@data[,i]),i] <-  extrct
  }
  return(spdf)
}



################ considerando meses no wintering range, weighted by distance from breeding record
# from file 2.3 Annual Variables Migratory batch
#### function to compute weighted mean across raster layers
wDfP <- function(r, w, na.rm=T){
  w.r <- integer(nlayers(r))
  for(i in 1:length(w.r)){
    w.r[i] <- weighted.mean(getValues(r[[i]]), getValues(w), na.rm=T)
  }
  return(w.r)
}

#### função weighted mean of raster
f.wr <- function(r, xy, sp.m.p, m.w) {
  require(rgeos)
  require(raster)

  # crop raster to sp.shp
  c.raster <- r[[m.w]]*sp.m.p

  df.r <- data.frame(matrix(ncol = length(m.w), nrow = length(xy)))

  for(i in 1:length(xy)){
    # create weight raster based on distance from point # inverse of distance from points
    DfP <- distanceFromPoints(sp.m.p, xy[i,])
    w.xy <- abs(DfP/max(getValues(DfP)) - 1)*sp.m.p
    df.r[i,] <- wDfP(c.raster, w.xy)
  }
  return(df.r)
}

f.wr(r, xy, sp.m.p, m.w)


