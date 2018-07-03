################### Função que extrai presença sazonal dos shapefiles de várias spp e
################### salva mapas de pontos sobre a distribuição
# from file 2.2 Annual Variables Migratory
a.var.seas <- function(r, OBS, sp.shp.path, basemap=NULL, Maps= "save") {
  require(sp)
  require(raster)
  OBS$SeasPres <- NA
  OBS$AproxSeasPres <- F

  # 1. extrair presença sazonal de todas as espécies para o OBS
  # 1.1 criar lista de espécies observadas
  spp.list <- sort(unique(OBS$species))
  # 1.1b criar lista de presença sazonal de acordo com os shapes
  migratory <- vector(length = length(spp.list))
  # 1.2 criar lista de espécies da pasta
  spp.folder <- list.files(sp.shp.path, pattern = "\\.shp$", full.names = FALSE)
  spp.folder.F <- list.files(sp.shp.path, pattern = "\\.shp$", full.names = T)

  spp.found <- grep(paste(sub(" ", "_", spp.list),collapse="|"), spp.folder, value = T)
  # 1.2b testar se todas as spp possuem shapefiles
  if(length(spp.found) < length(spp.list)) {
    #spp.found <- grep(paste(sub(" ", "_", spp.list),collapse="|"), spp.folder, value = T)
    print("Shapefiles for some species from dataset are not present in the shapefile folder. See below:")
    spp.not.found <- spp.list[is.na(pmatch(sub(" ", "_", spp.list), spp.found))]
    print(spp.not.found)
    return(spp.not.found)
    stop("Operation will not continue")
  }

  # 1.3 criar Spatial Data Frame
  if(class(OBS) != "SpatialPointsDataFrame"){
    coordinates(OBS) = cbind(OBS[,grep("^lon$|^long$|^longitude$", colnames(OBS), ignore.case = T)],
                             OBS[,grep("^lat$|^latitude$", colnames(OBS), ignore.case = T)])
  }

  # crs(OBS) <- crs(shapefile(paste(sp.shp.path, spp.folder[1], sep= "/")))
  crs(OBS) <- crs(shapefile(spp.folder.F[1]))

  if(Maps == "save") if(dir.exists("sppmaps")==F) dir.create("sppmaps")

  # 1.4 selecionar
  for(i in 1:length(spp.list)){
    # 1.4.1 cada espécie migratória na pasta e..
    sp <- grep(sub(" ", "_", spp.list[i]), spp.folder.F)
    #grep(paste(sub(" ", "_", spp.list),collapse="|"), spp.folder)
    print(paste(spp.list[i], " ", i, " of ", length(spp.list)))
    # sp.shp <- shapefile(paste(sp.shp.path, spp.folder[sp], sep= "/"))
    sp.shp <- shapefile(spp.folder.F[sp])

    # 1.4.2 separar OBS de cada sp
    OBS.spi <- subset(OBS, OBS$species == spp.list[i])

    if(Maps == "save"){
      # change plot extent, to show all points and sp distribution
      x.min <- min(bbox(sp.shp)[1,1],bbox(OBS.spi)[1,1])*1.05
      x.max <- max(bbox(sp.shp)[1,2],bbox(OBS.spi)[1,2])*.95
      y.min <- min(bbox(sp.shp)[2,1],bbox(OBS.spi)[2,1])*1.05
      y.max <- max(bbox(sp.shp)[2,2],bbox(OBS.spi)[2,2])*.95
      x = c(x.min, x.min, x.max, x.max, x.min)
      y = c(y.min, y.max, y.max, y.min, y.min)
      xy <- matrix(c(x,y), ncol=2, byrow=F)

      # cria um arquivo .pdf com nome igual a "i.pdf".
      pdf(paste("sppmaps", paste0(spp.list[i], ".pdf"), sep="/"))
      plot(as.SpatialPolygons.PolygonsList(list(
        Polygons(list(Polygon( xy ) ), "1")) , CRS(proj4string(sp.shp))), border="white", main = spp.list[i] )
      plot(sp.shp[sp.shp$SEASONAL<3,], col="greenyellow", add=T)
      plot(sp.shp[sp.shp$SEASONAL==3,], col="steelblue1", add=T)
      if(!is.null(basemap)) plot(basemap , add=T)
      plot(OBS.spi, add=T, col = "red")
      dev.off()
    } else if(Maps == "plot") {
      # change plot extent, to show all points and sp distribution
      x.min <- min(bbox(sp.shp)[1,1],bbox(OBS.spi)[1,1])*1.05
      x.max <- max(bbox(sp.shp)[1,2],bbox(OBS.spi)[1,2])*.95
      y.min <- min(bbox(sp.shp)[2,1],bbox(OBS.spi)[2,1])*1.05
      y.max <- max(bbox(sp.shp)[2,2],bbox(OBS.spi)[2,2])*.95
      x = c(x.min, x.min, x.max, x.max, x.min)
      y = c(y.min, y.max, y.max, y.min, y.min)
      xy <- matrix(c(x,y), ncol=2, byrow=F)

      plot(as.SpatialPolygons.PolygonsList(list(
        Polygons(list(Polygon( xy ) ), "1")) , CRS(proj4string(sp.shp))), border="white" )
      plot(sp.shp, main = spp.list[i], col="yellow")
      plot(OBS.spi, add=T, col = "red")
    }

    # 1.4.3 extrair presença sazonal de cada espécie para o OBS
    # OBS.spi$SeasPres <- extract(sp.shp[sp.shp$SEASONAL<3, ], OBS.spi)$SEASONAL
    OBS.spi$SeasPres <- over(OBS.spi, sp.shp[sp.shp$SEASONAL<3, ])$SEASONAL
    # 1.4.3b Checar se todos os pontos contém a presença sazonal
    if(any(is.na(OBS.spi$SeasPres) == T)) {
      sp.shp.B <- rasterize(sp.shp[sp.shp$SEASONAL<3, ], r, field="SEASONAL") # rasterize(sp.shp[sp.shp$SEASONAL<3, ], r)
      coord.rep <- matrix(coordinates(OBS.spi[is.na(OBS.spi@data$SeasPres),]), ncol=2)#matrix(coordinates(OBS.spi), ncol=2)[is.na(OBS.spi@data$SeasPres),]
      DfP <- distanceFromPoints(r, coord.rep)
      OBS.spi@data$AproxSeasPres[is.na(OBS.spi@data$SeasPres)] <- T
      OBS.spi@data$SeasPres[is.na(OBS.spi@data$SeasPres)] <- subs.na(sp.shp.B, coord.rep, DfP=DfP )
    }

    # 1.5 Salvar no DF principal
    OBS$SeasPres[OBS$species == spp.list[i]] <- OBS.spi@data$SeasPres
    OBS$AproxSeasPres[OBS$species == spp.list[i]] <- OBS.spi@data$AproxSeasPres

    # 2 spp migratórias - aproveitar que o shape já está carregado e
    # 2.1 checar quais sp são migratórias nos shapes
    migratory[i] <- any(sp.shp$SEASONAL==3)
    # print(any(sp.shp$SEASONAL==3))
    # print(migratory[i])
    #
  }
  migratory <- cbind(data.frame(spp.list), migratory)
  return(list(OBS, migratory))
}

Eggs_tyranni.geo.L <- a.var.seas(grade.r.NW, OBS=Eggs_tyranni.geo, sp.shp.path, Maps="not")





# Insert breeding season information
# DF.spp <- EggsT.Seas.test
# m.bw.sp <- data.frame(species=unique(EggsT.Seas.test$species), A.BG=4, D.BG=7, A.WG=10, D.WG=1)
m.bw.sp <- data.frame(species=unique(EggsT.Seas$species), A.BG=4, D.BG=7, A.WG=10, D.WG=1)

a.var.mov <- function(DF.spp, m.bw.spp){
  DF.spp.df <- as.data.frame(DF.spp)
  if(class(DF.spp) != "SpatialPointsDataFrame") {
    coordinates(DF.spp) <- cbind(DF.spp[,grep("^lon$|^long$|^longitude$", colnames(DF.spp), ignore.case = T)],
                                 DF.spp[,grep("^lat$|^latitude$", colnames(DF.spp), ignore.case = T)])
  }
  DF.spp.df$row <- row.names(DF.spp.df)
  DF.spp.df <- merge(DF.spp.df, m.bw.spp, by.x = "species", by.y = "species", all.x = T)
  row.names(DF.spp.df) <- DF.spp.df$row
  # NW.df <- NW.df[order(as.numeric(NW.df$row)),]
  # DF.spp.df$TyrnnBR[is.na(DF.spp.df$TyrnnBR)] <- 0
  DF.spp <- SpatialPointsDataFrame(DF.spp, DF.spp.df)

  if(any(is.na(DF.spp@data[DF.spp@data$SeasPres==2,c("A.BG", "D.BG", "A.WG", "D.WG")]))){
    print("Species without movement information:")
    print(paste(unique(DF.spp@data$species[which(is.na(DF.spp@data[DF.spp@data$SeasPres==2,c("A.BG")]))])))
  }
  return(DF.spp)
}


