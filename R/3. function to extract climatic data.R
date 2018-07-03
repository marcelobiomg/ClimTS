


# function extract climate data by hemisphere and month/year of breeding record
a.var.mig8 <- function(r.Tmax.LTM = NULL, r.Tmin.LTM = NULL, r.Tmax.TS = NULL, r.Tmin.TS = NULL, r.P.LTM = NULL, r.P.TS = NULL, m.bw.sp, OBS, sp.shp, dist.w = F, max.dist.km=200){ #, DfP = NA
  { require(raster)
    require(sp)
    require(rgeos)

    crs(OBS) <- crs(sp.shp)
    var.names.T <- c("Tmax", "Tmin", "TyrRange", "Tisot", "Tseas.Hussell") #, "seas.Ricklefs")
    var.names.LTM <- paste0("LTM.", var.names.T)
    var.names.TS <- paste0("TS.", var.names.T)
    var.names.Anom <- paste0("Anom.", var.names.T)
    all.var.names.T <- c(var.names.LTM, var.names.TS, var.names.Anom)

    var.names.P <- c("Pmax", "Pmin", "PyrRange", "Pseas.Hussell") #, "seas.Ricklefs")
    var.names.P.LTM <- paste0("LTM.", var.names.P)
    var.names.P.TS <- paste0("TS.", var.names.P)
    var.names.P.Anom <- paste0("Anom.", var.names.P)
    all.var.names.P <- c(var.names.P.LTM, var.names.P.TS, var.names.P.Anom)

    # avoid warnings when only NAs are extracted
    robust.min <- function(x) {if(sum(!is.na(x))>0) min(x, na.rm=T) else NA}
    robust.max <- function(x) {if(sum(!is.na(x))>0) max(x, na.rm=T) else NA}

    # Temperature variables
    #    # ~BIO2 - mean month range
    month.range <- function(mx, mn) apply(mx - mn, 1, function(x) mean(x, na.rm=T) )
    #    # BIO4 = Temperature Seasonality (standard deviation *100)
    # yr.seas <- function(mx, mn) apply(c(mx, mn), 1, sd(x, na.rm=T))*100  #function(avg) apply(avg, 1, sd)*100
    ## computes seasonality - must be done with resources, not temperature
    # yrSeasR <- function(mx, mn) Tmax.yr(mx)/Tmin.yr(mn) #computes seasonality (ratio Ricklefs 1980) between x->min(winter) and y->max(summer) values across RasterBrick
    yrSeasH <- function(mx, mn) (Tmax.yr(mx)-Tmin.yr(mn))/(Tmin.yr(mn)) #computes seasonality (Hussell 1985) based on x->min(winter) and y->max(summer) values across RasterBrick
    #    # BIO5 - Max Temperature of Warmest Month
    Tmax.yr <- function(mx) apply(mx, 1, function(x) robust.max(x))
    #    # BIO6 - Min Temperature of Coldest Month
    Tmin.yr <- function(mn) apply(mn, 1, function(x) robust.min(x))
    #    # BIO7 = Temperature Annual Range (P5-P6) # calc range # the most extreme temperatures
    yr.range <- function(mx, mn) Tmax.yr(mx) - Tmin.yr(mn)
    #    # BIO3 = Isothermality (B2/B7) (* 100)
    yr.month.isotherm <- function(mx, mn) month.range(mx, mn)/yr.range(mx, mn)*100


    # Precipitation variables
    #    # BIO4 = Temperature Seasonality (standard deviation *100)
    ## computes seasonality - must be done with resources, not temperature
    # PyrSeasR <- function(mp) Pmax.yr(mp)/Pmin.yr(mp) #computes seasonality (ratio Ricklefs 1980) between x->min(winter) and y->max(summer) values across RasterBrick
    PyrSeasH <- function(mp) (Pmax.yr(mp)-Pmin.yr(mp))/(1+Pmin.yr(mp)) #computes seasonality (Hussell 1985) based on x->min(winter) and y->max(summer) values across RasterBrick
    #    # BIO5 - Max Temperature of Warmest Month
    Pmax.yr <- function(mp) apply(mp, 1, function(x) robust.max(x, na.rm=T) )
    #    # BIO6 - Min Temperature of Coldest Month
    Pmin.yr <- function(mp) apply(mp, 1, function(x) robust.min(x, na.rm=T) )
    #    # BIO7 = Temperature Annual Range (P5-P6) # calc range # the most extreme temperatures
    Pyr.range <- function(mp) Tmax.yr(mp) - Tmin.yr(mp)
  }

  #  Make spdf with unique values
  m.col.pos <- grep("^month$|?month$|^mes$", colnames(OBS@data), ignore.case = T) #, value=T
  m.col.val <- grep("^month$|?month$|^mes$", colnames(OBS@data), ignore.case = T, value=T)
  colnames(OBS@data)[m.col.pos] <- "month"

  OBS.unq.xy <- data.frame(unique(cbind(coordinates(OBS), OBS$SeasPres, OBS@data$month ))) # OBS$month
  colnames(OBS.unq.xy) <- c("lon", "lat","SeasPres", "month")
  coordinates(OBS.unq.xy) <- ~lon+lat
  crs(OBS.unq.xy) <- crs(OBS)
  OBS.unq.xy@data <- cbind(OBS.unq.xy@data, as.data.frame(setNames(replicate(length(c(all.var.names.P, all.var.names.T)), NA, simplify = F), c(all.var.names.P, all.var.names.T))))

  # cat(paste("SeasPres:", sort(unique(OBS.unq.xy$SeasPres))), "\n")
  cat(paste(sort(unique(OBS.unq.xy$SeasPres))), sep=" & ", fill = T, labels = c("-> SeasPres:"))

  #----#########----######## resident
  OBS.resid <- OBS.unq.xy$SeasPres == 1

  if(sum(OBS.resid)>0) {
    cat("Extracting values of breeding points - RESIDENT - LTM & TS", "\n")

    df.Tmax.LTM <- data.frame(extract(r.Tmax.LTM, OBS.unq.xy[OBS.resid,]))
    df.Tmax.LTM <- data.frame(matrix(NA, nrow(df.Tmax.LTM), ncol(df.Tmax.LTM)))
    coordinates(df.Tmax.LTM) <- coordinates(OBS.unq.xy[OBS.resid,])
    crs(df.Tmax.LTM) <- crs(OBS)
    df.Tmin.LTM <- df.Tmax.LTM
    df.Tmax.TS <- df.Tmax.LTM
    df.Tmin.TS <- df.Tmax.LTM
    df.P.LTM <- df.Tmax.LTM
    df.P.TS <- df.Tmax.LTM

    #----#########----######## resident
    #----######## LTM & TS - res - Breeding : 12 meses
    monthR <- unique(OBS.unq.xy$month[OBS.resid])
    for(mr in monthR){ ###### extrair registros por mês
      m.e.LTM <- ((mr+1):(mr+12))
      m.e.LTM <- ifelse(m.e.LTM > 12, m.e.LTM-12, m.e.LTM)
      m.e.TS <- ((mr+1):(mr+12))+12 # months of resident species to extract
      obs.m <- OBS.unq.xy$month[OBS.resid] == mr # obs to extract

      ## Temperatura
      #----######## Extraindo valores LTM ao longo do ano para todos os pontos # resident
      #  # Long Term Mean  - LTM
      { df.Tmax.LTM@data[obs.m,] <- extract(r.Tmax.LTM[[m.e.LTM]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])
        df.Tmin.LTM@data[obs.m,] <- extract(r.Tmin.LTM[[m.e.LTM]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])

        to.fill <- apply(df.Tmax.LTM@data, 1, function(x)any(is.na(x))) & obs.m
        if(any(to.fill) == T){
          print("filling NAs of resident temp LTM df.Tmax & df.Tmin")
          print(nlayers(r.Tmax.LTM))
          df.Tmax.LTM@data[to.fill,] <- ClosestAV(df.Tmax.LTM[to.fill,], r.Tmax.LTM[[m.e.LTM]], max.dist.km)@data
          df.Tmin.LTM@data[to.fill,] <- ClosestAV(df.Tmin.LTM[to.fill,], r.Tmin.LTM[[m.e.LTM]], max.dist.km)@data
        }
      }

      #----######## Extraindo valores TS ao longo do ano para todos os pontos # resident
      #  # Time Series  - TS
      { df.Tmax.TS@data[obs.m,] <- extract(r.Tmax.TS[[m.e.TS]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])
        df.Tmin.TS@data[obs.m,] <- extract(r.Tmin.TS[[m.e.TS]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])

        to.fill <- apply(df.Tmax.TS@data, 1, function(x)any(is.na(x))) & obs.m
        if(any(to.fill) ==T ){
          print("filling NAs of resident temp TS df.Tmax & df.Tmin")
          print(nlayers(r.Tmax.TS[[m.e.TS]]))
          df.Tmax.TS@data[to.fill,] <- ClosestAV(df.Tmax.TS[to.fill,], r.Tmax.TS[[m.e.TS]], max.dist.km)@data
          df.Tmin.TS@data[to.fill,] <- ClosestAV(df.Tmin.TS[to.fill,], r.Tmin.TS[[m.e.TS]], max.dist.km)@data
        } #
      }

      ## Precipitação
      #----######## Extraindo valores LTM ao longo do ano para todos os pontos # resident
      #  # Long Term Mean  - LTM
      { df.P.LTM@data[obs.m,] <- extract(r.P.LTM[[m.e.LTM]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])

        to.fill <- apply(df.P.LTM@data, 1, function(x)any(is.na(x))) & obs.m
        if(any(to.fill) == T){
          print("filling NAs of resident LTM df.Tmax & df.Tmin")
          print(nlayers(r.P.LTM))
          df.P.LTM@data[to.fill,] <- ClosestAV(df.P.LTM[to.fill,], r.P.LTM[[m.e.LTM]], max.dist.km)@data
        }
      }

      #----######## Extraindo valores TS ao longo do ano para todos os pontos # resident
      #  # Time Series  - TS
      { df.P.TS@data[obs.m,] <- extract(r.P.TS[[m.e.TS]], OBS.unq.xy[OBS.resid & OBS.unq.xy$month==mr,])

        to.fill <- apply(df.P.TS@data, 1, function(x)any(is.na(x))) & obs.m
        if(any(to.fill) ==T ){
          print("filling NAs of resident TS df.Tmax & df.Tmin")
          print(nlayers(r.P.TS[[m.e.TS]]))
          df.P.TS@data[to.fill,] <- ClosestAV(df.P.TS[to.fill,], r.P.TS[[m.e.TS]], max.dist.km)@data
        } #
      }

      if(mr == monthR[length(monthR)]){
        #### Variação de TEMPERATURA experimentada pela espécie
        #### year round # breeding ground only
        {# LTM
          OBS.unq.xy@data[OBS.resid, var.names.LTM[1]] <- Tmax.yr(df.Tmax.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.LTM[2]] <- Tmin.yr(df.Tmin.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.LTM[3]] <- yr.range(df.Tmax.LTM@data, df.Tmin.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.LTM[4]] <- yr.month.isotherm(df.Tmax.LTM@data, df.Tmin.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.LTM[5]] <- yrSeasH(df.Tmax.LTM@data, df.Tmin.LTM@data)
          # OBS.unq.xy@data[OBS.resid,var.names.LTM[6]] <- yrSeasR(df.Tmax.LTM@data, df.Tmin.LTM@data)

          # TS
          TS.max <- df.Tmax.TS@data + df.Tmax.LTM@data
          TS.min <- df.Tmin.TS@data + df.Tmin.LTM@data
          OBS.unq.xy@data[OBS.resid, var.names.TS[1]] <- Tmax.yr(TS.max)
          OBS.unq.xy@data[OBS.resid, var.names.TS[2]] <- Tmin.yr(TS.min)
          OBS.unq.xy@data[OBS.resid, var.names.TS[3]] <- yr.range(TS.max, TS.min)
          OBS.unq.xy@data[OBS.resid, var.names.TS[4]] <- yr.month.isotherm(TS.max, TS.min)
          OBS.unq.xy@data[OBS.resid, var.names.TS[5]] <- yrSeasH(TS.max, TS.min)
          # OBS.unq.xy@data[, var.names.TS[5]] <- yrSeasR(TS.max, TS.min)

          # Anom
          OBS.unq.xy@data[OBS.resid, var.names.Anom[1]] <-  OBS.unq.xy@data[OBS.resid, var.names.TS[1]] - OBS.unq.xy@data[OBS.resid, var.names.LTM[1]]
          OBS.unq.xy@data[OBS.resid, var.names.Anom[2]] <-  OBS.unq.xy@data[OBS.resid, var.names.TS[2]] - OBS.unq.xy@data[OBS.resid, var.names.LTM[2]]
          OBS.unq.xy@data[OBS.resid, var.names.Anom[3]] <-  OBS.unq.xy@data[OBS.resid, var.names.TS[3]] - OBS.unq.xy@data[OBS.resid, var.names.LTM[3]]
          OBS.unq.xy@data[OBS.resid, var.names.Anom[4]] <-  OBS.unq.xy@data[OBS.resid, var.names.TS[4]] - OBS.unq.xy@data[OBS.resid, var.names.LTM[4]]
          OBS.unq.xy@data[OBS.resid, var.names.Anom[5]] <-  OBS.unq.xy@data[OBS.resid, var.names.TS[5]] - OBS.unq.xy@data[OBS.resid, var.names.LTM[5]]
          # OBS.unq.xy@data[,var.names.Anom[6]] <-  OBS.unq.xy@data[,var.names.TS[5]] - OBS.unq.xy@data[,var.names.LTM[5]]
        }


        #### Variação de PRECIPITAÇÃO experimentada pela espécie
        #### year round # breeding ground only
        {# LTM
          OBS.unq.xy@data[OBS.resid, var.names.P.LTM[1]] <- Pmax.yr(df.P.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.P.LTM[2]] <- Pmin.yr(df.P.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.P.LTM[3]] <- Pyr.range(df.P.LTM@data)
          OBS.unq.xy@data[OBS.resid, var.names.P.LTM[4]] <- PyrSeasH(df.P.LTM@data)
          # OBS.unq.xy@data[OBS.resid,var.names.LTM[6]] <- yrSeasR(df.Tmax.LTM@data, df.Tmin.LTM@data)

          # CHECAR!!!        # TS
          TS.P <- df.P.TS@data # + df.P.LTM@data
          # TS.Pmin <- df.P.TS@data + df.P.LTM@data
          OBS.unq.xy@data[OBS.resid, var.names.P.TS[1]] <- Pmax.yr(TS.P)
          OBS.unq.xy@data[OBS.resid, var.names.P.TS[2]] <- Pmin.yr(TS.P)
          OBS.unq.xy@data[OBS.resid, var.names.P.TS[3]] <- Pyr.range(TS.P)
          OBS.unq.xy@data[OBS.resid, var.names.P.TS[4]] <- PyrSeasH(TS.P)
          # OBS.unq.xy@data[, var.names.TS[5]] <- yrSeasR(TS.max, TS.min)

          # Anom
          OBS.unq.xy@data[OBS.resid, var.names.P.Anom[1]] <-  OBS.unq.xy@data[OBS.resid, var.names.P.TS[1]] - OBS.unq.xy@data[OBS.resid, var.names.P.LTM[1]]
          OBS.unq.xy@data[OBS.resid, var.names.P.Anom[2]] <-  OBS.unq.xy@data[OBS.resid, var.names.P.TS[2]] - OBS.unq.xy@data[OBS.resid, var.names.P.LTM[2]]
          OBS.unq.xy@data[OBS.resid, var.names.P.Anom[3]] <-  OBS.unq.xy@data[OBS.resid, var.names.P.TS[3]] - OBS.unq.xy@data[OBS.resid, var.names.P.LTM[3]]
          OBS.unq.xy@data[OBS.resid, var.names.P.Anom[4]] <-  OBS.unq.xy@data[OBS.resid, var.names.P.TS[4]] - OBS.unq.xy@data[OBS.resid, var.names.P.LTM[4]]
          # OBS.unq.xy@data[OBS.resid, var.names.P.Anom[5]] <-  OBS.unq.xy@data[OBS.resid, var.names.P.TS[5]] - OBS.unq.xy@data[OBS.resid, var.names.P.LTM[5]]
          # OBS.unq.xy@data[,var.names.P.Anom[6]] <-  OBS.unq.xy@data[,var.names.TS[5]] - OBS.unq.xy@data[,var.names.LTM[5]]
        }
      }
    }

  }

  #----#########----######## migratory
  OBS.mig.NS <- OBS.unq.xy$SeasPres != 1
  if(sum(OBS.mig.NS)>0) {
    #####--# Breeding and wintering seasons species specific
    # criar o vetor corretamente para os meses de invernada, que pegam os últimos meses antes do solstício
    # e primeiros meses após o solstício
    if(m.bw.sp$A.WG > m.bw.sp$D.WG){ # A - arrival, D - departure, WG - wintering ground, BG - breeding ground
      m.w.slstc <- m.bw.sp$A.WG:(m.bw.sp$D.WG+12)
      # m.w <- ifelse(m.w > 12, m.w-12, m.w)  # excluded because of taking months of previous year
      # m.b <- m.bw.sp$A.BG:m.bw.sp$D.BG      # excluded because of taking months of previous year
      m.b.slstc <- m.bw.sp$A.BG:m.bw.sp$D.BG   # included because of taking months of previous year
    } else {
      # m.w <- m.bw.sp$A.WG:m.bw.sp$D.WG      # excluded because of taking months of previous year
      m.w.slstc <- m.bw.sp$A.WG:m.bw.sp$D.WG   # included because of taking months of previous year
      # m.b <- m.bw.sp$A.BG:(m.bw.sp$D.BG+12)
      # m.b <- ifelse(m.b > 12, m.b-12, m.b)
      m.b.slstc <- m.bw.sp$A.BG:m.bw.sp$D.BG
    }


    #----######## separando obs por hemisf
    OBS.unq.S <- coordinates(OBS.unq.xy)[, 2] < 0 & OBS.mig.NS
    OBS.unq.N <- coordinates(OBS.unq.xy)[, 2] >= 0 & OBS.mig.NS

    # determinar qual(is) hemisfério(s) para o for()
    if(length(OBS.unq.xy[OBS.unq.S,]) > 0 & length(OBS.unq.xy[OBS.unq.N,]) == 0) {
      hems <- "S"
    } else if(length(OBS.unq.xy[OBS.unq.S,]) == 0 & length(OBS.unq.xy[OBS.unq.N,]) > 0) {
      hems <- "N"
    } else {
      hems <- c("N","S")
    }


    ###### Cálculos variáveis por hemisfério
    # print(hems)
    for(h in hems) {
      cat(paste("->", h, "Hemisphere"), "\n")
      if(h == "N"){
        OBS.mig.h <- OBS.unq.N
      } else {
        OBS.mig.h <- OBS.unq.S
      }

      ###### extrair registros de cada mês
      monthM <- unique(OBS.unq.xy$month[OBS.mig.h])
      for(m in monthM){
        OBS.mig.hm <- OBS.mig.h & OBS.unq.xy$month==m

        #----######## Extraindo valores nas áreas de reprodução e invernada # migratory
        ## só calcula para espécies migratórias
        if(length(OBS.mig.hm)>0 & length(sp.shp[sp.shp$SEASONAL==3,])>0){ ## só calcula para espécies migratórias

          # se for Hemisf. N, o m.w são os meses reais # se for Hemisf. S tem que modificar
          if(h == "N"){
            m.b <- ((m+11):(m+12)) + 12
            m.w <- m.w.slstc + 12
            m.w.LTM <- ifelse(m.w.slstc > 12, m.w.slstc-12, m.w.slstc)
          } else {
            m.b <- ((m+11):(m+12)) + 12 # m.b [month breeding] é o próprio valor, pq é o mês do registro e não o mês após o solsctício, como o m.w
            m.w <- m.w.slstc + 18
            m.w.LTM <- m.w.slstc - 6
          }
          m.w <- ifelse(m.w > m.b[1], m.w <- m.w-12, m.w)
          m.b.LTM <- (m-1):m
          m.b.LTM <- ifelse(m.b.LTM < 1, m.b.LTM + 12, m.b.LTM)


          #----######## Extraindo valores para a área de invernada
          #----######## LTM TS - mig - W
          cat("Extracting values of wintering range - MIGRATORY - LTM & TS", "\n")
          # sp.shp$SEASONAL
          # 1 = year-round areas (resident)
          # 2 = breeding areas (summer)
          # 3 = wintering areas (winter)
          if(dist.w == F) {
            # 1. extract mean values from the whole wintering range
            cat(paste("  Using mean values from the whole wintering range:", "dist.w == F"), "\n")
            if(m == monthM[1]){
              Tmax.w.LTM <- extract(r.Tmax.LTM[[m.w.LTM]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)
              Tmin.w.LTM <- extract(r.Tmin.LTM[[m.w.LTM]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)
              P.w.LTM <- extract(r.P.LTM[[m.w.LTM]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)
            }
            Tmax.w.TS <- extract(r.Tmax.TS[[m.w]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)
            Tmin.w.TS <- extract(r.Tmin.TS[[m.w]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)
            P.w.TS <- extract(r.P.TS[[m.w]], gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), fun=mean, na.rm=TRUE)

          } else {
            # 2. use function to extract weighted values
            cat(paste("  Using values from wintering range weighted by distance:", "dist.w == T"), "\n")
            # create spp seasonal presence based on polygon
            sp.m.p.LTM <- rasterize(gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), r.Tmax.LTM)
            sp.m.p.TS <- rasterize(gUnaryUnion(sp.shp[sp.shp$SEASONAL==3,]), r.Tmax.TS)

            # function to extract weighted values
            if(m == monthM[1]){
              Tmax.w.LTM <- f.wr(r.Tmax.LTM, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w.LTM)
              Tmin.w.LTM <- f.wr(r.Tmin.LTM, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w.LTM)
              P.w.LTM <- f.wr(r.P.LTM, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w.LTM)
            }
            Tmax.w.TS <- f.wr(r.Tmax.TS, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w)
            Tmin.w.TS <- f.wr(r.Tmin.TS, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w)
            P.w.TS <- f.wr(r.P.TS, xy = OBS.mig.hm, sp.m.p = sp.m.p.LTM, m.w = m.w)
          }



          #----######## Extraindo valores para pontos da área reprodutiva
          #----######## LTM TS - mig - B
          cat("Extracting values of breeding points - MIGRATORY - LTM & TS", "\n")
          { Tmax.b.LTM <- as.data.frame(extract(r.Tmax.LTM[[m.b.LTM]], OBS.unq.xy[OBS.mig.hm,]))
            Tmin.b.LTM <- as.data.frame(extract(r.Tmin.LTM[[m.b.LTM]], OBS.unq.xy[OBS.mig.hm,]))
            P.b.LTM <- as.data.frame(extract(r.P.LTM[[m.b.LTM]], OBS.unq.xy[OBS.mig.hm,]))
            Tmax.b.TS <- as.data.frame(extract(r.Tmax.TS[[m.b]], OBS.unq.xy[OBS.mig.hm,]))
            Tmin.b.TS <- as.data.frame(extract(r.Tmin.TS[[m.b]], OBS.unq.xy[OBS.mig.hm,]))
            P.b.TS <- as.data.frame(extract(r.P.TS[[m.b]], OBS.unq.xy[OBS.mig.hm,]))

            {coordinates(Tmax.b.LTM) <- coordinates(OBS.unq.xy[OBS.mig.hm,])
              coordinates(Tmin.b.LTM) <- coordinates(Tmax.b.LTM)
              coordinates(P.b.LTM) <- coordinates(Tmax.b.LTM)
              coordinates(Tmax.b.TS) <- coordinates(Tmax.b.LTM)
              coordinates(Tmin.b.TS) <- coordinates(Tmax.b.LTM)
              coordinates(P.b.TS) <- coordinates(Tmax.b.LTM)

              crs(Tmax.b.LTM) <- crs(OBS)
              crs(Tmin.b.LTM) <- crs(OBS)
              crs(P.b.LTM) <- crs(OBS)
              crs(Tmax.b.TS) <- crs(OBS)
              crs(Tmin.b.TS) <- crs(OBS)
              crs(P.b.TS) <- crs(OBS)}}

          # checking for NAs in the extraction DF
          # TEMPERATURE
          if(any(is.na(Tmax.b.LTM@data)) ==T ){
            print("filling NAs of migratory LTM df.Tmax & df.Tmin")
            print(nlayers(r.Tmax.LTM[[m.b.LTM]]))
            Tmax.b.LTM <- ClosestAV(Tmax.b.LTM, r.Tmax.LTM[[m.b.LTM]], max.dist.km)
            Tmin.b.LTM <- ClosestAV(Tmin.b.LTM, r.Tmin.LTM[[m.b.LTM]], max.dist.km)
          }
          if(any(is.na(Tmax.b.TS@data)) ==T ){
            print("filling NAs of migratory TS df.Tmax & df.Tmin")
            print(nlayers(r.Tmax.TS[[m.b]]))
            Tmax.b.TS <- ClosestAV(Tmax.b.TS, r.Tmax.TS[[m.b]], max.dist.km)
            Tmin.b.TS <- ClosestAV(Tmin.b.TS, r.Tmin.TS[[m.b]], max.dist.km)
          }

          # checking for NAs in the extraction DF
          # PRECIPITATION
          if(any(is.na(P.b.LTM@data)) ==T ){
            print("filling NAs of migratory LTM df.Tmax & df.Tmin")
            print(nlayers(r.P.LTM[[m.b.LTM]]))
            P.b.LTM <- ClosestAV(P.b.LTM, r.P.LTM[[m.b.LTM]], max.dist.km)
          }
          if(any(is.na(P.b.TS@data)) ==T ){
            print("filling NAs of migratory TS df.Tmax & df.Tmin")
            print(nlayers(r.P.TS[[m.b]]))
            P.b.TS <- ClosestAV(P.b.TS, r.P.TS[[m.b]], max.dist.km)
          }


          # joining W and B extractions
          df.Tmax.LTM <- data.frame(Tmax.b.LTM@data, Tmax.w.LTM)
          df.Tmin.LTM <- data.frame(Tmin.b.LTM@data, Tmin.w.LTM)
          df.P.LTM <- data.frame(P.b.LTM@data, P.w.LTM)
          df.Tmax.TS <- data.frame(Tmax.b.TS@data, Tmax.w.TS)
          df.Tmin.TS <- data.frame(Tmin.b.TS@data, Tmin.w.TS)
          df.P.TS <- data.frame(P.b.TS@data, P.w.TS)
          # }

          #### Variação de TEMPERATURA experimentada pela espécie
          #### year round # breeding ground only
          # LTM
          OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[1]] <- Tmax.yr(df.Tmax.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[2]] <- Tmin.yr(df.Tmin.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[3]] <- yr.range(df.Tmax.LTM, df.Tmin.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[4]] <- yr.month.isotherm(df.Tmax.LTM, df.Tmin.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[5]] <- yrSeasH(df.Tmax.LTM, df.Tmin.LTM)
          # OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[6]] <- yrSeasR(df.Tmax.LTM@data, df.Tmin.LTM@data)
          # TS
          TS.max <- df.Tmax.TS + df.Tmax.LTM
          TS.min <- df.Tmin.TS + df.Tmin.LTM
          OBS.unq.xy@data[OBS.mig.hm, var.names.TS[1]] <- Tmax.yr(TS.max)
          OBS.unq.xy@data[OBS.mig.hm, var.names.TS[2]] <- Tmin.yr(TS.min)
          OBS.unq.xy@data[OBS.mig.hm, var.names.TS[3]] <- yr.range(TS.max, TS.min)
          OBS.unq.xy@data[OBS.mig.hm, var.names.TS[4]] <- yr.month.isotherm(TS.max, TS.min)
          OBS.unq.xy@data[OBS.mig.hm, var.names.TS[5]] <- yrSeasH(TS.max, TS.min)
          # OBS.unq.xy@data[OBS.mig.hm, var.names.TS[6]] <- yrSeasR(TS.max, TS.min)
          # Anom
          OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[1]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[1]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[1]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[2]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[2]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[2]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[3]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[3]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[3]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[4]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[4]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[4]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[5]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[5]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[5]]
          # OBS.unq.xy@data[OBS.mig.hm, var.names.Anom[6]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.TS[6]] - OBS.unq.xy@data[OBS.mig.hm, var.names.LTM[6]]


          #### Variação de PRECIPITAÇÃO experimentada pela espécie
          #### year round # breeding ground only
          # LTM
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[1]] <- Pmax.yr(df.P.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[2]] <- Pmin.yr(df.P.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[3]] <- Pyr.range(df.P.LTM)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[4]] <- PyrSeasH(df.P.LTM)
          # OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[5]] <- yrSeasR(df.P.LTM)

          # CHECAR!!!          # TS
          TS.P <- df.P.TS #+ df.P.LTM
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[1]] <- Pmax.yr(TS.P)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[2]] <- Pmin.yr(TS.P)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[3]] <- Pyr.range(TS.P)
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[4]] <- PyrSeasH(TS.P)
          # OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[5]] <- PyrSeasR(TS.P)

          # Anom
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.Anom[1]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[1]] - OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[1]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.Anom[2]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[2]] - OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[2]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.Anom[3]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[3]] - OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[3]]
          OBS.unq.xy@data[OBS.mig.hm, var.names.P.Anom[4]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[4]] - OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[4]]
          # OBS.unq.xy@data[OBS.mig.hm, var.names.P.Anom[5]] <-  OBS.unq.xy@data[OBS.mig.hm, var.names.P.TS[5]] - OBS.unq.xy@data[OBS.mig.hm, var.names.P.LTM[5]]

          # if(length(hems) == 1) {OBS.unq.xy <- OBS.unq.N}
        } # if OBS.mig

      } # end for monthM
    } # end for hems
  } # end for OBS.mig.NS

  OBS.unq.xy@data$lat <- coordinates(OBS.unq.xy)[,2]
  OBS.unq.xy@data$lon <- coordinates(OBS.unq.xy)[,1]

  OBS <- merge(OBS@data, OBS.unq.xy@data, by.x= c("SeasPres", "lat", "lon", "month"), by.y=c("SeasPres", "lat", "lon", "month"),
               all.x=T)[, union(names(OBS@data), names(OBS.unq.xy@data))]
  colnames(OBS)[grep("^month$|?month$|^mes$", colnames(OBS), ignore.case = T)] <- m.col.val
  return(OBS)
}


# for several species
a.var.mig.batch7 <- function(r.Tmax.LTM = NULL, r.Tmin.LTM = NULL, r.Tmax.TS = NULL, r.Tmin.TS = NULL, r.P.LTM = NULL, r.P.TS = NULL, m.bw.spp = NULL, DF, sp.shp.path, dist.w = F, max.dist.km=150){ #, DfP = NA
  library(raster)
  spp.list <- sort(unique(DF$species))
  spp.mig <- sort(unique(DF$species[DF$SeasPres==2]))
  if(is.null(m.bw.spp)){
    m.bw.spp <- unique(data.frame(species=DF$species, A.BG=DF$A.BG, D.BG=DF$D.BG, A.WG=DF$A.WG, D.WG=DF$D.WG))
  }

  # 1.2 criar lista de espécies da pasta
  spp.folder <- list.files(sp.shp.path, pattern = "\\.shp$", full.names = FALSE)
  spp.folder.F <- list.files(sp.shp.path, pattern = "\\.shp$", full.names = T)

  yr <- sort(unique(as.numeric(paste(DF@data[,grep("^year$|?year$|^ano$", colnames(DF@data), ignore.case = T)]))))
  for(s in 1:length(spp.list)){
    cat("\n")
    cat(yr, paste(spp.list[s]), fill=1, labels = c("---------------- Year:","---------------- Species:"))
    # cat(yr, paste(spp.list[s]), fill=paste(T, T), labels = paste( "Year:","Species"))

    # print(paste0("Year: ", yr, " - Species: ", spp.list[s]))
    # OBS <- subset(DF, DF$species == spp.list[s])
    OBS <- DF[DF$species == spp.list[s],]

    # 1.4.1 cada espécie migratória na pasta e..
    sp.shp <- shapefile(spp.folder.F[grep(sub(" ", "_", spp.list[s]), spp.folder.F)])

    ## TODO - make a loop for each breeding season area
    #- separate OBS by area of m.bw.sp
    m.bw.sp <- m.bw.spp[m.bw.spp$species == paste(spp.list[s]),]


    OBSi <- a.var.mig8(r.Tmax.LTM=r.Tmax.LTM, r.Tmin.LTM=r.Tmin.LTM,
                       r.Tmax.TS=r.Tmax.TS, r.Tmin.TS=r.Tmin.TS,
                       r.P.LTM = r.P.LTM, r.P.TS = r.P.TS,
                       m.bw.sp=m.bw.sp, OBS=OBS, sp.shp=sp.shp, dist.w=dist.w, max.dist.km=max.dist.km) #, DfP=DfP

    if(s==1){OBS.batch <- OBSi} else {OBS.batch <- rbind(OBS.batch, OBSi)}
  }
  return(OBS.batch)
}




######
a.var.mig.TS2 <- function(r.Tmax.LTM = NULL, r.Tmin.LTM = NULL, r.Tmax.TS = NULL, r.Tmin.TS = NULL,
                          r.P.LTM = NULL, r.P.TS = NULL, m.bw.spp = NULL, DF, sp.shp.path,
                          dist.w = F, max.dist.km=150) {
  require(raster)
  require(sp)
  require(rgeos)
  # select all observations from resident? probably not

  #-------- mudar para year depois
  yTS <- sort(unique(as.numeric(paste(DF@data[,grep("^year$|?year$|^ano$", colnames(DF@data), ignore.case = T)]))))
  cat("================================\n", paste(length(yTS)), "year(s) to extract:\n", paste(yTS, collapse = ", "), "\n")

  # finding first and last years of each time-series raster stack
  lastYrTS <- trunc(max(as.numeric(unlist(r.Tmax.TS@z)))) # Temperature
  frstYrTS <- trunc(min(as.numeric(unlist(r.Tmax.TS@z))))
  lastYrPTS <- trunc(max(as.numeric(format(r.P.TS@z$time, format="%Y")))) # Precipitation
  frstYrPTS <- trunc(min(as.numeric(format(r.P.TS@z$time, format="%Y"))))

  # list of species
  spp.list <- sort(unique(DF$species))
  spp.mig <- sort(unique(DF$species[DF$SeasPres==2 ]))

  # criar lista de espécies com dados de movimentação
  # 1.2 criar lista de espécies da pasta
  spp.folder <- list.files(sp.shp.path, pattern = "\\.shp$", full.names = FALSE)

  ### corrigir! colocar opção dos dados na tabela
  if(is.null(m.bw.spp)){
    m.bw.spp <- unique(data.frame(species=DF$species, A.BG=DF$A.BG, D.BG=DF$D.BG, A.WG=DF$A.WG, D.WG=DF$D.WG))
  }
  # testar se todas as espécies do data.frame possuem dados de movimentação e estão nas pastas dos shapefiles
  if(length(spp.mig) > length(unique(m.bw.spp$species))) paste("Seasonal movement dates are lacking for", length(spp.mig) - length(m.bw.spp$species) ,"species")
  if(length(spp.mig) > length(pmatch(sub(" ", "_", spp.mig), spp.folder))) paste("Shapefiles are lacking for", length(spp.mig) - length(m.bw.spp$species) ,"species")

  # year loops
  for(y in yTS){
    # print("Year:")
    # print(y)
    # cat("================================", "", yr, sep = "",  fill = 1, labels = c(""," Year(s) to extract:", ""))
    # select all obs from a single year
    #-------- mudar para year depois
    # OBSy <- DF[DF$Year.Publ == y,] # DF$year
    OBSy <- DF[DF@data[,grep("^year$|?year$|^ano$", colnames(DF@data), ignore.case = T)] == y,]

    # select months from a single year
    # check if year is before (< than) first year in raster brick
    {
    # Temperature
    if(y >= frstYrTS){
      y_0 <- as.numeric(grep(y, names(r.Tmax.TS)))
    }else{
      y_0 <- 1:12
      print(paste("BE CAUTIOUS! No climatic data for year:", y, "Using first layers, of year:", frstYrTS))}
    if(y-1 >= frstYrTS) {
      y_1 <- as.numeric(grep(y-1, names(r.Tmax.TS)))
    }else{y_1 <- y_0}
    if(y-2 >= frstYrTS) {
      y_2 <- as.numeric(grep(y-2, names(r.Tmax.TS)))
    }else{y_2 <- y_1}

    # Precipitation
    if(y >= frstYrPTS){
      yP_0 <- as.numeric(grep(y, gsub("X", "", names(r.P.TS))))
    }else{
      yP_0 <- 1:12
      print(paste("BE CAUTIOUS! No climatic data for year:", y, "Using first layers, of year:", frstYrPTS))}
    if(y-1 >= frstYrPTS) {
      yP_1 <- as.numeric(grep(y-1, gsub("X", "", names(r.P.TS))))
    }else{yP_1 <- yP_0}
    if(y-2 >= frstYrPTS) {
      yP_2 <- as.numeric(grep(y-2, gsub("X", "", names(r.P.TS))))
    }else{yP_2 <- yP_1}
}

    # check if year is after (> than) last year in raster brick
    {
    # Temperature
    if(y > lastYrTS){
      y_0 <- as.numeric(grep(lastYrTS, names(r.Tmax.TS)))
      y_1 <- y_0 #as.numeric(grep(lastYrTS, names(r.Tmax.TS)))
      y_2 <- y_0 #as.numeric(grep(lastYrTS, names(r.Tmax.TS)))
      print(paste("BE CAUTIOUS! No temperature data for year:", y, "Using last layers, of year:", lastYrTS))
    }
    months.year.y <- c(y_2, y_1, y_0)
    # months.year.y <- c(as.numeric(grep(y-1, names(r.Tmax.TS))), as.numeric(grep(y-1, names(r.Tmax.TS))), as.numeric(grep(y, names(r.Tmax.TS))) )#,

    # Precipitation
    if(y > lastYrPTS){
      yP_0 <- as.numeric(grep(lastYrPTS, names(r.P.TS)))
      yP_1 <- yP_0 #as.numeric(grep(lastYrTS, names(r.Tmax.TS)))
      yP_2 <- yP_0 #as.numeric(grep(lastYrTS, names(r.Tmax.TS)))
      print(paste("BE CAUTIOUS! No precipitation data for year:", y, "Using last layers, of year:", lastYrPTS))
    }
    months.yearP.y <- c(yP_2, yP_1, yP_0)
    }

    r.Tmax.TSy <- r.Tmax.TS[[months.year.y]]
    r.Tmin.TSy <- r.Tmin.TS[[months.year.y]]
    r.P.TSy <- r.P.TS[[months.yearP.y]]

    # (r.Tmax.LTM,           r.Tmin.LTM,           r.Tmax.TS,          r.Tmin.TS,          m.bw.spp = NULL,   DF,      sp.shp.path,             dist.w = F){ #, DfP = NA
    OBSy <- a.var.mig.batch7(r.Tmax.LTM=r.Tmax.LTM, r.Tmin.LTM=r.Tmin.LTM,
                             r.Tmax.TS=r.Tmax.TSy, r.Tmin.TS=r.Tmin.TSy,
                             r.P.LTM = r.P.LTM, r.P.TS=r.P.TSy,
                             m.bw.spp=m.bw.spp, DF=OBSy, sp.shp.path=sp.shp.path, dist.w=dist.w)

    if(y == yTS[1]) {OBS.final <- OBSy} else {OBS.final <- rbind(OBS.final, OBSy)}
  }

  return(OBS.final)
}

system.time(EggsT.Seas.Clim.test <- a.var.mig.TS2(r.Tmax.LTM = LTM.Tmax.NW, r.Tmin.LTM = LTM.Tmin.NW,
                                                  r.Tmax.TS = BEST.Tmax.NW, r.Tmin.TS = BEST.Tmin.NW,
                                                  r.P.LTM = GPCC.prec.LTM.NW, r.P.TS = GPCC.prec.m.NW.1850.2013,
                                                  m.bw.spp=NULL, DF=EggsT.Seas.test2, sp.shp.path=sp.shp.path,
                                                  dist.w = F, max.dist.km=200))

system.time(EggsT.Seas.Clim <- a.var.mig.TS2(r.Tmax.LTM = LTM.Tmax.NW, r.Tmin.LTM = LTM.Tmin.NW,
                                             r.Tmax.TS = BEST.Tmax.NW, r.Tmin.TS = BEST.Tmin.NW,
                                             r.P.LTM = GPCC.prec.LTM.NW, r.P.TS = GPCC.prec.m.NW.1850.2013,
                                             m.bw.spp=NULL, DF=EggsT.Seas, sp.shp.path=sp.shp.path,
                                             dist.w = F, max.dist.km=150.5))


#########

### drop observations with not enough number of observations/(spp & Köppen−Geiger climate)
## make a function to subset the data
select <- function(x, y, v.x, v.y, df) df[df[[v.x]] == x & df[[v.y]] == y, ]

EggsT.Seas.Clim.sel.TP <- mapply(select, x=EggsT.spp.KGclim.S.TP$species, y=EggsT.spp.KGclim.S.TP$KGclimate.ID,
                                 MoreArgs=list(v.x="species", v.y="KGclimate.ID", df=EggsT.Seas.Clim.TP), SIMPLIFY=F)



#####

######## clear workspace
rm(m.b, m.b.LTM, m.b.slstc, m.bw.sp, m.bw.spi, m.bw.spp, seas.pres.shp, xy, col, cols, DF, DfP, h,
   hems, i, m, m.e, migratory, m.w, m.w.LTM, m.w.slstc, month, monthM, monthR, months.year.y, mr, OBS, obs.m,
   OBS.mig.h, OBS.mig.hm, OBS.mig.NS, OBS.resid, OBS.spi, OBS.unq.N, OBS.unq.S, OBS.unq.xy, OBSy,
   r.Tmax.LTM.t, r.Tmax.TS, r.Tmax.TSy,r.Tmin.TS, r.Tmin.TSy, remov, s, sp, spp.folder, spp.folder.F,
   spp.found, spp.list, spp.mig, spp.rec, var.names, var.names.Anom, var.names.LTM, var.names.TS,
   x, x.min, x.max, y, y.min, y.max, yTS, y_0, y_1, y_2, yP_0, yP_1, yP_2, spp.SP, TS.max, NW.df, Tmax.w.LTM, Tmax.w.TS,
   df.Tmax.LTM, df.Tmax.LTM2, df.Tmax.TS, df.Tmin.LTM, df.Tmin.TS, frstYrTS, lastYrTS, m.col, m.col.pos,
   m.col.val, max.dist.km, path.shp.Cot, path.shp.Pipr, path.shp.Tyr, Tmax.yr, Tmin.yr, yr.month.isotherm,
   Tmin.b.LTM, Tmin.b.TS, Tmin.w.LTM, Tmin.w.TS, TS.min, extrct, m.e.LTM, m.e.TS, r, spdf, to.fill,
   Tmax.b.LTM, Tmax.b.TS, Tmin.b.LTM, Tmin.b.TS, yr, validYrMth, OBSi, OBS.na, TS.P, all.var.names.P, all.var.names.T,
   df.P.LTM, df.P.TS, EggsT.Seas.test2, frstYrPTS, lastYrPTS, months.yearP.y, r.P.LTM, r.P.TSm, r.P.TS.1900,
   r.P.TSy, r.Tmax.LTM, r.Tmin.LTM, sp.shp, var.names.P, var.names.P.Anom, var.names.P.LTM, var.names.P.TS,
   var.names.T, Pmax.yr, Pmin.yr, Pyr.range, PyrSeasH, yr.range, yrSeasH, yrSeasR)
m.e.TS
rm(r.Tmax.LTM, r.Tmin.LTM, r.Tmax.TS, r.Tmin.TS, DF, dist.w, max.dist.km)
