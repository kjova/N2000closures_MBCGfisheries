#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
##  Rscript for the analyse underlying the article 
#      "Evaluating the efficacy of closing the current European network of marine protected areas 
#       to mobile bottom-contacting fishing gears"
#
##  Step 1: Processing of (publicly available) fisheries data
#
##  Code by Karin van der Reijden. (kjvdreijden@gmail.com)
#        published under the MIT-license
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Step 0. Load libraries, set data paths, load ICES Ecoregions ----
library(sf)
library(data.table)
library(terra)
library(RColorBrewer)
library(basemaps)
library(csquares)
library(grDevices)
sf_use_s2(FALSE)

##-- set data paths
datPath <- "D:/MS_Natura2000closures/Data/"
outPath <- "D:/MS_Natura2000closures/Routput/"

##-- Load ICES Ecoregions, and clip to the ones of interest
# Downloadable here: https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb
ICESEcoregs <- st_read(paste0(datPath, "ICES_ecoregions_20171207_erase_ESRI.shp"))
ICESEcoregs <- subset(ICESEcoregs, Ecoregion %in% c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast"))
ICESEcoregs <- ICESEcoregs[,c("Ecoregion")]

#-#-# Step 1. OSPAR-fisheries data ----
FishDat <- st_read(paste0(datPath, "ICES.2021.OSPAR_production_of_spatial_fishing_pressure_data_layers/simple_features/benthic_metiers.csv"))
FishDat <- data.table(st_drop_geometry(FishDat))
FishDat <- FishDat[,c(1:5, 7, 16:17)]
NumCols <- c("Year", "lat", "lon", "surface", "TotWt_low", "TotWt_upp")
FishDat <- FishDat[, (NumCols) := lapply(.SD, as.numeric), .SDcols = NumCols]
rm(NumCols); gc()

##-- Clip to ICES Ecoregions of interest -- only GNS, CS, BoB&IC
FD <- FishDat[, c("C.square", "lat", "lon")]
FD <- FD[!duplicated(FD),]
FD2 <- st_as_sf(FD, coords = c("lon","lat"))
st_crs(FD2) <- st_crs(ICESEcoregs)
FD2$lon = FD$lon
FD2$lat =FD$lat
OSPARarea <- subset(ICESEcoregs, Ecoregion %in% c("Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast"))
FD2 <- subset(FD2, lon >= -16.1 & lon <= 13.1)
FD2 <- subset(FD2, lat <= 63.9)
FD3 <- st_intersection (FD2, OSPARarea)
FishDat <- subset(FishDat, C.square %in% FD3$C.square)
FishDat$Ecoregion = FD3$Ecoregion [match(FishDat$C.square, FD3$C.square)]
rm(FD, FD2, FD3); gc()

##-- Assign 1 landings weight per Csquare (average of low and upp range, except for highest landings category which get the lowest value).
FD = data.table()
for(iYear in unique(FishDat$Year)){
  for(iMet in unique(FishDat$benthisMet)){
    FDsub = subset(FishDat, Year == iYear & benthisMet == iMet)
    mxlo = max(FDsub$TotWt_low, na.rm=T)
    print(paste(iYear, iMet, mxlo, sep=" _ "))
    FDsub$LandWt = ifelse(FDsub$TotWt_low == mxlo, mxlo, (FDsub$TotWt_low + FDsub$TotWt_upp)/2)
    FD = rbind(FD, FDsub)
  }
}
FishDat = FD ; rm(FD, FDsub, mxlo, iMet, iYear)
FishDat$TotWt_low <- NULL ; FishDat$TotWt_upp <- NULL

FD <- subset(FishDat, !benthisMet %in% c("OT_MIX_CRU_DMF", "OT_CRU", "OT_MIX_DMF_BEN", "OT_MIX"))

##-- Combine some metiers: "OT_MIX_CRU_DMF" + "OT_CRU" 
a <- subset(FishDat, benthisMet %in% c("OT_MIX_CRU_DMF", "OT_CRU"))
a <- data.table(st_drop_geometry(a))
SDcols <- c("surface", "LandWt")
ab <- a[,(lapply(.SD, sum)), by=.(Year, C.square, Ecoregion, lat, lon), .SDcols=SDcols]
ab$benthisMet <- "OT_CRU"
FD <- rbind(FD, ab)

##-- Combine some metiers: "OT_MIX_DMF_BEN" with "OT_MIX"
a <- subset(FishDat, benthisMet %in% c("OT_MIX_DMF_BEN", "OT_MIX"))
a <- data.table(st_drop_geometry(a))
SDcols <- c("surface", "LandWt")
ab <- a[,(lapply(.SD, sum)), by=.(Year, C.square, Ecoregion, lat, lon), .SDcols=SDcols]
ab$benthisMet <- "OT_MIX"
FD <- rbind(FD, ab)

FishDat <- FD
rm(a,ab, SDcols, FD)

##-- Change formatting: all years/gears combi's get separate column
FishDatNEW <- FishDat[, c("C.square", "Ecoregion")]
FishDatNEW <- FishDatNEW[!duplicated(FishDatNEW),]

for(iYear in unique(FishDat$Year)){
  for(iMet in unique(FishDat$benthisMet)){
    subdat <- subset(FishDat, Year == iYear & benthisMet == iMet)
    FishDatNEW[[paste(iMet, "surface_SA", iYear, sep="_")]] = subdat$surface [match(FishDatNEW$C.square, subdat$C.square)]
    FishDatNEW[[paste(iMet, "total_weight", iYear, sep="_")]] = subdat$LandWt [match(FishDatNEW$C.square, subdat$C.square)]
  }
}

names(FishDatNEW)[1] <- "csquares"
saveRDS(FishDatNEW, paste0(datPath, "OSPAR0920.rds"))

#-#-# Step 2. HELCOM-fisheries data ----
HELCOMFD <- st_read(paste0(datPath, "2022_helcom_spatial_effort_intensity/shapefiles/intensity.shp"))
HELCOMFD <- subset(HELCOMFD, gear_group %in% c("DRB_MOL", "OT_CRU", "OT_DMF", "OT_MIX_CRU_DMF", "OT_SPF", "SDN_DMF"))
HELCOMFD <- data.table(st_drop_geometry(HELCOMFD))
HELCOMFD <- HELCOMFD[,c(1:3, 6, 8, 12, 13, 18, 19, 26, 27)]

##-- Clip to ICES Baltic Sea Ecoregion
HC <- HELCOMFD[, c("c_square", "lat", "lon")]
HC <- HC[!duplicated(HC),]
HC2 <- st_as_sf(HC, coords = c("lon","lat"))
st_crs(HC2) <- st_crs(ICESEcoregs)
HC2$lon = HC$lon
HC2$lat = HC$lat
HELCOMarea <- subset(ICESEcoregs, Ecoregion == "Baltic Sea")
HC3 <- st_intersection (HC2, HELCOMarea)
HELCOMFD <- subset(HELCOMFD, c_square %in% HC3$c_square)
rm(HC, HC2, HC3); gc()

##-- Assign 1 landings weight per Csquare (average of low and upp range.)
FD = data.table()
for(iYear in unique(HELCOMFD$year)){
  for(iMet in unique(HELCOMFD$gear_group)){
    FDsub = subset(HELCOMFD, gear_group == iMet & year == iYear)
    mxlo = max(FDsub$totwt_cl)
    FDsub$LandWt = ifelse(is.na(FDsub$totwt) == F, FDsub$totwt,
                          ifelse(FDsub$totwt_cl == mxlo, FDsub$totwt_cl, (FDsub$totwt_cl + FDsub$totwt_ch)/2))
    FD = rbind(FD, FDsub)
  }}
HELCOMFD = FD
HELCOMFD$totwt_cl <- NULL ; HELCOMFD$totwt_ch <- NULL ; HELCOMFD$totwt <- NULL

##-- Determine number of Csquares with no anonymous data - 139 csquares without anonymous data
ancsqs = subset(HELCOMFD, anonymous == 0)
noancsq = subset(HELCOMFD, anonymous == 1)
noancsq = subset(noancsq, !c_square %in% ancsqs$c_square)
noancsq = unique(noancsq$c_square)
saveRDS(noancsq, file=paste0(outPath, "noancsqBS.rds"))

##-- Calculate annual values
HELCOMFD <- HELCOMFD[,(lapply(.SD, sum)), by=c("year", "c_square", "gear_group"), .SDcols = c("sur", "LandWt")]

##-- Combine some metiers: "OT_MIX_CRU_DMF" + "OT_CRU" 
FD <- subset(HELCOMFD, !gear_group %in% c("OT_MIX_CRU_DMF", "OT_CRU"))
a <- subset(HELCOMFD, gear_group %in% c("OT_MIX_CRU_DMF", "OT_CRU"))
SDcols <- c("sur", "LandWt")
ab <- a[,(lapply(.SD, sum)), by=.(year, c_square), .SDcols=SDcols]
ab$gear_group <- "OT_CRU"
FD <- rbind(FD, ab)

HELCOMFD <- FD
rm(a,ab, SDcols, FD)

##--  Change format with all years/metiers in columns
HELCOMNEW <- HELCOMFD[, c("c_square")]
HELCOMNEW <- HELCOMNEW[!duplicated(HELCOMNEW),]
HELCOMNEW$Ecoregion = "Baltic Sea"

for(iYear in unique(HELCOMFD$year)){
  for(iMet in unique(HELCOMFD$gear_group)){
    subdat <- subset(HELCOMFD, year == iYear & gear_group == iMet)
    HELCOMNEW[[paste(iMet, "surface_SA", iYear, sep="_")]] = subdat$sur [match(HELCOMNEW$c_square, subdat$c_square)]
    HELCOMNEW[[paste(iMet, "total_weight", iYear, sep="_")]] = subdat$LandWt [match(HELCOMNEW$c_square, subdat$c_square)]
  }
}
names(HELCOMNEW)[1] <- "csquares"
saveRDS(HELCOMNEW, paste0(datPath, "HELCOM1621.rds"))

#-#-# Step 3. Combine all effort data into 1 file ----
OSPAR <- readRDS(file=paste0(datPath, "OSPAR0920.rds"))
HC1621 <- readRDS(file = paste0(datPath, "HELCOM1621.rds"))

##-- Remove 2009-2015 in OSPAR data
OSPAR <- OSPAR[,c(1:2, 143:242)]

##-- Remove 2021 in HELCOM data
HC1621 <- HC1621[,c(1:52)]

##-- Check colnames and add column when needed
a <- setdiff(colnames(HC1621), colnames(OSPAR)) # empty
b <- setdiff(colnames(OSPAR), colnames(HC1621))

for(i in b){
  HC1621[[i]] <- as.numeric(NA)
}
setdiff(colnames(OSPAR), colnames(HC1621)) #check
FishData <- rbind(OSPAR, HC1621)
saveRDS(FishData, paste0(datPath, "FishData1620.rds"))

#-#-# Step 4. Landing-weights capping ----
FishDataALL = readRDS(paste0(datPath, "FishData1620.rds"))
check = data.table(YG = colnames(FishDataALL[,3:102]),
                   original = colSums(FishDataALL[,3:102], na.rm=T))

##-- Remove impossible landing weights (166736 kg removed)
for(iYear in 2016:2020){
  for(iMet in c("DRB_MOL", "OT_DMF", "OT_SPF", "SDN_DMF", "SSC_DMF", "TBB_CRU", "TBB_MOL", "TBB_DMF", "OT_CRU", "OT_MIX")){
    FishDataALL[[paste0(iMet, "_total_weight_", iYear)]] = ifelse(FishDataALL[[paste0(iMet, "_total_weight_", iYear)]]>0 & FishDataALL[[paste0(iMet, "_surface_SA_", iYear)]]==0,
                                                                  0, FishDataALL[[paste0(iMet, "_total_weight_", iYear)]])
}}
check$afterremoval = colSums(FishDataALL[,3:102], na.rm=T)
check$removed = ifelse(check$original==check$afterremoval, 0, 1) ## all years and almost all gears.

##-- Determine LpUE
for(iMet in c("DRB_MOL", "OT_DMF", "OT_SPF", "SDN_DMF", "SSC_DMF", "TBB_CRU", "TBB_MOL", "TBB_DMF", "OT_CRU", "OT_MIX")){
  for(iYear in 2016:2020){
    lpue_col = paste0("LPUE_", iMet, "_", iYear)
    ## calculate LpuE
    FishDataALL[[lpue_col]] = FishDataALL[[paste0(iMet, "_total_weight_", iYear)]] / FishDataALL[[paste0(iMet, "_surface_SA_", iYear)]]
  }}

##-- Perform LpUE-capping
capping = data.table()
for(iMet in c("DRB_MOL", "OT_DMF", "OT_SPF", "SDN_DMF", "SSC_DMF", "TBB_CRU", "TBB_MOL", "TBB_DMF", "OT_CRU", "OT_MIX")){
  for(iYear in 2016:2020){
    lpue_col = paste0("LPUE_", iMet, "_", iYear)
    TLW_col =  paste0(iMet, "_total_weight_", iYear)
    SA_col =   paste0(iMet, "_surface_SA_", iYear)
    
    for(iReg in unique(FishDataALL$Ecoregion)){
      ## Determine outlier LpUE values
      idx = which(FishDataALL$Ecoregion == iReg)
      x = FishDataALL[[lpue_col]][idx] 
      bp = boxplot.stats(x)
      upper_w = bp$stats[5]
      
      ## Cap outliers to max whisker value
      x_capped = ifelse(x > upper_w, upper_w, x)
      FishDataALL[[lpue_col]][idx] = x_capped
      
      ## Correct affected weights
      OrigTLW = FishDataALL[[TLW_col]][idx]
      FishDataALL[[TLW_col]][idx] = FishDataALL[[SA_col]][idx] * FishDataALL[[lpue_col]][idx]
      
      capreg = data.table(gear = iMet,
                          year = iYear,
                          reg = iReg,
                          orig_LPUE = sum(x, na.rm=T),
                          cap_LPUE = sum(x_capped, na.rm=T),
                          nr_cap = length(bp$out),
                          orig_TLW = sum(OrigTLW, na.rm=T),
                          cap_TLW = sum(FishDataALL[[TLW_col]][idx], na.rm=T),
                          TLW_diff = sum(OrigTLW, na.rm=T) - sum(FishDataALL[[TLW_col]][idx], na.rm=T))
      capping = rbind(capping, capreg)
    }}}

##-- Determine csquare-size
csquares = FishDataALL[,c(1:2)]
csquares = st_as_sf.csquares(csquares)
csquares$area_sqkm = as.numeric(st_area(csquares))/1E6

##-- Determine annual SA and TLW
FD = FishDataALL
for(iYear in 2016:2020){
  idx = grep(iYear, colnames(FD))
  yearsub = FD[,..idx]
  idx = grep("surface_SA", colnames(yearsub))
  idy = grep("total_weight", colnames(yearsub))
  saryear = rowSums(yearsub[,..idx], na.rm=T)
  tlwyear = rowSums(yearsub[,..idy], na.rm=T)
  FD[[paste0("SA", iYear)]] = saryear
  FD[[paste0("TLW", iYear)]] = tlwyear
}
SAcols = paste0("SA", 2016:2020) ; TLWcols = paste0("TLW", 2016:2020)
FD$AV_SA1620 = rowMeans(FD[,..SAcols], na.rm=T)
FD$AV_TLW1620 = rowMeans(FD[,..TLWcols], na.rm = T)
FD$area_sqkm = csquares$area_sqkm [match(FD$csquares, csquares$csquares)]
FD$AV_SAR1620 = FD$AV_SA1620 / FD$area_sqkm
saveRDS(FD, file=paste0(datPath, "FD_capped.rds"))

##-- Create FishdataANN (with spatial info)
csquares = csquares[,1:3]
FD_ann = merge(csquares, FD_capped[,c(1:2, 153:166)], by=c("csquares", "Ecoregion"))
saveRDS(FD_ann, file=paste0(datPath, "FD_ann.rds"))



