#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
##  Rscript for the analyse underlying the article 
#      "Evaluating the efficacy of closing the current European network of marine protected areas 
#       to mobile bottom-contacting fishing gears"
#
##  Step 2: Perform spatial overlay between csquares, MPAs, and MSFD BHT.
#
##  Code by Karin van der Reijden. (kjvdreijden@gmail.com)
#        published under the MIT-license
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Step 0. Load libraries, set data paths, load ICES Ecoregions ----
library(sf)
library(data.table)
library(sp)
library(terra)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(basemaps)
sf_use_s2(FALSE)

##-- set data paths
datPath <- "D:/MS_Natura2000closures/Data/"
outPath <- "D:/MS_Natura2000closures/Routput/"

##-- Load ICES Ecoregions, and clip to the ones of interest
# Downloadable here: https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb
ICESEcoregs <- st_read(paste0(datPath, "ICES_ecoregions_20171207_erase_ESRI.shp"))
ICESEcoregs <- subset(ICESEcoregs, Ecoregion %in% c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast"))
ICESEcoregs <- ICESEcoregs[,c("Ecoregion")]

#-#-# Step 1. Obtain habitat sensitivity data ----- 
##-- The file used here was recieved from PD van Denderen, but is available here: https://zenodo.org/records/15176198)
load(paste0(datPath, "250107_outputs_combined/Infauna_sens_fishing_state.Rdata"))
inf$Ecoregion <- as.character(inf$Ecoregion)
inf$EEZ = as.character(inf$EEZ)
inf <- subset(inf, Ecoregion %in% c("Baltic Sea", "Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast"))
load(paste0(datPath, "250107_outputs_combined/Epifauna_sens_fishing_state.Rdata"))
epi$Ecoregion <- as.character(epi$Ecoregion)
epi$EEZ = as.character(epi$EEZ)
epi <- subset(epi, Ecoregion %in% c("Celtic Seas", "Greater North Sea", "Bay of Biscay and the Iberian Coast"))

##-- Get unique csquares, with long, lat, depth, ecoregion
csquares1 = data.table(inf[,c(1:7)])
csquares2 = data.table(epi[,c(1:7)])
csquares = rbind(csquares1, csquares2)
csquares = csquares[!duplicated(csquares),]

##-- Remove Norwegian and Russian waters
csquares = subset(csquares, !EEZ %in% c("Russia", "Norway"))

##-- Combine epi- & infauna
csquares$inf_medlong = inf$medlong [match(csquares$csquares, inf$csquares)]
csquares$inf_intercept= inf$intercept [match(csquares$csquares, inf$csquares)]
csquares$inf_slope = inf$slope [match(csquares$csquares, inf$csquares)]
csquares$epi_medlong = epi$medlong [match(csquares$csquares, epi$csquares)]
csquares$epi_intercept= epi$intercept [match(csquares$csquares, epi$csquares)]
csquares$epi_slope = epi$slope [match(csquares$csquares, epi$csquares)]

saveRDS(csquares, file = paste0(datPath, "Regiontot.rds"))

##-- create a sf-dataset for regiontot as well.
##-- function to create polygons
create_square_sf <- function(long, lat, offset = 0.025) {
  square_coords <- matrix(c(
    long - offset, lat - offset,
    long + offset, lat - offset,
    long + offset, lat + offset,
    long - offset, lat + offset,
    long - offset, lat - offset
  ), ncol = 2, byrow = TRUE)
  
  square <- st_polygon(list(square_coords))
  square_sf <- st_sfc(square, crs = 4326)
  return(square_sf)
}

csquares_sf <- csquares %>%
  rowwise() %>%
  mutate(geometry = create_square_sf(long, lat)) %>%
  st_as_sf()
saveRDS(csquares_sf, file=paste0(datPath, "RegiontotSF.rds"))

#-#-# Step 1a. Determine swept area % in and outside assessment area ----
# FishDataALL = readRDS(file=paste0(datPath, "FD_capped.rds"))
# Regiontot = readRDS(file = paste0(datPath, "Regiontot.rds"))
# Regiontot$AV_SA = FishDataALL$AV_SA1620 [ match(Regiontot$csquares, FishDataALL$csquares)]
# 
# ecoreg = aggregate(AV_SA1620~ Ecoregion, data=FishDataALL, FUN="sum")
# names(ecoreg) = c("Ecoregion", "tot_SA")
# assare = aggregate(AV_SA~ Ecoregion, data=Regiontot, FUN= "sum")
# ecoreg$ass_SA = assare$AV_SA [ match(ecoreg$Ecoregion, assare$Ecoregion)]
# ecoreg$SAprop = round(ecoreg$ass_SA / ecoreg$tot_SA * 100, digits=2)
# 
# # Ecoregion                             tot_SA      ass_SA     SAprop
# # Baltic Sea                            74405.45     74254.45  99.80
# # Bay of Biscay and the Iberian Coast   358162.97   357254.76  99.75
# # Celtic Seas                           937720.29   852588.25  90.92
# # Greater North Sea                     1054185.70  912529.11  86.56

#-#-# Step 2. Load BHT data and merge with csquares ----
##-- Load BHT data - crop to ICES Ecoregions (available here: https://emodnet.ec.europa.eu/geonetwork/srv/eng/catalog.search#/metadata/0a1cb988-22de-48b2-8cda-d90947ef77d1)
EUSeaMap <- st_read(dsn=paste0(datPath, "EUSeaMap_2023.gdb"))
EUSeaMap <- EUSeaMap[,c("MSFD_BBHT", "Shape")]
EUSeaMap <- st_crop(EUSeaMap, st_bbox(ICESEcoregs)) # set to boundaries of ICESEcoregs
saveRDS(EUSeaMap, paste0(outPath, "EUSeaMap.rds"))

##-- Load BHT data and match to csquares from sensitivity data
Regiontot <- readRDS(paste0(datPath, "RegiontotSF.rds"))
csquaresMSFDBHT <- st_intersection(Regiontot, EUSeaMap)

##-- Determine area per CS and BHT
csquaresMSFDBHT <- csquaresMSFDBHT %>%
  group_by(csquares, MSFD_BBHT) %>%
  summarise()
csquaresMSFDBHT$BHTareakm2 <- as.numeric(st_area(st_make_valid(csquaresMSFDBHT))/1e6)
CSarea <- aggregate(BHTareakm2~csquares, data=csquaresMSFDBHT, FUN="sum")
csquaresMSFDBHT$CSareaNEW <- CSarea$BHTareakm2 [match(csquaresMSFDBHT$csquares, CSarea$csquares)]
csquaresMSFDBHT$Ecoregion <- Regiontot$Ecoregion [match(csquaresMSFDBHT$csquares, Regiontot$csquares)]

##-- Save as SF and data.frame
saveRDS(csquaresMSFDBHT, paste0(datPath, "csquaresMSFDBHT.rds"))
RegionMSFD <- st_drop_geometry(csquaresMSFDBHT)
saveRDS(RegionMSFD, file=paste0(datPath, "RegionMSFD.rds"))

#-#-# Step 3. Load N2000 and UK MPA sites and merge with csquares ----
##-- N2000 sites (available here: https://sdi.eea.europa.eu/data/be2142b0-7dc4-42e3-afac-f7cc5f1a9ac6)
N2000gp = paste0(datPath, "N2000_sites/Natura2000_end2023.gpkg")
N2000 = st_read(N2000gp, layer="NaturaSite_polygon")
N2000 <- st_transform(N2000, 4326)
N2000 <- N2000[,c("SITECODE", "SITENAME", "SITETYPE")]
st_geometry(N2000) = "geometry"
names(N2000) <- c("MPA_ID", "Name", "MPA_status", "geometry")
N2000$source = "EEA"

##-- JNCC UK MPA sites (available here: https://jncc.gov.uk/resources/ade43f34-54d6-4084-b66a-64f0b4a5ef27)
UKMPA = st_read(dsn = paste0(datPath, "UK_MPAs"), layer="c20230705_OffshoreMPAs_WGS84")
UKMPA <- st_transform(UKMPA, 4326)
UKMPA <- UKMPA[,c("SITE_CODE", "SITE_NAME")]
names(UKMPA) <- c("MPA_ID", "Name", "geometry")
UKMPA$source = "JNCC"

##-- UK sites: Determine feature type protected: Hab when either: habitat, geology, ecosystem, or specific (fish/benthic species)
# based on spreasheet UK Offshore MPA spreadsheet XLSX: UKOffshoreMPAs-20230705-LIVE-WEB.xlsx (accessed: 02/04/2025)
Specs = c("UKNCMPA021","UKMCZ0023","UKMCZ0043","UKNCMPA022","UKNCMPA018",
          "UKMCZ0046","UKMCZ0078","UKMCZ0024","UKNCMPA019","UKNCMPA026",
          "UKMCZ0082","UKMCZ0025","UKNCMPA028","UKNCMPA029","UKMPA0001")
Geo = c("UKNCMPA020","UKNCMPA022","UKNCMPA018","UKNCMPA023","UKMCZ0047",
        "UKNCMPA024","UKMCZ0078","UKNCMPA025","UKNCMPA019","UKMCZ0044",
        "UKMCZ0089","UKMCZ0025","UKMCZ0026","UKNCMPA028","UKMPA0001")
Ecosys = c("UKEHPMA002","UKEHPMA003")
Hab = c("UK0030387","UK0030368","UK0030357","UKMCZ0076","UKNCMPA020","UK0030381","UK0030317","UK0030352",
        "UKNCMPA021","UKMCZ0023","UKMCZ0077","UK0030389","UKMCZ0043","UKNCMPA022","UKNCMPA018","UKMCZ0046",
        "UKNCMPA023","UKMCZ0047","UK0030353","UK0030369","UK0030388","UKNCMPA024","UKMCZ0078","UKMCZ0079",
        "UK0030370","UKMCZ0080","UKMCZ0084","UKMCZ0024","UK0030358","UK0030363","UKNCMPA025","UKMCZ0085",
        "UKMCZ0048","UKMCZ0049","UKMCZ0044","UKMCZ0081","UK0030379","UK0030385","UKMCZ0086","UK0030354",
        "UK0030386","UKMCZ0022","UKMCZ0087","UKMCZ0082","UKMCZ0088","UKMCZ0083","UKMCZ0089","UKMCZ0025",
        "UK0030359","UKMCZ0026","UKNCMPA028","UKMCZ0027","UKMCZ0090","UKMPA0001","UKMCZ0045","UKMCZ0091",
        "UKNCMPA030","UKMCZ0050","UK0030380","UK0030355")
UKMPA$MPA_status = ifelse(UKMPA$MPA_ID %in% c(Specs, Geo, Ecosys, Hab), "Hab", "noHab")

##-- Combine all MPA sites into 1 file
MPAs = rbind(N2000, UKMPA)
saveRDS(MPAs, file=paste0(outPath, "MPAs_dataset.rds"))
rm(N2000, UKMPA)

##-- Select all csquares that intersect with MPAs
Regiontot   <- readRDS(paste0(datPath, "RegiontotSF.rds"))
grid_over   <- st_intersects(st_make_valid(Regiontot),st_make_valid(MPAs))
grid_over   <- as.data.frame(grid_over)
csquaresSUB <- Regiontot[grid_over$row.id,] # only keep csquares with MPA match
csquaresSUB <- csquaresSUB[!duplicated(csquaresSUB$csquares), ] # only keep unique csquare-MPA matches.

##-- Obtain actual intersection for all MPA-data
grid_over   <- st_intersection(csquaresSUB, st_make_valid(MPAs))
saveRDS(grid_over, file=paste0(datPath, "grid_over.rds"))

##-- Merge N2000 with BHT and csquares
csquaresMSFDBHT <- readRDS(paste0(datPath, "csquaresMSFDBHT.rds"))
MSFDpoly <- csquaresMSFDBHT[,c("csquares", "MSFD_BBHT", "geometry")]
MSFDpoly <- subset(MSFDpoly, csquares %in% grid_over$csquares) # only include csquares with MPAs.
CS_BHT_N2000<- st_intersection(st_make_valid(grid_over), st_make_valid(MSFDpoly))
CS_BHT_N2000$csquares.1 <- NULL
saveRDS(CS_BHT_N2000, file=paste0(datPath, "CS_BHT_N2000.rds"))

#-#-# Step 4. Determine % cover of MPAs per csquare; ALL & HPS separately ----
##-- Set 1: all MPAs
Regiontot <- readRDS(paste0(datPath, "Regiontot.rds"))
grid_over <- readRDS(file=paste0(datPath, "grid_over.rds"))
grid_overALL <- grid_over %>%
  group_by(csquares) %>%
  summarise()
grid_overALL$S1area <- as.numeric(st_area(st_make_valid(grid_overALL)) / 1e6)
Regiontot$S1area <- grid_overALL$S1area [match(Regiontot$csquares, grid_overALL$csquares)]
Regiontot$S1area[is.na(Regiontot$S1area)] <- 0
Regiontot$S1pct  <- Regiontot$S1area/Regiontot$area_sqkm

##-- Set 2: Habitat-protection sites only 
grid_overS2 <- subset(grid_over, MPA_status %in% c("B", "C", "Hab")) %>%
  group_by(csquares) %>%
  summarise()
grid_overS2$S2area <- as.numeric(st_area(st_make_valid(grid_overS2)) / 1e6)
Regiontot$S2area <- grid_overS2$S2area [match(Regiontot$csquares, grid_overS2$csquares)]
Regiontot$S2area[is.na(Regiontot$S2area)] <- 0
Regiontot$S2pct  <- Regiontot$S2area/Regiontot$area_sqkm

saveRDS(Regiontot, file=paste0(datPath, "Regiontot2.rds"))

#-#-# Step 5. Determine % cover of MPAs per csquare per MSFD-BHT; ALL & HPS separately ----
##-- Set 1: all MPAs
CS_BHT_N2000 <- readRDS(paste0(datPath, "CS_BHT_N2000.rds"))
RegionMSFD <- readRDS(paste0(datPath, "RegionMSFD.rds"))
grid_over2S1 <- CS_BHT_N2000 %>%
  group_by(csquares, MSFD_BBHT) %>%
  summarise()
grid_over2S1$S1area <- as.numeric(st_area(st_make_valid(grid_over2S1))/ 1e6)
grid_over2S1$CS_MSFD <- paste(grid_over2S1$csquares, grid_over2S1$MSFD_BBHT, sep="_")
RegionMSFD$CS_MSFD <- paste(RegionMSFD$csquares, RegionMSFD$MSFD_BBHT, sep="_")
RegionMSFD$S1area <- grid_over2S1$S1area [match(RegionMSFD$CS_MSFD, grid_over2S1$CS_MSFD)]
RegionMSFD$S1area[is.na(RegionMSFD$S1area)] <- 0
RegionMSFD$S1pct <- RegionMSFD$S1area / RegionMSFD$BHTareakm2

##-- Set 2: Habitat Directive only
grid_over2S2 <- subset(CS_BHT_N2000, MPA_status %in% c("B", "C", "Hab")) %>%
  group_by(csquares, MSFD_BBHT) %>%
  summarise()
grid_over2S2$S2area <- as.numeric(st_area(st_make_valid(grid_over2S2)) / 1e6)
grid_over2S2$S2_MSFD <- paste(grid_over2S2$csquares, grid_over2S2$MSFD_BBHT, sep="_")
RegionMSFD$S2area <- grid_over2S2$S2area [match(RegionMSFD$CS_MSFD, grid_over2S2$S2_MSFD)]
RegionMSFD$S2area[is.na(RegionMSFD$S2area)] <- 0
RegionMSFD$S2pct <- RegionMSFD$S2area / RegionMSFD$BHTareakm2

saveRDS(RegionMSFD, file=paste0(datPath, "RegionMSFD2.rds"))


#-#-# Step 6. Determine overlap Ecoregion with MSFD BHT and MPAs ----
##-- Load BHT data and match to ICES Ecoregions
EUSeaMap <- readRDS(paste0(outPath, "EUSeaMap.rds"))
MSFDBHTreg <- st_intersection(ICESEcoregs, EUSeaMap)
MSFDBHTreg$area_sqkm = as.numeric(st_area(MSFDBHTreg))/1E6
saveRDS(MSFDBHTreg, file=paste0(outPath, "MSFDBHTreg.rds"))

dt = data.table(aggregate(area_sqkm~ Ecoregion+MSFD_BBHT, data=MSFDBHTreg, FUN="sum"))
names(dt) = c("Ecoregion", "MSFD_BBHT", "Total_area")
dt$ID = paste(dt$Ecoregion, dt$MSFD_BBHT, sep="_")

MPAs = readRDS(file=paste0(outPath, "MPAs_dataset.rds"))
MSFDBHTregMPA = st_intersection(MSFDBHTreg, MPAs)
saveRDS(MSFDBHTregMPA, file=paste0(outPath, "MSFDBHTregMPA.rds"))

##-- Obtain all MPAs (correct for overlap)
MSFDBHTregMPA_ALL <- MSFDBHTregMPA %>%
  group_by(Ecoregion, MSFD_BBHT) %>%
  summarise()

MSFDBHTregMPA_ALL$area = as.numeric(st_area(st_make_valid(MSFDBHTregMPA_ALL)) / 1e6)
MSFDBHTregMPA_ALL$ID = paste(MSFDBHTregMPA_ALL$Ecoregion, MSFDBHTregMPA_ALL$MSFD_BBHT, sep="_")
dt$ALL_MPA_area = MSFDBHTregMPA_ALL$area [match(dt$ID, MSFDBHTregMPA_ALL$ID)]

##-- Add HPS MPAs (correct for overlap)
MSFDBHTregMPA_HPS <- subset(MSFDBHTregMPA, MPA_status %in% c("B", "C", "Hab")) %>%
  group_by(Ecoregion, MSFD_BBHT) %>%
  summarise()
MSFDBHTregMPA_HPS$area = as.numeric(st_area(st_make_valid(MSFDBHTregMPA_HPS)) / 1e6)
MSFDBHTregMPA_HPS$ID = paste(MSFDBHTregMPA_HPS$Ecoregion, MSFDBHTregMPA_HPS$MSFD_BBHT, sep="_")
dt$HPS_MPA_area = MSFDBHTregMPA_HPS$area [match(dt$ID, MSFDBHTregMPA_HPS$ID)]

##-- Calculate proportions
totext = aggregate(cbind(Total_area, ALL_MPA_area, HPS_MPA_area)~ Ecoregion, FUN="sum", data=dt)
totext$MSFD_BBHT = "Total"
totext$MPAHabprop = round(totext$ALL_MPA_area / totext$Total_area * 100, digits=2)
totext$HPSHabprop = round(totext$HPS_MPA_area / totext$Total_area * 100, digits=2)
totext$TotExt = totext$Total_area
totext$HabProp = 100
dt$TotExt = totext$Total_area [match(dt$Ecoregion, totext$Ecoregion)]

dt$HabProp = round(dt$Total_area / dt$TotExt * 100, digits=2)
dt$MPAHabprop = round(dt$ALL_MPA_area / dt$Total_area * 100, digits=2)
dt$HPSHabprop = round(dt$HPS_MPA_area / dt$Total_area * 100, digits=2)
dt$ID = NULL
dt = rbind(dt, totext)
saveRDS(dt, file=paste0(outPath, "reg_MPABHT.rds"))

##-- Create table S1 for all ecoregion extent
TabH <- data.table(MSFD = c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
                            "Infralittoral mixed sediment", "Infralittoral sand", 
                            "Infralittoral mud", "Infralittoral mud or Infralittoral sand", 
                            "Circalittoral rock and biogenic reef", "Circalittoral coarse sediment",
                            "Circalittoral mixed sediment", "Circalittoral sand", 
                            "Circalittoral mud", "Circalittoral mud or Circalittoral sand",
                            "Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                            "Offshore circalittoral mixed sediment", "Offshore circalittoral sand", 
                            "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand",
                            "Upper bathyal rock and biogenic reef", "Upper bathyal sediment", 
                            "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                            "Lower bathyal rock and biogenic reef", "Lower bathyal sediment", 
                            "Lower bathyal sediment or Lower bathyal rock and biogenic reef",
                            "Abyssal", "na", "Total"), 
                   habOrd = 1:27, 
                   MSFDName = c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
                                "Infralittoral mixed sediment", "Infralittoral sand", 
                                "Infralittoral mud", "Infralittoral mud or sand", 
                                "Circalittoral rock and biogenic reef", "Circalittoral coarse sediment",
                                "Circalittoral mixed sediment", "Circalittoral sand", 
                                "Circalittoral mud", "Circalittoral mud or sand",
                                "Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                "Offshore circalittoral mixed sediment", "Offshore circalittoral sand", 
                                "Offshore circalittoral mud", "Offshore circalittoral mud or sand",
                                "Upper bathyal rock and biogenic reef", "Upper bathyal sediment", 
                                "Upper bathyal sediment or rock and biogenic reef",
                                "Lower bathyal rock and biogenic reef", "Lower bathyal sediment", 
                                "Lower bathyal sediment or rock and biogenic reef",
                                "Abyssal", "Unknown", "Total"))
regs = data.frame(Ecoregion =c("Baltic Sea", "Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast"),
                  code = c("BS", "GNS", "CS", "BI"))

for(iReg in c("Baltic Sea", "Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
  subdat = subset(dt, Ecoregion==iReg)
  TabH[[paste0(subset(regs, Ecoregion==iReg)$code, "_Ext_km2")]] = round(subdat$Total_area [match(TabH$MSFD, subdat$MSFD_BBHT)] /1000, digits=2)
  TabH[[paste0(subset(regs, Ecoregion==iReg)$code, "_Ext_pct")]] = subdat$HabProp [match(TabH$MSFD, subdat$MSFD_BBHT)]
  TabH[[paste0(subset(regs, Ecoregion==iReg)$code, "_MPA_pct")]] = subdat$MPAHabprop [match(TabH$MSFD, subdat$MSFD_BBHT)]
  TabH[[paste0(subset(regs, Ecoregion==iReg)$code, "_HPS_pct")]] = subdat$HPSHabprop [match(TabH$MSFD, subdat$MSFD_BBHT)]
}
saveRDS(TabH, file=paste0(outPath, "tabS1allextreg.rds"))

#-#-# Step 7. Overlap fishing with BHT in Ecoregion ----
MSFDBHTreg = readRDS(file=paste0(outPath, "MSFDBHTreg.rds"))
FishDataANN = readRDS(file=paste0(datPath, "FD_ann.rds"))

sf1 <- st_sf(
  geometry = st_sfc(),  # Empty geometry column
  crs = 4326            # Specify CRS (WGS84)
)
for(iReg in unique(MSFDBHTreg$Ecoregion)){
  BHT = subset(MSFDBHTreg, Ecoregion == iReg)
  FD = subset(FishDataANN, Ecoregion == iReg)
  int = st_intersection(BHT, FD)
  sf1 = rbind(sf1, int)
}
sf1$Ecoregion.1 = NULL
sf1$BHT_area_sqkm = as.numeric(st_area(sf1))/1E6 
sf1$area_sqkm = NULL ; sf1$area_sqkm = sf1$area_sqkm.1 ; sf1$area_sqkm.1 = NULL
saveRDS(sf1, file=paste0(outPath, "sf1.rds"))

##-- Aggregate fished area and TLW per BHT per ecoregion
sf1$BHT_SA = sf1$AV_SA1620 * (sf1$BHT_area_sqkm/sf1$area_sqkm)
sf1$BHT_TLW = sf1$AV_TLW * (sf1$BHT_area_sqkm/sf1$area_sqkm)
BHTfished = aggregate(cbind(BHT_area_sqkm, BHT_SA, BHT_TLW)~ MSFD_BBHT+Ecoregion, data=sf1, FUN="sum")
totfished = aggregate(cbind(BHT_SA, BHT_TLW)~ Ecoregion, data=sf1, FUN="sum")
BHTfished$Tot_SA = totfished$BHT_SA [match(BHTfished$Ecoregion, totfished$Ecoregion)]
BHTfished$Tot_TLW = totfished$BHT_TLW [match(BHTfished$Ecoregion, totfished$Ecoregion)]
BHTfished$PropSA_BHT = round(BHTfished$BHT_SA / BHTfished$Tot_SA * 100, digits=2)
BHTfished$PropTLW_BHT = round(BHTfished$BHT_TLW / BHTfished$Tot_TLW * 100, digits=2)
saveRDS(BHTfished, file=paste0(outPath, "BHTfished.rds"))

##-- Aggregate fished area per ecoregion
FishedReg = aggregate(BHT_area_sqkm~Ecoregion, data=sf1, FUN="sum")
totext = aggregate(area_sqkm ~ Ecoregion, data=MSFDBHTreg, FUN="sum")
FishedReg$totExt = totext$area_sqkm [match(FishedReg$Ecoregion, totext$Ecoregion)]
FishedReg$Propfished = FishedReg$BHT_area_sqkm / FishedReg$totExt * 100
saveRDS(FishedReg, file=paste0(outPath, "fishedREG.rds"))

##-- Aggregate extent and swept area per sed and depth class per ecoregion
sf1$Sediment = ifelse(sf1$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment", "Circalittoral rock and biogenic reef",
                                           "Circalittoral coarse sediment", "Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                           "Upper bathyal rock and biogenic reef", "Lower bathyal rock and biogenic reef",
                                           "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                           "Lower bathyal sediment or Lower bathyal rock and biogenic reef"), "Gravel",
                      ifelse(sf1$MSFD_BBHT %in% c("Infralittoral sand",
                                                  "Circalittoral sand", 
                                                  "Offshore circalittoral sand"), "Sand", 
                             ifelse(sf1$MSFD_BBHT %in% c("Infralittoral mud", "Infralittoral mud or Infralittoral sand",
                                                         "Circalittoral mud", "Circalittoral mud or Circalittoral sand",
                                                         "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand",
                                                         "Upper bathyal sediment", "Lower bathyal sediment", "Abyssal"), "Mud", 
                                    ifelse(sf1$MSFD_BBHT %in% c("Infralittoral mixed sediment", "Circalittoral mixed sediment", "Offshore circalittoral mixed sediment"), "Mixed", "Unknown"))))
sf1$Depthclass = ifelse(sf1$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
                                             "Infralittoral mixed sediment", "Infralittoral sand", 
                                             "Infralittoral mud", "Infralittoral mud or Infralittoral sand"), "Infralittoral", 
                        ifelse(sf1$MSFD_BBHT %in% c("Circalittoral rock and biogenic reef", "Circalittoral coarse sediment",
                                                    "Circalittoral mixed sediment", "Circalittoral sand", 
                                                    "Circalittoral mud", "Circalittoral mud or Circalittoral sand"), "Circalittoral", 
                               ifelse(sf1$MSFD_BBHT %in% c("Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                                           "Offshore circalittoral mixed sediment", "Offshore circalittoral sand", 
                                                           "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand"), "Offshore circalittoral", 
                                      ifelse(sf1$MSFD_BBHT %in% c("Upper bathyal rock and biogenic reef", "Upper bathyal sediment", 
                                                                        "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                                                        "Lower bathyal rock and biogenic reef", "Lower bathyal sediment", 
                                                                        "Lower bathyal sediment or Lower bathyal rock and biogenic reef",
                                                                        "Abyssal"), "Bathyal & Abyssal", "Unknown"))))

sedfished = aggregate(cbind(BHT_area_sqkm, BHT_SA)~ Sediment, data=sf1, FUN="sum")
dptfished = aggregate(cbind(BHT_area_sqkm, BHT_SA) ~ Depthclass, data=sf1, FUN="sum")
sedfished$totSA = sum(totfished$BHT_SA)
dptfished$totSA = sum(totfished$BHT_SA)
sedfished$PropSA_sed = round(sedfished$BHT_SA / sedfished$totSA * 100, digits=2)
dptfished$PropSA_dpt = round(dptfished$BHT_SA / dptfished$totSA * 100, digits=2)

##-- Determine total extent per variable for ecoregion
MSFDBHTreg$Sediment = ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment", "Circalittoral rock and biogenic reef",
                                                         "Circalittoral coarse sediment", "Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                                         "Upper bathyal rock and biogenic reef", "Lower bathyal rock and biogenic reef",
                                                         "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                                         "Lower bathyal sediment or Lower bathyal rock and biogenic reef"), "Gravel",
                             ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Infralittoral sand",
                                                                "Circalittoral sand", 
                                                                "Offshore circalittoral sand"), "Sand", 
                                    ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Infralittoral mud", "Infralittoral mud or Infralittoral sand",
                                                                       "Circalittoral mud", "Circalittoral mud or Circalittoral sand",
                                                                       "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand",
                                                                       "Upper bathyal sediment", "Lower bathyal sediment", "Abyssal"), "Mud", 
                                           ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Infralittoral mixed sediment", "Circalittoral mixed sediment", "Offshore circalittoral mixed sediment"), "Mixed", "Unknown")))) 
MSFDBHTreg$Depthclass = ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
                                                           "Infralittoral mixed sediment", "Infralittoral sand", 
                                                           "Infralittoral mud", "Infralittoral mud or Infralittoral sand"), "Infralittoral", 
                               ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Circalittoral rock and biogenic reef", "Circalittoral coarse sediment",
                                                                  "Circalittoral mixed sediment", "Circalittoral sand", 
                                                                  "Circalittoral mud", "Circalittoral mud or Circalittoral sand"), "Circalittoral", 
                                      ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                                                         "Offshore circalittoral mixed sediment", "Offshore circalittoral sand", 
                                                                         "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand"), "Offshore circalittoral", 
                                             ifelse(MSFDBHTreg$MSFD_BBHT %in% c("Upper bathyal rock and biogenic reef", "Upper bathyal sediment", 
                                                                                "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                                                                "Lower bathyal rock and biogenic reef", "Lower bathyal sediment", 
                                                                                "Lower bathyal sediment or Lower bathyal rock and biogenic reef",
                                                                                "Abyssal"), "Bathyal & Abyssal", "Unknown"))))

dptreg = aggregate(area_sqkm ~ Depthclass, data=MSFDBHTreg, FUN="sum")
sedreg = aggregate(area_sqkm ~ Sediment, data=MSFDBHTreg, FUN="sum")
sedfished$sedarea = sedreg$area_sqkm [match(sedfished$Sediment, sedreg$Sediment)]
sedfished$totprop = sedfished$sedarea / sum(sedreg$area_sqkm) * 100
dptfished$dptarea = dptreg$area_sqkm [match(dptfished$Depthclass, dptreg$Depthclass)]
dptfished$totprop = dptfished$dptarea / sum(dptreg$area_sqkm) * 100

saveRDS(sedfished, file=paste0(outPath, "sedfishedREG.rds"))
saveRDS(dptfished, file=paste0(outPath, "dptfishedREG.rds"))

#-#-# Step 8. Overlap fishing with MPA in Ecoregion ----
##-- Select all csquares that intersect with MPAs
MPAs = readRDS(file=paste0(outPath, "MPAs_dataset.rds"))
FishDataANN = readRDS(file=paste0(datPath, "FD_ann.rds"))

grid_over   <- st_intersects(st_make_valid(FishDataANN), st_make_valid(MPAs))
grid_over   <- as.data.frame(grid_over)
FDA         <- FishDataANN[grid_over$row.id,] # only keep csquares with MPA match
FDA         <- FDA[!duplicated(FDA$csquares), ] # only keep unique csquare-MPA matches.

##-- Obtain actual intersection for all MPAVME-data
grid_over   <- st_intersection(FDA, st_make_valid(MPAs))
saveRDS(grid_over, file=paste0(datPath, "grid_overREG.rds")) 

##-- Determine MPA area per csquare, BHT and ecoregion (correct for potential overlap)
allMPA_reg <- subset(grid_over) %>%
  group_by(csquares) %>%
  summarise()
allMPA_reg$MPA_area = as.numeric(st_area(allMPA_reg))/1E6

##-- Determine HPS area per csquare, BHT and ecoregion (correct for potential overlap)
hpsMPA_reg <- subset(grid_over, MPA_status %in% c("B", "C", "Hab")) %>%
  group_by(csquares) %>%
  summarise()
hpsMPA_reg$HPS_area = as.numeric(st_area(hpsMPA_reg))/1E6

##-- Add to FishDataANN table
FE_MPA = as.data.table(st_drop_geometry(readRDS(file=paste0(datPath, "FD_ann.rds"))))
FE_MPA$MPA_area = allMPA_reg$MPA_area [match(FE_MPA$csquares, allMPA_reg$csquares)]
FE_MPA$MPA_pct = FE_MPA$MPA_area / FE_MPA$area_sqkm
FE_MPA$HPS_area = hpsMPA_reg$HPS_area [match(FE_MPA$csquares, hpsMPA_reg$csquares)]
FE_MPA$HPS_pct = FE_MPA$HPS_area / FE_MPA$area_sqkm
FE_MPA[is.na(FE_MPA)] = 0

##-- Determine SA and TLW per csquare in MPA and HPS areas
FE_MPA$MPA_SA = FE_MPA$AV_SA1620 * FE_MPA$MPA_pct
FE_MPA$MPA_TLW = FE_MPA$AV_TLW1620 * FE_MPA$MPA_pct
FE_MPA$HPS_SA = FE_MPA$AV_SA1620 * FE_MPA$HPS_pct
FE_MPA$HPS_TLW = FE_MPA$AV_TLW1620 * FE_MPA$HPS_pct

##-- Calculate total AV and TLW per ecoregion in MPA and HPS areas 
MPASATLW = aggregate(data=FE_MPA, cbind(AV_SA1620, AV_TLW1620, MPA_SA, MPA_TLW, HPS_SA, HPS_TLW)~ Ecoregion, FUN="sum")
MPASATLW$HPS = MPASATLW$HPS_SA / MPASATLW$AV_SA * 100
MPASATLW$MPA = (MPASATLW$MPA_SA - MPASATLW$HPS_SA) / MPASATLW$AV_SA * 100
MPASATLW$noMPA = (MPASATLW$AV_SA - MPASATLW$MPA_SA)/ MPASATLW$AV_SA * 100
MPASATLW$HPSL = MPASATLW$HPS_TLW / MPASATLW$AV_TLW * 100
MPASATLW$MPAL = (MPASATLW$MPA_TLW - MPASATLW$HPS_TLW) / MPASATLW$AV_TLW * 100
MPASATLW$noMPAL = (MPASATLW$AV_TLW - MPASATLW$MPA_TLW)/ MPASATLW$AV_TLW * 100
saveRDS(MPASATLW, file=paste0(outPath, "MPASATLWreg.rds"))

#-#-# Step 9. Create figure 1 ----
##-- Load data, settings etc
RegiontotSF = readRDS(paste0(datPath, "RegiontotSF.rds")) 
Regiontot = readRDS(paste0(datPath, "Regiontot.rds")) 
plotcols <- data.frame(colcode = c("#69696980", brewer.pal("Purples", n=6)),
                       valrange = 1:7)
FishDataANN = readRDS(paste0(datPath, "FD_ann.rds"))
FishDataANN$AVSARclass <- cut(FishDataANN$AV_SAR1620, breaks=c(-1, 0, 0.1, 0.5, 1, 5, 10, 100000), labels=FALSE)
FishDataANN$SARplotcol = plotcols$colcode [match(FishDataANN$AVSARclass, plotcols$valrange)]
RegMap <- basemap_terra(st_bbox(st_buffer(ICESEcoregs, 0.25)), map_service ="esri", map_type="world_ocean_base")
RegMap <- project(RegMap, "epsg:4326")
Regiontot$medlong = rowMeans(Regiontot[,c("inf_medlong", "epi_medlong")], na.rm=T)
RegiontotSF$Medlong = Regiontot$medlong [match(RegiontotSF$csquares, Regiontot$csquares)]
medlongcols = data.frame(valrange= 1:9,
                         colcode = brewer.pal(n=9, name="YlGnBu"))
RegiontotSF$MLclass = as.numeric(as.factor(cut(RegiontotSF$Medlong, breaks=c(0,2,3, 4, 5,6,7,8,10,20))))
RegiontotSF$MLplotcol = medlongcols$colcode [match(RegiontotSF$MLclass, medlongcols$valrange)]
MPAsICESmarine = readRDS(file=paste0(datPath, "MPAsICESmarine.rds")) 
MSFDBHTreg = readRDS(file=paste0(outPath, "MSFDBHTreg.rds")) 

BHTcols = data.frame(BHTcol = c("coral1", "firebrick3", "darkred", #mud 
                                "orange1", "orange3", "orange4", #mud/sand
                                "khaki1", "khaki3", "lightgoldenrod4", #sand 
                                "palegreen1", "palegreen3", "palegreen4", #mixed 
                                "paleturquoise1", "paleturquoise3", "paleturquoise4", #coarse 
                                "steelblue1",  "steelblue3", "steelblue4", #rock/biogenic 
                                "plum1", "plum1", "plum1", #upper bathyal
                                "lightpink1", "lightpink1", "lightpink1", # lower bathyal
                                "mediumpurple4", "dimgrey"),# abyssal
                     MSFD_BBHT = c("Infralittoral mud", "Circalittoral mud", "Offshore circalittoral mud", 
                                   "Infralittoral mud or Infralittoral sand", "Circalittoral mud or Circalittoral sand", "Offshore circalittoral mud or Offshore circalittoral sand",
                                   "Infralittoral sand", "Circalittoral sand", "Offshore circalittoral sand", 
                                   "Infralittoral mixed sediment", "Circalittoral mixed sediment", "Offshore circalittoral mixed sediment", 
                                   "Infralittoral coarse sediment", "Circalittoral coarse sediment", "Offshore circalittoral coarse sediment", 
                                   "Infralittoral rock and biogenic reef", "Circalittoral rock and biogenic reef", "Offshore circalittoral rock and biogenic reef", 
                                   "Upper bathyal sediment", "Upper bathyal sediment or Upper bathyal rock and biogenic reef", "Upper bathyal rock and biogenic reef",
                                   "Lower bathyal sediment", "Lower bathyal sediment or Lower bathyal rock and biogenic reef", "Lower bathyal rock and biogenic reef",
                                   "Abyssal", "Na"))
MSFDBHTreg$BHTcol = BHTcols$BHTcol [match(MSFDBHTreg$MSFD_BBHT, BHTcols$MSFD_BBHT)]

##-- Create Fig 1. Combination of assessment area, MPAs, MSFD, and fishing
tiff(filename=paste0(outPath, "Fig1.tiff"), width = 19, height = 20.3, res=500, units = "cm")
layout(mat=matrix(c(5,5,1,2,6,6,3,4), ncol=2, nrow=4, byrow=T), heights=c(0.02, 0.48, 0.02, 0.48))
par(oma=c(0.5, 0.5, 0.5, 0.5))
par(mar=c(1,1,1,1))
plotRGB(RegMap)
plot(st_geometry(ICESEcoregs), col=c("palegreen3", "coral2", "plum3", "steelblue2"), add=T)
plot(st_geometry(RegiontotSF), col=adjustcolor("dimgrey", alpha.f = 0.5), border=NA, add=T)
plot(st_geometry(ICESEcoregs), border="black", add=T, lwd=1)
legend("bottomright", legend=c("Baltic Sea", "Greater North Sea", "Celtic Seas", "Bay of Biscay &\nthe Iberian Coast", "", "Assessment area"), 
       fill=c("plum3", "steelblue2", "coral2", "palegreen3", "white", "#69696980"), cex=1, title="ICES Ecoregion", inset=0.02, 
       title.font=2, y.intersp = c(rep(0.8, 4), 0.3, 0.8), border=c(rep("black", 4), "white", "black"))
mtext(side=3, line=-0.5, adj=0.02, font=2, "a")
text(x=c(-7.2, -11, 18.7, 3), y=c(45, 56.5, 55.8, 55.6), labels=c("BoBIC", "CS", "BS", "GNS"), font=2)

par(mar=c(1,1,1,1))
plotRGB(RegMap)
plot(st_geometry(subset(MPAsICESmarine, MPA_status %in% c("A", "noHab"))), 
     col=subset(MPAsICESmarine, MPA_status %in% c("A", "noHab"))$plotcol, border=NA, add=T)
plot(st_geometry(subset(MPAsICESmarine, MPA_status %in% c("B", "Hab", "C"))), 
     col=subset(MPAsICESmarine, MPA_status %in% c("B", "Hab", "C"))$plotcol, border=NA, add=T)
plot(st_geometry(RegiontotSF), col=adjustcolor("dimgrey", alpha.f = 0.5), border=NA, add=T)
plot(st_geometry(ICESEcoregs), border="black", add=T, lwd=1)
legend("bottomright", legend=c("Non-habitat protection", "Habitat protection", "", "Assessment area"), border=c("black", "black", "white", "black"),
       fill=c("darkorange", "darkred","white", "#69696980"), cex=1, title="Natura 2000 sites", inset=0.02, title.font=2, y.intersp = c(0.8,0.8,0.3,0.8))
mtext(side=3, line=-0.5, adj=0.02, font=2, "b")

par(mar=c(1,1,1,1))
plotRGB(RegMap)
plot(st_geometry(MSFDBHTreg), col=MSFDBHTreg$BHTcol, border=NA, add=T)
plot(st_geometry(RegiontotSF), col=adjustcolor("dimgrey", alpha.f = 0.2), border=NA, add=T)
plot(st_geometry(ICESEcoregs), border="black", add=T, lwd=1)
l1 = legend("bottomright", fill=c("white", "coral1", "orange1", "khaki1", "palegreen1", "paleturquoise1", "steelblue1", "plum1", "white", "white"), border=c("white", rep("black", 7,), rep("white",2)),
            legend=c("", "Mud", "Mud or Sand", "Sand", "Mixed sediments", "Coarse sediments", "Rock or biogenic reef", "Upper bathyal, lower bathyal\n& Abyssal", "", "Assessment area"), 
            cex=1, title="MSFD Habitat types", inset=0.02, title.font=2, y.intersp = 0.6, x.intersp=c(0, rep(3, 9)), plot=T, text.col = "white", title.col = "black")
l2 = legend("bottomright", fill=c("white", "firebrick3", "orange3", "khaki3", "palegreen3", "paleturquoise3", "steelblue3","lightpink1", "white",  "#69696980"), border=c("white", rep("black", 7, "white", "black")),
            legend = c("", "Mud", "Mud or Sand", "Sand", "Mixed sediments", "Coarse sediments", "Rock or biogenic reef", "Upper bathyal, lower bathyal\n& Abyssal", "", "Assessment area"), 
            cex=1, title="", bty="n", plot=T, y.intersp = 0.6, x.intersp=0, text.col = "white", inset=c(0.08, 0.02))
l3 = legend("bottomright", fill=c("white", "darkred", "orange4", "khaki4", "palegreen4", "paleturquoise4", "steelblue4","mediumpurple4", "white", "white"), border=c("white", rep("black", 7), "white", "white"),
            legend = c("", "Mud", "Mud or Sand", "Sand", "Mixed sediments", "Coarse sediments", "Rock or biogenic reef", "Upper bathyal, lower bathyal\n& abyssal", "", "Assessment area"), cex=1, title="", 
            inset=c(0.042, 0.02), y.intersp = 0.6, bty="n", x.intersp=0.4)
text(x=5.7, y=48.8, "Infr. / Circ. / Off. Circ", pos=4)
mtext(side=3, line=-0.5, adj=0.02, font=2, "c")

par(mar=c(1,1,1,1))
plotRGB(RegMap)
plot(st_geometry(FishDataANN), col=FishDataANN$SARplotcol, border=NA, add=T)
plot(st_geometry(RegiontotSF), col=adjustcolor("dimgrey", alpha.f = 0.2), border=NA, add=T)
plot(st_geometry(ICESEcoregs), border="black", add=T, lwd=1)
l = legend("bottomright", legend=c("0 - 0.1", "0.1 - 0.5",  "0.5 - 1","", "", "1 - 5", "5 - 10", ">10", "", ""), 
           fill=c(plotcols$colcode[2:4], "white", "#69696980", plotcols$colcode[5:7], "white", "white"), border=c(rep("black", 3), "white", rep("black", 4), rep("white",2)), 
           cex=1, title="Fishing intensity", inset=0.02, title.font=2, y.intersp = c(0.8, 0.8, 0.8, 0.3, 0.8), ncol=2, plot=T)
text(x=l$text$x[5], y=l$text$y[5], "Assessment area", pos=4, offset=-0.09)
mtext(side=3, line=-0.5, adj=0.02, font=2, "d")

dev.off()

#-#-# Step 10. Determine spatial overlap ecoregion and assessed area (based on BHT-extent!) ----
RegionMSFD = readRDS(file=paste0(datPath, "RegionMSFD.rds"))
ICESEcoregs$area = st_area(ICESEcoregs)
ICESEcoregs$Ext_sqkm = as.numeric(ICESEcoregs$area)/1E6
a = aggregate(BHTareakm2~ Ecoregion, data=RegionMSFD, FUN="sum")
a$EcoregExt = ICESEcoregs$Ext_sqkm [match(a$Ecoregion, ICESEcoregs$Ecoregion)]
a$overlapPCT = a$BHTareakm2 / a$EcoregExt * 100

# a
# Ecoregion                             BHTareakm2 EcoregExt  overlapPCT
# Baltic Sea                            356544.5   399959.9   89.14506
# Bay of Biscay and the Iberian Coast   176452.1   755807.9   23.34615
# Celtic Seas                           472191.1   917770.8   51.44978
# Greater North Sea                     498816.3   670838.4   74.35715

