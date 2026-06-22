#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
##  Rscript for the analyse underlying the article 
#      "Evaluating the efficacy of closing the current European network of marine protected areas 
#       to mobile bottom-contacting fishing gears"
#
##  Step 3: Determine Status quo of fishing pressure + impact
#
##  Code by Karin van der Reijden. (kjvdreijden@gmail.com)
#        published under the MIT-license
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Step 0. Load libraries, set data paths ----
library(sf)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(grDevices)
library(plyr)
sf_use_s2(FALSE)

##-- set data paths
datPath <- "D:/MS_Natura2000closures/Data/"
outPath <- "D:/MS_Natura2000closures/Routput/"
scriptsPath <- "D:/MS_Natura2000closures/Rscripts/"

#-#-# Step 1. Determine state based on current fishing effort ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))

##-- Determine SAR from the swept area records
VMSdatMet = VMSdatMet[,c(1:102, 165)]
idx = grep("surface_SA", colnames(VMSdatMet))
SAdat = VMSdatMet[,..idx]
SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN = function(x){ x/VMSdatMet$area_sqkm}))
names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
idy = setdiff(1:ncol(VMSdatMet), idx)
VMSdatMet = cbind(VMSdatMet[,..idy], SAdat2)

source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo (https://github.com/ices-eg/FBIT)
Period = 2016:2020

source(paste0(scriptsPath, "Habitatstatefishing_inf.R")) # adjusted from FBIT repo (https://github.com/ices-eg/FBIT)
saveRDS(State_reg, paste0(outPath, "State_reg_S0_inf.rds"))

source(paste0(scriptsPath, "Habitatstatefishing_epi.R")) # adjusted from FBIT repo (https://github.com/ices-eg/FBIT)
saveRDS(State_reg, paste0(outPath, "State_reg_S0_epi.rds"))

#-#-# Step 2. Determine average PD per BHT & ecoregion ----
##-- Load data
RegionMSFD <- as.data.frame(readRDS(paste0(datPath, "RegionMSFD2.rds")))
State_regINF   = readRDS(paste0(outPath, "State_reg_S0_inf.rds"))
State_regEPI   = readRDS(paste0(outPath, "State_reg_S0_epi.rds"))

##-- Assign annual and average PD for both inf and epi.
RegionMSFD[,"PD_inf_S0_2016"] = State_regINF$state_2016 [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]
RegionMSFD[,"PD_inf_S0_2017"] = State_regINF$state_2017 [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]
RegionMSFD[,"PD_inf_S0_2018"] = State_regINF$state_2018 [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]
RegionMSFD[,"PD_inf_S0_2019"] = State_regINF$state_2019 [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]
RegionMSFD[,"PD_inf_S0_2020"] = State_regINF$state_2020 [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]
RegionMSFD[,"PD_inf_S0_1620"] = State_regINF$AV1620_PD [match(RegionMSFD$csquares, State_regINF$FisheriesMet.csquares)]

RegionMSFD[,"PD_epi_S0_2016"] = State_regEPI$state_2016 [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]
RegionMSFD[,"PD_epi_S0_2017"] = State_regEPI$state_2017 [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]
RegionMSFD[,"PD_epi_S0_2018"] = State_regEPI$state_2018 [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]
RegionMSFD[,"PD_epi_S0_2019"] = State_regEPI$state_2019 [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]
RegionMSFD[,"PD_epi_S0_2020"] = State_regEPI$state_2020 [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]
RegionMSFD[,"PD_epi_S0_1620"] = State_regEPI$AV1620_PD [match(RegionMSFD$csquares, State_regEPI$FisheriesMet.csquares)]

##-- Determine mean PD for epi and inf combined.
RegionMSFD$PD_comb_S0_2016 = rowMeans(RegionMSFD[,c(11, 17)], na.rm=T)
RegionMSFD$PD_comb_S0_2017 = rowMeans(RegionMSFD[,c(12, 18)], na.rm=T)
RegionMSFD$PD_comb_S0_2018 = rowMeans(RegionMSFD[,c(13, 19)], na.rm=T)
RegionMSFD$PD_comb_S0_2019 = rowMeans(RegionMSFD[,c(14, 20)], na.rm=T)
RegionMSFD$PD_comb_S0_2020 = rowMeans(RegionMSFD[,c(15, 21)], na.rm=T)
RegionMSFD$PD_comb_S0_1620 = rowMeans(RegionMSFD[,c(16, 22)], na.rm=T)

##-- Remove the epi and inf columns
epi = grep("epi", colnames(RegionMSFD))
RegionMSFD = RegionMSFD[,setdiff(1:ncol(RegionMSFD), epi)]
inf = grep("inf", colnames(RegionMSFD))
RegionMSFD = RegionMSFD[,setdiff(1:ncol(RegionMSFD), inf)]

##-- Save results
saveRDS(RegionMSFD, file=paste0(datPath, "RegionMSFD3.rds"))

#-#-# Step 3. Determine average PD per ecoregion ----
##-- Load data
Regiontot <- as.data.frame(readRDS(paste0(datPath, "Regiontot2.rds")))
State_regINF   = readRDS(paste0(outPath, "State_reg_S0_inf.rds"))
State_regEPI   = readRDS(paste0(outPath, "State_reg_S0_epi.rds"))

##-- Assign annual and mean PD for epi and inf
Regiontot[,"PD_inf_S0_2016"] = State_regINF$state_2016 [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]
Regiontot[,"PD_inf_S0_2017"] = State_regINF$state_2017 [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]
Regiontot[,"PD_inf_S0_2018"] = State_regINF$state_2018 [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]
Regiontot[,"PD_inf_S0_2019"] = State_regINF$state_2019 [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]
Regiontot[,"PD_inf_S0_2020"] = State_regINF$state_2020 [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]
Regiontot[,"PD_inf_S0_1620"] = State_regINF$AV1620_PD [match(Regiontot$csquares, State_regINF$FisheriesMet.csquares)]

Regiontot[,"PD_epi_S0_2016"] = State_regEPI$state_2016 [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]
Regiontot[,"PD_epi_S0_2017"] = State_regEPI$state_2017 [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]
Regiontot[,"PD_epi_S0_2018"] = State_regEPI$state_2018 [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]
Regiontot[,"PD_epi_S0_2019"] = State_regEPI$state_2019 [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]
Regiontot[,"PD_epi_S0_2020"] = State_regEPI$state_2020 [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]
Regiontot[,"PD_epi_S0_1620"] = State_regEPI$AV1620_PD [match(Regiontot$csquares, State_regEPI$FisheriesMet.csquares)]

##-- Determine annual and mean PD for epi and inf combined.
Regiontot$PD_comb_S0_2016 = rowMeans(Regiontot[,c(18, 24)], na.rm=T)
Regiontot$PD_comb_S0_2017 = rowMeans(Regiontot[,c(19, 25)], na.rm=T)
Regiontot$PD_comb_S0_2018 = rowMeans(Regiontot[,c(20, 26)], na.rm=T)
Regiontot$PD_comb_S0_2019 = rowMeans(Regiontot[,c(21, 27)], na.rm=T)
Regiontot$PD_comb_S0_2020 = rowMeans(Regiontot[,c(22, 28)], na.rm=T)
Regiontot$PD_comb_S0_1620 = rowMeans(Regiontot[,c(23, 29)], na.rm=T)

##-- Remove epi and inf data
epi = grep("epi", colnames(Regiontot))
Regiontot = Regiontot[,setdiff(1:ncol(Regiontot),epi)]
inf = grep("inf", colnames(Regiontot))
Regiontot = Regiontot[,setdiff(1:ncol(Regiontot),inf)]

##-- Save results
saveRDS(Regiontot, file=paste0(datPath, "Regiontot3.rds"))

#-#-# Step 4. Calculate quality-extent threshold trade-offs - per BHT per ecoregion ----
QET_BHT_table <- data.table()
RegionMSFD <- readRDS(paste0(datPath, "RegionMSFD3.rds"))
RegionMSFD <- subset(RegionMSFD, !MSFD_BBHT=="Na")

##-- Determine 5 most abundant BHT per Ecoregion
sub = subset(RegionMSFD, is.na(RegionMSFD$PD_comb_S0_1620)==F)
BHTssub <- aggregate(BHTareakm2~MSFD_BBHT+Ecoregion, data=sub, FUN=sum)
BHTssub <- BHTssub[order(BHTssub$Ecoregion, -BHTssub$BHTareakm2),]
BHTs = data.frame(Ecoregion = rep(c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast"), each=5),
                  MSFD_BBHT = c("Offshore circalittoral sand",
                                "Circalittoral sand", 
                                "Offshore circalittoral mud",
                                "Offshore circalittoral coarse sediment",
                                "Circalittoral coarse sediment", # end GNS
                                # start CS
                                "Offshore circalittoral sand", 
                                "Offshore circalittoral coarse sediment",
                                "Upper bathyal sediment",
                                "Offshore circalittoral mud",
                                "Circalittoral coarse sediment", # end CS
                                # start BS
                                "Circalittoral sand", "Infralittoral sand", 
                                "Offshore circalittoral mud",
                                "Offshore circalittoral mixed sediment",
                                "Offshore circalittoral mud or Offshore circalittoral sand", # end BS
                                # start BBIC
                                "Offshore circalittoral sand",
                                "Offshore circalittoral mud",
                                "Upper bathyal sediment",
                                "Offshore circalittoral coarse sediment",
                                "Circalittoral sand"))

##-- Create table to save PD values per ecoregion for 5 most abundant BHT
BHTReg <- data.table()
for(iReg in c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast")){
  BHTsd <- subset(BHTs, Ecoregion==iReg)
  BHTsd$ID <- 1:nrow(BHTsd)
  BHTReg <- rbind(BHTReg, BHTsd[1:5,])
}
iInd = "PD" ; iBtype = "comb"
for(iReg in c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast")){
  for(iYear in c(2016:2020, 1620)){
    MA_BHT <- subset(BHTReg, Ecoregion == iReg)$MSFD_BBHT
    subRegMSFD <- subset(RegionMSFD, Ecoregion == iReg)
    subRegMSFD2 = subset(subRegMSFD, is.na(PD_comb_S0_1620)==F)
    
    table1 <- data.frame(MSFD = rep(NA, 5), ID =1:5, Indicator = iInd, Ecoregion = iReg, Biotype = iBtype, Period = iYear)
    colID <- which(colnames(subRegMSFD2) == paste(iInd, iBtype, "S0", iYear, sep="_"))
    colnames(subRegMSFD2)[colID] <- "ColofInterest"
    
    for(iBHT in 1:length(MA_BHT)){
      table1$MSFD[iBHT] = MA_BHT[iBHT]
      Unit.data = subset(subRegMSFD2, MSFD_BBHT == MA_BHT[iBHT])
      if(nrow(Unit.data)>0){
        Unit.data = Unit.data[order(Unit.data$ColofInterest, -Unit.data$BHTareakm2),]
        TotBHText = sum(Unit.data$BHTareakm2)
        table1$q.100[iBHT] = sum(subset(Unit.data, ColofInterest >= 1)$BHTareakm2)/TotBHText
        table1$q.99[iBHT] = sum(subset(Unit.data, ColofInterest >= .99)$BHTareakm2)/TotBHText
        table1$q.95[iBHT] = sum(subset(Unit.data, ColofInterest >= .95)$BHTareakm2)/TotBHText
        table1$q.90[iBHT] = sum(subset(Unit.data, ColofInterest >= .90)$BHTareakm2)/TotBHText
        table1$q.85[iBHT] = sum(subset(Unit.data, ColofInterest >= .85)$BHTareakm2)/TotBHText
        table1$q.80[iBHT] = sum(subset(Unit.data, ColofInterest >= .80)$BHTareakm2)/TotBHText
        table1$q.75[iBHT] = sum(subset(Unit.data, ColofInterest >= .75)$BHTareakm2)/TotBHText
        table1$q.70[iBHT] = sum(subset(Unit.data, ColofInterest >= .70)$BHTareakm2)/TotBHText
        table1$q.65[iBHT] = sum(subset(Unit.data, ColofInterest >= .65)$BHTareakm2)/TotBHText
        table1$q.60[iBHT] = sum(subset(Unit.data, ColofInterest >= .60)$BHTareakm2)/TotBHText
        table1$q.55[iBHT] = sum(subset(Unit.data, ColofInterest >= .55)$BHTareakm2)/TotBHText
        table1$q.50[iBHT] = sum(subset(Unit.data, ColofInterest >= .50)$BHTareakm2)/TotBHText
        table1$q.45[iBHT] = sum(subset(Unit.data, ColofInterest >= .45)$BHTareakm2)/TotBHText
        table1$q.40[iBHT] = sum(subset(Unit.data, ColofInterest >= .40)$BHTareakm2)/TotBHText
        table1$q.35[iBHT] = sum(subset(Unit.data, ColofInterest >= .35)$BHTareakm2)/TotBHText
        table1$q.30[iBHT] = sum(subset(Unit.data, ColofInterest >= .30)$BHTareakm2)/TotBHText
        table1$q.25[iBHT] = sum(subset(Unit.data, ColofInterest >= .25)$BHTareakm2)/TotBHText
        table1$q.20[iBHT] = sum(subset(Unit.data, ColofInterest >= .20)$BHTareakm2)/TotBHText
        table1$q.15[iBHT] = sum(subset(Unit.data, ColofInterest >= .15)$BHTareakm2)/TotBHText
        table1$q.10[iBHT] = sum(subset(Unit.data, ColofInterest >= .10)$BHTareakm2)/TotBHText
        table1$q.05[iBHT] = sum(subset(Unit.data, ColofInterest >= .05)$BHTareakm2)/TotBHText
        table1$q.00[iBHT] = sum(subset(Unit.data, ColofInterest >= .00)$BHTareakm2)/TotBHText
      } else {
        table1$q.100[iBHT] = as.numeric(NA)
        table1$q.99[iBHT] = as.numeric(NA)
        table1$q.95[iBHT] = as.numeric(NA)
        table1$q.90[iBHT] = as.numeric(NA)
        table1$q.85[iBHT] = as.numeric(NA)
        table1$q.80[iBHT] = as.numeric(NA)
        table1$q.75[iBHT] = as.numeric(NA)
        table1$q.70[iBHT] = as.numeric(NA)
        table1$q.65[iBHT] = as.numeric(NA)
        table1$q.60[iBHT] = as.numeric(NA)
        table1$q.55[iBHT] = as.numeric(NA)
        table1$q.50[iBHT] = as.numeric(NA)
        table1$q.45[iBHT] = as.numeric(NA)
        table1$q.40[iBHT] = as.numeric(NA)
        table1$q.35[iBHT] = as.numeric(NA)
        table1$q.30[iBHT] = as.numeric(NA)
        table1$q.25[iBHT] = as.numeric(NA)
        table1$q.20[iBHT] = as.numeric(NA)
        table1$q.15[iBHT] = as.numeric(NA)
        table1$q.10[iBHT] = as.numeric(NA)
        table1$q.05[iBHT] = as.numeric(NA)
        table1$q.00[iBHT] = as.numeric(NA)
      }
    } # end iBHT
    
    ##-- Formatting table to nice output  
    table2 <- gather(table1, key = "quality", value = "extent",
                     q.100, q.99, q.95, q.90, q.85, q.80, q.75, q.70, q.65, q.60, q.55, q.50, q.45, q.40, q.35, q.30, q.25, q.20, q.15, q.10, q.05, q.00)
    table2 <- table2 %>%
      mutate(quality = recode(quality, "q.100" = "1.00",
                              "q.99" = "0.99", "q.95" = "0.95", "q.90" = "0.90", 
                              "q.85" = "0.85", "q.80" = "0.80", "q.75" = "0.75", 
                              "q.70" = "0.70", "q.65" = "0.65", "q.60" = "0.60", 
                              "q.55" = "0.55", "q.50" = "0.50", "q.45" = "0.45", 
                              "q.40" = "0.40", "q.35" = "0.35", "q.30" = "0.30", 
                              "q.25" = "0.25", "q.20" = "0.20", "q.15" = "0.15", 
                              "q.10" = "0.10", "q.05" = "0.05", "q.00" = "0.00"))
    
    table2$extent <- round(table2$extent, digits = 3)
    table2 <- table2[order(table2$ID),]

    ##-- Store results in QET_BHT_table
    QET_BHT_table <- rbind(QET_BHT_table, table2)
  } # end iYear
} # end iReg

saveRDS(QET_BHT_table, file=paste0(outPath, "QET_BHT_table_S0.rds"))

#-#-# Step 5. Calculate quality-extent threshold trade-offs - per ecoregion ----
QET_table <- data.table()
Regiontot <- readRDS(paste0(datPath, "Regiontot3.rds"))
iInd = "PD" ; iBtype = "comb"
for(iReg in c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast")){
  for(iYear in c(2016:2020, 1620)){
      subReg <- subset(Regiontot, Ecoregion == iReg)
      subReg2 = subset(subReg, is.na(PD_comb_S0_1620)==F)
      
      table1 <- data.frame(Indicator = iInd, Ecoregion = iReg, Biotype = iBtype, Period = iYear)
      colID  <- which(colnames(subReg2) == paste(iInd, iBtype, "S0", iYear, sep="_"))
      colnames(subReg2)[colID] <- "ColofInterest"
      
      Unit.data = subReg2
      Unit.data = Unit.data[order(Unit.data$ColofInterest, -Unit.data$area_sqkm),]
      TotExt = sum(Unit.data$area_sqkm)
      table1$q.100= sum(subset(Unit.data, ColofInterest >= 1)$area_sqkm)/TotExt
      table1$q.99 = sum(subset(Unit.data, ColofInterest >= .99)$area_sqkm)/TotExt
      table1$q.95 = sum(subset(Unit.data, ColofInterest >= .95)$area_sqkm)/TotExt
      table1$q.90 = sum(subset(Unit.data, ColofInterest >= .90)$area_sqkm)/TotExt
      table1$q.85 = sum(subset(Unit.data, ColofInterest >= .85)$area_sqkm)/TotExt
      table1$q.80 = sum(subset(Unit.data, ColofInterest >= .80)$area_sqkm)/TotExt
      table1$q.75 = sum(subset(Unit.data, ColofInterest >= .75)$area_sqkm)/TotExt
      table1$q.70 = sum(subset(Unit.data, ColofInterest >= .70)$area_sqkm)/TotExt
      table1$q.65 = sum(subset(Unit.data, ColofInterest >= .65)$area_sqkm)/TotExt
      table1$q.60 = sum(subset(Unit.data, ColofInterest >= .60)$area_sqkm)/TotExt
      table1$q.55 = sum(subset(Unit.data, ColofInterest >= .55)$area_sqkm)/TotExt
      table1$q.50 = sum(subset(Unit.data, ColofInterest >= .50)$area_sqkm)/TotExt
      table1$q.45 = sum(subset(Unit.data, ColofInterest >= .45)$area_sqkm)/TotExt
      table1$q.40 = sum(subset(Unit.data, ColofInterest >= .40)$area_sqkm)/TotExt
      table1$q.35 = sum(subset(Unit.data, ColofInterest >= .35)$area_sqkm)/TotExt
      table1$q.30 = sum(subset(Unit.data, ColofInterest >= .30)$area_sqkm)/TotExt
      table1$q.25 = sum(subset(Unit.data, ColofInterest >= .25)$area_sqkm)/TotExt
      table1$q.20 = sum(subset(Unit.data, ColofInterest >= .20)$area_sqkm)/TotExt
      table1$q.15 = sum(subset(Unit.data, ColofInterest >= .15)$area_sqkm)/TotExt
      table1$q.10 = sum(subset(Unit.data, ColofInterest >= .10)$area_sqkm)/TotExt
      table1$q.05 = sum(subset(Unit.data, ColofInterest >= .05)$area_sqkm)/TotExt
      table1$q.00 = sum(subset(Unit.data, ColofInterest >= .00)$area_sqkm)/TotExt
      
      ##-- Formatting table to nice output  
      table2 <- gather(table1, key = "quality", value = "extent",
                       q.100, q.99, q.95, q.90, q.85, q.80, q.75, q.70, q.65, q.60, q.55, q.50, q.45, q.40, q.35, q.30, q.25, q.20, q.15, q.10, q.05, q.00)
      table2 <- table2 %>%
        mutate(quality = recode(quality, "q.100" = "1.00",
                                "q.99" = "0.99", "q.95" = "0.95", "q.90" = "0.90", 
                                "q.85" = "0.85", "q.80" = "0.80", "q.75" = "0.75", 
                                "q.70" = "0.70", "q.65" = "0.65", "q.60" = "0.60", 
                                "q.55" = "0.55", "q.50" = "0.50", "q.45" = "0.45", 
                                "q.40" = "0.40", "q.35" = "0.35", "q.30" = "0.30", 
                                "q.25" = "0.25", "q.20" = "0.20", "q.15" = "0.15", 
                                "q.10" = "0.10", "q.05" = "0.05", "q.00" = "0.00"))
      
      table2$extent <- round(table2$extent, digits = 3)

      ##-- Store results in QET_BHT_table
      QET_table <- rbind(QET_table, table2)
  } # end iYear
} # end iReg

saveRDS(QET_table, file=paste0(outPath, "QET_table_S0.rds"))

#-#-# Step 6. Determine overlap BHT and MPAs per ecoregion (Table S1) ----
RegionMSFD              <- as.data.table(readRDS(file=paste0(datPath, "RegionMSFD3.rds")))
RegionMSFD$MSFD_BBHT = ifelse(RegionMSFD$MSFD_BBHT == "Na", "Unknown", RegionMSFD$MSFD_BBHT)
tabdat <- aggregate(cbind(S1area, S2area, BHTareakm2)~Ecoregion+MSFD_BBHT, data=RegionMSFD, FUN=sum)
tabdattots <- aggregate(cbind(S1area, S2area, BHTareakm2)~Ecoregion, data=RegionMSFD, FUN=sum)
tabdat2 <- aggregate(cbind(S1area, S2area, BHTareakm2)~1, data=RegionMSFD, FUN=sum)

regs = data.table(Reg = c("Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast", "Baltic Sea"),
                  RegCode = c("GNS", "CS", "BBIC", "BS"))

for(iReg in c("Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast", "Baltic Sea")) {
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
                              "Abyssal", "Unknown"), 
                     habOrd = 1:26, 
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
                                  "Abyssal", "Unknown"))
  subdattab <- subset(tabdat, Ecoregion== iReg)
  subdattab <- subdattab[,-1]
  subdattab$S1 <- (subdattab$S1area - subdattab$S2area) / subdattab$BHTareakm2 * 100
  subdattab$S2 <- subdattab$S2area / subdattab$BHTareakm2 * 100
  subdattab$noMPA <- (subdattab$BHTareakm2 - subdattab$S1area) / subdattab$BHTareakm2 * 100
  subdattab <- subdattab[,c(1,6:5, 7, 4)]
  names(subdattab)[1] <- c("MSFD")
  TabH <- merge(TabH, subdattab, by="MSFD", all.x=T)
  TabH <- TabH[order(TabH$habOrd, decreasing = T)]
  TabH$label <- ifelse(TabH$BHTareakm2 < 1000, paste0(round(TabH$BHTareakm2)),
                       ifelse(TabH$BHTareakm2 >= 1000 & TabH$BHTareakm2 < 10000, paste0(substr(round(TabH$BHTareakm2),1,1), ".", substr(round(TabH$BHTareakm2),2,nchar(round(TabH$BHTareakm2)))),
                              ifelse(TabH$BHTareakm2 >= 10000 & TabH$BHTareakm2 < 100000, paste0(substr(round(TabH$BHTareakm2),1,2), ".", substr(round(TabH$BHTareakm2),3,nchar(round(TabH$BHTareakm2)))), 
                                     ifelse(TabH$BHTareakm2 >= 100000 & TabH$BHTareakm2 < 1000000, paste0(substr(round(TabH$BHTareakm2),1,3), ".", substr(round(TabH$BHTareakm2),4,nchar(round(TabH$BHTareakm2)))), paste0("XX")))))
  Htab <- t(TabH[,4:6])
  } # end iReg

##-- Create Table S1: MPA and HDS coverage per MSFD-BHT, per region
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
                            "Abyssal", "Unknown", "Total"), 
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

tabdattots$MSFD_BBHT = "Total"
tabdat = rbind(tabdat, tabdattots)

a = subset(tabdat, Ecoregion == "Baltic Sea")
totext = subset(a, MSFD_BBHT == "Total")$BHTareakm2
TabH$BS_TOT = round(a$BHTareakm2/1000, digits=2) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BS_TOTp = round(a$BHTareakm2 / totext *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BS_ALL = round((a$S1area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BS_HDS= round((a$S2area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]

a = subset(tabdat, Ecoregion == "Greater North Sea")
totext = subset(a, MSFD_BBHT == "Total")$BHTareakm2
TabH$GNS_TOT = round(a$BHTareakm2/1000, digits=2) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$GNS_TOTp = round(a$BHTareakm2 / totext *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$GNS_ALL = round((a$S1area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$GNS_HDS= round((a$S2area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]

a = subset(tabdat, Ecoregion == "Celtic Seas")
totext = subset(a, MSFD_BBHT == "Total")$BHTareakm2
TabH$CS_TOT = round(a$BHTareakm2/1000, digits=2) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$CS_TOTp = round(a$BHTareakm2 / totext *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$CS_ALL = round((a$S1area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$CS_HDS= round((a$S2area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]

a = subset(tabdat, Ecoregion == "Bay of Biscay and the Iberian Coast")
totext = subset(a, MSFD_BBHT == "Total")$BHTareakm2
TabH$BoBIC_TOT = round(a$BHTareakm2/1000, digits=2) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BoBIC_TOTp = round(a$BHTareakm2 / totext *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BoBIC_ALL = round((a$S1area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]
TabH$BoBIC_HDS= round((a$S2area/a$BHTareakm2) *100, digits=1) [match(TabH$MSFD, a$MSFD_BBHT)]

TableS1 = TabH[,3:19]

write.table(TableS1, file=paste0(outPath, "TableS1.txt"), sep="\t", row.names = F)


#-#-# Step 7. Spatial overlap BHT and MPA + fisheries (Figure 2) ----
##-- Determine spatial overlap BHT-MPA per regio
RegionMSFD              <- as.data.table(readRDS(file=paste0(datPath, "RegionMSFD3.rds")))
RegionMSFD$Sediment = ifelse(RegionMSFD$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Circalittoral rock and biogenic reef",
                                                         "Offshore circalittoral rock and biogenic reef",
                                                         "Upper bathyal rock and biogenic reef", "Lower bathyal rock and biogenic reef",
                                                         "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                                         "Lower bathyal sediment or Lower bathyal rock and biogenic reef"), "Gravel",
                             ifelse(RegionMSFD$MSFD_BBHT %in% c("Infralittoral sand",
                                                                "Circalittoral sand", 
                                                                "Offshore circalittoral sand"), "Sand", 
                                    ifelse(RegionMSFD$MSFD_BBHT %in% c("Infralittoral mud", "Infralittoral mud or Infralittoral sand",
                                                                       "Circalittoral mud", "Circalittoral mud or Circalittoral sand",
                                                                       "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand",
                                                                       "Upper bathyal sediment", "Lower bathyal sediment", "Abyssal"), "Mud", 
                                           ifelse(RegionMSFD$MSFD_BBHT %in% c("Infralittoral mixed sediment", "Circalittoral mixed sediment", "Offshore circalittoral mixed sediment"), "Mixed", 
                                                  ifelse(RegionMSFD$MSFD_BBHT %in% c("Circalittoral coarse sediment", "Infralittoral coarse sediment", "Offshore circalittoral coarse sediment"), "Coarse", "Unknown")))))
RegionMSFD$Depthclass   = ifelse(RegionMSFD$MSFD_BBHT %in% c("Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
                                                             "Infralittoral mixed sediment", "Infralittoral sand", 
                                                             "Infralittoral mud", "Infralittoral mud or Infralittoral sand"), "Infralittoral", 
                                 ifelse(RegionMSFD$MSFD_BBHT %in% c("Circalittoral rock and biogenic reef", "Circalittoral coarse sediment",
                                                                    "Circalittoral mixed sediment", "Circalittoral sand", 
                                                                    "Circalittoral mud", "Circalittoral mud or Circalittoral sand"), "Circalittoral", 
                                        ifelse(RegionMSFD$MSFD_BBHT %in% c("Offshore circalittoral rock and biogenic reef", "Offshore circalittoral coarse sediment",
                                                                           "Offshore circalittoral mixed sediment", "Offshore circalittoral sand", 
                                                                           "Offshore circalittoral mud", "Offshore circalittoral mud or Offshore circalittoral sand"), "Offshore circalittoral", 
                                               ifelse(RegionMSFD$MSFD_BBHT %in% c("Upper bathyal rock and biogenic reef", "Upper bathyal sediment", 
                                                                                  "Upper bathyal sediment or Upper bathyal rock and biogenic reef",
                                                                                  "Lower bathyal rock and biogenic reef", "Lower bathyal sediment", 
                                                                                  "Lower bathyal sediment or Lower bathyal rock and biogenic reef",
                                                                                  "Abyssal"), "Bathyal", "Unknown"))))
tabdatD <- aggregate(cbind(S1area, S2area, BHTareakm2)~Depthclass+Ecoregion, data=RegionMSFD, FUN=sum)
tabdatS <- aggregate(cbind(S1area, S2area, BHTareakm2)~Sediment+Ecoregion, data=RegionMSFD, FUN=sum)
tabdatT <- aggregate(cbind(S1area, S2area, BHTareakm2)~Ecoregion, data=RegionMSFD, FUN=sum)

## Determine proportions
dats = c("tabdatD", "tabdatS")
for(iDat in dats){
  ds = get(iDat)
  names(ds) = c("var", "Reg", "MPA", "HDS", "TOT")
  ds$HDSfrac = ds$HDS / ds$TOT * 100
  ds$MPAfrac = ds$MPA / ds$TOT * 100
  ds$totext = tabdatT$BHTareakm2 [match(ds$Reg, tabdatT$Ecoregion)]
  ds$Extfrac = ds$TOT / ds$totext * 100
  ds$HDSratio = ds$HDSfrac / ds$Extfrac
  ds$MPAratio = ds$MPAfrac / ds$Extfrac
  assign(iDat, ds)
}

tabdatD$y = ifelse(tabdatD$var == "Infralittoral", 4, 
                   ifelse(tabdatD$var == "Circalittoral", 3,
                          ifelse(tabdatD$var == "Offshore circalittoral", 2, 1)))
tabdatS$y = ifelse(tabdatS$var == "Gravel", 5,
                   ifelse(tabdatS$var == "Coarse", 4,
                          ifelse(tabdatS$var == "Mixed", 3,
                                 ifelse(tabdatS$var == "Sand", 2, 1))))

##-- Determine spatial overlap fisheries per BHT per regio
FishDataANN = readRDS(file=paste0(datPath, "FD_capped.rds"))
RegionMSFD$AV_SA = FishDataANN$AV_SA1620 [match(RegionMSFD$csquares, FishDataANN$csquares)]
RegionMSFD$BHT_SA = RegionMSFD$AV_SA * (RegionMSFD$BHTareakm2 / RegionMSFD$CSareaNEW)

sedass = aggregate(BHT_SA~ Sediment+Ecoregion, data=RegionMSFD, FUN="sum")
dptass = aggregate(BHT_SA~ Depthclass+Ecoregion, data=RegionMSFD, FUN="sum")
totass = aggregate(BHT_SA ~ Ecoregion, data=RegionMSFD, FUN="sum")
sedass$TOT_BHT = totass$BHT_SA[match(sedass$Ecoregion, totass$Ecoregion)]
dptass$TOT_BHT = totass$BHT_SA[match(dptass$Ecoregion, totass$Ecoregion)]
sedass$PropSA_sed = round(sedass$BHT_SA / sedass$TOT_BHT * 100, digits=2)
dptass$PropSA_dpt = round(dptass$BHT_SA / dptass$TOT_BHT * 100, digits=2)

sedEass = aggregate(BHTareakm2~Sediment+Ecoregion, data=RegionMSFD, FUN="sum")
dptEass = aggregate(BHTareakm2~Depthclass+Ecoregion, data=RegionMSFD, FUN="sum")
totEass = aggregate(BHTareakm2~Ecoregion, data=RegionMSFD, FUN="sum")
sedEass$TOTe = totEass$BHTareakm2 [match(sedEass$Ecoregion, totEass$Ecoregion)]
dptEass$TOTe = totEass$BHTareakm2 [match(dptEass$Ecoregion, totEass$Ecoregion)]
sedEass$propE = round(sedEass$BHTareakm2 / sedEass$TOTe * 100, digits=2)
dptEass$propE = round(dptEass$BHTareakm2 / dptEass$TOTe * 100, digits=2)
dptass$totprop = dptEass$propE
sedass$totprop = sedEass$propE
dptass$FishRatio = dptass$PropSA_dpt / dptass$totprop
sedass$FishRatio = sedass$PropSA_sed / sedass$totprop

sedass$y = ifelse(sedass$Sediment == "Gravel", 5,
                  ifelse(sedass$Sediment == "Coarse", 4,
                         ifelse(sedass$Sediment == "Mixed", 3,
                                ifelse(sedass$Sediment == "Sand", 2, 1))))

dptass$y = ifelse(dptass$Depthclass == "Infralittoral", 4, 
                  ifelse(dptass$Depthclass == "Circalittoral", 3,
                         ifelse(dptass$Depthclass == "Offshore circalittoral", 2, 1)))

##-- Combine gears per region (all) to dataset
tabdatALLD$Reg = "All comb" ; tabdatALLS$Reg = "All comb"
dptassALL$Reg = "All comb" ; sedassALL$Reg = "All comb"
dptassALL$FishRatio = dptassALL$FishProp
sedassALL$FishRatio = sedassALL$FishProp
dptass$Reg = dptass$Ecoregion ; sedass$Reg = sedass$Ecoregion
tabdatALLD$Extfrac= tabdatALLD$Extent ; tabdatALLS$Extfrac =tabdatALLS$Extent

TabDepth = rbind(tabdatALLD[, c("var", "Reg", "Extfrac", "MPAratio", "HDSratio")],
                 tabdatD[,c("var", "Reg", "Extfrac", "MPAratio", "HDSratio")], 
                 data.table(var="Bathyal & Abyssal", Reg = "Baltic Sea",
                            Extfrac= 0, MPAratio=NA, HDSratio=NA))
TabSedim = rbind(tabdatALLS[, c("var", "Reg", "Extfrac", "MPAratio", "HDSratio")],
                 tabdatS[,c("var", "Reg", "Extfrac", "MPAratio", "HDSratio")])
dptTable = rbind(dptassALL[,c("Depthclass", "Reg", "totprop", "FishRatio")],
                 dptass[,c("Depthclass", "Reg", "totprop", "FishRatio")],
                 data.table(Depthclass = "Bathyal & Abyssal", Reg = "Baltic Sea",
                            totprop= 0, FishRatio=NA))
sedTable = rbind(sedassALL[,c("Sediment", "Reg", "totprop", "FishRatio")],
                 sedass[,c("Sediment", "Reg", "totprop", "FishRatio")])

TabDepth$y = c(26.5, 19.5, 12.5, 5.5,   18, 25, 11, NA,   1, 15, 22, 8, NA,  2, 16, 23, 9, NA,  3, 17, 24, 10, NA,  4)
TabDepth$MPAratio = ifelse(TabDepth$Extfrac < 2, NA, TabDepth$MPAratio)
TabDepth$HDSratio = ifelse(TabDepth$Extfrac < 2, NA, TabDepth$HDSratio)
dptTable$FishRatio = ifelse(dptTable$totprop < 2, NA, dptTable$FishRatio)
TabDepth$ylab = c(rep("All", 4), rep("BS", 4), rep("BoBIC", 5), rep("CS", 5), rep("GNS", 5), "BS")
dptTable$y = c(26.5, 19.5, 12.5, 5.5,   18, 25, 11, NA,   1, 15, 22, 8, NA,  2, 16, 23, 9, NA,  3, 17, 24, 10, NA,  4)

TabSedim$MPAratio = ifelse(TabSedim$Extfrac <2, NA, TabSedim$MPAratio)
TabSedim$HDSratio = ifelse(TabSedim$Extfrac <2, NA, TabSedim$HDSratio)
sedTable$FishRatio = ifelse(sedTable$totprop < 2, NA, sedTable$FishRatio)
TabSedim$ylab = c(rep("All", 5), rep("BS", 6), rep("BoBIC", 6), rep("CS", 6), rep("GNS", 6))
TabSedim$y = c(33.5, 26.5, 19.5, 12.5, 5.5,  25, 32, 18, 4, 11, NA,    
               22, 29, 15, 1, 8, NA,   23, 30, 16, 2, 9, NA,  24, 31, 17, 3, 10, NA)
sedTable$y = c(33.5, 26.5, 19.5, 12.5, 5.5,  25, 32, 18, 4, 11, NA,    
               22, 29, 15, 1, 8, NA,   23, 30, 16, 2, 9, NA,  24, 31, 17, 3, 10, NA)

##-- Create figure 2.
tiff(filename=paste0(outPath, "Figure2.tiff"),
     width=17, height=10, units="cm", res=500)
layout(mat=matrix(data=1:2, nrow=1, ncol=2))

## plot MPA and fisheries coverage per depthclass
par(mar=c(4,5,2,1))
plot(x=log(TabDepth$MPAratio), y=TabDepth$y, pch=16, col="grey20", axes=F, ann=F, cex=1)
abline(h=c(20.75, 13.75, 6.75), lwd=1.5, lty=3, col="grey80")
abline(v=0, col="red", lty=2)
points(x=log(TabDepth$HDSratio), y=TabDepth$y, pch=16, col=adjustcolor("grey60", alpha.f = 0.6), cex=1)
points(x=log(dptTable$FishRatio), y=dptTable$y, pch=18, col="cornflowerblue", cex=1)
box()
axis(1, cex.axis=0.8, tck=-0.01, labels=F)
axis(1, tick=F, line=-1, cex.axis=0.8)
axis(2, at=TabDepth$y, tick=F, labels=TabDepth$ylab, cex.axis=0.45, las=1, line=-0.75)
axis(2, at=TabDepth$y, tck=-0.01, labels=F)

mtext(side=2, at=c(5.5, 12.5, 19.5, 26.5), cex=0.6, line=4, las=1, adj=0,
      text=c("Bathyal", "Offshore\ncircalittoral", "Circalittoral", "Infralittoral"))
mtext("a", font=2, side=3, adj=0, cex=0.8)

par(mar=c(4,5,2,1))
plot(log(TabSedim$MPAratio), y=TabSedim$y, pch=16, col="grey20", axes=F, ann=F, xlim=c(-2.1, 2.8))
abline(h=c(27.75, 20.75, 13.75, 6.75), lwd=1.5, lty=3, col="grey80")
abline(v=0, col="red", lty=2)
points(x=log(TabSedim$HDSratio), y=TabSedim$y, pch=16, col=adjustcolor("grey60", alpha.f = 0.6))
points(x=log(sedTable$FishRatio), y=sedTable$y, pch=18, col="cornflowerblue")
box()
axis(1, cex.axis=0.8, tck=-0.01, labels=F)
axis(1, tick=F, line=-1, cex.axis=0.8)
axis(2, at=TabSedim$y, tick=F, labels=TabSedim$ylab, cex.axis=0.45, las=1, line=-0.75)
axis(2, at=TabSedim$y, tck=-0.01, labels=F)
mtext("b", font=2, side=3, adj=0)
mtext(side=2, at=c(5.5, 12.5, 19.5, 26.5, 33.5), cex=0.6, line=3, las=1, adj=0,
      text=c("Mud", "Sand", "Mixed", "Coarse", "Gravel"))

## Include legend
par(fig=c(0,1,0,1), new=T, mar=c(0,0,0,0))
plot(c(0,1), c(0,1), type="n", axes=F, ann=F)
mtext(side=1, adj=0.5, "Relative representation (ln[observed / expected])", line=-2.7, cex=0.8)
legend("bottom", legend= c("Natura 2000 sites", "Habitat-protection sites", "Fishing effort"),
       cex=0.7, col=c("grey20", "grey60", "cornflowerblue"), pch=c(16,16,18), pt.cex=1.5, bty="n", horiz=T)
dev.off()


#-#-# Step 8. Spatial overlap fisheries--MPAs (Figure 3) ----
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
Regiontot <- readRDS(paste0(datPath, "Regiontot3.rds"))
Regtot = Regiontot[,c(1:2, 4, 9, 11)]
VMSdat = merge(Regtot, VMSdatMet[,c(1, 3:102)], on="csquares", all.x=T)
csquares = as.data.table(VMSdat[,c(1:5)])
for(iMet in c("DRB_MOL", "TBB_MOL", "TBB_CRU", "TBB_DMF", "SDN_DMF", "SSC_DMF", "OT_CRU", "OT_DMF", "OT_SPF", "OT_MIX")) {
  period = 2016:2020
  nms1 = paste0(iMet, "_surface_SA_", period)
  nms2 = paste0(iMet, "_total_weight_", period)
  idx1 = which(colnames(VMSdat) %in% nms1)
  csquares[[paste0(iMet, "_SA")]] = rowSums(VMSdat[,idx1], na.rm=T)/5
  csquares[[paste0(iMet, "_SA_MPA")]] = (rowSums(VMSdat[,idx1], na.rm=T)/5) * VMSdat$S1pct
  csquares[[paste0(iMet, "_SA_HDS")]] = (rowSums(VMSdat[,idx1], na.rm=T)/5) * VMSdat$S2pct
  idx2 = which(colnames(VMSdat) %in% nms2)
  csquares[[paste0(iMet, "_TLW")]] = rowSums(VMSdat[,idx2], na.rm=T)/5
  csquares[[paste0(iMet, "_TLW_MPA")]] = (rowSums(VMSdat[,idx2], na.rm=T)/5) * VMSdat$S1pct
  csquares[[paste0(iMet, "_TLW_HDS")]] = (rowSums(VMSdat[,idx2], na.rm=T)/5) * VMSdat$S2pct
} # end iMet-loop

## Aggregate data
SATLWMPA = csquares[, lapply(.SD, sum), .SDcols = colnames(csquares[,6:65])]
test= data.table(t(SATLWMPA[,1:60]))
names(test)= "val"
test$var = colnames(SATLWMPA[,1:60])
test$gear = ifelse(substr(test$var, 1, 2)== "OT", substr(test$var, 1, 6), substr(test$var, 1, 7))
test$var2 = rep(c("SA", "SA", "SA", "TLW", "TLW", "TLW"),10)
test$var3 = rep(c("TOT", "ALL", "HDS"), 20)
df = data.frame(gear = rep(unique(test$gear), 2),
                var = rep(c("SA", "TLW"), each=10),
                TOT = as.numeric(NA),
                ALL = as.numeric(NA),
                HDS = as.numeric(NA))
for(iMet in unique(df$gear)){
  for(iVar in c("SA", "TLW")){
    a = subset(test, gear==iMet & var2 == iVar)
    df$TOT = ifelse(df$gear == iMet & df$var == iVar, a$val[1], df$TOT)
    df$ALL = ifelse(df$gear == iMet & df$var == iVar, a$val[2], df$ALL)
    df$HDS = ifelse(df$gear == iMet & df$var == iVar, a$val[3], df$HDS)
}} # end all loops

## Calculate proportions
dfcomb = aggregate(cbind(TOT, ALL, HDS)~var, data=df, FUN="sum")
dfcomb$gear = "All combined"
df = rbind(df, dfcomb)

df$Hab = df$HDS/df$TOT * 100
df$noHab = (df$ALL - df$HDS)/df$TOT * 100
df$noMPA = (df$TOT - df$ALL)/df$TOT * 100

dftot = aggregate(TOT~var, data=subset(df, !gear == "All combined"), FUN="sum")
names(dftot) = c("var", "overalltot")
df = merge(df, dftot, by=c("var"))    
df$Prop = df$TOT / df$overalltot * 100
df = df[order(df$var, df$gear),]
df$proplab = ifelse(df$Prop < 1, "(< 1%)", paste0("(", round(df$Prop), "%)"))
df$Hablab = ifelse(df$Hab < 6, "", paste0(round(df$Hab), "%"))
saveRDS(df, file=paste0(outPath, "df_gearMPA.rds"))

df$MPAratio = (100-df$noMPA)/df$Prop
df$HDSratio = df$Hab/df$Prop
df$label = ifelse(df$Prop == 0, "",
                  ifelse(round(df$Prop)== 0, "(<1%)", paste0("(", round(df$Prop), "%)")))
df$HabLabel = paste0(round(df$Hab), "%")

##-- Create Fig 3: metier-network importance
tiff(filename=paste0(outPath, "Fig_3.tiff"),  
     width=17, height=8, units="cm", res=500)
layout(mat=matrix(data=c(3,1,2), nrow=1, ncol=3), widths=c(0.1,0.45,0.45))
par(oma=c(1, 0.5, 0.5, 0.5))

# swept area
plotdatSAR = data.frame(t(df[1:11,6:8]))
names(plotdatSAR) = df$gear[1:11]
par(mar=c(2,0.5,0,0))
b = barplot(as.matrix(plotdatSAR), las=1, names.arg = rep("", 11), horiz=T, axes=F, 
            ann=F, xlim=c(0,119), ylim=c(0,15), col=rev(brewer.pal(n=3, "Greys")), space=c(0.2, 1, rep(0.2, 9)))
mtext(side=1, at=c(0,50, 100), line=-1, paste0(c(0,50,100), "%"), cex=0.6)
text(x=100-df$noMPA[1:11], y=b, paste0(round(100-df$noMPA[1:11]), "%"), pos=4, cex=0.6, offset=0.1)
text(x=df$Hab[1:11], y=b, df$HabLabel[1:11], pos=2, col="white", cex=0.6, offset=0.05)
text(x=107, y=b[2:11], df$label[2:11], font=3, cex=0.6)
mtext(side=3, at=107, line=-1.7, "Prop. (%) of\ntotal swept area", cex=0.5, font=3)
mtext(side=3, line=-1, adj=-0, "a", font=2)
mtext(side=1, at=50, "Swept area (%)", line=-0.1, cex=0.6)

# landings
plotdatTLW = data.frame(t(df[12:22,6:8]))
names(plotdatTLW) = df$gear[12:22]
par(mar=c(2,0.5,0,0))
b= barplot(as.matrix(plotdatTLW), las=1, names.arg = rep("", 11), horiz=T, axes=F, 
           ann=F, xlim=c(0,119), ylim=c(0,15), col=rev(brewer.pal(n=3, "Greys")), space=c(0.2, 1, rep(0.2, 9)))
mtext(side=1, at=c(0,50, 100), line=-1, paste0(c(0,50,100), "%"), cex=0.6)
text(x=100-df$noMPA[12:22], y=b, paste0(round(100-df$noMPA[12:22]), "%"), pos=4, cex=0.6, offset=0.1)
text(x=df$Hab[12:22], y=b, df$HabLabel[12:22], pos=2, cex=0.6, col="white", offset=0.05)
text(x=107, y=b[2:11], df$label[13:22], font=3, cex=0.6)
mtext(side=3, at=107, line=-1.7, "Fraction (%) of\ntotal landings", cex=0.5, font=3)
mtext(side=3, adj=-0, "b", font=2, line=-1)
mtext(side=1, at=50, "Landings (%)", line=-0.1, cex=0.6)

# gear names
plot(y=b, x=rep(1, 11), xlim=c(0,1), ylim=c(0,15), type="n", axes=F, ann=F)
axis(2, at=b, tick=F, labels=df$gear[12:22], las=1, line=-5.75)

# legend
par(fig=c(0,1,0,1), new=T, mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(c(0,1), c(0,1), type="n", axes=F, ann=F)
legend("bottom", horiz=T, legend= c("Habitat-protection sites", "Non-habitat protection sites", "Outside Natura 2000 network"), 
       cex=0.9, fill=(brewer.pal(n=3, "Greys")[c(3,2,1)]), bty="n", inset=0)

dev.off()

#-#-# Step 9. Plot LpUE-effort per region (Figure S1) ---- 
##-- Load data
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
Regiontot <- readRDS(paste0(datPath, "Regiontot3.rds"))
Regtot = Regiontot[,c(1:2, 4, 9, 11)]
VMSdat = merge(Regtot, VMSdatMet[,c(1, 3:102)], on="csquares", all.x=T)
csquares = as.data.table(VMSdat[,1:5])

##-- Determine metier-specific landings/activity
for(iMet in c("DRB_MOL", "TBB_MOL", "TBB_CRU", "TBB_DMF", "SDN_DMF", "SSC_DMF", "OT_CRU", "OT_DMF", "OT_SPF", "OT_MIX")) {
  period = 2016:2020
  nms1 = paste0(iMet, "_surface_SA_", period)
  idx1 = which(colnames(VMSdat) %in% nms1)
  csquares[[paste0(iMet, "_SA")]] = rowSums(VMSdat[,idx1], na.rm=T)/5
  nms2 = paste0(iMet, "_total_weight_", period)
  idx2 = which(colnames(VMSdat) %in% nms2)
  csquares[[paste0(iMet, "_TLW")]] = rowSums(VMSdat[,idx2], na.rm=T)/5
} # end iMet-loop

mets = c("DRB_MOL", "TBB_MOL", "TBB_CRU", "TBB_DMF", "SDN_DMF", "SSC_DMF", "OT_CRU", "OT_DMF", "OT_SPF", "OT_MIX")
clnms = c("csquares", "Ecoregion", paste0(mets, "_TLW"), paste0(mets, "_SA"))
subcsquares = csquares[,..clnms]
idx = grep("_SA", colnames(subcsquares))
subcsquares$TOT_SA = rowSums(subcsquares[,..idx], na.rm=T)
idx = grep("_TLW", colnames(subcsquares))
subcsquares$TOT_TLW = rowSums(subcsquares[,..idx], na.rm=T)
subcsquares$TOT_LPUE = subcsquares$TOT_TLW / subcsquares$TOT_SA
subcsquares2 = subset(subcsquares, TOT_SA > 0.5 & TOT_LPUE > 1)

tiff(filename=paste0(outPath, "Figure_S1.tiff"), height = 8, width = 8, units= "cm", res=300)
par(mar=c(2,2,1,1))
plot(TOT_LPUE/1000~TOT_SA, data=subcsquares2, axes=F, ann=F, type="n", ylim=c(0,10), xlim=c(0,100))
axis(1, cex.axis=0.5, tck=-0.01, at=seq(0,100,20), labels=rep("",6))
axis(1, cex.axis=0.5, tick=F, at=seq(0,100,20), labels=seq(0,100,20), line=-1.2)
axis(2, cex.axis=0.5, las=1, at=seq(0,10, 2), labels=rep("",6), tck=-0.01)
axis(2, cex.axis=0.5, las=1, at=seq(0,10, 2), labels=seq(0,10,2), tick=F, line=-0.75)

mtext(side=1, adj=0.5, line=0.75, expression('Fishing effort (Swept area - km'^2*')'), cex=0.6)
mtext(side=2, adj=0.5, line=0.9, expression('Landings per Unit Effort (LpUE - ×1000 kg/km'^2*')'), cex=0.6)
points(TOT_LPUE/1000~TOT_SA, data=subset(subcsquares2, Ecoregion == "Baltic Sea"), cex=0.6, pch=16, col=adjustcolor("plum3", alpha.f = 0.4))
points(TOT_LPUE/1000~TOT_SA, data=subset(subcsquares2, Ecoregion == "Celtic Seas"), cex=0.6, pch=16, col=adjustcolor("coral2", alpha.f = 0.3))
points(TOT_LPUE/1000~TOT_SA, data=subset(subcsquares2, Ecoregion == "Greater North Sea"), cex=0.6, pch=16, col=adjustcolor("steelblue2", alpha.f = 0.2))
points(TOT_LPUE/1000~TOT_SA, data=subset(subcsquares2, Ecoregion == "Bay of Biscay and the Iberian Coast"), cex=0.6, pch=16, col=adjustcolor("palegreen3", alpha.f = 0.1))
box()
legend("topright", bty="n", pch=16, col=c("plum3", "coral2", "steelblue2", "palegreen3"), cex=0.5, pt.cex=0.8,
       legend=c("Baltic Sea", "Celtic Seas", "Greather North Sea", "Bay of Biscay &\nthe Iberian Coast"))
dev.off()

