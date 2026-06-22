#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
##  Rscript for the analyse underlying the article 
#      "Evaluating the efficacy of closing the current European network of marine protected areas 
#       to mobile bottom-contacting fishing gears"
#
##  Step 5: Determine scenario implications
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
library(plotrix)
library(viridis)
sf_use_s2(FALSE)

##-- set data paths
datPath <- "D:/MS_Natura2000closures/Data/"
outPath <- "D:/MS_Natura2000closures/Routput/"
scriptsPath <- "D:/MS_Natura2000closures/Rscripts/"

#-#-# Step 1. Assign annual mean PD - per BHT & ecoregion ----
RegionMSFD  = as.data.frame(readRDS(paste0(datPath, "RegionMSFD3.rds")))
PD_ScAS1_INF   = readRDS(paste0(outPath, "State_reg_ScAS1_inf.rds"))

##-- Assign habitat state
Period = 2016:2020
state_names <- paste0("state_", Period)
idx = which(colnames(PD_ScAS1_INF) %in% state_names)

for(iSC in c("A", "B", "C", "D", "E")){
  for(iMPA in c("S1", "S2")){
    for(iBio in c("epi", "inf")){
      datset = readRDS(paste0(outPath, "State_reg_Sc", iSC, iMPA, "_", iBio, ".rds"))
      datset$AV1620_PD <- rowMeans(datset[,idx], na.rm=T)
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "2016", sep="_")] = datset$state_2016 [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "2017", sep="_")] = datset$state_2017 [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "2018", sep="_")] = datset$state_2018 [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "2019", sep="_")] = datset$state_2019 [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "2020", sep="_")] = datset$state_2020 [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
      RegionMSFD[,paste("PD", iBio, iSC, iMPA, "1620", sep="_")] = datset$AV1620_PD [match(RegionMSFD$csquares, datset$FisheriesMet.csquares)]
    }
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "2016", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "2016", sep="_"),paste("PD_epi", iSC, iMPA, "2016", sep="_")) ], na.rm=T)
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "2017", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "2017", sep="_"),paste("PD_epi", iSC, iMPA, "2017", sep="_")) ], na.rm=T)
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "2018", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "2018", sep="_"),paste("PD_epi", iSC, iMPA, "2018", sep="_")) ], na.rm=T)
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "2019", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "2019", sep="_"),paste("PD_epi", iSC, iMPA, "2019", sep="_")) ], na.rm=T)
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "2020", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "2020", sep="_"),paste("PD_epi", iSC, iMPA, "2020", sep="_")) ], na.rm=T)
    RegionMSFD[,paste("PD_comb", iSC, iMPA, "1620", sep="_")] = rowMeans(RegionMSFD[,c(paste("PD_inf", iSC, iMPA, "1620", sep="_"),paste("PD_epi", iSC, iMPA, "1620", sep="_")) ], na.rm=T)
  }
}  

##-- Remove "epi" and "inf" columns
epi = grep("epi", colnames(RegionMSFD))
RegionMSFD = RegionMSFD[,setdiff(1:ncol(RegionMSFD), epi)]
inf = grep("inf", colnames(RegionMSFD))
RegionMSFD = RegionMSFD[,setdiff(1:ncol(RegionMSFD), inf)]

saveRDS(RegionMSFD, file=paste0(datPath, "RegionMSFD4.rds"))

#-#-# Step 2. Assign annual mean PD - per ecoregion ----
Regiontot   = as.data.frame(readRDS(paste0(datPath, "Regiontot3.rds")))
PD_ScAS1_INF   = readRDS(paste0(outPath, "State_reg_ScAS1_inf.rds"))

## Assign habitat state
Period = 2016:2020
state_names <- paste0("state_", Period)
idx = which(colnames(PD_ScAS1_INF) %in% state_names)

for(iSC in c("A", "B", "C", "D", "E")){
  for(iMPA in c("S1", "S2")){
    for(iBio in c("epi", "inf")){
      datset = readRDS(paste0(outPath, "State_reg_Sc", iSC, iMPA, "_", iBio, ".rds"))
      datset$AV1620_PD <- rowMeans(datset[,idx], na.rm=T)
      Regiontot[,paste("PD", iBio, iSC, iMPA, "2016", sep="_")] = datset$state_2016 [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
      Regiontot[,paste("PD", iBio, iSC, iMPA, "2017", sep="_")] = datset$state_2017 [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
      Regiontot[,paste("PD", iBio, iSC, iMPA, "2018", sep="_")] = datset$state_2018 [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
      Regiontot[,paste("PD", iBio, iSC, iMPA, "2019", sep="_")] = datset$state_2019 [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
      Regiontot[,paste("PD", iBio, iSC, iMPA, "2020", sep="_")] = datset$state_2020 [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
      Regiontot[,paste("PD", iBio, iSC, iMPA, "1620", sep="_")] = datset$AV1620_PD [match(Regiontot$csquares, datset$FisheriesMet.csquares)]
    }
    Regiontot[,paste("PD_comb", iSC, iMPA, "2016", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "2016", sep="_"), paste("PD_inf", iSC, iMPA, "2016", sep="_"))], na.rm=T)
    Regiontot[,paste("PD_comb", iSC, iMPA, "2017", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "2017", sep="_"), paste("PD_inf", iSC, iMPA, "2017", sep="_"))], na.rm=T)
    Regiontot[,paste("PD_comb", iSC, iMPA, "2018", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "2018", sep="_"), paste("PD_inf", iSC, iMPA, "2018", sep="_"))], na.rm=T)
    Regiontot[,paste("PD_comb", iSC, iMPA, "2019", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "2019", sep="_"), paste("PD_inf", iSC, iMPA, "2019", sep="_"))], na.rm=T)
    Regiontot[,paste("PD_comb", iSC, iMPA, "2020", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "2020", sep="_"), paste("PD_inf", iSC, iMPA, "2020", sep="_"))], na.rm=T)
    Regiontot[,paste("PD_comb", iSC, iMPA, "1620", sep="_")] = rowMeans(Regiontot[,c(paste("PD_epi", iSC, iMPA, "1620", sep="_"), paste("PD_inf", iSC, iMPA, "1620", sep="_"))], na.rm=T)
  }
}
   
## remove "epi" and "inf" columns from dataset
epi = grep("epi", colnames(Regiontot))
Regiontot = Regiontot[,setdiff(1:ncol(Regiontot), epi)]
inf = grep("inf", colnames(Regiontot))
Regiontot = Regiontot[,setdiff(1:ncol(Regiontot), inf)]

saveRDS(Regiontot, file=paste0(datPath, "Regiontot4.rds"))

#-#-# Step 3. Calculate quality-extent-relation - per ecoregion (prep fig 4) ----
QET_table <- data.table()
Regiontot <- readRDS(paste0(datPath, "Regiontot4.rds"))

iInd = "PD" ; iBtype = "comb"
for(iReg in c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast")){
  for(iYear in c(2016:2020, 1620)){
    for(iSC in c("A", "B", "C", "D", "E")){
      for(iMPA in c("S1", "S2")){
          subReg <- subset(Regiontot, Ecoregion == iReg)
          subReg2 = subset(subReg, is.na(PD_comb_S0_1620)==F)
          
          table1 <- data.frame(Indicator = iInd,
                               Ecoregion = iReg, 
                               Biotype = iBtype,
                               Scenario = iSC,
                               MPAclose = iMPA, 
                               Period = iYear)
          colID <- which(colnames(subReg2) == paste(iInd, iBtype, iSC, iMPA, iYear, sep="_"))
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
          
          ## Formatting table to nice output  
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
        } # end iMPA
      } # end iSC
    } # end iYear
  } # end iReg
    
##-- add the status quo QET_table
QET_table_S0 = readRDS(paste0(outPath, "QET_table_S0.rds"))
QET_table_S0$Scenario = "SQ"
QET_table_S0$MPAclose = "S0"
QET_table_all = rbind(QET_table_S0, QET_table)

saveRDS(QET_table_all, file=paste0(outPath, "QET_table_all.rds"))


#-#-# Step 4. Calculate quality-extent-relation - BHT-dpt & sed classes overall (prep fig 5) ----
QET_BHT_table <- data.table()
RegionMSFD <- readRDS(paste0(datPath, "RegionMSFD4.rds"))
RegionMSFD <- subset(RegionMSFD, !MSFD_BBHT=="Na")
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



BHTdpts = c("Infralittoral", "Circalittoral", "Offshore circalittoral", "Bathyal")
BHTseds = c("Gravel", "Coarse", "Mixed", "Sand", "Mud")

##-- Create table to save PD values 
subRegMSFD = subset(RegionMSFD, is.na(PD_comb_S0_1620)==F)
for(iYear in c(2016:2020, 1620)){
  for(iVal in c(BHTdpts, BHTseds)){
    for(iSC in c("S0", "A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
      table1 <- data.frame(val = iVal,
                           Scenario = iSC, 
                           Period = iYear)
      
      if(iVal %in% BHTdpts){
        subRegMSFD2 = subset(subRegMSFD, Depthclass == iVal)}
      if(iVal %in% BHTseds){
        subRegMSFD2 = subset(subRegMSFD, Sediment == iVal)}
      
      colID <- which(colnames(subRegMSFD2) == paste("PD_comb", iSC, iYear, sep="_"))
      colnames(subRegMSFD2)[colID] <- "ColofInterest"
      
      Unit.data = subRegMSFD2
      Unit.data = Unit.data[order(Unit.data$ColofInterest, -Unit.data$BHTareakm2),]
      TotBHText = sum(Unit.data$BHTareakm2)
      table1$q.100 = sum(subset(Unit.data, ColofInterest >= 1)$BHTareakm2)/TotBHText
      table1$q.99 = sum(subset(Unit.data, ColofInterest >= .99)$BHTareakm2)/TotBHText
      table1$q.95 = sum(subset(Unit.data, ColofInterest >= .95)$BHTareakm2)/TotBHText
      table1$q.90 = sum(subset(Unit.data, ColofInterest >= .90)$BHTareakm2)/TotBHText
      table1$q.85 = sum(subset(Unit.data, ColofInterest >= .85)$BHTareakm2)/TotBHText
      table1$q.80 = sum(subset(Unit.data, ColofInterest >= .80)$BHTareakm2)/TotBHText
      table1$q.75 = sum(subset(Unit.data, ColofInterest >= .75)$BHTareakm2)/TotBHText
      table1$q.70 = sum(subset(Unit.data, ColofInterest >= .70)$BHTareakm2)/TotBHText
      table1$q.65 = sum(subset(Unit.data, ColofInterest >= .65)$BHTareakm2)/TotBHText
      table1$q.60 = sum(subset(Unit.data, ColofInterest >= .60)$BHTareakm2)/TotBHText
      table1$q.55 = sum(subset(Unit.data, ColofInterest >= .55)$BHTareakm2)/TotBHText
      table1$q.50 = sum(subset(Unit.data, ColofInterest >= .50)$BHTareakm2)/TotBHText
      
      ##-- Formatting table to nice output  
      table2 <- gather(table1, key = "quality", value = "extent",
                       q.100, q.99, q.95, q.90, q.85, q.80, q.75, q.70, q.65, q.60, q.55, q.50)
      table2 <- table2 %>%
        mutate(quality = recode(quality, "q.100" = "1.00",
                                "q.99" = "0.99", "q.95" = "0.95", 
                                "q.90" = "0.90", "q.85" = "0.85",
                                "q.80" = "0.80", "q.75" = "0.75",
                                "q.70" = "0.70", "q.65" = "0.65",
                                "q.60" = "0.60", "q.55" = "0.55",
                                "q.50" = "0.50"))
      
      table2$extent <- round(table2$extent, digits = 3)
      
      ##-- Store results in QET_BHT_table
      QET_BHT_table <- rbind(QET_BHT_table, table2)
      
    } # end iSC
  } # end iVal
} # end iYear
saveRDS(QET_BHT_table, file=paste0(outPath, "QET_dptsed_table.rds"))

#-#-# Step 5. Create figure 4 ----
QET_table =readRDS (file=paste0(outPath, "QET_table_all.rds"))
QET_table$Scen = ifelse(QET_table$Scenario == "SQ", "S0", paste0(QET_table$Scenario, "_", QET_table$MPAclose))
SC = c("A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")

##-- Determine difference with status quo
QET_SQPR = data.table()
for(iYear in c(2016:2020, 1620)){
  for(iReg in unique(QET_table$Ecoregion)){
    QET_SQ = subset(QET_table, Scen =="S0" & Period == iYear & Ecoregion == iReg)
    for(iSC in 1:10){
      subdt = subset(QET_table, Scen == SC[iSC] & Period == iYear & Ecoregion == iReg)
      QET_SQ[[paste0("SC_",iSC)]] = subdt$extent - QET_SQ$extent
    }
   QET_SQPR = rbind(QET_SQPR, QET_SQ)
  }}

##-- plotting settings
yeardt = data.table(year = 2016:2020,
                    shape = 21:25,
                    xmn = c(-0.4, -0.2, 0, 0.2, 0.4))
plotcols = c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
             "skyblue3","skyblue1", "pink3", "pink", 
             "darkolivegreen3","darkolivegreen1")

##-- Start plotting figure 4
tiff(filename=paste0(outPath, "Figure4.tiff"), width=19, height=12, units="cm", res=500)
layout(matrix(data = c(1:4, 5), ncol=1, byrow=T), heights = c(0.21, 0.21, 0.21, 0.21, 0.16))
par(mar=c(0.5,6,2,1))

for(iReg in c("Baltic Sea", "Greater North Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
  dtall = subset(QET_SQPR, Ecoregion == iReg & quality %in% c("0.80", "0.60"))
  ymn = round_any(min(dtall[,10:19]), 0.05, floor)
  ymx = round_any(max(dtall[,10:19]), 0.05, ceiling)
  
  dt1 = subset(QET_SQPR, Period==1620 & Ecoregion == iReg)
  dt_1620 = as.data.frame(t(dt1[c(6,10),c(10:19)]))
  colnames(dt_1620)= c("Q80", "Q60")
  
  b=barplot(as.matrix(dt_1620), beside=T, axes=F, names.arg=rep("",2),
            col=plotcols, ylim=c(ymn, ymx), space = c(0,1.5))
  box()
  abline(v=b[10,1]+1.25, col="grey60", lty=3)
  axis(2, las=1, cex.axis=0.9, at=pretty(x=c(ymn, ymx), n=5), labels=rep("", length(pretty(x=c(ymn, ymx)))), tcl=-0.25)
  axis(2, las=1, cex.axis=0.9, at=pretty(x=c(ymn, ymx), n=5), tick=F, line=-0.3)
  abline(h=0, col="red", lwd=2)
  mtext(side=3, adj=0, font=2, paste0(iReg), cex=0.8)
  
 
  for(iYear in 2016:2020){
    dt2 = subset(QET_SQPR, Period == iYear & Ecoregion == iReg)
    dt = as.data.frame(t(dt2[c(6,10),c(10:19)]))
    points(x=c(b[,1], b[,2])+subset(yeardt, year==iYear)$xmn, y=c(dt[,1], dt[,2]), pch=subset(yeardt, year == iYear)$shape, 
           bg=plotcols, col="black", cex=0.5)
  } # end iYear
  if(iReg == "Bay of Biscay and the Iberian Coast"){
    axis(1, at=b[5,]+0.5, labels=c("QT-80", "QT-60"), las=1)}
}# end iReg

plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
l = legend("bottom", bty="n", ncol=5, legend=rep(c("all", "hps"),5),text.width=0.13,inset=c(0.01,0),
           fill=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                  "skyblue3","skyblue1", "pink3", "pink", 
                  "darkolivegreen3","darkolivegreen1"),
           cex=0.8, y.intersp=0.75)
text(x=l$text$x[c(2,4,6,8,10)], y=0.8, cex=0.8, 
     c("NORED", "SA-EFF", "SA-LPUE","TLW-EFF", "TLW-LPUE"))

par(mar=c(0,5,0,0), mfrow=c(1,1), new=T)
plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
mtext(side=2, adj=0.5, line=3, text="Difference with status quo in\nhabitat extent meeting the QTV", font=3, cex=0.7)
legend("bottomright", legend=2016:2020, pch=21:25, bty="n", title="Year", cex=0.5, inset=c(0.01, 0.02))
dev.off()

##-- Determine difference with status quo for unfished extent
unfished = subset(QET_SQPR, quality %in% c("1.00") & Period == 1620)

#-#-# Step 6. Create figure 5 ----
QET_dptsed = readRDS(paste0(outPath, "QET_dptsed_table.rds"))

##-- Determine difference with status quo
QETdt = data.table()
for(iYear in c(2016:2020, 1620)){
  for(iVal in unique(QET_dptsed$val)){
    QET = subset(QET_dptsed, val == iVal & Period == iYear)
    QET_SQ = subset(QET, Scenario == "S0")
    for(iSC in c("A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
      subQET = subset(QET, Scenario == iSC)
      QET_SQ[[iSC]] = subQET$extent - QET_SQ$extent
    } # end iSC
    QETdt = rbind(QETdt, QET_SQ)
  } # end iVal
} # end iYear
  

yeardt = data.table(year = 2016:2020,
                    shape = 21:25,
                    xmn = c(-0.4, -0.2, 0, 0.2, 0.4))


## Plot figure 5
tiff(filename=paste0(outPath, "Figure5.tiff"), res=300, width=19, height=10, units = "cm")
layout(mat=matrix(data=c(11, rep(12,20), 10,11, rep(1,5), rep(2,5), rep(3,5), rep(4,5), 10, 11, rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 10), nrow=22, ncol=3), 
       heights=c(0.03,rep(0.044, 20), 0.09), widths=c(0.02, 0.48, 0.48))
par(mar=c(1,4,1,0), oma=c(0,0,0,0))
for(iDPT in c("Infralittoral", "Circalittoral", "Offshore circalittoral", "Bathyal")){
  
  subdt = subset(QETdt, val == iDPT & Period == 1620)
  dt_1620 = as.data.frame(t(subdt[c(6,10),c(6:15)]))
  colnames(dt_1620)= c("Q80", "Q60")
  
  b=barplot(as.matrix(dt_1620), beside=T, axes=F, names.arg=rep("",2),
            col=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                  "skyblue3","skyblue1", "pink3", "pink", 
                  "darkolivegreen3","darkolivegreen1"), 
            ylim=c(-0.1, 0.1), space = c(0,1.5))  
  
  for(iYear in 2016:2020){
    dt2 = subset(QETdt, val ==iDPT & Period == iYear)
    dt = as.data.frame(t(dt2[c(6,10),c(6:15)]))
    points(x=c(b[,1], b[,2])+subset(yeardt, year==iYear)$xmn, y=c(dt[,1], dt[,2]), pch=subset(yeardt, year == iYear)$shape, 
           bg=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                "skyblue3","skyblue1", "pink3", "pink", 
                "darkolivegreen3","darkolivegreen1"), col="black", cex=0.5)
  }
  
  abline(v=b[10,1]+1.25, col="grey80", lty=3)
  mtext(side=3, adj=0, text=iDPT, font=2, cex=0.8)
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.05), labels = rep("",5), tcl=-0.25)
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.1), tcl=-0.4, labels=rep("", 3))
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.1), tick=F, line=-0.3)
  abline(h=0, col="red", lwd=2)
  
  box(col="grey")
  axis(2, at=c(-3,3))
  axis(1, at=c(-5,50))
  if(iDPT== "Infralittoral"){
    mtext(side=3, line=0.6, adj=-0.15, text="a", font=2)
  }
  if(iDPT == "Bathyal"){
    axis(side=1, at=b[5,]+0.5, labels=c("QT-80", "QT-60"), tick=F, line=-1)
  }
}# end iDPT-loop

par(mar=c(1,3.5,1,0.5))
for(iSED in c("Gravel", "Coarse", "Mixed", "Sand", "Mud")){
  subdt1620 = subset(QETdt, val == iSED & Period == 1620)
  
  dt_1620 = as.data.frame(t(subdt1620[c(6,10),c(6:15)]))
  colnames(dt_1620)= c("Q80", "Q60")

  b=barplot(as.matrix(dt_1620), beside=T, axes=F, names.arg=rep("",2),
             col=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                   "skyblue3","skyblue1", "pink3", "pink", 
                   "darkolivegreen3","darkolivegreen1"), 
             ylim=c(-0.1, 0.1), space = c(0,1.5))  
  
  for(iYear in 2016:2020){
    dt2 = subset(QETdt, val == iSED & Period == iYear)
    dt = as.data.frame(t(dt2[c(6,10),c(6:15)]))
    points(x=c(b[,1], b[,2])+subset(yeardt, year==iYear)$xmn, y=c(dt[,1], dt[,2]), pch=subset(yeardt, year == iYear)$shape, 
           bg=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                "skyblue3","skyblue1", "pink3", "pink", 
                "darkolivegreen3","darkolivegreen1"), col="black", cex=0.5)
  }
  
  abline(v=b[10,1]+1.25, col="grey80", lty=3)
  mtext(side=3, adj=0, text=iSED, font=2, cex=0.8)
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.05), labels = rep("",5), tcl=-0.25)
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.1), tcl=-0.4, labels=rep("", 3))
  axis(2, las=1, cex.axis=0.9, at=seq(-0.1, 0.1, 0.1), tick=F, line=-0.3)
  abline(h=0, col="red", lwd=2)
  
  box(col="grey")
  axis(2, at=c(-3,3))
  axis(1, at=c(-5,50))
  if(iSED== "Gravel"){
    mtext(side=3, line=0.6, adj=-0.15, text="b", font=2)
  }
  if(iSED == "Mud"){
    axis(side=1, at=b[5,]+0.5, labels=c("QT-80", "QT-60"), tick=F, line=-1)
  }
}# end iSED-loop
par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
l = legend("bottom", bty="n", ncol=5, legend=rep(c("all", "hps"),5),text.width=0.1,
       fill=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
              "skyblue3","skyblue1", "pink3", "pink", 
              "darkolivegreen3","darkolivegreen1"),
       cex=0.8, y.intersp=0.75)
text(x=l$text$x[c(1,3,5,7,9)], y=0.75, cex=0.8, 
     c("NORED", "SA-EFF", "SA-LPUE","TLW-EFF", "TLW-LPUE"))
legend("bottomright", pch=yeardt$shape, legend=yeardt$year, ncol=3, cex=0.7, bty="n", title="Year", inset = c(0,0.02))

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1), axes=F, ann=F, type="n")

par(mar=c(0,0,0,0))
plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
mtext(side=2, adj=0.5, line=-2.5, text="Difference with status quo in\nhabitat extent meeting the QTV", font=3, cex=0.7)
dev.off()

#-#-# Step 7. Create Figs S2-5 (histogram comparison dpt+sed per region) ----
QET_dsr = readRDS(paste0(outPath, "QET_dsr.rds"))

##-- Determine difference with status quo
QETdt = data.table()
for(iYear in c(2016:2020, 1620)){
  for(iReg in unique(QET_dsr$reg)){
    for(iVal in unique(QET_dsr$val)){
      QET = subset(QET_dsr, val == iVal & reg == iReg & Period == iYear)
      QET_SQ = subset(QET, Scenario == "S0")
      for(iSC in c("A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
        subQET = subset(QET, Scenario == iSC)
        QET_SQ[[iSC]] = subQET$extent - QET_SQ$extent
      } # end iSC
      QETdt = rbind(QETdt, QET_SQ)
    } # end iVal
  }# end iReg
}# end iYear

##-- Plotting settings
yeardt = data.table(year = 2016:2020,
                    shape = 21:25,
                    xmn = c(-0.4, -0.2, 0, 0.2, 0.4))
regs = data.table(Ecoregion = c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay & the Iberian Coast"),
                  regcode = c("GNS", "CS", "BS", "BoBIC"))

##-- Start plotting
for(iReg in unique(QET_dsr$reg)){
  tiff(filename=paste0(outPath, "Scen_res_", iReg, ".tiff"), res=300, width=17, height=10, units = "cm")
  layout(mat=matrix(data=c(11, rep(1,5), rep(2,5), rep(3,5), rep(4,5), 10, 11, rep(5,4), rep(6,4), rep(7,4), rep(8,4), rep(9,4), 10), nrow=22, ncol=2), heights=c(0.03,rep(0.044, 20), 0.09))
  par(mar=c(1,6,1,0), oma=c(0,0,0,0))
  for(iDPT in c("Infralittoral", "Circalittoral", "Offshore circalittoral", "Bathyal")){
    
    subdt1620 = subset(QETdt, val == iDPT & reg == iReg & Period == 1620)
    dt_1620 = as.data.frame(t(subdt1620[c(6,10),c(7:16)]))
    colnames(dt_1620)= c("Q80", "Q60")
    ymn = round_any(min(dt_1620), 0.05, floor)
    ymx = round_any(max(dt_1620), 0.05, ceiling)
    ymn = ifelse(ymn == 0, -0.05, ymn) ; ymx = ifelse(ymx == 0, 0.05, ymx)
    
    if(iReg == "BS" & iDPT == "Bathyal"){
      plot(y=c(-1,1), x=c(0,22.5), type="n", axes=F, ann=F)
      box(col="grey")
      axis(1, at=c(-5,50))
      axis(2, at=c(-2,2))
      text("Habitat not present", col="black", font=3, x=12.25, y=0, adj=c(0.5, 0.5))
      mtext(side=3, adj=0, text=iDPT, font=2, cex=0.8)
    } else {
      if(iReg == "GNS" & iDPT == "Bathyal"){
        plot(y=c(-1,1), x=c(0,22.5), type="n", axes=F, ann=F)
        box(col="grey")
        axis(1, at=c(-5,50))
        axis(2, at=c(-2,2))
        text("Total extent < 2%", col="black", font=3, x=12.25, y=0, adj=c(0.5, 0.5))
        mtext(side=3, adj=0, text=iDPT, font=2, cex=0.8)
      } else {
        if(iReg == "CS" & iDPT == "Infralittoral"){
          plot(y=c(-1,1), x=c(0,47), type="n", axes=F, ann=F)
          box(col="grey")
          axis(1, at=c(-5,50))
          axis(2, at=c(-2,2))
          text("Total extent < 2%", col="black", font=3, x=14.5, y=0, adj=c(0.5, 0.5))
          mtext(side=3, adj=0, text=iDPT, font=2, cex=0.8)
        } else {
          b = barplot(as.matrix(dt_1620), beside=T, axes=F, names.arg=rep("",2),
                      col=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                            "skyblue3","skyblue1", "pink3", "pink", 
                            "darkolivegreen3","darkolivegreen1"), 
                      ylim = c(ymn, ymx), space = c(0,1.5))  
          for(iYear in 2016:2020){
            dt2 = subset(QETdt, val == iDPT & reg == iReg & Period == iYear)
            dt = as.data.frame(t(dt2[c(6,10),c(7:16)]))
            points(x=c(b[,1], b[,2])+subset(yeardt, year==iYear)$xmn, y=c(dt[,1], dt[,2]), pch=subset(yeardt, year == iYear)$shape, 
                   bg=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                        "skyblue3","skyblue1", "pink3", "pink", 
                        "darkolivegreen3","darkolivegreen1"), col="black", cex=0.5)
          }# end iYear
          
          abline(v=b[10,1]+1.25, col="grey80", lty=3)
          mtext(side=3, adj=0, text=iDPT, font=2, cex=0.8)
          axis(2, las=1, cex.axis=0.8, at=pretty(c(ymn, ymx), n=5))
          abline(h=0, col="red", lwd=2)
          box(col="grey")
          axis(2, at=c(-3,3))
          axis(1, at=c(-5,50))
        }}}
    if(iDPT== "Infralittoral"){
      mtext(side=3, line=0.6, adj=-0.15, text="a", font=2)
    }
    if(iDPT == "Bathyal"){
      axis(side=1, at=b[5,], labels=c("QT-80", "QT-60"), tick=F, line=-1)
    }
  }# end iDPT-loop
  par(mar=c(1,3.5,1,0.5))
  for(iSED in c("Gravel", "Coarse", "Mixed", "Sand", "Mud")){
    subdt1620 = subset(QETdt, val == iSED & reg == iReg & Period == 1620)
    dt_1620 = as.data.frame(t(subdt1620[c(6,10),c(7:16)]))
    colnames(dt_1620)= c("Q80", "Q60")
    ymn = round_any(min(dt_1620), 0.05, floor)
    ymx = round_any(max(dt_1620), 0.05, ceiling)
    ymn = ifelse(ymn == 0, -0.05, ymn) ; ymx = ifelse(ymx == 0, 0.05, ymx)
    
    if(iReg == "CS" & iSED == "Mixed"){
      plot(y=c(-1,1), x=c(0,47), type="n", axes=F, ann=F)
      box(col="grey")
      axis(1, at=c(-5,50))
      axis(2, at=c(-2,2))
      text("Total extent < 2%", col="black", font=3, x=14.5, y=0, adj=c(0.5, 0.5))
      mtext(side=3, adj=0, text=iSED, font=2, cex=0.8)
    } else {
      if(iReg == "GNS" & iSED == "Gravel") {
        plot(y=c(-1,1), x=c(0,47), type="n", axes=F, ann=F)
        box(col="grey")
        axis(1, at=c(-5,50))
        axis(2, at=c(-2,2))
        text("Total extent < 2%", col="black", font=3, x=14.5, y=0, adj=c(0.5, 0.5))
        mtext(side=3, adj=0, text=iSED, font=2, cex=0.8)
      } else {
        b = barplot(as.matrix(dt_1620), beside=T, axes=F, names.arg=rep("",2),
                    col=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                          "skyblue3","skyblue1", "pink3", "pink", 
                          "darkolivegreen3","darkolivegreen1"), 
                    ylim=c(ymn, ymx), space = c(0,1.5))  
        for(iYear in 2016:2020){
          dt2 = subset(QETdt, val == iSED & reg == iReg & Period == iYear)
          dt = as.data.frame(t(dt2[c(6,10),c(7:16)]))
          points(x=c(b[,1], b[,2])+subset(yeardt, year==iYear)$xmn, y=c(dt[,1], dt[,2]), pch=subset(yeardt, year == iYear)$shape, 
                 bg=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                      "skyblue3","skyblue1", "pink3", "pink", 
                      "darkolivegreen3","darkolivegreen1"), col="black", cex=0.5)
        }# end iYear
        
        abline(v=b[10,1]+1.25, col="grey80", lty=3)
        mtext(side=3, adj=0, text=iSED, font=2, cex=0.8)
        axis(2, las=1, cex.axis=0.8, at=pretty(c(ymn, ymx), n=5))
        abline(h=0,col="red", lwd=2)
        box(col="grey")
        axis(2, at=c(-3,3))
        axis(1, at=c(-5,50))
      }}
    if(iSED== "Gravel"){
      mtext(side=3, line=0.6, adj=-0.15, text="b", font=2)
      mtext(side=3, line=0.6, adj=1, text=subset(regs, regcode == iReg)$Ecoregion, cex=0.8, font=3)
    }
    if(iSED == "Mud"){
      axis(side=1, at=b[5,]+0.5, labels=c("QT-80", "QT-60"), tick=F, line=-1)
    }
  }# end iSED-loop
  par(mar=c(0,0,0,0))
  plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
  l = legend("bottom", bty="n", ncol=5, legend=rep(c("all", "hps"),5),text.width=0.1,
             fill=c("mediumseagreen", "mediumspringgreen", "tomato3", "tomato1", 
                    "skyblue3","skyblue1", "pink3", "pink", 
                    "darkolivegreen3","darkolivegreen1"),
             cex=0.8, y.intersp=0.75)
  text(x=l$text$x[c(1,3,5,7,9)], y=0.75, cex=0.8, 
       c("NORED", "SA-EFF", "SA-LPUE","TLW-EFF", "TLW-LPUE"))
  
  par(mar=c(0,5,0,0), new=T, mfrow=c(1,1), oma=c(0,0,0,0))
  plot(c(0,1), c(0,1), axes=F, ann=F, type="n")
  mtext(side=2, adj=0.5, line=3, text="Difference with status quo in\nhabitat extent meeting the QTV", font=3, cex=0.7)
  legend("bottomright", legend=2016:2020, pch=21:25, bty="n", title="Year", cex=0.5, inset=c(0, 0), ncol=3)
  dev.off()
} # end iReg

#-#-# Step 8. Create Figures S6-13 (continuous comparison dpt+sed per region) ----
QET_dsr =readRDS(paste0(outPath, "QET_dsr.rds"))
QET_dsr = subset(QET_dsr, Scenario %in% c("S0", "A_S1", "B_S1", "C_S1", "D_S1", "E_S1"))
QET_dsr$qual2 = as.numeric(QET_dsr$quality)
plotdt = data.table(scenario = unique(QET_dsr$Scenario),
                    LWD = c(1.5, 1,1,1,1,1),
                    label = c("Status quo",  "NORED", "SA-EFF", "SA-LPUE", "TLW-EFF", "TLW-LPUE"),
                    cols = c("black", "#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"))

regdt = data.table(reg = c("GNS", "BS", "CS", "BoBIC"),
                   name = c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay &\nthe Iberian Coast"))

for(iReg in c("GNS", "BS", "CS", "BoBIC")){
  dat = subset(QET_dsr, reg == iReg & Period == 1620)
  
  ##-- make dpt-class plot
  tiff(filename=paste0(outPath, iReg, "_dpt_QETplot.tiff"), width=17, height = 10, res=300, units="cm")
  layout(matrix(data=c(7,1,2,5,7,3,4,5,7,6,6,5), nrow=3, byrow=T), widths = c(0.03, 0.4, 0.4, 0.17), heights = c(0.48, 0.48, 0.04))
  par(mar=c(2,2,2,1))
  
  for(iDpt in c("Infralittoral","Circalittoral","Offshore circalittoral","Bathyal")){
    dat1 = subset(dat, val == iDpt)
    if(iReg == "BS" & iDpt == "Bathyal"){
      plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
      box()
      axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
      axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
      axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
      axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
      mtext(iDpt, side=3, adj=0)
      text(x=0.7, y=0.7, "Habitat not present", font=3)
    } else {
      if(iReg == "CS" & iDpt == "Infralittoral"){
        plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
        box()
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
        axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
        axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
        mtext(iDpt, side=3, adj=0)
        text(x=0.7, y=0.7, "Total extent <2%", font=3)
      } else {
        if(iReg == "GNS" & iDpt == "Bathyal"){
          plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
          box()
          axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
          axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
          axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
          axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
          mtext(iDpt, side=3, adj=0)
          text(x=0.7, y=0.7, "Total extent <2%", font=3)
        } else {
          plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
          box()
          axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
          axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
          axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
          axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
          abline(v=c(0.6, 0.8), lty=3, lwd=0.7, col="dimgrey")
          mtext(iDpt, side=3, adj=0)
          lines(x=subset(dat1, Scenario == "A_S1")$qual2, y=subset(dat1, Scenario == "A_S1")$extent, 
                col=subset(plotdt, scenario == "A_S1")$cols, lwd=1.5)
          lines(x=subset(dat1, Scenario == "B_S1")$qual2, y=subset(dat1, Scenario == "B_S1")$extent, 
                col=subset(plotdt, scenario == "B_S1")$cols, lwd=1.5)
          lines(x=subset(dat1, Scenario == "C_S1")$qual2, y=subset(dat1, Scenario == "C_S1")$extent, 
                col=subset(plotdt, scenario == "C_S1")$cols, lwd=1.5)
          lines(x=subset(dat1, Scenario == "D_S1")$qual2, y=subset(dat1, Scenario == "D_S1")$extent, 
                col=subset(plotdt, scenario == "D_S1")$cols, lwd=1.5)
          lines(x=subset(dat1, Scenario == "E_S1")$qual2, y=subset(dat1, Scenario == "E_S1")$extent, 
                col=subset(plotdt, scenario == "E_S1")$cols, lwd=1.5)
          lines(x=subset(dat1, Scenario == "S0")$qual2, y=subset(dat1, Scenario == "S0")$extent, 
                col=subset(plotdt, scenario == "S0")$cols, lwd=2)
        }}}
  }# end iDpt
  
  par(new=T, fig=c(0,1,0,1), mar=c(2,2,1,1))
  plot(x=c(0,1), y=c(0,1), type="n", axes=F, ann=F)
  legend("right", legend=plotdt$label, fill=plotdt$cols, bty="n",
         title="Scenario")
  mtext(side=1, adj=0.415, "Habitat quality (RBS)", line=0.7, cex=0.8)
  mtext(side = 2, adj=0.5, "Habitat extent (%)", line=0.6, cex=0.8)
  mtext(side = 1, adj=1, subset(regdt, reg== iReg)$name, cex=0.8, font=3)
  dev.off()
  
  ##-- make sed-class plot
  tiff(filename=paste0(outPath, iReg, "_sed_QETplot.tiff"), width=17, height = 10, res=300, units="cm")
  layout(matrix(data=c(7,1,2,3,7,4,5,6,7,8,8,8), nrow=3, byrow=T), widths = c(0.03, 0.31, 0.31, 0.31), heights = c(0.48, 0.48, 0.04))
  par(mar=c(2,2,2,1))
  
  for(iSed in c("Gravel","Coarse","Mixed","Sand","Mud")){
    dat2 = subset(dat, val == iSed)
    if(iReg == "CS" & iSed == "Mixed"){
      plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
      box()
      axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
      axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
      axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
      axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
      mtext(iDpt, side=3, adj=0)
      text(x=0.7, y=0.7, "Total extent <2%", font=3)
    } else {
      if(iReg == "GNS" & iSed == "Gravel"){
        plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
        box()
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
        axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
        axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
        mtext(iDpt, side=3, adj=0)
        text(x=0.7, y=0.7, "Total extent <2%", font=3)
      } else {
        plot(x=c(0.4,1), y=c(0.4,1), type="n", axes=F, ann=F)
        box()
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tck=-0.02, labels=rep("",7))
        axis(1, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.75, cex.axis=0.8)
        axis(2, at=seq(0.4, 1, 0.1), tck=-0.02, labels=rep("",7))
        axis(2, at=c(0.4,0.5,0.6,0.7,0.8,0.9,1), tick=F,line=-0.5, las=1, cex.axis=0.8)
        abline(v=c(0.6, 0.8), lty=3, lwd=0.7, col="dimgrey")
        mtext(iSed, side=3, adj=0)
        lines(x=subset(dat2, Scenario == "A_S1")$qual2, y=subset(dat2, Scenario == "A_S1")$extent, 
              col=subset(plotdt, scenario == "A_S1")$cols, lwd=1.5)
        lines(x=subset(dat2, Scenario == "B_S1")$qual2, y=subset(dat2, Scenario == "B_S1")$extent, 
              col=subset(plotdt, scenario == "B_S1")$cols, lwd=1.5)
        lines(x=subset(dat2, Scenario == "C_S1")$qual2, y=subset(dat2, Scenario == "C_S1")$extent, 
              col=subset(plotdt, scenario == "C_S1")$cols, lwd=1.5)
        lines(x=subset(dat2, Scenario == "D_S1")$qual2, y=subset(dat2, Scenario == "D_S1")$extent, 
              col=subset(plotdt, scenario == "D_S1")$cols, lwd=1.5)
        lines(x=subset(dat2, Scenario == "E_S1")$qual2, y=subset(dat2, Scenario == "E_S1")$extent, 
              col=subset(plotdt, scenario == "E_S1")$cols, lwd=1.5)
        lines(x=subset(dat2, Scenario == "S0")$qual2, y=subset(dat2, Scenario == "S0")$extent, 
              col=subset(plotdt, scenario == "S0")$cols, lwd=2)
      }}
  }# end iSed
  
  plot(x=c(0,1), y=c(0,1), type="n", axes=F, ann=F)
  legend("center", legend=plotdt$label, fill=plotdt$cols, bty="n",
         title="Scenario")
  
  par(new=T, fig=c(0,1,0,1), mar=c(2,2,1,1))
  plot(x=c(0,1), y=c(0,1), type="n", axes=F, ann=F)
  mtext(side=1, adj=0.5, "Habitat quality (RBS)", line=0.7, cex=0.8)
  mtext(side = 2, adj=0.5, "Habitat extent (%)", line=0.6, cex=0.8)
  mtext(side = 1, adj=1, subset(regdt, reg== iReg)$name, cex=0.8, font=3)
  dev.off()
}# end iReg


#-#-# Step 9. Create Tables S2-5 ----
RegionMSFD <- readRDS(paste0(datPath, "RegionMSFD4.rds"))
RegionMSFD = subset(RegionMSFD, is.na(PD_comb_S0_1620)==F)
QET_table_all = readRDS(paste0(outPath, "QET_table_all.rds"))

TabH <- data.table(MSFD = c("Total", "Infralittoral rock and biogenic reef", "Infralittoral coarse sediment",
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
                            "Abyssal", "Na"), 
                   habOrd = 1:27, 
                   MSFDName = c("Total", "Infr. rock and biog. reef", "Infr. coarse sed.",
                                "Infr. mixed sed.", "Infr. sand", 
                                "Infr. mud", "Infr. mud or sand", 
                                "Circ. rock and biog. reef", "Circ. coarse sed.",
                                "Circ. mixed sed.", "Circ. sand", 
                                "Circ. mud", "Circ. mud or sand",
                                "Off. circ. rock and biog. reef", "Off. circ. coarse sed.",
                                "Off. circ. mixed sed.", "Off. circ. sand", 
                                "Off. circ. mud", "Off circ. mud or sand",
                                "Up. bath. rock and biog. reef", "Up. bath. sed.", 
                                "Up. bath. sed. or rock and biog. reef",
                                "Low. bath. rock and biog. reef", "Low. bath. sed.", 
                                "Low. bath. sed. or rock and biog. reef",
                                "Abyssal", "Unknown"))

tabInd = data.table()
for(iReg in c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast")){
  for(iYear in c(2016:2020, 1620)){
    MA_BHT <- unique(subset(RegionMSFD, Ecoregion == iReg)$MSFD_BBHT)
    subRegMSFD <- subset(RegionMSFD, Ecoregion == iReg)
    table1 <- data.frame(MSFD = MA_BHT,
                         Ecoregion = iReg,
                         Period = iYear)
    tabEXT = aggregate(BHTareakm2~ MSFD_BBHT, data=subRegMSFD, FUN="sum")
    table1$BHT_EXT = round(tabEXT$BHTareakm2 / sum(tabEXT$BHTareakm2) *100, digits=0) [match(table1$MSFD, tabEXT$MSFD_BBHT)]
    
    for(iSC in c("S0", "A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
      subRegMSFD2 <- subRegMSFD
      colID <- which(colnames(subRegMSFD2) == paste("PD_comb", iSC, iYear, sep="_"))
      colnames(subRegMSFD2)[colID] <- "ColofInterest"
      
      for(iBHT in 1:length(MA_BHT)){
        Unit.data = subset(subRegMSFD2, MSFD_BBHT == MA_BHT[iBHT])
        if(nrow(Unit.data)>0){
          Unit.data = Unit.data[order(Unit.data$ColofInterest, -Unit.data$BHTareakm2),]
          TotBHText = sum(Unit.data$BHTareakm2)
          table1[[paste0(iSC, "_100")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= 1)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_99")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .99)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_95")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .95)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_90")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .90)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_85")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .85)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_80")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .80)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_75")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .75)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_70")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .70)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_65")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .65)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_60")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .60)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_55")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .55)$BHTareakm2)/TotBHText *100, digits=0)
          table1[[paste0(iSC, "_50")]][iBHT] = round(sum(subset(Unit.data, ColofInterest >= .50)$BHTareakm2)/TotBHText *100, digits=0)
        } else {
          table1[[paste0(iSC, "_100")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_99")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_95")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_90")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_85")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_80")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_75")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_70")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_65")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_60")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_55")]][iBHT] = as.numeric(NA)
          table1[[paste0(iSC, "_50")]][iBHT] = as.numeric(NA)
        }
      } # end iBHT
    } # end iSC
    
  ### merge in the total estimates per ecoregion (from QET_table_all)
  a = subset(QET_table_all, Biotype=="comb" & Ecoregion == iReg & Indicator == "PD" & Period == iYear)
  a$Scenario2 = ifelse(a$Scenario == "SQ", "S0", paste0(a$Scenario, "_", a$MPAclose))
  tabtot = data.table(MSFD = "Total",
                      Ecoregion = iReg,
                      BHT_EXT = as.numeric(100), 
                      Period = iYear)
  for(iSc in c("S0", "A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
    tabtot[[paste0(iSc, "_100")]] = round(subset(a, Scenario2 == iSc & quality == "1.00")$extent *100)
    tabtot[[paste0(iSc, "_99")]] = round(subset(a, Scenario2 == iSc & quality == "0.99")$extent *100)
    tabtot[[paste0(iSc, "_95")]] = round(subset(a, Scenario2 == iSc & quality == "0.95")$extent *100)
    tabtot[[paste0(iSc, "_90")]] = round(subset(a, Scenario2 == iSc & quality == "0.90")$extent *100)
    tabtot[[paste0(iSc, "_85")]] = round(subset(a, Scenario2 == iSc & quality == "0.85")$extent *100)
    tabtot[[paste0(iSc, "_80")]] = round(subset(a, Scenario2 == iSc & quality == "0.80")$extent *100)
    tabtot[[paste0(iSc, "_75")]] = round(subset(a, Scenario2 == iSc & quality == "0.75")$extent *100)
    tabtot[[paste0(iSc, "_70")]] = round(subset(a, Scenario2 == iSc & quality == "0.70")$extent *100)
    tabtot[[paste0(iSc, "_65")]] = round(subset(a, Scenario2 == iSc & quality == "0.65")$extent *100)
    tabtot[[paste0(iSc, "_60")]] = round(subset(a, Scenario2 == iSc & quality == "0.60")$extent *100)
    tabtot[[paste0(iSc, "_55")]] = round(subset(a, Scenario2 == iSc & quality == "0.55")$extent *100)
    tabtot[[paste0(iSc, "_50")]] = round(subset(a, Scenario2 == iSc & quality == "0.50")$extent *100)
  }
  table1 = rbind(table1, tabtot)
  table1 = merge(TabH, table1, by="MSFD", all=T)
  tabInd = rbind(tabInd, table1)
} # end iReg-loop 
}# end iYear-loop

## some formatting of the table
tabInd$Ecoreg2 =ifelse(tabInd$Ecoregion == "Bay of Biscay and the Iberian Coast", 4,
                       ifelse(tabInd$Ecoregion == "Baltic Sea", 1,
                              ifelse(tabInd$Ecoregion == "Greater North Sea", 2,
                                     ifelse(tabInd$Ecoregion == "Celtic Seas", 3, 5))))
tabInd = tabInd[order(tabInd$Ecoreg2, tabInd$habOrd)]
idx = c(2,6,139,5, grep("_100", colnames(tabInd)), grep("_80", colnames(tabInd)), grep("_60", colnames(tabInd)))
tabInd2 = tabInd[,..idx]

TI = data.table()
for(iYear in c(2016:2020, 1620)){
  TI1= subset(tabInd2, Period == iYear)
  TI2 = TI1[,1:4]
  
  for(iQT in c("_100", "_80", "_60")){
    for(iSC in c("A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
      TI2 [[paste0(iSC, iQT)]] = TI1[[paste0(iSC, iQT)]] - TI1[[paste0("S0", iQT)]]
      
    }}
  TI = rbind(TI, TI2)
}

TI2 = subset(TI,Period == 1620)
TI3 = subset(TI, Period %in% 2016:2020)
TImin = TI3[, lapply(.SD, min), by=.(habOrd, BHT_EXT, Ecoreg2), 
          .SDcols = colnames(TI3)[5:34]]
TImax = TI3[, lapply(.SD, max), by=.(habOrd, BHT_EXT, Ecoreg2), 
            .SDcols = colnames(TI3)[5:34]]


TabInd3 = TI2[,1:3]
for(iQT in c("_100", "_80", "_60")){
  for(iSC in c("A_S1", "A_S2", "B_S1", "B_S2", "C_S1", "C_S2", "D_S1", "D_S2", "E_S1", "E_S2")){
    TabInd3[[paste0(iSC,iQT)]] = paste0(TI2[[paste0(iSC,iQT)]], " [", TImin[[paste0(iSC, iQT)]], "-", TImax[[paste0(iSC, iQT)]], "]")
  }
}


write.csv(TI2, file=paste0(outPath, "Table_SupMat_S2S5.csv"), row.names = F)
 


#-#-# Step 10. Create Table 4 (scenario consequences SA & TLW) ----
##-- Load scenario-outputs with updated fishing activity
FishMet = readRDS(file=paste0(datPath, "FD_capped.rds"))
FishMetA1 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScA_S1.rds")))
FishMetA2 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScA_S2.rds")))
FishMetB1 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScB_S1.rds")))
FishMetB2 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScB_S2.rds")))
FishMetC1 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScC_S1.rds")))
FishMetC1$Ecoregion = FishMetA1$Ecoregion [match(FishMetC1$csquares, FishMetA1$csquares)]
FishMetC2 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScC_S2.rds")))
FishMetC2$Ecoregion = FishMetA1$Ecoregion [match(FishMetC2$csquares, FishMetA1$csquares)]
FishMetD1 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScD_S1.rds")))
FishMetD1$Ecoregion = FishMetA1$Ecoregion [match(FishMetD1$csquares, FishMetA1$csquares)]
FishMetD2 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScD_S2.rds")))
FishMetD2$Ecoregion = FishMetA1$Ecoregion [match(FishMetD2$csquares, FishMetA1$csquares)]
FishMetE1 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScE_S1.rds")))
FishMetE1$Ecoregion = FishMetA1$Ecoregion [match(FishMetE1$csquares, FishMetA1$csquares)]
FishMetE2 = as.data.table(readRDS(file=paste0(outPath, "FishMet4_ScE_S2.rds")))
FishMetE2$Ecoregion = FishMetA1$Ecoregion [match(FishMetE2$csquares, FishMetA1$csquares)]

dat = c("FishMet", "FishMetA1", "FishMetA2", "FishMetB1", "FishMetB2", "FishMetC1", "FishMetC2", "FishMetD1","FishMetD2", "FishMetE1", "FishMetE2")
dtall = data.table()
regs = data.table(Ecoregion = c("Greater North Sea", "Celtic Seas", "Baltic Sea", "Bay of Biscay and the Iberian Coast"),
                  regcode = c("GNS", "CS", "BS", "BoBIC"))

##-- Determine scenario effects
FishMet[is.na(FishMet),] = 0
for(iDat in 1:11){
  subdat = get(dat[iDat])
  dtregs = data.table(matrix(nrow=2, ncol=7, data=as.numeric(NA)))
  names(dtregs) = c("Val", "Scenario", regs$regcode, "Total")
  dtregs$Val = c("SA", "TLW")
  dtregs$Scenario = iDat
  
  ## where needed, turn SAR into Swept Area estimates
  if(iDat %in% 2:11){
    idx = grep("surface_sar", colnames(subdat))
    SAdat = as.data.table(subdat)[,..idx]
    SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= "*", STAT=subdat$area_sqkm))
    names(SAdat2) = gsub("surface_sar", replacement= "surface_SA", x=names(SAdat2))
    idy = setdiff(1:ncol(subdat), idx)
    subdat = cbind(as.data.table(subdat)[,..idy], SAdat2)} # end iDat-if
  
  for(iReg in 1:4){
    subdat2 = subset(subdat, Ecoregion  == regs$Ecoregion[iReg])
    idx = grep("surface_SA", colnames(subdat2))
    idy = grep("total_weight", colnames(subdat2))
    SA = colSums(subdat2[, ..idx])
    TLW = colSums(subdat2[,..idy])
    y1620_SA = c(grep("2016", names(SA)), grep("2017", names(SA)), grep("2018", names(SA)), grep("2019", names(SA)), grep("2020", names(SA)))
    y1620_TLW = c(grep("2016", names(TLW)), grep("2017", names(TLW)), grep("2018", names(TLW)), grep("2019", names(TLW)), grep("2020", names(TLW)))
    avanSA1620 = round((sum(SA[y1620_SA])/5)/1000, digits=2)
    avanTLW1620 = round((sum(TLW[y1620_TLW])/5)/1000000, digits=2)
    
    ##- store in table dtregs
    dtregs[1,iReg+2] = avanSA1620
    dtregs[2,iReg+2] = avanTLW1620
    } # end iReg-loop
  
  ##- Determine total
  dtregs[1,7] = sum(dtregs[1, 3:6])
  dtregs[2,7] = sum(dtregs[2, 3:6])
  
  ##-= store in dtall
  dtall = rbind(dtall, dtregs)
} # end iDat-loop
write.csv2(dtall, file=paste0(outPath, "SA&TLW_scenarios.csv"), row.names = FALSE)
