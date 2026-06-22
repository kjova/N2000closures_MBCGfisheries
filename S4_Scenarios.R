#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
##  Rscript for the analyse underlying the article 
#      "Evaluating the efficacy of closing the current European network of marine protected areas 
#       to mobile bottom-contacting fishing gears"
#
##  Step 4: Run the displacement scenarios
#
##  Code by Karin van der Reijden. (kjvdreijden@gmail.com)
#        published under the MIT-license
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-# Step 0. Load libraries, set data paths ----
library(sf)
library(data.table)
library(tidyverse)
sf_use_s2(FALSE)

##-- set data paths
datPath <- "D:/MS_Natura2000closures/Data/"
outPath <- "D:/MS_Natura2000closures/Routput/"
scriptsPath <- "D:/MS_Natura2000closures/Rscripts/"

#-#-# Scenario A: NORED ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
CsquareInfo$S1pct <- ifelse(CsquareInfo$S1pct > 0.95, 1, CsquareInfo$S1pct)
CsquareInfo$S2pct <- ifelse(CsquareInfo$S2pct > 0.95, 1, CsquareInfo$S2pct)

source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo (https://github.com/ices-eg/FBIT)

##-- remove ecoregion and annual SA estimates from VMSdatMet. change NA to 0
VMSdatMet = VMSdatMet[,c(1,3:102)]
VMSdatMet[is.na(VMSdatMet)] = 0

for(iSC in c("S1", "S2")){
  FishMet4 = data.frame()
  
  ##-- Determine remaining fishing effort per ecoregion
  for(iReg in c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
    print(paste0(iReg, ". Scenario ", iSC))
    CsquareInfoReg = subset(CsquareInfo, Ecoregion == iReg)
    
    FishMet = as.data.frame(merge(CsquareInfoReg, VMSdatMet, by="csquares"))
    SA = grep("surface_SA", colnames(FishMet))
    TLW = grep("total_weight", colnames(FishMet))
    FishMet = FishMet[,c(1:17, SA, TLW)]
    print(paste0("Total SA: ", round(sum(colSums(FishMet[,18:67])))))
    print(paste0("Total TLW: ", round(sum(colSums(FishMet[,68:117])))))
    FishMet[,18:67]  <- FishMet[,18:67]*(1-FishMet[,paste0(iSC, "pct")]) # remaining fishing
    FishMet[,18:67][is.na(FishMet[,18:67])] <- 0
    print(paste0("Total remaining SA: ", round(sum(colSums(FishMet[,18:67])))))
    FishMet[,68:117] <- FishMet[,68:117]*(1-FishMet[,paste0(iSC, "pct")]) # remaining landings
    FishMet[,68:117][is.na(FishMet[,68:117])] <- 0
    print(paste0("Total remaining TLW:", round(sum(colSums(FishMet[,68:117])))))

    FishMet4 = rbind(FishMet4, FishMet)
  } # end iReg-loop
  
  ##-- Turn remaining SA into SAR
  idx = grep("surface_SA", colnames(FishMet4))
  SAdat = as.data.table(FishMet4)[,..idx]
  SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= function(x){ x/FishMet4$area_sqkm}))
  names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
  idy = setdiff(1:ncol(FishMet4), idx)
  FishMet4 = cbind(as.data.table(FishMet4)[,..idy], SAdat2)
  
  ##-- Save final vmsdatfile used
  saveRDS(FishMet4, file=paste0(outPath, "FishMet4_ScA_", iSC, ".rds"))
  
  ## Determine state
  period = 2016:2020
  source(paste0(scriptsPath, "HSF_Scenarios_inf.R")) # adjusted from FBIT repo (https://github.com/ices-eg/FBIT)
  saveRDS(State_reg, paste0(outPath, "State_reg_ScA",iSC, "_inf.rds"))

  source(paste0(scriptsPath, "HSF_Scenarios_epi.R")) # adjusted from FBIT repo (https://github.com/ices-eg/FBIT)
  saveRDS(State_reg, paste0(outPath, "State_reg_ScA", iSC, "_epi.rds"))
}# end iSC-loop


#-#-# Scenario B: SA-EFF ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
CsquareInfo$S1pct <- ifelse(CsquareInfo$S1pct > 0.95, 1, CsquareInfo$S1pct)
CsquareInfo$S2pct <- ifelse(CsquareInfo$S2pct > 0.95, 1, CsquareInfo$S2pct)

##-- remove ecoregion and annual SA estimates from VMSdatMet
VMSdatMet = VMSdatMet[,c(1,3:102)]
VMSdatMet[is.na(VMSdatMet)] = 0

source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo 

for(iSC in c("S1", "S2")){
  FishMet4 = data.frame()
  
  ##-- Determine remaining fishing effort per ecoregion
  for(iReg in c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
    print(paste0(iReg, ". Scenario ", iSC))
    CsquareInfoReg = subset(CsquareInfo, Ecoregion == iReg)
    
    FishMet <- as.data.frame(merge(CsquareInfoReg, VMSdatMet,by="csquares")) # status quo
    SA = grep("surface_SA", colnames(FishMet))
    TLW = grep("total_weight", colnames(FishMet))
    FishMet = FishMet[,c(1:17, SA, TLW)]
    print(paste0("Total SA: ", round(sum(colSums(FishMet[,18:67])))))
    print(paste0("Total TLW: ", round(sum(colSums(FishMet[,68:117])))))
    
    FishMet2 <- FishMet
    FishMet2[,18:67] <- FishMet2[,18:67]*(FishMet2[,paste0(iSC, "pct")]) # removed effort
    FishMet2[,18:67][is.na(FishMet2[,18:67])] <- 0
    print(paste0("Total removed SA: ", round(sum(colSums(FishMet2[,18:67])))))
    
    FishMet3 <- FishMet
    FishMet3[,18:67] <- FishMet3[,18:67]*(1-FishMet3[,paste0(iSC, "pct")]) # remaining fishing
    FishMet3[,18:67][is.na(FishMet3[,18:67])] <- 0
    print(paste0("Total remaining SA: ", round(sum(colSums(FishMet3[,18:67])))))
    
    Sc0_FI <- colSums(FishMet[,18:67]) # This is what is done originally
    ScX_FIremoved <- colSums(FishMet2[,18:67]) # This is what is removed due to closures
    ScX_FIremain <- colSums(FishMet3[,18:67]) # This is what remains after closures (no displacement)
    ScX_Rmvd_pct <- 1 + ScX_FIremoved/ScX_FIremain # Determine multiplication factor per gear/year-combi
    
    FishMetD <- FishMet3
    FishMetD[,18:67] <- sweep(FishMetD[,18:67], MARGIN=2, ScX_Rmvd_pct, "*") # new effort after displacement
    FishMetD[,18:67][is.na(FishMetD[,18:67])] <- 0
    print(paste0("Total new SA: ", round(sum(colSums(FishMetD[,18:67])))))
    ScX_newFI <- colSums(FishMetD[,18:67])
    
    ##-- Determine (old) LpUE per c-square per gear/year
    LpUE <- cbind(FishMet$csquares, FishMet[,68:117] / FishMet[,18:67], FishMet[paste0(iSC, "pct")])
    LpUE[is.na(LpUE)] <- 0 # if no kg and no effort is registered, this yields a NaN, which is set to 0
    a <- c("csquares", gsub(x=colnames(LpUE)[2:51], pattern="total_weight", replacement="LpUE"), paste0(iSC, "pct"))
    colnames(LpUE) <- a
    
    ##-- Determine (new) landings based on new effort and old LpUE
    FishMetD[,68:117] = FishMetD[,18:67] * LpUE[,2:51]
    print(paste0("Total new TLW: ", round(sum(colSums(FishMetD[,68:117])))))
    
    FishMet4 = rbind(FishMet4, FishMetD)
  } # end iReg-loop
  
  ##-- Turn updated SA into SAR
  idx = grep("surface_SA", colnames(FishMet4))
  SAdat = as.data.table(FishMet4)[,..idx]
  SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= function(x){ x/FishMet4$area_sqkm}))
  names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
  idy = setdiff(1:ncol(FishMet4), idx)
  FishMet4 = cbind(as.data.table(FishMet4)[,..idy], SAdat2)
  
  ##-- Save final vmsdatfile used
  saveRDS(FishMet4, file=paste0(outPath, "FishMet4_ScB_", iSC, ".rds"))
  
  ## Determine state
  Period=2016:2020
  source(paste0(scriptsPath, "HSF_Scenarios_inf.R")) # adjusted from FBIT repo 
  saveRDS(State_reg, paste0(outPath, "State_reg_ScB",iSC, "_inf.rds"))
 
  source(paste0(scriptsPath, "HSF_Scenarios_epi.R")) # adjusted from FBIT repo 
  saveRDS(State_reg, paste0(outPath, "State_reg_ScB", iSC, "_epi.rds"))
}# end iSC-loop


#-#-# Scenario C: SA-LPUE ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
CsquareInfo$S1pct <- ifelse(CsquareInfo$S1pct > 0.95, 1, CsquareInfo$S1pct)
CsquareInfo$S2pct <- ifelse(CsquareInfo$S2pct > 0.95, 1, CsquareInfo$S2pct)

VMSdatMet = VMSdatMet[,c(1,3:102)]
VMSdatMet[is.na(VMSdatMet)] = 0

source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo 

for(iSC in c("S1", "S2")){
  FishMet4 = data.frame()
  
  ##-- Determine remaining fishing effort per ecoregion
  for(iReg in c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
    print(paste0(iReg, ". Scenario ", iSC))
    CsquareInfoReg = subset(CsquareInfo, Ecoregion == iReg)
    
    FishMet <- as.data.frame(merge(CsquareInfoReg, VMSdatMet,by="csquares")) # status quo
    SA = grep("surface_SA", colnames(FishMet))
    TLW = grep("total_weight", colnames(FishMet))
    FishMet = FishMet[,c(1:17, SA, TLW)]
    print(paste0("Total SA: ", round(sum(colSums(FishMet[,18:67])))))
    print(paste0("Total TLW: ", round(sum(colSums(FishMet[,68:117])))))
    
    FishMet2 <- FishMet
    FishMet2[,18:67] <- FishMet2[,18:67]*(FishMet2[,paste0(iSC, "pct")]) # removed effort
    FishMet2[,18:67][is.na(FishMet2[,18:67])] <- 0
    print(paste0("Total removed SA: ", round(sum(colSums(FishMet2[,18:67])))))
    
    FishMet3 <- FishMet
    FishMet3[,18:67] <- FishMet3[,18:67]*(1-FishMet3[,paste0(iSC, "pct")]) # remaining fishing
    FishMet3[,18:67][is.na(FishMet3[,18:67])] <- 0
    print(paste0("Total remaining SA: ", round(sum(colSums(FishMet3[,18:67])))))
    
    Sc0_FI <- colSums(FishMet[,18:67]) # This is what is done originally
    ScX_FIremoved <- colSums(FishMet2[,18:67]) # This is what is removed due to closures
    ScX_FIremain <- colSums(FishMet3[,18:67]) # This is what remains after closures (no displacement)
    
    ##-- Determine (old) LpUE per c-square per gear/year
    LpUE <- cbind(FishMet$csquares, FishMet[,68:117] / FishMet[,18:67], FishMet[paste0(iSC, "pct")])
    LpUE[is.na(LpUE)] <- 0 # if no kg and no effort is registered, this yields a NaN, which is set to 0
    a <- c("csquares", gsub(x=colnames(LpUE)[2:51], pattern="total_weight", replacement="LpUE"), paste0(iSC, "pct"))
    colnames(LpUE) <- a
    
    ##-- Determine for each met/year combi the LpUE weights to distribute removed effort
    LpUE2 <- LpUE
    LpUE2$MPA_PCT <- ifelse(LpUE2[[paste0(iSC, "pct")]] == 1, 0, 1) # Identify c-squares that are completely closed
    LpUE2[,2:51] <- LpUE2[,2:51] * LpUE2$MPA_PCT # remove LpUE data for completely closed c-squares 
    LpUE2[,2:51] <- sweep(x=LpUE2[,2:51], MARGIN = 2, STATS = colSums(LpUE2[,2:51]), FUN="/") # determine weighting factor as % of sum(LpUE)
    LpUE2[is.na(LpUE2)] <- 0
    
    ##-- Determine for each met/year combi the total FE that should be added to the remaining effort
    LpUE3 <- LpUE2
    LpUE3[,2:51] <- sweep(x=LpUE3[,2:51], MARGIN = 2, STATS = ScX_FIremoved, FUN="*")
    ScX_FIadded <- colSums(LpUE3[,2:51])
    
    ##-- Determine redistribution of fishing effort
    FishMetD <- cbind(FishMet3$csquares, (FishMet3[,18:67] + LpUE3[,2:51]))
    colnames(FishMetD)[1] <- "csquares"
    print(paste0("Total new SA: ", round(sum(colSums(FishMetD[,2:51])))))
    ScX_newFI <- colSums(FishMetD[,2:51])
    
    ##-- Determine new landings using old LpUE and new effort data.
    landings = FishMetD[,2:51] * LpUE[,2:51]
    a <- c(gsub(x=colnames(landings), pattern="surface_SA", replacement="total_weight"))
    colnames(landings) <- a
    FishMetD = cbind(FishMetD, landings)
    print(paste0("Total new TLW: ", round(sum(colSums(FishMetD[,52:101])))))
    
    ##-- add some other columns
    FishMetD$inf_intercept <- FishMet$inf_intercept
    FishMetD$inf_slope <- FishMet$inf_slope
    FishMetD$epi_intercept <- FishMet$epi_intercept
    FishMetD$epi_slope <- FishMet$epi_slope
    FishMetD$area_sqkm = FishMet$area_sqkm
    
    FishMet4 = rbind(FishMet4, FishMetD)
  } # end iReg-loop
  
  ##-- Turn updated SA into SAR
  idx = grep("surface_SA", colnames(FishMet4))
  SAdat = as.data.table(FishMet4)[,..idx]
  SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= function(x){ x/FishMet4$area_sqkm}))
  names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
  idy = setdiff(1:ncol(FishMet4), idx)
  FishMet4 = cbind(as.data.table(FishMet4)[,..idy], SAdat2)
  
  ##-- Save final vmsdatfile used
  saveRDS(FishMet4, file=paste0(outPath, "FishMet4_ScC_", iSC, ".rds"))
  
  ##-- Determine state
  Period = 2016:2020
  source(paste0(scriptsPath, "HSF_Scenarios_inf.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScC",iSC, "_inf.rds"))

  source(paste0(scriptsPath, "HSF_Scenarios_epi.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScC", iSC, "_epi.rds"))
}# end iSC-loop

#-#-# Scenario D: TLW-EFF ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
CsquareInfo$S1pct <- ifelse(CsquareInfo$S1pct > 0.95, 1, CsquareInfo$S1pct)
CsquareInfo$S2pct <- ifelse(CsquareInfo$S2pct > 0.95, 1, CsquareInfo$S2pct)

VMSdatMet = VMSdatMet[,c(1,3:102)]
VMSdatMet[is.na(VMSdatMet)] = 0
source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo

for(iSC in c("S1", "S2")){
  FishMet4 = data.frame()
  
  ##-- Determine remaining fishing effort per ecoregion
  for(iReg in c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
    print(paste0(iReg, ". Scenario ", iSC))
    CsquareInfoReg = subset(CsquareInfo, Ecoregion == iReg)
    
    FishMet <- as.data.frame(merge(CsquareInfoReg, VMSdatMet,by="csquares")) # status quo
    SA = grep("surface_SA", colnames(FishMet))
    TLW = grep("total_weight", colnames(FishMet))
    FishMet = FishMet[,c(1:17, SA, TLW)]
    print(paste0("Total SA: ", round(sum(colSums(FishMet[,18:67])))))
    print(paste0("Total TLW: ", round(sum(colSums(FishMet[,68:117])))))
    
    FishMet2 <- FishMet
    FishMet2[,68:117] <- FishMet2[,68:117]*(FishMet2[,paste0(iSC, "pct")]) # removed landings
    FishMet2[,68:117][is.na(FishMet2[,68:117])] <- 0
    print(paste0("Total removed TLW: ", round(sum(colSums(FishMet2[,68:117])))))
    
    FishMet3 <- FishMet
    FishMet3[,68:117] <- FishMet3[,68:117]*(1-FishMet3[,paste0(iSC, "pct")]) # remaining landings
    FishMet3[,68:117][is.na(FishMet3[,68:117])] <- 0
    print(paste0("Total remaining TLW: ", round(sum(colSums(FishMet3[,68:117])))))
    
    Sc0_TLW <- colSums(FishMet[,68:117]) # This is what is landed originally
    Sc0_FI <- colSums(FishMet[,18:67]) # This is effort originally
    ScX_TLWremoved <- colSums(FishMet2[,68:117]) # This is what is removed due to closures
    ScX_TLWremain <- colSums(FishMet3[,68:117]) # This is what remains after closures (no displacement)
    
    ##-- Determine LpUE per c-square per gear/year
    LpUE <- cbind(FishMet$csquares, FishMet[,68:117] / FishMet[,18:67], FishMet[paste0(iSC, "pct")])
    LpUE[is.na(LpUE)] <- 0 # if no kg and no effort is registered, this yields a NaN, which is set to 0
    a <- c("csquares", gsub(x=colnames(LpUE)[2:51], pattern="total_weight", replacement="LpUE"), paste0(iSC, "pct"))
    colnames(LpUE) <- a
    
    ##-- Determine for each met/year combi the effort-based weights to distribute removed landings
    EFF <- FishMet
    EFF[,18:67] = EFF[,18:67]*(1-EFF[,paste0(iSC, "pct")]) # determine remaining effort (after fisheries closures)
    EFF[,18:67][is.na(EFF[,18:67])] = 0
    EFF[,18:67] <- sweep(x=EFF[,18:67], MARGIN = 2, STATS = colSums(EFF[,18:67]), FUN="/") # determine weighting factor as % of sum(effort)
    EFF[is.na(EFF)] <- 0
    
    ##-- Determine for each met/year combi the new landings that should be added, relative to old effort
    addLand <- EFF
    addLand[,18:67] <- sweep(x=addLand[,18:67], MARGIN = 2, STATS = ScX_TLWremoved, FUN="*")
    ScX_TLWadded <- colSums(addLand[,18:67])
    
    ##-- Determine the new fishing effort, using old LpUE
    FishMetD <- cbind(FishMet3$csquares, (FishMet3[,68:117] + addLand[,18:67]))
    colnames(FishMetD)[1] <- "csquares"
    print(paste0("Total new TLW: ", round(sum(colSums(FishMetD[,2:51])))))
    
    FishMetD2 = cbind(FishMetD$csquares, (FishMetD[,2:51] / LpUE[,2:51]))
    FishMetD2[is.na(FishMetD2)] = 0
    #FishMetD2[sapply(FishMetD2, is.infinite)] <- 0
    colnames(FishMetD2)[1] <- "csquares"
    print(paste0("Total new SA: ", round(sum(colSums(FishMetD2[,2:51])))))
    a <- c("csquares", gsub(x=colnames(FishMetD2)[2:51], pattern="total_weight", replacement="surface_SA"))
    colnames(FishMetD2) <- a
    FishMetD2 = cbind(FishMetD2, FishMetD[,2:51])
    
    ScX_newFI <- colSums(FishMetD2[,2:51])
    FishMetD2$inf_intercept <- FishMet$inf_intercept
    FishMetD2$inf_slope <- FishMet$inf_slope
    FishMetD2$epi_intercept <- FishMet$epi_intercept
    FishMetD2$epi_slope <- FishMet$epi_slope
    FishMetD2$area_sqkm = FishMet$area_sqkm
    
    FishMet4 = rbind(FishMet4, FishMetD2)
  } # end iReg-loop
  
  ##-- Turn updated SA into SAR
  idx = grep("surface_SA", colnames(FishMet4))
  SAdat = as.data.table(FishMet4)[,..idx]
  SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= function(x){ x/FishMet4$area_sqkm}))
  names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
  idy = setdiff(1:ncol(FishMet4), idx)
  FishMet4 = cbind(as.data.table(FishMet4)[,..idy], SAdat2)
  
  ##-- Save final vmsdatfile used
  saveRDS(FishMet4, file=paste0(outPath, "FishMet4_ScD_", iSC, ".rds"))
  
  ##-- Determine state
  Period=2016:2020
  source(paste0(scriptsPath, "HSF_Scenarios_inf.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScD",iSC, "_inf.rds"))

  source(paste0(scriptsPath, "HSF_Scenarios_epi.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScD", iSC, "_epi.rds"))
}# end iSC-loop


#-#-# Scenario E: TLW-LPUE ----
CsquareInfo <- st_drop_geometry(readRDS(file=paste0(datPath, "Regiontot2.rds")))
VMSdatMet <- readRDS(file=paste0(datPath, "FD_capped.rds"))
CsquareInfo$S1pct <- ifelse(CsquareInfo$S1pct > 0.95, 1, CsquareInfo$S1pct)
CsquareInfo$S2pct <- ifelse(CsquareInfo$S2pct > 0.95, 1, CsquareInfo$S2pct)

VMSdatMet = VMSdatMet[,c(1,3:102)]
VMSdatMet[is.na(VMSdatMet)] = 0
source(paste0(scriptsPath, "Impact_continuous_longevity.R")) # get RBS-functions from FBIT repo

for(iSC in c("S1", "S2")){
  FishMet4 = data.frame()
  
  ##-- Determine remaining fishing effort per ecoregion
  for(iReg in c("Greater North Sea", "Baltic Sea", "Celtic Seas", "Bay of Biscay and the Iberian Coast")){
    print(paste0(iReg, ". Scenario ", iSC))
    CsquareInfoReg = subset(CsquareInfo, Ecoregion == iReg)
    
    FishMet <- as.data.frame(merge(CsquareInfoReg, VMSdatMet,by="csquares")) # status quo
    SA = grep("surface_SA", colnames(FishMet))
    TLW = grep("total_weight", colnames(FishMet))
    FishMet = FishMet[,c(1:17, SA, TLW)]
    print(paste0("Total SA: ", round(sum(colSums(FishMet[,18:67])))))
    print(paste0("Total TLW: ", round(sum(colSums(FishMet[,68:117])))))
    
    FishMet2 <- FishMet
    FishMet2[,68:117] <- FishMet2[,68:117]*(FishMet2[,paste0(iSC, "pct")]) # removed landings
    FishMet2[,68:117][is.na(FishMet2[,68:117])] <- 0
    print(paste0("Total removed TLW: ", round(sum(colSums(FishMet2[,68:117])))))
    
    FishMet3 <- FishMet
    FishMet3[,68:117] <- FishMet3[,68:117]*(1-FishMet3[,paste0(iSC, "pct")]) # remaining landings
    FishMet3[,68:117][is.na(FishMet3[,68:117])] <- 0
    print(paste0("Total remaining TLW: ", round(sum(colSums(FishMet3[,68:117])))))
    
    Sc0_TLW <- colSums(FishMet[,68:117]) # This is what is landed originally
    Sc0_FI <- colSums(FishMet[,18:67]) # This is effort originally
    ScX_TLWremoved <- colSums(FishMet2[,68:117]) # This is what is removed due to closures
    ScX_TLWremain <- colSums(FishMet3[,68:117]) # This is what remains after closures (no displacement)
    
    ##-- Determine LpUE per c-square per gear/year
    LpUE <- cbind(FishMet$csquares, FishMet[,68:117] / FishMet[,18:67], FishMet[paste0(iSC, "pct")])
    LpUE[is.na(LpUE)] <- 0 # if no kg and no effort is registered, this yields a NaN, which is set to 0
    #LpUE[sapply(LpUE, is.infinite)] <- 0 # If no effort is registered, but there are kg, this yields an inf, which is set to 0
    a <- c("csquares", gsub(x=colnames(LpUE)[2:51], pattern="total_weight", replacement="LpUE"), paste0(iSC, "pct"))
    colnames(LpUE) <- a
    
    ##-- Determine for each met/year combi the LpUE weights to distribute removed landings
    LpUE2 <- LpUE
    LpUE2$MPA_PCT <- ifelse(LpUE2[[paste0(iSC, "pct")]] == 1, 0, 1) # Identify c-squares that are completely closed
    LpUE2[,2:51] <- LpUE2[,2:51] * LpUE2$MPA_PCT # remove LpUE data for completely closed c-squares 
    LpUE2[,2:51] <- sweep(x=LpUE2[,2:51], MARGIN = 2, STATS = colSums(LpUE2[,2:51]), FUN="/") # determine weighting factor as % of sum(LpUE)
    LpUE2[is.na(LpUE2)] <- 0
    
    ##-- Determine for each met/year combi the new landings that should be added, relative to old LpUE.
    LpUE3 <- LpUE2
    LpUE3[,2:51] <- sweep(x=LpUE3[,2:51], MARGIN = 2, STATS = ScX_TLWremoved, FUN="*")
    ScX_TLWadded <- colSums(LpUE3[,2:51])
    
    ##-- Determine the new fishing effort, using old LpUE
    FishMetD <- cbind(FishMet3$csquares, (FishMet3[,68:117] + LpUE3[,2:51]))
    colnames(FishMetD)[1] <- "csquares"
    print(paste0("Total new TLW: ", round(sum(colSums(FishMetD[,2:51])))))
    
    FishMetD2 = cbind(FishMetD$csquares, (FishMetD[,2:51] / LpUE[,2:51]))
    FishMetD2[is.na(FishMetD2)] = 0
    #FishMetD2[sapply(FishMetD2, is.infinite)] <- 0
    colnames(FishMetD2)[1] <- "csquares"
    print(paste0("Total new SA: ", round(sum(colSums(FishMetD2[,2:51])))))
    a <- c("csquares", gsub(x=colnames(FishMetD2)[2:51], pattern="total_weight", replacement="surface_SA"))
    colnames(FishMetD2) <- a
    FishMetD2 = cbind(FishMetD2, FishMetD[,2:51])
    
    ScX_newFI <- colSums(FishMetD2[,2:51])
    FishMetD2$inf_intercept <- FishMet$inf_intercept
    FishMetD2$inf_slope <- FishMet$inf_slope
    FishMetD2$epi_intercept <- FishMet$epi_intercept
    FishMetD2$epi_slope <- FishMet$epi_slope
    FishMetD2$area_sqkm = FishMet$area_sqkm
    
    FishMet4 = rbind(FishMet4, FishMetD2)
  } # end iReg-loop
  
  ##-- Turn updated SA into SAR
  idx = grep("surface_SA", colnames(FishMet4))
  SAdat = as.data.table(FishMet4)[,..idx]
  SAdat2 = as.data.table(apply(SAdat, MARGIN=2, FUN= function(x){ x/FishMet4$area_sqkm}))
  names(SAdat2) = gsub("surface_SA", replacement= "surface_sar", x=names(SAdat2))
  idy = setdiff(1:ncol(FishMet4), idx)
  FishMet4 = cbind(as.data.table(FishMet4)[,..idy], SAdat2)
  
  ##-- Save final vmsdatfile used
  saveRDS(FishMet4, file=paste0(outPath, "FishMet4_ScE_", iSC, ".rds"))
  
  # ##-- Determine state 
  Period = 2016:2020
  source(paste0(scriptsPath, "HSF_Scenarios_inf.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScE",iSC, "_inf.rds"))

  source(paste0(scriptsPath, "HSF_Scenarios_epi.R")) # adjusted from FBIT repo
  saveRDS(State_reg, paste0(outPath, "State_reg_ScE", iSC, "_epi.rds"))
}# end iSC-loop

