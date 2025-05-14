#Script to Analyse Background Respirometry of PumpResp 
#Experiment in Laukaa LUKE facility Summer 2023 using pumpresp set up in 16 chambers (Shrek2) 
#Lucy Cotgrove, Jenni Prokkola, Evan Sala
#Using fishresp package 

# This Code calculates 1)Pre and post background respiration of each chamber, 
# and creates a file of the mean from three slops before and after SMR, then 2) SMR for each fish

#clear data & gc
rm(list=ls(all=TRUE));gc() 

#load packages
library(tidyverse)
library(stringi)
library(data.table)
library(respirometry)
library(respR)
library(FishResp)
library(dplyr)

#working directory
setwd("")
path <- getwd()

#list files

fishdat <- list.files(patt="PR[1-4]_fish_B[1-4]")
predat <- list.files(patt="PR[1-4]_preblank_B[1-4]")
postdat <- list.files(patt="PR[1-4]_postblank_B[1-4]")


#This file is created from the manual inspection of slopes further in this code, and the code is re-run with this file loaded after inspection and deciding to remove certain phases.
#have included as a manually made df
PhasesToRemove <- df <- data.frame(
  Ind = c(
    "PR1.B4.1.1", "PR1.B4.1.1", "PR1.B1.1.1", "PR1.B1.1.1", "PR1.B1.1.1", 
    "PR1.B1.1.3", "PR1.B1.1.3", "PR1.B2.1.2", "PR1.B2.1.2", "PR1.B2.1.2", 
    "PR1.B2.1.3", "PR1.B2.1.3", "PR1.B2.1.4", "PR2.B1.2.1", "PR2.B1.2.1", 
    "PR2.B1.2.3", "PR2.B1.2.4", "PR2.B2.2.2", "PR2.B2.2.2", "PR2.B2.2.4", 
    "PR4.B7.4.3", "PR4.B7.4.3", "PR4.B7.4.3", "PR4.B7.4.3", "PR4.B3.4.3", 
    "PR4.B3.4.3", "PR3.B7.3.3", "PR3.B7.3.3", "PR3.B7.3.3", "PR3.B7.3.4", 
    "PR3.B7.3.4", "PR3.B8.3.2", "PR3.B1.3.3", "PR3.B1.3.4", "PR3.B1.3.4", 
    "PR3.B1.3.2", "PR3.B4.3.2", "PR3.B2.3.2", "PR3.B2.3.2", "PR3.B2.3.2", 
    "PR3.B2.3.2", "PR3.B2.3.3", "PR3.B2.3.3", "PR3.B2.3.3"
  ),
  Phase = c(
    "M12", "M39", "M37", "M35", "M46", "M43", "M36", "M31", "M35", "M41", 
    "M27", "M39", "M38", "M39", "M45", "M50", "M50", "M34", "M35", "M41", 
    "M35", "M36", "M40", "M45", "M24", "M40", "M31", "M36", "M47", "M46", 
    "M37", "M28", "M35", "M33", "M35", "M33", "M25", "M11", "M13", "M31", 
    "M42", "M32", "M35", "M18"
  ),
  stringsAsFactors = FALSE
)

#Info r.e chambers needed for SMR calculation
chamberdat <- read.csv("combinedchamberdata.csv")
chamberdat$Chamber.No <- paste0("CH", chamberdat$chamber)

#making the Batch/pumpresp list, accounting for the lost data

batch <- c("B1","B2","B4")
pumpresp <- c("PR1", "PR2", "PR3", "PR4")
batchnumbers <- expand.grid( batch = batch, pumpresp = pumpresp)
batch <- c("B5","B6","B7","B8")
pumpresp <- c("PR3", "PR4")
batchnumbers1 <- expand.grid( batch = batch, pumpresp = pumpresp)
batch <- c("B3")
pumpresp <- c("PR3", "PR4")
batchnumbers2 <- expand.grid( batch = batch, pumpresp = pumpresp)
batchnumbers <- rbind(batchnumbers, batchnumbers1, batchnumbers2)

##############
##  PreBG   ##
##############

predatBG <-list.files(pattern = "BGpre_B")

listBG <- list()

rm(j)
for(j in 1:length(predatBG)){
  
  #load chamber data for batch and pumpresp
  chamber <- subset(chamberdat, chamberdat$Batch == (str_extract(predatBG[j], pattern = "B[1-8]")) & (chamberdat$pumpresp == str_extract(predatBG[j], pattern = "PR[1-4]")))
  
  #load info from fishresp function
  info <- input.info(ID = c(chamber$Fish_ID),
                     Mass = c(chamber$Weight_g),
                     Volume = c(chamber$Volume_g),
                     DO.unit = "mg/L")  
  
  #import test from fishresp, convert from tibble to df, generates oxygen data
  listBG[[j]] <- as.data.frame(import.test(file = predatBG[[j]],
                                           info.data = info,
                                           logger = "FishResp",
                                           n.chamber = 4,
                                           meas.to.wait = 210, #how long to wait post measure period start in seconds 
                                           plot.temperature = F,
                                           plot.oxygen = F))
  #add data from extract
  listBG[[j]]$file <-  predatBG[[j]]
  listBG[[j]]$Batch <- str_extract(predatBG[j], pattern = "B[1-8]") 
  listBG[[j]]$pumpresp <- str_extract(predatBG[j], pattern = "PR[1-4]")
  listBG[[j]]$measure <- str_extract(predatBG[j], pattern = "M[1-6]")
  
}

allpreBG <- rbindlist(listBG)


#remove bad pre values B5, 2.3 pre B3 2.3 pre
allpreBG <-  allpreBG %>% filter((Chamber.No !="CH3" & Batch != "3" & pumpresp != "PR2") & (Chamber.No !="CH3" & Batch != "5" & pumpresp != "PR2"))%>% filter(measure == "M1" | measure == "M2" | measure == "M3")

Bs <- unique(batchnumbers$batch)

#calculates mean of pre background respiration 
for(i in 1:length(Bs)){
  
  allBG <- subset(allpreBG, Batch == as.character(Bs[i]))
  
  
  meanpreBG1 <- allpreBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH1", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpreBG2 <- allpreBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH2", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpreBG3 <- allpreBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH3", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpreBG4 <- allpreBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH4", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  
meanpreBG <- rbind(meanpreBG1, meanpreBG2, meanpreBG3, meanpreBG4)
  
write.table(meanpreBG, "meanPreBG.txt", sep = '\t')

  }
  
##############
##  PostBG  ##
##############


postdatBG <- list.files(pattern = "BGpost_B")
listBG <- list()
rm(j)

for(j in 1:length(postdatBG)){
  tryCatch({
    #load chamber data for batch and pumpresp
    chamber <- subset(chamberdat, chamberdat$Batch == (str_extract(postdatBG[j], pattern = "B[1-8]")) & (chamberdat$pumpresp == str_extract(postdatBG[j], pattern = "PR[1-4]")))
  
    #load info from fishresp function
  info <- input.info(ID = c(chamber$Fish_ID),
                     Mass = c(chamber$Weight_g),
                     Volume = c(chamber$Volume_g),
                     DO.unit = "mg/L")  
  
  #import test from fishresp, convert from tibble to df, generates oxygen data
  listBG[[j]] <- as.data.frame(import.test(file = postdatBG[[j]],
                                           info.data = info,
                                           logger = "FishResp",
                                           n.chamber = 4,
                                           meas.to.wait = 210,#how long to wait post measure period start in seconds
                                           plot.temperature = F,
                                           plot.oxygen = F))
  #add data from extract
  listBG[[j]]$file <-  postdatBG[[j]]
  listBG[[j]]$Batch <- str_extract(postdatBG[j], pattern = "B[1-8]") 
  listBG[[j]]$pumpresp <- str_extract(postdatBG[j], pattern = "PR[1-4]")
  listBG[[j]]$measure <- str_extract(postdatBG[j], pattern = "M[1-6]")
  

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
}

allpostBG <- rbindlist(listBG)

#remove bad values B5, 2.3 pre B3 2.3 pre

allpostBG <-  allpostBG %>% filter((Chamber.No !="CH3" & Batch != "3" & pumpresp != "PR2") & (Chamber.No !="CH3" & Batch != "5" & pumpresp != "PR2"))%>% filter(measure == "M1" | measure == "M2" | measure == "M3")


Bs <- unique(batchnumbers$batch)
#calculates mean of POST background respiration 
for(i in 1:length(Bs)){
  
  allBG <- subset(allpostBG, Batch == as.character(Bs[i]))
  
  
  meanpostBG1 <- allpostBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH1", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpostBG2 <- allpostBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH2", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpostBG3 <- allpostBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH3", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  meanpostBG4 <- allpostBG %>% group_by(Batch, Time) %>% summarise(delta.O2 = mean(delta.O2, na.rm = T), Temp = mean(Temp, na.rm = T), Init.O2 = mean(Init.O2, na.rm = T), O2 = mean(O2, na.rm = T)) %>% mutate(Chamber.No = "CH4", measure = "M1", Test = "test") %>% dplyr::select(Chamber.No, Test, Time, Init.O2,Temp, O2, delta.O2, Batch)
  
  meanpostBG <- rbind(meanpostBG1, meanpostBG2, meanpostBG3, meanpostBG4)
  
  write.table(meanpostBG, "meanPostBG.txt", sep = '\t')

}



##############
##  SMR     ##
##############

R2list <- list()
SMRfull <- list()
SMRfullsum <- list()

rm(i)


# test: PR4.B1, PR1.B3, PR4.B2, PR4.B1, 
for (i in 1:nrow(batchnumbers)){
  tryCatch({
    #subset chamber data to match batchnumbers
   chamber <- subset(chamberdat, chamberdat$Batch == (as.character(batchnumbers[i,][,1])) & (chamberdat$pumpresp == (as.character(batchnumbers[i,][,2]))))
  
   #make dataframe of NAs to join for those with missing chambers
   missingchamber <- data.frame(Chamber_ID = c(paste0(unique(chamber$FishResp_ID),".1"),
                                                              paste0(unique(chamber$FishResp_ID),".2"),
                                                              paste0(unique(chamber$FishResp_ID),".3"),
                                                              paste0(unique(chamber$FishResp_ID),".4")),
                                               Volume_g = rep(NA, 4),
                                               pumpresp = rep(unique(chamber$pumpresp), 4),
                                               chamber = c(1, 2, 3, 4),
                                               FishResp_ID = rep(unique(chamber$FishResp_ID), 4),
                                               Sampling.date = rep(NA, 4),
                                               Temperature = rep(unique(chamber$Temperature), 4),
                                               Batch = rep(unique(chamber$Batch), 4),
                                               Sample_ID = rep(NA, 4),
                                               DNA_well = rep(NA, 4),
                                               Weight_g = rep(NA, 4),
                                               Total_length_mm = rep(NA, 4),
                                               Fish_ID = rep(NA, 4),
                                               Chamber.No = c("CH1", "CH2", "CH3", "CH4"))
   #fill chamber df with NAs if missing data 
   chamber <- plyr::rbind.fill(chamber, missingchamber)
    chamber <- chamber %>% distinct(Chamber_ID, .keep_all = T) %>% arrange(Chamber_ID)
    
    #ID which chmabers are missing
    missing <- chamber[is.na(chamber$Volume_g), ]$Chamber.No
    
    #input info for fishresp package
    info <- input.info(ID = c(chamber$Fish_ID),
                                     Mass = c(chamber$Weight_g),
                                     Volume = c(chamber$Volume_g),
                                     DO.unit = "mg/L")
    # Import pre-test data
    pre <- subset(meanpreBG, Batch == as.character(batchnumbers[i, ][,1]))
    # Import post-test data
    post <- subset(meanpostBG, Batch == as.character(batchnumbers[i, ][,1]))
    # get fish data
    fishname <- paste0(unique(chamber$pumpresp), "_fish_", unique(chamber$Batch), ".txt")
    # Import all actual SMR data.
    SMR.raw <- import.meas(file = fishname,
                            info.data = info,
                            logger = "FishResp",
                            n.chamber = 4,
                            meas.to.wait = 210,
                            meas.to.flush = 120,
                            date.format = "YMD",
                            start.measure = "15:00:00",
                            stop.measure = "07:00:00",
                            plot.temperature = F,
                            plot.oxygen = F)
    #############################################################
    #-----------------------------------------------------------#
    #------------ STEP 3: BACK RESP CORRECTION -----------------#
    #-----------------------------------------------------------#
    #############################################################
    SMR.clean <- correct.meas(info.data = info,
                              pre.data =  pre,
                              post.data = post,
                              meas.data = SMR.raw,
                              method = "linear")
    #########################################################
    #-------------------------------------------------------#
      #-------- STEP 4: FILTERING AND SMOOTHING --------------#
    #-------------------------------------------------------#
    #########################################################
    #from visual inspection, if Ind is in list, then smooth using filter
      #PR1 B2 & B4 smooth then reanalyse, PR3 B1, B2, B3 and B4 smooth: noisy data but trends look good.
      #CHabot et al 2021 recommend 29 frame moving average
      
smoothlist <- c("PR1.B2.1.3", "PR1.B2.1.4", "PR1.B4.1.1", "PR3.B1.3.2",
                "PR3.B1.3.3","PR3.B2.3.2", "PR3.B2.3.3", "PR3.B3.3.1", 
                "PR3.B4.3.1", "PR3.B4.3.2", #smoothing because low number of slopes
 "PR3.B7.3.3", "PR3.B5.3.4", "PR2.B1.2.4", "PR2.B2.2.4",
 "PR1.B4.1.2", "PR1.B1.1.3")

dat <- list()

rm(k)

#smooth data slopes where needed
for(k in 1:length(unique(SMR.clean$Ind))){
  dat[[k]] <- as.data.frame(SMR.clean) %>% filter(Ind == as.character(unique(SMR.clean$Ind)[k]))
  
  if(nrow(dat[[k]]) == 0){
    dat[[k]] <- dat[[k]]
  } else if(unique(dat[[k]]$Ind) %in% smoothlist & nrow(dat[[k]]) != 0){
    
  dat[[k]]$O2.correct <- zoo::rollapplyr(dat[[k]]$O2.correct, 29, mean, na.rm = TRUE, partial = T, fill = NA)
  dat[[k]]$O2 <- zoo::rollapplyr(dat[[k]]$O2, 29, mean, na.rm = TRUE, partial = T, fill = NA)
  dat[[k]]$BR <- zoo::rollapplyr(dat[[k]]$BR, 29, mean, na.rm = TRUE, partial = T, fill = NA)
  } else  {
  dat[[k]] <- dat[[k]]
  }
  
  }

#tidy data up for imaging
dat <- dat[sapply(dat, nrow)>0]
SMR.clean <- rbindlist(dat)
SMR.clean <- na.omit(SMR.clean)
SMR.clean <- droplevels(SMR.clean)

SMR.clean <- SMR.clean %>% group_by(Chamber.No, Phase) %>% slice(120:(n()-30))


#########################################################
#-------------------------------------------------------#
#------------------- STEP 5: IMAGES --------------------#
#-------------------------------------------------------#
#########################################################


#can use these in the SMR loop to generate images for quality control - RENAME WITH YOUR FILEPATH


png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_temp.png"),
    width=1800, height=1050)
temp <- QC.meas(SMR.clean, "Temperature")
temp$main <- fishname
print(temp)
dev.off()

png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_phas.png"),
    width=1800, height=1050)
phas <- QC.meas(SMR.clean, "Total.O2.phases")
phas$main <- fishname
print(phas)
dev.off()

png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_corrphas.png"),
    width=1800, height=1050)
corr.phas <- QC.meas(SMR.clean, "Corrected.O2.phases")
corr.phas$main <- fishname
print(corr.phas)
dev.off()

png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_chamb.png"),
    width=1800, height=1050)
chamb <- QC.meas(SMR.clean, "Total.O2.chambers")
chamb$main <- fishname
print(chamb)
dev.off()

png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_corrchamb.png"),
    width=1800, height=1050)
corr.chamb <- QC.meas(SMR.clean, "Corrected.O2.chambers")
corr.chamb$main <- fishname
print(corr.chamb)  
dev.off()

png(file=paste0("qc_",unique(chamber$pumpresp), "_fish_", unique(chamber$Batch),"_activty.png"),
    width=1800, height=1050)
act<- QC.activity(SMR.clean, compare = T) # MR + BG vs just MR
dev.off()



#########################################################
#-------------------------------------------------------#
#------------ STEP 4: SLOPE EXTRACTION -----------------#
#-------------------------------------------------------#
#########################################################


#"all" means individual slopes for each Measurement phase, SMRslope is ordered from low to high in slope, so when using MLND from CHabot 2016, thats why cl is 1,1,1,1,2,2,2,2,3,3,3 etc because it's numbering which slope it is using (not in order of slope measurement)
SMR.slope <- extract.slope(SMR.clean,
                           method = "all",
                           r2 = 0.95)

#remove slopes manually
SMR.slope$Phase<-as.character(SMR.slope$Phase)
SMR.slope <- anti_join(SMR.slope, PhasesToRemove, by = c("Ind"="Ind", "Phase"="Phase")) #removes phases that are "bad" listed in the csv at the beginning.
SMR.slope$Phase<-ordered(SMR.slope$Phase)

#needs all channels for calculate MR from fishresp to work
SMR <- calculate.MR(SMR.slope,
                    density = 1000,
                    plot.BR = F,
                    plot.MR.abs = F,
                    plot.MR.mass = F)
rm(m)

SMRlist <- list()
SMRCH <- list()
fs <- list()



######### FUNCTION FROM CHABOT 2016 #########
require(mclust)

calcSMR = function(Y, q=c(0.1,0.15,0.2,0.25,0.3), G=1:4){
  u = sort(Y)
  the.Mclust <- Mclust(Y,  G=G)
  cl <- the.Mclust$classification
  # sometimes, the class containing SMR is not called 1
  # the following presumes that when class 1 contains > 10% of cases, 
  # it contains SMR, otherwise we take class 2
  cl2 <- as.data.frame(table(cl))
  cl2$cl <- as.numeric(levels(cl2$cl))
  valid <- cl2$Freq>=0.1*length(cl)  
  the.cl <- min(cl2$cl[valid])
  left.distr <- Y[the.Mclust$classification==the.cl]
  mlnd = the.Mclust$parameters$mean[the.cl]
  CVmlnd = sd(left.distr)/mlnd * 100
  quant=quantile(Y, q)
  low10=mean(u[1:10])
  low10pc = mean(u[6:(5 + round(0.1*(length(u)-5)))])
  # remove 5 outliers, keep lowest 10% of the rest, average
  # Herrmann & Enders 2000
  return(list(mlnd=mlnd, quant=quant, low10=low10, low10pc=low10pc,
              cl=cl, CVmlnd=CVmlnd))
}
 ######### this function returns all the types of SMR for comparison i.e. mlnd, low10, quants ######

#uses the chabot function to calculate SMR and extract info from chabot fucntion which is needed to check data (i.e. MLND features, total phases used in SMR calc) 
for(m in 1:length(unique(SMR$Chamber.No))){
  SMRCH[[m]] <- SMR %>% filter(Chamber.No == unique(SMR$Chamber.No)[[m]])
  SMRchab <- calcSMR(SMRCH[[m]]$MR.abs, q=c(0.1,0.15,0.2,0.25,0.3), G=4)
  SMRCH[[m]]$cl <- SMRchab$cl
  SMRCH[[m]]$No.mlnd =   as.numeric(names(SMRchab$mlnd[1]))
  SMRlist[[m]] <- data.frame(Chamber.No = unique(SMRCH[[m]]$Chamber.No),
                                Ind = unique(SMRCH[[m]]$Ind),
                                PhaseUsed = as.numeric(length(which(SMRchab$cl ==as.numeric(names(SMRchab$mlnd[1]))))),
                              No.mlnd =   as.numeric(names(SMRchab$mlnd[1])),
                             TotalPhase = as.numeric(length(SMRchab$cl)),
                                MLND = SMRchab$mlnd,
                                MLND.cv = SMRchab$CVmlnd,
                                low10 = SMRchab$low10,
                                low10pc = SMRchab$low10pc,
                                quant10 = SMRchab$quant[1],
                                quant15 = SMRchab$quant[2],
                                quant20 = SMRchab$quant[3],
                                quant25 = SMRchab$quant[4],
                                quant30 = SMRchab$quant[5])
  
  # run MLND then check which slopes are being used - only linear 
  fs[[m]] <- SMRCH[[m]] %>% filter(cl == as.numeric(names(SMRchab$mlnd[1])))
  
}

SMRsummary <- data.table::rbindlist(SMRlist)
SMR <- rbindlist(SMRCH)

SMRfull[[i]] <- SMR

#all slopes used in MLND saved in csv for analysis
write.csv(SMRfull[[i]], paste0("AllSlopesMLND_",unique(chamber$pumpresp),"_", unique(chamber$Batch),".csv"))
SMRfullsum[[i]] <- SMRsummary #all the summary values


# run MLND then check which slopes are being used - we only want linear. 
FinalSlopes <- rbindlist(fs)

FinalSlopeSUM <- FinalSlopes %>% dplyr::select(Chamber.No, Phase)
FinalPlotDat <- SMR.clean %>% filter(interaction(Phase, Chamber.No) %in% interaction(FinalSlopeSUM$Phase, FinalSlopeSUM$Chamber.No))
FinalPlotDat <- anti_join(FinalPlotDat, PhasesToRemove, by = c("Ind"="Ind", "Phase"="Phase"))

#check visually
png(file=paste0("FinalSlopesMLND_",unique(chamber$pumpresp), unique(chamber$Batch),"_corrphas.png"),
    width=1800, height=1050)
corr.phas <- QC.meas(FinalPlotDat, "Corrected.O2.phases")
corr.phas$main <- fishname
print(corr.phas)
dev.off()

#after MLND calculated using chabot function, filter SMR.clean for all used phases so can be plotted etc. DO NOT ANALYSE THIS DATA FOR 
#NEW MR CALCS we have the data in the chabot function. 
write.csv(FinalSlopes, paste0("FinalSlopesMLND_",unique(chamber$pumpresp),"_", unique(chamber$Batch),".csv"))


  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

#save SMR data
 FullSMR <- rbindlist(SMRfull)
write.csv(FullSMR, "Laukaa_Summer_2023_SMR_data_FULL.csv")

FullSMRsum <- rbindlist(SMRfullsum)
write.csv(FullSMRsum, "Laukaa_Summer_2023_SMR_data_SUMMARY.csv")

dev.off()



