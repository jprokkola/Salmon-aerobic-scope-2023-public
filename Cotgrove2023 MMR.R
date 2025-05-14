#Script to join Respirometry data and PumpResp times 
#Experiment in Laukaa LUKE facility Summer 2023 using pumpresp set up in 16 chambers (Shrek2) 
#Lucy Cotgrove, Jenni Prokkola, Evan Sala

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
library(purrr)
library(dplyr)
library(rlist)
library(lubridate)

#working directory
setwd("")
path <- getwd()

#post and pre background respiration
meanpostBG <- read.delim("meanPostBG.txt", stringsAsFactors=TRUE)
meanpreBG <- read.delim("meanPreBG.txt", stringsAsFactors=TRUE)

# Load chamber data
chamberdat <- read.csv("combinedchamberdata.csv")
chamberdat$Chamber.No <- paste0("CH", chamberdat$chamber)
#
#

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



mmrfiles <- list.files(patt="MMR_final_")
mmrfiles <- mmrfiles[-grep(".pdf", mmrfiles, fixed=T)]
mmrfinal <- list()
mmr.clean <- list()
MRdat <- list()


rm(j)
for (j in 1:length(mmrfiles)){
  tryCatch({  
    
    
    mmrfinal[[j]] <- read.delim(mmrfiles[[j]], stringsAsFactors=TRUE)
    
    
    ##################################
    # recalculate mr from this new df #
    #read and filter BGdata
    chamber <- subset(chamberdat, chamberdat$Batch == unique(mmrfinal[[j]]$Batch) & (chamberdat$pumpresp == unique(mmrfinal[[j]]$pumpresp)) & (chamberdat$Chamber.No == unique(mmrfinal[[j]]$Chamber_ID)))
    
    #make dummy NA data for missing info
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
    chamber <- plyr::rbind.fill(chamber, missingchamber)
    chamber <- chamber %>% distinct(Chamber_ID, .keep_all = T) %>% arrange(Chamber_ID)
    missing <- chamber[is.na(chamber$Volume_g), ]$Chamber.No
    
    #fishresp input info 
    info <- input.info(ID = c(chamber$Fish_ID),
                       Mass = c(chamber$Weight_g),
                       Volume = c(chamber$Volume_g),
                       DO.unit = "mg/L")
    
    
    
    # Import pre-test data
    pre <- subset(meanpreBG, Batch == as.character(unique(mmrfinal[[j]]$Batch)))
    # Import post-test data
    post <- subset(meanpostBG, Batch == as.character(unique(mmrfinal[[j]]$Batch)))
    fishname <- paste0(unique(chamber$pumpresp), "_MMR_final_", unique(mmrfinal[[j]]$Chamber_ID), "_", unique(chamber$Batch), ".txt")
    info <- info %>% drop_na()
    
    
    #this was because of a 10sec misalignment with times (see meas.to.wait)
    endlist <- list("PR1_MMR_final_CH2_B2.txt, PR1_MMR_final_CH3_B2.txt, PR3_MMR_final_CH1_B4.txt,
       PR3_MMR_final_CH2_B1.txt, PR3_MMR_final_CH3_B1.txt, PR3_MMR_final_CH3_B2.txt, PR3_MMR_final_CH4_B4.txt")
    
    
    mmr.raw <- if(fishname %in% endlist){
      # Import all actual SMR data. We only have one chamber per file here.
      import.meas(file = fishname,
                  info.data = info,
                  logger = "FishResp",
                  n.chamber = 1,
                  meas.to.wait = 30,
                  meas.to.flush = 60,
                  date.format = "YMD",
                  start.measure = "11:00:00",
                  stop.measure = "15:00:00",
                  plot.temperature = F,
                  plot.oxygen = F)
    } else {
      # Import all actual MR data.
      import.meas(file = fishname,
                  info.data = info,
                  logger = "FishResp",
                  n.chamber = 1,
                  meas.to.wait = 20,
                  meas.to.flush = 60,
                  date.format = "YMD",
                  start.measure = "11:00:00",
                  stop.measure = "15:00:00",
                  plot.temperature = F,
                  plot.oxygen = F)
    }
    
    
    #############################################################
    #-----------------------------------------------------------#
    #------------ STEP 3: BACK RESP CORRECTION -----------------#
    #-----------------------------------------------------------#
    #############################################################
    mmr.clean[[j]] <- correct.meas(info.data = info,
                                   pre.data =  pre,
                                   post.data = post,
                                   meas.data = mmr.raw,
                                   method = "linear")
    mmr.clean[[j]]$Chamber.No <- unique(mmrfinal[[j]]$Chamber_ID)
    
    
    #visualise  
    #RENAME HERE WITH YOUR FILE PATH
    png(file=paste0("/qc_",unique(chamber$pumpresp),unique(mmrfinal[[j]]$Chamber_ID) , unique(chamber$Batch),"_corrphas.png"),
        width=1800, height=1050)
    corr.phas <- QC.meas(mmr.clean[[j]], "Corrected.O2.phases")
    corr.phas$main <- fishname
    print(corr.phas)
    dev.off()
    
    
    }
    
    # Function for plot:
    rawdat.multiplot <- function(data) {
      lapply(data, function(x){
        plotdat = data
        plotdat <-  plotdat[!is.na(plotdat$O2.correct), ]
        plotdat$O2.correct <- is.finite(plotdat$O2.correct)
        DD = 1:length ( plotdat$O2.correct )
        MM = plotdat$O2.correct 
        spline <- smooth.spline(plotdat$O2.correct  ~ 1:length ( plotdat$O2.correct)  , df = 10  )
        plotdat$spl <- spline$y
        p <- ggplot(data = plotdat, 
                    aes(x = Time, y = O2.correct ))+
          geom_point(cex = 0.4)+
          geom_line(aes(x = Time, y = spl ), color = "green")+ 
          ylab("O2 mg")+
          theme_bw()+
          theme(axis.text.x = element_text(angle = 90, size = 8),
                axis.title.x = element_blank())
        
        print(p)
        rm(plotdat)
        rm(i)
      })}
    
      #  
    if(unique(mmrfinal[[j]]$Batch) %in% batchnumbers$batch & unique(mmrfinal[[j]]$pumpresp) %in% batchnumbers$pumpresp){
      
      # Plot and check plots for slope quality (if too much noise, exclude. Some can have pump malfunctions or bubbles)
      pdf(file = paste0(unique(chamber$pumpresp), "_MMR_final_", unique(mmrfinal[[j]]$Chamber_ID), "_", unique(chamber$Batch), ".pdf"), width = 5, height = 5)
      rawdat.multiplot(mmr.clean[[j]])
      dev.off()
      
      
      ### Spline method to get tangents at 1 and 20 seconds (T1 and T2). This also plots all slopes.
      mmr.proc = t(sapply ( 1: length(mmr.clean[[j]]) , function(i) {
        mmr.clean[[j]] <-  mmr.clean[[j]][!is.na(mmr.clean[[j]]$O2.correct), ]
        mmr.clean[[j]]$O2.correct <- is.finite(mmr.clean[[j]]$O2.correct)
        DD = 1:length ( mmr.clean[[j]]$O2.correct )
        MM = mmr.clean[[j]]$O2.correct
        spl <- smooth.spline(MM ~ DD , df = 10)
        plot ( MM~DD, ylab= "cool", cex = 0.4, xaxp = c(0, 900, 10))
        lines(spl, col="red" , lwd=3)
        T1 = predict(spl, x=1, deriv=1)$y #Compare tangents here T1
        T2 = predict(spl, x=20, deriv=1)$y #and here T2
        c(T1,T2)
      }))
      
      plot(mmr.proc[,1]) ## slope b. all look reasonable. 
      plot(mmr.proc[,2]) 
      
      
      # Compare slopes from T1 and T2.
      plot(mmr.proc[,1],mmr.proc[,2])
      abline(0,1) 
      #Almost identical.
      colnames(mmr.proc)[1:2] = c("T1slp","T20slp")
      
      
      ## Exclude poor data identified earlier.
      
      # Make a new data frame with slopes and original cols (this only takes the first row for each fish,
      # i.e. the BR data will not be relevant). Slopes from T1 and T2 in 1 and 2.
      MMR_data = cbind (as.data.frame(mmr.clean[[j]][1,]) , unique(mmr.proc)) 
      
      
      # Add index to each fish
      MMR_data$index = unique(as.numeric(gsub("\\D", "", mmr.clean[[j]]$Chamber.No)))
      
      ## Continue with FishResp for MMR calc
      
      # Change format into Fishresp slope format:
      # Cols that need to be in the data (empty or not) "Chamber.No"    "Ind"           "Mass"          "Volume"        "Date.Time"     "Phase"  "Temp"          "Slope.with.BR" "Slope"         "SE"            "R2"            "DO.unit"     
      
      
      MMR_T1 <- MMR_data %>%
        dplyr::select("Chamber.No",   "Ind","Mass","Volume",
                      "Date.Time","Phase", "Temp","T1slp", 
                      "index",  "DO.unit") %>%
        mutate(SE = NA,
               R2 = NA,
               BR = NA,
               Slope.with.BR =NA) %>%
        dplyr::rename("Slope" = T1slp) %>%
        mutate(Ind = (as.factor(Ind)))
      
      # Calculate MMR from slope
      MRdat[[j]] <- calculate.MR(MMR_T1, 
                                 density=1000, 
                                 plot.BR = F,
                                 plot.MR.abs = F,
                                 plot.MR.mass = F)
      
      
      #Remove NA cols
      MRdat[[j]] <- MRdat[[j]] %>%
        dplyr::select(-c(MR.abs.with.BR, Slope.with.BR, BR, SE, Ind)) %>%
        droplevels() %>%
        mutate(pumpresp = unique(chamber$pumpresp),
               Batch = unique(chamber$Batch),
               Chamber_ID = unique(mmrfinal[[j]]$Chamber_ID))
    }else{NULL}
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), print(j), "\n")})
}

MMRdata <- rbindlist(MRdat, fill = T)
### Combine with sampling and SMR data in another script

write.table(MMRdata, "Laukaa_2023_MMR.txt", col.names =T, row.names = F, sep = '\t')

#######################################
