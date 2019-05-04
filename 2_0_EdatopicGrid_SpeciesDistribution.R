
######################################################################
#####CREATE EDATOPIC GRIDS WITH VEG STATS IN EACH CELL################
###############################################################
##Kiri Daust, July 2018
##MacKenzie, August 2018 extensive updates
# Builds edatopic grids by subzone variant showing species abundance in three categories
#by edatopic position from BECMaster Plot data records

.libPaths("E:/R packages351")
#install.packages("Hmisc")
require(reshape)
require(reshape2)
require(vegan)
require(caret)
require(tcltk)
require(randomForest)
require(Matrix)
require(labdsv)
require(gdata)
require(MASS)
require(openxlsx)
require (C50)
require(tidyr)
require(stringr)
require(rpart)
require(tree)
require(rattle)
require(rpart.plot)
require(partykit)
require(vegclust)
require(standardize)
require(dplyr)
require(tictoc)
require(plyr)
require(Hmisc)
require(ggplot2)
require(ggdendro)
require(pvclust)
require(dendextend)
require(ape)
rm(list=ls())
wd=tk_choose.dir(); setwd(wd)

###Function used above to combine data into formatted string
combineSpp <- function(x){ x <- x[order(-x$Order),]
  if(any(x$Pres == 1)){
    dom <- paste(x$Species[x$Pres == 1], collapse = ",")
    dom <- paste("*",dom,"*", sep = "") }else{ dom <- ""  }
  if(any(x$Pres == 2)){sub <- x$Species[x$Pres == 2]
    if(length(sub) > 5){
      sec <- paste("(", paste(sub[1:4], collapse = ","),"\n",paste(sub[5:length(sub)], collapse = ","),")",sep = "")
    }else{sec <- paste(sub, collapse = ",")    }
  }else{    sec <- ""  }
  if(any(x$Pres == 3)){sub <- x$Species[x$Pres == 3]
    if(length(sub) > 5){un <- paste("(", paste(sub[1:4], collapse = ","),"\n",paste(sub[5:length(sub)], collapse = ","),")",sep = "")
    }else{      un <- paste(sub, collapse = ",")
      un <- paste("(",un,")", sep = "")    }
  }else{    un <- ""  }
    return(paste(dom,"\n", sec, "\n", un, sep = ""))}

output.folder = ("../../../Output") # location of model outputs 

####Import Raw vegetation data
vegAll <- read.table("BECMasterVeg_Oct26_2018.txt", header = TRUE) # R export from Vpro
codeCross <- read.csv("CodeCrosswalk.csv", stringsAsFactors = FALSE) #tree species codes

vegAll <- separate(vegAll, Species, c("Species","Type"), "-", remove = TRUE)
vegAll <- vegAll[vegAll$Type %in% c(1,2),]

BGCLookup <- read.csv("BGCLookup.csv", stringsAsFactors = FALSE)###GIS plot location BGC
BGCLookup <- BGCLookup[,c(3,12)]
BGCLookup <- BGCLookup[BGCLookup$BGC_LABEL != "",]
colnames(BGCLookup)[1] <- "PlotNumber"

##import edatopic data
plotEnv <- read.csv("KiriEnvDat.csv", stringsAsFactors = FALSE)# plot, smr, snr, BGC
plotEnv <- plotEnv[plotEnv$NutrientRegime %in% c("A","B","C","D","E"),]
plotEnv <- plotEnv[plotEnv$MoistureRegime %in% c(0,1,2,3,4,5,6,7,8),]
plotEnv <- plotEnv[,-2]
plotEnv <- merge(plotEnv,BGCLookup, by = "PlotNumber", all.x = TRUE)
plotEnv$BGC_LABEL <- gsub("[[:space:]]","",plotEnv$BGC_LABEL)
plotEnv <- plotEnv[plotEnv$BGC_LABEL != "",]
colnames(plotEnv)[4] <- "Unit"

modBGC <- read.csv("ModelledBGC_Forested_err_removed.csv", stringsAsFactors = FALSE)
##problems with missing info. ESSFmmw, MHun, MHunp,MSun, SWBvk, SWBvks

for(i in 1:length(modBGC$BGC)){
  Unit <- modBGC$BGC[i]
  envSub <- plotEnv[plotEnv$Unit == Unit,]
  vegData <- merge(vegAll, envSub, by = "PlotNumber", all.y = TRUE)
  
  if(length(vegData$PlotNumber) > 2){
    vegData$Group <- paste(vegData$MoistureRegime,"-",vegData$NutrientRegime, sep = "")
    vegData <- vegData[!is.na(vegData$Species),]
    vegData <- vegData[,-c(3,5:7)]
    constCut <- 0.2
    ###roll up
    temp <- foreach(SS = unique(vegData$Group), .combine = rbind, .packages = "foreach") %dopar% {
      sub <- vegData[vegData$Group == SS,]
      num <- length(unique(sub$PlotNumber))
      foreach(Spp = unique(sub$Species), .combine = rbind) %do% {
        sub2 <- sub[sub$Species == Spp,]
        numSpp <- dim(unique(sub2[,1:2]))[1]
        mean <- mean(sub2$Cover)
        const <- numSpp/num
        if(const >= constCut){
          out <- data.frame(Group = SS, Species = Spp, MeanCov = mean, Constancy = const*100, Num = num)
        }
        
      }
    }
    
    vegGrid <- temp
    ###classify as 1,2 or 3
    vegGrid$Pres <- ifelse(vegGrid$MeanCov > 20 & vegGrid$Constancy > 50, 1,
                           ifelse(vegGrid$MeanCov > 10 & vegGrid$Constancy > 25,2,3))
    vegGrid$Order <- vegGrid$MeanCov*vegGrid$Constancy
    vegGrid <- merge(vegGrid,codeCross, by.x = "Species", by.y = "Code", all.x = TRUE)
    vegGrid <- vegGrid[,c(2,5:8)]
    colnames(vegGrid)[5] <- "Species"
    
    Lab <- ave(vegGrid[,c("Group","Pres","Species","Order")], vegGrid$Group, FUN = combineSpp)
    vegGrid <- separate(vegGrid, Group, c("Numeric","Alph"), "-", remove = TRUE)
    vegGrid$Lab <- Lab$Species
    vegGrid <- unique(vegGrid[,c(1:3,7)])
    
    ##plot
   pdf(file = paste("Output/", "EdaGrid_",Unit,".pdf",sep = ""), height = 10.5, paper = "letter")
    #pdf(file = paste("EdaGrid_",Unit,".pdf",sep = ""), height = 10.5, paper = "letter")
    
        print(ggplot(data = vegGrid)+
            geom_tile(aes(x= Alph, y = Numeric), color = "black", fill = "white")+
            geom_text(aes(x = Alph, y = Numeric, label = Lab), size = 3)+
            geom_text(aes(x = Alph, y = Numeric, label = Num,hjust = -4, vjust = -6), size = 2, color = "red")+
            scale_y_discrete(limits = c("8","7","6","5","4","3","2","1","0"))+
            scale_x_discrete(limits = c("A","B","C","D","E"))+
            labs(x = "Relative Soil Nutrient Regime", y = "Relative Soil Moisture Regime", cex = 1.5, title = paste(Unit," (",sum(vegGrid$Num)," plots)", sep = ""))+
            theme_bw(base_size = 10)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            coord_fixed())
    dev.off()
    
  }
}
