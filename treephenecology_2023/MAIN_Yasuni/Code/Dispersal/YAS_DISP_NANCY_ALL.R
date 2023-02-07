##################################################################
#### CODE THAT INCORPORATES THE NEW DISPERSAL GROUP FOR NANCY ###
##################################################################

#note: should be COI modified!

#####################
###package required##
#####################

library(here)
library(reshape2)
library(mvcwt) #not on CRAN, get the updated version
library(data.table)
library(dplyr)
library(biwavelet)
library(stringr)
library(ggplot2)
library(patchwork)

###########################
###########################

source("Permutation_Function_N.R")
MVCWT_YAS_WMR_COI0 <- readRDS(here("YASUNI","Output","MVCWT_YAS_WMR_COI0.RDS"))

#Most of the code is the making of the data files-
#if you're just interested in the graphs-
#run this
#
#Actual:
#load(here("YASUNI","Output","DISP_ALL","WMR_WIND_ACTUAL_NONLOG.RData"))
#load(here("YASUNI","Output","DISP_ALL","WMR_WIND_ACTUAL_LOG.RData"))
#load(here("YASUNI","Output","DISP_ALL","WMR_ANIM_ACTUAL_NONLOG.RData"))
#load(here("YASUNI","Output","DISP_ALL","WMR_ANIM_ACTUAL_LOG.RData"))

#saveRDS(WMR_WIND_ACTUAL_NONLOG,file = "WMR_WIND_ACTUAL_NONLOG.RDS", compress=FALSE)
#saveRDS(WMR_WIND_ACTUAL_LOG,file = "WMR_WIND_ACTUAL_LOG.RDS", compress=FALSE)
#saveRDS(WMR_ANIM_ACTUAL_NONLOG,file = "WMR_ANIM_ACTUAL_NONLOG.RDS", compress=FALSE)
#saveRDS(WMR_ANIM_ACTUAL_LOG,file = "WMR_ANIM_ACTUAL_LOG.RDS", compress=FALSE)

#Perm:
YAS_ANIMAL_ALL_LOG<- 
  readRDS(here("YASUNI","Output","DISP_ALL", "PERM_LOG",
            "YAS_ANIMAL_ALL_LOG.rds"))

YAS_WIND_ALL_LOG<- 
  readRDS(here("YASUNI","Output","DISP_ALL", "PERM_LOG",
               "YAS_WIND_ALL_LOG.rds"))

YAS_ANIMAL_ALL_NONLOG<- 
  readRDS(here("YASUNI","Output","DISP_ALL", "Perm_NONLOG",
               "YAS_ANIMAL_ALL_NONLOG.rds"))

YAS_WIND_ALL_NONLOG<- 
  readRDS(here("YASUNI","Output","DISP_ALL", "Perm_NONLOG",
               "YAS_WIND_ALL_NONLOG.rds"))
#########################################################
#########################################################
#########################################################

###The file needed [-1] because the rows are included
YASUNI_F <- fread(here("YASUNI","Data","YASUNI_ALL_FINAL.csv"))[,-1]


########################################################################
#############THE DISPERSAL CODE ########################################
########################################################################

###THERE ARE 1067 rows
YAS_DISPERSAL_CODE  <- read.csv(here("YASUNI",
                                     "Data",
                                     "YasuniAnimalWindDispersalSpecies.csv"))
nrow(YAS_DISPERSAL_CODE)

############
###ANIMAL###
############
#901 species 
YAS_ANIMAL <- subset(YAS_DISPERSAL_CODE,
                     YAS_DISPERSAL_CODE$Dispersal == 'Animal')


###Only 650 

ACTUAL_PERM_ANIMAL <- subset(YASUNI_F, YASUNI_F$variable%in% 
                    YAS_ANIMAL$CODIGO)

#########################################################################
########################################################################

###This is analysis just for the animal 

YAS_NANCY_ANIM_WMR_ACT <- function(log){
  
  ANIMAL_WIDE_DAT<-reshape2::dcast(ACTUAL_PERM_ANIMAL,
                 Day_Year_Y~variable,value.var='value',
                      direction="wide")

x_A = ANIMAL_WIDE_DAT$Day_Year_Y
#nonlogged
if(log == 'FALSE'){
  y_A =  ANIMAL_WIDE_DAT[,2:ncol( 
    ANIMAL_WIDE_DAT)]
}

#logged
if(log == 'TRUE'){

y_A =  log(ANIMAL_WIDE_DAT[,2:ncol( 
              ANIMAL_WIDE_DAT)]+1)
}
###This does the wavelet transformation for each column/species
MVWCT_FAM_ANIM = mvcwt(x_A, y_A, min.scale = get.min.scale(x_A), 
                       max.scale = get.max.scale(x_A))
#####################################################################
###After that is accomplished, we now do the wavelet modulus ratio##
####################################################################
WMR_ANIM_ACTUAL = wmr(MVWCT_FAM_ANIM)

return(WMR_ANIM_ACTUAL)
}

WMR_ANIM_ACTUAL_LOG <- YAS_NANCY_ANIM_WMR_ACT (log="TRUE")
save(WMR_ANIM_ACTUAL_LOG , file=here("YASUNI","Output","DISP_ALL",
                                     "WMR_ANIM_ACTUAL_LOG.RData"))

WMR_ANIM_ACTUAL_NONLOG <- YAS_NANCY_ANIM_WMR_ACT (log="FALSE")
save(WMR_ANIM_ACTUAL_NONLOG , file=here("YASUNI","Output","DISP_ALL",
                                     "WMR_ANIM_ACTUAL_NONLOG.RData"))

############
###WIND###
############

###There are 166 species
YAS_WIND <- subset(YAS_DISPERSAL_CODE, YAS_DISPERSAL_CODE$Dispersal == 'Wind')
#nrow(YAS_WIND)


###131 accounted for 
ACTUAL_PERM_WIND <- subset(
  YASUNI_F, YASUNI_F$variable %in% YAS_WIND$CODIGO)

WIND_WIDE_DAT<- reshape(ACTUAL_PERM_WIND , idvar = c("Day_Year_Y"), 
                          timevar = "variable", direction = "wide")
######
YAS_NANCY_WIND_WMR_ACT <- function(log){
  
  WIND_WIDE_DAT<-reshape2::dcast(ACTUAL_PERM_WIND,
                                   Day_Year_Y~variable,
                                  value.var='value',
                                   direction="wide")
  
  x_W = WIND_WIDE_DAT$Day_Year_Y
  
  if(log == 'FALSE'){
    y_W =  WIND_WIDE_DAT[,2:ncol( 
      WIND_WIDE_DAT)]
  }
  
  if(log == 'TRUE'){
    
    y_W =  log(WIND_WIDE_DAT[,2:ncol( 
      WIND_WIDE_DAT)]+1)
  }
  ###This does the wavelet transformation for each column/species
  MVWCT_FAM_WIND = mvcwt(x_W, y_W, min.scale = get.min.scale(x_W), 
                         max.scale = get.max.scale(x_W))
  ###After that is accomplished, we now do the wavelet modulus ratio
  ###
  WMR_WIND_ACTUAL = wmr(MVWCT_FAM_WIND)
  
  return(WMR_WIND_ACTUAL)
}
WMR_WIND_ACTUAL_LOG <- YAS_NANCY_WIND_WMR_ACT (log="TRUE")


save(WMR_WIND_ACTUAL_LOG , file=here("YASUNI","Output","DISP_ALL",
                                     "WMR_WIND_ACTUAL_LOG.RData"))


WMR_WIND_ACTUAL_NONLOG <- YAS_NANCY_WIND_WMR_ACT (log="FALSE")
save(WMR_WIND_ACTUAL_NONLOG , file=here("YASUNI","Output","DISP_ALL",
                                     "WMR_WIND_ACTUAL_NONLOG.RData"))

####################################################################
###################################################################
YAS_ALL_DISPERSALCODE <- unique(rbind(ACTUAL_PERM_ANIMAL,ACTUAL_PERM_WIND)$variable)

save(YAS_ALL_DISPERSALCODE, file = "YAS_ALL_DISPERSALCODE.RData")

#########################################################
################This makes the data-frame################
#########################################################


YAS_ANIMAL_ALL_LOG<-
  c(YAS_ANIMAL_ALL_LOG_A,
  YAS_ANIMAL_ALL_LOG_B,
  YAS_ANIMAL_ALL_LOG_C,
  YAS_ANIMAL_ALL_LOG_D)


YAS_ANIMAL_ALL_NONLOG<-
  c(YAS_ANIMAL_ALL_NONLOG_A,
    YAS_ANIMAL_ALL_NONLOG_B,
    YAS_ANIMAL_ALL_NONLOG_C,
    YAS_ANIMAL_ALL_NONLOG_D)


###ANIM
DF_WMR_ANIM_LOG <- Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,YAS_ANMAL_ALL_LOG ,WMR_ANIM_ACTUAL_LOG)
DF_WMR_ANIM_LOG$id = "animal"

DF_WMR_ANIM_NONLOG <- Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,YAS_ANIMAL_ALL_NONLOG ,WMR_ANIM_ACTUAL_NONLOG)
DF_WMR_ANIM_NONLOG$id = "animal"

###WIND
DF_WMR_WIND_LOG <- Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,YAS_WIND_ALL_LOG, WMR_WIND_ACTUAL_LOG)
DF_WMR_WIND_LOG$id = "wind"

DF_WMR_WIND_NONLOG <- Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,YAS_WIND_ALL_NONLOG, WMR_WIND_ACTUAL_NONLOG)
DF_WMR_WIND_NONLOG$id = "wind"

DF_ANIM_WIND_LOG <- rbind(DF_WMR_ANIM_LOG,DF_WMR_WIND_LOG)
DF_ANIM_WIND_NONLOG <- rbind(DF_WMR_ANIM_NONLOG,DF_WMR_WIND_NONLOG)

#########################################################
################This makes the figures################
#########################################################


####LOGGED
SUB_ANIM_WIND_LOG <- subset(DF_ANIM_WIND_LOG ,DF_ANIM_WIND_LOG$time=='sub' )
                           
SUB_ANIM_WIND_LOG$Cut <- cut(SUB_ANIM_WIND_LOG$scale , breaks = c(seq(0.1,1.2,0.02)),
                         dig.lab =4)
SUB_ANIM_WIND_LOG$Cut <- gsub("[()]", "",SUB_ANIM_WIND_LOG$Cut)
SUB_ANIM_WIND_LOG $Cut <- gsub("\\[|\\]", "", SUB_ANIM_WIND_LOG$Cut)
SUB_ANIM_WIND_LOG$Cut <- factor(gsub("[,]","-",SUB_ANIM_WIND_LOG$Cut))

################################

SUB_ANIM_WIND_LOG_SPLIT <- split(SUB_ANIM_WIND_LOG, 
                            list(SUB_ANIM_WIND_LOG$Cut,
                                 SUB_ANIM_WIND_LOG$id))

###This gets rid of empty entries
SUB_ANIM_WIND_LOG_SPLIT <-SUB_ANIM_WIND_LOG_SPLIT [sapply(SUB_ANIM_WIND_LOG_SPLIT, nrow) > 0]

for (k in seq(1,length(SUB_ANIM_WIND_LOG_SPLIT))){
  tmp <-  SUB_ANIM_WIND_LOG_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
           unique(tmp$synch_compen)[unique(tmp$synch_compen)!=("0")],
           tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }
  
  SUB_ANIM_WIND_LOG_SPLIT[[k]] = tmp
}

SUB_ANIM_WIND_LOG_DF <- do.call(rbind,  SUB_ANIM_WIND_LOG_SPLIT)


##############
###ANNNUAL####
##############

ANN_DF_ANIM_WIND_LOG <- subset(DF_ANIM_WIND_LOG, 
                       DF_ANIM_WIND_LOG$time=='ann' &
                         DF_ANIM_WIND_LOG$scale < 10)


ANN_DF_ANIM_WIND_LOG$Cut <- cut(ANN_DF_ANIM_WIND_LOG$scale , breaks = c(seq(1,8.2,0.2)))
ANN_DF_ANIM_WIND_LOG$Cut <- gsub("[()]", "",ANN_DF_ANIM_WIND_LOG$Cut)
ANN_DF_ANIM_WIND_LOG$Cut <- gsub("\\[|\\]", "", ANN_DF_ANIM_WIND_LOG$Cut)
ANN_DF_ANIM_WIND_LOG$Cut <- factor(gsub("[,]","-",ANN_DF_ANIM_WIND_LOG$Cut) )



ANN_DF_ANIM_WIND_LOG_SPLIT <- split(ANN_DF_ANIM_WIND_LOG, 
                            list(ANN_DF_ANIM_WIND_LOG$Cut,ANN_DF_ANIM_WIND_LOG$id))


for (k in seq(1,length(ANN_DF_ANIM_WIND_LOG_SPLIT))){
  tmp = ANN_DF_ANIM_WIND_LOG_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
           unique(tmp$synch_compen)[ unique(tmp$synch_compen)!=("0")],
           tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }
  
  ANN_DF_ANIM_WIND_LOG_SPLIT[[k]] = tmp
}

ANN_ANIM_WIND_LOG_DF <- do.call(rbind, ANN_DF_ANIM_WIND_LOG_SPLIT)

ANIM_LOG<-rbind.data.frame(
  subset(SUB_ANIM_WIND_LOG_DF ,
         SUB_ANIM_WIND_LOG_DF $id=='animal'),
  subset(ANN_ANIM_WIND_LOG_DF,
         ANN_ANIM_WIND_LOG_DF$id=='animal'))
WIND_LOG<-rbind.data.frame(
  subset(SUB_ANIM_WIND_LOG_DF ,
         SUB_ANIM_WIND_LOG_DF $id=='wind'),
  subset(ANN_ANIM_WIND_LOG_DF,
         ANN_ANIM_WIND_LOG_DF$id=='wind'))
#############################################

###
ANIM_LOG_LINE_GG <- Grapher_Dispersal_Line(ANIM_LOG)+ggtitle("Animal (N = 653)")
WIND_LOG_LINE_GG <- Grapher_Dispersal_Line(WIND_LOG)+ggtitle("Wind (N = 653)")

ANIM_LOG_LINE_GG/WIND_LOG_LINE_GG 



################
####NONLOGGED###
################
SUB_ANIM_WIND_NONLOG <- subset(DF_ANIM_WIND_NONLOG ,DF_ANIM_WIND_NONLOG$time=='sub')
SUB_ANIM_WIND_NONLOG$Cut <- cut(SUB_ANIM_WIND_NONLOG$scale , breaks = c(seq(0.1,1.2,0.02)),
                             dig.lab =4)
SUB_ANIM_WIND_NONLOG$Cut <- gsub("[()]", "",SUB_ANIM_WIND_NONLOG$Cut)
SUB_ANIM_WIND_NONLOG $Cut <- gsub("\\[|\\]", "", SUB_ANIM_WIND_NONLOG$Cut)
SUB_ANIM_WIND_NONLOG$Cut <- factor(gsub("[,]","-",SUB_ANIM_WIND_NONLOG$Cut))

################################

SUB_ANIM_WIND_NONLOG_SPLIT <- split(SUB_ANIM_WIND_NONLOG, 
                                 list(SUB_ANIM_WIND_NONLOG$Cut,
                                      SUB_ANIM_WIND_NONLOG$id))

###This gets rid of empty entries
SUB_ANIM_WIND_NONLOG_SPLIT <-SUB_ANIM_WIND_NONLOG_SPLIT [sapply(SUB_ANIM_WIND_NONLOG_SPLIT, nrow) > 0]

for (k in seq(1,length(SUB_ANIM_WIND_NONLOG_SPLIT))){
  tmp <-  SUB_ANIM_WIND_NONLOG_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
           unique(tmp$synch_compen)[unique(tmp$synch_compen)!=("0")],
           tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }
  
  SUB_ANIM_WIND_NONLOG_SPLIT[[k]] = tmp
}

SUB_ANIM_WIND_NONLOG_DF <- do.call(rbind,  SUB_ANIM_WIND_NONLOG_SPLIT)


##############
###ANNNUAL####
##############

ANN_DF_ANIM_WIND_NONLOG <- subset(DF_ANIM_WIND_NONLOG, 
                               DF_ANIM_WIND_NONLOG$time=='ann')


ANN_DF_ANIM_WIND_NONLOG$Cut <- cut(ANN_DF_ANIM_WIND_NONLOG$scale , breaks = c(seq(1,8.4,0.2)))
ANN_DF_ANIM_WIND_NONLOG$Cut <- gsub("[()]", "",ANN_DF_ANIM_WIND_NONLOG$Cut)
ANN_DF_ANIM_WIND_NONLOG$Cut <- gsub("\\[|\\]", "", ANN_DF_ANIM_WIND_NONLOG$Cut)
ANN_DF_ANIM_WIND_NONLOG$Cut <- factor(gsub("[,]","-",ANN_DF_ANIM_WIND_NONLOG$Cut) )



ANN_DF_ANIM_WIND_NONLOG_SPLIT <- split(ANN_DF_ANIM_WIND_NONLOG, 
                                    list(ANN_DF_ANIM_WIND_NONLOG$Cut,ANN_DF_ANIM_WIND_NONLOG$id))


for (k in seq(1,length(ANN_DF_ANIM_WIND_NONLOG_SPLIT))){
  tmp = ANN_DF_ANIM_WIND_NONLOG_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
           unique(tmp$synch_compen)[ unique(tmp$synch_compen)!=("0")],
           tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }
  
  ANN_DF_ANIM_WIND_NONLOG_SPLIT[[k]] = tmp
}

ANN_ANIM_WIND_NONLOG_DF <- do.call(rbind, ANN_DF_ANIM_WIND_NONLOG_SPLIT)

###

ANIM_NONLOG<-rbind.data.frame(
  subset(SUB_ANIM_WIND_NONLOG_DF ,
         SUB_ANIM_WIND_NONLOG_DF $id=='animal'),
  subset(ANN_ANIM_WIND_NONLOG_DF,
         ANN_ANIM_WIND_NONLOG_DF$id=='animal'))
WIND_NONLOG<-rbind.data.frame(
  subset(SUB_ANIM_WIND_NONLOG_DF ,
         SUB_ANIM_WIND_NONLOG_DF $id=='wind'),
  subset(ANN_ANIM_WIND_NONLOG_DF,
         ANN_ANIM_WIND_NONLOG_DF$id=='wind'))
#############################################

###
ANIM_NONLOG_LINE_GG <- Grapher_Dispersal_Line(ANIM_NONLOG)+ggtitle("Animal (N = 653)")
WIND_NONLOG_LINE_GG <- Grapher_Dispersal_Line(WIND_NONLOG)+ggtitle("Wind (N = 653)")

ANIM_NONLOG_LINE_GG/WIND_NONLOG_LINE_GG 

#Size: 8/9 


