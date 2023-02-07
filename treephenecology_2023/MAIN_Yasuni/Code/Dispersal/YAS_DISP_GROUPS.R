##################################################################
#### CODE THAT INCORPORATES THE NEW DISPERSAL GROUP FROM SIMON ###
##################################################################
#####################
###package required##
#####################

library(here) # To deal with the working directory
library(reshape2) # To deal with the melting long to wide/wide to long
library(mvcwt) # not on CRAN, get the updated version from the graveyard
library(data.table)
library(dplyr)
library(biwavelet)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tidyverse)

###
###
source("Permutation_Function_N.R")
MVCWT_YAS_WMR_COI0 <- readRDS(here("YASUNI","Output","MVCWT_YAS_WMR_COI0.RDS"))

#Most of the code is the making of the data files-
#if you're just interested in the graphs-
#run this
#
#########################################################
#########################################################
#########################################################

###The file needed [-1] because the rows names are included
###This is the cleaned-up data file and should be used for all analysis 
YASUNI_F <- fread(here("YASUNI","Data","YASUNI_ALL_FINAL.csv"))[,-1]

###Cast this into wide format- because this is in long format right now
YASUNI_WIDE_F_ALL<-reshape2::dcast(YASUNI_F,
                                   Day_Year_Y~variable,value.var='value',
                                   direction="wide")
############################
#Used for actual table######
#############################
YASUNI_DT <- data.table(YASUNI_WIDE_F_ALL)


#####################################
###NEW DISPERSAL CODE FOR YASUNI ###
####################################
################
###SIMON CODE###
################
# Read in SSP data from Nancy (I converted it to .CSV and removed all the rows with extra info at the top and the unidents at the bottom)
ssp <-read.csv(here("YASUNI","Data",'YASUNI_TRAP_CODE_SUMMARY_20180725_CSV.csv'))
# Read in Dispersal data from Simon
yas <- read.csv(here("YASUNI","Data",'YASUNI_DISPERSE_NEW.csv'))

## Create new column in 'yas' of the binomial
yas$Name <- paste(yas$genus, yas$species, sep = ' ')

# merge the two data sets by the latin binomial, keeping all rows of both
MERGED_CODIGO_DISPERSE <- merge(ssp, yas, by.x = 'NombreActual', by.y = 'Name', all = TRUE)

###############################################################
######################END HERE #################################
################################################################

##First get the unique species names and omit any that is Na
SPECIES_YAS <- na.omit(data.frame(CODIGO = unique(YASUNI_F$variable)))

###I'm getting the code and the Code and Dispersal Syndromes 
CD_YAS <-na.omit(data.frame(CODIGO = MERGED_CODIGO_DISPERSE $CODIGO, 
                     Syndrome  = MERGED_CODIGO_DISPERSE$Dsynd))
CD_YAS_F<- subset(CD_YAS, CD_YAS$Syndrome != '')

###Unnecessary row 
CD_YAS_F <- CD_YAS_F[-328,]

###Merging the species that are in the full time series as well as 
###Simon's code that include the dispersal syndrome 
UNIQUE_CODE <-na.omit(left_join(SPECIES_YAS, CD_YAS_F))


###Here are all 6 dispersal syndromes and all the specices 
###associted with each group
LESS_2CM <- subset(UNIQUE_CODE, UNIQUE_CODE$Syndrome == '<2cm')
WIND <- subset(UNIQUE_CODE, UNIQUE_CODE$Syndrome == 'wind')
BALLISTIC <- subset(UNIQUE_CODE,UNIQUE_CODE$Syndrome == 'ballistic')
TERRESTRIAL <- subset(UNIQUE_CODE, UNIQUE_CODE$Syndrome == 'terrestrial')
BET_2_5CM <- subset(UNIQUE_CODE,UNIQUE_CODE$Syndrome == "2 - 5 cm")
GREATER_5CM<- subset(UNIQUE_CODE,UNIQUE_CODE$Syndrome == '>5cm')


###################################################
###################################################
###################################################

###Calculating the actual WMR fr each dispersal group


###This is the function that is specific to this 
###script where we can calculate the actual WMR
###(not the permutated )- so this function basically
###ask what group you're interested in and if you want
###it logged or not (TRUE or False)
Cal_WMR_Dispersal_Group <- function(Group, Log){
  tmp <- subset(YASUNI_F,
                YASUNI_F$variable %in% 
                Group$CODIGO)
  
  tmp_WIDE<-dcast(tmp ,Day_Year_Y~variable,value.var='value',
                        direction="wide")
  
  tmp_X = tmp_WIDE$Day_Year_Y 
  
  if(Log == TRUE){
  tmp_Y  =log(tmp_WIDE[,2:ncol(tmp_WIDE)]+1)
  
  tmp_WMR = wmr(mvcwt(tmp_X,  tmp_Y , min.scale =get.min.scale(tmp_X), 
                   max.scale = get.max.scale(tmp_X)))
  
  }
  if(Log == FALSE){
  tmp_Y  =tmp_WIDE[,2:ncol(tmp_WIDE)]
    
  tmp_WMR = wmr(mvcwt(tmp_X,  tmp_Y , min.scale =get.min.scale(tmp_X), 
                        max.scale = get.max.scale(tmp_X)))
    
  }
  return(tmp_WMR)
}
###################
###Less than 2CM ##
###################
##211 SPECIES
LESS_2CM_FIND_Log<-Cal_WMR_Dispersal_Group(LESS_2CM,Log=TRUE)
LESS_2CM_FIND_NonLog<-Cal_WMR_Dispersal_Group(LESS_2CM,Log=FALSE)

#save(LESS_2CM_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "LESS_2CM_FIND_Log.RData"))
#save(LESS_2CM_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "LESS_2CM_FIND_NonLog.RData"))

#############################################################################
######
#WIND#
######
####THERE ARE 30 SPECIES 

WIND_FIND_Log<-Cal_WMR_Dispersal_Group(WIND,Log=TRUE)
WIND_FIND_NonLog<-Cal_WMR_Dispersal_Group(WIND,Log=FALSE)

#save(WIND_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "WIND_FIND_Log.RData"))
#save(WIND_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "WIND_FIND_NonLog.RData"))

##################################################
##################################################
###########
#BALLISTIC#
###########
########### 
####THERE ARE 16 SPECIES 

BALL_FIND_Log<-Cal_WMR_Dispersal_Group(BALLISTIC,Log=TRUE)
BALL_FIND_NonLog<-Cal_WMR_Dispersal_Group(BALLISTIC,Log=FALSE)

#save(BALL_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "BALL_FIND_Log.RData"))
#save(BALL_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "BALL_FIND_NonLog.RData"))

##################################################
##################################################
#TERRESTRIAL 
####THERE ARE 25 SPECIES 

TERR_FIND_Log<-Cal_WMR_Dispersal_Group(TERRESTRIAL,Log=TRUE)
TERR_FIND_NonLog<-Cal_WMR_Dispersal_Group(TERRESTRIAL,Log=FALSE)

#save(TERR_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "TERR_FIND_Log.RData"))
#save(TERR_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "TERR_FIND_NonLog.RData"))

################################################
#####################
#BETWEEN 2 and 5 cm ##
######################
####THERE ARE 63 SPECIES 

BET_2_5CM_FIND_Log<-Cal_WMR_Dispersal_Group(BET_2_5CM,Log=TRUE)
BET_2_5CM_FIND_NonLog<-Cal_WMR_Dispersal_Group( BET_2_5CM,Log=FALSE)

#save(TERR_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "BET_2_5CM_FIND_Log.RData"))
#save(TERR_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "BET_2_5CM_FIND_NonLog.RData"))

#################################################
##############
#GREATER 5CM##
##############
####THERE ARE 10 SPECIES 

GREATER_5CM_FIND_Log<-Cal_WMR_Dispersal_Group( GREATER_5CM,Log=TRUE)
GREATER_5CM_FIND_NonLog<-Cal_WMR_Dispersal_Group( GREATER_5CM,Log=FALSE)

#save(GREATER_5CM_FIND_Log, file=here("YASUNI","Output","DISP_GROUP", "GREATER_5CM_FIND_Log.RData"))
#save(GREATER_5CM_FIND_NonLog, file=here("YASUNI","Output","DISP_GROUP", "GREATER_5CM_FIND_NonLog.RData"))

##################END SECTION OF NON PERMUTATED ANALYSIS
#######################################################################################
############################################PERMUTATION#################################
#######################################################################################

#LABELLER INCLDUDES- ALL THE SPECIES
UNIQUE_LABEL <-  UNIQUE_CODE$CODIGO


length((UNIQUE_LABEL)) 

###Unique Numbers 
unique_n<- c(211,30,16,25,63,10)

###LOG PERM
for (k in seq(1,6)){
  
  tmp <- Permutation_Boot_Strap(
    actualdata=YASUNI_DT,
    UNIQUE_LABEL, 
    disp_n= unique_n[[k]], 
    log=TRUE)
  
  assign(paste("LOGPERM", unique_n[[k]], sep = ""), tmp)
  write.csv(get(paste("LOGPERM", unique_n[[k]], sep = "")),
            file =here("YASUNI","Output",
            paste("LOGPERM", unique_n[[k]],".csv", sep = "")))
            rm(list=paste("LOGPERM", unique_n[[k]], sep = ""))
  
}

####################################################
####################################################
###REMEMBER TO REMOVE AS YOU GO OR IT WILL CHUG#####
####################################################
##########################
###I manually loaded it###
##########################
#######################
###The NonTransformed #
#######################
LESS_2CM_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM211, LESS_2CM_FIND_NonLog)

WIND_NONLOG_FULL <-   
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM30, WIND_FIND_NonLog)

BALLISTIC_NONLOG_FULL <-  
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM16,BALL_FIND_NonLog)

TERRESTRIAL_NONLOG_FULL <-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM25, TERR_FIND_NonLog)

BET_2_5CM_NONLOG_FULL <-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM63,BET_2_5CM_FIND_NonLog)

GREATER_5CM_NONLOG_FULL<-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,NONPERM10, GREATER_5CM_FIND_NonLog)

#######################
###The Transformed #
#######################

LESS_2CM_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM211, LESS_2CM_FIND_Log)

WIND_LOG_FULL <-   
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM30, WIND_FIND_Log)

BALLISTIC_LOG_FULL <-  
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM16,  BALL_FIND_Log)

TERRESTRIAL_LOG_FULL <-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM25, TERR_FIND_Log)

BET_2_5CM_LOG_FULL <-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM63,BET_2_5CM_FIND_Log)

GREATER_5CM_LOG_FULL<-
  Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,LOGPERM10, GREATER_5CM_FIND_Log)


#########################
#########Figures#########
########################
##########################
###NON LOG TRANSFORMED ###
##########################
LESS_2CM_NONLOG_GG <- Grapher_Dispersal_Line(LESS_2CM_NONLOG_FULL)+ggtitle("Less than 2 cm (N=211)")
WIND_NONLOG_GG <- Grapher_Dispersal_Line(WIND_NONLOG_FULL)+ggtitle("Wind (N=30)")
BALLISTIC_NONLOG_GG <- Grapher_Dispersal_Line(BALLISTIC_NONLOG_FULL)+ggtitle("Ballistic (N=16)")
TERR_NONLOG_GG <- Grapher_Dispersal_Line(TERRESTRIAL_NONLOG_FULL)+ggtitle("Terrestrial (N=25)")
BET_2_5CM_NONLOG_GG <- Grapher_Dispersal_Line(BET_2_5CM_NONLOG_FULL)+ggtitle("Between 2-5 cm (N=63)")
GREATER_5CM_NONLOG_GG <- Grapher_Dispersal_Line(GREATER_5CM_NONLOG_FULL)+ggtitle("Greater than 5 cm(N=10)")

(LESS_2CM_NONLOG_GG +BET_2_5CM_NONLOG_GG )/
(GREATER_5CM_NONLOG_GG +BALLISTIC_NONLOG_GG)/
(WIND_NONLOG_GG + TERR_NONLOG_GG)




#####################
###LOG TRANSFORMED###
#####################
LESS_2CM_LOG_GG <- Grapher_Dispersal_Line(LESS_2CM_LOG_FULL)+ggtitle("Less than 2 cm (N=211)")
WIND_LOG_GG <- Grapher_Dispersal_Line(WIND_LOG_FULL)+ggtitle("Wind (N=30)")
BALLISTIC_LOG_GG <- Grapher_Dispersal_Line(BALLISTIC_LOG_FULL)+ggtitle("Ballistic (N=16)")
TERR_LOG_GG <- Grapher_Dispersal_Line(TERRESTRIAL_LOG_FULL)+ggtitle("Terrestrial (N=25)")
BET_2_5CM_LOG_GG <- Grapher_Dispersal_Line(BET_2_5CM_LOG_FULL)+ggtitle("Between 2-5 cm (N=63)")
GREATER_5CM_LOG_GG <- Grapher_Dispersal_Line(GREATER_5CM_LOG_FULL)+ggtitle("Greater than 5 cm(N=10)")


(LESS_2CM_LOG_GG +BET_2_5CM_LOG_GG )/
  (GREATER_5CM_LOG_GG +BALLISTIC_LOG_GG)/
  (WIND_LOG_GG + TERR_LOG_GG)

