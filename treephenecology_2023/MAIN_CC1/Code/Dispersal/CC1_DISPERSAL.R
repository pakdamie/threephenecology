#############################
#This is the analyses for   #
#CC1 dispersal groups       #
#############################

###All new code
#######################
###Packages required###
#######################
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



#########################
###THE DATA #############
#########################
###This is the cleaned up data
CC1_Agg <- read.csv(here('CC1','Data',"CC1_ALL_FINAL.csv"))

##################################
##################################

species_name_CC1 <- colnames(CC1_Agg)

###Here I'm cleaning up a lot of the species name
###Getting rid of dashes, periods, and spaces and such
species_name_CC1<- gsub("_", "",
                        species_name_CC1)
species_name_CC1<- sub(" ", "",species_name_CC1)
species_name_CC1 <- gsub("[.]","",species_name_CC1 )
species_name_CC1 <- gsub("\\(","",species_name_CC1)
species_name_CC1<- gsub("\\)","",species_name_CC1)
species_name_CC1 <- gsub("\"","",species_name_CC1)

species_name_CC1 <- data.frame(name=(species_name_CC1),
                               found =0)
species_name_CC1$name <- as.character(species_name_CC1$name)

colnames(CC1_Agg)[2:655] <-  species_name_CC1$name[2:655]

###Cleaned up the columnn names

###############################
###############################
###############################

DISPERSAL1 <- read.csv(here("CC1","Data","MDDCC1clust_.csv")) #Dispersal group 1 information
DISPERSAL2 <- read.csv(here("CC1","Data","CC1_dispsynMiss_23Apr2020_varun.csv")) #Dispersal group 2 information

###CLEANING UP DISPERSAL 1 characters
DISPERSAL1$Spp <- tolower(DISPERSAL1$Spp)
DISPERSAL1$Spp<- gsub("_", "",
                          DISPERSAL1$Spp)
DISPERSAL1$Spp<- sub(" ", "",DISPERSAL1$Spp)
DISPERSAL1$Spp <- gsub("[.]","",DISPERSAL1$Spp )
DISPERSAL1$Spp <- gsub("\\(","",DISPERSAL1$Spp)
DISPERSAL1$Spp<- gsub("\\)","",DISPERSAL1$Spp)
DISPERSAL1$Spp <- gsub("\"","",DISPERSAL1$Spp)
DISPERSAL1$Spp <- gsub("#","",DISPERSAL1$Spp)

###CLEANING UP DISPERSAL 2 characters
DISPERSAL2$Species <- tolower(DISPERSAL2$Species)
DISPERSAL2$Species<- gsub("_", "",
                        DISPERSAL2$Species)
DISPERSAL2$Species<- sub(" ", "",DISPERSAL2$Species)
DISPERSAL2$Species <- gsub("[.]","",DISPERSAL2$Species )
DISPERSAL2$Species <- gsub("\\(","",DISPERSAL2$Species)
DISPERSAL2$Species<- gsub("\\)","",DISPERSAL2$Species)
DISPERSAL2$Species <- gsub("\"","",DISPERSAL2$Species)
DISPERSAL2$Species <- gsub("#","",DISPERSAL2$Species)

CC1_SPECIES <-colnames(CC1_Agg[,c(2:655)])

#####################################################
#####################################################

############################
###FOR DISPERSAL GROUP 1###
###########################

#########################################################################
###THIS JUST SEES IF THE SPECIES IN THE DISPERSAL GROUP LIST IS FOUND####
###IN THE TIME SEREIS - IF FOUND, THEN TRUE.                         ####
##IF NOT FOUND, THEN FALSE                                           #### 
#########################################################################

CC1_SPECIES_T_F_D1 <- NULL
###There are 655 species in the time series and I"m checking if it's found
###in the first dispersal group 
for (s in seq(1,654)){
  CC1_SPECIES_T_F_D1 [[s]]= data.frame(
    Species =  CC1_SPECIES[[s]],
    Found= as.character(any(CC1_SPECIES[[s]]==
                        DISPERSAL1$Spp)))}

CC1_SPECIES_T_F_D1_2 <- do.call(rbind, CC1_SPECIES_T_F_D1)

###SPECIES FROM THE TIME SERIES THAT WERE FOUND IN THE DISPERSAL GROUP 1
CC1_SPECIES_DISPERSAL_1 <- subset(CC1_SPECIES_T_F_D1_2,
  CC1_SPECIES_T_F_D1_2$Found==TRUE)


###654 but only 113 are found.


#########################
#########################
###FOR DISPERSAL 2#######
#########################
#########################
#######################################################################
###THIS JUST SEES IF THE SPECIES are found in the time THEN TRUE.                          #
##IF NOT FOUND, THEN FALSE                                            #
#######################################################################
CC1_SPECIES_T_F_D2 <- NULL
for (s in seq(1, 654)) {
  CC1_SPECIES_T_F_D2[[s]] = data.frame(Species =
                                        as.character( CC1_SPECIES[[s]]), Found =
                                         as.character(any(CC1_SPECIES[[s]] == DISPERSAL2$Species)))
  
}
CC1_SPECIES_T_F_D2_2 <- do.call(rbind,CC1_SPECIES_T_F_D2)

###SPECIES FOUND
CC1_SPECIES_DISPERSAL_2 <- subset(CC1_SPECIES_T_F_D2_2,
                                  CC1_SPECIES_T_F_D2_2$Found==TRUE)

###654 but only 119 are found

#####################
###ASSIGNING GROUP###
#####################
CC1_SPECIES_DISPERSAL_1$Dispersal <- 0

for (n in seq(1, 113)) {
  temp =  DISPERSAL1[(which(CC1_SPECIES_DISPERSAL_1$Species[[n]] ==
                              DISPERSAL1$Spp)), ]
  dispersal = temp$group
  CC1_SPECIES_DISPERSAL_1$Dispersal[n] = dispersal
}
#############################################################
############DISPERSAL GROUP 2 ################################
#############################################################

CC1_SPECIES_DISPERSAL_2$Dispersal <- 0


for (n in seq(1,119)){
  temp=  DISPERSAL2[(which(CC1_SPECIES_DISPERSAL_2$Species[[n]]==
                                DISPERSAL2$Species)),]
      dispersal=colnames(temp[which((temp==1)==TRUE)])
      CC1_SPECIES_DISPERSAL_2$Dispersal[n]=dispersal
}

##################################################
############GROUPS##############################
#################################################

CC1_GROUPS_DIS_1 <- split(CC1_SPECIES_DISPERSAL_1, 
                    CC1_SPECIES_DISPERSAL_1$Dispersal)

###THE SAME FOR ALL DISPERSAL GROUPS
x = CC1_Agg$Day_Year_Y

CC1_GROUPS_DIS_1_G1 <- CC1_GROUPS_DIS_1[[1]]$Species; length(CC1_GROUPS_DIS_1_G1)
CC1_GROUPS_DIS_1_G2 <- CC1_GROUPS_DIS_1[[2]]$Species; length(CC1_GROUPS_DIS_1_G2)
CC1_GROUPS_DIS_1_G3 <- CC1_GROUPS_DIS_1[[3]]$Species; length(CC1_GROUPS_DIS_1_G3)
CC1_GROUPS_DIS_1_G5 <- CC1_GROUPS_DIS_1[[5]]$Species; length(CC1_GROUPS_DIS_1_G5)
CC1_GROUPS_DIS_1_G6 <- CC1_GROUPS_DIS_1[[6]]$Species; length(CC1_GROUPS_DIS_1_G6 )


###This is the function that is specific to this 
###script where we can calculate the actual WMR
###(not the permutated )- so this function basically
###ask what group you're interested in and if you want
###it logged or not (TRUE or False)
Cal_WMR_Dispersal_Group <- function(Group, Log){
  tmp <- CC1_Agg[,Group]
  
  tmp_X = CC1_Agg$Day_Year_Y
  
  if(Log == TRUE){
    tmp_Y  =log(tmp+1)
    
    tmp_WMR = wmr(mvcwt(tmp_X,  tmp_Y , min.scale =get.min.scale(tmp_X), 
                        max.scale = get.max.scale(tmp_X)))
    
  }
  if(Log == FALSE){
    tmp_Y  =tmp
    
    tmp_WMR = wmr(mvcwt(tmp_X,  tmp_Y , min.scale =get.min.scale(tmp_X), 
                        max.scale = get.max.scale(tmp_X)))
    
  }
  return(tmp_WMR)
}


G1_WMR_DF_NONLOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G1,FALSE)
G1_WMR_DF_LOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G1,TRUE)

G2_WMR_DF_NONLOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G2,FALSE)
G2_WMR_DF_LOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G2,TRUE)

G3_WMR_DF_NONLOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G3,FALSE)
G3_WMR_DF_LOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G3,TRUE)

G5_WMR_DF_NONLOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G5,FALSE)
G5_WMR_DF_LOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G5,TRUE)

G6_WMR_DF_NONLOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G6,FALSE)
G6_WMR_DF_LOG <- Cal_WMR_Dispersal_Group(CC1_GROUPS_DIS_1_G6,TRUE)



##############################################################
##############################################################
########################################################
###COMBINING THE DISPERSAL GROUP 1 and GROUP 2 together#
########################################################
##################################################################
###DISPERSAL GROUP 1 is 113 species and DISPERSAL GROUP 2 is 119 #
### 232 species to look at                                       #
##################################################################
DISPERSAL_TOTAL_SP <- c(as.character(CC1_SPECIES_DISPERSAL_1$Species),
 as.character( CC1_SPECIES_DISPERSAL_2$Species))

#####################################################
#####################################################
#########################PERMUTATING FUNCTION########
#####################################################
Perm25_NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,25,FALSE)
Perm25_LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,25,TRUE)

saveRDS(Perm25_NONLOG, file = 'Perm25NONLOG.RDS')
saveRDS(Perm25_LOG, file = 'Perm25_LOG.RDS')


Perm20NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,20,FALSE)
Perm20_LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,20,TRUE)

saveRDS(Perm20NONLOG, file = 'Perm20NONLOG.RDS')
saveRDS(Perm20_LOG, file = 'Perm20_LOG.RDS')


Perm52NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,52,FALSE)
Perm52_LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,52,TRUE)


saveRDS(Perm52NONLOG, file = 'Perm52NONLOG.RDS')
saveRDS(Perm52_LOG, file = 'Perm52_LOG.RDS')


Perm8NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,8,FALSE)
Perm8_LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,8,TRUE)


saveRDS(Perm8NONLOG, file = 'Perm8NONLOG.RDS')
saveRDS(Perm8_LOG, file = 'Perm8_LOG.RDS')

Perm5NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,5,FALSE)
Perm5_LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_SPECIES_DISPERSAL_1$Species,5,TRUE)


saveRDS(Perm5NONLOG, file = 'Perm5NONLOG.RDS')
saveRDS(Perm5_LOG, file = 'Perm5_LOG.RDS')

#####################################################################
#####################################################################
#######################GROUP 1 #####################################
#####################################################################

G1_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm25NONLOG, G1_WMR_DF_NONLOG )
G1_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm25_LOG, G1_WMR_DF_LOG)


G1_NONLOG_FULL_GG <- Grapher_Dispersal_Line(G1_NONLOG_FULL)+ggtitle("Small birds (N=25)")
G1_LOG_FULL_GG <- Grapher_Dispersal_Line(G1_LOG_FULL)+ggtitle("Small birds (N=25)")

#######################################################################
#######################################################################

#####################################################################
#######################GROUP 2 #####################################
#####################################################################

G2_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm20NONLOG, G2_WMR_DF_NONLOG )
G2_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm20_LOG, G2_WMR_DF_LOG)


G2_NONLOG_FULL_GG <- Grapher_Dispersal_Line(G2_NONLOG_FULL)+ggtitle("Small vertebrate (N = 20)")
G2_LOG_FULL_GG <- Grapher_Dispersal_Line(G2_LOG_FULL)+ggtitle("Small vertebrate (N = 20)")

#####################################################################
#######################GROUP 3 #####################################
#####################################################################
G3_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm52NONLOG, G3_WMR_DF_NONLOG )
G3_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm52_LOG, G3_WMR_DF_LOG)


G3_NONLOG_FULL_GG <- Grapher_Dispersal_Line(G3_NONLOG_FULL)+ggtitle("Large Vertebrate (N = 52)")
G3_LOG_FULL_GG <- Grapher_Dispersal_Line(G3_LOG_FULL)+ggtitle("Large Vertebrate (N = 52)")



#####################################################################
#######################GROUP 5 #####################################
#####################################################################
G5_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm8NONLOG, G5_WMR_DF_NONLOG )
G5_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm8_LOG, G5_WMR_DF_LOG) 


G5_NONLOG_FULL_GG <- Grapher_Dispersal_Line(G5_NONLOG_FULL)+ggtitle("Wind (N = 8)")
G5_LOG_FULL_GG <- Grapher_Dispersal_Line(G5_LOG_FULL)+ggtitle("Wind (N = 8)")




#####################################################################
#######################GROUP 6 #####################################
#####################################################################
G6_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,PERM5NONLOG, G6_WMR_DF_NONLOG )
G6_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm5_LOG, G6_WMR_DF_LOG) 


G6_NONLOG_FULL_GG <- Grapher_Dispersal_Line(G6_NONLOG_FULL)+ggtitle("Bats (N = 5)")
G6_LOG_FULL_GG <- Grapher_Dispersal_Line(G6_LOG_FULL)+ggtitle("Bats (N = 5)")

###NONLOG###
(G1_NONLOG_FULL_GG + G6_NONLOG_FULL_GG)/
  (G2_NONLOG_FULL_GG + G5_NONLOG_FULL_GG)/
  (G3_NONLOG_FULL_GG +G3_NONLOG_FULL_GG)


###LOG###
(G1_LOG_FULL_GG + G6_LOG_FULL_GG)/
  (G2_LOG_FULL_GG + G5_LOG_FULL_GG)/
  (G3_LOG_FULL_GG +G3_LOG_FULL_GG)


##########################################################################
##########################################################################
##########################COMBINING BOTH #################################
#############################################################################
############################################################################

#####################
###ZOOPHORIC ########
####################
ZOO <- subset(CC1_SPECIES_DISPERSAL_2, CC1_SPECIES_DISPERSAL_2$Dispersal=='Zoochoric')

CC1_GROUP1_ZOO <- subset(CC1_SPECIES_DISPERSAL_1, 
                         CC1_SPECIES_DISPERSAL_1$Dispersal 
                         %in% c(1,2,3,6))

###THis should have 174 species 
CC1_GROUP12_ZOO_SP <- c(as.character(CC1_GROUP1_ZOO$Species),
                        as.character(ZOO $Species))

CC1_Agg_ZOO_SP_TS <- CC1_Agg[,CC1_GROUP12_ZOO_SP ]

ZOO_Y_NONLOG = CC1_Agg_ZOO_SP_TS 
ZOO_Y_LOG = log(CC1_Agg_ZOO_SP_TS + 1 ) 

ZOO_NONLOG_WMR<- wmr(mvcwt(x,  ZOO_Y_NONLOG , min.scale =
        get.min.scale(x), max.scale = 
        get.max.scale(x)))

ZOO_LOG_WMR<- wmr(mvcwt(x,  ZOO_Y_LOG , min.scale =
            get.min.scale(x), max.scale = 
            get.max.scale(x)))
#####################
###WIND ########
####################

WIND<- subset(CC1_SPECIES_DISPERSAL_2, CC1_SPECIES_DISPERSAL_2$Dispersal=='Wind')


CC1_GROUP1_WIND <- subset(CC1_SPECIES_DISPERSAL_1, 
                          CC1_SPECIES_DISPERSAL_1$Dispersal 
                          %in% c(5))

CC1_GROUP12_WIND_SP <- c(as.character(CC1_GROUP1_WIND$Species),
                         as.character(WIND$Species))

CC1_Agg_WIND_SP_TS <- CC1_Agg[,CC1_GROUP12_WIND_SP ]


WIND_Y_NONLOG = CC1_Agg_WIND_SP_TS 
WIND_Y_LOG = log(CC1_Agg_WIND_SP_TS + 1 ) 

WIND_NONLOG_WMR<- wmr(mvcwt(x,  WIND_Y_NONLOG , min.scale =
                             get.min.scale(x), max.scale = 
                             get.max.scale(x)))

WIND_LOG_WMR<- wmr(mvcwt(x,  WIND_Y_LOG , min.scale =
                          get.min.scale(x), max.scale = 
                          get.max.scale(x)))
###TOTAL 



CC1_DISP_CODES <- c(CC1_GROUP12_ZOO_SP,CC1_GROUP12_WIND_SP)


Perm174NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_DISP_CODES,174,FALSE)
Perm174LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_DISP_CODES,174,TRUE)

saveRDS(Perm174NONLOG, file = "Perm174NONLOG.RDS")
saveRDS(Perm174LOG, file = "Perm174LOG.RDS")

Perm48NONLOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_DISP_CODES,48,FALSE)
Perm48LOG<- Permutation_Boot_Strap(data.table(CC1_Agg),CC1_DISP_CODES,48,TRUE)

saveRDS(Perm48NONLOG, file = "Perm48NONLOG.RDS")
saveRDS(Perm48LOG, file = "Perm48LOG.RDS")

###################################################

CC1_Agg_WIND_SP_TS <- CC1_Agg[,CC1_GROUP12_WIND_SP ]
CC1_Agg_WIND_SP_TS_Y <- CC1_Agg_WIND_SP_TS

CC1_WIND12_MVCWT = mvcwt(x,  CC1_Agg_WIND_SP_TS_Y , min.scale =
                          get.min.scale(x), max.scale = 
                          get.max.scale(x))
WIND12_WMR = wmr(CC1_WIND12_MVCWT)
WIND12_WMR_DF <- data.frame(Scale =WIND12_WMR$y, WMR = colMeans(WIND12_WMR$z),
                           Group="WIND12")

ZOO_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm174NONLOG, ZOO_NONLOG_WMR )
ZOO_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm174LOG,ZOO_LOG_WMR) 

WIND_NONLOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm48NONLOG, WIND_NONLOG_WMR )
WIND_LOG_FULL <- 
  Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI0,Perm48LOG,WIND_LOG_WMR) 



ZOO_NONLOG_FULL_GG <- Grapher_Dispersal_Line(ZOO_NONLOG_FULL)+ggtitle("Animal (N = 174)")
ZOO_LOG_FULL_GG <- Grapher_Dispersal_Line(ZOO_LOG_FULL)+ggtitle("Animal (N = 174)")


WIND_NONLOG_FULL_GG <- Grapher_Dispersal_Line(WIND_NONLOG_FULL)+ggtitle("Wind (N = 48)")
WIND_LOG_FULL_GG <- Grapher_Dispersal_Line(WIND_LOG_FULL)+ggtitle("Wind (N = 48)")

ZOO_NONLOG_FULL_GG/WIND_NONLOG_FULL_GG

ggsave("NON_LOG_CC1_ALL.pdf", width =9, height =9,
       units='in')

ZOO_LOG_FULL_GG/WIND_LOG_FULL_GG


ggsave("LOG_CC1_ALL.pdf", width =9, height =9,
       units='in')
