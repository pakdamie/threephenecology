
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
library(patchwork)

###This is a function needed
###

#########################
###THE DATA #############
#########################
###This is the cleaned up data
CC1_Agg <- read.csv(here('CC1','Data',"CC1_ALL_FINAL.csv"))


#######################
###The family codes# 
#######################
###This is the family code 
FAMILY_CC1<- read.csv(here("CC1","Data","SPECIES_NAME_2020.csv"))

###Fixing the strings- no _, #, empty space, period, dashes, ?
FAMILY_CC1$Species<- gsub("_", "",FAMILY_CC1$Species)
FAMILY_CC1$Species<- gsub("#", "",FAMILY_CC1$Species)
FAMILY_CC1$Species<- gsub(" ", "",FAMILY_CC1$Species)
FAMILY_CC1$Species <- gsub("\\.","",FAMILY_CC1$Species)
FAMILY_CC1$Species <- gsub("\\(","",FAMILY_CC1$Species)
FAMILY_CC1$Species <- gsub("\\)","",FAMILY_CC1$Species)
FAMILY_CC1$Species <- gsub("\"","",FAMILY_CC1$Species)
FAMILY_CC1$Species <- gsub("?","",FAMILY_CC1$Species,fixed = TRUE)

#lower case it
FAMILY_CC1$Species<- tolower(FAMILY_CC1$Species)

#don't make it into a factor
FAMILY_CC1$Family <- as.character(FAMILY_CC1$Family)


###Assigining family to the species 
species_name_CC1 <- colnames(CC1_Agg)

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

colnames(CC1_Agg)[2:656] <-  species_name_CC1$name[2:656]



for (n in seq(2,656)){
  species_name_CC1$found[[n]]= as.character(any(species_name_CC1$name[[n]]==FAMILY_CC1$Species))
}

###MANUALLY CHECK THROGH EXCEL

species_not_found <- subset(
  species_name_CC1,
  species_name_CC1$found == 'FALSE')

#hypocratheaceae is a family, but probably not going to use either
#because it has only one species ?

###So only thing to include:
#ingaspsmall1 is a Fabaceae (include that in manually)

species_name_CC1$family = 0

species_name_CC1_Found <- subset(species_name_CC1,
                                 species_name_CC1$found=='TRUE')


for (n in seq(1,633)){
  species_name_CC1_Found$family[n]=  FAMILY_CC1$Family[[(which(species_name_CC1_Found$name[[n]]==
                                                                 FAMILY_CC1$Species)[[1]])
  ]]
}



##Look for Urticaceae
BIG_Code <- subset(species_name_CC1_Found,
                   species_name_CC1_Found$family == "Bignoniaceae")$name

CC1_BIG_SP_TS <- CC1_Agg[,BIG_Code]

CC1_Sum_Seed <- data.frame(sum=rowSums(CC1_BIG_SP_TS))
CC1_BIG_SP_TS$Day_Year_Y <- CC1_Agg$Day_Year_Y
CC1_Sum_Seed$Day_Year_Y <-  CC1_Agg$Day_Year_Y
CC1_BIG_SP_TS_melt <- reshape2::melt(CC1_BIG_SP_TS, id.vars=c("Day_Year_Y"))

ggplot(CC1_BIG_SP_TS_melt, aes (x = Day_Year_Y, y = log(value+1), color=variable,group=variable))+
  geom_line()+scale_colour_viridis(option='turbo',discrete=TRUE)+
  xlab("Years")+ylab("Abundance (log + 1)")+ggtitle("Cocha Cashu - Bignoniaceae") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y=element_text(size = 14,color='black'),
        axis.title.x =element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



###This is for the total seed count
 ggplot(CC1_Sum_Seed, aes(x = Day_Year_Y, y= log(sum+1)))+
  geom_line()+
  scale_x_continuous(breaks = c(2.5,5,7.5),
                     expand = c(0, 0),
                     labels= c(2004,2007,2009))+
  xlab("Year")+
  ylab("Total seed count (log10)" )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=20))
