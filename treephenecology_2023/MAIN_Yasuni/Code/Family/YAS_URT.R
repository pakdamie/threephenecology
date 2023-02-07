
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
library(viridis)
###This is a function needed
###

###The file needed [-1] because the rows are included
YASUNI_F <- fread(here("YASUNI","Data","YASUNI_ALL_FINAL.csv"))[,-1]
###USE THE FAMILY_CODE.csv- This is the file that has the unique species name #
### and the family it belongs to   
FAMILY_CODE <-  read.csv(here("YASUNI","Data",'CODE_YAS_Lasky_ncg.csv'))


###Cast this into wide format 
YASUNI_WIDE_F_ALL<-reshape2::dcast(YASUNI_F,
                                   Day_Year_Y~variable,value.var='value',
                                   direction="wide")

#Used for actual table
YASUNI_DT <- data.table(YASUNI_WIDE_F_ALL)


##Look for Urticaceae
URT_Code <- subset(FAMILY_CODE,
                   FAMILY_CODE $Family.APG.IV. == "Urticaceae")$CODIGO

YAS_URT_SP_TS <- subset(YASUNI_F[YASUNI_F$variable %in% URT_Code, ])



ggplot(YAS_URT_SP_TS, aes (x = Day_Year_Y, y = log(value+1), color=variable,group=variable))+
  geom_line()+scale_colour_viridis(option='turbo',discrete=TRUE)+
  scale_x_continuous(breaks = seq(5,17.5,5), labels= c(2005,2010,2015))+
  xlab("Years")+ylab("Abundance (log + 1)")+ggtitle("Yasuni - Urticaceae") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y=element_text(size = 14,color='black'),
        axis.title.x =element_text(size = 16,color='black'),
        axis.title.y = element_text(size = 16,color='black'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



URT_Count <- cbind.data.frame(Day_Year_Y = 
                              YASUNI_WIDE_F_ALL$Day_Year_Y, 
                              Seed_Count  = 
                                aggregate(YAS_URT_SP_TS[,c("value")],
                                          YAS_URT_SP_TS[,c("Day_Year_Y")],sum)
)

ggplot(URT_Count, aes( x= Day_Year_Y, y=log(Seed_Count.value+1)
))+geom_line()+
  scale_x_continuous(breaks = seq(5,17.5,5),
                     expand = c(0, 0),
                     labels= c(2005,2010,2015))+
  xlab("Year")+
  ylab("Total seed count (log10)" )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=20))
