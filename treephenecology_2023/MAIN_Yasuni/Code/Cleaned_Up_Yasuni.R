
#######################################
#This function is to clean up the     #
#main Yasuni time series file for     #
#further analyses- specifically,      #
#getting rid of the censored species  #
#as well as certain species that      #
#Nancy suggested to remove from       #
#analyses. This also ensures that     #
#the time scale is yearly.            #
#All analyses for Yasuni will use this#
#time series                          #
#######################################



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

###Files needed
CODE <- read.csv(here("YASUNI", "Data", "CODE_YAS_Lasky_ncg.csv"))
YAS_REMOVED_NANCY <-read.csv(here("YASUNI","Data","11_16_NANCY_REMOVED.csv"))

# Taxon	Level	End Date for censoring
#--------------------------------------
# Ficus	Genus	1/1/2008
# Moraceae	Family 	1/1/2008
# Solanum	Genus	1/1/2007
# Solanaceae	Family 	1/1/2007
# Clusia	Genus	1/1/2007
# Clusiaceae	Family	1/1/2007
# Asteraceae	Family	Exclude
# Araceae	Family	Exclude

###############################################
# These are the censored species that we are  #
# excluding from all analyses                 #
###############################################
CENSORED_SP <- c("Moraceae","Solanaceae","Clusiaceae",
                 "Clusiaceae?","Asteraceae","Araceae")
###################################################################
###This gives me the code for the species that should be excluded##
###################################################################
CENSORED_SP_DF <- CODE %>%
  filter(Family.APG.IV. %in% CENSORED_SP) %>%
  select(CODIGO)

##############################################################################
### The Original Yasuni Time Series ##########################################
##############################################################################
YASUNI <- read.csv(here("YASUNI", "Data", "YasuniTS.csv"))
#Ensure that the first colum is titled Date
colnames(YASUNI)[1] <- "Date"

# Make into actual date format
YASUNI$Date <- as.Date(YASUNI$Date, format = "%m/%d/%Y")

### Seems like each date is unique and in order!
### 2000-02-28 is the beginning of the time series such that (59) is 
### the Julian day. So we're starting at day 59.

###Putting everything into a yearly scale
###So 365 days = 1, 365*2 = 2, and so on.
YASUNI <- YASUNI %>% 
  mutate(Day_Length = c(59, diff(Date)), 
         Day_Year = cumsum(Day_Length),
         Day_Year_Y = (Day_Year) / 365.25)
###In total should have 1745 columns, but only
### 1741 should be the species columns (aka
### excluding the date columns)
###  
### Gets rid of any species where there are no entries -
###  should be 1198 total columns (1194 species)
YASUNI <- YASUNI[, colSums(YASUNI != 0) != 0]

###This will melt the YASUNI dataframe 
###EXCLUDING the Date (Because we're using
###Day_Year_Y as the variable to calculate the scales)
###This will create a dataframe, with the Day_Year_Y, with the variable name
###(The Species name) and the value.

YASUNI_MELT <- reshape2::melt(YASUNI[, -which(names(YASUNI) %in% 
                                      c("Day_Length","Day_Year"))],
                       id.vars = c("Day_Year_Y", "Date")
)
### This removes all species that have the code that pertains to the 
###censored family 
###
YASUNI_MELT <- YASUNI_MELT %>%
  filter(!(variable %in% CENSORED_SP_DF$CODIGO)) 

###This more convoluted way was because previous analyses tried
###incorporate censored species unlike our current method where we just
###remove them entirely
###########################################
###1072 species EXCLUDING censored species#
###########################################

##############################################################
##Also there are some removed species that Nancy provided#####
##later down the line                                    #####
##############################################################

YASUNI_F <- YASUNI_MELT %>%  
  filter(!(variable %in% YAS_REMOVED_NANCY$CODIGO))


length(unique(YASUNI_F $variable))

##########################################
###This should now bring it down to 1054##
###species for all analyses             ##
##########################################

###This function only takes the 
###time scale we are using, the species name,
### and the value. 
YASUNI_F <- YASUNI_F %>%
  select(Day_Year_Y, variable,value)
###################################
###YASUNI_F should be the analysis# 
###################################

write.csv(YASUNI_F, file = 'YASUNI_ALL_FINAL.csv')
