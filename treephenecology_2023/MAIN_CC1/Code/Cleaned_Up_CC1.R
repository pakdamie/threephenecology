#######################################
#This function is to clean up the     #
#main CC1 time series file for        #
#further analyses. This also ensures  #
#the time scale is yearly.            #
#All analyses for CC1 will use this   #
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

###Following the same example as the lrlake from mvcwt package

main.dat_CC1 <-read.table(here("CC1", "Data", "PE_CC1_1_active_seedtraps_12May2020.txt"),header = TRUE)

###Seem to have the same dates for many of the entries- aggregate by date
# Create a date column with the month day and year provided

###Combine the month, day, year together to make a date 
###data column
main.dat_CC1$Date<-paste(main.dat_CC1$month,
                         main.dat_CC1$day, 
                         main.dat_CC1$year, sep='-') 

main.dat_CC1$Date<-as.Date(main.dat_CC1$Date, format = "%m-%d-%Y") 

#Not everything is in order so order everything by date because it is not 
#in order
main.dat_CC1<- main.dat_CC1[order(main.dat_CC1$Date),]

# head(main.dat_CC1[,1:4])
#trapnum month day seed.yr
#  138     9   9       1
#   276     9   9       1
#   134     9   9       1
#3      27     9   9       1

#When we take a look, similar sampling day, but trap number differ.
#For right now, ignoring trap num and aggregating it based on date,
#function is to sum everything up

CC1_Agg <-  aggregate(main.dat_CC1[,6:727], by= list(main.dat_CC1$Date),'sum') 

#CC1_Agg <- aggregate(main.dat_CC1[,6:727], by= list(main.dat_CC1$Date),'sum') 
#I need to find the day length 
#The Julian Day of 9/9/2002 is 252 
#We can start counting from there to convert
#this into a yearly scale

CC1_Agg <- CC1_Agg %>%
  mutate(Day_Length =c(252,diff(Group.1)),
         Day_Year = cumsum(Day_Length),
         Day_Year_Y = (Day_Year)/365.25)

###This removes any species where all the species entries are 0
CC1_Agg<-CC1_Agg[,colSums(CC1_Agg!= 0)!=0]

###This should have 658 columns or 654 species excluding
###date columns

CC1_F <- CC1_Agg

write.csv(CC1_F, file = here('CC1','Data',"CC1_ALL_FINAL.csv"))
