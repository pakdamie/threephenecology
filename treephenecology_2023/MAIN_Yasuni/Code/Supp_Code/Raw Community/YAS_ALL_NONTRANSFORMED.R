#####################################
### Yasuni- ALL COMMUNITY ANALYSIS###
#####################################

##############################################
#This function is for the total              #
#community analysis: Use YASUNI_ALL_FINAL.csv#
##############################################

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

###The file needed [-1] because the rows are included
YASUNI_F <- fread(here("YASUNI","Data","YASUNI_ALL_FINAL.csv"))[,-1]

###Cast this into wide format 
YASUNI_WIDE_F_ALL<-reshape2::dcast(YASUNI_F,
                          Day_Year_Y~variable,value.var='value',
                          direction="wide")

tmp <- data.table(YASUNI_WIDE_F_ALL)

###The x (time)
x = YASUNI_WIDE_F_ALL$Day_Year_Y
###The y (abundance that is not log transformed)
y =YASUNI_WIDE_F_ALL[,2:ncol(YASUNI_WIDE_F_ALL)]

###Run the multivariate continuous wavelet transformation-
###Let the package choose the minimum scale and the maximum scale
MVCWT_YAS_ALL <-mvcwt(x,y, min.scale=get.min.scale(x),
                      max.scale=get.max.scale(x))
###Get the WMR 
MVCWT_YAS_WMR<- wmr(MVCWT_YAS_ALL)
###This is the bootstrapping-phase randomized
MVCWT_YAS_WMR_BOOT<- wmr.boot(MVCWT_YAS_ALL, smoothing = 1, reps = 1000, mr.func = "wmr")

MVCWT_YAS_WMR_BOOT_NONLOG_2<- getBootWMR(MVCWT_YAS_WMR  , reps = 1000)

##############################################
###Make the figures ###
#######################
image.mvcwt2(MVCWT_YAS_WMR_BOOT,
             z.fun = "Re", 
             bound = 1,reset.par=FALSE,
           legend.only = TRUE,ylim =c(0.08,8.5))
contour(MVCWT_YAS_WMR,add=TRUE,ylim=c(0.08,8.5))

#contour(MVCWT_YAS_WMR_BOOT$x,
      #  MVCWT_YAS_WMR_BOOT$y, MVCWT_YAS_WMR_BOOT$z.boot[,,1], 
      #  levels = c(0.025, 0.975), lwd = 3, add = T,drawlabels = F)
#ontour(MVCWT_YAS_WMR,add=TRUE)

##############################################
##############################################

##################################################
###Total species count and seed count figures ####
##################################################

YAS_Count <- cbind.data.frame(Day_Year_Y = 
                              
                                
                                YASUNI_WIDE_F_ALL$Day_Year_Y, 
                              Species_Count = rowSums(YASUNI_WIDE_F_ALL[,-1]  !=0),
                              Seed_Count  =  (rowSums(YASUNI_WIDE_F_ALL)))



YAS_SC_GG <- ggplot(YAS_Count, aes(x= Day_Year_Y, y =Species_Count ))+
  geom_line(color='darkgreen')+theme_bw()+
  scale_x_continuous(breaks = seq(5,17.5,5),expand = c(0, 0),  labels= c(2005,2010,2015))+
  ylim(0,150)+
  xlab("Year")+
  ylab("Total species count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=20))


YAS_PC_GG <- ggplot(YAS_Count, aes(x = Day_Year_Y, y= log(Seed_Count+1)))+geom_line()+theme_bw()+
  ylim(0, 15)+
  
  scale_x_continuous(breaks = seq(5,17.5,5),expand = c(0, 0),
                     labels= c(2005,2010,2015))+
  xlab("Year")+
  ylab("Total seed count (log)" )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=20))


(YAS_SC_GG +CC1_SP_COUNT_GG_2) /(YAS_PC_GG+CC1_PHEN_SUM_GG)

##########################################
###AVERAGE WMR ACROSS TIME SQUARE PLOT ###
##########################################

YAS_AVG_WMR =data.frame (time=MVCWT_YAS_WMR$y,WMR=colMeans(MVCWT_YAS_WMR$z[,,1]))
View(MVCWT_YAS_WMR_BOOT_2_nonlog[[1]])

###So each list element represents the period
###each element's row represent a time point and the 
###each element's column represent a replicate

scale_window <- NULL
for (i in seq(1, length(MVCWT_YAS_WMR_BOOT_2_nonlog))){
  tmp = MVCWT_YAS_WMR_BOOT_2_nonlog[[i]]
  mean <- rowMeans(tmp)

  lb <-quantile(mean, probs= c(0.025),na.rm=TRUE)
  ub <-quantile(mean, probs= c(0.975),na.rm=TRUE) 

  scale_window[[i]]=data.frame(scale=MVCWT_YAS_WMR$y[[i]],lb,ub)

  }
scale_window <-do.call(rbind, scale_window)


YAS_GG <- ggplot(YAS_AVG_WMR, aes(x=time,y=WMR))+
  geom_line(size=1.2)+ylim(0,1)+theme_bw()+
  xlab("Period (Years)") + ylab("Avg. WMR")+
  geom_ribbon(data=scale_window,aes(ymin=lb,ymax=ub))
CC1_GG<- ggplot(CC1_AVG_WMR, aes(x=time,y=WMR))+geom_line(size=1.2)+ylim(0,1)+theme_bw()+
  xlab("Period (Years)") + ylab("Avg. WMR")



##########################################
### AVERAGE WMR ACROSS TIME SQUARE PLOT ###
##########################################

###This ensures that anything out of the cone of influence is 
###not included for the anaysis.

MVCWT_YAS_WMR_COI <- matrix(NA, ncol = ncol(MVCWT_YAS_WMR$z[, , 1]), nrow = nrow(MVCWT_YAS_WMR$z[, , 1]))

for (i in 1:ncol(MVCWT_YAS_WMR_COI)) MVCWT_YAS_WMR_COI[which(MVCWT_YAS_WMR$x > (min(MVCWT_YAS_WMR$x) + 1 * MVCWT_YAS_WMR$y[i]) & MVCWT_YAS_WMR$x < (max(MVCWT_YAS_WMR$x) - 1* MVCWT_YAS_WMR$y[i])), i] <- 1
MVCWT_YAS_WMR_COI0 <- MVCWT_YAS_WMR_COI



yas_nontransformed_CI<-NULL
for (k in seq(1, length(MVCWT_YAS_WMR_BOOT_2_nonlog))){
  tmp <-  MVCWT_YAS_WMR_BOOT_2_nonlog[[k]] 
  new_tmp <- colMeans(tmp * MVCWT_YAS_WMR_COI0[,k],na.rm=TRUE)
  
  
  lowerbound <-quantile(new_tmp, probs= c(0.025),na.rm=TRUE)
  upperbound <- quantile(new_tmp, probs= c( 0.975),na.rm=TRUE)
  
  yas_nontransformed_CI[[k]]<- cbind(lb = lowerbound, up = upperbound)
}

yas_nontransformed_CI<- do.call(rbind,yas_nontransformed_CI)



###################
###The average WMR#
###################
YAS_AVG_WMR_nontransformed <- data.frame(Scale = MVCWT_YAS_WMR$y,
                          WMR = colMeans(MVCWT_YAS_WMR_COI0 * MVCWT_YAS_WMR$z[,,1],na.rm=TRUE),
                          lb = yas_nontransformed_CI[,1],
                          ub = yas_nontransformed_CI[,2])


YAS_AVG_WMR_nontransformed$sig <- ifelse(YAS_AVG_WMR_nontransformed$WMR >YAS_AVG_WMR_nontransformed$ub
                          |YAS_AVG_WMR_nontransformed$WMR <YAS_AVG_WMR_nontransformed$lb, 'sig','notsig')



#######################################
### Creating the confidence interval###
#######################################

YAS_GG_nontransformed <- ggplot(YAS_AVG_WMR_nontransformed,
                 aes(x = Scale, y = WMR)) +
  geom_line(size = 1.1,aes(color=sig,group=1)) +
  geom_ribbon(aes(ymin= lb, ymax=ub),alpha = 0.1)+
  scale_color_manual(values=c('black','#FFB61D'))+
  ylim(0, 1) +
  theme_bw() +
  xlab("Period (Years)") +
  ylab("Avg. WMR")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
