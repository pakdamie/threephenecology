################################################################
### CC1- Community all                                        #
################################################################


##############################################
#This function is for the total              #
#community analysis: Use CC1_ALL_FINAL.csv#
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
CC1_F <- fread(here("CC1","Data","CC1_ALL_FINAL.csv"))[,-1]
source(here("SUPP","image_2.R"))

###FOR THE MVCWT PACKAGE
x_cc1 = CC1_F$Day_Year_Y
y_cc1 = CC1_F[,2:655]
ylog_cc1 = log(CC1_F[,2:655]+1)

CC1_MVCWT_nonlog = mvcwt(x_cc1, y_cc1, min.scale = get.min.scale(x_cc1), 
                      max.scale = get.max.scale(x_cc1))

CC1_MVCWT_log = mvcwt(x_cc1, ylog_cc1, min.scale = get.min.scale(x_cc1), 
                  max.scale = get.max.scale(x_cc1))
CC1_WMR_nonlog = wmr(CC1_MVCWT_nonlog)

CC1_WMR_log = wmr(CC1_MVCWT_log)


###BOOTSTRAP###
CC1_WMR_BOOT<- wmr.boot(CC1_MVCWT, 
                        reps = 1000,smoothing = 1,  mr.func = "wmr")

save(CC1_WMR_BOOT, file = here("CC1","Output", "CC1_TRANSFORMED","CC1_WMR_BOOT_TRANSFORMED.RData"))

MVCWT_CC1_WMR_BOOT_NONLOG_2<- getBootWMR(CC1_WMR_nonlog  , reps = 1000)
MVCWT_CC1_WMR_BOOT_LOG_2<- getBootWMR(CC1_WMR_log  , reps = 1000)

#######################
###Making the plots ###
#######################

image.mvcwt2(CC1_WMR_BOOT,
             z.fun = "Re", 
             bound = 1,
             reset.par=FALSE,
             ylim =c(0.08,8.5))

contour(CC1_WMR_log,
        add=TRUE, 
        ylim =c(0.08,8.5))


image.mvcwt2(CC1_WMR_BOOT,
             z.fun = "Re", 
             bound = 1,
             reset.par=FALSE,
             ylim =c(0.08,8.5))

contour(CC1_WMR_log,
        add=TRUE, 
        ylim =c(0.08,8.5))


#####################################
###TOTAL SPECIES- ACROSS THE YEARS###
#####################################

CC1_FIN <- CC1_F %>%
          select(!c(Group.1,
                    Day_Length,
                    Day_Year
                  ))

CC1_Count <- cbind.data.frame(Day_Year_Y = 
                                CC1_FIN$Day_Year_Y, 
                              Species_Count = rowSums(CC1_FIN[,-1]  !=0),
                              Seed_Count  =  (rowSums(CC1_FIN)))


high_species_name <- names(sort(colSums(CC1_FIN[,-c("Day_Year_Y")]), decreasing=TRUE)[1:5])

cc1_fin_01<- melt(data.frame(log(CC1_FIN [,..high_species_name ]+0.1)))
qqnorm(cc1_fin_01$value)
qqline(cc1_fin_01$value,col='red')


cc1_fin_1<- melt(data.frame(log(CC1_FIN [,..high_species_name ]+1)))
qqnorm(cc1_fin_1$value)
qqline(cc1_fin_1$value,col='red')

cc1_fin_10<- melt(data.frame(log(CC1_FIN [,..high_species_name ]+10)))
qqnorm(cc1_fin_10$value)
qqline(cc1_fin_10$value,col='red')

###This is for the total Species Count
YAS_SC_GG <- ggplot(YAS_Count, 
                    aes(x= Day_Year_Y, y =Species_Count ))+
  geom_line(color='darkgreen')+
  scale_x_continuous(breaks = seq(5,17.5,5),expand = c(0, 0),  
                     labels= c(2005,2010,2015))+
  xlab("Year")+
  ylab("Total species count")+theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=20),
        axis.text.x = element_text(color='black'),
        axis.text.y = element_text(color='black'))

###This is for the total seed count
YAS_PC_GG <- ggplot(YAS_Count, aes(x = Day_Year_Y, y= log10(Seed_Count+1)))+
  geom_line()+
  scale_x_continuous(breaks = seq(5,17.5,5),
                     expand = c(0, 0),
                     labels= c(2005,2010,2015))+
  xlab("Year")+
  ylab("Total seed count (log10)" )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=20))

###Average across scale with null window

MVCWT_CC1_WMR_COI <- matrix(NA, ncol = ncol(CC1_WMR_log$z[, , 1]), nrow = nrow(CC1_WMR_log$z[, , 1]))

for (i in 1:ncol(MVCWT_CC1_WMR_COI)) MVCWT_CC1_WMR_COI[which(CC1_WMR_log$x > (min(CC1_WMR_log$x) + 1 * CC1_WMR_log$y[i]) & CC1_WMR_log$x < (max(CC1_WMR_log$x) - 1 * CC1_WMR_log$y[i])), i] <- 1
MVCWT_CC1_WMR_COI0 <- MVCWT_CC1_WMR_COI




CC1_transformed_CI<-NULL
CC1_untransformed_CI<-NULL

for (k in seq(1, length(MVCWT_CC1_WMR_BOOT_2_log))){
  tmp <-  MVCWT_CC1_WMR_BOOT_2_log[[k]] 
  new_tmp <- colMeans(tmp * MVCWT_CC1_WMR_COI0[,k],na.rm=TRUE)
  
  tmp2 <- MVCWT_CC1_WMR_BOOT_2_nonlog[[k]] 
  new_tmp_2 <- colMeans(tmp2 * MVCWT_CC1_WMR_COI0[,k],na.rm=TRUE)
  
  
  lowerbound <-quantile(new_tmp, probs= c(0.025),na.rm=TRUE)
  upperbound <- quantile(new_tmp, probs= c( 0.975),na.rm=TRUE)
  
  lowerbound_2 <-quantile(new_tmp_2, probs= c(0.025),na.rm=TRUE)
  upperbound_2 <- quantile(new_tmp_2, probs= c( 0.975),na.rm=TRUE)
  
  
 CC1_transformed_CI[[k]]<- cbind(lb = lowerbound, up = upperbound)
 CC1_untransformed_CI[[k]]<- cbind(lb = lowerbound_2, up = upperbound_2)
 
 }

CC1_transformed_CI<- do.call(rbind, CC1_transformed_CI)
CC1_untransformed_CI<- do.call(rbind, CC1_untransformed_CI)


###The average WMR 
CC1_AVG_WMR <- data.frame(Scale = CC1_WMR_log$y,
                          WMR = colMeans(MVCWT_CC1_WMR_COI0 * CC1_WMR_log$z[,,1],na.rm=TRUE),
                          lb = CC1_transformed_CI[,1],
                          ub = CC1_transformed_CI[,2])

CC1_AVG_WMR_nonlog <- data.frame(Scale = CC1_WMR_nonlog$y,
                          WMR = colMeans(MVCWT_CC1_WMR_COI0 * CC1_WMR_nonlog$z[,,1],na.rm=TRUE),
                          lb = CC1_untransformed_CI[,1],
                          ub = CC1_untransformed_CI[,2])



CC1_AVG_WMR$sig <- ifelse(CC1_AVG_WMR$WMR >CC1_AVG_WMR$ub
                          |CC1_AVG_WMR$WMR <CC1_AVG_WMR$lb, 'sig','notsig')

CC1_AVG_WMR_nonlog$sig <- ifelse(CC1_AVG_WMR_nonlog$WMR >CC1_AVG_WMR_nonlog$ub
                          |CC1_AVG_WMR_nonlog$WMR <CC1_AVG_WMR_nonlog$lb, 'sig','notsig')

### Creating the confidence interval

CC1_GG <- ggplot(na.omit(CC1_AVG_WMR),
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


CC1_nonlog_GG <- ggplot(na.omit(CC1_AVG_WMR_nonlog),
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
