####################################################
###Environmental analysis of the CC1             ###
###- looking at the minimum temperature          ###
###and precipitation- this was last updated      ###
###in 5/5/2022 and now includes the log/unlogged ###
####################################################

###All new code

###Packages required###
library(here) # To deal with the working directory
library(reshape2) # To deal with the melting long to wide/wide to long
library(mvcwt) # not on CRAN, get the updated version from the graveyard
library(data.table)
library(dplyr)
library(WaveletComp)
library(ggplot2)
library(gridExtra)
library(stringr)
library(tidyverse)
library(patchwork)
###Following the same example as the lrlake from mvcwt package

main.dat_CC1 <-read.table(here("CC1", "Data", "PE_CC1_1_active_seedtraps_12May2020.txt"),header = TRUE)

###Seem to have the same dates- aggregate by date
# Create a date column with the month day and year provided

main.dat_CC1$Date<-paste(main.dat_CC1$month,
                         main.dat_CC1$day, 
                         main.dat_CC1$year, sep='-') 

main.dat_CC1$Date<-as.Date(main.dat_CC1$Date, format = "%m-%d-%Y") 

#Not everything is in order so order everything by date
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
#The Julian Day of 9/9/2002 is at the 252 
#We can start counting from there?

CC1_Agg <- CC1_Agg %>%
  mutate(Day_Length =c(252,diff(Group.1)),
         Day_Year = cumsum(Day_Length),
         Day_Year_Y = (Day_Year)/365.25)


CC1_Agg<-CC1_Agg[,colSums(CC1_Agg!= 0)!=0]

#I need to find the day length 
#The Julian Day of 9/9/2002 is at the 252 
#We can start counting from there?

CC1_Agg$Month <- as.numeric(format(CC1_Agg$Group.1, '%m'))
CC1_Agg$Year <- as.numeric(format(CC1_Agg$Group.1, '%Y'))

CC1_Agg <- aggregate(CC1_Agg[,-1], by = list(CC1_Agg$Month, 
                                             CC1_Agg$Year),'mean')
top<- CC1_Agg[1:63,]
missing<- NA
bottom<- CC1_Agg[64:99,]

CC1_Agg_2 <- rbind(top,missing)
CC1_Agg_3 <- rbind(CC1_Agg_2, bottom)

CC1_Agg_3[,3:ncol(CC1_Agg_3 )] <- na.approx(CC1_Agg_3 [,3:ncol(CC1_Agg_3 )])


Actual_Date <- c(seq(1,100))

CC1_Agg_3$Date = Actual_Date

CC1_FINAL = CC1_Agg_3[,!(colnames(CC1_Agg_3) %in% c("Month","Year","Group.1","Group.2",
                                                    "Day_Year","Day_Year_Y","Day_Length"
                                                ))]


x= CC1_FINAL$Date
### For the y, you must rely on the species data from the 2nd column to the end of the
###column- LOG TRANSFORMED
y_nonlog=  CC1_FINAL[,1:654]
y_log=  log(CC1_FINAL[,1:654]+1)

###This does the wavelet transformation for each column/species
MVWCT_FAM_nonlog = mvcwt(x, y_nonlog, scales =seq(2,50))
MVWCT_FAM_log = mvcwt(x, y_log, scales =seq(2,50))

###After that is accomplished, we now do the wavelet modulus ratio
###
WMR_FAM_PERU_MONTH_nonlog= wmr(MVWCT_FAM_nonlog)
WMR_FAM_PERU_MONTH_log= wmr(MVWCT_FAM_log)


#save(WMR_FAM_UNCENSORED_PERU_MONTH, file=  "WMR_MONTH_PERU.RData")
#####################################################
#####################################################

TMIN_CC1 <- read.table(here("CC1","Data","ECRtmin_PE_CC1.txt"),header=TRUE)
PRECIP_CC1 <- read.table(here("CC1","Data","CHIRP_PE_CC1.txt"),header=TRUE)

TMIN_CC1_F <- NULL
for (i in seq(1,15)){
  df <- data.frame(tmin = t(TMIN_CC1 [i,])[-1],
                   year = TMIN_CC1[i,][1],
                   month = as.numeric(seq(1,12)))
  
  TMIN_CC1_F[[i]] <- df  
  
}
TMIN_CC1_F<-do.call(rbind, TMIN_CC1_F)
TMIN_CC1_F <- na.omit(TMIN_CC1_F)

TMIN_CC1_F<-TMIN_CC1_F[33:133,]
#######################################################


###PRECIP WIDE TO LONG
PRECIP_CC1_F <- NULL
for (i in seq(1,34)){
  df <-data.frame(pre = t(PRECIP_CC1[i,])[-1], year = PRECIP_CC1[i,][1],month = seq(1,12))
  
  PRECIP_CC1_F[[i]] <- df  
  
}
PRECIP_CC1_F<-do.call(rbind, PRECIP_CC1_F)
PRECIP_CC1_F <- na.omit(PRECIP_CC1_F)

PRECIP_CC1_F <-PRECIP_CC1_F [261:361,]

############################################################


###COMBINED ENVIRONMENTAL VARIABLES

ENVIR_CC1 <- cbind.data.frame(TMIN_CC1_F, precip = PRECIP_CC1_F$pre, id = 'CC1')
ENVIR_CC1<- ENVIR_CC1[-1,]

ENVIR_CC1$jul <- seq(1, nrow(ENVIR_CC1))

###############################################################


TMIN = data.frame(tmin =ENVIR_CC1$tmin)
TEMP_WAV <- analyze.wavelet(TMIN,'tmin',
                            loess.span = 0,dt = 1, lowerPeriod =2, upperPeriod = 50)
wt.image(TEMP_WAV)

Period = seq(2,50)



M_CC1_COI <-  matrix(NA, ncol = ncol(WMR_FAM_PERU_MONTH_log$z[, , 1]), nrow = nrow(WMR_FAM_PERU_MONTH_log$z[, , 1]))

for (i in 1:ncol(M_CC1_COI )) M_CC1_COI [which(WMR_FAM_PERU_MONTH_log$x > (min(WMR_FAM_PERU_MONTH_log$x) + 1 * WMR_FAM_PERU_MONTH_log$y[i]) & WMR_FAM_PERU_MONTH_log$x < (max(WMR_FAM_PERU_MONTH_log$x) - 1 * WMR_FAM_PERU_MONTH_log$y[i])), i] <- 1
M_CC1_COI0<-M_CC1_COI  

WMR_FAM_PERU_MONTH_nonlog$z[,,1]<- WMR_FAM_PERU_MONTH_nonlog$z[,,1]*M_CC1_COI0
WMR_FAM_PERU_MONTH_log$z[,,1]<- WMR_FAM_PERU_MONTH_log$z[,,1]*M_CC1_COI0

#################This is where the fun begin##########################
TEST_RESULTS_TMIN_NONLOG=NULL
TEST_RESULTS_TMIN_LOG=NULL

for (i in seq(1,49)){
  RECONSTRUCTED =reconstruct(TEMP_WAV,sel.period = Period[i],
                             only.sig=FALSE)
  recon_tmin =RECONSTRUCTED$series$tmin.r
  
  CC1_WMR_SCALE_NONLOG=WMR_FAM_PERU_MONTH_nonlog$z[1:100,i,1]
  CC1_WMR_SCALE_LOG=WMR_FAM_PERU_MONTH_log$z[1:100,i,1]
  
  tempo_NONLOG = cor.test(recon_tmin, CC1_WMR_SCALE_NONLOG,method='pearson',na.rm=TRUE)
  tempo_LOG = cor.test(recon_tmin, CC1_WMR_SCALE_LOG,method='pearson',na.rm=TRUE)
  
  df_temp_NONLOG <- cbind.data.frame(Scale = Period[i],t= tempo_NONLOG$estimate) 
  df_temp_LOG <- cbind.data.frame(Scale = Period[i],t= tempo_LOG$estimate) 
  
  TEST_RESULTS_TMIN_NONLOG[[i]] = df_temp_NONLOG
  TEST_RESULTS_TMIN_LOG[[i]] = df_temp_LOG
  
  }


########
########
TEST_RESULT_TMIN_2_NONLOG <- do.call(rbind, TEST_RESULTS_TMIN_NONLOG)
TEST_RESULT_TMIN_2_LOG <- do.call(rbind, TEST_RESULTS_TMIN_LOG)





PHASE_RANDOMIZED_ALL_NONLOG = NULL
PHASE_RANDOMIZED_ALL_LOG = NULL

for (j in seq(1,46)){
  #Pull out a specific scale 
  RECONSTRUCTED =reconstruct(TEMP_WAV,sel.period = Period[j],
                             only.sig=FALSE)
  #Pulls out the reconstructed time series 
  recon_tmin =data.frame(jul=seq(1,100),
                         tmin=RECONSTRUCTED$series$tmin.r)
  #Sample where the time series would start (1 to 169)
  random_starter = sample(1:100,1000,replace =T)
  
  CC1_WMR_SCALE_NONLOG=WMR_FAM_PERU_MONTH_nonlog$z[1:100,j,1]
  CC1_WMR_SCALE_LOG=WMR_FAM_PERU_MONTH_log$z[1:100,j,1]
  
  
  cor_1000_NONLOG = NULL
  cor_1000_LOG = NULL
  
  for (i in seq(random_starter)){
    
    start = recon_tmin[recon_tmin$jul>=random_starter[i],]
    end = recon_tmin[recon_tmin$jul < random_starter[i],]
    connected = rbind(start,end)
    
    
    tmp_NONLOG<- cor.test(connected$tmin, 
                   CC1_WMR_SCALE_NONLOG,method='kendall',na.rm=TRUE)$estimate
    
    tmp_LOG<- cor.test(connected$tmin, 
                          CC1_WMR_SCALE_LOG,method='kendall',na.rm=TRUE)$estimate
    
    cor_1000_NONLOG[[i]]= tmp_NONLOG
    cor_1000_LOG[[i]]= tmp_LOG
    
  }
  PHASE_RANDOMIZED_ALL_NONLOG[[j]]=unlist(cor_1000_NONLOG)
  PHASE_RANDOMIZED_ALL_LOG[[j]]=unlist(cor_1000_LOG)
  
}

Summary_Stats_PHASE_NONLOG= NULL
Summary_Stats_PHASE_LOG= NULL

for (k in seq(1,46)){
  tmp_NONLOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_NONLOG[[k]]),
                       lower = quantile(PHASE_RANDOMIZED_ALL_NONLOG[[k]], probs= c(0.025),na.rm=TRUE),
                       upper = quantile(PHASE_RANDOMIZED_ALL_NONLOG[[k]], probs= c(0.975),na.rm=TRUE),
                       scale = Period[k])
 
   Summary_Stats_PHASE_NONLOG[[k]]=tmp_NONLOG
  
   tmp_LOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_LOG[[k]]),
                               lower = quantile(PHASE_RANDOMIZED_ALL_LOG[[k]], probs= c(0.025),na.rm=TRUE),
                               upper = quantile(PHASE_RANDOMIZED_ALL_LOG[[k]], probs= c(0.975),na.rm=TRUE),
                               scale = Period[k])
   
   Summary_Stats_PHASE_NONLOG[[k]]=tmp_NONLOG
   Summary_Stats_PHASE_LOG[[k]]=tmp_LOG
   
   
   
  #Make the p-value 
  tmp_NONLOG$p_value <- 0
  tmp_LOG$p_value <- 0
  
  #Count 
  count1_NONLOG =table( PHASE_RANDOMIZED_ALL_NONLOG[[k]] >= TEST_RESULT_TMIN_2_NONLOG $t[[k]])
  count2_NONLOG =table(  PHASE_RANDOMIZED_ALL_NONLOG[[k]] <= TEST_RESULT_TMIN_2_NONLOG $t[[k]])
  
  count1_LOG =table( PHASE_RANDOMIZED_ALL_LOG[[k]] >= TEST_RESULT_TMIN_2_LOG $t[[k]])
  count2_LOG =table( PHASE_RANDOMIZED_ALL_LOG[[k]] <= TEST_RESULT_TMIN_2_LOG $t[[k]])
  
  
  
  
  p_1_NONLOG <- try(c(2*(count1_NONLOG[[2]]/1000)))
  p_2_NONLOG <- try(c(2 *(count2_NONLOG[[2]]/1000)))
  
  p_1_LOG <- try(c(2*(count1_LOG[[2]]/1000)))
  p_2_LOG <- try(c(2*(count2_LOG[[2]]/1000)))
  
  p_1_NONLOG=ifelse(is.character(p_1_NONLOG)==TRUE,0,p_1_NONLOG)
  p_2_NONLOG=ifelse(is.character(p_2_NONLOG)==TRUE,0,p_2_NONLOG)

  p_1_LOG=ifelse(is.character(p_1_LOG)==TRUE,0,p_1_LOG)
  p_2_LOG=ifelse(is.character(p_2_LOG)==TRUE,0,p_2_LOG)
  
  
  tmp_NONLOG$p_value= min(p_1_NONLOG,p_2_NONLOG)
  tmp_LOG$p_value= min(p_1_LOG,p_2_LOG)
  
  

  
  
  Summary_Stats_PHASE_NONLOG[[k]]=tmp_NONLOG
  Summary_Stats_PHASE_LOG[[k]]=tmp_LOG
  
  
}







Summary_Stats_PHASE_2_NONLOG <- do.call(rbind, Summary_Stats_PHASE_NONLOG)
Summary_Stats_PHASE_2_NONLOG$actual <- TEST_RESULT_TMIN_2_NONLOG$t

Summary_Stats_PHASE_2_LOG <- do.call(rbind, Summary_Stats_PHASE_LOG)
Summary_Stats_PHASE_2_LOG$actual <- TEST_RESULT_TMIN_2_LOG$t


Summary_Stats_PHASE_2_NONLOG$p_value_adjusted <-  p.adjust( Summary_Stats_PHASE_2_NONLOG $p_value, 
                                    method ="BY")
Summary_Stats_PHASE_2_LOG$p_value_adjusted <-  p.adjust( Summary_Stats_PHASE_2_LOG $p_value, 
                                                            method ="BY")
Summary_Stats_PHASE_2_NONLOG$psig <- ifelse(Summary_Stats_PHASE_2_NONLOG$p_value_adjusted  <= 0.05,
                   1,0.5)   
Summary_Stats_PHASE_2_LOG$psig <- ifelse(Summary_Stats_PHASE_2_LOG$p_value_adjusted  <= 0.05,
                                            1,0.5) 
Summary_Stats_PHASE_2_NONLOG$sig <- ifelse( Summary_Stats_PHASE_2_NONLOG$actual>= 
                                       Summary_Stats_PHASE_2_NONLOG$lower &
                                       Summary_Stats_PHASE_2_NONLOG$actual<=Summary_Stats_PHASE_2_NONLOG$upper,
                                     'insig','sig')
Summary_Stats_PHASE_2_LOG$sig <- ifelse( Summary_Stats_PHASE_2_LOG$actual>= 
                                              Summary_Stats_PHASE_2_LOG$lower &
                                              Summary_Stats_PHASE_2_LOG$actual<=Summary_Stats_PHASE_2_LOG$upper,
                                            'insig','sig')
Summary_Stats_PHASE_2_NONLOG$yscale <- as.numeric(round(Summary_Stats_PHASE_2_NONLOG$scale/12,3))
Summary_Stats_PHASE_2_LOG$yscale <- as.numeric(round(Summary_Stats_PHASE_2_NONLOG$scale/12,3))




###########
###PLOTS###
###########

tmin_NONLOG<- 
  ggplot(na.omit(Summary_Stats_PHASE_2_NONLOG), 
         aes(x = (yscale), y = mean))+
  geom_point()+geom_segment(aes(x=(yscale), y=lower, xend=yscale,
                                yend=upper),alpha =0.3,size = 0.5)+
  geom_point(data=na.omit(Summary_Stats_PHASE_2_NONLOG),
             aes(x =yscale, y= actual,
                 fill=as.factor(sig),
                 color = as.factor(psig),
                 size = as.factor(sig),
                 shape = as.factor(psig)))+
  scale_fill_manual(values = c('grey','#3ffc97'))+
  scale_shape_manual(values=c(21,24))+
  scale_size_manual(values=c(2,3))+
  scale_color_manual(values = c('lightgrey','black'))+

  xlab("Scale (Years)") +   
  ylim(-1,1)+
  ylab("WMR")+ggtitle("Min. Temperature")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))

tmin_LOG<- 
  ggplot(na.omit(Summary_Stats_PHASE_2_LOG), 
         aes(x = (yscale), y = mean))+
  geom_point()+geom_segment(aes(x=(yscale), y=lower, xend=yscale,
                                yend=upper),alpha =0.3,size = 0.5)+
  geom_point(data=na.omit(Summary_Stats_PHASE_2_LOG),
             aes(x =yscale, y= actual,
                 fill=as.factor(sig),
                 color = as.factor(psig),
                 size = as.factor(sig),
                 shape = as.factor(psig)))+
  scale_fill_manual(values = c('grey','#3ffc97'))+
  scale_shape_manual(values=c(21,24))+
  scale_size_manual(values=c(2,3))+
  scale_color_manual(values = c('lightgrey','black'))+
 
  xlab("Scale (Years)") +   
  ylim(-1,1)+
  ylab("WMR")+ggtitle("Min. Temperature")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))

##############################################################
##############################################################
########################
###PRECIPITATION DATA###
########################

PRECIP = data.frame(precip =ENVIR_CC1$precip)
PRECIP_WAV <- analyze.wavelet(PRECIP,'precip',
                              loess.span = 0,dt = 1,
                              lowerPeriod =2, upperPeriod = 50)
wt.image(PRECIP_WAV)




TEST_RESULTS_PRECIP_NONLOG=NULL
TEST_RESULTS_PRECIP_LOG=NULL

for (i in seq(1,49)){
  RECONSTRUCTED =reconstruct(PRECIP_WAV,sel.period = Period[i],
                             only.sig=FALSE)
  recon_precip =RECONSTRUCTED$series$precip.r
  
  CC1_WMR_SCALE_NONLOG=WMR_FAM_PERU_MONTH_nonlog$z[1:100,i,1]
  CC1_WMR_SCALE_LOG=WMR_FAM_PERU_MONTH_log$z[1:100,i,1]
  
  tempo_NONLOG = cor.test(recon_precip, CC1_WMR_SCALE_NONLOG,method='pearson',na.rm=TRUE)$estimate
  tempo_LOG = cor.test(recon_precip, CC1_WMR_SCALE_LOG,method='pearson',na.rm=TRUE)$estimate
  
  df_precip_NONLOG <- cbind.data.frame(Scale = Period[i],t(tempo_NONLOG))
  df_precip_LOG <- cbind.data.frame(Scale = Period[i],t(tempo_LOG))
  
  TEST_RESULTS_PRECIP_NONLOG[[i]] = df_precip_NONLOG
  TEST_RESULTS_PRECIP_LOG[[i]] = df_precip_LOG
  
}


########
########



TEST_RESULT_PRECIP_2_NONLOG <- do.call(rbind, TEST_RESULTS_PRECIP_NONLOG)
TEST_RESULT_PRECIP_2_LOG <- do.call(rbind, TEST_RESULTS_PRECIP_LOG)


PHASE_RANDOMIZED_ALL_PRECIP_NONLOG= NULL
PHASE_RANDOMIZED_ALL_PRECIP_LOG= NULL

for (j in seq(1,49)){
  #Pull out a specific scale 
  RECONSTRUCTED =reconstruct(PRECIP_WAV,sel.period = Period[j],
                             only.sig=FALSE)
  #Pulls out the reconstructed time series 
  recon_precip =data.frame(jul=seq(1,100),
                           precip=RECONSTRUCTED$series$precip.r)
  #Sample where the time series would start (1 to 169)
  random_starter = sample(1:100,1000,replace =T)
  
  
  CC1_WMR_SCALE_NONLOG=WMR_FAM_PERU_MONTH_nonlog$z[1:100,j,1]
  CC1_WMR_SCALE_LOG=WMR_FAM_PERU_MONTH_log$z[1:100,j,1]
  
  cor_1000_NONLOG= NULL
  cor_1000_LOG= NULL
  
  for (i in seq(random_starter)){
    
    start = recon_precip[recon_precip$jul>=random_starter[i],]
    end = recon_precip[recon_precip$jul < random_starter[i],]
    connected = rbind(start,end)
    
    
    tmp_NONLOG<- cor.test(connected$precip, 
                   CC1_WMR_SCALE_NONLOG,method='pearson',na.rm=TRUE)$estimate
    tmp_LOG<- cor.test(connected$precip, 
                          CC1_WMR_SCALE_LOG,method='pearson',na.rm=TRUE)$estimate
    
    
    cor_1000_NONLOG[[i]]= tmp_NONLOG
    cor_1000_LOG[[i]] = tmp_LOG 
  }
  PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[j]]=unlist(cor_1000_NONLOG)
  PHASE_RANDOMIZED_ALL_PRECIP_LOG[[j]]=unlist(cor_1000_LOG)
  
}

Summary_Stats_PHASE_PRECIP_NONLOG= NULL
Summary_Stats_PHASE_PRECIP_LOG= NULL

for (k in seq(1,49)){
  tmp_NONLOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]]),
                       lower = quantile(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]], probs= c(0.025),na.rm=TRUE),
                       upper = quantile(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]], probs= c(0.975),na.rm=TRUE),
                       scale = Period[k])
  Summary_Stats_PHASE_PRECIP_NONLOG[[k]]=tmp_NONLOG
  
  tmp_LOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]]),
                              lower = quantile(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]], probs= c(0.025),na.rm=TRUE),
                              upper = quantile(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]], probs= c(0.975),na.rm=TRUE),
                              scale = Period[k])
  Summary_Stats_PHASE_PRECIP_LOG[[k]]=tmp_LOG
  
  
  #Make the p-value 
  tmp_NONLOG$p_value <- 0
  tmp_LOG$p_value <- 0
  
  #Count 
  count1_NONLOG =table( PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]] >=  TEST_RESULT_PRECIP_2_NONLOG$cor[[k]])
  count2_NONLOG =table( PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]] <=  TEST_RESULT_PRECIP_2_NONLOG$cor[[k]])
  
  #Count 
  count1_LOG =table( PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]] >=  TEST_RESULT_PRECIP_2_LOG$cor[[k]])
  count2_LOG =table( PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]] <=  TEST_RESULT_PRECIP_2_LOG$cor[[k]])
  
  
  p_1_NONLOG <- try(c(2*(count1_NONLOG[[2]]/1000)))
  p_2_NONLOG <- try(c(2 *(count2_NONLOG[[2]]/1000)))
  
  p_1_LOG <- try(c(2*(count1_LOG[[2]]/1000)))
  p_2_LOG <- try(c(2 *(count2_LOG[[2]]/1000)))
  
  p_1_NONLOG=ifelse(is.character(p_1_NONLOG)==TRUE,0,p_1_NONLOG)
  p_2_NONLOG=ifelse(is.character(p_2_NONLOG)==TRUE,0,p_2_NONLOG)
  
  p_1_LOG=ifelse(is.character(p_1_NONLOG)==TRUE,0,p_1_LOG)
  p_2_LOG=ifelse(is.character(p_2_NONLOG)==TRUE,0,p_2_LOG)
  
  tmp_NONLOG$p_value= min(p_1_NONLOG,p_2_NONLOG)
  tmp_LOG$p_value= min(p_1_LOG,p_2_LOG)
  
  Summary_Stats_PHASE_PRECIP_NONLOG[[k]]=tmp_NONLOG
  Summary_Stats_PHASE_PRECIP_LOG[[k]]=tmp_LOG
  
  
  
  
  
}

Summary_Stats_PHASE_PRECIP_2_NONLOG <- do.call(rbind, Summary_Stats_PHASE_PRECIP_NONLOG)
Summary_Stats_PHASE_PRECIP_2_LOG <- do.call(rbind, Summary_Stats_PHASE_PRECIP_LOG)

Summary_Stats_PHASE_PRECIP_2_NONLOG$actual <- TEST_RESULT_PRECIP_2_NONLOG$cor
Summary_Stats_PHASE_PRECIP_2_LOG$actual <- TEST_RESULT_PRECIP_2_LOG$cor



Summary_Stats_PHASE_PRECIP_2_NONLOG$sig <- ifelse(Summary_Stats_PHASE_PRECIP_2_NONLOG$actual>= 
                                             Summary_Stats_PHASE_PRECIP_2_NONLOG$lower &
                                             Summary_Stats_PHASE_PRECIP_2_NONLOG$actual<=
                                             Summary_Stats_PHASE_PRECIP_2_NONLOG$upper,
                                           'insig','sig')
Summary_Stats_PHASE_PRECIP_2_LOG$sig <- ifelse(Summary_Stats_PHASE_PRECIP_2_LOG$actual>= 
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$lower &
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$actual<=
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$upper,
                                                  'insig','sig')

Summary_Stats_PHASE_PRECIP_2_NONLOG$p_value_adjusted <-  p.adjust( Summary_Stats_PHASE_PRECIP_2_NONLOG $p_value, 
                                                            method ="BY")
Summary_Stats_PHASE_PRECIP_2_LOG$p_value_adjusted <-  p.adjust( Summary_Stats_PHASE_PRECIP_2_LOG $p_value, 
                                                         method ="BY")
Summary_Stats_PHASE_PRECIP_2_NONLOG$psig <- ifelse(Summary_Stats_PHASE_PRECIP_2_NONLOG$p_value_adjusted  <= 0.05,
                                            1,0.5)   
Summary_Stats_PHASE_PRECIP_2_LOG$psig <- ifelse(Summary_Stats_PHASE_PRECIP_2_LOG$p_value_adjusted  <= 0.05,
                                         1,0.5) 

Summary_Stats_PHASE_PRECIP_2_NONLOG$yscale <- as.numeric(round(Summary_Stats_PHASE_PRECIP_2_NONLOG$scale/12,3))
Summary_Stats_PHASE_PRECIP_2_LOG$yscale <- as.numeric(round(Summary_Stats_PHASE_PRECIP_2_LOG$scale/12,3))


precip_NONLOG <- ggplot(na.omit(Summary_Stats_PHASE_PRECIP_2_NONLOG), 
                            aes(x = (yscale), y = mean))+
  geom_point()+
  geom_segment(aes(x=(yscale), y=lower, xend=yscale,
                   yend=upper),alpha =0.3,size = 0.5)+
  geom_point(data=na.omit(Summary_Stats_PHASE_PRECIP_2_NONLOG),
             aes(x =yscale, y= actual,
                 fill=as.factor(sig),
                 color = as.factor(psig),
                 size = as.factor(sig),
                 shape = as.factor(psig)))+
  scale_fill_manual(values = c('grey','#3ffc97'))+
  scale_shape_manual(values=c(21,24))+
  scale_size_manual(values=c(2,3))+
  scale_color_manual(values = c('lightgrey','black'))+
  xlab("Scale (Years)") +   
  ylim(-1,1)+
  ylab("WMR")+ggtitle("Precipitation")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))


precip_LOG <- ggplot(na.omit(Summary_Stats_PHASE_PRECIP_2_LOG), 
                        aes(x = (yscale), y = mean))+
  geom_point()+
  geom_segment(aes(x=(yscale), y=lower, xend=yscale,
                   yend=upper),alpha =0.3,size = 0.5)+
  geom_point(data=na.omit(Summary_Stats_PHASE_PRECIP_2_LOG),
             aes(x =yscale, y= actual,
                 fill=as.factor(sig),
                 color = as.factor(psig),
                 size = as.factor(sig),
                 shape = as.factor(psig)))+
  scale_fill_manual(values = c('grey','#3ffc97'))+
  scale_shape_manual(values=c(21,24))+
  scale_size_manual(values=c(2,3))+
  scale_color_manual(values = c('lightgrey','black'))+

  xlab("Scale (Years)") +   
  ylim(-1,1)+
  ylab("WMR")+ggtitle("Precipitation")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))





####ENVIRONMENTAL DATA BY ITSELF ONLY


TMIN_CC1 <- read.table(here("CC1","Data","ECRtmin_PE_CC1.txt"),header=TRUE)
PRECIP_CC1 <- read.table(here("CC1","Data","CHIRP_PE_CC1.txt"),header=TRUE)

TMIN_CC1_F <- NULL
for (i in seq(1,15)){
  df <- data.frame(tmin = t(TMIN_CC1 [i,])[-1],
                   year = TMIN_CC1[i,][1],
                   month = as.numeric(seq(1,12)))
  
  TMIN_CC1_F[[i]] <- df  
  
}
TMIN_CC1_F<-do.call(rbind, TMIN_CC1_F)
TMIN_CC1_F <- na.omit(TMIN_CC1_F)

#######################################################

###PRECIP WIDE TO LONG
PRECIP_CC1_F <- NULL
for (i in seq(1,34)){
  df <- data.frame(pre = t(PRECIP_CC1[i,])[-1], year = PRECIP_CC1[i,][1],month = seq(1,12))
  
  PRECIP_CC1_F[[i]] <- df  
  
}
PRECIP_CC1_F<-do.call(rbind, PRECIP_CC1_F)
PRECIP_CC1_F <- na.omit(PRECIP_CC1_F)


############################################################
###COMBINED ENVIRONMENTAL VARIABLES 
tminz <- wt(TMIN_CC1_F$tmin , upperPeriod = 50)
tmin_gg_cc1<- as.ggplot(~plot((tminz $Scale/12), log(tminz $Power.avg+1),  type = 'l',ylim=c(0,0.8)))
tmin_gg_cc1
precipz <- wt(PRECIP_CC1_F$pre , upperPeriod = 50)
precip_gg_cc1 <- as.ggplot(~plot((precipz$Scale/12), log(precipz $Power.avg+1), type = 'l',ylim=c(0,0.8)))

tmin_gg + tmin_NONLOG
tmin_gg + tmin_LOG
precip_gg + precip_NONLOG
precip_gg + precip_LOG
