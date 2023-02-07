####################################################
###Environmental analysis of the Yasuni          ###
###- looking at the minimum temperature          ###
###and precipitation- this was last updated      ###
###in 5/5/2022 and now includes the log/unlogged ###
####################################################

library(here)
library(reshape2)
library(ggplot2)
library(WaveletComp)
library(mvcwt)
library(zoo)
library(ggplotify)


###Minimum temperature
TMIN_YAS <- read.table(here("YASUNI","Data","Environment","ECRtmin_EC_YAS.txt"),header=TRUE)
###Precipitation temperature
PRECIP_YAS <- read.table(here("YASUNI","Data","Environment","CHIRP_EC_YAS.txt"))


###TMIN WIDE TO LONG FORMAT
TMIN_YAS_F <- NULL
for (i in seq(1,15)){
  df <- data.frame(tmin = t(TMIN_YAS[i,])[-1], year = TMIN_YAS[i,][1],month = as.numeric(seq(1,12)))
  
TMIN_YAS_F[[i]] <- df  
  
}


TMIN_YAS_F<-do.call(rbind, TMIN_YAS_F)
TMIN_YAS_F <- na.omit(TMIN_YAS_F)

###PRECIP WIDE TO LONG FORMAT
PRECIP_YAS_F <- NULL
for (i in seq(1,34)){
  df <- data.frame(pre = t(PRECIP_YAS[i,])[-1], year = PRECIP_YAS[i,][1],month = seq(1,12))

  PRECIP_YAS_F[[i]] <- df  
  
}

PRECIP_YAS_F<-do.call(rbind, PRECIP_YAS_F)
PRECIP_YAS_F <- na.omit(PRECIP_YAS_F)

###Combining means that we are subsetting
PRECIP_YAS_F <- subset(PRECIP_YAS_F, PRECIP_YAS_F$V1 > 1999)
PRECIP_YAS_F <- PRECIP_YAS_F[-c(171:174),]

###COMBINED ENVIRONMENTAL VARIABLES

ENVIR_YAS <- cbind.data.frame(TMIN_YAS_F, precip = PRECIP_YAS_F$pre, id = 'YAS')
ENVIR_YAS <- ENVIR_YAS[-1,]

ENVIR_YAS$jul <- seq(1, nrow(ENVIR_YAS))


###############################################################
################# NEED THE CODE DATA FILE ######################
###############################################################

CODE <- read.csv(here("YASUNI", "Data", "CODE_YAS_Lasky_ncg.csv"))

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

# THESE ARE THE CENSORED

CENSORED_SP <- c("Moraceae","Solanaceae","Clusiaceae",
                 "Clusiaceae?","Asteraceae","Araceae")

###This gives me the code for the species that should be excluded 
CENSORED_SP_DF <- CODE %>%
  filter(Family.APG.IV. %in% CENSORED_SP) %>%
  select(CODIGO)

##############################################################################
### The Original Yasuni Time Series ##########################################
##############################################################################
YASUNI <- read.csv(here("YASUNI", "Data", "YasuniTS.csv"))
colnames(YASUNI)[1] <- "Date"

# Make into dates
YASUNI$Date <- as.Date(YASUNI$Date, format = "%m/%d/%Y")

### Seems like each date is unique and in order!
### IT SEEMS THAT 2000-02-28 is the beginning (59) is 
###the Julian day

YASUNI <- YASUNI %>% 
  mutate(Day_Length = c(59, diff(Date)),
         Day_Year = cumsum(Day_Length),
         Day_Year_Y = (Day_Year) / 365.25)

###Should have 1741 species/columns 
### Gets rid of any species where there is no entries - should be 1194 species
YASUNI <- YASUNI[, colSums(YASUNI != 0) != 0]


###This will melt the YASUNI dataframe 
###EXCLUDING the Date (Because we're using)
###Day_Year_Y as the variable to calculate the scales.
###This will create a dataframe, with the Day_Year_Y, with the variable name
###(The Species name) and the value.

YASUNI_MELT <- melt(YASUNI[, -which(names(YASUNI) ==
                                      "Day_Length" |
                                      names(YASUNI) == "Day_Year")],
                    id.vars = c("Day_Year_Y", "Date")
)
### This removes all species that have the code that pertains to the 
###censored family 
YASUNI_MELT <- YASUNI_MELT%>%
  filter(!(variable %in% CENSORED_SP_DF$CODIGO)) 

###This more convoluted way was because previous analyses tried
###incorporate censored species unlike our current method where we just
###remove them entirely

###1072 species EXCLUDING censored species 

##############################################################
##Also there are some removed species that Nancy provided#####
##############################################################
YAS_REMOVED_NANCY <-read.csv(here("YASUNI","Data","11_16_NANCY_REMOVED.csv"))

YASUNI_F <- YASUNI_MELT %>%  
  filter(!(variable %in% YAS_REMOVED_NANCY $CODIGO))

length(unique(YASUNI_F $variable))

YASUNI_F <- YASUNI_F %>%
  select(Date, variable,value)

#Should equal to 1059
YASUNI_F <- dcast(YASUNI_F, Date~variable)
YASUNI_F$Month <- as.numeric(format(YASUNI_F$Date,"%m"))
YASUNI_F$Year<- as.numeric(format(YASUNI_F$Date,"%Y"))

#####################################
###We're going to aggregate by month and year and we're going to average
###the amount of seeds
YASUNI_AGG_M_Y<- aggregate(YASUNI_F[,2:1055], by=list(YASUNI_F$Month, YASUNI_F$Year), 'mean',na.action=na.pass())

######################################################
###I checked and the months are in sequential order
###We're missing the 5 from 2008 and the 6 from 2009
#####################################################
First_Fifth <- YASUNI_AGG_M_Y[1:99,]
Second_Fifth <- YASUNI_AGG_M_Y[100:110,]
Third_Fifth <- YASUNI_AGG_M_Y[111:203,]

First_Part <- rbind(First_Fifth, NA)
Second_Part <- rbind(First_Part,Second_Fifth)
Third_Part <- rbind(Second_Part, NA)
Fourth_Part<- rbind(Third_Part, Third_Fifth)
Fourth_Part[,1]<- seq(1,205)
Fourth_Part[,3:ncol(Fourth_Part)] <- na.approx(Fourth_Part[,3:ncol(Fourth_Part)])


YASUNI_MELT_M <- melt(Fourth_Part[,-2], id.vars = c('Group.1'))

#####################################################
YAS_ALL_MONTH <-dcast(YASUNI_MELT_M,
                   Group.1~variable,value.var='value',
                   direction="wide")

YAS_ALL_MONTH_2<- YAS_ALL_MONTH[,colSums(is.na(YAS_ALL_MONTH)) == 0]
YAS_ALL_MONTH_F <- subset(YAS_ALL_MONTH_2,YAS_ALL_MONTH_2$Group.1 < 170)

x=YAS_ALL_MONTH_F$Group.1
y_nonlog= YAS_ALL_MONTH_F[,2:ncol(YAS_ALL_MONTH_F)]
y_log= log(YAS_ALL_MONTH_F[,2:ncol(YAS_ALL_MONTH_F)]+1)

YAS_WMR_NONLOG= wmr(mvcwt(x,y_nonlog, scales = seq(2,70)))
YAS_WMR_LOG= wmr(mvcwt(x,y_log, scales = seq(2,70)))

image(YAS_WMR_NONLOG)


M_YAS_COI <-  matrix(NA, ncol = ncol(YAS_WMR_LOG$z[, , 1]), nrow = nrow(YAS_WMR_LOG$z[, , 1]))

for (i in 1:ncol(M_YAS_COI )) M_YAS_COI [which(YAS_WMR_LOG$x > (min(YAS_WMR_LOG$x) + 1 * YAS_WMR_LOG$y[i]) & YAS_WMR_LOG$x < (max(YAS_WMR_LOG$x) -1* YAS_WMR_LOG$y[i])), i] <- 1
M_YAS_COI0 <- M_YAS_COI 

saveRDS(M_YAS_COI0 , file ='M_YAS_COI0.RDS')


###########################
###TMIN Wavelet analysis###
###########################
TMIN = data.frame(tmin =ENVIR_YAS$tmin)

TEMP_WAV <- analyze.wavelet(TMIN,'tmin',
                            loess.span = 0,dt = 1, lowerPeriod =2, upperPeriod = 70)
wt.image(TEMP_WAV)

Period = seq(2,70)

##########################################
#Reconstructing the temperature data at ##
#specific periods                       ##
##########################################
TEST_RESULTS_TMIN_LOG=NULL
TEST_RESULTS_TMIN_NONLOG=NULL

YAS_WMR_NONLOG$z[,,1] <- YAS_WMR_NONLOG$z[,,1]* M_YAS_COI0
YAS_WMR_LOG$z[,,1] <- YAS_WMR_LOG$z[,,1]* M_YAS_COI0

for (i in seq(1,69)){
RECONSTRUCTED =reconstruct(TEMP_WAV,
                           sel.period = Period[i],
                           only.sig=FALSE)

recon_tmin =RECONSTRUCTED$series$tmin.r

YAS_WMR_SCALE_NONLOG=YAS_WMR_NONLOG$z[1:169,i,1]
YAS_WMR_SCALE_LOG=YAS_WMR_LOG$z[1:169,i,1]

tempo_NONLOG = cor.test(recon_tmin, YAS_WMR_SCALE_NONLOG,method='pearson',na.rm=TRUE)
tempo_LOG = cor.test(recon_tmin, YAS_WMR_SCALE_LOG,method='pearson',na.rm=TRUE)

df_temp_NONLOG <- cbind.data.frame(Scale = Period[i],t= tempo_NONLOG$estimate) 
df_temp_LOG <- cbind.data.frame(Scale = Period[i],t= tempo_LOG$estimate) 

TEST_RESULTS_TMIN_NONLOG[[i]] = df_temp_NONLOG
TEST_RESULTS_TMIN_LOG[[i]] = df_temp_LOG
}


########
########
TEST_RESULT_TMIN_2_NONLOG <- do.call(rbind, TEST_RESULTS_TMIN_NONLOG)
TEST_RESULT_TMIN_2_LOG <- do.call(rbind, TEST_RESULTS_TMIN_LOG)

PHASE_RANDOMIZED_ALL_LOG = NULL
PHASE_RANDOMIZED_ALL_NONLOG = NULL

for (j in seq(1,69)){
#Pull out a specific scale 
RECONSTRUCTED =reconstruct(TEMP_WAV,sel.period = Period[j],
                           only.sig=FALSE)
 #Pulls out the reconstructed time series 
recon_tmin =data.frame(jul=seq(1,169),
                         tmin=
                         RECONSTRUCTED$series$tmin.r)

#Sample where the time series would start (1 to 169)
random_starter = sample(1:169,1000,replace =T)

YAS_WMR_SCALE_NONLOG=YAS_WMR_NONLOG$z[1:169,j,1]
YAS_WMR_SCALE_LOG=YAS_WMR_LOG$z[1:169,j,1]


cor_1000_NONLOG = NULL
cor_1000_LOG = NULL

for (i in seq(random_starter)){
  
  start = recon_tmin[recon_tmin$jul>=random_starter[i],]
  end = recon_tmin[recon_tmin$jul < random_starter[i],]
  connected = rbind(start,end)
  
  
 tmp_NONLOG<- cor.test(connected$tmin, 
                YAS_WMR_SCALE_NONLOG,method='kendall',na.rm=TRUE)$estimate
 
 
 tmp_LOG<- cor.test(connected$tmin, 
                       YAS_WMR_SCALE_LOG,method='kendall',na.rm=TRUE)$estimate
 
 cor_1000_NONLOG[[i]]= tmp_NONLOG
 cor_1000_LOG[[i]]= tmp_LOG
 
}

PHASE_RANDOMIZED_ALL_NONLOG[[j]]=unlist(cor_1000_NONLOG)
PHASE_RANDOMIZED_ALL_LOG[[j]]=unlist(cor_1000_LOG)

}



Summary_Stats_PHASE_NONLOG= NULL
for (k in seq(1,69)){
  tmp_NONLOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_NONLOG[[k]]),
  lower_nonlog = quantile(PHASE_RANDOMIZED_ALL_NONLOG[[k]], probs= c(0.025),na.rm=TRUE),
  upper_nonlog = quantile(PHASE_RANDOMIZED_ALL_NONLOG[[k]], probs= c(0.975),na.rm=TRUE),
  scale = Period[k])

  
  #Make the p-value 
  tmp_NONLOG$p_value <- 0
  
   #Count 
  count1_NONLOG =table( PHASE_RANDOMIZED_ALL_NONLOG[[k]] >=  TEST_RESULT_TMIN_2_NONLOG$t[[k]])
  count2_NONLOG =table(  PHASE_RANDOMIZED_ALL_NONLOG[[k]] <=  TEST_RESULT_TMIN_2_NONLOG$t[[k]])
    
    
    p_1_NONLOG <- try(c(2*(count1_NONLOG[[2]]/1000)))
    p_2_NONLOG <- try(c(2 *(count2_NONLOG[[2]]/1000)))
    
    p_1_NONLOG=ifelse(is.character(p_1_NONLOG)==TRUE,0,p_1_NONLOG)
    p_2_NONLOG=ifelse(is.character(p_2_NONLOG)==TRUE,0,p_2_NONLOG)
    
    
    tmp_NONLOG$p_value= min(p_1_NONLOG,p_2_NONLOG)
  
  


 
  Summary_Stats_PHASE_NONLOG[[k]]=tmp_NONLOG
}  


Summary_Stats_PHASE_LOG= NULL
for (k in seq(1,69)){
  tmp_LOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_LOG[[k]]),
                              lower_log = quantile(PHASE_RANDOMIZED_ALL_LOG[[k]], probs= c(0.025),na.rm=TRUE),
                              upper_log = quantile(PHASE_RANDOMIZED_ALL_LOG[[k]], probs= c(0.975),na.rm=TRUE),
                              scale = Period[k])
  
  
  #Make the p-value 
  tmp_LOG$p_value <- 0
  
  #Count 
  count1_LOG =table( PHASE_RANDOMIZED_ALL_LOG[[k]] >=  TEST_RESULT_TMIN_2_LOG$t[[k]])
  count2_LOG =table(  PHASE_RANDOMIZED_ALL_LOG[[k]] <=  TEST_RESULT_TMIN_2_LOG$t[[k]])
  
  
  p_1_LOG <- try(c(2*(count1_LOG[[2]]/1000)))
  p_2_LOG <- try(c(2 *(count2_LOG[[2]]/1000)))
  
  p_1_LOG=ifelse(is.character(p_1_LOG)==TRUE,0,p_1_LOG)
  p_2_LOG=ifelse(is.character(p_2_LOG)==TRUE,0,p_2_LOG)
  
  
  tmp_LOG$p_value= min(p_1_LOG,p_2_LOG)

  
  Summary_Stats_PHASE_LOG[[k]]=tmp_LOG
}  

Summary_Stats_PHASE_2_NONLOG <- do.call(rbind, Summary_Stats_PHASE_NONLOG)
Summary_Stats_PHASE_2_LOG <- do.call(rbind, Summary_Stats_PHASE_LOG)

Summary_Stats_PHASE_2_NONLOG$p.adjust <- p.adjust(Summary_Stats_PHASE_2_NONLOG$p_value,
                                           method = 'BY')
Summary_Stats_PHASE_2_LOG$p.adjust <- p.adjust(Summary_Stats_PHASE_2_LOG$p_value,
                                                  method = 'BY')
Summary_Stats_PHASE_2_NONLOG$psig <- ifelse(Summary_Stats_PHASE_2_NONLOG$p.adjust  <= 0.05,
                                     1,0.5)
Summary_Stats_PHASE_2_LOG$psig <- ifelse(Summary_Stats_PHASE_2_LOG$p.adjust  <= 0.05,
                                            1,0.5)

Summary_Stats_PHASE_2_NONLOG$actual <- TEST_RESULT_TMIN_2_NONLOG$t
Summary_Stats_PHASE_2_LOG$actual <- TEST_RESULT_TMIN_2_LOG$t

Summary_Stats_PHASE_2_NONLOG$sig <- ifelse( Summary_Stats_PHASE_2_NONLOG$actual>= 
          Summary_Stats_PHASE_2_NONLOG$lower &
          Summary_Stats_PHASE_2_NONLOG$actual<=Summary_Stats_PHASE_2_NONLOG$upper,
        'insig','sig')

Summary_Stats_PHASE_2_LOG$sig <- ifelse( Summary_Stats_PHASE_2_LOG$actual>= 
                                       Summary_Stats_PHASE_2_LOG$lower &
                                       Summary_Stats_PHASE_2_LOG$actual<=Summary_Stats_PHASE_2_LOG$upper,
                                     'insig','sig')


Summary_Stats_PHASE_2_NONLOG$yscale <- as.numeric(round(Summary_Stats_PHASE_2_NONLOG$scale/12,3))
Summary_Stats_PHASE_2_LOG$yscale <- as.numeric(round(Summary_Stats_PHASE_2_LOG$scale/12,3))


###The Minimum temperature Plots###
tmin_yas_NONLOG<- 
  ggplot(na.omit(Summary_Stats_PHASE_2_NONLOG), 
         aes(x = (yscale), y = mean))+
  geom_point()+geom_segment(aes(x=(yscale), y=lower_nonlog, xend=yscale,
                              yend=upper_nonlog),alpha =0.3,size = 0.5)+
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

tmin_yas_LOG<- 
  ggplot(na.omit(Summary_Stats_PHASE_2_LOG), 
         aes(x = (yscale), y = mean))+
  geom_point()+geom_segment(aes(x=(yscale), y=lower_log, xend=yscale,
                                yend=upper_log),alpha =0.3,size = 0.5)+
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
  ylim(-1,1) +
  ylab("WMR")+ggtitle("Min. Temperature")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))
##############################################################
##############################################################

##############################
#########PRECIPITATION#########
##############################

PRECIP = data.frame(precip =ENVIR_YAS$precip)
PRECIP_WAV <- analyze.wavelet(PRECIP,'precip',
                            loess.span = 0,dt = 1,
                            lowerPeriod =2, upperPeriod = 70)
wt.image(PRECIP_WAV)




TEST_RESULTS_PRECIP_NONLOG=NULL
TEST_RESULTS_PRECIP_LOG=NULL

for (i in seq(1,69)){
  RECONSTRUCTED = reconstruct(PRECIP_WAV,sel.period = Period[i],
                             only.sig=FALSE)
  recon_precip = RECONSTRUCTED$series$precip.r
  
  YAS_WMR_SCALE_LOG = YAS_WMR_LOG$z[1:169,i,1]
  YAS_WMR_SCALE_NONLOG = YAS_WMR_NONLOG$z[1:169,i,1]
  
  tempo_NONLOG = cor.test(recon_precip, YAS_WMR_SCALE_NONLOG,method='pearson')$estimate
  tempo_LOG = cor.test(recon_precip, YAS_WMR_SCALE_LOG,method='pearson')$estimate
  
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

for (j in seq(1,69)){
  #Pull out a specific scale 
  RECONSTRUCTED =reconstruct(PRECIP_WAV,sel.period = Period[j],
                             only.sig=FALSE)
  #Pulls out the reconstructed time series 
  recon_precip =data.frame(jul=seq(1,169),
                         precip=RECONSTRUCTED$series$precip.r)
  #Sample where the time series would start (1 to 169)
  random_starter = sample(1:169,1000,replace =T)
  
  
  YAS_WMR_SCALE_NONLOG=YAS_WMR_NONLOG$z[1:169,j,1]
  YAS_WMR_SCALE_LOG=YAS_WMR_LOG$z[1:169,j,1]
  
  cor_1000_NONLOG = NULL
  cor_1000_LOG = NULL
  
  for (i in seq(random_starter)){
    
    start = recon_precip[recon_precip$jul>=random_starter[i],]
    end = recon_precip[recon_precip$jul < random_starter[i],]
    connected = rbind(start,end)
    
    
    tmp_NONLOG<- cor.test(connected$precip, 
              YAS_WMR_SCALE_NONLOG,method='pearson')$estimate
    tmp_LOG<- cor.test(connected$precip, 
              YAS_WMR_SCALE_LOG,method='pearson')$estimate
    

cor_1000_NONLOG[[i]]= tmp_NONLOG
cor_1000_LOG[[i]]= tmp_LOG

  }
  PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[j]]=unlist(cor_1000_NONLOG)
  PHASE_RANDOMIZED_ALL_PRECIP_LOG[[j]]=unlist(cor_1000_LOG)
  
}

Summary_Stats_PHASE_PRECIP_NONLOG= NULL
Summary_Stats_PHASE_PRECIP_LOG= NULL

for (k in seq(1,69)){
  tmp_NONLOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]]),
                       lower = quantile(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]], probs= c(0.025),na.rm=TRUE),
                       upper = quantile(PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]], probs= c(0.975),na.rm=TRUE),
                       scale = Period[k])
  tmp_LOG=cbind.data.frame(mean = mean(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]]),
                              lower = quantile(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]], probs= c(0.025),na.rm=TRUE),
                              upper = quantile(PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]], probs= c(0.975),na.rm=TRUE),
                              scale = Period[k])
  #Make the p-value 
  tmp_NONLOG$p_value <- 0
  tmp_LOG$p_value <- 0
  
  #Count 
  count1_NONLOG =table(  PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]] >=  TEST_RESULT_PRECIP_2_NONLOG$cor[[k]])
  count2_NONLOG =table(   PHASE_RANDOMIZED_ALL_PRECIP_NONLOG[[k]] <=  TEST_RESULT_PRECIP_2_NONLOG$cor[[k]])
  
  count1_LOG =table(  PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]] >=  TEST_RESULT_PRECIP_2_LOG$cor[[k]])
  count2_LOG =table(   PHASE_RANDOMIZED_ALL_PRECIP_LOG[[k]] <=  TEST_RESULT_PRECIP_2_LOG$cor[[k]])
  
  p_1_NONLOG<- try(c(2*(count1_NONLOG[[2]]/1000)))
  p_2_NONLOG <- try(c(2 *(count2_NONLOG[[2]]/1000)))
  
  p_1_NONLOG=ifelse(is.character(p_1_NONLOG)==TRUE,0,p_1_NONLOG)
  p_2_NONLOG=ifelse(is.character(p_2_NONLOG)==TRUE,0,p_2_NONLOG)
  
  p_1_LOG<- try(c(2*(count1_LOG[[2]]/1000)))
  p_2_LOG <- try(c(2 *(count2_LOG[[2]]/1000)))
  
  p_1_LOG=ifelse(is.character(p_1_LOG)==TRUE,0,p_1_LOG)
  p_2_LOG=ifelse(is.character(p_2_LOG)==TRUE,0,p_2_LOG)
  
  
  tmp_NONLOG$p_value= min(p_1_NONLOG,p_2_NONLOG)
  tmp_LOG$p_value= min(p_1_LOG,p_2_LOG)
  
  
  Summary_Stats_PHASE_PRECIP_NONLOG[[k]]=tmp_NONLOG
  Summary_Stats_PHASE_PRECIP_LOG[[k]]=tmp_LOG
  
  }

Summary_Stats_PHASE_PRECIP_2_NONLOG<- do.call(rbind, Summary_Stats_PHASE_PRECIP_NONLOG)
Summary_Stats_PHASE_PRECIP_2_NONLOG$actual <- TEST_RESULT_PRECIP_2_NONLOG$cor

Summary_Stats_PHASE_PRECIP_2_LOG<- do.call(rbind, Summary_Stats_PHASE_PRECIP_LOG)
Summary_Stats_PHASE_PRECIP_2_LOG$actual <- TEST_RESULT_PRECIP_2_LOG$cor



Summary_Stats_PHASE_PRECIP_2_NONLOG$p.adjust <- p.adjust(Summary_Stats_PHASE_PRECIP_2_NONLOG$p_value,
                                           method = 'BY')

Summary_Stats_PHASE_PRECIP_2_LOG$p.adjust <- p.adjust(Summary_Stats_PHASE_PRECIP_2_LOG$p_value,
                                                         method = 'BY')



Summary_Stats_PHASE_PRECIP_2_NONLOG$psig <- ifelse(Summary_Stats_PHASE_PRECIP_2_NONLOG$p.adjust  <= 0.05,
                                     1,0.5)
Summary_Stats_PHASE_PRECIP_2_LOG$psig <- ifelse(Summary_Stats_PHASE_PRECIP_2_LOG$p.adjust  <= 0.05,
                                                   1,0.5)
Summary_Stats_PHASE_PRECIP_2_NONLOG$sig <- 
                                ifelse(Summary_Stats_PHASE_PRECIP_2_NONLOG$actual>= 
                                      Summary_Stats_PHASE_PRECIP_2_NONLOG$lower &
                                      Summary_Stats_PHASE_PRECIP_2_NONLOG$actual<=
                                        Summary_Stats_PHASE_PRECIP_2_NONLOG$upper,
                                     'insig','sig')
Summary_Stats_PHASE_PRECIP_2_LOG$sig <- ifelse(Summary_Stats_PHASE_PRECIP_2_LOG$actual>= 
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$lower &
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$actual<=
                                                    Summary_Stats_PHASE_PRECIP_2_LOG$upper,
                                                  'insig','sig')
Summary_Stats_PHASE_PRECIP_2_NONLOG$yscale <- as.numeric(round(Summary_Stats_PHASE_PRECIP_2_NONLOG$scale/12,3))
Summary_Stats_PHASE_PRECIP_2_LOG$yscale <- as.numeric(round(Summary_Stats_PHASE_PRECIP_2_LOG$scale/12,3))



precip_yas_NONLOG <- ggplot(na.omit(Summary_Stats_PHASE_PRECIP_2_NONLOG), 
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
  ylim(-1.0,1.0)+
  ylab("WMR")+ggtitle("Precipitation")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))



precip_yas_LOG <- ggplot(na.omit(Summary_Stats_PHASE_PRECIP_2_LOG), 
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
  ylim(-1.0,1.0)+
  ylab("WMR")+ggtitle("Precipitation")+
  theme_bw()+theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size=13,color='black'))

##############################################################
###DOING ALL THE YEARS ENVIRONMENT FOR THE WAVELET ANALYSIS### 
##############################################################


###Minimum temperature
TMIN_YAS <- read.table(here("YASUNI","Data","Environment","ECRtmin_EC_YAS.txt"),header=TRUE)
###Precipitation temperature
PRECIP_YAS <- read.table(here("YASUNI","Data","Environment","CHIRP_EC_YAS.txt"))




###TMIN WIDE TO LONG
TMIN_YAS_F <- NULL
for (i in seq(1,15)){
  df <- data.frame(tmin = t(TMIN_YAS[i,])[-1], year = TMIN_YAS[i,][1],month = as.numeric(seq(1,12)))
  
  TMIN_YAS_F[[i]] <- df  
  
}


TMIN_YAS_F<-do.call(rbind, TMIN_YAS_F)
TMIN_YAS_F <- na.omit(TMIN_YAS_F)


tminz <- wt(TMIN_YAS_F$tmin, upperPeriod = 70)
tmin_GG_yas<- as.ggplot(~plot((tminz $Scale/12),log(tminz$Power.avg + 1), ylim=c(0,0.8),
     main='Minimum temperature', type='l',
     xlab= 'Scale (Years)',
     ylab ='Average power'))




###PRECIP WIDE TO LONG
PRECIP_YAS_F <- NULL
for (i in seq(1,34)){
  df <- data.frame(pre = t(PRECIP_YAS[i,])[-1], year = PRECIP_YAS[i,][1],month = seq(1,12))
  
  PRECIP_YAS_F[[i]] <- df  
  
}
PRECIP_YAS_F<-do.call(rbind, PRECIP_YAS_F)
PRECIP_YAS_F <- na.omit(PRECIP_YAS_F)


precipz<- wt(PRECIP_YAS_F $pre, upperPeriod = 70)
precip_GG_yas <- as.ggplot(~plot((precipz $Scale/12),log(precipz $Power.avg+1), type = 'l',ylim=c(0,0.8)))


###non-transformed minimum temperature
(tmin_GG_yas+tmin_yas_NONLOG)/
 (tmin_gg_cc1+ tmin_NONLOG)

(tmin_GG_yas+tmin_yas_LOG)/
  (tmin_gg_cc1 + tmin_LOG)

(precip_GG_yas+precip_yas_NONLOG)/
  (precip_gg_cc1 + precip_NONLOG)

(precip_GG_yas+precip_yas_LOG)/
  (precip_gg_cc1 + precip_LOG)

