####################################################
### Yasuni- ALL COMMUNITY ANALYSIS (Transformed) ###
####################################################
library(here)
library(reshape2)
library(mvcwt) #not on CRAN, get the updated version
library(data.table)
library(dplyr)
library(biwavelet)
library(stringr)
library(ggplot2)
library(patchwork)

##############################################
# This function is for the total              #
# community analysis: Use YASUNI_ALL_FINAL.csv#
##############################################
### The file needed [-1] because the rows are included-
###this is the data file that is cleaned up and used for all
###analyses
YASUNI_F <- fread(here("YASUNI", "Data", 
                       "YASUNI_ALL_FINAL.csv"))[, -1]
source(here("SUPP", "image_2.R"))
source(here("SUPP", "modified_boot_wmr.R"))

################
### reshape it##
################
YASUNI_WIDE_F_ALL <- reshape2::dcast(YASUNI_F,
  Day_Year_Y ~ variable,
  value.var = "value",
  direction = "wide"
)

#########################################################
#########THIS IS THE MAIN WMR PLOT ######################
#########################################################

### The x (time)
x <- YASUNI_WIDE_F_ALL$Day_Year_Y
### The y (abundance that is log transformed)
y <- log(YASUNI_WIDE_F_ALL[, 2:ncol(YASUNI_WIDE_F_ALL)] + 1)

### Run the multivariate continuous wavelet transformation-
### Let the package choose the minimum scale and the maximum scale
MVCWT_YAS_ALL <- mvcwt(x, y,
  min.scale = get.min.scale(x),
  max.scale = get.max.scale(x)
)
### Get the WMR
MVCWT_YAS_WMR <- wmr(MVCWT_YAS_ALL)


### This is the bootstrapping-phase randomized
MVCWT_YAS_WMR_BOOT <- wmr.boot(MVCWT_YAS_ALL, smoothing = 1, reps = 1000, mr.func = "wmr")
MVCWT_YAS_WMR_BOOT_2_TRANSFORMED <- getBootWMR(MVCWT_YAS_ALL, reps = 1000)

####################################################
### Make the figures that include both the Yasuni ###
### and CC1                                       ###
#####################################################

load(here( "YASUNI", "Output", "YASUNI_TRANSFORMED",
  "MVCWT_YAS_WMR_BOOT_TRANSFORMED.RData"
))


image.mvcwt2(MVCWT_YAS_WMR_BOOT,
  z.fun = "Re",
  bound = 1, reset.par = FALSE,
  ylim = c(0.08, 8.5)
)

contour(MVCWT_YAS_WMR, add = TRUE, ylim = c(0.08, 8.5))


# contour(MVCWT_YAS_WMR_BOOT$x,
#  MVCWT_YAS_WMR_BOOT$y, MVCWT_YAS_WMR_BOOT$z.boot[,,1],
#  levels = c(0.025, 0.975), lwd = 3, add = T,drawlabels = F)
# ontour(MVCWT_YAS_WMR,add=TRUE)

##################################################
##################################################
##################################################
### Total species count and seed count figures ####
##################################################

YAS_Count <- cbind.data.frame(
  Day_Year_Y =
    YASUNI_WIDE_F_ALL$Day_Year_Y,
  Species_Count = rowSums(YASUNI_WIDE_F_ALL[, -1] != 0),
  Seed_Count = rowSums(YASUNI_WIDE_F_ALL)
)
high_species_name <- names(sort(colSums( YASUNI_WIDE_F_ALL), decreasing=TRUE)[1:5])


### This is for the total Species Count
YAS_SC_GG <- ggplot(
  YAS_Count,
  aes(x = Day_Year_Y, y = Species_Count)
) +
  geom_line(color = "darkgreen") +
  scale_x_continuous(
    breaks = seq(5, 17.5, 5), expand = c(0, 0),
    labels = c(2005, 2010, 2015)
  ) +
  xlab("Year") +
  ylab("Total species count") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 20),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

### This is for the total seed count
YAS_PC_GG <- 
  ggplot(YAS_Count, aes(x = Day_Year_Y, y = log10(Seed_Count + 1))) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(5, 17.5, 5),
    expand = c(0, 0),
    labels = c(2005, 2010, 2015)
  ) +
  xlab("Year") +
  ylab("Total seed count (log10)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 20)
  )

###For making the large multiplot
###(YAS_SC_GG + CC1_SP_COUNT_GG_2) / (YAS_PC_GG + CC1_PHEN_SUM_GG)

##########################################
### AVERAGE WMR ACROSS TIME SQUARE PLOT ###
##########################################

###This ensures that anything out of the cone of influence is 
###not included for the anaysis.

MVCWT_YAS_WMR_COI <- matrix(NA, ncol = ncol(MVCWT_YAS_WMR$z[, , 1]), nrow = nrow(MVCWT_YAS_WMR$z[, , 1]))

for (i in 1:ncol(MVCWT_YAS_WMR_COI)) MVCWT_YAS_WMR_COI[which(MVCWT_YAS_WMR$x > (min(MVCWT_YAS_WMR$x) + 1 * MVCWT_YAS_WMR$y[i]) & MVCWT_YAS_WMR$x < (max(MVCWT_YAS_WMR$x) - 1 * MVCWT_YAS_WMR$y[i])), i] <- 1
MVCWT_YAS_WMR_COI0 <- MVCWT_YAS_WMR_COI

#saveRDS(MVCWT_YAS_WMR_COI0, file ='MVCWT_YAS_WMR_COI0.RDS')


yas_transformed_CI<-NULL
 for (k in seq(1, length(MVCWT_YAS_WMR_BOOT_2_log))){
 tmp <-  MVCWT_YAS_WMR_BOOT_2_log[[k]] 
 new_tmp <- colMeans(tmp * MVCWT_YAS_WMR_COI0[,k],na.rm=TRUE)
 
 lowerbound <-quantile(new_tmp, probs= c(0.025),na.rm=TRUE)
 upperbound <- quantile(new_tmp, probs= c( 0.975),na.rm=TRUE)

 yas_transformed_CI[[k]]<- cbind(lb = lowerbound, up = upperbound)
 }

yas_transformed_CI<- do.call(rbind,yas_transformed_CI)
 

###################
###The average WMR#
###################
YAS_AVG_WMR <- data.frame(Scale = MVCWT_YAS_WMR$y,
                          WMR = colMeans(MVCWT_YAS_WMR_COI0 * MVCWT_YAS_WMR$z[,,1],na.rm=TRUE),
                          lb = yas_transformed_CI[,1],
                          ub = yas_transformed_CI[,2])


YAS_AVG_WMR$sig <- ifelse(YAS_AVG_WMR$WMR >YAS_AVG_WMR$ub
                        |YAS_AVG_WMR$WMR <YAS_AVG_WMR$lb, 'sig','notsig')

#######################################
### Creating the confidence interval###
#######################################
YAS_GG <- 
  ggplot(YAS_AVG_WMR,
         aes(x = Scale, y = WMR)) +
  geom_line(size = 1.1,aes(color=sig)) +
  geom_ribbon(aes(ymin= lb, ymax=ub),alpha = 0.1)+
  scale_color_manual(values=c('black','#FFB61D'))+
  ylim(0, 1) +
  xlab("Period (Years)") +
  ylab("Avg. WMR")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  coord_flip()
