##################################################################
#### CODE THAT INCORPORATES THE NEW DISPERSAL GROUP FROM SIMON ###
##################################################################
#####################
###package required##
#####################


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
library(readr)
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



###This removes any species code from the FAMILY that DOES NOT
###APPEAR IN THE YASUNI DATA- this removes about 31 species-code
FAMILY_CODE    <- FAMILY_CODE  %>%
                  filter(CODIGO %in% 
                   as.character(YASUNI_F$variable))

FAMILY_LENGTH <- table(FAMILY_CODE$Family.APG.IV.)
FAMILY_5 =data.frame(t(FAMILY_LENGTH[FAMILY_LENGTH >= 5][-1]))[,-1]

colnames(FAMILY_5) <- c("Family","Species") 

NUMBER_NEEDED <- sort(unique(FAMILY_5$Species))

Yas_Dispersal_Code <- FAMILY_CODE$CODIGO

#####################
#####################
#####################
#####################
###PERMUTATED PART###
#####################

for (k in seq(1,25)){
  
tmp <- Permutation_Boot_Strap(
       actualdata=YASUNI_DT,
       Yas_Dispersal_Code, 
       disp_n= 54,#NUMBER_NEEDED[[k]], 
       log=TRUE)

assign(paste("PERM", 54, sep = ""), tmp)
saveRDS(get(paste("PERM", 54, sep = "")),
     file =here("YASUNI","Output",
     "FAMILY_PERM_LOGTRANSFORMED",
     paste("PERM", 54,".RDS", sep = "")))
rm(list=paste("PERM", NUMBER_NEEDED[[k]], sep = ""))

}
###This should be saved in a file already

################################################
################################################

#SET WORKING DIRECTORY TO 
lofs <-list.files( here("YASUNI","Output",
                "FAMILY_PERM_LOGTRANSFORMED"),pattern = ".RData")
############################################################
############################################################


###36 is gone 


###Notice that's it's not in order
NUMBER_ORDER <- parse_number(lofs)


YAS_WMR_PERM_DF <- NULL

for (r in seq(37,nrow(FAMILY_5))){
  FAMILY_NAME <-   FAMILY_5[r,1]
  FAMILY_NUMBER <- FAMILY_5[r,2]
  
  ######################
  ###The actual part###
  #####################
  
  FAM_CODE <- FAMILY_CODE[FAMILY_CODE$Family.APG.IV. == FAMILY_NAME,]$CODIGO
  
  FAM_TS <- YASUNI_DT[,..FAM_CODE]
  
  
  x = YASUNI_DT$Day_Year_Y
  ### For the y, you must rely on the species data from the 2nd column to the end of the
  ###column- LOG TRANSFORMED
  y =   log(FAM_TS+1)
  
  WMR_FAM = wmr(mvcwt(x, y, min.scale = get.min.scale(x), 
                    max.scale = get.max.scale(x)))

  
  
  
  ####################################################################
  ###After that is accomplished, we now do the wavelet modulus ration#
  ####################################################################


   index =which( NUMBER_ORDER == FAMILY_NUMBER)
  
   tmp = load(here("YASUNI","Output",
              "FAMILY_PERM_LOGTRANSFORMED",
              lofs[[index]]))
   
   PERM = get(tmp)

   
   YAS_WMR_PERM_DF[[r]]<-
     Permutation_Actual_Combiner(MVCWT_YAS_WMR_COI0,PERM,WMR_FAM )
   
   
  rm(list=str_remove(lofs[[index]], ".RData"))
}

for (f in seq(1, nrow(FAMILY_5))){
  
  YAS_WMR_PERM_DF[[f]]$fam = FAMILY_5$Family[[f]]
}

YAS_WMR_PERM_ACTUAL_DF<- do.call(rbind,YAS_WMR_PERM_DF)

saveRDS(YAS_WMR_PERM_ACTUAL_DF, file=here("YASUNI","Output",
          "YAS_FAMILY_WMR_PERM_LOG_ACTUAL_DF.RDS"))


########################################################
###WE NEED INFORMATION ON THE NUMBER FOR EACH FAMILY###
#######################################################

SPLIT_FAM <- split(YAS_WMR_PERM_ACTUAL_DF,
                   YAS_WMR_PERM_ACTUAL_DF$fam)

###This excludes any family that doesn't have any significance
###MAKING THE DOT GRAPH
###
family_k = data.frame(fam=FAMILY_5$Family, sig = 0)
for (k in seq(1,43)){
  k1= ifelse(sum(str_count(SPLIT_FAM[[k]]$factor,"in"))==283,'none','some')
  family_k [k,2] <-k1 
  
}

none_sig <- subset(family_k, family_k$sig == 'none' )


FAMILY_YAS_PERM_PVAL_F <-       subset(YAS_WMR_PERM_ACTUAL_DF,
                                 !(YAS_WMR_PERM_ACTUAL_DF$fam %in% none_sig$fam))
FAMILY_YAS_PERM_PVAL_F$family =  droplevels(FAMILY_YAS_PERM_PVAL_F$fam)

############################################
###This makes the labeler for the number of#
############################################

###Makes a copy
Family_5_Nspecies <- FAMILY_5

###Makes the new labeler
Family_5_Nspecies$Label <- paste0(
        Family_5_Nspecies$Family," (",
        Family_5_Nspecies$Species, 
        ")", sep=" ")

###Merge with the original data.frame
FAMILY_YAS_PERM_PVAL_F<- 
          left_join(FAMILY_YAS_PERM_PVAL_F,
          Family_5_Nspecies[,c(1,3)],
          by = c("family"=c("Family")))

###SUBANNUAL PART
SUB_YAS_PERM <- subset(FAMILY_YAS_PERM_PVAL_F, 
                       FAMILY_YAS_PERM_PVAL_F$time=='sub')



SUB_YAS_PERM $Cut <- cut(SUB_YAS_PERM $scale , breaks = c(seq(0.1,1.2,0.02)),
                         dig.lab =4)
SUB_YAS_PERM $Cut <- gsub("[()]", "",SUB_YAS_PERM $Cut)
SUB_YAS_PERM $Cut <- gsub("\\[|\\]", "", SUB_YAS_PERM $Cut)
SUB_YAS_PERM$Cut <- factor(gsub("[,]","-",SUB_YAS_PERM$Cut))

#############################################
#############################################

###this makes it so that point color shows up if any
###scale in the interval is significant it makes it sigificant.
### Otherwise, the points overlap too much


SUB_YAS_PERM_SPLIT <- split(SUB_YAS_PERM, 
                      list(SUB_YAS_PERM$Cut,
                           SUB_YAS_PERM$fam))

###This gets rid of empty entries
SUB_YAS_PERM_SPLIT <- SUB_YAS_PERM_SPLIT [sapply(SUB_YAS_PERM_SPLIT , nrow) > 0]

for (k in seq(1,length(SUB_YAS_PERM_SPLIT))){
  tmp <-  SUB_YAS_PERM_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
    unique(tmp$synch_compen)[unique(tmp$synch_compen)!=("0")],
         tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }

  SUB_YAS_PERM_SPLIT[[k]] = tmp
}

SUB_YAS_PERM_DF <- do.call(rbind, SUB_YAS_PERM_SPLIT)

#############################################
#############################################

#############################################

SUB_YAS_GG <- 
  ggplot(SUB_YAS_PERM_DF, aes( x= factor(Cut), y= as.factor(Label)))+
  geom_point(aes(fill= synch_compen,size= as.factor(psig),
                           stroke = as.factor(psig),
                           color = as.factor(factor),
                  shape = as.factor(psig)))+
  ###For Point Size 
  scale_size_manual(values = c(3,2.3),
                    guide = "none")+
  # For Shape Size
  scale_shape_manual(values = c(21,23),
                     guide = "none")+
  
  ###For Stroke 
  scale_discrete_manual(
    aesthetics = "stroke",
    values = c(0,1.5),
    guide = "none")+
  ###For color
  scale_color_manual(values=c('white','black'),
                     guide = "none")+
  ###For color fill 
  scale_fill_manual(values = c('#f0f1f2','#24bcf0',"#fc4780"))+
  scale_y_discrete(limits =rev(levels(as.factor(SUB_YAS_PERM_DF$Label))))+
  scale_x_discrete(breaks =unique(SUB_YAS_PERM_DF$Cut)[seq(1, length(SUB_YAS_PERM_DF$Cut), by = 2)])+
  xlab("Scale (Year)")+
  ylab("Family")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.40,
                                   size=12),
        axis.text.y = element_text(color='black',size=12,vjust = 0.20),
        axis.title.x = element_text(size=14,color='black'),
        axis.title.y = element_text(size=14,color='black'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")
##############
###ANNNUAL####
##############

ANN_YAS_PERM <- subset(FAMILY_YAS_PERM_PVAL_F, 
                       FAMILY_YAS_PERM_PVAL_F$time=='ann' )

ANN_YAS_PERM $Cut <- cut(ANN_YAS_PERM $scale , breaks = c(seq(1,8.5,0.2)))
ANN_YAS_PERM $Cut <- gsub("[()]", "",ANN_YAS_PERM $Cut)
ANN_YAS_PERM $Cut <- gsub("\\[|\\]", "", ANN_YAS_PERM $Cut)
ANN_YAS_PERM$Cut <- factor(gsub("[,]","-",ANN_YAS_PERM$Cut) )



ANN_YAS_PERM_SPLIT <- split(ANN_YAS_PERM, 
                            list(ANN_YAS_PERM$Cut, ANN_YAS_PERM$family))


for (k in seq(1,length(ANN_YAS_PERM_SPLIT))){
  tmp = ANN_YAS_PERM_SPLIT[[k]]
  if(any(tmp$factor=='sig')==TRUE){
    tmp$factor = 'sig'
  }
  tmp$synch_compen<-
    ifelse(length(unique(tmp$synch_compen))>1,
           unique(tmp$synch_compen)[ unique(tmp$synch_compen)!=("0")],
           tmp$synch_compen)  
  
  
  if(any(tmp$psig == 1)==TRUE){
    tmp$psig = 1
  }
  
  ANN_YAS_PERM_SPLIT[[k]] = tmp
}

ANN_YAS_PERM_DF <- do.call(rbind, ANN_YAS_PERM_SPLIT)


#############################################


ANN_YAS_GG <-
  ggplot(ANN_YAS_PERM_DF, aes( x= factor(Cut), y= as.factor(Label)))+
  geom_point(aes(fill= synch_compen,size= as.factor(psig),
                 stroke = as.factor(psig),
                 color = as.factor(factor),
                 shape = as.factor(psig)))+
  ###For Point Size 
  scale_size_manual(values =  c(3,2.3),
                    guide = "none")+
  # For Shape Size
  scale_shape_manual(values = c(21,23),
                     guide = "none")+
  
  ###For Stroke 
  scale_discrete_manual(
    aesthetics = "stroke",
    values = c(0,1.1),
    guide = "none")+
  ###For color
  scale_color_manual(values=c('white','black'),
                     guide = "none")+
  ###For color fill 
  scale_fill_manual(values = c('#f0f1f2','#24bcf0',"#fc4780"))+
  scale_y_discrete(limits =rev(levels(as.factor(ANN_YAS_PERM_DF$Label))))+
  scale_x_discrete(breaks =unique(ANN_YAS_PERM_DF$Cut)[seq(1, length(ANN_YAS_PERM_DF$Cut), by = 2)])+
  
  xlab("Scale (Year)")+ylab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.40,
                                   size=12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14,color='black'),
        axis.title.y = element_text(size=14,color='black'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="bottom")

SUB_YAS_GG +ANN_YAS_GG+ plot_layout(
  guides = 'collect')& theme(legend.position = 'bottom')

#6.2 x 14




