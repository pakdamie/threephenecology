################################
####REQUIRED PACKAGES ##########
################################
require(plyr)
require(ggplot2)
require(magrittr)
require(reshape2)
require(dplyr)
require(mvcwt)
require(here)
require(patchwork)
require(readr)
require(stringr)
###############################
###FAMILY FOR CC1 ##########
###############################

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

#####################################################################
#####################################################################
family_count_CC1<- dplyr::count(species_name_CC1_Found, family)

###Only looking at families with more than five species

family_five_CC1 <- subset(family_count_CC1, family_count_CC1$n >=5)

###There's 32 families
 ######################################################################
######################################################################

FAM_SPECIES_FOUND <- 
  subset(species_name_CC1_Found,
         species_name_CC1_Found$family %in% family_five_CC1$family)

CC1_Species_Code <- FAM_SPECIES_FOUND$name

CC1_Agg[,1] <- CC1_Agg$Day_Year_Y

CC1_Agg <- data.table(CC1_Agg)
NUMBER_NEEDED_CC1 <- unique(family_five_CC1 $n)


for (k in seq(1,19)){
  
  tmp <- Permutation_Boot_Strap(
    actualdata=CC1_Agg,
    CC1_Species_Code, 
    disp_n= NUMBER_NEEDED_CC1[[k]], 
    log=TRUE)
  
  assign(paste("PERM", NUMBER_NEEDED_CC1[[k]], sep = ""), tmp)
  
  save(list=paste("PERM", NUMBER_NEEDED_CC1[[k]], sep = ""),
       file =paste("PERM", NUMBER_NEEDED_CC1[[k]],".RData", sep = ""))
  rm(list=paste("PERM", NUMBER_NEEDED_CC1[[k]], sep = ""))

####################################################
#####################################################
#####################################################


#SET WORKING DIRECTORY TO P:\treephenecology\CC1\Output\FAMILY_PERM_LOG
lofs <-list.files( here("CC1","Output",
                        "FAMILY_PERM_LOG"),pattern = ".RData")

  NUMBER_ORDER <- parse_number(lofs)
  
  #If you want the transformed file
#lofs <-list.files( here("CC1","Output",
# FAMILY_PERM","Perm"),pattern = ".RData")

###Notice that's it's not in order
NUMBER_ORDER_CC1 <- parse_number(lofs)

WMR_CC1_PERM_DF<- NULL

for (r in seq(1,nrow(family_five_CC1))){
  
  FAMILY_NAME <-   family_five_CC1 [r,1]
  FAMILY_NUMBER <- family_five_CC1 [r,2]
  
  FAM_CODE <-subset(species_name_CC1_Found, 
                    species_name_CC1_Found$family ==   
                    FAMILY_NAME)
  
  fam_df <- cbind.data.frame(time=CC1_Agg$Day_Year_Y, 
                             as.data.frame(CC1_Agg)[FAM_CODE $name])
  
  
  x = fam_df$time
  ### For the y, you must rely on the species data from the 2nd column to the end of the
  ###column- LOG TRANSFORMED
  y =   log(fam_df[,2:ncol(fam_df)]+1)
  
  WMR_FAM = wmr(mvcwt(x, y, min.scale = get.min.scale(x), 
                      max.scale = get.max.scale(x)))
  
  index =which( NUMBER_ORDER == FAMILY_NUMBER)
  
  tmp = load(here("CC1","Output",
                  "FAMILY_PERM_LOG",
                  lofs[[index]]))
  
  PERM = get(tmp)
  
  
  WMR_CC1_PERM_DF[[r]]<-
    Permutation_Actual_Combiner(MVCWT_CC1_WMR_COI,PERM,WMR_FAM )
  
}


for (f in seq(1, nrow( family_five_CC1 ))){
  
  WMR_CC1_PERM_DF[[f]]$fam = family_five_CC1$family[[f]]
}

WMR_CC1_PERM_DF<-do.call(rbind, WMR_CC1_PERM_DF)

saveRDS(WMR_CC1_PERM_DF,file ='WMR_CC1_PERM_DF.RDS')

SPLIT_FAM_CC1 <- split( WMR_CC1_PERM_DF,
                        WMR_CC1_PERM_DF$fam)
###MAKING THE DOT GRAPH
family_k = data.frame(fam=family_five_CC1$family, sig = 0)
for (k in seq(1,32)){
  k1= ifelse(sum(str_count(SPLIT_FAM_CC1[[k]]$factor,"in"))==196,'none','some')
  family_k [k,2] <-k1 
  
}
###Non significant families 
none_sig <- subset(family_k, family_k$sig == 'none' )

###############################
###This makes the label part###
###############################

############################################
###This makes the labeler for the number of#
############################################

###Makes a copy
Family_5_Nspecies <- family_five_CC1

###Makes the new labeler
Family_5_Nspecies$Label <- paste0(
  Family_5_Nspecies$family," (",
  Family_5_Nspecies$n, 
  ")", sep=" ")

###Merge with the original data.frame
WMR_CC1_PERM_DF_F<- 
  left_join(WMR_CC1_PERM_DF,
            Family_5_Nspecies[,c(1,3)],
            by = c("fam"=c("family")))

###########################################
###########################################


FAMILY_CC1_PERM_F <- subset(WMR_CC1_PERM_DF_F,
                                 !(WMR_CC1_PERM_DF_F$fam %in% none_sig$fam))


FAMILY_CC1_PERM_F <- na.omit(FAMILY_CC1_PERM_F )
#FAMILY_CC1_PERM_PVAL_F$family =droplevels(FAMILY_CC1_PERM_PVAL_F$family)

SUB_CC1_PERM <- subset(FAMILY_CC1_PERM_F, 
                       FAMILY_CC1_PERM_F$time=='sub')



SUB_CC1_PERM $Cut <- cut(SUB_CC1_PERM  $scale , breaks = c(seq(0.08,1.2,0.02)))
SUB_CC1_PERM $Cut <- gsub("[()]", "",SUB_CC1_PERM  $Cut)
SUB_CC1_PERM $Cut <- gsub("\\[|\\]", "", SUB_CC1_PERM $Cut)
SUB_CC1_PERM $Cut <- factor(gsub("[,]","-",SUB_CC1_PERM $Cut) )

SUB_CC1_PERM$fam <- as.factor(SUB_CC1_PERM$fam)
#############################################
#############################################

###this makes it so that point color shows up if any
###scale in the interval is significant it makes it sigificant.
### Otherwise, the points overlap too much


SUB_CC1_PERM_SPLIT <- split(SUB_CC1_PERM, 
                            list(SUB_CC1_PERM$Cut,
                                 SUB_CC1_PERM$fam))

###This gets rid of empty entries
SUB_CC1_PERM_SPLIT <- SUB_CC1_PERM_SPLIT [sapply(SUB_CC1_PERM_SPLIT , nrow) > 0]

for (k in seq(1,length(SUB_CC1_PERM_SPLIT))){
  tmp <-  SUB_CC1_PERM_SPLIT[[k]]
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
  
  SUB_CC1_PERM_SPLIT[[k]] = tmp
}

SUB_CC1_PERM_DF <- do.call(rbind, SUB_CC1_PERM_SPLIT)

#############################################
#############################################

#############################################

SUB_CC1_GG <- 
  ggplot(SUB_CC1_PERM_DF, aes( x= factor(Cut), y= as.factor(Label)))+
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
  scale_y_discrete(limits =rev(levels(as.factor(SUB_CC1_PERM_DF$Label))))+
  scale_x_discrete(breaks =unique(SUB_CC1_PERM_DF$Cut)[seq(1, length(SUB_CC1_PERM_DF$Cut), by = 2)])+
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

###ANNNUAL

ANN_CC1_PERM <- subset(FAMILY_CC1_PERM_F, 
                       FAMILY_CC1_PERM_F$time=='ann')


ANN_CC1_PERM $Cut <- cut(ANN_CC1_PERM $scale , breaks = c(seq(1,2.2,0.2)))
ANN_CC1_PERM $Cut <- gsub("[()]", "",ANN_CC1_PERM $Cut)
ANN_CC1_PERM $Cut <- gsub("\\[|\\]", "", ANN_CC1_PERM $Cut)
ANN_CC1_PERM$Cut <- factor(gsub("[,]","-",ANN_CC1_PERM$Cut) )


#############################################


ANN_CC1_PERM_SPLIT<- split(ANN_CC1_PERM, 
                            list(ANN_CC1_PERM$Cut, ANN_CC1_PERM$fam))


for (k in seq(1,length(ANN_CC1_PERM_SPLIT))){
  tmp = ANN_CC1_PERM_SPLIT[[k]]
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
  
  ANN_CC1_PERM_SPLIT[[k]] = tmp
}

ANN_CC1_PERM_DF <- do.call(rbind, ANN_CC1_PERM_SPLIT)


#############################################


ANN_CC1_GG <-
  ggplot(ANN_CC1_PERM_DF, aes( x= factor(Cut), y= as.factor(Label)))+
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
  scale_y_discrete(limits =rev(levels(as.factor(ANN_CC1_PERM_DF$Label))))+

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

SUB_CC1_GG +ANN_CC1_GG+ plot_layout(
  guides = 'collect')& theme(legend.position = 'bottom')

