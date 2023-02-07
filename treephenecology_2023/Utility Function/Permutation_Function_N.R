################################
#This is the function          #
#that does the permutation     #
#for all analyses- you sample  #
#species based on a group-     #
#should work for family        #
#and dispersal groups for both##
#CC1 and Yasuni                #
################################

###The actual data should be of the YASUNI_F
###and of CC1_F 

Permutation_Boot_Strap <- function(actualdata,dispersal,disp_n,log){
  
  #this inputs the actual data, the dispersal type,
  #and how many 
  
  ###Through the cluster-we'redoing 50 perms
  ###each and then combining them at the end
  ###to make it faster-this samples for the
  ###permutation

  RANDOM_SP <- replicate(1000,
              sample(dispersal,
              disp_n,
              replace = FALSE),
              simplify = FALSE
              )
  
  ###Pull out the important data
  SUBSETTED_DAT<- lapply(RANDOM_SP, function(Data) actualdata[,..Data])
  
  ###This pulls out the Date
  Date <- actualdata[,1]

  if(log==FALSE){
  
  ###This runs the WMR (Untransformed)
  WMR_PERM <- lapply(SUBSETTED_DAT, 
         function (sp) wmr(mvcwt(Date, sp,
                   min.scale = get.min.scale(Date),
                   max.scale = get.max.scale(Date)
                     )))
  }
  ###This runs the WMR (Transformed)
  
  if(log==TRUE){
    WMR_PERM <- lapply(SUBSETTED_DAT, 
                       function (sp) wmr(mvcwt(Date, log(sp+1),
                        min.scale = get.min.scale(Date),
                       max.scale = get.max.scale(Date)
                       )))
    
  }
  
  return(WMR_PERM )
  
}


Permutation_Actual_Combiner <- function(COI, perm, actual){
  
  #We get each of the perm element and find the z (WMR) and then
  #get the colmnean (this averages across time)
  
  
  z_1000_COI <- lapply(perm, function(x) x$z[,,1] * COI)
  
  z_1000 <-lapply(  z_1000_COI , function(x) colMeans(x,na.rm=TRUE)) #
  
  ###We then take make sure that it is now a data.frame with scale
  z_1000 <- lapply(z_1000, function(x) cbind.data.frame(Scale=perm[[1]]$y,
                                                        mean = x))
  ###We merge all 1000 into a singel data.frame 
  z_1000_MERGED <- Reduce(function(x,y) merge(x,y,by="Scale"), z_1000)
  
  ###We get the data.frame we're interested in and then we find the 
  ###mean of each row- so we now average across the 1000 to get a mean
  ###for each scale 
  WMR_1000_DF<-  cbind.data.frame(Scale =  
                                    perm[[1]]$y,
                                  mean = rowMeans((z_1000_MERGED[,-1])))
  
  ###Find the lowerbound and the upperbound
  WMR_1000_DF$lowerbound <-apply( z_1000_MERGED ,1,quantile, probs= c(0.025),na.rm=TRUE)
  WMR_1000_DF$upperbound <- apply( z_1000_MERGED ,1,quantile, probs= c(0.975),na.rm=TRUE)
  
  ###This is the data.frame we're interested in 
  DF_WMR_PA<- cbind.data.frame(scale = actual$y,
                                 actual= colMeans(COI * actual$z[,,1],na.rm=TRUE), #actual WMR
                                 perm_mean =  WMR_1000_DF$mean,
                                 lb = WMR_1000_DF$lowerbound,
                                 ub = WMR_1000_DF$upperbound)
  
  DF_WMR_PA $factor <- ifelse( DF_WMR_PA $actual>= 
                                 DF_WMR_PA $lb &
                                 DF_WMR_PA $actual <=  DF_WMR_PA $ub,
                               'insig','sig')
  
  ###This just allows me to seperate the whole data into being subannual (less than a year)
  ###or annual (more than a year)
  DF_WMR_PA $time <- ifelse( DF_WMR_PA $scale <= 1, 'sub', 'ann')
  
  
  ###This is for figuring out if the significant values are either
  ###compensatory dynamics or synchronous
  DF_WMR_PA $synch_compen <- 0
  DF_WMR_PA $synch_compen[   DF_WMR_PA $factor=='sig'&    DF_WMR_PA $actual  > 
                               DF_WMR_PA $ub] <- 'synch'
  DF_WMR_PA $synch_compen[   DF_WMR_PA $factor=='sig'&   DF_WMR_PA $actual
                             <    DF_WMR_PA $lb] <- 'com'  
  
  
  
  DF_WMR_PA$p_value <- 0
  
  ###This is for calculating the pvalue 
  for (c in seq(1,nrow(DF_WMR_PA))){
    count1 =table( z_1000_MERGED[c,-1] >=  DF_WMR_PA$actual[c])
    count2 =table( z_1000_MERGED[c,-1] <=   DF_WMR_PA$actual[c])
    
    
    p_1 <- try(c(2*(count1[[2]]/1000)))
    p_2 <- try(c(2 *(count2[[2]]/1000)))
    
    p_1=ifelse(is.character(p_1)==TRUE,0,p_1)
    p_2=ifelse(is.character(p_2)==TRUE,0,p_2)
    
    
    DF_WMR_PA$p_value[c]= min(p_1,p_2)
  }
  
  
  DF_WMR_PA$p_value_adjusted <-  p.adjust(  DF_WMR_PA$p_value, 
                                            method ="BY")
  
  #Significant or not?
  DF_WMR_PA$psig <- ifelse(DF_WMR_PA$p_value_adjusted  <= 0.05,
                           1,0.5)   
  return(DF_WMR_PA)
}
