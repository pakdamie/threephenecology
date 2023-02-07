### The file needed [-1] because the rows are included
YASUNI_F <- fread(here("YASUNI", "Data", 
                       "YASUNI_ALL_FINAL.csv"))[,-1]


### reshape it
YASUNI_WIDE_F_ALL <- reshape2::dcast(YASUNI_F,
                                     Date ~ variable,
                                     value.var = "value",
                                     direction = "wide"
)


pdf("yas_ts.pdf",height = 12, width = 8, onefile = TRUE)
par(mfrow=c(5,2))

colnames(CodeJesse)[1] <- "CODIGO"

codejess_NULL= NULL
for (k in seq(2,1055)){
  
if ((colnames(YASUNI_WIDE_F_ALL)[k] %in% CodeJesse$CODIGO)==TRUE){
 name = CodeJesse[which(colnames(YASUNI_WIDE_F_ALL)[k] == CodeJesse$CODIGO),]$NombreActual

 if(is.na(name)==TRUE){
   name =colnames(YASUNI_WIDE_F_ALL)[k]
 }
 }
  else{
  name=  colnames(YASUNI_WIDE_F_ALL)[k]

  }
  codejess_NULL[[k]]=name
}

codejess_DF <- do.call(rbind, 
        codejess_NULL)

for (k in seq(2,1055)){


  plot(YASUNI_WIDE_F_ALL[,'Date'],
     log(YASUNI_WIDE_F_ALL[,k]+1),main=
       
       codejess_DF[k-1], 
     xlab = "Time",
     ylab = "Abundance (Log)",type='l')

}
dev.off()



#########################CC1############################
#########################
CC1_F <- data.frame(fread(here("CC1","Data","CC1_ALL_FINAL.csv")))[,-1]


pdf("cc1_ts.pdf",height = 12, width = 8, onefile = TRUE)
par(mfrow=c(5,2))

for (k in seq(2,655)){
  
  plot(CC1_F[,'Group.1'],
       log(CC1_F[,k]+1),main=
         colnames(CC1_F)[k], 
       xlab = "Time",
       ylab = "Abundance (Log)",type='l')
  
}
dev.off()
