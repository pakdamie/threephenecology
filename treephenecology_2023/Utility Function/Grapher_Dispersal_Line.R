###this makes the graph for the dispersal plots for yasuni and cc1



Grapher_Dispersal_Line <- function(dat){
  
  dat$time <- factor(dat$time, levels = c('sub','ann'))
  

  dat2 = dat
  #[seq(1,nrow(dat),2),]
  
  significant_points <- dat2[dat2$psig ==1,]
  significant_points_2 <- subset(significant_points,
                                 significant_points$actual!=0)
  
  
  dat2<-subset(dat2, dat2$scale < 2.1)
  
  plot<- 
    ggplot(dat2, aes(x= scale, y= perm_mean))+
    geom_ribbon(aes(x=scale,ymin= lb,ymax=ub),alpha = 0.4,fill='#73f5db') +
    #actual data- this should plotted
    geom_line(data=dat2, aes(x= scale, y= actual,color = synch_compen,group=1),
              size = 2, lineend = "round")+
    scale_color_manual(values = c('0'='darkgrey','com'='#24bcf0','synch'="#fc4780"))+
    geom_point(data=  significant_points_2,aes(x= scale, y= actual,fill=synch_compen),shape=21,size=3,
               )+
    scale_fill_manual(values = c('0'='darkgrey','com'='#24bcf0','synch'="#fc4780"))+
    
    ylim(0,1)+
    
    theme_classic()+
    facet_grid(.~time,scales = 'free')+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      panel.spacing = unit(0.2, "lines"), 
      panel.border=element_rect(colour="black",size=1,fill=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text=element_text(size=14, color='black'),
      axis.title=element_text(size=16, color='black'),
      plot.title = element_text(size = 18))+
    ylab("Average WMR")+xlab("Scale (Years)")
  
  
  return(plot)
}
