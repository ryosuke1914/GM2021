library(tidyverse)
library(Momocs)
#call source code to define new essential functions
source("morphospace for ggplot2.R")

#using sample data "bot" for outline data
outlines<-bot
stack(outlines)  
dev.off()
#================================================================================
#Elliptic fourier analysis
E<-efourier(outlines,12,norm=T)              
#PCA
P<-PCA(E)   
#define factors(groups)
Type<-fac_dispatcher(P,factor(outlines$fac$type))                  

PC1<-P$x[,"PC1"]
PC2<-P$x[,"PC2"]
PC3<-P$x[,"PC3"]
PC4<-P$x[,"PC4"]
da<-data.frame(PC1,PC2,PC3,PC4,Type)

#================================================================================
#using ggplot2 for PC plot with morphospaces
#install.pakages("ggThemeAssist") 
#to ease theme setting interactively. 
library(ggThemeAssist)
result<-morphospacePCA_custom(P,xax="PC1",yax="PC2",pos.shp="range",rotate.shp=pi/2,nb.shp=18, nr.shp = 6,nc.shp = 3)
result2<-result$gg+
  geom_point(cex=2,data=da, aes(x =PC1, y =PC2,colour=Type,shape=Type))+
  stat_ellipse(data=da, aes(x =PC1, y =PC2,colour=Type),level=0.8,lwd=1.0)+
  theme(panel.grid.major = element_line(colour = "gray90"), 
                panel.grid.minor = element_line(colour = "gray90"), 
                axis.title = element_text(size = 10), 
                plot.title = element_text(size = 10), 
                legend.text = element_text(size = 10), 
                panel.background = element_rect(fill = NA, 
                linetype = "twodash"), plot.background = element_rect(linetype = "twodash"), 
                legend.key = element_rect(fill = NA)) +labs(x = "PC1", y = "PC2")

#choose "result2" on row30 and click "Addins"→"ggThemeAssist"  
result2

#=================================================================================
#for interactive plot 
#install.packages("plotly")
library(plotly)
ggplotly(result2)

#save as...
ggsave("plot.sample.pdf",
       width = 400,
       height = 250,
       dpi = 300,
       units = "mm")












