## ======================cleaning==========================================
library(ggplot2)
library(maptools)
library(grid)
library(ggmap)
library(dplyr)
library(tidyr)
library(maptools)
rm(list=ls())

## ===========================data readln==================================
##========================================================================
conveySpCoord <- function(x) {
  tmp <- strsplit(as.character(x),"°")
  out <- NULL
  for(i in 1:length(tmp)) {
    out <- c(out,as.numeric(tmp[[i]][1]) + as.numeric(substr(tmp[[i]][2],1,nchar(tmp[[i]][2])-1))/60)
  }
  out
}



cname <- c("Longitude","Latitude")

coo.land <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Land.kml", ignoreAltitude = T)[[1]])
colnames(coo.land) <- cname

coo.sm1 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh1.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm1) <- cname

coo.sm2 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh2.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm2) <- cname

coo.sm3 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh3.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm3) <- cname

coo.sm4 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh4.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm4) <- cname

coo.sm5 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh5.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm5) <- cname

coo.sm6 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh6.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm6) <- cname

coo.sm7 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Saltmarsh7.kml", ignoreAltitude = T)[[1]])
colnames(coo.sm7) <- cname

coo.sw1 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall1.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw1) <- cname

coo.sw2 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall2.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw2) <- cname

coo.sw3 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall3.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw3) <- cname

coo.sw4 <- as.data.frame(
  getKMLcoordinates("Data/map/RudongSiteMap_Seawall4.kml", ignoreAltitude = T)[[1]])
colnames(coo.sw4) <- cname

coo.site <- read.csv("Data/map/sitecor.csv")
#coo.site$col <- as.character(coo.site$col)
coo.site$Longitude <- conveySpCoord(as.character(coo.site$Longitude))
coo.site$Latitude <- conveySpCoord(as.character(coo.site$Latitude))

bkmap <- readShapePoly("data/bou2_4p.shp")

##================= mapping============================
##============================================================
ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_point(data = coo.site , 
             aes(x = Longitude, y = Latitude, col = col ,size = 8, shape = group))+
  scale_x_continuous(breaks = 121.25 + 0:4 * 0.005, 
                     labels = paste("121°",c(15 + 0:4 * 0.3),"'E",sep=""))+
  geom_vline(xintercept = 121.25 + 0:4 * 0.005, alpha = 0.5)+
  scale_y_continuous(breaks = 32.45 + 0:4 * 0.005, 
                     labels = paste("32°",c(27 + 0:4 * 0.3),"'N",sep=""))+
  geom_hline(yintercept = 32.45 + 0:4 * 0.005, alpha = 0.5) +
  scale_size_identity()+
  scale_fill_identity() +
  scale_colour_identity()+
  scale_shape_manual(values = c(16,17,18,15))+
  coord_map(xlim = c(121.2475, 121.2725), ylim = c(32.4475, 32.4725),
            projection = "azequidistant") +
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "none") 
ggsave(filename = "map/line_1.png", dpi = 600, width = 6, height = 6)


ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#C1FFC1"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "black"))+
  geom_point(data = coo.site, 
             aes(x = Longitude, y = Latitude, col = col ,size = 8, shape = group))+
  scale_x_continuous(breaks = 121.4125 + 0:2 * 0.005, 
                     labels = paste("121°",(0.4125 + 0:2 * 0.005) * 60 ,"'E",sep=""))+
  geom_vline(xintercept = 121.4125 + 0:2 * 0.005, alpha = 0.5) +
  scale_y_continuous(breaks = 32.365 + 0:2 * 0.005, 
                     labels = paste("32°",(0.365 + 0:2 * 0.005) * 60,"'N",sep=""))+
  geom_hline(yintercept = 32.365 + 0:2 * 0.005, alpha = 0.5) +
  scale_size_identity()+
  scale_fill_identity()+
  scale_colour_identity()+
  scale_shape_manual(values = c(16,17,18,15))+
  coord_map(xlim = c(121.41, 121.425), ylim = c(32.3625, 32.3775),
            projection = "azequidistant")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "none") 
ggsave(filename = "map/line_2.png", dpi = 600, width = 6, height = 6)


ggplot()+
  geom_polygon(data = coo.land,
               aes(x = Longitude, y = Latitude, fill = "#B7B7B7"))+
  geom_polygon(data = coo.sm1,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm2,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm3,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm4,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm5,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm6,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sm7,
               aes(x = Longitude, y = Latitude, fill = "#698B22"))+
  geom_polygon(data = coo.sw1,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw2,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw3,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_polygon(data = coo.sw4,
               aes(x = Longitude, y = Latitude, fill = "#555555"))+
  geom_point(data = coo.site , 
             aes(x = Longitude, y = Latitude, col = col ,size = 2))+
  geom_path(aes(x = c(121.23+0.01,121.23+0.01+5/(pi*6378.140*cos(32/180*pi))*180), 
                y = c(32.36,32.36)), size = 1.7, color = "black")+
  geom_text(aes(x = 121.23+0.01+2.5/(pi*6378.140*cos(32/180*pi))*180, y = 32.365), 
            label = "5 km", size = 7, color = "black") + 
  geom_text(aes(x = 121.43, y = 32.47), 
            label = "N", size = 9, angle = 0, color = "black") + 
  geom_path(aes(x = c(121.43,121.43), y = c(32.43,32.46)),size = 1.3, arrow = arrow(), color = "black") + 
  geom_path(aes(x = c(121.42,121.44), y = c(32.445,32.445)),size = 1.3, color = "black") + 
  scale_x_continuous(breaks = 121.235 + 0:6 * 0.035, 
                     labels = paste("121°",(0.235 + 0:6 * 0.035) * 60 ,"'E",sep=""))+
  #geom_vline(xintercept = 121.235 + 0:6 * 0.035, alpha = 0.5) +
  scale_y_continuous(breaks = 32.355 + 0:3 * 0.040, 
                     labels = paste("32°",(0.355 + 0:3 * 0.040) * 60,"'N",sep=""))+
  #geom_hline(yintercept = 32.355 + 0:3 * 0.040, alpha = 0.5) +
  coord_map(xlim =  c(121.23, 121.45), ylim = c(32.35, 32.48),
            projection = "azequidistant")+
  scale_size_identity()+
  scale_fill_identity() +
  scale_colour_identity()+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(),
        panel.background = element_rect(fill = "#B9D3EE"),
        panel.grid = element_blank(),
        legend.position = "bottom")  
ggsave(filename = "map/base.png", dpi = 600, width = 12, height = 8)




ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey90", data = fortify(bkmap)) + 
  geom_point(aes(x = 121.25, y = 32.5), size = 4, color = "red") + 
  geom_text(aes(x = 121.25, y = 32.5), nudge_x = 0.8, nudge_y = 0.7,
            label = "Rudong", size = 5, angle = -45, color = "black") + 
  geom_text(aes(x = 125, y = 35), 
            label = "YELLOW SEA", size = 6, angle = -20, color = "black") + 
  geom_text(aes(x = 125, y = 29), 
            label = "EAST CHINA SEA", size = 6, angle = 47, color = "black") + 
  geom_text(aes(x = 120.3, y = 39), 
            label = "BOHAI SEA", size = 5, angle = 60, color = "black") + 
  geom_text(aes(x = 115, y = 31), 
            label = "CHINA", size = 17, angle = 20, color = "grey30") + 
  geom_text(aes(x = 129, y = 40), 
            label = "N", size = 9, angle = 0, color = "black") + 
  geom_path(aes(x = c(129,129), y = c(35,39)),size = 1.3, arrow = arrow(), color = "black") + 
  geom_path(aes(x = c(128,130), y = c(37,37)),size = 1.3, color = "black") + 
  coord_quickmap(xlim = c(110,130),ylim = c(25,40)) +
  scale_x_continuous("",breaks = 0:4 * 5 + 110, 
                     labels = paste(0:4 * 5 + 110,"°E", sep = "")) + 
  scale_y_continuous("",breaks = 0:3 * 5 + 25, 
                     labels = paste(0:3 * 5 + 25,"°N", sep = "")) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill = blues9[3]))
ggsave(filename = "map/rudong.png",height = 5, width = 6)
