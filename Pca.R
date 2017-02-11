## clean ----
rm(list = ls())
library(dplyr);library(tidyr);library(ggplot2);library(gridExtra);library(ade4):library(psych)

## Functions ----
datareadln <- function(sel) { ## data readln ----
  read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
    dplyr::select(N:Pb,SiteID,SplMonth,group)
}

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

## grouped pca
## ref...
dat <- datareadln() %>%
  filter(group != "EA") %>%
  filter(!is.na(N))%>%
  select(N:Pb)

fa.parallel(dat[,1:15],fa="pc", n.iter=10, ,show.legend=FALSE,main="Screen plot with parallel analysis")
tmp.pca <- principal(dat[,1:15],nfactors=3, scores = T, rotate="varimax")
tmp.pca 
tmp.pca <- dudi.pca(dat[,1:15],scann=F,nf = 3)
s.arrow(tmp.pca$c1, lab = names(dat[,1:15]),xax = 1, yax = 2, boxes = F,
        addaxes = T, grid = F, clabel = 1) 

## ----- EA
dat <- data %>% 
  inner_join(group) %>%
  filter(group == "EA") %>%
  filter(!is.na(N))%>%
  select(N:Pb)

fa.parallel(dat[,1:15],fa="pc", n.iter=10, ,show.legend=FALSE,main="Screen plot with parallel analysis")
tmp.pca <- principal(dat[,1:15],nfactors=2, scores = T, rotate="varimax")
tmp.pca 
tmp.pca <- dudi.pca(dat[,1:15],scann=F,nf = 3)
s.arrow(tmp.pca$c1, lab = names(dat[,1:15]),xax = 1, yax = 2, boxes = F,
        addaxes = T, grid = F, clabel = 1) 

##plot with cca
dat2 <- read.csv("elementPcaLoading.csv")
fix2pc <- as.data.frame(cbind(c("Al","As","C","Cd","Cr","Cu","Fe","Mn","N","Ni","orgC","P","Pb","S","Zn"),rep(0,15),rep(0,15)))
fix3pc <- as.data.frame(cbind(c("Al","As","C","Cd","Cr","Cu","Fe","Mn","N","Ni","orgC","P","Pb","S","Zn"),rep(0,15),rep(0,15),rep(0,15)))
colnames(fix2pc) <- c("Element","PC1","PC2")
colnames(fix3pc) <- c("Element","PC1","PC2","PC3")
fix2pc <- fix2pc %>%
  mutate(PC1 = as.numeric(0),PC2 = as.numeric(0))
fix3pc <- fix3pc %>%
  mutate(PC1 = as.numeric(0),PC2 = as.numeric(0),PC3 = as.numeric(0))

circle <- circleFun(c(0,0),2,npoints = 100)

## ----- ALL
datmp <- dat2%>%
  filter(group == "(CL+WE+NV+EA)") %>%
  spread(key = Pcaxe, value = Loading) %>%
  select(-group)
datmp1 <- rbind(fix2pc,datmp)
p1 <- ggplot() +
  geom_path(aes(x = PC1, y = PC2, group = Element),size = 0.7,alpha = 0.5,data = datmp1,
            arrow = arrow(angle = 15, length = unit(0.15, "inches"),
                          ends = "last", type = "open")) +
  geom_path(aes(x = x, y = y),col = "grey50", size = 0.7,linetype = 2, data = circle )+
  geom_text(aes(x = 1.1*PC1, y = 1.1*PC2, label = Element),check_overlap = F,data = datmp) +
  geom_text(aes(x = -1.0 , y = 1.0, label = "(a)"), size = 5) +
  xlim(-1.1,1.1) +ylim(-1.1,1.1) + 
  theme_bw() + theme(panel.grid = element_blank())

## ----- ref 
datmp <- dat2%>%
  filter(group == "(CL+WE+NV)") %>%
  spread(key = Pcaxe, value = Loading) %>%
  select(-group)
datmp1 <- rbind(fix3pc,datmp)
p2 <- ggplot() +
  geom_path(aes(x = PC1, y = PC2, group = Element),size = 0.7,alpha = 0.5,data = datmp1,
            arrow = arrow(angle = 15, length = unit(0.15, "inches"),
                          ends = "last", type = "open")) +
  geom_path(aes(x = x, y = y),col = "grey50", size = 0.7,linetype = 2, data = circle )+
  geom_text(aes(x = 1.1*PC1, y = 1.1*PC2, label = Element),check_overlap = F,data = datmp) +
  geom_text(aes(x = -1.0 , y = 1.0, label = "(b)"), size = 5) +
  xlim(-1.1,1.1) +ylim(-1.1,1.1) + 
  theme_bw() + theme(panel.grid = element_blank())

p3 <- ggplot() +
  geom_path(aes(x = PC1, y = PC3, group = Element),size = 0.7,alpha = 0.5,data = datmp1,
            arrow = arrow(angle = 15, length = unit(0.15, "inches"),
                          ends = "last", type = "open")) +
  geom_path(aes(x = x, y = y),col = "grey50", size = 0.7,linetype = 2, data = circle )+
  geom_text(aes(x = 1.1*PC1, y = 1.1*PC3, label = Element),check_overlap = F,data = datmp) +
  geom_text(aes(x = -1.0 , y = 1.0, label = "(c)"), size = 5) +
  xlim(-1.1,1.1) +ylim(-1.1,1.1) + 
  theme_bw() + theme(panel.grid = element_blank())

## ----- EA 
datmp <- dat2%>%
  filter(group == "(EA)") %>%
  spread(key = Pcaxe, value = Loading) %>%
  select(-group)
datmp1 <- rbind(fix2pc,datmp)
p4 <- ggplot() +
  geom_path(aes(x = PC1, y = PC2, group = Element),size = 0.7,alpha = 0.5,data = datmp1,
            arrow = arrow(angle = 15, length = unit(0.15, "inches"),
                          ends = "last", type = "open")) +
  geom_path(aes(x = x, y = y),col = "grey50", size = 0.7,linetype = 2, data = circle )+
  geom_text(aes(x = 1.1*PC1, y = 1.1*PC2, label = Element),check_overlap = F,data = datmp) +
  geom_text(aes(x = -1.0 , y = 1.0, label = "(d)"), size = 5) +
  xlim(-1.1,1.1) +ylim(-1.1,1.1) + 
  theme_bw() + theme(panel.grid = element_blank())

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2, widths=c(4,4), heights=c(4,4))