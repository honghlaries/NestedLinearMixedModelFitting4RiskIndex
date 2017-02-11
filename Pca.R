## clean ----
rm(list = ls())
folderCreating <- function(dirs) {
  for(i in 1:length(dirs)) {
    if(!dir.exists(dirs[i])) dir.create(dirs[i])
  }
}
folderCreating(dirs = c("pca","pca/plot","pca/log"))

## Functions ----
datareadln <- function() { ## data readln
  library(dplyr);library(tidyr)
  read.csv("Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("Data/meta_Quadrat4pca.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
    dplyr::select(N:Pb,SiteID,SplMonth,group)
}

pcaLoadingCal <- function(dat, grouped = T, log = T) { ## cal and log the pca loading
  library(dplyr);library(tidyr);library(psych)
  dat1 <- dat %>%
    filter(!is.na(N))%>%
    select_("N:Pb")
  png(paste("pca/plot/screenPlot_","all",".png",sep = ""))
  fa.parallel(dat1,fa="pc", n.iter=10, show.legend=T,
              main = paste("Screen plot with parallel analysis\n, group =","all"))
  dev.off()
  pca <- principal(dat1, nfactors=2, scores = T, rotate="varimax")
  tag <- dimnames(pca$loadings)[[1]]
  loadings <- matrix(as.numeric(pca$loadings), ncol = 2, 
                     dimnames = list(tag, c("RC1","RC2"))) 
  rst.tot <- data.frame(loadings,tag, group = "all")
  if (grouped) {
    grouptag <- levels(dat$group)
    for(i in 1:length(grouptag)) {
      dat1 <- dat %>%
        filter(group == grouptag[i]) %>%
        filter(!is.na(N))%>% 
        select_("N:Pb")
      png(paste("pca/plot/screenPlot_",grouptag[i],".png",sep = ""))
      fa.parallel(dat1,fa="pc", n.iter=10, show.legend=T,
                  main = paste("Screen plot with parallel analysis\n, group =",grouptag[i]))
      dev.off()
      pca <- principal(dat1, nfactors=2, scores = T, rotate="varimax")
      tag <- dimnames(pca$loadings)[[1]]
      loadings <- matrix(as.numeric(pca$loadings), ncol = 2, 
                         dimnames = list(tag, c("RC1","RC2"))) 
      rst.tot <- rbind(rst.tot, data.frame(loadings,tag, group = grouptag[i]))
    }
  }
  if(log) {
    write.csv(pca/log/pcaLoading.csv)
  }
  rst.tot
}

pcaLoadingPlot <- function(dat, grouped = T, themeset, suffix = ""){ ## output the pca loading plot
  library(dplyr);library(tidyr);library(ggplot2)
  circleFun <- function(center = c(0,0), r = 1, npoints = 100){
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  if (grouped) {
    dat <- dat%>%
      select(group != "all")
  } else {
    select(group == "all")
  }
  zeros <- data.frame(RC1 = 0, RC2 = 0, tag = dat$tag, group = dat$group)
  dat2 <- rbind(zeros,dat)
  p <- 
  ggsave(plot = ggplot() +
           geom_path(aes(x = RC1, y = RC2, group = tag),
                     arrow = arrow(angle = 15, length = unit(0.15, "inches"),
                                   ends = "last", type = "open"),
                     size = 0.7, data = dat2) +
           geom_path(aes(x = x, y = y),col = "black", size = 0.7, linetype = 2, data = circleFun())+
           geom_text(aes(x = 1.1 * RC1, y = 1.1 * RC2, label = tag),
                     check_overlap = F,data = dat) +
           facet_wrap(~ group) +
           xlim(-1.1, 1.1) +ylim(-1.1, 1.1) + 
           themeset,
         filename = paste("pca/plot/pcaLoading", suffix, ".png", sep = ""),
         width = 6, height = 6, dpi = 600)
}

## Example
pcaLoading <- pcaLoadingcal(datareadln())
library(ggolot2)
themeset <- theme_bw() + theme(aspect.ratio = 1, plot.background = element_blank())
pcaLoadingPlot(pcaLoading, themeset = themeset, suffix = "_g_tp")
pcaLoadingPlot(pcaLoading, grouped = F, themeset = themeset, suffix = "_a_tp")
themeset <- theme_bw() + theme(aspect.ratio = 1, plot.background = element_blank())
pcaLoadingPlot(pcaLoading, themeset = themeset, suffix = "_g")
pcaLoadingPlot(pcaLoading, grouped = F, themeset = themeset, suffix = "_g")