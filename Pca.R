## clean ----
rm(list = ls())
pkgInitialization(c("ggplot2","dplyr","tidyr","gridExtra"))
source("uniTls_presetPaths.R");source("uniTls_pkgInstall.R");source("const.R");

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
    write.csv("pca/log/pcaLoading.csv")
  }
  rst.tot
}

pcaLoadingPlot <- function(dat, grouped = T, themeset, suffix = "", emphtag = NA){ ## output the pca loading plot
  library(dplyr);library(tidyr);library(ggplot2)
  circleFun <- function(center = c(0,0), r = 1, npoints = 100){
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  if (grouped) {
    dat <- dat%>% filter(group != "all")
  } else {
    dat <- dat%>% filter(group == "all")
  }
  zeros <- data.frame(RC1 = 0, RC2 = 0, tag = dat$tag, group = dat$group)
  dat2 <- rbind(zeros,dat)
  p <- ggplot() +
    geom_path(aes(x = RC1, y = RC2, group = tag, col = tag, alpha = tag, size = tag),
              arrow = arrow(angle = 15, length = unit(0.10, "inches"),
                            ends = "last", type = "open"),
              data = dat2) +
    geom_path(aes(x = x, y = y),col = "black", size = 0.7, linetype = 2, data = circleFun())+
    geom_text(aes(x = 1.1 * RC1, y = 1.1 * RC2, label = tag, col = tag, alpha = tag),
              check_overlap = F,data = dat) +
    facet_wrap(~ group) +
    xlim(-1.1, 1.1) +ylim(-1.1, 1.1) + 
    themeset
  taglv <- levels(dat$tag); ntaglv <- length(taglv)
  alphaSet <- rep(if(grouped) 0.3 else 1, ntaglv); sizeSet <- rep(0.5, ntaglv); colSet <- rep("black", ntaglv)
  if(grouped) {if(!is.na(emphtag)) {
    for(i in 1:length(emphtag)) {
      alphaSet[taglv == emphtag[i]] <- 1
      sizeSet[taglv == emphtag[i]] <- 0.75
      colSet[taglv == emphtag[i]] <- "red"
    }
  }}
  p <- p + scale_alpha_manual(breaks = taglv, values = alphaSet) +
    scale_size_manual(breaks = taglv, values = sizeSet) +
    scale_color_manual(breaks = taglv, values = colSet) 
  ggsave(plot = p,
         filename = paste("pca/plot/pcaLoading", suffix, ".png", sep = ""),
         width = 6, height = 6, dpi = 600)
}

pairPlot <- function(dat, xtag, ytag, gtag = "group") {
  library(ggplot2);library(dplyr);library(tidyr)
  x <- dat[xtag][,1]; y <- dat[ytag][,1]; g<- dat[gtag][,1] 
  g.lv <- as.character(unique(g))
  
  lmfitness <- NULL
  for(i in 1: length(g.lv)) {
    dat1 <- data.frame(x = x, y = y)[g == g.lv[i],]
    #plot(y ~ x, dat = dat1)
    mod <- lm(y ~ x, data = dat1); aov <- anova(mod)
    rsq <- aov$`Sum Sq`
    rsq.fmt<- paste("R square = ",format(100*rsq[1]/(sum(rsq)), nsmall  = 1, digits = 1),"%", sep = "")
    p <- aov$`Pr(>F)`[1]
    p.fmt <- if(p < 0.001) "p < 0.001" else paste("p = ",format(p, nsmall = 3, digits = 1), sep = "")
    coeff <- format(mod$coefficients, scientific = T, nsmall = 3, digits = 3)
    coeff.fmt <- paste("[",ytag,"] = ", coeff[2], "[",xtag,"] + ",coeff[1], sep = "")
    lmfitness <- rbind(lmfitness, 
                       data.frame(x.l = 0.5*max(dat1$x) + 0.5*min(dat1$x), 
                                  y.l = 0.9*max(dat1$y) + 0.1*min(dat1$y),
                                  g = g.lv[i], r.fmt = rsq.fmt, p = p, p.fmt = p.fmt, 
                       sig = (p < 0.05), coeff.fmt = if(p < 0.05) coeff.fmt else ""))
  }
  
  dat <- data.frame(x = x, y = y, g = g) %>% 
    inner_join(lmfitness[c("g","sig")], by = c("g" = "g"))
  
  ggplot(aes(x = x, y = y, col = g), data = dat) + 
    geom_point() + 
    geom_smooth(aes(fill = g, alpha = sig, linetype = sig), method = "lm", col = "black") + 
    geom_text(aes(label = paste(r.fmt,", ",p.fmt,"\n",coeff.fmt, sep = ""), 
                  x = x.l, y = y.l),
              size = 3, col = "black", data = lmfitness) +
    facet_wrap(~ g, scales = "free") + 
    scale_alpha_manual(values = c(0,0.4), breaks = c(1,0)) +
    scale_linetype_manual(values = c(2,1), breaks = c(1,0)) 
}

## Example ----
pcaLoading <- pcaLoadingCal(datareadln())

pcaLoadingPlot(pcaLoading, themeset = themeset, emphtag = c("S","Cd"), suffix = "_g")
pcaLoadingPlot(pcaLoading, grouped = F, themeset = themeset, suffix = "_a")

pairPlot(dat = datareadln(), xtag = "orgC", ytag = "Cd") +
  scale_x_continuous("organic carbon (mg/kg)") + 
  scale_y_continuous("Cd (mg/kg)") + 
  scale_color_manual("Location", breaks = c("CL","EA","NV","WE"), 
                     values = c("#B45F04","#31B404","grey50","#013ADF")) +
  scale_fill_manual("Location", breaks = c("CL","EA","NV","WE"), 
                    values = c("#B45F04","#31B404","grey50","#013ADF")) +
  themeset -> p
ggsave(filename = "pca/pair_Cd.png", plot = p, dpi = 600)

pairPlot(dat = datareadln(), xtag = "orgC", ytag = "S") +
  scale_x_continuous("organic carbon (mg/kg)") + 
  scale_y_continuous("S (mg/kg)") + 
  scale_color_manual("Location", breaks = c("CL","EA","NV","WE"), 
                     values = c("#B45F04","#31B404","grey50","#013ADF")) +
  scale_fill_manual("Location", breaks = c("CL","EA","NV","WE"), 
                    values = c("#B45F04","#31B404","grey50","#013ADF")) +
  themeset -> p 
ggsave(filename = "pca/pair_S.png", plot = p, dpi = 600)




  