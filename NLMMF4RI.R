# clean ----
rm(list = ls())

# package loading ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(nlme)
library(merTools)

# functions ----
datareadln <- function() { ## data readln ----
    read.csv("./Data/Result_Sediment.csv") %>%
    inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID"))%>%
    mutate(siteXmonth = paste(SplMonth,SiteID, sep = ""))
}

lmefit <- function(dat, indx, tag) { ## Linear mix model fiting ----
  dat <- dat %>% filter(elem == tag)
  if (sd(dat$value)==0) return(NA)
  mod.nlme <- lme(value ~ group + month, random = ~ 1|siteXmonth, data = dat) 
  est <-c("intercept","group","month"); element <- tag
  cbind(anova(mod.nlme), element, est, indx)
}

mainCal <- function(dat,indx) { ## main calculationg output
  result <- NULL
  for(i in 1:length(indx)) {
    indx1 <- indx[i]; dat1 <- dat
    dat1 <- dat1 %>%
      dplyr::select(SiteID,SplMonth,group,siteXmonth,contains(indx1)) %>%
      gather(elem, value, contains(indx1)) %>%
      mutate(elem = gsub(pattern = paste("_",indx1,sep=""),
                         replacement = "", x=elem)) %>%
      dplyr::rename(site = SiteID, month = SplMonth)
    elemcol <- unique(dat1$elem)
    for(j in 1:length(elemcol)) {
      result <-rbind(result, lmefit(dat = dat1, indx = indx1, tag = elemcol[j]))
    }
  }
  result
}

lmeFECal <- function(dat, indx, tag, archiplot) { ## Fixed effect
  dat <- dat %>% filter(elem == tag)
  mod.lme4 <- lmer(value ~ group + month + (1|siteXmonth), data = dat)
  #shinyMer(mod.lme4, simData = dat)
  fe <- FEsim(mod.lme4)
  if (archiplot) {
    ggsave(plot = plotFEsim(fe), 
           paste("Plot/",indx,"/",tag,"_Fixeff.png",sep=""),
           width = 6, height = 4)
  }
  fe
}

lmeRECal <- function(dat, indx, tag, archiplot) { ## Random effect
  dat <- dat %>% filter(elem == tag)
  mod.lme4 <- lmer(value ~ group + month + (1|siteXmonth), data = dat)
  re <- REsim(mod.lme4)
  if(archiplot) {
    ggsave(plot = plotREsim(re), 
           paste("Plot/",indx,"/",tag,"_Raneff.png",sep=""),
           width = 6, height = 4)
  }
  re
}

plotREsim2 <- function (data, level = 0.95, stat = "median", sd = TRUE, sigmaScale = NULL, 
                        oddsRatio = FALSE, labs = FALSE, taglv = NA, ncol) {
  ## plotREsim2, modified from merTools::plotREsim
  if (!missing(sigmaScale)) {
    data[, "sd"] <- data[, "sd"]/sigmaScale
    data[, stat] <- data[, stat]/sigmaScale
  }
  data[, "sd"] <- data[, "sd"] * qnorm(1 - ((1 - level)/2))
  data[, "ymax"] <- data[, stat] + data[, "sd"]
  data[, "ymin"] <- data[, stat] - data[, "sd"]
  data[, "sig"] <- data[, "ymin"] > 0 | data[, "ymax"] < 0
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    data[, "ymax"] <- exp(data[, "ymax"])
    data[, stat] <- exp(data[, stat])
    data[, "ymin"] <- exp(data[, "ymin"])
    hlineInt <- 1
  }
  data <- data[order(data[, "groupFctr"], data[, "term"], data[, 
                                                               stat]), ]
  rownames(data) <- 1:nrow(data)
  data[, "xvar"] <- factor(paste(data$groupFctr, data$groupID, 
                                 sep = ""), levels = unique(paste(data$groupFctr, data$groupID, 
                                                                  sep = "")), ordered = TRUE)
  if (labs == TRUE) {
    xlabs.tmp <- element_text(face = "bold", angle = 90, 
                              vjust = 0.5)
  }
  else {
    data[, "xvar"] <- as.numeric(data[, "xvar"])
    xlabs.tmp <- element_blank()
  }
  if (!is.na(taglv)) { 
    tmp <- NULL
    for (k in 1: length(taglv)) {
      tmp <- rbind(tmp, dplyr::filter(data, tag == taglv[k]))
    }
    data <- tmp
    data$tag <- factor(data$tag, levels = taglv)
  }
  p <- ggplot(data, aes_string(x = "xvar", y = stat, ymax = "ymax", 
                               ymin = "ymin")) +
    geom_hline(yintercept = hlineInt, color = I("red"), size = I(1.1)) +
    geom_point(color = "gray75", alpha = 1/(nrow(data)^0.33),  size = I(0.5)) +
    geom_point(data = subset(data, sig ==  TRUE), size = I(3)) + 
    labs(x = "Group", y = "Effect Range") + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.text.x = xlabs.tmp, 
                       axis.ticks.x = element_blank())
  if (sd) {
    p <- p + geom_pointrange(alpha = 1/(nrow(data)^0.33)) + 
      geom_pointrange(data = subset(data, sig == TRUE), 
                      alpha = 0.25)
  }
  p + facet_wrap(~tag + groupFctr, ncol = ncol, scales = "free")
}

plotFEsim2 <- function (data, level = 0.95, stat = "median", sd = TRUE, intercept = FALSE, 
                        sigmaScale = NULL, oddsRatio = FALSE, taglv = NA, glv, gcode, ncol) {
  ## plotFEsim2, modified from merTools::plotFEsim
  if (!missing(sigmaScale)) {
    data[, "sd"] <- data[, "sd"]/sigmaScale
    data[, stat] <- data[, stat]/sigmaScale
  }
  if (intercept == FALSE) {
    data <- data[data$term != "(Intercept)", ]
  }
  data[, "sd"] <- data[, "sd"] * qnorm(1 - ((1 - level)/2))
  data[, "ymax"] <- data[, stat] + data[, "sd"]
  data[, "ymin"] <- data[, stat] - data[, "sd"]
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    data[, "ymax"] <- exp(data[, "ymax"])
    data[, stat] <- exp(data[, stat])
    data[, "ymin"] <- exp(data[, "ymin"])
    hlineInt <- 1
  }
  xvar <- "term"
  data$term <- as.character(data$term)
  data$term <- factor(data$term, levels = data[order(data[, stat]), 1])
  if (!is.na(taglv)) { 
    tmp <- NULL
    for (k in 1: length(taglv)) {
      tmp <- rbind(tmp, dplyr::filter(data, tag == taglv[k]))
    }
    data <- tmp
    data$tag <- factor(data$tag, levels = taglv)
  }
  gcol <- rep("black",length(data$term)/length(taglv))
  gsize <- rep(1,length(data$term)/length(taglv))
  for (m in 1: length(glv)) {
    gcol[levels(data$term) == glv[m]] <- gcode[m]
    gsize[levels(data$term) == glv[m]] <- 3
  }
  p <- ggplot(aes_string(x = xvar, y = stat, ymax = "ymax", ymin = "ymin"), data = data) + 
    geom_hline(yintercept = hlineInt, color = I("red")) +
    geom_point(aes(col = term, size = term)) + 
    labs(x = "Group", y = "Fixed Effect") + 
    scale_color_manual("Group", breaks = levels(data$term), values = gcol) + 
    scale_size_manual("Group", breaks = levels(data$term), values = gsize) + 
    facet_wrap(~ tag, ncol = ncol) +
    coord_flip() + theme_bw() + theme(legend.position = "none")
  if (sd) {
    p <- p + geom_errorbar(width = 0.2)
  }
  p
}

mainPlot <- function(dat,indx,ncolfe,ncolre, glv, gcode, taglv = NA, archiplot = FALSE, ranecal = TRUE, suffix ="") { ## main plot output
  for(i in 1:length(indx)) {
    indx1 <- indx[i]; dat1 <- dat; fe <- NULL; re <- NULL;
    dat1 <- dat1 %>%
      dplyr::select(SiteID,SplMonth,group,siteXmonth,contains(indx1)) %>%
      gather(elem, value, contains(indx1)) %>%
      mutate(elem = gsub(pattern = paste("_",indx1,sep=""),
                         replacement = "", x=elem)) %>%
      dplyr::rename(site = SiteID, month = SplMonth)
    elemcol <- unique(dat1$elem)
    for(j in 1:length(elemcol)) {
      fe <- rbind(fe,cbind(lmeFECal(dat = dat1, indx = indx1, tag = elemcol[j], archiplot = archiplot),tag = elemcol[j]))
    }
    fe <- fe %>% filter(term != "(Intercept)") %>% mutate(term = gsub("group","",term)) %>% mutate(term = gsub("month","",term))
    p.fe <- plotFEsim2(data = fe, ncol = ncolfe[i], taglv = taglv, glv = glv, gcode = gcode)
    ggsave(filename = paste("Plot/",indx1,"_FixEff",suffix,".png", sep = ""),
           plot = p.fe, dpi = 600, width = 8, height = 6)
    if (ranecal) {
      for(j in 1:length(elemcol)) {
        re <- rbind(re,cbind(lmeRECal(dat = dat1, indx = indx1, tag = elemcol[j], archiplot = archiplot),tag = elemcol[j]))
      }
      p.re <- plotREsim2(data = re, ncol = ncolre[i], taglv = taglv)
      ggsave(filename = paste("Plot/",indx1,"_RanEff.png", sep = ""),
             plot = p.re, dpi = 600, width = 12, height = 6)
    }
  }
  NULL
}

## example ----
mainCal(dat = datareadln()[-84:-85,], indx = c("Igeo","EnrichmentFator")) %>%
  write.csv("log/ModelFiting.csv")
mainPlot(dat = datareadln()[-84:-85,], indx = "Igeo", ncolfe = 5, ncolre = 5, archiplot = T, ranecal = T, suffix = "_G",
         taglv = c("N","C","orgC","S","P","Al","Fe","Mn","Cu","Zn","Ni","Cr","Pb","As","Cd"),
         glv = c("EA","NV","WE"), gcode = c("#31B404","grey50","#013ADF"))
mainPlot(dat = datareadln()[-84:-85,], indx = "EnrichmentFator",ncolfe = 4,ncolre = 4, archiplot = T, ranecal = T,  suffix = "_G",
         taglv = c("N","C","S","Fe","Mn","Cu","Zn","Ni","Cr","Pb","As","Cd"),
         glv = c("EA","NV","WE"), gcode = c("#31B404","grey50","#013ADF"))
mainPlot(dat = datareadln()[-84:-85,], indx = "Igeo", ncolfe = 5, ncolre = 5, archiplot = F, ranecal = F,  suffix = "_M",
         taglv = c("N","C","orgC","S","P","Al","Fe","Mn","Cu","Zn","Ni","Cr","Pb","As","Cd"),
         glv = c("Nov","Jul"), gcode = c("brown","green"))
mainPlot(dat = datareadln()[-84:-85,], indx = "EnrichmentFator",ncolfe = 4,ncolre = 4, archiplot = F, ranecal = F,  suffix = "_M",
         taglv = c("N","C","S","Fe","Mn","Cu","Zn","Ni","Cr","Pb","As","Cd"),
         glv = c("Nov","Jul"), gcode = c("brown","green"))
