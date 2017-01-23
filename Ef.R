## clean ----
rm(list = ls())
library(dplyr);library(tidyr)
folderCreating <- function(dirs) {
  for(i in 1:length(dirs)) {
    if(!file.exists(dirs[i])) dir.create(dirs[i])
  }
}
folderCreating(dirs = c("ef","ef/log","ef/plot","ef/plot/Raneff","ef/plot/Fixeff","ef/plot/diag","ef/plot/den"))

## Functions ----
datareadln <- function() { ## data readln ----
  read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::right_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
    dplyr::select(N_EnrichmentFator:Pb_EnrichmentFator,SiteID,SplMonth,group) %>%
    tidyr::gather(elem, resp, N_EnrichmentFator:Pb_EnrichmentFator) %>%
    dplyr::mutate(elem = gsub("_EnrichmentFator", "", elem, fixed = T)) 
}

meanseCal <- function(dat) { ## calulate and output mean and se
  se <- function(v) {
    if(length(v) < 2) {NA} else {sd(v, na.rm = T)/sqrt(length(v)-sum(is.na(v)))}
  }
  dat %>%
    dplyr::group_by(SiteID, group, SplMonth, elem) %>%
    dplyr::summarise(avg = mean(resp,na.rm = T),se = se(resp)) %>%
    tidyr::gather(stat, value, avg:se) %>% 
    dplyr::mutate(elemstat = paste(elem,stat,sep = "_")) %>%
    dplyr::select(-elem,-stat) %>%
    tidyr::spread(elemstat,value) %>% 
    write.csv("ef/log/IndexStat.csv")
}

shapiroTest <- function(dat, tag, grouplv, SplMonthlv) {
  dat <- dat %>% filter(elem == tag) %>% filter(!is.na(resp))
  shapiro <- NULL
  for(i in 1:length(grouplv)) {
    for(j in 1:length(SplMonthlv)) {
      subdat <- dat %>% filter(group == grouplv[i]) %>% filter(SplMonth == SplMonthlv[j])
      result <- stats::shapiro.test(subdat$resp)
      shapiro <- rbind(shapiro, data.frame(tag = tag,
                                           group = grouplv[i],
                                           SplMonth = SplMonthlv[j],
                                           W = result$statistic,
                                           Pvalue = result$p.value))
    }
  }
  shapiro
}

modelCompare <- function(dat, tag, fact, log = TRUE) {
  library(lme4)
  dat <- dat %>% filter(elem == tag) %>% filter(!is.na(resp))
  formulaList <- NULL; modlog <- NULL
  for (j in 1:length(fact)) {
    formulaList <- c(formulaList, as.formula(paste("resp","~",fact[j])))
  }
  mod.champion <- lmer(formulaList[[1]], data = dat)
  for (i in 2:length(formulaList)) {
    mod.challenger <- lmer(formulaList[[i]], data = dat)
    comp <- anova(mod.champion, mod.challenger)
    mod <- as.character(formula(mod.champion@call))[3]
    modp <- as.character(formula(mod.challenger@call))[3]
    modlog <- rbind(modlog, cbind(comp,
                                  model = c(paste(tag,"~",mod),paste(tag,"~",modp)),
                                  ref = c("champion", "challenger"),
                                  result = if (comp$AIC[2] < comp$AIC[1]) c("", "Win") else c("Win", ""),
                                  time = date())) 
    if (comp$AIC[2] < comp$AIC[1]) mod.champion <- mod.challenger
  }
  if (!file.exists("ef/log/ModelComparing.csv")) {
    write.table(t(c("Df","AIC","BIC","Loglik","deviance","Chisq","ChiDf","Pr(>Chisq)","model","ref","result","time")),
                file = "ef/log/ModelComparing.csv", sep = ",", row.names = F, col.names = F)
  } 
  if(log) {
    write.table(modlog, append = T, sep = ",", row.names = F, col.names = F,
                file = "ef/log/ModelComparing.csv")
  }
  mod.champion
}

plotFEsim2 <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                        sigmaScale = NULL, oddsRatio = FALSE, glv, gcode, theme) {
  ## plotFEsim2, modified from merTools::plotFEsim
  if (!missing(sigmaScale)) {
    fe[, "sd"] <- fe[, "sd"]/sigmaScale
    fe[, stat] <- fe[, stat]/sigmaScale
  }
  fe[, "sd"] <- fe[, "sd"] * qnorm(1 - ((1 - level)/2))
  fe[, "ymax"] <- fe[, stat] + fe[, "sd"]
  fe[, "ymin"] <- fe[, stat] - fe[, "sd"]
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    fe[, "ymax"] <- exp(fe[, "ymax"])
    fe[, stat] <- exp(fe[, stat])
    fe[, "ymin"] <- exp(fe[, "ymin"])
    hlineInt <- 1
  }
  xvar <- "term"
  fe$term <- as.character(fe$term)
  fe$term <- factor(fe$term, levels = fe[order(fe[, stat]), 1])
  
  gcol <- rep("black",length(levels(fe$term)))
  gsize <- rep(1,length(levels(fe$term)))
  gltp <- rep(2,length(levels(fe$term)))
  for (m in 1: length(glv)) {
    gcol[levels(fe$term) == glv[m]] <- gcode[m]
    gsize[levels(fe$term) == glv[m]] <- 4
  }
  llv <- c("Nov","Jul","EA","NV","WE")
  for (m in 1: length(llv)) {
    gltp[levels(fe$term) == llv[m]] <- 1
  }
  p <- ggplot(aes_string(x = xvar, y = stat, ymax = "ymax", ymin = "ymin"), data = fe) + 
    geom_hline(yintercept = hlineInt, color = I("red")) +
    geom_point(aes(col = term, size = term)) + 
    labs(x = "Group", y = "Fixed Effect") + 
    scale_color_manual("Group", breaks = levels(fe$term), values = gcol) + 
    scale_size_manual("Group", breaks = levels(fe$term), values = gsize) + 
    coord_flip() + 
    theme
  if (sd) {
    p <- p + geom_errorbar(aes(linetype = term), width = 0.2) +
      scale_linetype_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  p
}

plotFEsim2facet <- function (fe, modinf = NULL, level = 0.95, stat = "mean", sd = TRUE,  
                             sigmaScale = NULL, oddsRatio = FALSE, taglv = NA, glv, gcode, ncol, theme) {
  ## plotFEsim2facet, modified from merTools::plotFEsim
  if (!missing(sigmaScale)) {
    fe[, "sd"] <- fe[, "sd"]/sigmaScale
    fe[, stat] <- fe[, stat]/sigmaScale
  }
  fe[, "sd"] <- fe[, "sd"] * qnorm(1 - ((1 - level)/2))
  fe[, "ymax"] <- fe[, stat] + fe[, "sd"]
  fe[, "ymin"] <- fe[, stat] - fe[, "sd"]
  hlineInt <- 0
  if (oddsRatio == TRUE) {
    fe[, "ymax"] <- exp(fe[, "ymax"])
    fe[, stat] <- exp(fe[, stat])
    fe[, "ymin"] <- exp(fe[, "ymin"])
    hlineInt <- 1
  }
  xvar <- "term"
  fe$term <- factor(fe$term, levels = c("Jul","Nov",
                                        "WE","WE:Jul","WE:Nov",
                                        "EA","EA:Jul","EA:Nov",
                                        "NV","NV:Jul","NV:Nov"))
  if (!is.na(taglv)) { 
    tmp <- NULL
    for (k in 1: length(taglv)) {
      tmp <- rbind(tmp, dplyr::filter(fe, tag == taglv[k]))
    }
    fe <- tmp
    fe$tag <- factor(fe$tag, levels = taglv)
  }
  gcol <- rep("black",length(levels(fe$term)))
  gsize <- rep(1,length(levels(fe$term)))
  gltp <- rep(0.5,length(levels(fe$term)))
  for (m in 1: length(glv)) {
    gcol[levels(fe$term) == glv[m]] <- gcode[m]
    gsize[levels(fe$term) == glv[m]] <- 2
  }
  llv <- c("Nov","Jul","EA","NV","WE")
  for (m in 1: length(llv)) {
    gltp[levels(fe$term) == llv[m]] <- 1
  }
  p <- ggplot(aes_string(x = xvar, y = stat, ymax = "ymax", ymin = "ymin"), data = fe) + 
    geom_hline(yintercept = hlineInt, color = I("red")) +
    geom_point(aes(col = term, size = term)) + 
    labs(x = "Group", y = "Fixed Effect") + 
    scale_color_manual("Group", breaks = levels(fe$term), values = gcol) + 
    scale_size_manual("Group", breaks = levels(fe$term), values = gsize) + 
    facet_wrap(~ tag, ncol = ncol, scales = "free") +
    # coord_flip() + 
    theme
  if (sd) {
    p <- p + geom_errorbar(aes(alpha = term), width = 0.2) +
      scale_alpha_manual("Group", breaks = levels(fe$term), values = gltp)
  }
  if(!missing(modinf)) {
    p <- p + geom_text(x = 2, y = 5, label = "abc")
  }
  p
}

multiElementMod <- function(dat, fact, SplMonthlv, grouplv, glv, gcode, 
                            tag = NULL, suffix = "", archiplot = TRUE, log = TRUE) {
  library(merTools);library(lsmeans);library(multcompView);library(car);library(ggplot2);library(lattice)
  fe.g <- NULL; re.g <- NULL; modinf.g <- NULL; shapiro.g <- NULL; posthoc.g <- NULL; modavo.g <- NULL; 
  resid.g <- NULL;
  theme_HL <- theme_bw() + theme(legend.position = "none",axis.text = element_text(angle = 30))
  if(missing(tag)) tag <- unique(dat$elem)
  for (i in 1:length(tag)) {
    
    #shapiro
    shapiroTest(dat = dat, tag = tag[i], grouplv = grouplv, SplMonthlv = SplMonthlv) -> shapiro
    print(shapiro); 
    shapiro.g <- rbind(shapiro.g, cbind(shapiro, tag[i]))
    
    #mod choose
    mod <- modelCompare(dat = dat, tag = tag[i], fact = fact,log = log)
    paste(tag[i],"~",as.character(formula(mod@call))[3]) -> modinf
    print(modinf); 
    modinf.g <- rbind(modinf.g, cbind(modinf, tag[i]))
    
    #anova
    nyma <- anova(mod)
    nymb <- Anova(mod)
    cbind(nyma[,c(1,2,3)],Chisq = nymb[,1], Fvalue = nyma[,4],Pr = nymb[,3], 
          FixedFact = row.names(nyma), tag = tag[i]) -> modavo; print(modavo); 
    modavo.g <- rbind(modavo.g, cbind(modavo, tag[i]))
    
    #posthoc
    modfact <- strsplit(as.character(formula(mod@call))[3],"+", fixed = TRUE)[[1]]
    if (grepl('*', modfact[1], fixed = TRUE)) {
      modfact <- strsplit(modfact[1],"*", fixed = TRUE)[[1]]
    }
    modfxfact <-  modfact[!grepl('(', modfact, fixed = TRUE)]
    modfxfm <- if(length(modfxfact) == 1) {modfxfact[1]} else {paste(modfxfact[1],modfxfact[2],sep = "+")}
    posthoc <-lsmeans::cld(lsmeans(object = mod, adjust = "tukey",
                                   specs = as.formula(paste("pairwise~",modfxfm))),
                           Letters = LETTERS) 
    print(posthoc) 
    if(length(modfxfact) == 1) {
      if (modfxfact == "group ") {
        posthoc <- cbind(posthoc, SplMonth = NA)
      } else {
        if (modfxfact == "SplMonth ") {
          posthoc <- cbind(posthoc, group = NA)
        } else {stop("Kidding me?")}
      }
    }
    posthoc.g <- rbind(posthoc.g,cbind(posthoc, tag[i]))
    
    # fe
    fe <- FEsim(mod)
    fe.g <- rbind(fe.g, cbind(tag = tag[i],fe))
    ggsave(plot = plotFEsim2(fe%>% filter(term != "(Intercept)") %>%
                               mutate(term = gsub("group","",term)) %>% 
                               mutate(term = gsub("SplMonth","",term)),
                             glv = glv, gcode = gcode, theme = theme_HL), 
           paste("ef/plot/Fixeff/",tag[i],".png",sep=""),
           width = 6, height = 4)
    #re
    re <- REsim(mod)
    re.g <- rbind(re.g, cbind(tag = tag[i],re))
    ggsave(plot = merTools::plotREsim(re), 
           paste("ef/plot/Raneff/",tag[i],".png",sep=""),
           width = 6, height = 4)
    #resid
    png(paste("ef/plot/diag/",tag[i],"resid.png",sep=""), 
        width = 30, height = 20, units = "cm", res = 600)
    print(plot(mod, type = c("p", "smooth")))
    dev.off()
    png(paste("ef/plot/diag/",tag[i],"residQQ.png",sep=""), 
        width = 30, height = 20, units = "cm", res = 600)
    print(qqmath(mod, id = 0.05))
    dev.off()
    #png(paste("ef/plot/diag/",tag[i],"_diag_zeta.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using xyplot
    #dev.off()
    #png(paste("ef/plot/diag/",tag[i],"_diag_dens.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using densityplot
    #dev.off()
    #png(paste("ef/plot/diag/",tag[i],"_diag_pair.png",sep=""), 
    #    width = 30, height = 20, units = "cm", res = 600)
    #something using splom
    #dev.off()
    stats::shapiro.test(resid(mod)) -> resid
    resid.g <- rbind(resid.g, 
                     data.frame(tag = tag[i],
                                W = resid$statistic,
                                Pvalue = resid$p.value))
  }
  ## output
  if (log) {write.csv(x = fe.g, file = paste("ef/log/FixedEff",suffix,".csv",sep = ""), row.names = F)
    write.csv(x = re.g, file = paste("ef/log/RandomEff",suffix,".csv",sep = ""), row.names = F)
    write.csv(x = modinf.g, file = paste("ef/log/ModelChoice",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = shapiro.g, file = paste("ef/log/ShapiroRawData",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = posthoc.g, file = paste("ef/log/Posthoc",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = modavo.g, file = paste("ef/log/ModelAnova",suffix,".csv",sep = ""), row.names = F) 
    write.csv(x = resid.g, file = paste("ef/log/ShapiroResid",suffix,".csv",sep = ""), row.names = F)
  }
  "DONE!"
}

denPlot <- function(dat,tag) {
  for(i in 1:length(tag)) {
    ggsave(plot = ggplot(data = dat %>% dplyr::filter(elem == tag[i])) + 
             geom_density(aes(x = resp, fill = group, na.rm = FALSE, stat = "density")) +
             geom_vline(xintercept = c(0.6667,1.5), linetype = I(2), size = I(0.5), col = I("black")) +
             geom_vline(xintercept = c(1), linetype = I(1), size = I(0.75), col = I("black")) +
             facet_grid(group~., scales = "free_y") +
             scale_x_continuous("Enricmment Factor", limits = c(0,2.5)) +
             scale_fill_manual("Location", breaks = c("CL","EA","NV","WE"), 
                               values = c("#B45F04","#31B404","grey50","#013ADF")) +
             theme_bw() + theme(legend.position = "right"),
           filename = paste("ef/plot/den/",tag[i],".png"),
           width = 6, height = 4, dpi = 600)
    
  }
}

## Basic Stat information ----
meanseCal(datareadln())

## LMM fiting and ploting ----
multiElementMod(dat = datareadln() %>% filter(elem != "Al") %>% filter(elem != "orgC"),
                fact = c("group + (1|SiteID)",
                         "SplMonth + (1|SiteID)",
                         "group + (1|SiteID:SplMonth)",
                         "SplMonth + (1|SiteID:SplMonth)",
                         "group + (1|SiteID) + (1|SiteID:SplMonth)",
                         "SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                         "group + SplMonth + (1|SiteID)",
                         "group + SplMonth + (1|SiteID:SplMonth)",
                         "group + SplMonth + (1|SiteID) + (1|SiteID:SplMonth)",
                         "group * SplMonth + (1|SiteID)",
                         "group * SplMonth + (1|SiteID:SplMonth)",
                         "group * SplMonth + (1|SiteID) + (1|SiteID:SplMonth)"),
                SplMonthlv = c("Apr","Jul","Nov"), grouplv = c("EA","CL","WE","NV"), 
                glv = c("EA","WE","NV"), gcode = c("#31B404","#013ADF","grey50"))

## Gather ploting ----
ggsave(plot = plotFEsim2facet(read.csv("ef/log/FixedEff.csv") %>% 
                                filter(term != "(Intercept)") %>% 
                                mutate(term = gsub("group","",term)) %>% 
                                mutate(term = gsub("SplMonth","",term)),
                              glv = c("EA","WE", "NV"), gcode = c("#31B404","#013ADF","grey50"), ncol = 4,
                              theme = theme_bw() + theme(legend.position = "none",
                                                         axis.text = element_text(size= 4,angle = 30),
                                                         axis.title = element_text(size= 6),
                                                         strip.text = element_text(size= 6))), 
       "ef/plot/Fixeff.png",
       width = 6, height = 4, dpi = 600)

## density ploting ----
denPlot(dat = datareadln(), tag = c("N","C","S","N",
                                    "N","N","N","N",
                                    "N","N","N","N"))

  
  