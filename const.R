source("uniTls_presetPaths.R");source("uniTls_pkgInstall.R");

pkgLoad("ggplot2")

# Theme
themeset <- 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        legend.position = "none",
        panel.grid = element_blank())
  
  