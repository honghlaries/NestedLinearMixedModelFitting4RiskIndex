# Creating dataset for shareing

# 

# data creating
## Element content
dat <- read.csv("Data/Result_Sediment.csv") %>%
  dplyr::inner_join(read.csv("Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
  dplyr::inner_join(read.csv("Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
  dplyr::inner_join(read.csv("Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
  dplyr::select(group,SiteID,SplMonth,N:Pb)
write.csv(dat%>%select(group,SiteID,SplMonth,N:S,Cr:Pb),"C:/Users/hongh/Desktop/tmp/ElementContent.csv")

## Geo-accumulation Index
dat <- read.csv("./Data/Result_Sediment.csv") %>%
  dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
  dplyr::inner_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
  dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
  dplyr::select(group,SiteID,SplMonth,N_Igeo:Pb_Igeo)
colnames(dat) <- gsub("_Igeo","",colnames(dat),fixed=T)
write.csv(dat%>%select(group,SiteID,SplMonth,N:Pb),"C:/Users/hongh/Desktop/tmp/Igeo.csv")

## Enrichment factor 
dat <- read.csv("./Data/Result_Sediment.csv") %>%
    dplyr::inner_join(read.csv("./Data/meta_SedimentSampleList.csv"), by = c("SplNo" = "SplNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_Quadrat.csv"), by = c("QudNo" = "QudNo")) %>%
    dplyr::inner_join(read.csv("./Data/meta_SiteGroup.csv"), by = c("SiteID" = "SiteID")) %>% 
    dplyr::select(group,SiteID,SplMonth,N_EnrichmentFator:Pb_EnrichmentFator)
colnames(dat) <- gsub("_EnrichmentFator","",colnames(dat),fixed=T)
write.csv(dat%>%select(group,SiteID,SplMonth,N:S,Cr:Pb),"C:/Users/hongh/Desktop/tmp/EF.csv")