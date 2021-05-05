rm(list = ls())

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(SYNCSA) ## For calculating functional diversity metric



PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")
traits <- readRDS("Outputs/MeanTraits.rds")


#######################################
#### Functional Diversity of sites ####
#######################################

### RAO's Qaudratic Entropy an abundance-weighted measure of diversity - in our case functional diversity 

## first need to combine a couple of studies 

combine_stud <- PREDICTS_Aves_Am %>% filter(grepl(SS, pattern = "Lasky")) %>% group_by(Site_name, Jetz_Name) %>%
  dplyr::mutate(Effort_Corrected_Measurement = sum(Effort_Corrected_Measurement)) %>% ungroup() %>%
  dplyr::distinct(Site_name, Jetz_Name, .keep_all = TRUE) %>% dplyr::mutate(SSBS = ifelse(SSBS == "HP1_2010__Lasky 2  1", "HP1_2010__Lasky 1  1", paste(SSBS)),
                                                                            SSBS = ifelse(SSBS == "HP1_2010__Lasky 2  2", "HP1_2010__Lasky 1  2", paste(SSBS)),
                                                                            SSBS = ifelse(SSBS == "HP1_2010__Lasky 2  3", "HP1_2010__Lasky 1  3", paste(SSBS)),
                                                                            SSBS = ifelse(SSBS == "HP1_2010__Lasky 2  4", "HP1_2010__Lasky 1  4", paste(SSBS)),
                                                                            SSBS = ifelse(SSBS == "HP1_2010__Lasky 2  5", "HP1_2010__Lasky 1  5", paste(SSBS)),
                                                                            SSBS = factor(SSBS),
                                                                            SS = ifelse(SS == "HP1_2010__Lasky 2", "HP1_2010__Lasky 1", paste(SS)),
                                                                            SS = factor(SS),
                                                                            SSB = ifelse(SSB == "HP1_2010__Lasky 2", "HP1_2010__Lasky 1", paste(SS)),
                                                                            SSB = factor(SSB))

PREDICTS_Aves_Am <- PREDICTS_Aves_Am %>% filter(!(grepl(SS, pattern = "Lasky"))) %>% rbind(combine_stud)


## get abundance data for each species with each site  

abundance_data <- PREDICTS_Aves_Am %>% dplyr::filter(Diversity_metric == "abundance") %>% dplyr::filter(Effort_Corrected_Measurement != 0) %>%
  
  ### filter out studies of just a single species 
  
  dplyr::group_by(SS) %>% dplyr::mutate(study_n_species = n_distinct(Jetz_Name)) %>% dplyr::filter(study_n_species > 1) %>% ungroup() %>%  
  
  #### group by Site and Species get abundance if there are some sites that the same species is recorded multiple times 
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ungroup() %>%   droplevels() %>%  
  
  ### group by just site to get Total site abundance
  group_by(SSBS) %>% dplyr::mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance), site_species = n_distinct(Jetz_Name)) %>% 
  filter(site_species > 1 ) %>%
  
  ungroup() %>%
  
  ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  dplyr::filter(SS != "GN1_2010__Hvenegaard 1") %>%
  
  droplevels()

#### neeed to merge a couple of studies that are replicates of the same sites, eithr in different seasns (wet and dry in Lasky) 
### or point count has been conducted at different radii ( as in hevengaard)



write_rds(abundance_data, file = "Outputs/Rao_abundance_data.rds")



#####################################################
#### Function to calculate Rao's Q for each site ####
####################################################

source("Functions/Site_Rao_Q.R")


PREDICTS_Site_Rao <- data.frame(abundance_data) %>% distinct(SSBS, .keep_all = TRUE) %>% 
  dplyr::select(SS, SSB, SSBS,UN_subregion, LandUse, Use_intensity, LandUse_Intensity, Longitude,Latitude, site_species)

PREDICTS_Site_Rao$SSBS <- as.character(PREDICTS_Site_Rao$SSBS)
PREDICTS_Site_Rao$Bias_Rao <- NA
PREDICTS_Site_Rao$Unbias_Rao <- NA

Rao_data <- Rao_Q_Func(abundance_data,traits[["morpho_traits"]])

PREDICTS_Site_Rao$Bias_Rao <- Rao_data$Bias
PREDICTS_Site_Rao$Unbias_Rao <- Rao_data$Unbias

plot(Rao_data$Unbias ~ Rao_data$Bias)
abline(a=0,b=1)
## save for modelling

write_rds(PREDICTS_Site_Rao, file = "Outputs/PREDICTS_Site_Rao.rds")



########### Might want to try and look for outliers where and how we are going to treat them, it looks like there are a few sites that are dominated 
######### with thousands of the same species which is dragging down the measue of Rao's Q - This would also be a good point to possibly caluclate some
####### other indices of FD e.g FRic, FDis etc etc 

#### Can also calculate a few other metric for FD and will do a little later -- FD, Hvol etc


