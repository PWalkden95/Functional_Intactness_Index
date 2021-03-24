rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping

#### Load in the relevant datasets that will be collated to get all relevant information for FII calculations

Forage <- readRDS("../Datasets/GBD/Trophic_Foraging_Niche.rds")
PREDICTS <- readRDS("../Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds")
PREDICTS_Taxonomy <- readRDS("../Datasets/PREDICTS/PREDICTS_AVES_Updated_Taxonomy.rds")
Jetz_Traits <- read.csv("../Datasets/GBD/GBD_Jetz_averages_11_Nov_2020.csv")

colnames(Jetz_Traits)[1] <- "Jetz_Name" 

#### Subset the PREDICTS database to just include the Class Aves resolved to the species level and add in species level
#### Trophic niche and morphometric trait measurements 


PREDICTS_Aves <- PREDICTS %>% filter(Rank == "Species", Class == "Aves") %>% left_join(PREDICTS_Taxonomy[,c("Taxon", "Jetz_Name")], by = "Taxon") %>%
  left_join(Forage[,c(34:46)], by = c("Jetz_Name" = "Taxon")) %>%
  left_join(Jetz_Traits[,c(1,8:18)], by = "Jetz_Name")

### So as we are just going to practise the workflow for FII I will subset just those sites in the Americas


PREDICTS_Aves_Am <- PREDICTS_Aves %>% filter(UN_region == "Americas")

## Lets take a look at the sites in Americas 

table(PREDICTS_Aves_Am$Predominant_habitat, PREDICTS_Aves_Am$Use_intensity)

### so looking at this it would be good to collapse secondary vegetation into a single
### Land Use type and since there was not enough combinations of all Land_use type and 
### intensity create a new factor of LandUse Intensity 

PREDICTS_Aves_Am <- PREDICTS_Aves_Am %>% dplyr::mutate(LandUse = ifelse(grepl(pattern = "secondary", tolower(Predominant_habitat)), "Secondary Vegetation", paste(Predominant_habitat)),
                                                       LandUse = ifelse(grepl(LandUse, pattern = "Primary"), "Primary", paste(LandUse)),
                                                       
                                                       LandUse = ifelse(LandUse == "Cannot decide", NA, paste(LandUse)),
                                                       
                                                       Intensity = ifelse(Use_intensity == "Cannot decide", NA, paste(Use_intensity)),
                                                       
                                                       LandUse_Intensity = ifelse(is.na(Use_intensity), NA, paste(LandUse, Intensity, sep = "_")),
                                                       LandUse_Intensity = ifelse(grepl("NA", LandUse_Intensity), NA, paste(LandUse_Intensity)),
                                                       
                                                       LandUse_Intensity = factor(LandUse_Intensity),
                                                       
                                                       LandUse_Intensity = relevel(LandUse_Intensity, "Primary_Minimal use")) %>%
  dplyr::filter(!is.na(LandUse_Intensity)) %>%
  
  droplevels()


#### Visulise the sites geographically

wm<-map_data("world") %>% filter(region != "Antartica" ) %>% fortify()

## site coords

site_points <- PREDICTS_Aves_Am %>% distinct(SSBS,Longitude,Latitude)

# generate and plot map

site_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+
  geom_point(data = fortify(site_points), aes(Longitude, Latitude),
             colour = "blue", size = 1)+
  theme_classic()

plot(site_plot)


### Correct for sampling effort at each site assuming that measurement increases linearly with sampling effort
### Rescale sampling effort so that its between 0 and 1 by dividing samplin effort at site by the max sampling effort within the study


PREDICTS_Aves_Am <- PREDICTS_Aves_Am %>% 
  
  
  dplyr::group_by(SS) %>% ## group by study 
  
  dplyr::mutate(Max_Sampling_Effort = max(Sampling_effort)) %>%  # max samplin effort in each site
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(Rescaled_Sampling_Effort = Sampling_effort/Max_Sampling_Effort) %>%  ## rescale sampling effort so that the maximum effort with a study is 1
  
  dplyr::filter(Diversity_metric_type =="Abundance" | Diversity_metric_type ==  "Occurrence" | Diversity_metric_type == "effort_corrected_abundance") %>%
  
  dplyr::mutate(Effort_Corrected_Measurement = ifelse(Diversity_metric_is_effort_sensitive == TRUE,
                                                      Measurement/Rescaled_Sampling_Effort,
                                                      Measurement)) %>%
  
  droplevels()

## Now the Biodiversity data is ready to be used 

write_rds(PREDICTS_Aves_Am, "../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")



