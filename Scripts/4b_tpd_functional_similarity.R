rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse)
require(magrittr)
require(vegan)
require(hypervolume)
require(TPD)
require(geosphere)
require(doParallel)
require(foreach)

#### Load in the relevant datasets that will be collated to get all relevant information for FII calculations


Jetz_Traits <- read.csv("../Datasets/GBD/GBD_BiometricsRaw_combined_2020_Dec_07.csv")
PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")


colnames(Jetz_Traits)[6] <- "Jetz_Name"


#### Subset the PREDICTS database to just include the Class Aves resolved to the species level and add in species level
#### Trophic niche and morphometric trait measurements 



Jetz_Species <- PREDICTS_Aves_Am %>% filter(Diversity_metric == "abundance" | Diversity_metric == "effort_corrected_abundance" ) %>% distinct(Jetz_Name)


Jetz_Traits <- Jetz_Traits %>% dplyr::filter(Jetz_Name %in% Jetz_Species$Jetz_Name) %>% 
  dplyr::select(Jetz_Name,Bill_TotalCulmen,Bill_Depth,Bill_Nares,Bill_Width,
                Wing_Chord,Secondary1,Tail_Length,Tarsus_Length) %>% na.omit() %>% group_by(Jetz_Name) %>%
  dplyr::mutate(num = n()) %>% ungroup()
################# Fewer than 4 triats

Jetz_Traits_Am <- Jetz_Traits %>% filter(num >= 4) %>% pull(Jetz_Name) %>% unique()

Jetz_Traits_under <- Jetz_Traits %>% filter(num < 4 & num != 1) %>% pull(Jetz_Name) %>% unique()
single_Jetz_traits <- Jetz_Traits %>% filter(num == 1) %>% pull(Jetz_Name) %>% unique()




Similarity_data <- data.frame(PREDICTS_Aves_Am) %>% 
  
  filter(any(Diversity_metric == "abundance" | Diversity_metric == "effort_corrected_abundance" ) & Effort_Corrected_Measurement > 0) %>%

  filter(Jetz_Name %in% c(Jetz_Traits_Am, Jetz_Traits_under, single_Jetz_traits)) %>%
  
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ### calculating a hypervolume in 3 dimensions with fewer than 21 species on result in inaccurcies 
  group_by(SSBS) %>% dplyr::mutate(Site_spp = n_distinct(Jetz_Name),TotalSiteAbundance = sum(SpeciesSiteAbundance)) %>%
  
  filter(Site_spp > 1) %>% ungroup() %>%
  
  group_by(SS) %>% dplyr::mutate(n_primary_min = sum(LandUse_Intensity == "Primary_Minimal use" ), n_sites_within_studies = n_distinct(SSBS)) %>% ungroup() %>%
  
  filter(n_primary_min > 0 ) %>%  filter(n_sites_within_studies > 1) %>% dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  droplevels() %>%
  
  data.frame()


write_rds(Similarity_data, file = "Outputs/tpd_similarity_data.rds")




registerDoParallel(cores = 5)


studies <- levels(Similarity_data$SS)


TPD_data <- foreach(study = studies,
                        .combine = "rbind",
                        .packages = c("geosphere","tidyverse","TPD", "magrittr"),
                        .inorder = FALSE) %dopar% {



data <- Similarity_data %>% filter(SS == study) %>% droplevels()

### coordinates going to be used to calculate distance between sites

LatLong <- data.frame(data) %>% distinct(SSBS,Latitude,Longitude)

### land use and intensity at each site

LandUse <- data.frame(data) %>% distinct(SSBS, LandUse_Intensity)

### which sites are of primary minimal landuse type and intensity  

Primary_sites <- LandUse %>% filter(LandUse_Intensity == "Primary_Minimal use") %>% pull(SSBS) %>% as.character()

#### all other sites 

all_sites <- LandUse %>% pull(SSBS) %>% as.character()



###### get the comparisons with the primary minimal site and all other sites within the study 

site_comparisons <- expand.grid(Primary_sites, all_sites) %>% 
  dplyr::rename(site1 = Var1, site2 = Var2) %>% dplyr::mutate(site1 = as.character(site1), site2 = as.character(site2)) %>%
  filter(site1 != site2) %>%
  left_join(LatLong, by = c("site1" = "SSBS")) %>% dplyr::rename(site1Lat = Latitude, site1Long = Longitude) %>%
  left_join(LatLong, by = c("site2" = "SSBS")) %>% dplyr::rename(site2Lat = Latitude, site2Long = Longitude) %>%
  left_join(LandUse, by = c("site1" = "SSBS")) %>% dplyr::rename(site1LUI = LandUse_Intensity) %>%
  left_join(LandUse, by = c("site2" = "SSBS")) %>% dplyr::rename(site2LUI = LandUse_Intensity) %>%
  dplyr::left_join(distinct(Similarity_data[,c("SSBS","Site_spp", "Rescaled_Sampling_Effort")]), by = c("site1" = "SSBS")) %>%
  dplyr::rename(site1_spp = Site_spp) %>%
  dplyr::rename(site1_sampling_effort = Rescaled_Sampling_Effort) %>%
  left_join(distinct(Similarity_data[,c("SSBS","Site_spp","Rescaled_Sampling_Effort")]), by = c("site2" = "SSBS")) %>%
  dplyr::rename(site2_spp = Site_spp) %>%
  dplyr::rename(site2_sampling_effort = Rescaled_Sampling_Effort) %>%
  dplyr::mutate(Contrast = paste(site1LUI,site2LUI, sep = "-")) %>%
  filter(site1_sampling_effort == site2_sampling_effort)


  if(nrow(site_comparisons) != 0){



### calculate geographic distance between the sites 

site_one <- as.matrix(site_comparisons[,c("site1Long", "site1Lat")])
site_two <- as.matrix(site_comparisons[,c("site2Long", "site2Lat")])


dist <- data.frame(distHaversine(site_one,site_two))

## Collate together the variables for the dataset

study_data <- site_comparisons %>% dplyr::select(site1,site2, Contrast, site1Long, site1Lat, site2Long,site2Lat, site1_spp, site2_spp)
study_data$min_site_spp <- ifelse(study_data$site1_spp > study_data$site2_spp, study_data$site2_spp, study_data$site1_spp) 
study_data <- cbind(study_data, distance = dist[,1])
study_data$SS <- study
study_data$TPD_Overlap <- NA


for(i in 1:nrow(site_comparisons)){


  ##### get the species in both sites being compared
  site1 <- site_comparisons[i,"site1"]
  site1_spp <- Similarity_data %>% dplyr::filter(SSBS == site1) %>% dplyr::select(Jetz_Name, RelativeAbundance)



  ### site 2

  site2 <- site_comparisons[i,"site2"]
  site2_spp <- Similarity_data %>% dplyr::filter(SSBS == site2) %>% dplyr::select(Jetz_Name, RelativeAbundance)


species <- data.frame(Jetz_Name = unique(c(site1_spp$Jetz_Name,site2_spp$Jetz_Name)))



trait_ranges <- list(c(min(Full_PC_Scores[,2]) -(0.1 * min(Full_PC_Scores[,2])),max(Full_PC_Scores[,2]) + (0.1 * max(Full_PC_Scores[,2]))),
                     c(min(Full_PC_Scores[,3]) -(0.1 * min(Full_PC_Scores[,3])),max(Full_PC_Scores[,3]) + (0.1 * max(Full_PC_Scores[,3]))),
                     c(min(Full_PC_Scores[,4]) -(0.1 * min(Full_PC_Scores[,4])),max(Full_PC_Scores[,4]) + (0.1 * max(Full_PC_Scores[,4]))))


mean_TPD <- c()

if(any(species$Jetz_Name %in% c(single_Jetz_traits, Jetz_Traits_under))){
  under_spp <- species %>% filter(Jetz_Name %in% c(single_Jetz_traits, Jetz_Traits_under))
  under_spp_traits <- under_spp %>% dplyr::left_join(Full_PC_Scores, by = "Jetz_Name")

  spp <- unique(under_spp_traits$Jetz_Name)
  k <- spp[2]
  means <- c()
  sd <- c()
  for(k in spp){
  spp_mean <- data.frame(Jetz_Name = k)
  spp_sd <- data.frame(Jetz_Name = k)
  spp_data <- under_spp_traits %>% filter(Jetz_Name == k)
    if(nrow(spp_data) == 1){
      spp_mean <- spp_data
      spp_sd <- cbind(spp_sd, Foraging.PCA =  c(0.5*sd(Full_PC_Scores[,2])),
                      Loco.PCA =  c(0.5*sd(Full_PC_Scores[,3])),
                      Body.PCA =  c(0.5*sd(Full_PC_Scores[,3])))
    } else{

  for(j in 2:ncol(under_spp_traits)){
  mean <-  mean(under_spp_traits[,j])
  stan <-  sd(under_spp_traits[,j])

  spp_mean <- cbind(spp_mean,mean)
  colnames(spp_mean)[j] <- colnames(under_spp_traits)[j]
  spp_sd <- cbind(spp_sd, stan)
  colnames(spp_sd)[j] <- colnames(under_spp_traits)[j]
    }}
  means <- rbind(means,spp_mean)
  sd <- rbind(sd,spp_sd)
  }


  mean_TPD <- TPDsMean(species = spp, means = means[,c(2:4)], sds = sd[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 61, alpha = 0.95 )
}


Traits <- species %>% dplyr::left_join(Full_PC_Scores, by = "Jetz_Name") %>% filter(Jetz_Name %in% Jetz_Traits_Am)

species$Jetz_Name[which(!(species$Jetz_Name %in% unique(names(Trait_density$TPDs))))]
Trait_density <- TPDs(Traits[,1], Traits[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 61, alpha = 0.95)

if(!is.null(mean_TPD)){
  Trait_density$TPDs <- c(mean_TPD$TPDs, Trait_density$TPDs)
  Trait_density$data$species <- c(Trait_density$data$species, under_spp_traits$Jetz_Name)
  Trait_density$data$traits <- rbind(Trait_density$data$traits, under_spp_traits[,c(2:4)])
}


Community <- species %>% dplyr::left_join(site1_spp, by = "Jetz_Name") %>%
  dplyr::rename(site1_abun = RelativeAbundance) %>%
  dplyr::left_join(site2_spp, by = "Jetz_Name") %>%
  dplyr::rename(site2_abun = RelativeAbundance) %>%
  dplyr::filter(Jetz_Name %in% c(Jetz_Traits_Am, Jetz_Traits_under, single_Jetz_traits))



Community[is.na(Community)] <- 0

rownames(Community) <- Community$Jetz_Name

Community <- Community[,-1]

Overlap <- TPDc(Trait_density, t(Community))


TPDOverlap <- dissim(Overlap)




study_data[i, "TPD_Overlap"] <- 1 - TPDOverlap$communities$dissimilarity[2]


}



study_data <- study_data

}
}

registerDoSEQ()



write_rds(file = "Outputs/tpd_overlap_data.rds", TPD_data)
