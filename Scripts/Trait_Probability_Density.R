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


Jetz_Traits <- read.csv("Datasets/GBD/GBD_BiometricsRaw_combined_2020_Dec_07.csv")
PREDICTS_Aves_Am <- readRDS("Functional_Intactness_Index/Outputs/PREDICTS_Americas_Aves.rds")


colnames(Jetz_Traits)[6] <- "Jetz_Name"

#### Subset the PREDICTS database to just include the Class Aves resolved to the species level and add in species level
#### Trophic niche and morphometric trait measurements 



Jetz_Species <- PREDICTS_Aves_Am %>% filter(Diversity_metric == "abundance" ) %>%distinct(Jetz_Name)


Jetz_Traits_Am <- Jetz_Traits %>% filter(Jetz_Name %in% Jetz_Species$Jetz_Name) %>% 
  dplyr::select(Jetz_Name,Bill_TotalCulmen,Bill_Depth,Bill_Nares,Bill_Width,
                Wing_Chord,Secondary1,Tail_Length,Tarsus_Length) %>% na.omit() %>% group_by(Jetz_Name) %>%
  dplyr::mutate(num = n()) %>% ungroup() %>% filter(num >= 4)


PCA_Data <- data.frame(Jetz_Name = Jetz_Traits_Am[,1], scale(log(Jetz_Traits_Am[,c(2:9)])))



For.pca.data <- PCA_Data[,c(2:5)]

For.pca <- prcomp(For.pca.data, center = TRUE, scale. = TRUE)
For.pca
summary(For.pca)

### PCA on the Locomotory traits - Tarsus, tail and wing dimensions

Loco.pca.data <- PCA_Data[,c(6:9)]

Loco.pca <- prcomp(Loco.pca.data, center = TRUE, scale. = TRUE)
Loco.pca
summary(Loco.pca)


#### Final PCA on the first Principal component of each of the first PCAs to derive an axis of body size 

Body.pca.data <- data.frame(LocoPC1 = Loco.pca$x[,1], ForPC1 = For.pca$x[,1])
Body.pca <- prcomp(Body.pca.data, center = TRUE, scale. = TRUE)

Body.pca
summary(Body.pca)



### Match the independent axes of trait variation to species in PREDICTS 

Full_PC_Scores <- data.frame(Jetz_Name = PCA_Data[,1], Foraging.PCA = For.pca$x[,2], Loco.PCA = Loco.pca$x[,2], Body.PCA = Body.pca$x[,1])


### standardize the PC scores so that the maximum value is 1

for(col in colnames(Full_PC_Scores[,-1])){
  Full_PC_Scores[,col] <- Full_PC_Scores[,col]/max(Full_PC_Scores[,col])
}

write_rds(file = "Functional_Intactness_Index/Outputs/Full_PC_scores.rds", Full_PC_Scores)




Similarity_data <- data.frame(test_2) %>% filter(Diversity_metric == "abundance" & Effort_Corrected_Measurement > 0) %>% filter(Jetz_Name %in% Jetz_Traits_Am$Jetz_Name) %>%
  
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ### calculating a hypervolume in 3 dimensions with fewer than 21 species on result in inaccurcies 
  group_by(SSBS) %>% dplyr::mutate(Site_spp = n_distinct(Jetz_Name),TotalSiteAbundance = sum(SpeciesSiteAbundance)) %>%
  
  filter(Site_spp > 21) %>% ungroup() %>%
  
  group_by(SS) %>% dplyr::mutate(n_primary_min = sum(LandUse_Intensity == "Primary_Minimal use" ), n_sites_within_studies = n_distinct(SSBS)) %>% ungroup() %>%
  
  filter(n_primary_min > 0 ) %>%  filter(n_sites_within_studies > 1) %>% dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  droplevels() %>%
  
  data.frame()


write_rds(Similarity_data, file = "Functional_Intactness_Index/Outputs/abundance_similarity_data.rds")


test <- PREDICTS_Aves_Am %>% filter(SS == "DL1_2011__Latta 2")
test_2 <- PREDICTS_Aves_Am %>% filter(SS == "SE1_2011__Rosselli 1")

remove_dupl_comparisons <- function(data){
  
  data$drop <- NA
  data$comparison <- paste(data$site1,data$site2, sep = "_")
  data[1,"drop"] <- FALSE
  
  
  if(NROW(data) > 1){
  for(i in 2:NROW(data)){
    
    site1 <- data[i,"site1"]
    site2 <- data[i,"site2"]
    
    match_1 <- as.numeric(grep(data[1:i-1,"comparison"], pattern = site1))
    match_2 <- as.numeric(grep(data[1:i-1,"comparison"], pattern = site2))
    
    if(any(match_1 %in% match_2)){
      data[i,"drop"] <- TRUE
    } else {
      data[i,"drop"] <- FALSE
    }
    
  }
    }
  
  if(any(data$drop)){
    data <- data[-which(data$drop == TRUE & data$Contrast == "Primary_Minimal use-Primary_Minimal use"),-which(colnames(data) == "comparison" | colnames(data) == "drop")]
  } else {
    data <- data[,-which(colnames(data) == "comparison" | colnames(data) == "drop")]
  }
  
  
  return(data)
}



registerDoParallel(cores = 5)


studies <- levels(Similarity_data$SS)


Overlap_data <- foreach(study = studies,
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
  dplyr::left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site1" = "SSBS")) %>%
  dplyr::rename(site1_spp = Site_spp) %>%
  left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site2" = "SSBS")) %>%
  dplyr::rename(site2_spp = Site_spp) %>%
  dplyr::mutate(Contrast = paste(site1LUI,site2LUI, sep = "-"))


site_comparisons <- remove_dupl_comparisons(site_comparisons)


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



for(i in 1:NROW(site_comparisons)){
  
  
  ##### get the species in both sites being compared 
  site1 <- site_comparisons[i,"site1"]
  site1_spp <- Similarity_data %>% dplyr::filter(SSBS == site1) %>% dplyr::select(Jetz_Name, RelativeAbundance) 
  
  
  
  ### site 2
  
  site2 <- site_comparisons[i,"site2"]
  site2_spp <- Similarity_data %>% dplyr::filter(SSBS == site2) %>% dplyr::select(Jetz_Name, RelativeAbundance)
  
  
species <-data.frame(Jetz_Name = unique(c(site1_spp$Jetz_Name,site2_spp$Jetz_Name)))

Traits <- species %>% dplyr::left_join(Full_PC_Scores, by = "Jetz_Name")



trait_ranges <- list(c(min(Full_PC_Scores[,2]) -(0.1 * min(Full_PC_Scores[,2])),max(Full_PC_Scores[,2]) + (0.1 * max(Full_PC_Scores[,2]))),
                     c(min(Full_PC_Scores[,3]) -(0.1 * min(Full_PC_Scores[,3])),max(Full_PC_Scores[,3]) + (0.1 * max(Full_PC_Scores[,3]))),
                     c(min(Full_PC_Scores[,4]) -(0.1 * min(Full_PC_Scores[,4])),max(Full_PC_Scores[,4]) + (0.1 * max(Full_PC_Scores[,4]))))


Trait_density <- TPDs(Traits[,1], Traits[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 61, alpha = 0.95)




Community <- species %>% dplyr::left_join(site1_spp, by = "Jetz_Name") %>%
  dplyr::rename(site1_abun = RelativeAbundance) %>%
  dplyr::left_join(site2_spp, by = "Jetz_Name") %>%
  dplyr::rename(site2_abun = RelativeAbundance)



Community[is.na(Community)] <- 0 

rownames(Community) <- Community$Jetz_Name

Community <- Community[,-1]

Overlap <- TPDc(Trait_density, t(Community))



TPDOverlap <- dissim(Overlap)




study_data[i, "TPD_Overlap"] <- 1 - TPDOverlap$communities$dissimilarity[2]


}



study_data <- study_data

}

registerDoSEQ()



hist(Overlap_data$TPD_Overlap)


write_rds(file = "Functional_Intactness_Index/Outputs/Trait_Prob_den.rds", Overlap_data)

table(Overlap_data$Contrast)

##### so a trait density which works is trait ranges = 0.1 and n_divisions =- 61 thank you!