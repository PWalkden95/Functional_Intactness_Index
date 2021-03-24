rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(geosphere) ## calculating geographic distances
require(hypervolume) ## calculating functional similarity/ diversity metrics 
require(foreach)## for parallelising loops
require(doParallel) ## ditto


PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")
traits <- readRDS("Outputs/MeanTraits.rds")


#############################################
#### Functional Similarity between sites ####
#############################################

### Essentially we want to know how much functional overlap there is between Primary minimal sites and sites of other land-use types and intensities.
### This will be calculated by working out the functional beta diversity between sites and how much of this beta-diversity similarity is due to 
### nestedness of functional diversity - as opposed to turnover. This is done by first estimating the multidimensional functional space for the morpho
### -logical and foraging traits for species at each site. -- Actually to get a measure of similarity I think it's just best to get a measure of overlap
### looking at the nestedness is just say what proportion of the dissimilarity is due to the smaller site containing a subset of function of the larger
### site and nestedness is just the replacement of function. However, that can be a post analysis as we are mainly concerned with the similarity for FII


### Join trait values filter out species that are marked as absent from sites 


Similarity_data <- data.frame(PREDICTS_Aves_Am) %>% filter(Effort_Corrected_Measurement > 0) %>%
  
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ### calculating a hypervolume in 3 dimensions with fewer than 21 species on result in inaccurcies 
  group_by(SSBS) %>% dplyr::mutate(Site_spp = n_distinct(Jetz_Name),TotalSiteAbundance = sum(SpeciesSiteAbundance)) %>%
  
  filter(Site_spp > 21) %>% ungroup() %>%
  
  group_by(SS) %>% dplyr::mutate(n_primary_min = sum(LandUse_Intensity == "Primary_Minimal use" ), n_sites_within_studies = n_distinct(SSBS)) %>% ungroup() %>%
  
  filter(n_primary_min > 0 ) %>%  filter(n_sites_within_studies > 1) %>% dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  droplevels() %>%
  
  dplyr::filter(SS != "GN1_2010__Hvenegaard 1") %>%
  
  data.frame()

write_rds(Similarity_data, file = "Outputs/hyper_similarity_data.rds")




########################################################
###### Loop to generate functional overlap of sites ####
########################################################

### With studies with multiple pri-min sites it compares with each other twice -- prim1 v prim2 & prim2 v prim 1 etc etc should only really compare once therefore i have create a function that will remove the duplicated comparisons.

### What studies do we have in the dataset 

studies <- levels(Similarity_data$SS) 

#### empty data frame for results to go into 


registerDoParallel(cores = 4)

Overlap_data <- foreach(study = studies,
                        .combine = "rbind",
                        .packages = c("hypervolume","betapart","tidyverse","magrittr","geosphere","gower","raster")) %dopar% {
                          
                          
                          ### data for study 
                          Overlap_data <- c()
                          
                          
                          
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
                          
                          if(nrow(site_comparisons) != 0 ){
                            
                            
                            
                            
                            
                            ### calculate geographic distance between the sites 
                            
                            site_one <- as.matrix(site_comparisons[,c("site1Long", "site1Lat")])
                            site_two <- as.matrix(site_comparisons[,c("site2Long", "site2Lat")])
                            
                            dist <- data.frame(distHaversine(site_one,site_two))
                            
                            ## Collate together the variables for the dataset
                            
                            study_data <- site_comparisons %>% dplyr::select(site1,site2, Contrast, site1Long, site1Lat, site2Long,site2Lat, site1_spp, site2_spp)
                            study_data$min_site_spp <- ifelse(study_data$site1_spp > study_data$site2_spp, study_data$site2_spp, study_data$site1_spp) 
                            study_data <- cbind(study_data, distance = dist[,1])
                            study_data$SS <- study
                            #study_data$convex_overlap <- NA
                            study_data$hyper_overlap <- NA
                            
                            
                            ## Also want to add in a variable for the minimum number of species in either of the sites used to construct the hypervolumes. This will be used as weights in the models as there may be greater uncertainty in the hypervolume overlaps when fewer species have been recorded at either site.
                            

                            for(i in 1:NROW(site_comparisons)){


                              ##### get the species in both sites being compared

                              site1 <- site_comparisons[i,"site1"]
                              site1_spp <- Similarity_data %>% filter(SSBS == site1) %>% distinct(Jetz_Name, .keep_all = FALSE)


                              ### join the traits and drop species names

                              site1_data <- site1_spp %>% left_join(traits[["morpho_traits"]], by = "Jetz_Name")
                              rownames(site1_data) <- site1_data$Jetz_Name
                              site1_data <- as.matrix(site1_data[,-1])


                              ### calculate support vector machine and minimum convex hull hypervolumes

                              hypersvm_1 <- hypervolume(site1_data, method = "svm")
                              #convex_1 <- expectation_convex(site1_data, check.memory = FALSE)


                              ### site 2

                              site2 <- site_comparisons[i,"site2"]
                              site2_spp <- Similarity_data %>% filter(SSBS == site2) %>% distinct(Jetz_Name, .keep_all = FALSE)

                              site2_data <- site2_spp %>% left_join(xz, by = "Jetz_Name")
                              rownames(site2_data) <- site2_data$Jetz_Name
                              site2_data <- as.matrix(site2_data[,-1])

                              hypersvm_2 <- hypervolume(site2_data, method = "svm")
                              #convex_2 <- expectation_convex(site2_data, check.memory = FALSE)

                              svm_list <- hypervolume_set(hypersvm_1,hypersvm_2, check.memory = FALSE)
                              #convex_list <- hypervolume_set(convex_1, convex_2, check.memory = FALSE)

                              svm_overlap <- hypervolume_overlap_statistics(svm_list)
                              #convex_overlap <- hypervolume_overlap_statistics(convex_list)

                              study_data[i,"hyper_overlap"] <- svm_overlap[1]
                              #study_data[i,"convex_overlap"] <- convex_overlap[1]

                            }
                            return(study_data)
                          }
                          
                        }

table(Overlap_data$Contrast)


registerDoSEQ()

write_rds(file = "Outputs/hyper_overlap_data.rds", Overlap_data)