rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(SYNCSA) ## For calculating functional diversity metric
require(vegan) ## for performing ordination (PCA/PCoA)
require(gower) ## calculating gower distances
require(raster) ## for loading in rasters of environmental variables and roads
require(geosphere) ## calculating geographic distances
require(ade4)  ### i can't remember right now
require(FD) ## calculating functional diversity metrics
require(hypervolume) ## calculating functional similarity/ diversity metrics 
require(betapart) ## for calculating functional similarity and for decomposing functional diversity into turnover and nestedness components 
require(foreach)## for parallelising loops
require(doParallel) ## ditto
require(sf)   ### for calculate density of roads
require(rgeos) ## ditto 
require(lwgeom) ### ditto


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
  
  filter(Site_spp > 1) %>% ungroup() %>%
  
  group_by(SS) %>% dplyr::mutate(n_primary_min = sum(LandUse_Intensity == "Primary_Minimal use" ), n_sites_within_studies = n_distinct(SSBS)) %>% ungroup() %>%
  
  filter(n_primary_min > 0 ) %>%  filter(n_sites_within_studies > 1) %>% dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  droplevels() %>%
  
  dplyr::filter(SS != "GN1_2010__Hvenegaard 1") %>%
  
  data.frame()

write_rds(Similarity_data, file = "Outputs/similarity_data.rds")




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
                          
                          
                          ### Also want to add in a variable for the minimum number of species in either of the sites used to construct the hypervolumes. This will be used as weights in the models as there may be greater uncertainty in the hypervolume overlaps when fewer species have been recorded at either site.
                          
                          
                          # for(i in 1:NROW(site_comparisons)){
                          #   
                          #   
                          #   ##### get the species in both sites being compared 
                          #   
                          #   site1 <- site_comparisons[i,"site1"]
                          #   site1_spp <- Similarity_data %>% filter(SSBS == site1) %>% distinct(Jetz_Name, .keep_all = FALSE)
                          #   
                          #   
                          #   ### join the traits and drop species names
                          #   
                          #   site1_data <- site1_spp %>% left_join(PC_Scores, by = "Jetz_Name") 
                          #   rownames(site1_data) <- site1_data$Jetz_Name
                          #   site1_data <- as.matrix(site1_data[,-1])
                          #   
                          #   
                          #   ### calculate support vector machine and minimum convex hull hypervolumes 
                          #   
                          #   hypersvm_1 <- hypervolume(site1_data, method = "svm")
                          #   #convex_1 <- expectation_convex(site1_data, check.memory = FALSE)
                          #   
                          #   
                          #   ### site 2
                          #   
                          #   site2 <- site_comparisons[i,"site2"]
                          #   site2_spp <- Similarity_data %>% filter(SSBS == site2) %>% distinct(Jetz_Name, .keep_all = FALSE)
                          #   
                          #   site2_data <- site2_spp %>% left_join(PC_Scores, by = "Jetz_Name") 
                          #   rownames(site2_data) <- site2_data$Jetz_Name
                          #   site2_data <- as.matrix(site2_data[,-1])
                          #   
                          #   hypersvm_2 <- hypervolume(site2_data, method = "svm")
                          #   #convex_2 <- expectation_convex(site2_data, check.memory = FALSE)
                          #   
                          #   svm_list <- hypervolume_set(hypersvm_1,hypersvm_2, check.memory = FALSE)
                          #   #convex_list <- hypervolume_set(convex_1, convex_2, check.memory = FALSE)
                          #   
                          #   svm_overlap <- hypervolume_overlap_statistics(svm_list)
                          #   #convex_overlap <- hypervolume_overlap_statistics(convex_list)
                          #   
                          #   study_data[i,"hyper_overlap"] <- svm_overlap[1]
                          #   #study_data[i,"convex_overlap"] <- convex_overlap[1]
                          #   
                          # }
                          
                          return(study_data)
                          }
                          
                        }

table(Overlap_data$Contrast)


registerDoSEQ()

write_rds(file = "Outputs/Functional_Overlap_data.rds", Overlap_data)



###########################
###Human Population Density
###########################

rm(list = ls())

PREDICTS_Site_Rao <- readRDS("Outputs/PREDICTS_Site_Rao.rds")
Overlap_data <- readRDS("Outputs/Functional_Overlap_data.rds")



hpd <- raster("../Datasets/PREDICTS_variables/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")


### calculate human population density for Functional diversity sites 

hpd_values <- raster::extract(hpd,PREDICTS_Site_Rao[,c("Longitude","Latitude")])

## human population density trandformed with the log +1 transformation
PREDICTS_Site_Rao$logHPD <- log(hpd_values + 1) 
PREDICTS_Site_Rao <- PREDICTS_Site_Rao %>% group_by(SS) %>% dplyr::mutate(CNTRLlogHPD = mean(logHPD)) %>% ungroup()


### calculate for functional similarity/overlap data
control_hpd <- PREDICTS_Aves_Am %>% dplyr::distinct(SS,SSBS,Longitude,Latitude)
hpd_values_control <- raster::extract(hpd,control_hpd[,c("Longitude","Latitude")])
control_hpd$HPD <- log(hpd_values_control + 1)
control_hpd <- control_hpd  %>% group_by(SS) %>% dplyr::mutate(CNTRLlogHPD = mean(HPD))

hpd_values_site1 <- raster::extract(hpd,Overlap_data[,c("site1Long","site1Lat")])
hpd_values_site2 <- raster::extract(hpd,Overlap_data[,c("site2Long","site2Lat")])

Overlap_data$S1logHPD <- log(hpd_values_site1 + 1)
Overlap_data$S2logHPD <- log(hpd_values_site2 + 1)

Overlap_data$logHPDdiff <- Overlap_data$S2logHPD - Overlap_data$S1logHPD
Overlap_data <- Overlap_data %>% dplyr::left_join(unique(control_hpd[,c("SS","CNTRLlogHPD")], by = "SS"))

########################
#Environmental distance
########################


Bioclim_5 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_5.tif")
Bioclim_6 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_6.tif")
Bioclim_13 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_13.tif")
Bioclim_14 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_14.tif")
Elevation <- raster("../Datasets/Environmental_Variables/wc2.1_30s_elev.tif")


### we only have to calculate the environmental distance when comparing two sites.

site_one <- as.matrix(Overlap_data[,c("site1Long","site1Lat")])
site_two <- as.matrix(Overlap_data[,c("site2Long","site2Lat")])

environ_1 <- data.frame(Ele  = raster::extract(Elevation,site_one), B5 = raster::extract(Bioclim_5,site_one), B6 = raster::extract(Bioclim_6, site_one),
                        B13 = raster::extract(Bioclim_13, site_one), B14 = raster::extract(Bioclim_14, site_one))

environ_2 <- data.frame(Ele  = raster::extract(Elevation,site_two), B5 = raster::extract(Bioclim_5,site_two), B6 = raster::extract(Bioclim_6, site_two),
                        B13 = raster::extract(Bioclim_13, site_two), B14 = raster::extract(Bioclim_14, site_two))


## calculate the gowers distance between environmental variables 

environ_dist <- gower_dist(x = environ_1, y = environ_2)
#### I dont know how to best transform it now so Im just going to keep it as is

Overlap_data$env_distance <- environ_dist




########################
## Density of Roads
########################

### load in the roads shape file for the americas
Roads <- st_read("../Datasets/PREDICTS_variables/groads-v1-americas-shp/groads-v1-americas-shp/gROADS-v1-americas.shp")
### combine into a single shapefile
Roads <- st_combine(Roads)
### transform to be projected on the mercator projection that deals in meters rather than latlong 
Roads <- st_transform(Roads, crs = "+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Both datasets are going to need information on road densities and will be some overlap of sites so make a vector of all unique sites across both datasets

sites <- unique(c(Overlap_data$site1,Overlap_data$site2, PREDICTS_Site_Rao$SSBS))

road_densities <- c()

registerDoParallel(cores = 4)


road_densities <- foreach(site = sites,
                          .combine = "rbind",
                          .packages = c("tidyverse", "sf", "rgeos", "lwgeom")) %dopar%{
                            
                            
                            site_LongLat <- PREDICTS_Aves_Am %>% filter(SSBS %in% site) %>% distinct(SSBS,Longitude,Latitude)
                            
                            
                            point <- as.matrix(site_LongLat[,c("Longitude","Latitude")])
                            
                            point <- st_point(point)
                            point <- st_sfc(point)
                            
                            ### first crs needs to be in longlat format as the points are coordinates then transformed into a meters based projection such as mercator to calculate distance in km
                            
                            st_crs(point) <- "+proj=longlat +datum=WGS84 +no_defs"
                            point <- st_transform(point, crs = "+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
                            
                            ## 1k and 50k radius around points
                            
                            buffer_1km <- st_buffer(point, 1000)
                            buffer_50km <- st_buffer(point, 50000)
                            
                            
                            ## overlap buffer with roads shape file to just get roads within each radius
                            
                            intersect_1km <- st_intersection(Roads, buffer_1km)
                            
                            #### if there are any roads caluclate the density as length of roads/area of radius
                            
                            if(length(intersect_1km) != 0 ){
                              density_1km <- (st_length(intersect_1km)/1000)/(st_area(buffer_1km)/1000000)
                            } else {
                              density_1km <- 0
                            }
                            
                            intersect_50km <- st_intersection(Roads, buffer_50km)
                            
                            if(length(intersect_50km) != 0) {
                              density_50km <- (st_length(intersect_50km)/1000)/(st_area(buffer_50km)/1000000)
                            } else {
                              density_50km <- 0
                            }
                            
                            
                            densities <- data.frame(site = paste(site), density_1km = as.numeric(density_1km), density_50km = as.numeric(density_50km))
                            
                            
                          }

registerDoSEQ()

## simply join road densities to functional diversity dataset

PREDICTS_Site_Rao <- PREDICTS_Site_Rao %>% dplyr::left_join(road_densities, by = c("SSBS" = "site"))

### join to similarity dataset and calculate difference in density of roads

Overlap_data <- Overlap_data %>% 
  dplyr::left_join(road_densities, by = c("site1" = "site")) %>%
  dplyr::rename(S1RD1K = density_1km, S1RD50K = density_50km) %>%
  dplyr::left_join(road_densities, by = c("site2" = "site")) %>%
  dplyr::rename(S2RD1K = density_1km, S2RD50K = density_50km)


write_rds(file = "../Datasets/PREDICTS_variables/Road_densities1_and_50k.rds", road_densities)




write_rds(file = "Outputs/Functional_Overlap_data.rds", Overlap_data)
write_rds(file = "Outputs/PREDICTS_Site_Rao.rds", PREDICTS_Site_Rao)


##############
##############
#############
