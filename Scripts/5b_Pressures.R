rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(gower) ## calculating gower distances
require(raster) ## for loading in rasters of environmental variables and roads
require(geosphere) ## calculating geographic distances



###########################
###Human Population Density
###########################

PREDICTS_Site_Rao <- readRDS("Outputs/PREDICTS_Site_Rao.rds")
hyper_data <- readRDS("Outputs/hyper_overlap_data.rds")
tpd_data <- readRDS("Outputs/tpd_overlap_data.rds")
PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")

hpd <- raster("../Datasets/PREDICTS_variables/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")


### calculate human population density for Functional diversity sites 

hpd_values <- raster::extract(hpd,PREDICTS_Site_Rao[,c("Longitude","Latitude")])

## human population density trandformed with the log +1 transformation
PREDICTS_Site_Rao$logHPD <- log(hpd_values + 1) 
PREDICTS_Site_Rao <- PREDICTS_Site_Rao %>% group_by(SS) %>% dplyr::mutate(CNTRLlogHPD = mean(logHPD)) %>% ungroup()

####################################################
### calculate for functional similarity/overlap data
#####################################################

#### first for hypervolume


control_hpd <- PREDICTS_Aves_Am %>% dplyr::distinct(SS,SSBS,Longitude,Latitude)
hpd_values_control <- raster::extract(hpd,control_hpd[,c("Longitude","Latitude")])
control_hpd$HPD <- log(hpd_values_control + 1)
control_hpd <- control_hpd  %>% group_by(SS) %>% dplyr::mutate(CNTRLlogHPD = mean(HPD))

hpd_values_site1 <- raster::extract(hpd,hyper_data[,c("site1Long","site1Lat")])
hpd_values_site2 <- raster::extract(hpd,hyper_data[,c("site2Long","site2Lat")])

hyper_data$S1logHPD <- log(hpd_values_site1 + 1)
hyper_data$S2logHPD <- log(hpd_values_site2 + 1)

hyper_data$logHPDdiff <- hyper_data$S2logHPD - hyper_data$S1logHPD
hyper_data <- hyper_data %>% dplyr::left_join(unique(control_hpd[,c("SS","CNTRLlogHPD")], by = "SS"))


##### then for tpd



hpd_values_site1 <- raster::extract(hpd,tpd_data[,c("site1Long","site1Lat")])
hpd_values_site2 <- raster::extract(hpd,tpd_data[,c("site2Long","site2Lat")])

tpd_data$S1logHPD <- log(hpd_values_site1 + 1)
tpd_data$S2logHPD <- log(hpd_values_site2 + 1)

tpd_data$logHPDdiff <- tpd_data$S2logHPD - tpd_data$S1logHPD
tpd_data <- tpd_data %>% dplyr::left_join(unique(control_hpd[,c("SS","CNTRLlogHPD")], by = "SS"))




########################
#Environmental distance
########################


Bioclim_5 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_5.tif")
Bioclim_6 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_6.tif")
Bioclim_13 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_13.tif")
Bioclim_14 <- raster("../Datasets/Environmental_Variables/wc2.1_30s_bio_14.tif")
Elevation <- raster("../Datasets/Environmental_Variables/wc2.1_30s_elev.tif")


### we only have to calculate the environmental distance when comparing two sites.

site_one <- as.matrix(hyper_data[,c("site1Long","site1Lat")])
site_two <- as.matrix(hyper_data[,c("site2Long","site2Lat")])

environ_1 <- data.frame(Ele  = raster::extract(Elevation,site_one), B5 = raster::extract(Bioclim_5,site_one), B6 = raster::extract(Bioclim_6, site_one),
                        B13 = raster::extract(Bioclim_13, site_one), B14 = raster::extract(Bioclim_14, site_one))

environ_2 <- data.frame(Ele  = raster::extract(Elevation,site_two), B5 = raster::extract(Bioclim_5,site_two), B6 = raster::extract(Bioclim_6, site_two),
                        B13 = raster::extract(Bioclim_13, site_two), B14 = raster::extract(Bioclim_14, site_two))


## calculate the gowers distance between environmental variables 

environ_dist <- gower_dist(x = environ_1, y = environ_2)
#### I dont know how to best transform it now so Im just going to keep it as is

hyper_data$env_distance <- environ_dist


#### then for tpd


site_one <- as.matrix(tpd_data[,c("site1Long","site1Lat")])
site_two <- as.matrix(tpd_data[,c("site2Long","site2Lat")])

environ_1 <- data.frame(Ele  = raster::extract(Elevation,site_one), B5 = raster::extract(Bioclim_5,site_one), B6 = raster::extract(Bioclim_6, site_one),
                        B13 = raster::extract(Bioclim_13, site_one), B14 = raster::extract(Bioclim_14, site_one))

environ_2 <- data.frame(Ele  = raster::extract(Elevation,site_two), B5 = raster::extract(Bioclim_5,site_two), B6 = raster::extract(Bioclim_6, site_two),
                        B13 = raster::extract(Bioclim_13, site_two), B14 = raster::extract(Bioclim_14, site_two))


## calculate the gowers distance between environmental variables 

environ_dist <- gower_dist(x = environ_1, y = environ_2)
#### I dont know how to best transform it now so Im just going to keep it as is

tpd_data$env_distance <- environ_dist




########################
## Density of Roads
########################

### Load Road Desnities

road_densities <- readRDS("../Datasets/PREDICTS_variables/ave_site_road_densities.rds")

## simply join road densities to functional diversity dataset

PREDICTS_Site_Rao <- PREDICTS_Site_Rao %>% dplyr::left_join(road_densities, by = c("SSBS" = "site"))

### join to similarity dataset and calculate difference in density of roads

hyper_data <- hyper_data %>% 
  dplyr::left_join(road_densities, by = c("site1" = "site")) %>%
  dplyr::rename(S1RD1K = density_1km, S1RD50K = density_50km) %>%
  dplyr::left_join(road_densities, by = c("site2" = "site")) %>%
  dplyr::rename(S2RD1K = density_1km, S2RD50K = density_50km)

tpd_data <- tpd_data %>% 
  dplyr::left_join(road_densities, by = c("site1" = "site")) %>%
  dplyr::rename(S1RD1K = density_1km, S1RD50K = density_50km) %>%
  dplyr::left_join(road_densities, by = c("site2" = "site")) %>%
  dplyr::rename(S2RD1K = density_1km, S2RD50K = density_50km)



write_rds(file = "Outputs/hyper_overlap_data.rds", hyper_data)
write_rds(file = "Outputs/tpd_overlap_data.rds", tpd_data)
write_rds(file = "Outputs/PREDICTS_Site_Rao.rds", PREDICTS_Site_Rao)


##############
##############
#############
