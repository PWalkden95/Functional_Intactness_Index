rm(list = ls())

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(foreach)## for parallelising loops
require(doParallel) ## ditto
require(sf)   ### for calculate density of roads
require(rgeos) ## ditto 
require(lwgeom) ### ditto


PREDICTS <- readRDS("../Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds")
PREDICTS_Aves <- PREDICTS %>% filter(Class == "Aves") %>% droplevels()

sites <- levels(PREDICTS_Aves$SSBS)

### load in the roads shape file for the americas
Roads <- st_read("../Datasets/PREDICTS_variables/groads-v1-global-gdb/gROADS_v1.gdb")
### combine into a single shapefile
Roads <- st_combine(Roads)
### transform to be projected on the mercator projection that deals in meters rather than latlong 
Roads <- st_transform(Roads, crs = "+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Both datasets are going to need information on road densities and will be some overlap of sites so make a vector of all unique sites across both datasets


road_densities <- c()

registerDoParallel(cores =10)

memory.limit(1200000)
road_densities <- foreach(site = sites,
                          .combine = "rbind",
                          .packages = c("tidyverse", "sf", "rgeos", "lwgeom")) %dopar%{
                            
                            
                            site_LongLat <- PREDICTS %>% filter(SSBS %in% site) %>% distinct(SSBS,Longitude,Latitude)
                            
                            
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
closeAllConnections()

write_rds(file = "../Datasets/PREDICTS_variables/ave_site_road_densities.rds", road_densities)