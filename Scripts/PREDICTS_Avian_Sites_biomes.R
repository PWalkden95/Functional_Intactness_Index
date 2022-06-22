rm(list = ls())

require(tidyverse)
require(sf)


World_biomes <- sf::read_sf("../../Datasets/WWF_Terrestrial_Biomes/official/wwf_terr_ecos.shp")
World_biomes$BIOME <- factor(World_biomes$BIOME)

realm_data <- World_biomes %>% dplyr::filter(REALM %in% c("NT","PA","NA","IM","AT","AA"))
realm_data$REALM <- factor(realm_data$REALM)

levels(realm_data$REALM) <- c("Australasia","Afrotropic","Indo-Malay","Nearctic","Neotropic","Palearctic") 

realm_colours <- c("khaki4","springgreen4","brown4", "dodgerblue3","olivedrab4","steelblue4")


wm<-map_data("world")


realm_list <- list()

for(realm in levels(realm_data$REALM)){
  data <- realm_data %>% dplyr::filter(REALM == realm)
  
  polygon <- st_combine(data)
  
  realm_list[[realm]] <- polygon
}


boundry_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  lims(x = c(-155,175),  ### to capture the points
       y = c(-55,80))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line( colour = "white"),
        panel.grid.minor = element_line(colour = "white"))


for(i in 1:length(realm_list)){
  boundry_plot <- boundry_plot +
    geom_sf(data = realm_list[[i]], fill = realm_colours[i], colour = realm_colours[i])
  
}


###############################################
#################################################
#################################################

levels(World_biomes$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests",
                                "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
                                "Boreal Forests/Taiga","Tropical & Subtropical Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas & Shrublands",
                                "Flooded Grasslands & Savannas","Montane Grasslands & Shrublands","Tundra","Mediterranean Forests, Woodlands & Scrub",
                                "Deserts & Xeric Shrublands","Mangroves","Lakes","Rock & Ice")


colours <- c("darkolivegreen4","darkseagreen","lightgreen","darkolivegreen1","turquoise","darkslategray4","darkseagreen2","khaki1","powderblue",
             "cadetblue1","darkslategray3","salmon","lightyellow1","tomato1","steelblue4","white")


biome_list <- c(rep(list(NA),length(levels(World_biomes$BIOME))))


i <- 1

for(biome in levels(World_biomes$BIOME)){
  

  
  data <- World_biomes %>% dplyr::filter(BIOME == biome)
  
  polygon <- st_combine(data)

  biome_list[[i]] <- polygon
  names(biome_list)[i] <- biome
  
  i <- i + 1
  }



boundry_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  lims(x = c(-180,180),  ### to capture the points
       y = c(-90,90))+
    theme_classic() +
  theme(panel.background = element_rect(fill = "skyblue3", colour = "skyblue4"))


for(i in 1:length(biome_list)){
  boundry_plot <- boundry_plot +
    geom_sf(data = biome_list[[i]], fill = colours[i], colour = colours[i])
    
}

plot(boundry_plot)

############################################################
############################################################
############################################################


PREDICTS_full <- readRDS("Outputs/refined_predicts.rds")

PREDICTS <- PREDICTS_full %>%  group_by(SSBS) %>% dplyr::mutate(num_spp = n()) %>% ungroup()


PREDICTS_Sites <- PREDICTS %>%
  dplyr::distinct(SSBS,Longitude,Latitude,num_spp)




boundry_plot <- boundry_plot +
geom_point(data = fortify(PREDICTS_Sites), aes(Longitude, Latitude),
           colour = "black", size = 10, alpha = I(.5))


plot(boundry_plot)

ggsave(filename = "C:/Users/patri/Desktop/map_plot.png",boundry_plot,device = "png", height = 10.08, width = 18.46,dpi=500)

