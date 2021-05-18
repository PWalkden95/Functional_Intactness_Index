rm(list = ls())

require(tidyverse)
require(sf)


World_biomes <- sf::read_sf("../Datasets/WWF_Terrestrial_Biomes/official/wwf_terr_ecos.shp")
World_biomes$BIOME <- factor(World_biomes$BIOME)

wm<-map_data("world")

levels(World_biomes$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests",
                                "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests",
                                "Boreal Forests/Taiga","Tropical & Subtropical Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas & Shrublands",
                                "Flooded Grasslands & Savannas","Montane Grasslands & Shrublands","Tundra","Mediterranean Forests, Woodlands & Scrub",
                                "Deserts & Xeric Shrublands","Mangroves","Lakes","Rock & Ice")


colours <- c("darkolivegreen4","darkseagreen","lightgreen","darkolivegreen1","turquoise","darkslategray4","palegoldenrod","lightyellow1","powderblue",
             "cadetblue1","darkslategray3","salmon","khaki1","tomato1","steelblue4","white")



biome_list <- rep(list(NA),length(levels(World_biomes$BIOME)))
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
    geom_sf(data = biome_list[[i]], fill = colours[i], colour = colours[i], alpha = 0.8)
    
}


plot(boundry_plot)


PREDICTS <- readRDS("../Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds") 
  
  
PREDICTS_Aves <- PREDICTS %>% filter(Rank == "Species", Class == "Aves", Measurement > 0) %>% group_by(SSBS) %>% mutate(num_spp = n()) %>% ungroup() 

PREDICTS_Sites <- PREDICTS_Aves %>%
  dplyr::distinct(SSBS,Longitude,Latitude,num_spp)

PREDICTS_Species <- PREDICTS_Aves %>% dplyr::distinct(Country)



boundry_plot <- boundry_plot +
geom_point(data = fortify(PREDICTS_Aves), aes(Longitude, Latitude),
           colour = "black", size = log10(PREDICTS_Aves$num_spp)*5, alpha = I(.1))

ggsave(filename = "../test.png",boundry_plot,device = "png", height = 10, width = 20,dpi=300)

