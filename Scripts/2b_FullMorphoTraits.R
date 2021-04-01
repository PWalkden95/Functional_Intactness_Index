rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse)
require(magrittr)
require(vegan)


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

miss_spp <- Jetz_Species$Jetz_Name[which(!(Jetz_Species$Jetz_Name %in% unique(Jetz_Traits$Jetz_Name)))]

mean_Jetz_traits <-PREDICTS_Aves_Am %>% filter(Jetz_Name %in% miss_spp) %>%
  dplyr::select(Jetz_Name,Bill.TotalCulmen,Bill.Depth,Bill.Nares,Bill.Width,
  Wing.Chord,Secondary1,Tail.Length,Tarsus.Length) %>% dplyr::distinct() %>% group_by(Jetz_Name) %>%
  dplyr::mutate(num = n()) %>% ungroup()

colnames(mean_Jetz_traits) <- sub(x = colnames(mean_Jetz_traits), pattern = "\\.", replacement = "_")

Jetz_Traits <- rbind(Jetz_Traits,mean_Jetz_traits)



PCA_Data <- data.frame(Jetz_Name = Jetz_Traits[,1], scale(log(Jetz_Traits[,c(2:9)])))



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

write_rds(file = "Outputs/FullMorphTraits.rds", Full_PC_Scores)
