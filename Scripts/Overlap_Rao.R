rm(list = ls())

require(tidyverse)
require(TPD)


Full_PC_Scores <- read_rds("Functional_Intactness_Index/Full_PC_scores.rds") 
abundance_data <- read_rds("Functional_Intactness_Index/abundance_data.rds")

data <- abundance_data 

Species_abundance <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% filter(Jetz_Name %in% Full_PC_Scores$Jetz_Name) %>% droplevels()

### loop over every site in the dataset to collate the relative abundance of each species

for(site in levels(data$SSBS)){
  
  Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
  
  
  ### Join species withining site to dataframe and rename column as site name
  
  Species_abundance <- Species_abundance %>% left_join(Spp_abd[,c("Jetz_Name", "RelativeAbundance")], by = "Jetz_Name")
  colnames(Species_abundance)[which(colnames(Species_abundance) == "RelativeAbundance")] <- paste(site)
}

## rename rows as species and drop column from dataset 

rownames(Species_abundance) <- Species_abundance$Jetz_Name
Species_abundance <- as.matrix(Species_abundance[,-1])

## Nas to zeros

Species_abundance[is.na(Species_abundance)] <- 0  

### Join all species in datasets traits scores

traits <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels() %>% left_join(Full_PC_Scores, by = "Jetz_Name") %>% na.omit

trait_ranges <- list(c(min(Full_PC_Scores[,2]) -(0.3 * min(Full_PC_Scores[,2])),max(Full_PC_Scores[,2]) + (0.3 * max(Full_PC_Scores[,2]))),
                     c(min(Full_PC_Scores[,3]) -(0.3 * min(Full_PC_Scores[,3])),max(Full_PC_Scores[,3]) + (0.3 * max(Full_PC_Scores[,3]))),
                     c(min(Full_PC_Scores[,4]) -(0.3 * min(Full_PC_Scores[,4])),max(Full_PC_Scores[,4]) + (0.3 * max(Full_PC_Scores[,4]))))


TPD_spp <- TPD::TPDs(traits[,1], traits[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 66, alpha = 0.95)



comm <- TPDc(TPD_spp, t(Species_abundance))
dissim <- dissim(TPD_spp)


Rao <- TPD::Rao(diss = dissim, TPDc = comm)

hist(Rao$alpha_eqv)


write_rds(file = "Functional_Intactness_Index/Overlap_Rao.rds", Rao)


dissimilarity <- dissim$populations$dissimilarity
