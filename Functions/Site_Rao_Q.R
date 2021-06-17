require(tidyverse)
require(magrittr)
require(SYNCSA)



Rao_Q_Func_bias <- function(data, traits){
  
  ### get the list of uncique species across teh whole dataset
  
  Species_abundance <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels()
  
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
  
  spp_traits <- traits %>% filter(Jetz_Name %in% rownames(Species_abundance))
  rownames(spp_traits) <- spp_traits$Jetz_Name
  spp_traits <- spp_traits[,-1]
  
  
  
  Rao_Bias <- rao.diversity(comm = t(Species_abundance),traits =  spp_traits)#THIS is using the package SYNCSA that calcuates Rao's using gowdis 

  return(Rao_Bias)
  }
  ######################################################################
  ########## Here we are also going to calculate an "unbiased" Raos Q###
  ######################################################################
  
Rao_Q_Func_unbias <- function(data,traits){  

  Species_abundance_2 <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels()
  
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  for(site in levels(data$SSBS)){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance_2 <- Species_abundance_2 %>% left_join(Spp_abd[,c("Jetz_Name", "SpeciesSiteAbundance")], by = "Jetz_Name")
    colnames(Species_abundance_2)[which(colnames(Species_abundance_2) == "SpeciesSiteAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance_2) <- Species_abundance_2$Jetz_Name
  Species_abundance_2 <- as.matrix(Species_abundance_2[,-1])
  
  ## Nas to zeros
  
  Species_abundance_2[is.na(Species_abundance_2)] <- 0  
  
  
  comm <- t(Species_abundance_2)
  
  #### species traits 
  
  spp_traits <- traits %>% filter(Jetz_Name %in% rownames(Species_abundance_2))
  rownames(spp_traits) <- spp_traits$Jetz_Name
  spp_traits <- spp_traits[,-1]
  
  
  ##Load in altered SYNCSA function to calculate Rao's Q unbias 
  source("Functions/Rao_Diversity_2.R")
  
  
  Rao_Unbias <- rao_diversity_2(comm = comm, traits = spp_traits)
  Rao_Unbias <- Rao_Unbias$FunRao
  
  
  
  
  
  
  return(Rao_Unbias)
}
