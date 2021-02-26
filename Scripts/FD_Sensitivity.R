rm(list = ls())

require(tidyverse)
require(magrittr)
require(gtools)
require(SYNCSA)


PREDICTS_Aves <- readRDS("Functional_Intactness_Index/PREDICTS_Americas_Aves.rds")


PCA_Data <- data.frame(PREDICTS_Aves[,c(90,103:113)])
PCA_Data <- distinct(PCA_Data, Jetz_Name, .keep_all = TRUE)


#### Log transform - then standardise and centre for a mean of zero and a standard deviation of 1

for(col in colnames(PCA_Data)[c(2:12)]){
  PCA_Data[,col] <- log(PCA_Data[,col] + 1)
  PCA_Data[,col] <- (PCA_Data[,col] - mean(PCA_Data[,col])) / sd(PCA_Data[,col])
}

### The different combinations of the full traits dropping a single trait out then recalculating the PCA

Full.pca.data <- PCA_Data[,c(2:6,8:9,11)]

Full_Combinations <- combinations(NCOL(Full.pca.data), (NCOL(Full.pca.data) - 1))


Full.pca.list <-c() 
for(i in 1:NROW(Full_Combinations)){
  PCA <- prcomp(Full.pca.data[,c(Full_Combinations[i,])], , scale. = TRUE, center = TRUE)
  Full.pca.list <- c(Full.pca.list, list(PCA))
  names(Full.pca.list)[i] <- paste("Combination_", i, sep = "")
}


## Do the same for the Two-step PCA


#### Loco combinations

Loco.pca.data <- PCA_Data[,c(6,8,9,11)]

Loco_Combinations <- combinations(NCOL(Loco.pca.data), (NCOL(Loco.pca.data) - 1))


Loco.pca.list <-c() 
for(i in 1:NROW(Loco_Combinations)){
  PCA <- prcomp(Loco.pca.data[,c(Loco_Combinations[i,])], , scale. = TRUE, center = TRUE)
  Loco.pca.list <- c(Loco.pca.list, list(PCA))
  names(Loco.pca.list)[i] <- paste("Loco_Combination_", i, sep = "")
}


Loco.full.PCA <- prcomp(Loco.pca.data, scale. = TRUE, center = TRUE)
Loco.pca.list <- c(Loco.pca.list, list(Loco.full.PCA))
names(Loco.pca.list)[NROW(Loco_Combinations) + 1] <- "Loco_Full"


### Foraging combinations

For.pca.data <- PCA_Data[,c(2:5)]

For_Combinations <- combinations(NCOL(For.pca.data), (NCOL(For.pca.data) - 1))


For.pca.list <-c() 
for(i in 1:NROW(For_Combinations)){
  PCA <- prcomp(For.pca.data[,c(For_Combinations[i,])], , scale. = TRUE, center = TRUE)
  For.pca.list <- c(For.pca.list, list(PCA))
  names(For.pca.list)[i] <- paste("For_Combination_", i, sep = "")
}

For.full.PCA <- prcomp(For.pca.data, scale. = TRUE, center = TRUE)
For.pca.list <- c(For.pca.list, list(For.full.PCA))
names(For.pca.list)[NROW(For_Combinations) + 1] <- "For_Full"

#### Get the two-step permutations of the Loco and Foraging 


Two_Step_Combs <- combinations(n = length(names(c(For.pca.list,Loco.pca.list))),r = 2, v= names(c(For.pca.list,Loco.pca.list)))
Two_Step_Combs <- Two_Step_Combs[-c(which(substr(Two_Step_Combs[,1],1,4) == substr(Two_Step_Combs[,2],1,4))),]

Two_step_list <- c()
for(i in 1:NROW(Two_Step_Combs)){
  For.pca2 <- For.pca.list[[Two_Step_Combs[i,grep(Two_Step_Combs[i,], pattern = "For")]]]
  Loco.pca2 <- Loco.pca.list[[Two_Step_Combs[i,grep(Two_Step_Combs[i,], pattern = "Loco")]]]

  Body.pca.data2 <- data.frame(ForagePC1 = For.pca2$x[,1], LocoPC1 = Loco.pca2$x[,1])
  Body.pca2 <- prcomp(Body.pca.data2, scale. = TRUE, center = TRUE)
  
  two_step_PCdata <- data.frame(Jetz_Name = PCA_Data$Jetz_Name, Foraging.PCA = For.pca2$x[,2], Loco.PCA = Loco.pca2$x[,2], Body.PCA = Body.pca2$x[,1])
  
  #### Can here add the foraging trait data so that each item in the lst contains all the trait data necessary to calculate Raos Q and other functional
  ### Diversity metricss
  
  two_step_PCdata <- two_step_PCdata %>% left_join(PREDICTS_Aves[,c("Jetz_Name", "Trophic_Niche")], by = "Jetz_Name") %>% 
    distinct(Jetz_Name, .keep_all = TRUE)
  
  
  Two_step_list <- c(Two_step_list, list(two_step_PCdata))
  names(Two_step_list)[i] <- paste(Two_Step_Combs[i,1],Two_Step_Combs[i,2], sep = "*")
  }


#### Now we can calculate the Raos Q for each of the 25 different datasets 


abundance_data <- PREDICTS_Aves %>%
  
  #### group by Site and Species get abundance
  group_by(SSBS,Jetz_Name) %>% mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement)) %>%
  
  ungroup() %>%
  ### group by just site to get Total site abundance
  group_by(SSBS) %>% mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance)) %>%
  
  ### filter sites that have no abundance
  ungroup() %>% filter(TotalSiteAbundance != 0, Effort_Corrected_Measurement != 0) %>%
  
  ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()


#### Need a dataframe named "abundance_data" and "trait_data"  --- Abundance data doesnt change it'll be the trait data that changes. 

Rao_Q_Func <- function(site){
  
  SiteSpecies <- abundance_data %>% filter(SSBS == site) %>% droplevels()
  
  Spp_abd <- distinct(data.frame(SiteSpecies[,c("Jetz_Name", "RelativeAbundance")]),Jetz_Name, .keep_all = TRUE)
  Spp_abd <- Spp_abd[order(Spp_abd$Jetz_Name),]
  rownames(Spp_abd) <- Spp_abd$Jetz_Name
  Spp_abd <- Spp_abd %>% select(RelativeAbundance)
  colnames(Spp_abd) <- site
  
  
  traits <- data.frame(trait_data) %>% filter(Jetz_Name %in% rownames(Spp_abd))
  traits <- traits[order(traits$Jetz_Name),]
  rownames(traits) <- traits$Jetz_Name
  traits <- traits %>% select(Trophic_Niche,Foraging.PCA,Loco.PCA,Body.PCA)
  traits$Trophic_Niche <- factor(traits$Trophic_Niche)
  
  Rao <- rao.diversity(comm = t(Spp_abd),traits =  traits)
  Rao
  
  return(Rao)
}


PREDICTS_Rao_sensitivty <- abundance_data %>% distinct(SSBS, .keep_all = TRUE)
PREDICTS_Rao_sensitivty <- PREDICTS_Rao_sensitivty[,c("Jetz_Name", "SSBS","RelativeAbundance")]
PREDICTS_Rao_sensitivty$SSBS <- as.character(PREDICTS_Rao_sensitivty$SSBS)

for(i in 1:length(Two_step_list)){
  
  trait_data <- Two_step_list[[i]]
  
  
  
  PREDICTS_Rao_sensitivty <- data.frame(PREDICTS_Rao_sensitivty) %>% mutate(col = NA)
  colnames(PREDICTS_Rao_sensitivty)[3 + i] <- names(Two_step_list)[i]
  
  for(site in PREDICTS_Rao_sensitivty$SSBS){
    RaoQ <- Rao_Q_Func(site)
    PREDICTS_Rao_sensitivty[PREDICTS_Rao_sensitivty$SSBS == site,(3 + i)] <- as.numeric(RaoQ$FunRao) 
  }
  
}


for(i in 4:NCOL(PREDICTS_Rao_sensitivty)){
hist(PREDICTS_Rao_sensitivty[,i], breaks = 30, xlim = c(0,0.7), main =  colnames(PREDICTS_Rao_sensitivty)[i])
}

