rm(list = ls())

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(SYNCSA) ## For calculating functional diversity metric
require(TPD)

Jetz_Traits <- read.csv("../Datasets/GBD/GBD_BiometricsRaw_combined_2020_Dec_07.csv")
colnames(Jetz_Traits)[6] <- "Jetz_Name"

PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")
morpho_traits <- readRDS("Outputs/FullMorphTraits.rds")




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

abundance_data <- PREDICTS_Aves_Am %>% dplyr::filter(Diversity_metric == "abundance" |Diversity_metric == "effort-corrected-abundance" ) %>%
  
  dplyr::filter(Effort_Corrected_Measurement != 0) %>%
  
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

######################################################################
#---------------------------------------------------------------------
#---------------------------------------------------------------------

Jetz_Species <- PREDICTS_Aves_Am %>% filter(Diversity_metric == "abundance" | Diversity_metric == "effort_corrected_abundance" ) %>% distinct(Jetz_Name)


Jetz_Traits <- Jetz_Traits %>% dplyr::filter(Jetz_Name %in% Jetz_Species$Jetz_Name) %>% 
  dplyr::select(Jetz_Name,Bill_TotalCulmen,Bill_Depth,Bill_Nares,Bill_Width,
                Wing_Chord,Secondary1,Tail_Length,Tarsus_Length) %>% na.omit() %>% group_by(Jetz_Name) %>%
  dplyr::mutate(num = n()) %>% ungroup()
################# Fewer than 4 triats

Jetz_Traits_Am <- Jetz_Traits %>% filter(num >= 4) %>% pull(Jetz_Name) %>% unique()
Jetz_Traits_under <- Jetz_Traits %>% filter(num < 4 & num != 1) %>% pull(Jetz_Name) %>% unique()
single_Jetz_traits <- Jetz_Traits %>% filter(num == 1) %>% pull(Jetz_Name) %>% unique()


PRED_spp <- unique(abundance_data$Jetz_Name)
miss_spp <- PRED_spp[which(!(PRED_spp %in% c(Jetz_Traits_Am,Jetz_Traits_under,single_Jetz_traits)))]


data <- abundance_data
traits <- morpho_traits

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


trait_ranges <- list(c(min(traits[,2]) -(0.1 * min(traits[,2])),max(traits[,2]) + (0.1 * max(traits[,2]))),
                     c(min(traits[,3]) -(0.1 * min(traits[,3])),max(traits[,3]) + (0.1 * max(traits[,3]))),
                     c(min(traits[,4]) -(0.1 * min(traits[,4])),max(traits[,4]) + (0.1 * max(traits[,4]))))

species <- data.frame(Jetz_Name = rownames(Species_abundance))

mean_TPD <- c()

if(any(species$Jetz_Name %in% c(single_Jetz_traits, Jetz_Traits_under))){
  under_spp <- species %>% filter(Jetz_Name %in% c(single_Jetz_traits, Jetz_Traits_under,miss_spp))
  under_spp_traits <- under_spp %>% dplyr::left_join(traits, by = "Jetz_Name")
  
  spp <- unique(under_spp_traits$Jetz_Name)
  k <- spp[2]
  means <- c()
  sd <- c()
  for(k in spp){
    spp_mean <- data.frame(Jetz_Name = k)
    spp_sd <- data.frame(Jetz_Name = k)
    spp_data <- under_spp_traits %>% filter(Jetz_Name == k)
    if(nrow(spp_data) == 1){
      spp_mean <- spp_data
      spp_sd <- cbind(spp_sd, Foraging.PCA =  c(0.5*sd(traits[,2])),
                      Loco.PCA =  c(0.5*sd(traits[,3])),
                      Body.PCA =  c(0.5*sd(traits[,3])))
    } else{
      
      for(j in 2:ncol(under_spp_traits)){
        mean <-  mean(under_spp_traits[,j])
        stan <-  sd(under_spp_traits[,j])
        
        spp_mean <- cbind(spp_mean,mean)
        colnames(spp_mean)[j] <- colnames(under_spp_traits)[j]
        spp_sd <- cbind(spp_sd, stan)
        colnames(spp_sd)[j] <- colnames(under_spp_traits)[j]
      }}
    means <- rbind(means,spp_mean)
    sd <- rbind(sd,spp_sd)
  }
  
  
  mean_TPD <- TPDsMean(species = spp, means = means[,c(2:4)], sds = sd[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 61, alpha = 0.95 )
}


over_traits <- spp_traits %>% filter(Jetz_Name %in% Jetz_Traits_Am)

sites <- as.character(unique(abundance_data$SSBS))

Trait_density <- TPDs(over_traits[,1], over_traits[,c(2:4)], trait_ranges = trait_ranges, n_divisions = 61, alpha = 0.95)


if(!is.null(mean_TPD)){
  Trait_density$TPDs <- c(mean_TPD$TPDs, Trait_density$TPDs)
  Trait_density$data$species <- c(Trait_density$data$species, under_spp_traits$Jetz_Name)
  Trait_density$data$traits <- rbind(Trait_density$data$traits, under_spp_traits[,c(2:4)])
}

memory.limit(120000)

Overlap <- TPDc(Trait_density, t(Species_abundance))

dissim <- TPD::dissim(Trait_density)

Rao <- TPD::Rao(diss = dissim, TPDc = Overlap)


hist(Rao$alpha_rao)


write_rds(dissim, file = "Outputs/tpd_species_dissimilarity.rds")
write_rds(Rao, file = "Outputs/tpd_rao.rds")
