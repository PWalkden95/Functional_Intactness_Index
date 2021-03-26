rm(list = ls())

require(phytools)
require(picante)
require(phyr)
require(tidyverse)
require(poolr)
require(motmot)

PREDICTS_site <- readRDS("Outputs/PREDICTS_Site_Rao.rds")
PREDICTS_abundance <- readRDS("Outputs/Rao_abundance_data.rds")

#### Load in Phylogeny - just tae the first one when running the inital models 


all_bird_tree <- read.tree("../Datasets/AllBirdsHackett1.tre")
all_bird_tree <- all_bird_tree[[1]]



#############################################
#### Unifrac for all sites among studies ####
#############################################


source("Functions/create_vcv.R")

site_vcv <- create_vcv(data = PREDICTS_abundance, level = "SSBS", tree = all_bird_tree)

write_rds(file = "Outputs/site_vcv.rds", site_vcv)


########################################################################################################################
##### Identify whether there is phylogenetic signal of in the responses of sites of LU change  #########################
########################################################################################################################

### first we are going to create a cluster dendrogram based on the phylogenetic distances between sites as calculated by unifrac
### Using an non-ultmetric tree we will be able to see just how heirarchial the resulting tree is, and whether this needs to be resolved
### in the modelling.   


PREDICTS_site <- PREDICTS_site %>% droplevels()

among_site_vcv <- readRDS("Outputs/site_vcv.rds")

studies <- as.character(unique(PREDICTS_site$SS))

## get the first site from each study as it is very difficult to observe the full tree with all sites

first_site <- c()
for(i in 1:length(studies)){
  first <- which(grepl(studies[i],colnames(among_site_vcv)))[1]
  first_site <- c(first_site, first)
}

first_sites <- among_site_vcv[first_site,first_site]

### performing a clustering algorithm based on the distances between sites 

site_dendro <- hclust(as.dist(1 - first_sites), method = "median")

plot(site_dendro)

### Can actually see that there is not a large amount of heirarchial structure within the sites with only studies conducted simultaneously
### showing a greater amount of similarity compared with the rest of the studies 

#### This initially does not suggest that there is a large amount of phylogenetic similarity between sites resulting in hierarcial and
#### correlated responses to land-use change - but another check I am going to conduct is to see whether the residuals of the random
#### slopes of landuse within studies are more similar in studies that are more phylogenetically similar compared to the overall distribution
#### of slope differences. 



#### First calculate study level phylogenetic similarity 


Study_sim <- create_vcv(PREDICTS_abundance, level = "SS", tree = all_bird_tree)

### Next run a model without accounting for phylogenetic signal of sites and extract the random slopes of response to land use across studies

#1. fit model with no phylogeny in error term 


test_mod <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k  + RD50k + CNTRLlogHPD +
                   (1|SS) + (1|SSB) + (1 + LUI|SS), data = PREDICTS_site)


summary(test_mod)
car::Anova(test_mod, type = "III")

#### extract the random effects 

Randoms <- lme4::ranef(test_mod)

study_slope <- Randoms$SS





significance <- c()

for(position in 1:c(length(studies)-1)){
  
  position_diff <- c()
  relatives <- c()
  
  for(study in studies){
    
    position_rel <- colnames(Study_sim)[order(Study_sim[study,],decreasing = TRUE)[position + 1]]
    rel_diff <- abs(study_slope[study,] - study_slope[position_rel,])
    rownames(rel_diff) <- paste(position_rel)
    
    
    position_diff <- rbind(position_diff, rel_diff)
    relatives <- rbind(relatives,data.frame(stud1 = study, stud2 = position_rel) )
    
  }
  
  stud_com <- data.frame(gtools::combinations(v = studies, r = 2, n = c(length(studies)-1)))
  
  drops <- c()
  for(k in 1:nrow(relatives)){
    drop_row <- which(stud_com$X1 == relatives[k,1] & stud_com$X2 == relatives[k,2] | stud_com$X1 == relatives[k,2] & stud_com$X2 == relatives[k,1] )
    drops <- c(drops,drop_row)
  }
  
  stud_com <- stud_com[-unique(drops),]
  
  global_differences <- c()
  for(i in 1:NROW(stud_com)){
    abs_diff <- abs(study_slope[stud_com[i,1],] - study_slope[stud_com[i,2],])
    rownames(abs_diff) <- paste(stud_com[i,2]) 
    
    global_differences <- rbind(global_differences, abs_diff)
  }
  
  
  
  p_val <- wilcox.test(position_diff[,2], global_differences[,2], alternative = "less")
  p_val <- data.frame(Position = position, significance =  p_val$p.value)
  
  significance <- rbind(significance, p_val)
  
}

#### plot doesn't clearly indicate that the responses of sites within studies are significantly more similar in more phylogenetically
#### similar sites 

plot(significance$significance ~ significance$Position)

### fishers method - combined p values testing whether the hypothesis test (are the difference in slopes significantly more similar
### to each other in studies that are more phylogenetically related?)


### No they are not 

poolr::fisher(significance$significance)
