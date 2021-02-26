rm(list = ls())

require(phytools)
require(picante)
require(phyr)
require(tidyverse)
require(raster)
require(plyr)
require(ggridges)
require(motmot)
require(lme4)

## Load in the datasets

PREDICTS_site <- readRDS("Functional_Intactness_Index/Outputs/PREDICTS_Site_Rao.rds")
PREDICTS_abundance <- readRDS("Functional_Intactness_Index/Outputs/abundance_data.rds")


hist(PREDICTS_site$Bias_Rao)
hist(PREDICTS_site$Unbias_Rao)


### apply the square root transformation on the Euclidean Rao measure

PREDICTS_site$rt3Euc_Rao <- PREDICTS_site$Euc_Rao^(1/3)
hist(PREDICTS_site$rt3Euc_Rao)

PREDICTS_site$cubeRD1k <- PREDICTS_site$density_1km^(1/3)
PREDICTS_site$cubeRD50k <- PREDICTS_site$density_50km^(1/3)

## have a look at how he sites are distributed across landuse types and intenstiy

table(PREDICTS_site$LandUse, PREDICTS_site$Use_intensity)

#### collapse Plantation forest into secondary vegetation and relevel the LUI variable to have Primary Minimal use as its reference 

PREDICTS_site$SSBS <- factor(PREDICTS_site$SSBS) %>% droplevels()

PREDICTS_site$LandUse_Intensity <- ifelse(grepl(PREDICTS_site$LandUse_Intensity, pattern = "Plantation forest_Minimal use"),
                                          "Secondary Vegetation_Light use",
                                          paste(PREDICTS_site$LandUse_Intensity))

PREDICTS_site$LandUse <- ifelse(PREDICTS_site$LandUse_Intensity == "Primary_Minimal use", "Primary_Minimal use", paste(PREDICTS_site$LandUse))
PREDICTS_site$LandUse <- ifelse(PREDICTS_site$LandUse == "Plantation forest", "Secondary Vegetation", paste(PREDICTS_site$LandUse))


PREDICTS_site$LandUse_Intensity <- relevel(factor(PREDICTS_site$LandUse_Intensity), ref = "Primary_Minimal use")
PREDICTS_site$LandUse <- relevel(factor(PREDICTS_site$LandUse), ref = "Primary_Minimal use")


table(PREDICTS_site$LandUse_Intensity)
table(PREDICTS_site$LandUse)


PREDICTS_site$logHPD <- scale(PREDICTS_site$logHPD)
PREDICTS_site$Road_Density1k <- scale(sqrt(PREDICTS_site$density_1km))
PREDICTS_site$Road_Density50k <- scale(sqrt(PREDICTS_site$density_50km))

#### UN subregion and Landuse intensit were very colinear and therefore I dropped UN subregion as a factor in the modelling


source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(PREDICTS_site[,c( "LandUse_Intensity", "logHPD", "Road_Density1k","Road_Density50k" ,"UN_subregion")])


corvif(PREDICTS_site[,c( "LandUse_Intensity", "logHPD", "Road_Density1k","Road_Density50k")])





#### Load in Phylogeny - just tae the first one when running the inital models 


all_bird_tree <- read.tree("Datasets/AllBirdsHackett1.tre")
all_bird_tree <- all_bird_tree[[1]]



#############################################
#### Unifrac for all sites among studies ####
#############################################

data <- PREDICTS_abundance

among_site_unifrac <- function(data) {
  
  data <- droplevels(data)
  species <- sub(unique(data$Jetz_Name),pattern = " " ,replacement = "_")
  drop.species <- all_bird_tree$tip.label[which(!(all_bird_tree$tip.label %in% species))]

  
  study_tree <- drop.tip(all_bird_tree, drop.species)
  #study_tree <- transformPhylo(phy = study_tree, model = "lambda", lambda =  lambda)
  
  
  sites <- as.character(unique(PREDICTS_site[, "SSBS"]))
  
  comm_data <- t(species)
  colnames(comm_data) <- species
  comm_data <- data.frame(comm_data[-1,])
  
  
  
  for(i in 1:length(sites)){
    
    site_data <- data.frame(data[data$SSBS == sites[i],c("Jetz_Name", "Effort_Corrected_Measurement")])
    
    for(spp in species){
      comm_data[i,paste(spp)] <- ifelse(any(site_data[site_data$Jetz_Name == sub(spp,pattern = "_", replacement = " "),"Effort_Corrected_Measurement"] > 0),1,0)
    }
    
    
    rownames(comm_data)[i] <- sites[i]
  }
  
  
  for(i in 1:ncol(comm_data)){
    comm_data[,i] <- as.numeric(comm_data[,i])
  }
  
  comm_data <- as.matrix(comm_data)
  
  memory.limit(120000)
  
  among_site_vcv <- unifrac(comm = comm_data, tree = study_tree)

  among_site_vcv <- 1- as.matrix(among_site_vcv)
  
  return(among_site_vcv)  
}
 

among_site_vcv <- among_site_unifrac(PREDICTS_abundance)

  write_rds(file = "Functional_Intactness_Index/Outputs/Site_Among_Study_Vcv.rds", among_site_vcv)
  
  
###############################################
##### MODELLLLLLLLING #########################
###############################################

PREDICTS_site <- PREDICTS_site %>% droplevels()
  
  
among_site_vcv_1 <- readRDS("Functional_Intactness_Index/Outputs/Site_Among_Study_Vcv.rds")


first_site <- c()
for(i in 1:length(studies)){
  first <- which(grepl(studies[i],colnames(among_site_vcv_1)))[1]
  first_site <- c(first_site, first)
}

first_sites <- among_site_vcv_1[first_site,first_site]


lala <- hclust(as.dist(1 - first_sites), method = "median")

plot(lala)




Phylo_check <- c()
mod <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                              LandUse:logHPD + LandUse:density_1km +
                              (1|SS) + (1|SSB) , data = PREDICTS_site, bayes = TRUE)
mod

phy_check <- data.frame(lamdba = "No_Phylo", DIC = mod$DIC) 
Phylo_check <- rbind(Phylo_check, phy_check)

for(i in seq(0.1,1,0.1)){
  cov_mat <- among_site_unifrac(PREDICTS_abundance, lambda = i)
  mod <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                LandUse:logHPD + LandUse:density_1km +
                                (1|SS) + (1|SSB) + (1|SSBS__) ,cov_ranef = list(SSBS = cov_mat), data = PREDICTS_site, bayes = TRUE)
  
  phy_check <- data.frame(lamdba = i, DIC = mod$DIC)
  Phylo_check <- rbind(Phylo_check, phy_check)
}


?betapart::functional.betapart.core
  
Rao_Model1 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                     LandUse:logHPD + LandUse:density_1km +
                             (1|SS) + (1|SSB) + (1|SSBS__), cov_ranef = list(SSBS = among_site_vcv_1), data = PREDICTS_site)

summary(Rao_Model1)
Rao_Model1 

Rao_Model25 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                     LandUse:logHPD + LandUse:density_1km +
                                     (1|SS) + (1|SSB) + (1|SSBS__) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv_25), bayes = TRUE)

summary(Rao_Model25)



Rao_Model50 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                      LandUse:logHPD + LandUse:density_1km +
                                      (1|SS) + (1|SSB) + (1|SSBS__) , data = PREDICTS_site,
                                    cov_ranef = list(SSBS = among_site_vcv_50), bayes = TRUE)

summary(Rao_Model50)

Rao_Model75 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                      LandUse:logHPD + LandUse:density_1km +
                                      (1|SS) + (1|SSB) + (1|SSBS__) , data = PREDICTS_site,
                                    cov_ranef = list(SSBS = among_site_vcv_75), bayes = TRUE)

summary(Rao_Model75)

Rao_Model0 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse + logHPD + density_1km  +
                                      LandUse:logHPD + LandUse:density_1km +
                                      (1|SS) + (1|SSB) , data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model0)

#### adding random slopes 

Rao_Model2 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + logHPD + density_1km +
                                     LandUse_Intensity:logHPD + LandUse_Intensity:density_1km +
                                       (1|SS) + (1|SSB) + (1|SSBS__) + (logHPD|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model2)





Rao_Model3 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + logHPD + density_1km +  
                                     LandUse_Intensity:logHPD  + LandUse_Intensity:density_1km +
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (LandUse_Intensity|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model3)

Rao_Model4 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + logHPD + density_1km + 
                                     LandUse_Intensity:logHPD + LandUse_Intensity:density_1km +
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (density_1km|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model4)

#####




Random_DIC <- data.frame(Rao_Model1 = Rao_Model1$DIC, Rao_Model2 = Rao_Model2$DIC, 
                         Rao_Model3 = Rao_Model3$DIC, Rao_Model4 = Rao_Model4$DIC)

#### second set Rao Model 4 looks to be the best 


#### Remove each interaction to see if thsis improves the model at all.

## LandUSE Intensity:RoadD1

Rao_Model5 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + logHPD + density_1km  + 
                                     LandUse_Intensity:logHPD +
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (density_1km|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE )

summary(Rao_Model5)  ### didn't significantly change the DIC of the model ~ 0.1 difference

## log HPd so lets see what happens if I remove the loghpd interaction

Rao_Model6 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + logHPD + density_1km +  
                                      LandUse_Intensity:density_1km +
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (density_1km|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model6)


####

ModelDiffs <- data.frame(RoadDint =  Rao_Model5$DIC  - Rao_Model4$DIC, 
                         LogHPD =  Rao_Model6$DIC - Rao_Model4$DIC)



## Model significantly improved if i removed the interaction between landUse)intensity and lopHPD
## so lets try to remove log HPD from the model

Rao_Model7 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + density_1km +  
                                     LandUse_Intensity:density_1km +
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (density_1km|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model7)


(Rao_Model7$DIC) - Rao_Model6$DIC


### DIC goes up model is not better

## remove the other interaction

Rao_Model8 <- phyr::communityPGLMM(rt3Euc_Rao ~ LandUse_Intensity + density_1km + logHPD  
                                     (1|SS) + (1|SSB) + (1|SSBS__) + (density_1km|SS) , data = PREDICTS_site,
                                   cov_ranef = list(SSBS = among_site_vcv), bayes = TRUE)

summary(Rao_Model8)


(Rao_Model8$DIC) - Rao_Model6$  ### again does not improve model therefore Rao 6 model is retained


#### fcuntion to r=visualise the residuals from the PGLMM and check whether they are normally distributed 

residuals.communityPGLMM <- function(
  object, 
  type = if(object$family %in% c("binomial","poisson")) "deviance" else "response",
  scaled = FALSE, ...){
  if(object$family == "gaussian"){
    y <- object$Y
    mu <- pglmm_predicted_values(object)$Y_hat
    res <- switch(type,
                  deviance = stop("no deviance residuals for gaussian model", call. = FALSE),
                  response = y - mu
    )
    if(scaled) res/sqrt(object$s2resid)
  }
  
  if(object$family %in% c("binomial","poisson")){
    y <- as.numeric(object$Y)
    mu <- unname(object$mu[, 1])
    if(object$family == "binomial") dres <- sqrt(binomial()$dev.resids(y, mu, 1))
    if(object$family == "poisson") dres <- sqrt(poisson()$dev.resids(y, mu, 1))
    res <- switch(type,
                  deviance = {
                    dres
                    ifelse(y > mu, dres, - dres)
                  },
                  response = y - mu
    )
  }
  if(object$family %nin% c("gaussian", "binomial", "poisson"))
    stop("no residual methods for family other than gaussian, binomial and poisson, yet", call. = FALSE)
  
  unname(res)
}


Resid_check <- residuals.communityPGLMM(Rao_model)

plot(Resid_check)
  
qqnorm(Resid_check)
qqline(Resid_check)

write_rds(Rao_Model6, file = "Functional_Intactness_Index/Rao_model.rds")


#######################################
######## Phylo autocorrelation ########
#######################################

##### We are looking at the phylogenetic signal of the residuals of the model - first we need to decide whther there is any phylogenetic signal wwithin the 
### model that we need to account for. We are going to do this by 

################################
#### Mthod 2 ##################
###############################



#1. fit model with no phylogeny in error term 
PREDICTS_test <- PREDICTS_site %>% filter(SS != "HP1_2010__Lasky 2")



test_mod <- lmer(Unbias_Rao ~ LandUse + logHPD + Road_Density1k  + Road_Density50k +
                   LandUse:logHPD + LandUse:Road_Density1k + LandUse:Road_Density50k +
                   (1|SS) + (1|SSB) + (1 + LandUse|SS), data = PREDICTS_site)


summary(test_mod)
car::Anova(test_mod, type = "III")

#### extract the random effects 

Randoms <- lme4::ranef(test_mod)

study_slope <- Randoms$SS


#### Do more closely related studies have similar random effects 

data <- PREDICTS_abundance

study_phylosim <- function(data) {
  
  data <- droplevels(data)
  species <- sub(unique(data$Jetz_Name),pattern = " " ,replacement = "_")
  drop.species <- all_bird_tree$tip.label[which(!(all_bird_tree$tip.label %in% species))]
  
  
  study_tree <- drop.tip(all_bird_tree, drop.species)

  
  
  studies <- as.character(unique(PREDICTS_site[, "SS"]))
  
  
  
  comm_data <- t(species)
  colnames(comm_data) <- species
  comm_data <- data.frame(comm_data[-1,])
  
  
  
  for(i in 1:length(studies)){
    
    study_data <- data.frame(data[data$SS == studies[i],c("Jetz_Name", "Effort_Corrected_Measurement")])
    
    for(spp in species){
      comm_data[i,paste(spp)] <- ifelse(any(study_data[study_data$Jetz_Name == sub(spp,pattern = "_", replacement = " "),"Effort_Corrected_Measurement"] > 0),1,0)
    }
    
    
    rownames(comm_data)[i] <- studies[i]
  }
  
  
  for(i in 1:ncol(comm_data)){
    comm_data[,i] <- as.numeric(comm_data[,i])
  }
  
  comm_data <- as.matrix(comm_data)
  
  

  study_sim <- betapart::phylo.beta.pair(comm_data, study_tree, index.family = "jaccard")
  
  #study_sim <- pcd(comm_data, study_tree)
  #study_sim <- as.matrix(study_sim$PCD)
  
  
  study_sim <- 1 - as.matrix(study_sim$phylo.beta.jac)
  
  
  
  
  #among_site_vcv <- 1 - as.matrix(unifrac(comm = comm_data, tree = study_tree))
  
  
  return(study_sim)  
}




Study_sim <- study_phylosim(PREDICTS_abundance)

studies <- as.character(unique(PREDICTS_site[, "SS"]))


significance <- c()

for(position in 1:c(length(studies)-1)){

position_diff <- c()
relatives <- c()

for(study in studies){
  
  position_rel <- colnames(Study_sim)[order(Study_sim[study,],decreasing = TRUE)[position + 1]]
  rel_diff <- abs(study_slope[study,-7] - study_slope[position_rel,-7])
  rownames(rel_diff) <- paste(position_rel)

  
  position_diff <- rbind(position_diff, rel_diff)
  relatives <- rbind(relatives,data.frame(stud1 = study, stud2 = position_rel) )
  
}

stud_com <- data.frame(gtools::combinations(v = studies, r = 2, n = length(studies) - 1))

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



p_val <- wilcox.test(position_diff$`LandUseSecondary Vegetation`, global_differences$`LandUseSecondary Vegetation`, alternative = "less")
p_val <- data.frame(Position = position, significance =  p_val$p.value)

significance <- rbind(significance, p_val)

}


plot(significance$significance ~ significance$Position)

poolr::fisher(significance$significance)


ggResidpanel::resid_panel(test_mod)
