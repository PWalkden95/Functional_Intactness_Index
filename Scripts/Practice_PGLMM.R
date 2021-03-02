rm(list = ls())

require(phytools)
require(picante)
require(phyr)
require(tidyverse)
require(raster)
require(ggridges)
require(motmot)
require(lme4)
require(poolr)
require(car)
require(robustlmm)

## Load in the datasets

PREDICTS_site <- readRDS("Outputs/PREDICTS_Site_Rao.rds")
PREDICTS_abundance <- readRDS("Outputs/abundance_data.rds")




hist(PREDICTS_site$Bias_Rao, breaks = 20)
hist(PREDICTS_site$Unbias_Rao, breaks = 20)

plot(PREDICTS_site$Bias_Rao ~ PREDICTS_site$Unbias_Rao)


###### First we will want outliers 

##Outleir check 

qqPlot(PREDICTS_site$Unbias_Rao)
### apply the square root transformation on the Euclidean Rao measure

PREDICTS_site$sqrtRD1k <- PREDICTS_site$density_1km^(1/2)
PREDICTS_site$sqrtRD50k <- PREDICTS_site$density_50km^(1/2)

## have a look at how he sites are distributed across landuse types and intenstiy

table(PREDICTS_site$LandUse, PREDICTS_site$Use_intensity)

#### collapse Plantation forest into secondary vegetation and relevel the LUI variable to have Primary Minimal use as its reference 

PREDICTS_site$SSBS <- factor(PREDICTS_site$SSBS) %>% droplevels()

PREDICTS_site <- PREDICTS_site %>% dplyr::mutate(LandUse_Intensity = ifelse(grepl(LandUse_Intensity, pattern = "Plantation forest_Minimal use"),
                                          "Secondary Vegetation_Light use",
                                          paste(LandUse_Intensity)),
                                          LandUse_Intensity = ifelse(grepl(LandUse_Intensity, pattern = "Plantation forest_Light use"),
                                                                     "Secondary Vegetation_Intense use",
                                                                      paste(LandUse_Intensity)),
                                          LandUse_Intensity = relevel(factor(LandUse_Intensity), 
                                                                      ref = "Primary_Minimal use"))

table(PREDICTS_site$LandUse_Intensity)

########## Going to rename for some ease of outputs
PREDICTS_site <- PREDICTS_site %>% dplyr::rename(LUI = LandUse_Intensity)
levels(PREDICTS_site$LUI) <- c("PriMin", "CrpLig", "CrpMin", "PasIn", "PasLig", 
                               "PasMin","PriIn","PriLig","SecIn","SecLig","SecMin","UrbLig","UrbMin")

PREDICTS_site <- PREDICTS_site %>% dplyr::mutate(LandUse = ifelse(LUI == "PriMin", "PriMin", paste(LandUse)),
                                                 LandUse = ifelse(LandUse == "Plantation forest", "Secondary Vegetation", paste(LandUse)),
                                                 LandUse = relevel(factor(LandUse), ref = "PriMin"))



table(PREDICTS_site$LUI)
table(PREDICTS_site$LandUse)


PREDICTS_site$logHPD <- scale(PREDICTS_site$logHPD)
PREDICTS_site$RD1k <- scale(PREDICTS_site$sqrtRD1k)
PREDICTS_site$RD50k <- scale(PREDICTS_site$sqrtRD50k)
PREDICTS_site$CNTRLlogHPD <- scale(PREDICTS_site$CNTRLlogHPD)

#### UN subregion and Landuse intensit were very colinear and therefore I dropped UN subregion as a factor in the modelling


source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(PREDICTS_site[,c( "LUI", "logHPD", "RD1k","RD50k" ,"UN_subregion")])


corvif(PREDICTS_site[,c( "LUI", "logHPD",  "RD1k","RD50k")])


#### Road densitity at 50 km radius is colinear witht the other variables and is dropped 


#### Load in Phylogeny - just tae the first one when running the inital models 


all_bird_tree <- read.tree("Datasets/AllBirdsHackett1.tre")
all_bird_tree <- all_bird_tree[[1]]



#############################################
#### Unifrac for all sites among studies ####
#############################################


create_vcv <- function(data, level) {
  
  data <- droplevels(data)
  species <- sub(unique(data$Jetz_Name),pattern = " " ,replacement = "_")
  drop.species <- all_bird_tree$tip.label[which(!(all_bird_tree$tip.label %in% species))]

  
  full_tree <- drop.tip(all_bird_tree, drop.species)
  
  
  ID <- as.character(unique(PREDICTS_site[, level]))
  
  comm_data <- t(species)
  colnames(comm_data) <- species
  comm_data <- data.frame(comm_data[-1,])
  
  
  
  for(i in 1:length(ID)){
    
    ID_data <- data.frame(data[data[,level] == ID[i],c("Jetz_Name", "Effort_Corrected_Measurement")])
    
    for(spp in species){
      comm_data[i,paste(spp)] <- ifelse(any(ID_data[ID_data$Jetz_Name == sub(spp,pattern = "_", replacement = " "),"Effort_Corrected_Measurement"] > 0),1,0)
    }
    
    
    rownames(comm_data)[i] <- ID[i]
  }
  
  
  for(i in 1:ncol(comm_data)){
    comm_data[,i] <- as.numeric(comm_data[,i])
  }
  
  comm_data <- as.matrix(comm_data)
  
  suppressWarnings(memory.limit(120000))
  
  vcv <- 1 - as.matrix(unifrac(comm = comm_data, tree = full_tree))
 
  

  
  return(vcv)  
}
 

among_site_vcv <- create_vcv(PREDICTS_abundance, level = "SSBS")

  write_rds(file = "Functional_Intactness_Index/Site_Among_Study_Vcv.rds", among_site_vcv)
  
  
########################################################################################################################
##### Identify whether there is phylogenetic signal of in the responses of sites of LU change  #########################
########################################################################################################################

### first we are going to create a cluster dendrogram based on the phylogenetic distances between sites as calculated by unifrac
### Using an non-ultmetric tree we will be able to see just how heirarchial the resulting tree is, and whether this needs to be resolved
### in the modelling.   
  
  
PREDICTS_site <- PREDICTS_site %>% droplevels()
  
among_site_vcv <- readRDS("Functional_Intactness_Index/Outputs/Site_Among_Study_Vcv.rds")

studies <- as.character(unique(PREDICTS_site$SS))

## get the first site from each study as it is very difficult to observe the full tree with all sites

first_site <- c()
for(i in 1:length(studies)){
  first <- which(grepl(studies[i],colnames(among_site_vcv)))[1]
  first_site <- c(first_site, first)
}

first_sites <- among_site_vcv_1[first_site,first_site]

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


Study_sim <- create_vcv(PREDICTS_abundance, level = "SS")

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


####################################
# With this we can proceed with modelling using GLMMs as opposed to PGLMMs 
###################################

Rao_Model_1 <- lmer(Unbias_Rao ~ LUI+ logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB), data = PREDICTS_site)

summary(Rao_Model_1)

Rao_Model_1b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)


summary(Rao_Model_1b)

#### lets see what random effect structure gives the best AIC 

### adding random slopes 1) LUI within study 

Rao_Model_2 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + LUI|SS), data = PREDICTS_site)

summary(Rao_Model_2)

Rao_Model_2b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_2b)

#2) logHPD within study 

Rao_Model_3 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + logHPD|SS), data = PREDICTS_site)

summary(Rao_Model_3)

Rao_Model_3b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (logHPD|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_3b)

#3) Road Density_1km within study  

Rao_Model_4 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + RD1k|SS), data = PREDICTS_site)

summary(Rao_Model_4)

Rao_Model_4b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (RD1k|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_4b)

# 4) RD50km within study 

Rao_Model_5 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + RD50k|SS), data = PREDICTS_site)

summary(Rao_Model_5)


Rao_Model_5b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (RD50k|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_5b)


Random_AIC <- data.frame(mod1 = AIC(Rao_Model_1), mod2 = AIC(Rao_Model_2),mod3 = AIC(Rao_Model_3),
                         mod4 = AIC(Rao_Model_4), mod5 = AIC(Rao_Model_5))

Random_DIC <- data.frame(mod1 = Rao_Model_1b$DIC, mod2 = Rao_Model_2b$DIC,mod3 = Rao_Model_3b$DIC,
                         mod4 = Rao_Model_4b$DIC, mod5 = Rao_Model_5b$DIC)


#### the DIC and AIC are both lowest with the rabndom effect structure of random slope of logHPD within study (MODEL 3)

### Model 3 gives the best result which has logHPD as a random slope within study 
## so now to test the best fixed effect structure

car::Anova(Rao_Model_1, type = "III") ## desnity of roads 50k is the least significant interaction so will test its exclusion


### GLMM fails to converge removing the interaction between Road_density50k and landuse so I will continue with Bayes GLMM (INLA)

#### Remove each interaction to see if thsis improves the model at all. 

#LUI:RD50k

Rao_Model_6b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



#LUI:RD1k

Rao_Model_7b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD50k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)


# LUI:logHPD

Rao_Model_8b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:RD50k + LUI:RD1k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)




Model_differences <- data.frame(RD50k = Rao_Model_6b$DIC - Rao_Model_1b$DIC,
                                RD1k = Rao_Model_7b$DIC - Rao_Model_1b$DIC,
                                logHPD = Rao_Model_8b$DIC - Rao_Model_1b$DIC)   


### Model improves the most when I removed the interaction between LUI:RD1k

### Next round of removals -- 

## fixed effect RD1k

Rao_Model_9b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD50k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)




# LUI:RD50k

Rao_Model_10b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                                        LUI:logHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



# LUI:logHPD

Rao_Model_11b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                                        LUI:RD50k +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)





Model_differences <- data.frame(RD1k = Rao_Model_9b$DIC - Rao_Model_7b$DIC,
                                RD50k = Rao_Model_10b$DIC - Rao_Model_7b$DIC,
                                logHPD = Rao_Model_11b$DIC - Rao_Model_7b$DIC)   

Model_differences

#### Model DIC is reduced in all but most when removing the interaction between LUI:RD50k so model 10b is best 

## RD50k fixed

Rao_Model_12b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                                        LUI:logHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



### RD1k fixed

Rao_Model_13b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + CNTRLlogHPD +
                                        LUI:logHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



#### LUI:logHPD


Rao_Model_14b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



Model_differences <- data.frame(RD50k = Rao_Model_12b$DIC - Rao_Model_10b$DIC,
                                RD1k = Rao_Model_13b$DIC - Rao_Model_10b$DIC,
                                logHPD = Rao_Model_14b$DIC - Rao_Model_10b$DIC)   

Model_differences


#### Model DIC is reduced wehn removing the interaction between LUI:logHPD and is increased in the other cases so model 14 is carried forward.

#### RD50k


Rao_Model_15b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



#### RD1k

Rao_Model_16b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD + RD50k + CNTRLlogHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



### logHPD


Rao_Model_17b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + RD1k + RD50k + CNTRLlogHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)




Model_differences <- data.frame(RD50k = Rao_Model_15b$DIC - Rao_Model_14b$DIC,
                                RD1k = Rao_Model_16b$DIC - Rao_Model_14b$DIC,
                                logHPD = Rao_Model_17b$DIC - Rao_Model_14b$DIC)  

####


Model_differences

### the model wasn't improved with the removal of any of the fixed effects therefore the best model is Model_14b

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

   resid_check <- residuals.communityPGLMM(Rao_Model_14b)
plot(resid_check)   
hist(resid_check)


qqPlot(resid_check)
### looking at the residuals it seems that there are some influential points that are caused by flocks of birds dominating the site in terms 
## or relative abundance resulting in very low estimates of functional diversity. Because these values are true values we want to still include
### them in the model but we would like to down-weight their contribution to the model therefore I am going to perform a Robustlmm with 
## the final selected model 

Robust_mod <- rlmer(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                      (1|SS) + (1|SSB), data = PREDICTS_site)

summary(Robust_mod)
summary(Rao_Model_14b)

plot(Robust_mod)

###### so lets have a look at whats going on.

utils::sessionInfo()

### Functional diversity is significantly lower in Cropland vs Primary minimal habitat
