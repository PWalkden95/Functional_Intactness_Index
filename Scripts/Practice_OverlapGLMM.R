rm(list = ls())

require(tidyverse)
require(car)
require(lme4)
require(ggResidpanel)



#### Load in Similarity data

Overlap_data <- readRDS("Functional_Intactness_Index/Functional_Overlap_data.rds")
Similarity_data <- readRDS("Functional_Intactness_Index/similarity_data.rds")

## extract the lowest number of species between the two sites 


Overlap_data <- Overlap_data %>% left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site1" = "SSBS")) %>%
  dplyr::rename(site1_spp = Site_spp) %>%
left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site2" = "SSBS")) %>%
  dplyr::rename(site2_spp = Site_spp)

Overlap_data$min_site_spp <- ifelse(Overlap_data$site1_spp > Overlap_data$site2_spp, Overlap_data$site2_spp, Overlap_data$site1_spp) 

  
  
### Want to work out human population density at 


Overlap_data <- Overlap_data %>% dplyr::mutate(logitOver = car::logit(hyper_overlap, adjust = 0.001, percents = FALSE)) %>%
  dplyr::mutate(logdist = log(distance + 1)) %>%
  dplyr::mutate(Contrast = factor(Contrast),
         Contrast = relevel(Contrast, ref = "Primary_Minimal use-Primary_Minimal use"))





table(Overlap_data$Contrast)

Overlap_data$Contrast <- ifelse(grepl(Overlap_data$Contrast, pattern = "Plantation forest"), 
                                   "Primary_Minimal use-Secondary Vegetation_Light use",
                                   paste(Overlap_data$Contrast))

table(Overlap_data$Contrast)

Overlap_data$Contrast <- relevel(factor(Overlap_data$Contrast), ref = "Primary_Minimal use-Primary_Minimal use")

Overlap_data$env_distance <- scale(Overlap_data$env_distance)
Overlap_data$s2logHPD <- scale(Overlap_data$S2logHPD)
Overlap_data$logHPDdiff <- scale(Overlap_data$logHPDdiff)
Overlap_data$logdist <- scale(Overlap_data$logdist)
Overlap_data$S2RD1K <- scale(Overlap_data$S2RD1K)
Overlap_data$S2RD50K <- scale(Overlap_data$S2RD50K)
Overlap_data$RD1Kdiff <- scale(Overlap_data$RD1Kdiff)
Overlap_data$RD50Kdiff <- scale(Overlap_data$RD50Kdiff)

source("https://highstat.com/Books/Book2/HighstatLibV10.R")

corvif(Overlap_data[,c("Contrast", "logdist", "env_distance","S2logHPD", "logHPDdiff", "S2RD1K", "S2RD50K", "RD1Kdiff", "RD50Kdiff")])

### some colineararity problmems lets drop site2HPD

corvif(Overlap_data[,c("Contrast", "logdist", "env_distance", "logHPDdiff", "S2RD1K", "S2RD50K", "RD1Kdiff", "RD50Kdiff")])


### This resolved the issues

write_rds(file = "Functional_Intactness_Index/Overlap_modelling_data.rds", Overlap_data)


################################
##### Modelling ################
################################

Overlap_data <- readRDS("Functional_Intactness_Index/Overlap_modelling_data.rds")

Model_1 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                 (1|SS),
              data = Overlap_data, weights = min_site_spp )

summary(Model_1)

AIC(Model_1)

#### random slopes

## logdist

Model_2 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + logdist|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_2)

AIC(Model_2)

## env_distance

Model_3 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + env_distance|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_3)

AIC(Model_3)


## logHPDdiff

Model_4 <- lmer(logitOver ~ Contrast + logdist + env_distance  + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + logHPDdiff|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_4)

AIC(Model_4)

## S2RD1K

Model_5 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                   (1|SS) + (1 + S2RD1K|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_5)

AIC(Model_5)

## S2RD50K

Model_6 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + S2RD50K|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_6)

AIC(Model_6)

## RoadDdiff1k

Model_7 <- lmer(logitOver ~ Contrast + logdist + env_distance  + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance + Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + RD1Kdiff|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_7)

AIC(Model_7)

## RoadDdiff50k

Model_8 <- lmer(logitOver ~ Contrast + logdist + env_distance  + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff + Contrast:RD50Kdiff +
                  (1|SS) + (1 + RD50Kdiff|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_8)

AIC(Model_8)



MOD_AIC <- data.frame(Mod1 = AIC(Model_1),Mod2 = AIC(Model_2),Mod3 = AIC(Model_3),
                      Mod4 = AIC(Model_4),Mod5 = AIC(Model_5),Mod6 = AIC(Model_6),
                      Mod7 = AIC(Model_7),Mod8 = AIC(Model_8))

#### Model simplification - no indication for the need of a random slope for any of the continuou variables within site



###### Permutated dataset x 1000

Permuted_data <- rep(list(NA),1000)

for(i in 1:1000){
  
  sample_data <-c()
  
  for(study in unique(Overlap_data$SS)){
    data <- Overlap_data %>% filter(SS == study)
    
    data$logitOver <- data[sample(NROW(data)),"logitOver"]
    
    sample_data <- rbind(sample_data,data)
     
  }
  
  Permuted_data[[i]] <- sample_data
  
}

##### Liklihood ratio function

Permuted_model_simplification <- function(data, model1, model2){

  LRT_dist <- c()
  
  for(i in 1:length(data)){
  
  mod1 <- lmer(model1@call, data = data[[i]], weights = min_site_spp, REML = FALSE)
  mod2 <- lmer(model2@call, data = data[[i]], weights = min_site_spp, REML = FALSE)
  
  LRT <- anova(mod1,mod2)
  LRT <- LRT[which(!is.na(LRT$Chisq)),"Chisq"]
  
  LRT_dist <- rbind(LRT_dist,LRT)
  
  }

  return(LRT_dist)  
  
}


Anova(Model_1, type = "II")


######## dropping interaction between Contrast:RoadDdiff50k

Model_9 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff + RD50Kdiff + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff  +
                  (1|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_9)

AIC(Model_9)


Like_ratio <- anova(Model_1, Model_9)


LRT_dist1 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, model2 = Model_9)

ninety_five <- quantile(LRT_dist1,0.95)

Like_ratio[2,"Chisq"] > ninety_five


##### LRT of models using the observed data is not significantly different from the distribution of LRT using the permuted data therefore the
##### interaction can be dropped from the model.

Anova(Model_9,type = "II")



### dropping roadDdiff50k

Model_10 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                  S2RD1K + S2RD50K + RD1Kdiff  + 
                  Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                  Contrast:S2RD1K + Contrast:S2RD50K + Contrast:RD1Kdiff  +
                  (1|SS),
                data = Overlap_data, weights = min_site_spp )

summary(Model_10)

AIC(Model_10)

Like_ratio <- anova(Model_9, Model_10)


LRT_dist2 <- Permuted_model_simplification(Permuted_data,model1 = Model_9, model2 = Model_10)

ninety_five <- quantile(LRT_dist2, 0.95)


Like_ratio[2,"Chisq"] > ninety_five

###### Likelihood ratio is not significantly greater than the distribution therefore roadDdiff50k is also dropped from the model

Anova(Model_10, type = "II")

#### try removing s2RD1k interaction


Model_11 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                   S2RD1K + S2RD50K + RD1Kdiff  + 
                   Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                   Contrast:S2RD50K + Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_11)

AIC(Model_11)


Like_ratio <- anova(Model_10, Model_11)


LRT_dist3 <- Permuted_model_simplification(Permuted_data,model1 = Model_10, model2 = Model_11)

ninety_five <- quantile(LRT_dist3, 0.95)


Like_ratio[2,"Chisq"] > ninety_five


##### the liklihood ratio is not sigificantly greater than the distribution therefore the simplified model can proceed

###

Anova(Model_11, type = "II")

#### Lets try to remove the interaction between Contrast and S2RD50k

Model_12 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                   S2RD1K + S2RD50K + RD1Kdiff  + 
                   Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                   Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_12)

AIC(Model_12)

Like_ratio <- anova(Model_11, Model_12)


LRT_dist4 <- Permuted_model_simplification(Permuted_data,model1 = Model_11, model2 = Model_12)

ninety_five <- quantile(LRT_dist4, 0.95)


Like_ratio[2,"Chisq"] > ninety_five

##### the liklihood ratio is not sigificantly greater than the distribution therefore the simplified model can proceed

Anova(Model_12, type = "II")

### lets try to remove S2RD50K

Model_13 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                   S2RD1K + RD1Kdiff  + 
                   Contrast:logdist + Contrast:env_distance +  Contrast:logHPDdiff +
                   Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_13)

AIC(Model_13)

Like_ratio <- anova(Model_12, Model_13)


LRT_dist5 <- Permuted_model_simplification(Permuted_data,model1 = Model_12, model2 = Model_13)

ninety_five <- quantile(LRT_dist5, 0.95)


Like_ratio[2,"Chisq"] > ninety_five


### The liklihood ratio s not sigificantly greater than the distribution therefore the simplified model can proceeed


Anova(Model_13, type = "III")


### so lets remove contrast:logdist

Model_14 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                   S2RD1K + RD1Kdiff  + 
                   Contrast:env_distance +  Contrast:logHPDdiff +
                   Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_14)

AIC(Model_14)

Like_ratio <- anova(Model_13, Model_14)


LRT_dist6 <- Permuted_model_simplification(Permuted_data,model1 = Model_13, model2 = Model_14)

ninety_five <- quantile(LRT_dist6, 0.95)


Like_ratio[2,"Chisq"] > ninety_five

### The liklihood ratio s not sigificantly greater than the distribution therefore the simplified model can proceeed


Anova(Model_14, type = "II")


### so lets remove the interaction between Contrast and logdist

Model_15 <- lmer(logitOver ~ Contrast + logdist + env_distance + logHPDdiff + 
                   S2RD1K + RD1Kdiff  + 
                    Contrast:logHPDdiff +
                   Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_1)

AIC(Model_15)

Like_ratio <- anova(Model_14, Model_15)


LRT_dist7 <- Permuted_model_simplification(Permuted_data,model1 = Model_14, model2 = Model_15)

ninety_five <- quantile(LRT_dist7, 0.95)


Like_ratio[2,"Chisq"] > ninety_five


### The liklihood ratio s not sigificantly greater than the distribution therefore the simplified model can proceeed


Anova(Model_15, type = "II")


### so lets remove env_dist

Model_16 <- lmer(logitOver ~ Contrast + logdist + logHPDdiff + 
                   S2RD1K + RD1Kdiff  + 
                   Contrast:logHPDdiff +
                   Contrast:RD1Kdiff  +
                   (1|SS),
                 data = Overlap_data, weights = min_site_spp )

summary(Model_16)

AIC(Model_16)

Like_ratio <- anova(Model_15, Model_16)


LRT_dist8 <- Permuted_model_simplification(Permuted_data,model1 = Model_15, model2 = Model_16)

ninety_five <- quantile(LRT_dist8, 0.95)


Like_ratio[2,"Chisq"] > ninety_five

### The liklihood ratio is sigificantly greater than the distribution therefore the previous model is retained!


Anova(Model_15, type = "II") ### Model looks good but still want to have a look if any other variables can be dropped

summary(Model_15)

ggResidpanel::resid_panel(Model_15)


saveRDS(Model_15, "Functional_Intactness_Index/Overlap_GLMM.rds")



