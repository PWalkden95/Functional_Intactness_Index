rm(list = ls())

require(tidyverse)
require(car)
require(lme4)
require(ggResidpanel)



#### Load in Similarity data


Similarity_data <- readRDS("Functional_Intactness_Index/Outputs/similarity_data.rds")
TPD_data <- readRDS("Functional_Intactness_Index/Outputs/Trait_Prob_den.rds")


### Want to work out human population density at 


TPD_data <- TPD_data %>% dplyr::mutate(logitOver = car::logit(TPD_Overlap, adjust = 0.001, percents = FALSE)) %>%
  dplyr::mutate(logdist = log(distance + 1),
                rt3env = env_distance^(1/3),
                sqrtS2RD1K = S2RD1K^(1/2),
                sqrtS2RD50K = S2RD50K^(1/2),
                RD1Kdiff = sqrtS2RD1K - (S1RD1K^(1/2)),
                RD50Kdiff = sqrtS2RD50K - (S1RD50K^(1/2)))




table(TPD_data$Contrast)



TPD_data <- TPD_data %>% dplyr::mutate(Cont = ifelse(grepl(Contrast, pattern = "Plantation forest"), 
                                "PriMin-SecLig",
                                paste(Contrast)),
                                Cont = ifelse(grepl(Cont, pattern = "Primary_Light")|grepl(Cont, pattern = "Primary_Intense"),
                                "PriMin-Primary",
                                paste(Cont)),
                                Cont = ifelse(grepl(Cont, pattern = "Cropland"),
                                              "PriMin-Cropland", 
                                              paste(Cont)),
                              Cont = ifelse(grepl(Cont, pattern = "Secondary Vegetation_Light use"),
                                            "PriMin-SecLig", 
                                            paste(Cont)),
                              Cont  = relevel(factor(Cont), ref = "Primary_Minimal use-Primary_Minimal use"))


levels(TPD_data$Cont) <- c("PriMin-PriMin","PriMin-SecMin", "PriMin-UrbMin","PriMin-Cropland", "PriMin-Primary","PriMin-SecLig")


TPD_data$rt3env <- scale(TPD_data$rt3env)
TPD_data$s2logHPD <- scale(TPD_data$S2logHPD)
TPD_data$logHPDdiff <- scale(TPD_data$logHPDdiff)
TPD_data$logdist <- scale(TPD_data$logdist)
TPD_data$sqrtS2RD1K <- scale(TPD_data$sqrtS2RD1K)
TPD_data$sqrtS2RD50K <- scale(TPD_data$sqrtS2RD50K)
TPD_data$RD1Kdiff <- scale(TPD_data$RD1Kdiff)
TPD_data$RD50Kdiff <- scale(TPD_data$RD50Kdiff)


source("https://highstat.com/Books/Book2/HighstatLibV10.R")

corvif(TPD_data[,c("Cont", "logdist", "rt3env","S2logHPD", "logHPDdiff", "sqrtS2RD1K", "sqrtS2RD50K", "RD1Kdiff", "RD50Kdiff")])

### some colinearity problems lets drop sqrtS2RD50k

corvif(TPD_data[,c("Cont", "logdist", "rt3env", "S2logHPD","logHPDdiff", "sqrtS2RD1K", "RD1Kdiff", "RD50Kdiff")])




################################
##### Modelling ################
################################



Model_1 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS),
                data = TPD_data)

summary(Model_1)

AIC(Model_1)

#### random slopes

## logdist

Model_2 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + logdist|SS),
                data = TPD_data )

summary(Model_2)

AIC(Model_2)

## rt3env

Model_3 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + rt3env|SS),
                data = TPD_data)

summary(Model_3)

AIC(Model_3)


## logHPDdiff

Model_4 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + logHPDdiff|SS),
                data = TPD_data )

summary(Model_4)

AIC(Model_4)

## sqrtS2RD1K

Model_5 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + sqrtS2RD1K|SS),
                data = TPD_data)

summary(Model_5)

AIC(Model_5)

## S2logHPD

Model_6 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + S2logHPD|SS),
                data = TPD_data )

summary(Model_6)

AIC(Model_6)

## RoadDdiff1k

Model_7 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env + Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + RD1Kdiff|SS),
                data = TPD_data)

summary(Model_7)

AIC(Model_7)

## RoadDdiff50k

Model_8 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + RD50Kdiff|SS),
                data = TPD_data)

summary(Model_8)

AIC(Model_8)


## Cont 

Model_9 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont:sqrtS2RD1K + Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + Cont|SS),
                data = TPD_data)

summary(Model_9)

AIC(Model_9)


MOD_AIC <- data.frame(Mod1 = AIC(Model_1),Mod2 = AIC(Model_2),Mod3 = AIC(Model_3),
                      Mod4 = AIC(Model_4),Mod5 = AIC(Model_5),Mod6 = AIC(Model_6),
                      Mod7 = AIC(Model_7),Mod8 = AIC(Model_8), Mod9 = AIC(Model_9))







Permuted_data <- rep(list(NA),1000)

for(i in 1:1000){
  
  sample_data <-c()
  
  for(study in unique(TPD_data$SS)){
    data <- TPD_data %>% filter(SS == study)
    
    data$logitOver <- data[sample(NROW(data)),"logitOver"]
    
    sample_data <- rbind(sample_data,data)
    
  }
  
  Permuted_data[[i]] <- sample_data
  
}

#### function to generate the LR distribution across the 100 datasets

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



### so lets look at our best maximal model

Anova(Model_1, type = "II")

### with the lowest p-value we should remove the interaction between Contrast and S2/rd1K

Model_10 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + 
                  Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                  Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) ,
                data = TPD_data)

summary(Model_10)

AIC(Model_10)


Like_ratio <- anova(Model_1, Model_10)


LRT_dist1 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, model2 = Model_10)

ninety_five <- quantile(LRT_dist1,0.95)

Like_ratio[2,"Chisq"] > ninety_five

##### LRT of models using the observed data is not significantly different from the distribution of LRT using the permuted data therefore the
##### interaction can be dropped from the model.

Anova(Model_10,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove  sqrtS2RD1k

Model_11 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                    S2logHPD + RD1Kdiff + RD50Kdiff + 
                   Cont:logdist + Cont:rt3env +  Cont:logHPDdiff +
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_11)

AIC(Model_11)



Like_ratio <- anova(Model_10, Model_11)


LRT_dist2 <- Permuted_model_simplification(Permuted_data,model1 = Model_10, model2 = Model_11)

ninety_five <- quantile(LRT_dist2,0.95)

Like_ratio[2,"Chisq"] > ninety_five


##### LRT of models using the observed data is not significantly different from the distribution of LRT using the permuted data therefore the
##### interaction can be dropped from the model.

Anova(Model_11,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove Cont:logHPDdiff

Model_12 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff + RD50Kdiff + 
                   Cont:logdist + Cont:rt3env +
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_12)

AIC(Model_12)



Like_ratio <- anova(Model_11, Model_12)


LRT_dist3 <- Permuted_model_simplification(Permuted_data,model1 = Model_11, model2 = Model_12)

ninety_five <- quantile(LRT_dist3,0.95)

Like_ratio[2,"Chisq"] > ninety_five

########## LRT test comes back negative and therefore the reduced model is not significantly reduced 

Anova(Model_12,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove interaction between contrast:rt3env

Model_13 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff + RD50Kdiff + 
                   Cont:logdist + 
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:RD50Kdiff +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_13)

AIC(Model_13)


Like_ratio <- anova(Model_12, Model_13)


LRT_dist4 <- Permuted_model_simplification(Permuted_data,model1 = Model_12, model2 = Model_13)

ninety_five <- quantile(LRT_dist4,0.95)

Like_ratio[2,"Chisq"] > ninety_five

########## LRT test comes back negative and therefore the reduced model is not significantly reduced 

Anova(Model_13,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove interaction between contrast:RD50Kdiff

Model_14 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff + RD50Kdiff + 
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:logdist +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_14)

AIC(Model_14)


Like_ratio <- anova(Model_13, Model_14)


LRT_dist5 <- Permuted_model_simplification(Permuted_data,model1 = Model_13, model2 = Model_14)

ninety_five <- quantile(LRT_dist5,0.95)

Like_ratio[2,"Chisq"] > ninety_five

########## LRT test comes back negative and therefore the reduced model is not significantly reduced 

Anova(Model_14,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove RD50Kdiff

Model_15 <- lmer(logitOver ~ Cont  + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff +  logdist +
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff + Cont:logdist +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_15)

AIC(Model_15)


Like_ratio <- anova(Model_14, Model_15)


LRT_dist6 <- Permuted_model_simplification(Permuted_data,model1 = Model_14, model2 = Model_15)

ninety_five <- quantile(LRT_dist6,0.95)

Like_ratio[2,"Chisq"] > ninety_five

########## LRT test comes back negative and therefore the reduced model is not significantly reduced 

Anova(Model_15,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove rt3env

Model_16 <- lmer(logitOver ~ Cont  + logHPDdiff + 
                   S2logHPD + RD1Kdiff + logdist + 
                    Cont:S2logHPD + Cont:logdist + Cont:RD1Kdiff +
                   (1|SS) ,
                 data = TPD_data)

summary(Model_16)

AIC(Model_16)


Like_ratio <- anova(Model_15, Model_16)


LRT_dist7 <- Permuted_model_simplification(Permuted_data,model1 = Model_15, model2 = Model_16)

ninety_five <- quantile(LRT_dist7,0.95)

Like_ratio[2,"Chisq"] > ninety_five


########## LRT test comes back positive therefore the model loses a significant amount of explanatiry power and shuold not continue 

Anova(Model_15,type = "II")

## Most are significant and fixed effects that aren't their interactions are therefore I will try to remove interaction between Cont:RD50kdiff

Model_17 <- lmer(logitOver ~ Cont  + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff +  logdist +
                   Cont: Cont:S2logHPD + Cont:RD1Kdiff +
                   (1|SS) ,
                 data = TPD_data)
summary(Model_17)

AIC(Model_17)


Like_ratio <- anova(Model_15, Model_17)


LRT_dist8 <- Permuted_model_simplification(Permuted_data,model1 = Model_15, model2 = Model_17)

ninety_five <- quantile(LRT_dist8,0.95)

Like_ratio[2,"Chisq"] > ninety_five

########## LRT test comes back negative and therefore the reduced model is not significantly reduced 

Anova(Model_15,type = "II")

summary(Model_15)

test <- lmer(logitOver ~ Cont  + rt3env + logHPDdiff + 
                   S2logHPD + RD1Kdiff +  logdist +
                   (1|SS) ,
                 data = TPD_data)
summary(test)

