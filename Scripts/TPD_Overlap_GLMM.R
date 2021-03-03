rm(list = ls())

require(tidyverse)
require(car)
require(lme4)
require(ggResidpanel)



#### Load in Similarity data


TPD_data <- readRDS("Outputs/Trait_Prob_den.rds")


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
table(TPD_data$Cont)
TPD_data <- TPD_data[-118,]

TPD_data$rt3env <- scale(TPD_data$rt3env)
TPD_data$s2logHPD <- scale(TPD_data$S2logHPD)
TPD_data$logHPDdiff <- scale(TPD_data$logHPDdiff)
TPD_data$logdist <- scale(TPD_data$logdist)
TPD_data$sqrtS2RD1K <- scale(TPD_data$sqrtS2RD1K)
TPD_data$sqrtS2RD50K <- scale(TPD_data$sqrtS2RD50K)
TPD_data$RD1Kdiff <- scale(TPD_data$RD1Kdiff)
TPD_data$RD50Kdiff <- scale(TPD_data$RD50Kdiff)
TPD_data$CNTRLlogHPD <- scale(TPD_data$CNTRLlogHPD)


##### Because data is being compared multiple times - the same primary minimal site compared multiple times you are not able to assess
#### colinearity the usual way using VIF or otherwise so backwards stepwise model selection is going to take a while

##BUT first lets assess the optimal random effect structure in the maximal model


################################
##### Modelling ################
################################



Model_1 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff +
                  + Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS),
                data = TPD_data)


#### random slopes

## logdist

Model_2 <- lmer(logitOver ~ Cont + logdist + rt3env + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff + 
                    Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + logdist|SS),
                data = TPD_data )



## logHPDdiff

Model_3 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff + 
                    Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + logHPDdiff|SS),
                data = TPD_data )


## RoadDdiff1k

Model_4 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff + 
                    Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + RD1Kdiff|SS),
                data = TPD_data)


## RoadDdiff50k

Model_5 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff + 
                  Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + RD50Kdiff|SS),
                data = TPD_data)


## Cont 

Model_6 <- lmer(logitOver ~ Cont + logdist + rt3env  + logHPDdiff + 
                  sqrtS2RD1K + S2logHPD + RD1Kdiff + RD50Kdiff + sqrtS2RD50K + CNTRLlogHPD +
                  Cont:logdist +  Cont:logHPDdiff + 
                    Cont:RD1Kdiff + Cont:RD50Kdiff +
                  (1|SS) + (1 + Cont|SS),
                data = TPD_data)


MOD_AIC <- data.frame(Mod1 = AIC(Model_1),Mod2 = AIC(Model_2),Mod3 = AIC(Model_3),
                      Mod4 = AIC(Model_4),Mod5 = AIC(Model_5),Mod6 = AIC(Model_6))







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

#### function to generate the LR distribution across the 1000 datasets

##### Liklihood ratio function

Permuted_model_simplification <- function(data, model1, remove){
  
  formula <- as.formula(paste("~.-",remove,sep = ""))
  model2 <- update(model1,formula)
  
  LRT_dist <- c()
  
  for(i in 1:length(data)){
    
    mod1 <- lmer(model1@call, data = data[[i]], REML = FALSE)
    mod2 <- lmer(model2@call, data = data[[i]], REML = FALSE)
    
    LRT <- anova(mod1,mod2)
    LRT <- LRT[which(!is.na(LRT$Chisq)),"Chisq"]
    
    LRT_dist <- rbind(LRT_dist,LRT)
    
  }
  
  mod_LRT <- anova(model1, model2)
  ChiSq <- mod_LRT[2,"Chisq"]
  
  dist_quant <- quantile(LRT_dist, 0.95)
  
  DROP <- ChiSq < dist_quant
  
  percentile <- 0.01
  test <- TRUE
  while(test & percentile < 1.01){
    dq <- quantile(LRT_dist, percentile)
    test <- ChiSq > dq
    percentile <- percentile + 0.01
  }
perecentile <- percentile - 0.01
  
res <- data.frame(DROP = DROP, Percentile = perecentile)
rownames(res) <- paste(remove)

return(res)
}
### so lets look at our best maximal model

Anova(Model_1, type = "II")
### Because there may be some colinearity issues in the explanatory variables it is not reliable to pick the term to remove in the simplification
### using the highest p-value, therefore I will have to try removing all possibilities and proceeding with the model that imroves the most

###### int Cont:RD50Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, remove = "Cont:RD50Kdiff")

###### int Cont:RD1Kdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, remove = "Cont:logHPDdiff")


###### int Cont:logdist
mod_sim4 <- Permuted_model_simplification(Permuted_data,model1 = Model_1, remove = "Cont:logdist")


Model_simp <- rbind(mod_sim1,mod_sim2,mod_sim3,mod_sim4)


### Could drop all but Cont:logdist and Cont:S2logHPD shows the lowest probability of significantly lowering the explanatory
### power of the model

Model_7 <- update(Model_1, ~.-Cont:RD50Kdiff)

Anova(Model_7, type = "II")

###### RD50Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_7, remove = "RD50Kdiff")

###### int Cont:RD1Kdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_7, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_7, remove = "Cont:logHPDdiff")

###### int Cont:logdist
mod_sim4 <- Permuted_model_simplification(Permuted_data,model1 = Model_7, remove = "Cont:logdist")


Model_simp2 <- rbind(mod_sim1,mod_sim2,mod_sim3,mod_sim4)

### remove Cont:sqrtS2RD1k

Model_8 <- update(Model_7, ~.-RD50Kdiff)


Anova(Model_8, type = "III")


###### int Cont:RD1Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_8, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_8, remove = "Cont:logHPDdiff")

###### int Cont:logdist
mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_8, remove = "Cont:logdist")

##### S2 RD50k

mod_sim4 <- Permuted_model_simplification(Permuted_data,model1 = Model_8, remove = "sqrtS2RD50K")


Model_simp3 <- rbind(mod_sim1,mod_sim2,mod_sim3,mod_sim4)



### remove interaction Cont:RD50Kdiff

Model_9 <- update(Model_8, ~.-sqrtS2RD50K)


Anova(Model_9, type = "III")

###### int Cont:RD1Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_9, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_9, remove = "Cont:logHPDdiff")

###### int Cont:logdist
mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_9, remove = "Cont:logdist")

##### S2 RD50k
Model_simp4 <- rbind(mod_sim1,mod_sim2,mod_sim3)


#### Can remove RD50Kdiff

Model_10 <- update(Model_9, ~.-Cont:logHPDdiff)


Anova(Model_10, type = "III")

####### int Cont:RD1Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_10, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_10, remove = "logHPDdiff")

###### int Cont:logdist
mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_10, remove = "Cont:logdist")

##### S2 RD50k
Model_simp5 <- rbind(mod_sim1,mod_sim2,mod_sim3)

## Can't remove either of the S2 effects while the diff variable is in the model therefore Cont:rt3env interaction is removed

Model_11 <- update(Model_10, ~.-logHPDdiff)

Anova(Model_11, type = "III")

####### int Cont:RD1Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_11, remove = "Cont:RD1Kdiff")

###### int Cont:logHPDdiff

mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_11, remove = "S2logHPD")

###### int Cont:logdist
mod_sim3 <- Permuted_model_simplification(Permuted_data,model1 = Model_11, remove = "Cont:logdist")

##### S2 RD50k
Model_simp6 <- rbind(mod_sim1,mod_sim2,mod_sim3)




Model_12 <- update(Model_11, ~.-S2logHPD)



Anova(Model_12, type = "III")

####### int Cont:RD1Kdiff

mod_sim1 <- Permuted_model_simplification(Permuted_data,model1 = Model_12, remove = "Cont:RD1Kdiff")


###### int Cont:logdist
mod_sim2 <- Permuted_model_simplification(Permuted_data,model1 = Model_12, remove = "Cont:logdist")


Model_simp7 <- rbind(mod_sim1,mod_sim2)


summary(Model_12)


