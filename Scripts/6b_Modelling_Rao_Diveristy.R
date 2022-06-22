rm(list = ls())

require(phyr)
require(tidyverse)
require(raster)
require(ggridges)
require(lme4)
require(car)
require(robustlmm)

## Load in the datasets

PREDICTS_site <- readRDS("../Functional_Intactness_Index/Outputs/PREDICTS_Site_Rao.rds")
PREDICTS_abundance <- readRDS("Outputs/Rao_abundance_data.rds")



hist(PREDICTS_site$Bias_Rao, breaks = 20)
hist(PREDICTS_site$Unbias_Rao, breaks = 20)

plot(PREDICTS_site$Unbias_Rao ~ PREDICTS_site$Bias_Rao)


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
                                          LandUse_Intensity = ifelse(grepl(LandUse_Intensity, pattern = "Urban"),
                                                                     "Urban",
                                                                     paste(LandUse_Intensity)),
                                          LandUse_Intensity = ifelse(grepl(LandUse_Intensity, pattern = "Cropland"),
                                                                     "Cropland",
                                                                     paste(LandUse_Intensity)),
                                          LandUse_Intensity = relevel(factor(LandUse_Intensity), 
                                                                      ref = "Primary_Minimal use"))

table(PREDICTS_site$LandUse_Intensity)

########## Going to rename for some ease of outputs
PREDICTS_site <- PREDICTS_site %>% dplyr::rename(LUI = LandUse_Intensity)
levels(PREDICTS_site$LUI) <- c("PriMin", "Cropland", "PasIn", "PasLig", 
                               "PasMin","PriIn","PriLig","SecIn","SecLig","SecMin","Urban")

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

Model_differences
### Model improves the most when I removed the interaction between LUI:RD50k

### Next round of removals -- 

## fixed effect RD50k

Rao_Model_9b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD1k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)




# LUI:RD1k

Rao_Model_10b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                                        LUI:logHPD +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)



# LUI:logHPD

Rao_Model_11b <- phyr::communityPGLMM(Unbias_Rao ~ LUI + logHPD  + RD50k + RD1k + CNTRLlogHPD +
                                        LUI:RD1k +
                                        (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)





Model_differences <- data.frame(RD50k = Rao_Model_9b$DIC - Rao_Model_7b$DIC,
                                RD1k = Rao_Model_10b$DIC - Rao_Model_7b$DIC,
                                logHPD = Rao_Model_11b$DIC - Rao_Model_7b$DIC)   

Model_differences

#### Model DIC is reduced in all but most when removing the interaction between LUI:RD1k so model 10b is best 

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


Robust_mod
###### so lets have a look at whats going on.

utils::sessionInfo()

### Functional diversity is significantly lower in Cropland vs Primary minimal habitat

phyr::fixef(Rao_Model_14b)


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#-----------------------------------------------------------
#----------------------------------------------------------


new_data <- phyr::fixef(Rao_Model_14b) 
new_data <- new_data[-(12:15),]

intercept <- new_data$Value[1]

alpha2 <- new_data %>% mutate(pred = Value + intercept)
alpha2$pred[1] <- alpha2$pred[1] - alpha2$Value[1] 



rownames(alpha2) <- c("PriMin", "Cropland", "PasIn", "PasLig", "PasMin", "PriIn", "PriLig","SecIn","SecLig","SecMin","Urban")

alpha2 <- alpha2 %>% dplyr::mutate(landuse = c("PRM","CRP","PAI","PAL","PAM","PRI","PRL","SEI","SEL","SEM","URB"),
                                   intensity = c("Minimal","All","Intense","Light","Minimal","Intense","Light","Intense","Light", "Minimal", "All"),
                                   type = c("primary","crop","pasture","pasture","pasture","primary","primary","secondary","secondary","secondary","urban"),
                                   Gow.Rao = pred)



