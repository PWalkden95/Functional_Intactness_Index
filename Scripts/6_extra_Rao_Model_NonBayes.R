rm(list = ls())

require(phyr)
require(tidyverse)
require(raster)
require(ggridges)
require(lme4)
require(car)
require(robustlmm)

## Load in the datasets

PREDICTS_site <- readRDS("Outputs/PREDICTS_Site_Rao.rds")
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


#### lets see what random effect structure gives the best AIC 

### adding random slopes 1) LUI within study 

Rao_Model_2 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + LUI|SS), data = PREDICTS_site)

summary(Rao_Model_2)


#2) logHPD within study 

Rao_Model_3 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + logHPD|SS), data = PREDICTS_site)

summary(Rao_Model_3)


#3) Road Density_1km within study  

Rao_Model_4 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + RD1k|SS), data = PREDICTS_site)

summary(Rao_Model_4)


# 4) RD50km within study 

Rao_Model_5 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD1k + LUI:RD50k +
                      (1|SS) + (1|SSB) + (1 + RD50k|SS), data = PREDICTS_site)

summary(Rao_Model_5)




Random_AIC <- data.frame(mod1 = AIC(Rao_Model_1), mod2 = AIC(Rao_Model_2),mod3 = AIC(Rao_Model_3),
                         mod4 = AIC(Rao_Model_4), mod5 = AIC(Rao_Model_5))




### Model 3 gives the best result which has logHPD as a random slope within study 
## so now to test the best fixed effect structure

car::Anova(Rao_Model_1, type = "III") ## desnity of roads 50k is the least significant interaction so will test its exclusion




### GLMM fails to converge removing the interaction between Road_density50k and landuse so I will continue with Bayes GLMM (INLA)

#### Remove each interaction to see if thsis improves the model at all. 

#LUI:RD50k

Rao_Model_6 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site)



#LUI:RD1k

Rao_Model_7 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD50k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site)


# LUI:logHPD

Rao_Model_8 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:RD50k + LUI:RD1k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site)


AIC(Rao_Model_6)

Model_differences <- data.frame(RD50k = AIC(Rao_Model_6) - AIC(Rao_Model_1),
                                RD1k = AIC(Rao_Model_7) - AIC(Rao_Model_1),
                                logHPD = AIC(Rao_Model_8) - AIC(Rao_Model_1))   

Model_differences


#### Remove each interaction to see if thsis improves the model at all. 

#LUI:RD50k

Rao_Model_9 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                      LUI:logHPD  +
                      (1|SS) + (1|SSB), data = PREDICTS_site)



#RD1k

Rao_Model_10 <- lmer(Unbias_Rao ~ LUI + logHPD + RD50k + CNTRLlogHPD +
                      LUI:logHPD + LUI:RD50k +
                      (1|SS) + (1|SSB), data = PREDICTS_site)


# LUI:logHPD

Rao_Model_11 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                       LUI:RD50k +
                      (1|SS) + (1|SSB), data = PREDICTS_site)




Model_differences <- data.frame(RD50k = AIC(Rao_Model_9) - AIC(Rao_Model_7),
                                RD1k = AIC(Rao_Model_10) - AIC(Rao_Model_7),
                                logHPD = AIC(Rao_Model_11) - AIC(Rao_Model_7))   

Model_differences


#### Remove each interaction to see if thsis improves the model at all. 

#RD50k

Rao_Model_12 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                      LUI:logHPD  +
                      (1|SS) + (1|SSB), data = PREDICTS_site)



#RD1k

Rao_Model_13 <- lmer(Unbias_Rao ~ LUI + logHPD + RD50k + CNTRLlogHPD +
                       LUI:logHPD  +
                       (1|SS) + (1|SSB), data = PREDICTS_site)


# LUI:logHPD

Rao_Model_14 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)




Model_differences <- data.frame(RD50k = AIC(Rao_Model_12) - AIC(Rao_Model_9),
                                RD1k = AIC(Rao_Model_13) - AIC(Rao_Model_9),
                                logHPD = AIC(Rao_Model_14) - AIC(Rao_Model_9))   

Model_differences


#### Remove each interaction to see if thsis improves the model at all. 

#RD50k

Rao_Model_15 <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)



#RD1k

Rao_Model_16 <- lmer(Unbias_Rao ~ LUI + logHPD + RD50k + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)


# logHPD

Rao_Model_17 <- lmer(Unbias_Rao ~ LUI + RD1k + RD50k + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)




Model_differences <- data.frame(RD50k = AIC(Rao_Model_15) - AIC(Rao_Model_14),
                                RD1k = AIC(Rao_Model_16) - AIC(Rao_Model_14),
                                logHPD = AIC(Rao_Model_17) - AIC(Rao_Model_14))   

Model_differences


#### Remove each interaction to see if thsis improves the model at all. 

#RD1k

Rao_Model_18 <- lmer(Unbias_Rao ~ LUI + logHPD + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)


# logHPD

Rao_Model_19 <- lmer(Unbias_Rao ~ LUI + RD1k + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)




Model_differences <- data.frame(RD1k = AIC(Rao_Model_18) - AIC(Rao_Model_15),
                                logHPD = AIC(Rao_Model_19) - AIC(Rao_Model_15))   

Model_differences


##### remove the final other pressure#


#lpgHPD

Rao_Model_20 <- lmer(Unbias_Rao ~ LUI + CNTRLlogHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)

Model_differences <- AIC(Rao_Model_20) - AIC(Rao_Model_18)

Model_differences


Rao_Model_20


require(ggResidpanel)


ggResidpanel::resid_panel(Rao_Model_20)


summary(Rao_Model_20)


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

