PREDICTS_site <- PREDICTS_site %>% droplevels()




Rao_Model_1b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km + LUI:Road_Density_50km +
                                       (1|SS) + (1|SSB), data = PREDICTS_site, bayes = TRUE)


summary(Rao_Model_1b)

#### lets see what random effect structure gives the best AIC 

### adding random slopes 1) LUI within study 


Rao_Model_2b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km + LUI:Road_Density_50km +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_2b)

#2) logHPD within study 



Rao_Model_3b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km + LUI:Road_Density_50km +
                                       (1|SS) + (1|SSB) + (logHPD|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_3b)

#3) Road Density_1km within study  


Rao_Model_4b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km + LUI:Road_Density_50km +
                                       (1|SS) + (1|SSB) + (Road_Density_1km|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_4b)

# 4) Road_Density_50kmm within study 


Rao_Model_5b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km + LUI:Road_Density_50km +
                                       (1|SS) + (1|SSB) + (Road_Density_50km|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_5b)




Random_DIC <- data.frame(mod1 = Rao_Model_1b$DIC, mod2 = Rao_Model_2b$DIC,mod3 = Rao_Model_3b$DIC,
                         mod4 = Rao_Model_4b$DIC, mod5 = Rao_Model_5b$DIC)



Rao_Model_6b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD + LUI:Road_Density_1km  +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_6b)

Rao_Model_6b$DIC - Rao_Model_2b$DIC

########
Rao_Model_7b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km +
                                       LUI:logHPD + LUI:Road_Density_1km  +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_7b)

Rao_Model_7b$DIC - Rao_Model_6b$DIC #### didnt reduce enough

#### LUI: RD1k

Rao_Model_8b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       LUI:logHPD  +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_8b)

Rao_Model_8b$DIC - Rao_Model_6b$DIC  ## Much reduce

### 

Rao_Model_9b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km +
                                       LUI:logHPD  +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_9b)

Rao_Model_9b$DIC - Rao_Model_8b$DIC  ## nope reduce

## RD1k

Rao_Model_10b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_50km +
                                       LUI:logHPD  +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_10b)

Rao_Model_10b$DIC - Rao_Model_8b$DIC  ## Much increase

### LUI:loghpd

Rao_Model_11b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km + Road_Density_50km +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_11b)

Rao_Model_11b$DIC - Rao_Model_8b$DIC  ## Much reduce

####      RD50k


Rao_Model_12b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_1km +
                                        (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_12b)

Rao_Model_12b$DIC - Rao_Model_11b$DIC  ## not enough



#### RD1k

Rao_Model_13b <- phyr::communityPGLMM(Bias_Rao ~ LUI + logHPD + Road_Density_50km +
                                        (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_13b)

Rao_Model_13b$DIC - Rao_Model_11b$DIC  ## Much increase

##logHPD

Rao_Model_14b <- phyr::communityPGLMM(Bias_Rao ~ LUI  + Road_Density_1km + Road_Density_50km +
                                        (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site, bayes = TRUE)

summary(Rao_Model_14b)

Rao_Model_14b$DIC - Rao_Model_11b$DIC  ## nope modell 11b best again


summary(Rao_Model_11b)



