PREDICTS_site <- PREDICTS_site %>% droplevels()




Rao_Model_1b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB), data = PREDICTS_site)




#### lets see what random effect structure gives the best AIC 

### adding random slopes 1) LUI within study 


Rao_Model_2b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB) + (LUI|SS), data = PREDICTS_site)


#2) logHPD within study 



Rao_Model_3b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB) + (logHPD|SS), data = PREDICTS_site)


#3) Road Density_1km within study  


Rao_Model_4b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB) + (RD1k|SS), data = PREDICTS_site)


# 4) RD50km within study 


Rao_Model_5b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                                       LUI:logHPD + LUI:RD1k + LUI:RD50k +
                                       (1|SS) + (1|SSB) + (RD50k|SS), data = PREDICTS_site)



Random_DIC <- data.frame(mod1 = AIC(Rao_Model_1b), mod2 = AIC(Rao_Model_2b),mod3 = AIC(Rao_Model_3b),
                         mod4 = AIC(Rao_Model_4b), mod5 = AIC(Rao_Model_5b))

Anova(Rao_Model_1b)


## remove LUI:RD50k

Rao_Model_6b <-  lmer(Unbias_Rao ~ LUI + logHPD + RD1k + RD50k + CNTRLlogHPD +
                         LUI:logHPD + LUI:RD1k +
                         (1|SS) + (1|SSB), data = PREDICTS_site)

anova(Rao_Model_6b, Rao_Model_1b)


### proceed

Anova(Rao_Model_6b)

Rao_Model_7b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                       LUI:logHPD + LUI:RD1k +
                       (1|SS) + (1|SSB), data = PREDICTS_site)

anova(Rao_Model_6b, Rao_Model_7b)


Anova(Rao_Model_7b)

Rao_Model_8b <- lmer(Unbias_Rao ~ LUI + logHPD + RD1k + CNTRLlogHPD +
                       LUI:logHPD +
                       (1|SS) + (1|SSB), data = PREDICTS_site)

anova(Rao_Model_7b, Rao_Model_8b)


### proceed 

Anova(Rao_Model_8b)

summary(Rao_Model_8b)

resid_panel(Rao_Model_8b)
