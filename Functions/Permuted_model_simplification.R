require(tidyverse)
require(car)
require(lme4)


#### function to generate the LR distribution across the 1000 datasets

##### Liklihood ratio function

Permuted_model_simplification <- function(data, model1, remove){
  
  remove_variables <- function(data,model1,variable){ 
    
    
    ##### make the variable to drop as formula 
    formula <- as.formula(paste("~.-",variable,sep = ""))
    
    #### this drops the variable from the model 
    model2 <- update(model1,formula)
    
    
    #### create blank object to then hold LRT of the random comparisons of models 
    LRT_dist <- c()
    
    #### i in number of randomisations
    
    for(i in 1:length(data)){
      
      
      ### "Full" model
      mod1 <- lmer(model1@call, data = data[[i]], REML = FALSE)
      
      ### reduced model 
      mod2 <- lmer(model2@call, data = data[[i]], REML = FALSE)
      
      ### calculate LRT
      LRT <- anova(mod1,mod2)
      
      #### extract LR
      LRT <- LRT[which(!is.na(LRT$Chisq)),"Chisq"]
      
      ## add to distribution
      LRT_dist <- rbind(LRT_dist,LRT)
      
    }
    
    
    ### same again but with the observed data
    mod_LRT <- anova(model1, model2)
    ChiSq <- mod_LRT[2,"Chisq"]
    
    ## get the 95 percentile of the distribution
    
    dist_quant <- quantile(LRT_dist, 0.95)
    
    ### if the LR of the observed model is less than the 95 percentile then the variable can be dropped from the model as the model is not significantly 
    ### impacted by the removal of the variable 
    
    DROP <- ChiSq < dist_quant
    
    #### now to see which variable has the least impact on the model and can therefore be dropped  
    
    percentile <- 0.01
    test <- TRUE
    while(test & percentile < 1.01){
      dq <- quantile(LRT_dist, percentile)
      test <- ChiSq > dq
      percentile <- percentile + 0.01
    }
    perecentile <- percentile - 0.01
    
    res <- data.frame(DROP = DROP, Percentile = perecentile)
    rownames(res) <- paste(variable)
    
    return(res)
  }
  
  ### cycle through the designated variables to drop
  
  frame <- data.frame(DROP = NULL, Percentile = NULL)
  for(var in remove){
    mod <- remove_variables(data = data, model1 = model1, variable = var) 
    frame <- rbind(frame,mod)
  }
  
  ### return frame indicating which variables to drop 
  return(frame)
}