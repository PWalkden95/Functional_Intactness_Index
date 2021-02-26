rm(list = ls())

require(tidyverse)
require(geosphere)
require(hypervolume)
require(doParallel)
require(foreach)
require(TPD)

Similarity_data <- read_rds("Functional_Intactness_Index/similarity_data.rds")
Similarity_data <- Similarity_data %>% filter(Diversity_metric == "abundance") %>% droplevels()
PC_Scores <- read_rds("Functional_Intactness_Index/PC_Scores.rds")
Full_PC_Scores <- read_rds(file ="Functional_Intactness_Index/Full_PC_scores.rds")


remove_dupl_comparisons <- function(data){
  
  data$drop <- NA
  data$comparison <- paste(data$site1,data$site2, sep = "_")
  data[1,"drop"] <- FALSE
  
  
  if(NROW(data)> 1){
    
  for(i in 2:NROW(data)){
    
    site1 <- data[i,"site1"]
    site2 <- data[i,"site2"]
    
    match_1 <- as.numeric(grep(data[1:i-1,"comparison"], pattern = site1))
    match_2 <- as.numeric(grep(data[1:i-1,"comparison"], pattern = site2))
    
    if(any(match_1 %in% match_2)){
      data[i,"drop"] <- TRUE
    } else {
      data[i,"drop"] <- FALSE
    }
    
  }
  
  }
  
  if(any(data$drop)){
    data <- data[-which(data$drop == TRUE & data$Contrast == "Primary_Minimal use-Primary_Minimal use"),-which(colnames(data) == "comparison" | colnames(data) == "drop")]
  } else {
    data <- data[,-which(colnames(data) == "comparison" | colnames(data) == "drop")]
  }
  
  
  return(data)
}

##### #### ##########################
#### Function for effort check ######
####################################

sampling_effort_check <- function(site1,site2, trait_ranges){
  
  individuals_1 <- data.frame(individuals = rep(site1$Jetz_Name, site1$Measurement))
  individuals_2 <- data.frame(individuals = rep(site2$Jetz_Name, site2$Measurement))
  
  
  proportion <- seq(0.05,1, by = 0.05)
  
  Overlap_stats <- c()
  
  for(prop in proportion){
    
    
    length_1 <- as.numeric(NROW(individuals_1))
    sample_1 <- round(length_1 * prop)
    
    extract_1 <- sample_n(individuals_1, sample_1, replace = FALSE)
    
    
    spp_trait_1 <- data.frame(Jetz_Name = unique(extract_1$individuals)) %>% dplyr::left_join(PC_Scores, by = "Jetz_Name")
    rownames(spp_trait_1) <- spp_trait_1$Jetz_Name
    spp_trait_1 <- spp_trait_1[,-1]
    
    if(NROW(spp_trait_1) < 10) {
      hyper_1 <- NULL
      volume_1 <- 0
    } else {
      
      hyper_1 <- hypervolume(spp_trait_1, method = "svm")
      
      volume_1 <- hyper_1@Volume
    }
    
    
    ###
    
    length_2 <- as.numeric(NROW(individuals_2))
    sample_2 <- round(length_2 * prop)
    
    extract_2 <- sample_n(individuals_2, sample_2,replace = FALSE)
    
    
    spp_trait_2 <- data.frame(Jetz_Name = unique(extract_2$individuals)) %>% dplyr::left_join(PC_Scores, by = "Jetz_Name")
    rownames(spp_trait_2) <- spp_trait_2$Jetz_Name
    spp_trait_2 <- spp_trait_2[,-1]
    
    if(is.null(hyper_1)| NROW(spp_trait_2) < 10){
      Overlap <- NULL
      volume_2 <- 0
    } else {
      
      hyper_2 <- hypervolume(spp_trait_2, method = "svm")
      
      volume_2 <- hyper_2@Volume
      
      
      hyper_set <- hypervolume_set(hyper_1, hyper_2, check.memory = FALSE)
      
      
      Overlap <- hypervolume_overlap_statistics(hyper_set)
      
    }
    ####### TPD calc
    
    over_spp <- data.frame(Jetz_Name = unique(c(extract_1$individuals,extract_2$individuals)))
    
    spp_TPD <- over_spp %>% dplyr::left_join(Full_PC_Scores, by = "Jetz_Name") %>% na.omit
    
    T_D <- TPD::TPDs(spp_TPD[,1],spp_TPD[,c(2:4)], n_divisions = 66, alpha = 0.95,  trait_ranges = trait_ranges)
    
    
    
  
    S1_spp <- extract_1 %>% group_by(individuals) %>% count() %>% ungroup() %>% dplyr::mutate(RelativeAbundance = n/sum(n))
    S2_spp <- extract_2 %>% group_by(individuals) %>% count() %>% ungroup() %>% dplyr::mutate(RelativeAbundance = n/sum(n))
    
    
    
    community <- data.frame(Jetz_Name = unique(spp_TPD[,1])) %>% dplyr::left_join(S1_spp[,c("individuals", "RelativeAbundance")], by = c("Jetz_Name" = "individuals")) %>%
      dplyr::rename(site1_abun = RelativeAbundance) %>%
      dplyr::left_join(S2_spp[,c("individuals", "RelativeAbundance")], by = c("Jetz_Name" = "individuals")) %>%
      dplyr::rename(site2_abun = RelativeAbundance)
    
    community[is.na(community)] <- 0 
    rownames(community) <- community$Jetz_Name
    community <- community[,-1]
    
    TPD_over <- TPDc(T_D, t(community))
    
    TPD_diss <- dissim(TPD_over)
    
    #########
    
    prop_stats <- data.frame(site1 = paste(unique(as.character(site1$SSBS))), site2 = paste(unique(as.character(site2$SSBS))),Proportion = paste(prop),nspp1 = sample_1,nspp2 = sample_2 ,vol1 = volume_1, vol2 = volume_2,
                             overlap = ifelse(is.null(Overlap),0,Overlap[1]), TPD_overlap = 1 - TPD_diss$communities$dissimilarity[2])
    
    
    
    Overlap_stats <- rbind(Overlap_stats, prop_stats)
    
  }
  
  
  overplot <- ggplot(Overlap_stats, aes(x = proportion, y = overlap)) +
    geom_line() +
    labs(title = paste("Site1 Range: ",min(Overlap_stats$nspp1),"-",max(Overlap_stats$nspp1),"Site2 Range: ", min(Overlap_stats$nspp2),"-",max(Overlap_stats$nspp2)))
  
  TPD_overplot <- ggplot(Overlap_stats, aes(x = proportion, y = TPD_overlap)) +
    geom_line() +
    labs(title = paste("Site1 Range: ",min(Overlap_stats$nspp1),"-",max(Overlap_stats$nspp1),"Site2 Range: ", min(Overlap_stats$nspp2),"-",max(Overlap_stats$nspp2)))
  
  vol1 <- ggplot(Overlap_stats, aes(x = proportion, y = vol1)) +
    geom_line()+
    labs(title = paste("Site1 Range: ",min(Overlap_stats$nspp1),"-",max(Overlap_stats$nspp1),"Site2 Range: ", min(Overlap_stats$nspp2),"-",max(Overlap_stats$nspp2)))
  
  vol2 <- ggplot(Overlap_stats, aes(x = proportion, y = vol2)) +
    geom_line()+
    labs(title = paste("Site1 Range: ",min(Overlap_stats$nspp1),"-",max(Overlap_stats$nspp1),"Site2 Range: ", min(Overlap_stats$nspp2),"-",max(Overlap_stats$nspp2)))
  
  plots <- list(data = Overlap_stats, overplot = overplot, TPD_overplot = TPD_overplot,
                vol1 = vol1, vol2 = vol2)
  
  
  return(plots)
}






### What studies do we have in the dataset 

studies <- levels(Similarity_data$SS) 

#### empty data frame for results to go into 

trait_ranges <- list(c(min(Full_PC_Scores[,2]) -(0.3 * min(Full_PC_Scores[,2])),max(Full_PC_Scores[,2]) + (0.3 * max(Full_PC_Scores[,2]))),
                     c(min(Full_PC_Scores[,3]) -(0.3 * min(Full_PC_Scores[,3])),max(Full_PC_Scores[,3]) + (0.3 * max(Full_PC_Scores[,3]))),
                     c(min(Full_PC_Scores[,4]) -(0.3 * min(Full_PC_Scores[,4])),max(Full_PC_Scores[,4]) + (0.3 * max(Full_PC_Scores[,4]))))



registerDoParallel(cores = 5)

Effort_Check <- foreach(study = studies,
                        .combine = "c",
                        .packages = c("hypervolume","betapart","tidyverse","magrittr","geosphere","gower","raster", "TPD")) %dopar% {
                          
                          
                          ### data for study 
                          
                          
                          
                          
                          data <- Similarity_data %>% filter(SS == study) %>% droplevels()
                          
                          ### coordinates going to be used to calculate distance between sites
                          
                          
                          LandUse <- data.frame(data) %>% distinct(SSBS, LandUse_Intensity)
                          
                          
                          Primary_sites <- LandUse %>% filter(LandUse_Intensity == "Primary_Minimal use") %>% pull(SSBS) %>% as.character()
                          
                          #### all other sites 
                          
                          all_sites <- LandUse %>% pull(SSBS) %>% as.character()
                          
                          
                          
                          
                          ###### get the comparisons with the primary minimal site and all other sites within the study 
                          
                          site_comparisons <- expand.grid(Primary_sites, all_sites) %>% 
                            dplyr::rename(site1 = Var1, site2 = Var2) %>% dplyr::mutate(site1 = as.character(site1), site2 = as.character(site2)) %>%
                            filter(site1 != site2) %>%
                            left_join(LandUse, by = c("site1" = "SSBS")) %>% dplyr::rename(site1LUI = LandUse_Intensity) %>%
                            left_join(LandUse, by = c("site2" = "SSBS")) %>% dplyr::rename(site2LUI = LandUse_Intensity) %>%
                            dplyr::mutate(Contrast = paste(site1LUI,site2LUI, sep = "-"))
                          
                          
                          site_comparisons <- remove_dupl_comparisons(site_comparisons)
                          
                          
                          ### calculate geographic distance between the sites 
                          
                          
                          ## Collate together the variables for the dataset
                          
                          study_data <- site_comparisons %>% dplyr::select(site1,site2, Contrast)
                          study_data$SS <- study
                          study_data$convex_overlap <- NA
                          study_data$hyper_overlap <- NA
                          
                          
                          if(NROW(study_data) > 10 ){
                            
                            study_data <- sample_n(tbl = study_data,size =  10, replace = FALSE)
                          }
                          
                          ### Also want to add in a variable for the minimum number of species in either of the sites used to construct the hypervolumes. This will be used as weights in the models as there may be greater uncertainty in the hypervolume overlaps when fewer species have been recorded at either site.
                          
                          sample_plot <- list()
                          
                          for(i in 1:NROW(study_data)){
                            
   
           ##### get the species in both sites being compared 
                            
                            site1 <- study_data[i,"site1"]
                            site1_spp <- Similarity_data %>% filter(SSBS == site1) 
                            
                            
                          
                            ### site 2
                            
                            site2 <- study_data[i,"site2"]
                            site2_spp <- Similarity_data %>% filter(SSBS == site2)
                            
                            plot <- sampling_effort_check(site1 = site1_spp, site2 = site2_spp, trait_ranges = trait_ranges)
                            
                            sample_plot[[i]] <- plot
                            
                          }
                          
                          
                        plots <- sample_plot
                          
                        }
registerDoSEQ()


write_rds(file = "Functional_Intactness_Index/Effort_test.rds", Effort_Check)


plot(Effort_Check[[48]][[3]])


