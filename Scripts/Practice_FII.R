rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(SYNCSA) ## For calculating functional diversity metric
require(vegan) ## for performing ordination (PCA/PCoA)
require(gower) ## calculating gower distances
require(raster) ## for loading in rasters of environmental variables and roads
require(geosphere) ## calculating geographic distances
require(ade4)  ### i can't remember right now
require(FD) ## calculating functional diversity metrics
require(hypervolume) ## calculating functional similarity/ diversity metrics 
require(betapart) ## for calculating functional similarity and for decomposing functional diversity into turnover and nestedness components 
require(foreach)## for parallelising loops
require(doParallel) ## ditto
require(sf)   ### for calculate density of roads
require(rgeos) ## ditto 
require(lwgeom) ### ditto

#### Load in the relevant datasets that will be collated to get all relevant information for FII calculations

Forage <- readRDS("Datasets/GBD/Trophic_Foraging_Niche.rds")
PREDICTS <- readRDS("Datasets/PREDICTS/diversity-2021-02-24-03-32-59.rds")
PREDICTS_Taxonomy <- readRDS("Datasets/PREDICTS/PREDICTS_AVES_Updated_Taxonomy.rds")
Jetz_Traits <- read.csv("Datasets/GBD/GBD_Jetz_averages_11_Nov_2020.csv")

colnames(Jetz_Traits)[1] <- "Jetz_Name" 

#### Subset the PREDICTS database to just include the Class Aves resolved to the species level and add in species level
#### Trophic niche and morphometric trait measurements 

PREDICTS_Aves <- PREDICTS %>% filter(Rank == "Species", Class == "Aves") %>% left_join(PREDICTS_Taxonomy[,c("Taxon", "Jetz_Name")], by = "Taxon") %>%
  left_join(Forage[,c(34:46)], by = c("Jetz_Name" = "Taxon")) %>%
  left_join(Jetz_Traits[,c(1,8:18)], by = "Jetz_Name")

### So as we are just going to practise the workflow for FII I will subset just those sites in the Americas


PREDICTS_Aves_Am <- PREDICTS_Aves %>% filter(UN_region == "Americas")

## Lets take a look at the sites in Americas 

table(PREDICTS_Aves_Am$Predominant_habitat, PREDICTS_Aves_Am$Use_intensity)

### so looking at this it would be good to collapse secondary vegetation into a single
### Land Use type and since there was not enough combinations of all Land_use type and 
### intensity create a new factor of LandUse Intensity 

PREDICTS_Aves_Am <- PREDICTS_Aves_Am %>% dplyr::mutate(LandUse = ifelse(grepl(pattern = "secondary", tolower(Predominant_habitat)), "Secondary Vegetation", paste(Predominant_habitat)),
                                                       LandUse = ifelse(grepl(LandUse, pattern = "Primary"), "Primary", paste(LandUse)),
                                                       
                                                       LandUse = ifelse(LandUse == "Cannot decide", NA, paste(LandUse)),
                                                       
                                                       Intensity = ifelse(Use_intensity == "Cannot decide", NA, paste(Use_intensity)),
                                                       
                                                       LandUse_Intensity = ifelse(is.na(Use_intensity), NA, paste(LandUse, Intensity, sep = "_")),
                                                       LandUse_Intensity = ifelse(grepl("NA", LandUse_Intensity), NA, paste(LandUse_Intensity)),
                                                       
                                                       LandUse_Intensity = factor(LandUse_Intensity),
                                                       
                                                       LandUse_Intensity = relevel(LandUse_Intensity, "Primary_Minimal use")) %>%
  dplyr::filter(!is.na(LandUse_Intensity)) %>%
  
  droplevels()

table( PREDICTS_Aves_Am$LandUse, PREDICTS_Aves_Am$Intensity)

#### Visulise the sites geographically

wm<-map_data("world") %>% filter(region != "Antartica" ) %>% fortify()

## site coords

site_points <- PREDICTS_Aves_Am %>% distinct(SSBS,Longitude,Latitude)

# generate and plot map

site_plot<-ggplot()+ coord_fixed()+
  geom_map(data =wm, map = wm,
           aes(group = group, map_id= region),
           fill = "darkgrey")+
  geom_point(data = fortify(site_points), aes(Longitude, Latitude),
             colour = "blue", size = 1)+
  theme_classic()

plot(site_plot)


### Correct for sampling effort at each site assuming that measurement increases linearly with sampling effort
### Rescale sampling effort so that its between 0 and 1 by dividing samplin effort at site by the max sampling effort within the study


PREDICTS_Aves_Am <- PREDICTS_Aves_Am %>% 
  
  
  dplyr::group_by(SS) %>% ## group by study 
  
  dplyr::mutate(Max_Sampling_Effort = max(Sampling_effort)) %>%  # max samplin effort in each site
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(Rescaled_Sampling_Effort = Sampling_effort/Max_Sampling_Effort) %>%  ## rescale sampling effort so that the maximum effort with a study is 1
  
  dplyr::filter(Diversity_metric_type =="Abundance" | Diversity_metric_type ==  "Occurrence") %>%
  
  dplyr::mutate(Effort_Corrected_Measurement = ifelse(Diversity_metric_is_effort_sensitive == TRUE,
                                                      Measurement/Rescaled_Sampling_Effort,
                                                      Measurement)) %>%
  
  droplevels()

## Now the Biodiversity data is ready to be used 

write_rds(PREDICTS_Aves_Am, "Functional_Intactness_Index/Outputs/PREDICTS_Americas_Aves.rds")

### To be able to calculate our metric of functional diversity at the site level we need to have 
### abundance of each species in each site (rescaled), fucntional distances between PCAs of morphometric traits
### and PCoA of species foraging and trophic niches - we have rescaled and relative abundance of species at each site within study 



###### We are going to generate major axes of variation with the morphometric traits using a two-step PCA proposed by Trisos et al 2014
#####  Most traits were positively correlated due to their positive association with body size, therefore to derive independant axes of trait
### variation I performed two PCAs. First, on Locomotory traits (Tarsus, Wing and tail) and second on Foraging traits (Beak Dimensions) - The first 
### PC of each PCA would represent an index of body size and another PCA on these scores would resolve the axes to one dimension.The two subsequent
### PCs in each of the first PCAs would represent an axis of locomotory and Foraging traits respectively. 



## Just get the columns with the morphometric trait data
PCA_Data <- data.frame(PREDICTS_Aves_Am[,c("Jetz_Name", "Bill.TotalCulmen", "Bill.Nares", "Bill.Width", "Bill.Depth", "Tarsus.Length",    "Kipp.s.Distance", "Secondary1", "Wing.Chord", "Hand.Wing.Index", "Tail.Length", "Mass")])

PCA_Data <- distinct(PCA_Data, Jetz_Name, .keep_all = TRUE)


#### Log transform - then standardise and centre for a mean of zero and a standard deviation of 1

PCA_Data <- data.frame(Jetz_Name = PCA_Data[,1], scale(log(PCA_Data[c(2:12)])))

### Perform a PCA on all the traits - This is just used to compare the utility of the Two-step PCA as opposed to just a single full PCA

Full.pca.data <- PCA_Data[,c(2:6,8:9,11)]

Full.pca <- prcomp(Full.pca.data, center = TRUE, scale. = TRUE)

summary(Full.pca)
Full.pca


### PCA on the Foraging traits - Beak Dimensions

For.pca.data <- PCA_Data[,c(2:5)]

For.pca <- prcomp(For.pca.data, center = TRUE, scale. = TRUE)
For.pca
summary(For.pca)

### PCA on the Locomotory traits - Tarsus, tail and wing dimensions

Loco.pca.data <- PCA_Data[,c(6,8,9,11)]

Loco.pca <- prcomp(Loco.pca.data, center = TRUE, scale. = TRUE)
Loco.pca
summary(Loco.pca)


#### Final PCA on the first Principal component of each of the first PCAs to derive an axis of body size 

Body.pca.data <- data.frame(LocoPC1 = Loco.pca$x[,1], ForPC1 = For.pca$x[,1])
Body.pca <- prcomp(Body.pca.data, center = TRUE, scale. = TRUE)

Body.pca
summary(Body.pca)



### Match the independent axes of trait variation to species in PREDICTS 

PC_Scores <- data.frame(Jetz_Name = PCA_Data[,1], Foraging.PCA = For.pca$x[,2], Loco.PCA = Loco.pca$x[,2], Body.PCA = Body.pca$x[,1])


### standardize the PC scores so that the maximum value is 1

for(col in colnames(PC_Scores[,-1])){
  PC_Scores[,col] <- PC_Scores[,col]/max(PC_Scores[,col])
}


write_rds(file = "Functional_Intactness_Index/Outputs/PC_Scores.rds", PC_Scores)
###### To determine variation in dietary and foraging niches I performed an Principal Coordinate analysis

PCoA_Data <- data.frame(PREDICTS_Aves_Am[,c(90:99)])
PCoA_Data <- distinct(PCoA_Data, Jetz_Name, .keep_all = TRUE)


rownames(PCoA_Data) <- PCoA_Data$Jetz_Name
PCoA_Data <- PCoA_Data[,(2:10)]

### calculate manly distances between species based on the proportion of diet obattain from different foraging or trophic guild. 

Manly <- dist.prop(PCoA_Data,1)

## the function cmdscale performs the principal coordinate analysis. 

ForPCoA <- cmdscale(Manly, eig = TRUE, x.ret = TRUE, k = 8)

### extract the proportion of variance explain by each PCoA axis and then the point of each species no the first two axes. 

PCoAVar <- round(ForPCoA$eig/sum(ForPCoA$eig)*100,1)
mds.values <- ForPCoA$points
mds.data <- data.frame(X = mds.values[,1], Y = mds.values[,2])

## Visualise
 
ggplot(data = mds.data, aes(x= X, y = Y)) +
  geom_point() +
  theme_bw() +
  xlab(paste("PCoA1 - ", PCoAVar[1], "%", sep = "")) +
  ylab(paste("PCoA2 - ", PCoAVar[2], "%", sep = ""))


PCoA_Scores <- data.frame(Jetz_Name = rownames(PCoA_Data), PCoA1 = mds.values[,1], PCoA2 = mds.values[,2], PCoA3 = mds.values[,3], PCoA4 = mds.values[,4])

Manly_t <- dist.prop(data.frame(t(PCoA_Data)),1)

Axes_PCo <- cmdscale(Manly_t, eig = TRUE, x.ret = TRUE, k = 8)

Axes_PCo$points



#########################################
############# Trait Data ################
#########################################

trait_data <- data.frame(PREDICTS_Aves_Am) %>% dplyr::select(Jetz_Name) %>%
  left_join(PC_Scores, by = "Jetz_Name") %>% 
  left_join(PCoA_Scores, by = "Jetz_Name") %>%
  distinct(Jetz_Name, .keep_all = TRUE)

### standardise the PCA and PcoA scores

for(col in colnames(trait_data[,-1])){
  trait_data[,col] <- trait_data[,col]/max(trait_data[,col])
}


#######################################
#### Functional Diversity of sites ####
#######################################

### RAO's Qaudratic Entropy an abundance-weighted measure of diversity - in our case functional diversity 


## get abundance data for each species with each site  
  
abundance_data <- PREDICTS_Aves_Am %>% dplyr::filter(Diversity_metric == "abundance") %>% dplyr::filter(Effort_Corrected_Measurement != 0) %>%
  
  ### filter out studies of just a single species 
  
  dplyr::group_by(SS) %>% dplyr::mutate(study_n_species = n_distinct(Jetz_Name)) %>% dplyr::filter(study_n_species > 1) %>% ungroup() %>%  
  
  #### group by Site and Species get abundance if there are some sites that the same species is recorded multiple times 
  dplyr::group_by(SSBS,Jetz_Name) %>% dplyr::mutate(SpeciesSiteAbundance = sum(Effort_Corrected_Measurement), n_spp = n()) %>%
  
  filter(!duplicated(n_spp) | n_spp == 1) %>%
  
  ungroup() %>%   droplevels() %>%  
  
  ### group by just site to get Total site abundance
  group_by(SSBS) %>% dplyr::mutate(TotalSiteAbundance = sum(SpeciesSiteAbundance), site_species = n_distinct(Jetz_Name)) %>% 
  filter(site_species > 1 ) %>%
  
  ungroup() %>%
  
  ### relative abundance of each species at each site SpeciesSitelevel abundance/TotalSite abundance
  dplyr::mutate(RelativeAbundance = SpeciesSiteAbundance/TotalSiteAbundance) %>%
  
  droplevels()

write_rds(abundance_data, file = "Functional_Intactness_Index/Outputs/abundance_data.rds")


morpho_traits <- readRDS("Functional_Intactness_Index/Outputs/PC_Scores.rds")

#####################################################
#### Function to calculate Rao's Q for each site ####
####################################################


Rao_Q_Func <- function(data, traits){
  
  ### get the list of uncique species across teh whole dataset
  
  Species_abundance <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels()
  
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  for(site in levels(data$SSBS)){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance <- Species_abundance %>% left_join(Spp_abd[,c("Jetz_Name", "RelativeAbundance")], by = "Jetz_Name")
    colnames(Species_abundance)[which(colnames(Species_abundance) == "RelativeAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance) <- Species_abundance$Jetz_Name
  Species_abundance <- as.matrix(Species_abundance[,-1])
  
  ## Nas to zeros
  
  Species_abundance[is.na(Species_abundance)] <- 0  
  
  ### Join all species in datasets traits scores
  
  spp_traits <- traits %>% filter(Jetz_Name %in% rownames(Species_abundance))
  rownames(spp_traits) <- spp_traits$Jetz_Name
  spp_traits <- spp_traits[,-1]
  

  
  Rao_Bias <- rao.diversity(comm = t(Species_abundance),traits =  spp_traits)    #THIS is using the package SYNCSA that calcuates Rao's using gowdis 
  
  ######################################################################
  ########## Here we are also going to calculate an "unbiased" Raos Q###
  ######################################################################
  
  
  Species_abundance_2 <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels()
  
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  for(site in levels(data$SSBS)){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance_2 <- Species_abundance_2 %>% left_join(Spp_abd[,c("Jetz_Name", "SpeciesSiteAbundance")], by = "Jetz_Name")
    colnames(Species_abundance_2)[which(colnames(Species_abundance_2) == "SpeciesSiteAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance_2) <- Species_abundance_2$Jetz_Name
  Species_abundance_2 <- as.matrix(Species_abundance_2[,-1])
  
  ## Nas to zeros
  
  Species_abundance_2[is.na(Species_abundance_2)] <- 0  
  
  
  comm <- t(Species_abundance_2)
 
   Rao_Unbias <- rao_diversity_2(comm = comm, traits = spp_traits)
   Rao_Unbias <- Rao_Unbias$FunRao
  
  
  Rao <- data.frame(Bias = as.numeric(Rao_Bias$FunRao),
                    Unbias = Rao_Unbias)
  
  
  
  return(Rao)
}


PREDICTS_Site_Rao <- data.frame(abundance_data) %>% distinct(SSBS, .keep_all = TRUE)%>% dplyr::select(SS, SSB, SSBS,UN_subregion, LandUse, Use_intensity, LandUse_Intensity, Longitude,Latitude, site_species)

PREDICTS_Site_Rao$SSBS <- as.character(PREDICTS_Site_Rao$SSBS)
PREDICTS_Site_Rao$Bias_Rao <- NA
PREDICTS_Site_Rao$Unbias_Rao <- NA

Rao_data <- Rao_Q_Func(abundance_data,morpho_traits)

PREDICTS_Site_Rao$Bias_Rao <- Rao_data$Bias
PREDICTS_Site_Rao$Unbias_Rao <- Rao_data$Unbias

plot(Rao_data$Unbias ~ Rao_data$Bias)
abline(a=0,b=1)
## save for modelling

write_rds(PREDICTS_Site_Rao, file = "Functional_Intactness_Index/Outputs/PREDICTS_Site_Rao.rds")



########### Might want to try and look for outliers where and how we are going to treat them, it looks like there are a few sites that are dominated 
######### with thousands of the same species which is dragging down the measue of Rao's Q - This would also be a good point to possibly caluclate some
####### other indices of FD e.g FRic, FDis etc etc 

#### Can also calculate a few other metric for FD and will do a little later -- FD, Hvol etc




#############################################
#### Functional Similarity between sites ####
#############################################

### Essentially we want to know how much functional overlap there is between Primary minimal sites and sites of other land-use types and intensities.
### This will be calculated by working out the functional beta diversity between sites and how much of this beta-diversity similarity is due to 
### nestedness of functional diversity - as opposed to turnover. This is done by first estimating the multidimensional functional space for the morpho
### -logical and foraging traits for species at each site. -- Actually to get a measure of similarity I think it's just best to get a measure of overlap
### looking at the nestedness is just say what proportion of the dissimilarity is due to the smaller site containing a subset of function of the larger
### site and nestedness is just the replacement of function. However, that can be a post analysis as we are mainly concerned with the similarity for FII


### Join trait values filter out species that are marked as absent from sites 


Similarity_data <- data.frame(PREDICTS_Aves_Am) %>% left_join(PC_Scores, by = "Jetz_Name") %>% filter(Effort_Corrected_Measurement > 0) %>%
  
  ### calculating a hypervolume in 3 dimensions with fewer than 21 species on result in inaccurcies 
  group_by(SSBS) %>% dplyr::mutate(Site_spp = n_distinct(Jetz_Name)) %>% filter(Site_spp > 21) %>% ungroup() %>%
  
  ## how many studies have at least one site of primary minimal for comparisons to be made and make sure the sites have more than a single site. 
  
  group_by(SS) %>% dplyr::mutate(n_primary_min = sum(LandUse_Intensity == "Primary_Minimal use" ), n_sites_within_studies = n_distinct(SSBS)) %>%
  
  
  ungroup() %>% filter(n_primary_min > 0 ) %>%  filter(n_sites_within_studies > 1) %>%
  
  droplevels() %>%
  
  data.frame()

write_rds(Similarity_data, file = "Functional_Intactness_Index/Outputs/similarity_data.rds")




########################################################
###### Loop to generate functional overlap of sites ####
########################################################

### With studies with multiple pri-min sites it compares with each other twice -- prim1 v prim2 & prim2 v prim 1 etc etc should only really compare once therefore i have create a function that will remove the duplicated comparisons.

remove_dupl_comparisons <- function(data){
  
  data$drop <- NA
  data$comparison <- paste(data$site1,data$site2, sep = "_")
  data[1,"drop"] <- FALSE
  
  
  if(NROW(data) > 1){
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


### What studies do we have in the dataset 

studies <- levels(Similarity_data$SS) 

#### empty data frame for results to go into 





registerDoParallel(cores = 4)

Overlap_data <- foreach(study = studies,
                        .combine = "rbind",
                        .packages = c("hypervolume","betapart","tidyverse","magrittr","geosphere","gower","raster")) %dopar% {
                          
                          
                          ### data for study 
                          Overlap_data <- c()
                          
                          
                          
                          data <- Similarity_data %>% filter(SS == study) %>% droplevels()
                          
                          ### coordinates going to be used to calculate distance between sites
                          
                          LatLong <- data.frame(data) %>% distinct(SSBS,Latitude,Longitude)
                          
                          ### land use and intensity at each site
                          
                          LandUse <- data.frame(data) %>% distinct(SSBS, LandUse_Intensity)
                          
                          ### which sites are of primary minimal landuse type and intensity  
                          
                          Primary_sites <- LandUse %>% filter(LandUse_Intensity == "Primary_Minimal use") %>% pull(SSBS) %>% as.character()
                          
                          #### all other sites 
                          
                          all_sites <- LandUse %>% pull(SSBS) %>% as.character()
                          
                          
                          
                          
                          ###### get the comparisons with the primary minimal site and all other sites within the study 
                          
                          site_comparisons <- expand.grid(Primary_sites, all_sites) %>% 
                            dplyr::rename(site1 = Var1, site2 = Var2) %>% dplyr::mutate(site1 = as.character(site1), site2 = as.character(site2)) %>%
                            filter(site1 != site2) %>%
                            left_join(LatLong, by = c("site1" = "SSBS")) %>% dplyr::rename(site1Lat = Latitude, site1Long = Longitude) %>%
                            left_join(LatLong, by = c("site2" = "SSBS")) %>% dplyr::rename(site2Lat = Latitude, site2Long = Longitude) %>%
                            left_join(LandUse, by = c("site1" = "SSBS")) %>% dplyr::rename(site1LUI = LandUse_Intensity) %>%
                            left_join(LandUse, by = c("site2" = "SSBS")) %>% dplyr::rename(site2LUI = LandUse_Intensity) %>%
                            dplyr::left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site1" = "SSBS")) %>%
                            dplyr::rename(site1_spp = Site_spp) %>%
                            left_join(distinct(Similarity_data[,c("SSBS","Site_spp")]), by = c("site2" = "SSBS")) %>%
                            dplyr::rename(site2_spp = Site_spp) %>%
                            dplyr::mutate(Contrast = paste(site1LUI,site2LUI, sep = "-"))
                          
                          
                          site_comparisons <- remove_dupl_comparisons(site_comparisons)
                          
                          
                          ### calculate geographic distance between the sites 
                          
                          site_one <- as.matrix(site_comparisons[,c("site1Long", "site1Lat")])
                          site_two <- as.matrix(site_comparisons[,c("site2Long", "site2Lat")])
                          
                          dist <- data.frame(distHaversine(site_one,site_two))
                          
                          ## Collate together the variables for the dataset
                          
                          study_data <- site_comparisons %>% dplyr::select(site1,site2, Contrast, site1Long, site1Lat, site2Long,site2Lat, site1_spp, site2_spp)
                          study_data$min_site_spp <- ifelse(study_data$site1_spp > study_data$site2_spp, study_data$site2_spp, study_data$site1_spp) 
                          study_data <- cbind(study_data, distance = dist[,1])
                          study_data$SS <- study
                          #study_data$convex_overlap <- NA
                          study_data$hyper_overlap <- NA
                          
                          
                          ### Also want to add in a variable for the minimum number of species in either of the sites used to construct the hypervolumes. This will be used as weights in the models as there may be greater uncertainty in the hypervolume overlaps when fewer species have been recorded at either site.
                          
                          
                          for(i in 1:NROW(site_comparisons)){
                            
                            
                            ##### get the species in both sites being compared 
                            
                            site1 <- site_comparisons[i,"site1"]
                            site1_spp <- Similarity_data %>% filter(SSBS == site1) %>% distinct(Jetz_Name, .keep_all = FALSE)
                            
                            
                            ### join the traits and drop species names
                            
                            site1_data <- site1_spp %>% left_join(PC_Scores, by = "Jetz_Name") 
                            rownames(site1_data) <- site1_data$Jetz_Name
                            site1_data <- as.matrix(site1_data[,-1])
                            
                            
                            ### calculate support vector machine and minimum convex hull hypervolumes 
                            
                            hypersvm_1 <- hypervolume(site1_data, method = "svm")
                            convex_1 <- expectation_convex(site1_data, check.memory = FALSE)
                            
                            
                            ### site 2
                            
                            site2 <- site_comparisons[i,"site2"]
                            site2_spp <- Similarity_data %>% filter(SSBS == site2) %>% distinct(Jetz_Name, .keep_all = FALSE)
                            
                            site2_data <- site2_spp %>% left_join(PC_Scores, by = "Jetz_Name") 
                            rownames(site2_data) <- site2_data$Jetz_Name
                            site2_data <- as.matrix(site2_data[,-1])
                            
                            hypersvm_2 <- hypervolume(site2_data, method = "svm")
                            #convex_2 <- expectation_convex(site2_data, check.memory = FALSE)
                            
                            svm_list <- hypervolume_set(hypersvm_1,hypersvm_2, check.memory = FALSE)
                            #convex_list <- hypervolume_set(convex_1, convex_2, check.memory = FALSE)
                            
                            svm_overlap <- hypervolume_overlap_statistics(svm_list)
                            #convex_overlap <- hypervolume_overlap_statistics(convex_list)
                            
                            study_data[i,"hyper_overlap"] <- svm_overlap[1]
                            #study_data[i,"convex_overlap"] <- convex_overlap[1]
                            
                          }
                          
                          Overlap_data <- rbind(Overlap_data,study_data)
                          
                          
                          
                        }



write_rds(file = "Functional_Intactness_Index/Outputs/Functional_Overlap_data.rds", Overlap_data)


registerDoSEQ()




###########################
###Human Population Density
###########################

PREDICTS_Site_Rao <- readRDS("Functional_Intactness_Index/Outputs/PREDICTS_Site_Rao.rds")
Overlap_data <- readRDS("Functional_Intactness_Index/Outputs/Functional_Overlap_data.rds")



hpd <- raster("Datasets/PREDICTS_variables/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")


### calculate human population density for Functional diversity sites 

hpd_values <- raster::extract(hpd,PREDICTS_Site_Rao[,c("Longitude","Latitude")])

## human population density trandformed with the log +1 transformation
PREDICTS_Site_Rao$logHPD <- log(hpd_values + 1) 



### calculate for functional similarity/overlap data

hpd_values_site1 <- raster::extract(hpd,Overlap_data[,c("site1Long","site1Lat")])
hpd_values_site2 <- raster::extract(hpd,Overlap_data[,c("site2Long","site2Lat")])

Overlap_data$S1logHPD <- log(hpd_values_site1 + 1)
Overlap_data$S2logHPD <- log(hpd_values_site2 + 1)

Overlap_data$logHPDdiff <- Overlap_data$S2logHPD - Overlap_data$S1logHPD


########################
#Environmental distance
########################


Bioclim_5 <- raster("Datasets/Environmental_Variables/wc2.1_30s_bio_5.tif")
Bioclim_6 <- raster("Datasets/Environmental_Variables/wc2.1_30s_bio_6.tif")
Bioclim_13 <- raster("Datasets/Environmental_Variables/wc2.1_30s_bio_13.tif")
Bioclim_14 <- raster("Datasets/Environmental_Variables/wc2.1_30s_bio_14.tif")
Elevation <- raster("Datasets/Environmental_Variables/wc2.1_30s_elev.tif")


### we only have to calculate the environmental distance when comparing two sites.

site_one <- as.matrix(Overlap_data[,c("site1Long","site1Lat")])
site_two <- as.matrix(Overlap_data[,c("site2Long","site2Lat")])

environ_1 <- data.frame(Ele  = raster::extract(Elevation,site_one), B5 = raster::extract(Bioclim_5,site_one), B6 = raster::extract(Bioclim_6, site_one),
                        B13 = raster::extract(Bioclim_13, site_one), B14 = raster::extract(Bioclim_14, site_one))

environ_2 <- data.frame(Ele  = raster::extract(Elevation,site_two), B5 = raster::extract(Bioclim_5,site_two), B6 = raster::extract(Bioclim_6, site_two),
                        B13 = raster::extract(Bioclim_13, site_two), B14 = raster::extract(Bioclim_14, site_two))


## calculate the gowers distance between environmental variables 

environ_dist <- gower_dist(x = environ_1, y = environ_2)
#### I dont know how to best transform it now so Im just going to keep it as is

Overlap_data$env_distance <- environ_dist




########################
## Density of Roads
########################

### load in the roads shape file for the americas
Roads <- st_read("Datasets/PREDICTS_variables/groads-v1-americas-shp/groads-v1-americas-shp/gROADS-v1-americas.shp")
### combine into a single shapefile
Roads <- st_combine(Roads)
### transform to be projected on the mercator projection that deals in meters rather than latlong 
Roads <- st_transform(Roads, crs = "+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Both datasets are going to need information on road densities and will be some overlap of sites so make a vector of all unique sites across both datasets

sites <- unique(c(Overlap_data$site1,Overlap_data$site2, PREDICTS_Site_Rao$SSBS))

road_densities <- c()

registerDoParallel(cores = 4)


road_densities <- foreach(site = sites,
                          .combine = "rbind",
                          .packages = c("tidyverse", "sf", "rgeos", "lwgeom")) %dopar%{
                            
                            
                            site_LongLat <- PREDICTS_Aves_Am %>% filter(SSBS %in% site) %>% distinct(SSBS,Longitude,Latitude)
                            
                            
                            point <- as.matrix(site_LongLat[,c("Longitude","Latitude")])
                            
                            point <- st_point(point)
                            point <- st_sfc(point)
                            
                            ### first crs needs to be in longlat format as the points are coordinates then transformed into a meters based projection such as mercator to calculate distance in km
                            
                            st_crs(point) <- "+proj=longlat +datum=WGS84 +no_defs"
                            point <- st_transform(point, crs = "+proj=merc +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
                            
                            ## 1k and 50k radius around points
                            
                            buffer_1km <- st_buffer(point, 1000)
                            buffer_50km <- st_buffer(point, 50000)
                            
                            
                            ## overlap buffer with roads shape file to just get roads within each radius
                            
                            intersect_1km <- st_intersection(Roads, buffer_1km)
                            
                            #### if there are any roads caluclate the density as length of roads/area of radius
                            
                            if(length(intersect_1km) != 0 ){
                              density_1km <- (st_length(intersect_1km)/1000)/(st_area(buffer_1km)/1000000)
                            } else {
                              density_1km <- 0
                            }
                            
                            intersect_50km <- st_intersection(Roads, buffer_50km)
                            
                            if(length(intersect_50km) != 0) {
                              density_50km <- (st_length(intersect_50km)/1000)/(st_area(buffer_50km)/1000000)
                            } else {
                              density_50km <- 0
                            }
                            
                            
                            densities <- data.frame(site = paste(site), density_1km = as.numeric(density_1km), density_50km = as.numeric(density_50km))
                            
                            
                          }

registerDoSEQ()

## simply join road densities to functional diversity dataset

PREDICTS_Site_Rao <- PREDICTS_Site_Rao %>% dplyr::left_join(road_densities, by = c("SSBS" = "site"))

### join to similarity dataset and calculate difference in density of roads

Overlap_data <- Overlap_data %>% 
  dplyr::left_join(road_densities, by = c("site1" = "site")) %>%
  dplyr::rename(S1RD1K = density_1km, S1RD50K = density_50km) %>%
  dplyr::left_join(road_densities, by = c("site2" = "site")) %>%
  dplyr::rename(S2RD1K = density_1km, S2RD50K = density_50km)


write_rds(file = "Datasets/PREDICTS_variables/Road_densities1_and_50k.rds", road_densities)




write_rds(file = "Functional_Intactness_Index/Outputs/Functional_Overlap_data.rds", Overlap_data)
write_rds(file = "Functional_Intactness_Index/Outputs/PREDICTS_Site_Rao.rds", PREDICTS_Site_Rao)




rao_diversity_2 <- function (comm, traits = NULL, phylodist = NULL, checkdata = TRUE, 
          ord = "metric", put.together = NULL, standardize = TRUE, 
          ...) 
{
  diver.internal <- function(community1,community2, distance) {
    if (any(is.na(distance))) {
      distance.na <- ifelse(is.na(distance), 0, 1)
      inter.na <- community %*% distance.na
      adjustment <- rowSums(sweep(community, 1, inter.na, 
                                  "*", check.margin = FALSE))
      distance[is.na(distance)] <- 0
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", check.margin = FALSE))
      res <- ifelse(adjustment > 0, res/adjustment, res)
    }
    else {
      inter <- community1 %*% distance
      res <- rowSums(sweep(community2, 1, inter, "*", check.margin = FALSE))
    }
    return(res)
  }
  res <- list(call = match.call())
  if (inherits(comm, "metacommunity.data")) {
    if (!is.null(traits) | !is.null(phylodist) | !is.null(put.together)) {
      stop("\n When you use an object of class metacommunity.data the arguments traits, phylodist and put.together must be null. \n")
    }
    traits <- comm$traits
    phylodist <- comm$phylodist
    put.together <- comm$put.together
    comm <- comm$community
  }
  list.warning <- list()
  if (checkdata) {
    organize.temp <- organize.syncsa(comm, traits = traits, 
                                     phylodist = phylodist, check.comm = TRUE)
    if (!is.null(organize.temp$stop)) {
      organize.temp$call <- match.call()
      return(organize.temp)
    }
    list.warning <- organize.temp$list.warning
    comm <- organize.temp$community
    traits <- organize.temp$traits
    phylodist <- organize.temp$phylodist
  }
  if (length(list.warning) > 0) {
    res$list.warning <- list.warning
  }
  if (any(is.na(comm))) {
    stop("\n community data with NA\n")
  }
  comm <- as.matrix(comm)
  N <- nrow(comm)
  S <- ncol(comm)
  dist.1 <- 1 - diag(x = rep(1, S))
  if (!is.null(traits)) {
    traits <- as.data.frame(traits)
    m <- ncol(traits)
    weights <- rep(1, m)
    make.names <- is.null(colnames(traits))
    colnames(traits) <- colnames(traits, do.NULL = FALSE, 
                                 prefix = "T")
    names(weights) <- colnames(traits)
    if (!is.null(put.together)) {
      if (!inherits(put.together, "list")) {
        stop("\n put.together must be a object of class list\n")
      }
      if (make.names) {
        for (k in 1:length(put.together)) {
          put.together[[k]] <- paste("T", put.together[[k]], 
                                     sep = "")
        }
      }
      if (max(table(unlist(put.together))) > 1) {
        stop("\n The same trait appears more than once in put.together\n")
      }
      if (length(setdiff(unlist(put.together), colnames(traits))) > 
          0) {
        stop("\n Check traits names in put.together\n")
      }
      for (k in 1:length(put.together)) {
        weights[put.together[[k]]] <- 1/length(put.together[[k]])
      }
    }
    dist.functional <- sqrt(as.matrix(FD::gowdis(x = traits, 
                                                 asym.bin = NULL, ord = ord, w = weights, ...)))
    if (checkdata) {
      if (any(is.na(dist.functional))) {
        warning("Warning: NA in distance between species", 
                call. = FALSE)
      }
    }
  }
  if (!is.null(phylodist)) {
    dist.phylogenetic <- as.matrix(phylodist)
    if (checkdata) {
      if (any(is.na(dist.phylogenetic))) {
        warning("Warning: NA in phylodist", call. = FALSE)
      }
    }
    if (standardize) {
      dist.phylogenetic <- dist.phylogenetic/max(dist.phylogenetic, 
                                                 na.rm = TRUE)
    }
  }
  comm1 <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
  comm2 <- sweep(comm, 1, c(rowSums(comm, na.rm = TRUE)-1), "/")
  SD <- diver.internal(community2 = comm1,community1 = comm1, distance =  dist.1)
  res$Simpson <- SD
  if (!is.null(traits)) {
    FD <- diver.internal(community1 = comm1,community2 = comm2, distance = dist.functional)
    res$FunRao <- FD
    res$FunRedundancy <- SD - FD
  }
  if (!is.null(phylodist)) {
    PD <- diver.internal(community2 = comm1,community1 = comm1, distance = dist.phylogenetic)
    res$PhyRao <- PD
    res$PhyRedundancy <- SD - PD
  }
  return(res)
}

