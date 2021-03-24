rm(list = ls())

#dir.create("Functional_Intactness_Index")

require(tidyverse) ## For data manipulation
require(magrittr) ## for piping
require(vegan) ## for performing ordination (PCA/PCoA)
require(ade4)


### To be able to calculate our metric of functional diversity at the site level we need to have 
### abundance of each species in each site (rescaled), fucntional distances between PCAs of morphometric traits
### and PCoA of species foraging and trophic niches - we have rescaled and relative abundance of species at each site within study 

###### We are going to generate major axes of variation with the morphometric traits using a two-step PCA proposed by Trisos et al 2014
#####  Most traits were positively correlated due to their positive association with body size, therefore to derive independant axes of trait
### variation I performed two PCAs. First, on Locomotory traits (Tarsus, Wing and tail) and second on Foraging traits (Beak Dimensions) - The first 
### PC of each PCA would represent an index of body size and another PCA on these scores would resolve the axes to one dimension.The two subsequent
### PCs in each of the first PCAs would represent an axis of locomotory and Foraging traits respectively. 


PREDICTS_Aves_Am <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")


## Just get the columns with the morphometric trait data
PCA_Data <- data.frame(PREDICTS_Aves_Am[,c("Jetz_Name", "Bill.TotalCulmen", "Bill.Nares", "Bill.Width", "Bill.Depth", "Tarsus.Length",    "Kipp.s.Distance", "Secondary1", "Wing.Chord", "Hand.Wing.Index", "Tail.Length", "Mass")])

PCA_Data <- distinct(PCA_Data, Jetz_Name, .keep_all = TRUE)

rownames(PCA_Data) <- PCA_Data$Jetz_Name
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






###########################################
###### Foraging Traits ####################
###########################################



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



#########################################
############# Trait Data ################
#########################################


### standardise the PCA and PcoA scores

for(col in colnames(PCoA_Scores[,-1])){
  PCoA_Scores[,col] <- PCoA_Scores[,col]/max(PCoA_Scores[,col])
}


trait_list <- list(morpho_traits = PC_Scores, forag_traits = PCoA_Scores)


write_rds(trait_list, "Outputs/MeanTraits.rds")
