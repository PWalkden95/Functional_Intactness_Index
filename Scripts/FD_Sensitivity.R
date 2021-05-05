rm(list = ls())

require(tidyverse)
require(magrittr)
require(gtools)
require(SYNCSA)


PREDICTS_Aves <- readRDS("../Datasets/PREDICTS/PREDICTS_Americas_Aves.rds")


PCA_Data <- data.frame(PREDICTS_Aves[,c("Jetz_Name", "Bill.TotalCulmen", "Bill.Nares", "Bill.Width", "Bill.Depth", "Tarsus.Length", "Secondary1", "Wing.Chord", "Tail.Length")])
PCA_Data <- distinct(PCA_Data, Jetz_Name, .keep_all = TRUE)


#### Log transform - then standardise and centre for a mean of zero and a standard deviation of 1

for(col in colnames(PCA_Data)[c(2:9)]){
  PCA_Data[,col] <- log(PCA_Data[,col])
  PCA_Data[,col] <- (PCA_Data[,col] - mean(PCA_Data[,col])) / sd(PCA_Data[,col])
}

rownames(PCA_Data) <- PCA_Data$Jetz_Name
PCA_Data <- PCA_Data[,-1]

### The different combinations of the full traits dropping a single trait out then recalculating the PCA


Full_Combinations <- combinations(NCOL(PCA_Data), (NCOL(PCA_Data) - 1))

i <- 1
Full.pca.list <-c() 
for(i in 1:NROW(Full_Combinations)){
  PCA <- prcomp(PCA_Data[,c(Full_Combinations[i,])], scale. = TRUE, center = TRUE)
  Full.pca.list <- c(Full.pca.list, list(PCA))
  names(Full.pca.list)[i] <- paste("Combination_", i, sep = "")
}

PCA_Data[,Full_Combinations[i,]]
## Do the same for the Two-step PCA


#### Loco combinations

Loco.pca.data <- PCA_Data[,c("Tarsus.Length","Secondary1","Wing.Chord","Tail.Length")]

Loco_Combinations <- combinations(NCOL(Loco.pca.data), (NCOL(Loco.pca.data) - 1))


Loco.pca.list <-c() 
for(i in 1:NROW(Loco_Combinations)){
  PCA <- prcomp(Loco.pca.data[,c(Loco_Combinations[i,])], scale. = TRUE, center = TRUE)
  Loco.pca.list <- c(Loco.pca.list, list(PCA))
  names(Loco.pca.list)[i] <- paste("Loco_Combination_", i, sep = "")
}


Loco.full.PCA <- prcomp(Loco.pca.data, scale. = TRUE, center = TRUE)
Loco.pca.list <- c(Loco.pca.list, list(Loco.full.PCA))
names(Loco.pca.list)[NROW(Loco_Combinations) + 1] <- "Loco_Full"


### Foraging combinations

For.pca.data <- PCA_Data[,c("Bill.TotalCulmen","Bill.Nares","Bill.Width","Bill.Depth")]

For_Combinations <- combinations(NCOL(For.pca.data), (NCOL(For.pca.data) - 1))


For.pca.list <-c() 
for(i in 1:NROW(For_Combinations)){
  PCA <- prcomp(For.pca.data[,c(For_Combinations[i,])], scale. = TRUE, center = TRUE)
  For.pca.list <- c(For.pca.list, list(PCA))
  names(For.pca.list)[i] <- paste("For_Combination_", i, sep = "")
}

For.full.PCA <- prcomp(For.pca.data, scale. = TRUE, center = TRUE)
For.pca.list <- c(For.pca.list, list(For.full.PCA))
names(For.pca.list)[NROW(For_Combinations) + 1] <- "For_Full"

#### Get the two-step permutations of the Loco and Foraging 


Two_Step_Combs <- combinations(n = length(names(c(For.pca.list,Loco.pca.list))),r = 2, v= names(c(For.pca.list,Loco.pca.list)))
Two_Step_Combs <- Two_Step_Combs[-c(which(substr(Two_Step_Combs[,1],1,4) == substr(Two_Step_Combs[,2],1,4))),]

Two_step_list <- c()
for(i in 1:NROW(Two_Step_Combs)){
  For.pca2 <- For.pca.list[[Two_Step_Combs[i,grep(Two_Step_Combs[i,], pattern = "For")]]]
  Loco.pca2 <- Loco.pca.list[[Two_Step_Combs[i,grep(Two_Step_Combs[i,], pattern = "Loco")]]]

  Body.pca.data2 <- data.frame(ForagePC1 = For.pca2$x[,1], LocoPC1 = Loco.pca2$x[,1])
  Body.pca2 <- prcomp(Body.pca.data2, scale. = TRUE, center = TRUE)
  
  two_step_PCdata <- data.frame(Jetz_Name = rownames(PCA_Data), Foraging.PCA = For.pca2$x[,2], Loco.PCA = Loco.pca2$x[,2], Body.PCA = Body.pca2$x[,1])
  
  #### Can here add the foraging trait data so that each item in the lst contains all the trait data necessary to calculate Raos Q and other functional
  ### Diversity metricss
  

  
  
  Two_step_list <- c(Two_step_list, list(two_step_PCdata))
  names(Two_step_list)[i] <- paste(Two_Step_Combs[i,1],Two_Step_Combs[i,2], sep = "*")
  }

require(coRanking)


Full_PCA <- FD::gowdis(Two_step_list[[25]])

reduced_dist <- FD::gowdis(Two_step_list[[24]])



quality <- as.matrix(coRanking::coranking(Full_PCA,reduced_dist, input_Xi = "dist", input_X = "dist"))

imageplot(quality)

AUC <- coRanking::R_NX(quality)
score <- coRanking::AUC_ln_K(AUC)


abundance_data <- readRDS("Outputs/Rao_abundance_data.rds")



source("Functions/Site_Rao_Q.R")




Rao_sensitivity <- data.frame(matrix(rep(NA,1401),nrow = 1401,ncol = 1))
for(comb in 1:length(Two_step_list)){
Rao_data <- Rao_Q_Func(abundance_data,Two_step_list[[comb]])
Rao_sensitivity <- cbind(Rao_sensitivity,Rao_data)
}

Rao_sensitivity <- Rao_sensitivity[,-1]

Rao_unbias <- Rao_sensitivity[,seq(2,50,2)]

for(i in 1:ncol(Rao_unbias)){
hist(Rao_unbias[,i])
}
t.test(x = Rao_unbias[,1],y = Rao_unbias[,2])
## save for modelling

write_rds(PREDICTS_Site_Rao, file = "Outputs/PREDICTS_Site_Rao.rds")

