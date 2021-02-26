#########################################################################
#####  Estimation of functional alpha diversity across communities ######
#########################################################################

# warnings: check the two species that are lacking for diet and foraging behaviour and the NAs of one species
# warnings: adiv has some incopatibilities with other packages; if it does not work, start again, upload the package and use only the packages subsequently requested 

### GOALS: ### 

# Estimate for each community:

# Quadratic Entropy (QE) for all morphological and behavioral traits + phylogeny
# Uniqueness* for all morphological and behavioral traits + phylogeny
# Redundancies (CR = 1-Uniqueness) for all morphological and behavioral traits + phylogeny
# Simpson index (QE taxonomy or HGS)
# Species richness
# The meanD
# The balance factor
# Correlation coefficient of the balance factor


# *A community containning species that are functionally different will achieve a high Uniqueness value, and this will increase if the relative abundance of the species is even (the absolute abundance has no effect)
# *conversely, the same community will exhibit low redundancies

## We will make the estimations for

# all species
# all native species (excluding exotics)


### INPUTS: ###

# Full species*community matrix (comm) of relative abundances
# 9 Morphological traits
# Morphological axes derived from the traits defining body size, beak shape, locmotory (tarsus) shape and wing shape
# 33 foraging behavioural traits
# 8 diet categories
# Two full phylogenies


### OUTPUT: ### 

# Data.frame containning all metrics ("Morphological diversity metrics for communities.txt")
# Data.frame containning all metrics for natives only ("Morphological diversity metrics for communities natives.txt")


### ANALYSES START HERE ### 


## There are some incompatibatilities among packages (adivd vs BAT) so we will download only the essential libraries (or unclick BAT)

library(reshape2)
library(adiv)
library(plyr)
library(ecodist)

## Community data preparation (if you run this part, you get communities*species abundances of natives "comm")

{

dat<-read.table(paste0(workingData,"/Urban global data April 25 2018 for R.txt"), header=TRUE)
dat$community <- factor(dat$community)
dat <- subset(dat,status=="native")  # if we want to exclude exotics
dat[] <- lapply(dat, function(x) if(is.factor(x)) factor(x) else x)

dat01 <- subset(dat, urban.analyses=="yes") # if we focus on communities within cities or along urbanisation gradients
dat01[] <- lapply(dat01, function(x) if(is.factor(x)) factor(x) else x)

dat02 <- dat01[dat01$relative.abundance>0,] # we exclude species absent
# dat02 <- dat01  # if we want to include all occurrence data
dat02[] <- lapply(dat02, function(x) if(is.factor(x)) factor(x) else x)

comm <- acast(na.omit(dat02[,c(5,8,17)]), community~animal, value.var="relative.abundance", fun.aggregate=mean)   # we transform it to a matrix species*community    # fun.aggregate=mean is used as otherwise it gives the length
# comm <- acast(na.omit(dat02[,c(5,8,15)]), community~animal, value.var="occurrence", fun.aggregate=mean)   # if we use presence/absence instead of abundances
comm[comm=="NaN"] <- 0



}


########### Analyses for morphological data ##############
##########################################################

{
## Functional data preparation

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)

funcdat<-func[,c(4,6,7,15, 8:14, 16)]
names(funcdat)
rownames(funcdat) <- func$animal

funcdat <- funcdat[order(rownames(funcdat)),] # we need to order the species
funcdat <- funcdat[labels(comm[1,]),]

beakshape <- as.data.frame(funcdat$beak.shape)  # dataset for beak shape
rownames(beakshape) <- rownames(funcdat)

locomshape <- as.data.frame(funcdat$locom.shape)  # dataset for locomotory system
rownames(locomshape) <- rownames(funcdat)

bodysize <- as.data.frame(funcdat$body.size)  # dataset for body size
rownames(bodysize) <- rownames(funcdat)

hand.wing <- as.data.frame(funcdat$hand_wing)  # dataset for wing shape
rownames(hand.wing) <- rownames(funcdat)

all <-funcdat[,5:12]  # dataset for all traits when analysed separately
rownames(all) <- rownames(funcdat)


## estimation of euclidean distances among morphological traits

distallmorphology <- distance(all, "euclidean")   # all 8 traits
distallmorphology <- distallmorphology/max(distallmorphology)  # we standardize to 0-1 range

distbeak <- distance(beakshape, "euclidean")    # beak shape
distbeak <- distbeak/max(distbeak)   

distlocom <- distance(locomshape, "euclidean")   # locomotory shape
distlocom <- distlocom/max(distlocom)

distsize <- distance(bodysize, "euclidean")    # body size PCA all 8 traits
distsize <- distsize/max(distsize)

distwinghand<- distance(hand.wing, "euclidean")   # wing hand index
distwinghand<- distwinghand/max(distwinghand)   


## estimation of phylogenetic distances among species
library(ape)
library(picante)

ctree1 <- read.nexus(paste0(workingData,"/AllBirdsEricson1_summary.tre"))    # This is Ericson concensus tree
ctree2 <- read.nexus(paste0(workingData,"/AllBirdsHackett1_summary.tre"))    # This is Hackett concensus tree

combined <- match.phylo.comm(ctree1, comm)
ctree.Eric <- combined$phy
comm.Eric <- combined$comm

combined <- match.phylo.comm(ctree2, comm)
ctree.Hack <- combined$phy
comm.Hack <- combined$comm

phydisE <- as.dist(cophenetic(ctree.Eric))
phydisE <- phydisE/max(phydisE)

phydisH <- as.dist(cophenetic(ctree.Hack))
phydisH <- phydisH/max(phydisH)




## We will first use adiv to estimate:

  # N (species richness, used for mistakes control)
  # Q (quadratic diversity)
  # D (Simpson diversity, i.e taxonomic Q)
  # Community uniqueness: Q/D
  # Community redundancy (CR = 1-uniqueness)

all.morph <- uniqueness(comm, distallmorphology)
beak <- uniqueness(comm, distbeak)
locom <- uniqueness(comm, distlocom)
size <- uniqueness(comm, distsize)
winghand <- uniqueness(comm, distwinghand)
phyE <- uniqueness(comm.Eric, dis = phydisE)
phyH <- uniqueness(comm.Hack, dis = phydisH)


## We will next use the function QEpart.R to estimate:

  # The meanD
  # The balance factor
  # The balance factor, described as correlation coefficient

all.morph.2 <- QEpartition(comm, distallmorphology)
beak.2 <- QEpartition(comm, distbeak)
locom.2 <- QEpartition(comm, distlocom)
size.2 <- QEpartition(comm, distsize)
winghand.2 <- QEpartition(comm, distwinghand)
phyE.2 <- QEpartition(comm.Eric, dis = phydisE)
phyH.2 <- QEpartition(comm.Hack, dis = phydisH)

## we will estimate community-weighted mean traits 

library(FD)

Wbeak <- dbFD(beakshape, comm, calc.CWM=TRUE)
Wlocomshape <- dbFD(locomshape, comm, calc.CWM=TRUE)
Wbodysize <- dbFD(bodysize, comm, calc.CWM=TRUE)
Whand.wing <- dbFD(hand.wing, comm, calc.CWM=TRUE)



## Preparing data for subsequent analyses

FDmorphology<-as.data.frame(cbind(labels(comm[,2]),all.morph$red$N,all.morph$red$D,all.morph$red$Q,all.morph$red$U,1-all.morph$red$U,beak$red$Q,beak$red$U,1-beak$red$U,locom$red$Q,locom$red$U,1-locom$red$U,size$red$Q,size$red$U,1-size$red$U,winghand$red$Q,winghand$red$U,1-winghand$red$U,phyE$red$Q,phyE$red$U,1-phyE$red$U,phyH$red$Q,phyH$red$U,1-phyH$red$U,all.morph.2$meanD,beak.2$meanD,locom.2$meanD,size.2$meanD,winghand.2$meanD,phyE.2$meanD,phyH.2$meanD,all.morph.2$Balance,beak.2$Balance,locom.2$Balance,size.2$Balance,winghand.2$Balance,phyE.2$Balance,phyH.2$Balance,Wbeak$CWM, Wlocomshape$CWM, Wbodysize$CWM, Whand.wing$CWM, all.morph.2$Balance_Cor,beak.2$Balance_Cor,locom.2$Balance_Cor,size.2$Balance_Cor,winghand.2$Balance_Cor,phyE.2$Balance_Cor,phyH.2$Balance_Cor))

colnames(FDmorphology)<-c("community","Species.richness","QE.taxonomy","QE.all.morph","Uniqueness.all.morph","CR.all.morph","QE.beak","Uniqueness.beak","CR.beak","QE.locom","Uniqueness.locom","CR.locom","QE.size","Uniqueness.size","CR.size","QE.winghand","Uniqueness.winghand","CR.winghand","QE.phyE","Uniqueness.phyE","CR.phyE","QE.phyH","Uniqueness.phyH","CR.phyH","all.morph.meanD","beak.meanD","locom.meanD","size.meanD","winghand.meanD","phyE.meanD","phyH.meanD","all.morph.Balance","beak.Balance","locom.Balance","size.Balance","winghand.Balance","phyE.Balance","phyH.Balance","CWM.beak.shape","CWM.locom.shape","CWM.body.size","CWM.hand.wing","all.morph.Balance.cor","beak.Balance.cor","locom.Balance.cor","size.Balance.cor","winghand.Balance.cor","phyE.Balance.cor","phyH.Balance.cor")


# We add habitat and study site information

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
               Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FDmorphology,tmp, by="community")
      
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities.txt"))
write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities natives.txt"))
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities ocurrences.txt"))
# write.table(tmp2, paste0(workingData,"/Morphological diversity metrics for communities ocurrences natives.txt"))


}


################ Analyses for diet  ######################
##########################################################

{
### Functional diet data

## Need to estimate first "comm" running the firts part of this code

diet <-read.table(paste0(workingData,"/Foraging Data June 13 2018 for R.txt"), header=TRUE)

names(diet)
diet <- diet[order(diet$animal),] # wee need to order the species
diet<-diet[,c(1,4:11)]
names(diet)
rownames(diet) <- diet$animal
diet <- diet[,-1]

## Estimation of distances

dietdat <- na.omit(diet[labels(comm[1,]),])  # we take only species in communities

# comm1 <- comm[,rownames(dietdat)]  # if one species is missing in diet, we will need to update community

distdiet<- dist.prop(dietdat, method = 1) # method 1 is Manly
distdiet <- distdiet/max(distdiet) 

all.diet <- QEpartition(comm, distdiet)
all.diet1 <- uniqueness(comm, distdiet)

# redundancies are estimated as 1-QE/Simpson (Ricotta et al. 2016)

redundancy <- 1-(all.diet$QE/all.diet$Simpson)



## Diet axes separately

diet.axes <-read.table(paste0(workingData,"/diet.axes.txt"), header=TRUE)

diet.axes <- na.omit(diet.axes[labels(comm[1,]),])  # we take only species in communities

comm1 <- comm[,rownames(diet.axes)]  # as one species is missing in diet, we need to update community

Inv_SeedFruit <- as.data.frame(diet.axes$Inv_SeedFruit)
rownames(Inv_SeedFruit) <- rownames(diet.axes)
DInv_SeedFruit <- distance(Inv_SeedFruit, "euclidean")   
DInv_SeedFruit <- DInv_SeedFruit/max(DInv_SeedFruit)
Inv_SeedFruit.diet <- QEpartition(comm1, DInv_SeedFruit)
Inv_SeedFruit.redundancy <- 1-(Inv_SeedFruit.diet$QE/Inv_SeedFruit.diet$Simpson)
CWMInv_SeedFruit <- dbFD(Inv_SeedFruit, comm1, calc.CWM=TRUE)


Inv_Seed <- as.data.frame(diet.axes$Inv_Seed)
rownames(Inv_Seed) <- rownames(diet.axes)
DInv_Seed <- distance(Inv_Seed, "euclidean")   
DInv_Seed <- DInv_Seed/max(DInv_Seed)
Inv_Seed.diet <- QEpartition(comm1, DInv_Seed)
Inv_Seed.redundancy <- 1-(Inv_Seed.diet$QE/Inv_Seed.diet$Simpson)
CWMInv_Seed <- dbFD(Inv_Seed, comm1, calc.CWM=TRUE)


Seed <- as.data.frame(diet.axes$Seed)
rownames(Seed) <- rownames(diet.axes)
DSeed <- distance(Seed, "euclidean")   
DSeed <- DSeed/max(DSeed)
Seed.diet <- QEpartition(comm1, DSeed)
Seed.redundancy <- 1-(Seed.diet$QE/Seed.diet$Simpson)
CWMSeed <- dbFD(Seed, comm1, calc.CWM=TRUE)



## Preparing data for subsequent analyses

FDdiet <-as.data.frame(cbind(labels(comm1[,2]), all.diet$QE, all.diet$meanD, all.diet$Balance, redundancy, Inv_SeedFruit.diet$QE, Inv_SeedFruit.diet$meanD, Inv_SeedFruit.diet$Balance, Inv_SeedFruit.redundancy, Inv_Seed.diet$QE, Inv_Seed.diet$meanD, Inv_Seed.diet$Balance, Inv_Seed.redundancy, Seed.diet$QE, Seed.diet$meanD, Seed.diet$Balance, Seed.redundancy, CWMInv_SeedFruit$CWM, CWMInv_Seed$CWM,CWMSeed$CWM, all.diet$Balance_Cor,Inv_SeedFruit.diet$Balance_Cor,Inv_Seed.diet$Balance_Cor,Seed.diet$Balance_Cor))

colnames(FDdiet)<-c("community","QE.diet","diet.meanD","diet.Balance","CR.diet","QE.PCoA1","PCoA1.meanD","PCoA1.Balance","CR.PCoA1","QE.PCoA2","PCoA2.meanD","PCoA2.Balance","CR.PCoA2","QE.PCoA3","PCoA3.meanD","PCoA3.Balance","CR.PCoA3","CWMInv_SeedFruit", "CWMInv_Seed", "CWMSeed","diet.Balance.cor","PCoA1.diet.Balance.cor","PCoA2.diet.Balance.cor","PCoA3.diet.Balance.cor")



# We add habitat and study site information

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FDdiet,tmp, by="community")

# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities.txt"))
write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities natives.txt"))
# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities ocurrences.txt"))
# write.table(tmp2, paste0(workingData,"/Diet diversity metrics for communities ocurrences natives.txt"))

}


################ Analyses for foraging  ######################
##########################################################

{
  ### Functional foraging data
  
  ## Need to estimate first "comm" running the first part of this code
  
  forag <-read.table(paste0(workingData,"/Foraging Data June 13 2018 for R.txt"), header=TRUE)
  
  names(forag)
  forag <- forag[order(forag$animal),] # wee need to order the species
  forag<-forag[,c(1,12:41)]
  names(forag)
  rownames(forag) <- forag$animal
  forag <- forag[,-1]
  
  foragdat <- na.omit(forag[labels(comm[1,]),])  # we take only species in communities
  
  distforag<- dist.prop(foragdat, method = 1) # method 1 is Manly
  distforag <- distforag/max(distforag) ## Estimation of distances
  
  all.forag <- QEpartition(comm, distforag)
  all.forag1 <- uniqueness(comm, distforag)
  
  # redundancies are estimated as 1-QE/Simpson (Ricotta et al. 2016)
  
  redundancy <- 1-(all.forag$QE/all.forag$Simpson)
  
  
  
  ## Foraging axes separately
  
  forag.axes <-read.table(paste0(workingData,"/forag.axes.txt"), header=TRUE)
  
  forag.axes <- na.omit(forag.axes [labels(comm[1,]),])  # we take only species in communities
  
  
  PCo1 <- as.data.frame(forag.axes$PCo1)
  rownames(PCo1) <- rownames(forag.axes)
  DPCo1 <- distance(PCo1, "euclidean")   
  DPCo1 <- DPCo1/max(DPCo1)
  PCo1.forag <- QEpartition(comm, DPCo1)
  PCo1.redundancy <- 1-(PCo1.forag$QE/PCo1.forag$Simpson)
  CWMPCo1 <- dbFD(PCo1, comm, calc.CWM=TRUE)
  
  
  PCo2 <- as.data.frame(forag.axes$PCo2)
  rownames(PCo2) <- rownames(forag.axes)
  DPCo2 <- distance(PCo2, "euclidean")   
  DPCo2 <- DPCo2/max(DPCo2)
  PCo2.forag <- QEpartition(comm, DPCo2)
  PCo2.redundancy <- 1-(PCo2.forag$QE/PCo2.forag$Simpson)
  CWMPCo2 <- dbFD(PCo2, comm, calc.CWM=TRUE)
  
  
  PCo3 <- as.data.frame(forag.axes$PCo3)
  rownames(PCo3) <- rownames(forag.axes)
  DPCo3 <- distance(PCo3, "euclidean")   
  DPCo3 <- DPCo3/max(DPCo3)
  PCo3.forag <- QEpartition(comm, DPCo3)
  PCo3.redundancy <- 1-(PCo3.forag$QE/PCo3.forag$Simpson)
  CWMPCo3 <- dbFD(PCo3, comm, calc.CWM=TRUE)
  
  PCo4 <- as.data.frame(forag.axes$PCo4)
  rownames(PCo4) <- rownames(forag.axes)
  DPCo4 <- distance(PCo4, "euclidean")   
  DPCo4 <- DPCo4/max(DPCo4)
  PCo4.forag <- QEpartition(comm, DPCo4)
  PCo4.redundancy <- 1-(PCo4.forag$QE/PCo4.forag$Simpson)
  CWMPCo4 <- dbFD(PCo4, comm, calc.CWM=TRUE)
  
  
  ## Preparing data for subsequent analyses
  
  FDforag <-as.data.frame(cbind(labels(comm[,2]), all.forag$QE, all.forag$meanD, all.forag$Balance, redundancy, PCo1.forag$QE, PCo1.forag$meanD, PCo1.forag$Balance, PCo1.redundancy, PCo2.forag$QE, PCo2.forag$meanD, PCo2.forag$Balance, PCo2.redundancy, PCo3.forag$QE, PCo3.forag$meanD, PCo3.forag$Balance, PCo3.redundancy, CWMPCo1$CWM, CWMPCo2$CWM,CWMPCo3$CWM, CWMPCo3$CWM, all.forag$Balance_Cor,PCo1.forag$Balance_Cor,PCo2.forag$Balance_Cor, PCo3.forag$Balance_Cor))
  
  colnames(FDforag)<-c("community","QE.forag","forag.meanD","forag.Balance","CR.forag","QE.PCoA1","PCoA1.meanD","PCoA1.Balance","CR.PCoA1","QE.PCoA2","PCoA2.meanD","PCoA2.Balance","CR.PCoA2","QE.PCoA3","PCoA3.meanD","PCoA3.Balance","CR.PCoA3","CWM_PCo1", "CWM_PCo2","CWM_PCo3","CWM_PCo4","forag.Balance.cor","PCo1.forag.Balance.cor","PCo2.forag.Balance.cor","PCo3.forag.Balance.cor")
  
  
  
  # We add habitat and study site information
  
  tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
               Regional.spp.richness = length(relative.abundance))
  
  tmp2 <- merge(FDforag,tmp, by="community")
  
  # write.table(tmp2, paste0(workingData,"/Forag diversity metrics for communities.txt"))
  write.table(tmp2, paste0(workingData,"/Forag diversity metrics for communities natives.txt"))
  # write.table(tmp2, paste0(workingData,"/Forag diversity metrics for communities ocurrences.txt"))
  # write.table(tmp2, paste0(workingData,"/Forag diversity metrics for communities ocurrences natives.txt"))
  
}



################ Species-level uniqueness and vulnerability ####################
################################################################################

{## Extracted from adiv results above

## Morphology

# Kbar is the mean functional dissimilarity of each species from the rest of species from the assemblage
all.morph$Kbar[1:10,1:10]
length(all.morph$Kbar[,1])
Kbar<- as.data.frame(t(all.morph$Kbar))
Kbar[1:10,1:10]
Kbar$community <- rownames(Kbar)
Kbar$community <- as.factor(Kbar$community)
res.Kbar <- melt(Kbar,id.bars="Kbar$community")

# Vulnerability (V) is a function of species functional relevance and its extinction risk
all.morph$V[1:10,1:10]
length(all.morph$V[,1])
V<- as.data.frame(t(all.morph$V))
V[1:10,1:10]
V$community <- rownames(V)
V$community <- as.factor(V$community)
res.V <- melt(V,id.bars="Kbar$community")

res <- na.omit(cbind(res.Kbar, res.V))
res <- res[,-c(4,5)]
colnames(res) <- c("community", "animal", "Kbar.morphology", "V.morphology")
write.table(res,paste0(workingData,"/vulnerability_species_morphology_natives.txt"))


## Diet

# Kbar is the mean functional dissimilarity of each species from the rest of species from the assemblage
all.diet1$Kbar[1:10,1:10]
length(all.diet1$Kbar[,1])
Kbar<- as.data.frame(t(all.diet1$Kbar))
Kbar[1:10,1:10]
Kbar$community <- rownames(Kbar)
Kbar$community <- as.factor(Kbar$community)
res.Kbar <- melt(Kbar,id.bars="Kbar$community")

# Vulnerability (V) is a function of species functional relevance and its extinction risk
all.diet1$V[1:10,1:10]
length(all.diet1$V[,1])
V<- as.data.frame(t(all.diet1$V))
V[1:10,1:10]
V$community <- rownames(V)
V$community <- as.factor(V$community)
res.V <- melt(Kbar,id.bars="Kbar$community")

res <- na.omit(cbind(res.Kbar, res.V))
res <- res[,-c(4,5)]
colnames(res) <- c("community", "animal", "Kbar.diet", "V.diet")
write.table(res,paste0(workingData,"/vulnerability_species_diet_natives.txt"))



## Foraging

# Kbar is the mean functional dissimilarity of each species from the rest of species from the assemblage
all.forag1$Kbar[1:10,1:10]
length(all.forag1$Kbar[,1])
Kbar<- as.data.frame(t(all.forag1$Kbar))
Kbar[1:10,1:10]
Kbar$community <- rownames(Kbar)
Kbar$community <- as.factor(Kbar$community)
res.Kbar <- melt(Kbar,id.bars="Kbar$community")

# Vulnerability (V) is a function of species functional relevance and its extinction risk
all.forag1$V[1:10,1:10]
length(all.forag1$V[,1])
V<- as.data.frame(t(all.forag1$V))
V[1:10,1:10]
V$community <- rownames(V)
V$community <- as.factor(V$community)
res.V <- melt(V,id.bars="Kbar$community")


res <- na.omit(cbind(res.Kbar, res.V))
res <- res[,-c(4,5)]
colnames(res) <- c("community", "animal", "Kbar.foraging", "V.foraging")
write.table(res,paste0(workingData,"/vulnerability_species_foraging_natives.txt"))




}


################ Analyses for morphology by diet for natives ####################
#################################################################################

{
### We will examine FD for all traits, subsetted by insectivory, nectarivory, ...

## We need to estimate first comm for natives running the beginning of the code


diet <-read.table(paste0(workingData,"/Diet urban birds 28 April 2018 for R.txt"), header=TRUE)
names(diet)
diet<-diet[,c(5,6:16)]

func<-read.table(paste0(workingData,"/morphological.axes.txt"),header=TRUE)
names(func)
all <-func[,c(1,9:16)]  # dataset for all traits when analysed separately

diet.morph <- merge(diet,all, by="animal")

diet.morph <- diet.morph[order(diet.morph$animal),] # wee need to order the species
names(diet.morph)
rownames(diet.morph) <- diet.morph$animal
diet.morph <- diet.morph[,-c(1,2)]

diet.morph <- na.omit(diet.morph[labels(comm[1,]),])  # we take only species in communities


## Insectivorous, insects >= 50% in diet

insect.50 <- diet.morph[diet.morph$Diet.Inv>40,]
comm.insect.50 <- comm[,rownames(insect.50)]
distinsect.50 <- distance(insect.50[,c(11:14,15,16,17)], "euclidean")   # all 8 traits
morph.insect.50 <- QEpart(comm.insect.50, distinsect.50)
morph.insect.50.redundancy <- 1-(morph.insect.50$QE/morph.insect.50$Simpson)


## Granivorous, seeds >= 50% in diet

seed.50 <- diet.morph[diet.morph$Diet.Seed>40,]
comm.seed.50 <- comm[,rownames(seed.50)]
distseed.50 <- distance(seed.50[,c(11:14,15,16,17)], "euclidean")   # all 8 traits
morph.seed.50 <- QEpart(comm.seed.50, distseed.50)
morph.seed.50.redundancy <- 1-(morph.seed.50$QE/morph.seed.50$Simpson)


## Frugivorous, fruits >= 50% in diet

fruit.50 <- diet.morph[diet.morph$Diet.Fruit>40,]
comm.fruit.50 <- comm[,rownames(fruit.50)]
distfruit.50 <- distance(fruit.50[,c(11:14,15,16,17)], "euclidean")   # all 8 traits
morph.fruit.50 <- QEpart(comm.fruit.50, distfruit.50)
morph.fruit.50.redundancy <- 1-(morph.fruit.50$QE/morph.fruit.50$Simpson)




## Preparing data for subsequent analyses

FD.morphol.diet <-as.data.frame(cbind(labels(comm.insect.50[,2]), morph.insect.50$QE, morph.insect.50$meanD, morph.insect.50$Balance, morph.insect.50.redundancy, morph.seed.50$QE, morph.seed.50$meanD, morph.seed.50$Balance, morph.seed.50.redundancy, morph.fruit.50$QE, morph.fruit.50$meanD, morph.fruit.50$Balance, morph.fruit.50.redundancy))

colnames(FD.morphol.diet)<-c("community","QE.insectiv","insectiv.meanD","insectiv.Balance","CR.insectiv","QE.seeds","seeds.meanD","seeds.Balance","CR.seeds","QE.fruits","fruits.meanD","fruits.Balance","CR.fruits")



# We add habitat and study site information

tmp <- ddply(dat, c("country", "city", "community", "habitat", "used.urban.nonurban"), summarise,
             Regional.spp.richness = length(relative.abundance))

tmp2 <- merge(FD.morphol.diet,tmp, by="community")

write.table(tmp2, paste0(workingData,"/Morphology-Diet diversity metrics for communities natives.txt"))

}




