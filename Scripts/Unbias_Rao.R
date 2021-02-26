Rao_Q_Func <- function(data){
  
  ### get the list of uncique species across teh whole dataset
  
  Species_abundance <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels()
  
  ### loop over every site in the dataset to collate the relative abundance of each species
  
  for(site in levels(data$SSBS)){
    
    Spp_abd <- data %>% filter(SSBS == site) %>% droplevels() %>% data.frame()
    
    
    ### Join species withining site to dataframe and rename column as site name
    
    Species_abundance <- Species_abundance %>% left_join(Spp_abd[,c("Jetz_Name", "SpeciesSiteAbundance")], by = "Jetz_Name")
    colnames(Species_abundance)[which(colnames(Species_abundance) == "SpeciesSiteAbundance")] <- paste(site)
  }
  
  ## rename rows as species and drop column from dataset 
  
  rownames(Species_abundance) <- Species_abundance$Jetz_Name
  Species_abundance <- as.matrix(Species_abundance[,-1])
  
  ## Nas to zeros
  
  Species_abundance[is.na(Species_abundance)] <- 0  
  
  ### Join all species in datasets traits scores
  
  traits <- data.frame(data) %>% dplyr::distinct(Jetz_Name) %>% droplevels() %>% left_join(PC_Scores, by = "Jetz_Name")
  rownames(traits) <- traits$Jetz_Name
  traits <- traits[,-1]
  
  
  rao.diversity
  Rao_Gow <- rao.diversity(comm = t(Species_abundance),traits =  traits)    #THIS is using the package SYNCSA that calcuates Rao's using gowdis 
  
  
  Rao_Var <- FD::dbFD(x = traits, a = t(Species_abundance), scale.RaoQ = TRUE) ## dbFD calculates Rao's Q as a measure of variance
  Rao <- data.frame(Var = as.numeric(Rao_Var$RaoQ), Gow = Rao_Gow$FunRao)
  
  
  
  return(Rao)
}

comm <- t(Species_abundance)



function (comm, traits = NULL, phylodist = NULL, checkdata = TRUE, 
          ord = "metric", put.together = NULL, standardize = TRUE, 
          ...) 
{
  diver.internal <- function(community, distance) {
   
      inter <- community %*% distance
      res <- rowSums(sweep(community, 1, inter, "*", 
                           check.margin = FALSE))
    
    return(res)
  }
  
  
  res <- list(call = match.call())



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
    
    dist.functional <- sqrt(as.matrix(FD::gowdis(x = traits)))

  }
  
  comm_unbias <- sweep(comm, 1, c(rowSums(comm, na.rm = TRUE) - 1), "/")
  comm_2 <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")
  
  
  
  
  #### community rows * trait distance columns 
  
  inter <- comm_unbias %*% distance
  res <- rowSums(sweep(comm_2, 1, inter, "*", 
                       check.margin = FALSE))
  
  
  
test <-     rao.diversity(comm = t(Species_abundance),traits =  traits)
rao2 <- test$FunRao


adiv <- adiv::QE(comm = t(Species_abundance), dis = as.dist(test))
    



    