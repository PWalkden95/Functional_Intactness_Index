create_vcv <- function(data, level, tree) {
  
  data <- data.frame(droplevels(data))
  species <- sub(unique(data$Jetz_Name),pattern = " " ,replacement = "_")
  drop.species <- tree$tip.label[which(!(tree$tip.label %in% species))]
  
  
  full_tree <- drop.tip(tree, drop.species)
  
  
  ID <- unique(PREDICTS_site[, level]) %>% droplevels() %>% data.frame()
  ID <- unique(as.character(ID[,level]))
  
  
  comm_data <- t(species)
  colnames(comm_data) <- species
  comm_data <- data.frame(comm_data[-1,])
  
  
  
  for(i in 1:length(ID)){
    
    ID_data <- data.frame(data[data[,level] == ID[i],c("Jetz_Name", "Effort_Corrected_Measurement")])
    
    for(spp in species){
      comm_data[i,paste(spp)] <- ifelse(any(ID_data[ID_data$Jetz_Name == sub(spp,pattern = "_", replacement = " "),"Effort_Corrected_Measurement"] > 0),1,0)
    }
    
    
    rownames(comm_data)[i] <- ID[i]
  }
  
  
  for(i in 1:ncol(comm_data)){
    comm_data[,i] <- as.numeric(comm_data[,i])
  }
  
  comm_data <- as.matrix(comm_data)
  
  suppressWarnings(memory.limit(120000))
  
  vcv <- 1 - as.matrix(unifrac(comm = comm_data, tree = full_tree))
  
  
  
  
  return(vcv)  
}
