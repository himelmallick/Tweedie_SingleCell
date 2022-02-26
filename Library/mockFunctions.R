################################################################################
# Given features and metadata tables generate <nreps> shuffled (mock) datasets #
################################################################################

generate_mock_datasets<-function(features, 
                                 metadata, 
                                 ID, 
                                 libSize, 
                                 nreps = 100, 
                                 subset.features = TRUE,
                                 subset.features.size = 100,
                                 subset.samples = TRUE,
                                 subset.samples.size = 100,
                                 seed = 1234,
                                 append.original.data = FALSE){
  
                                     
  ############
  # Set seed #
  ############
  
  set.seed(seed)
  
  #######################
  # Quick sanity checks #
  #######################
  
  if(nrow(features)!=nrow(metadata) || ncol(features)<subset.features.size || nrow(features)<subset.samples.size){
    stop('Dimension mismatch detected! Check your input data or change the subset size!')
  }

  #####################################
  # Placeholder for shuffled datasets #
  #####################################
  
  if(append.original.data)
  {
    simlist<-vector("list", length = nreps + 1)
    names(simlist)<-c(paste('Shuffled', 1:nreps, sep = '_'), 'Unshuffled')
    
    #######################
    # Add unshuffled data #
    #######################
    
    simlist[[nreps + 1]]$features<-features
    simlist[[nreps + 1]]$metadata<-metadata
    simlist[[nreps + 1]]$ID<-ID
    simlist[[nreps + 1]]$libSize<-libSize
    
  } else{
    
    simlist<-vector("list", length = nreps)
    names(simlist)<-paste('Shuffled', 1:nreps, sep = '_')
  }
  
  #####################
  # Iteration counter #
  #####################
  
  features_full<-features
  metadata_full<-metadata
  
  for (i in 1:nreps){
    
    ############################
    # Create Shuffled Datasets #
    ############################
    
    metadata_shuffled<-as.data.frame(metadata[sample(nrow(metadata_full), replace=FALSE),])
    colnames(metadata_shuffled)<-colnames(metadata_full)
    rownames(metadata_shuffled)<-rownames(features_full)
    colnames(metadata_shuffled)<-colnames(metadata_full)
    features_reduced<-features_full
    
    #########################################
    # Create a smaller dataset if specified #
    #########################################
    
    if(subset.features) 
      {
      subset_features<-sample(ncol(features_reduced), subset.features.size, replace = FALSE)
      features_reduced<-features_reduced[, subset_features]
      }
    
    if(subset.samples){
      subset_samples<-sample(nrow(features_reduced), subset.samples.size, replace = FALSE)
      features_reduced<-features_reduced[subset_samples,]
      metadata_shuffled<-as.data.frame(metadata_shuffled[subset_samples,])
      colnames(metadata_shuffled)<-colnames(metadata)
      rownames(metadata_shuffled)<-rownames(features_reduced)
      libSize<- features_reduced %>% rowSums()
      ID<-features_reduced %>% rownames()
      }
    
    ###################
    # Create PCL file #
    ###################
    
    simlist[[i]]$features<-features_reduced
    simlist[[i]]$metadata<-metadata_shuffled
    simlist[[i]]$ID<-ID
    simlist[[i]]$libSize<-libSize
    
  }
  
  # Return 
  return(simlist)
}
