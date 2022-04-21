###################
# Clear workspace #
###################

rm(list = ls())

#########################
# Set working directory #
#########################

workingDirectory<-"<YOUR WORKING DIRECTORY>"
setwd(workingDirectory)

##################
# Load Libraries #
##################

library(tidyverse)
library(Seurat)
library(Matrix)
library(pkgmaker)
library(foreach)
library(doParallel)
library(data.table)
library(scDatasets) # library(devtools); devtools::install_github("gongx030/scDatasets")
library(SummarizedExperiment)
library(Biobase)
library(SC2P) # library(devtools); install_github("haowulab/SC2P")
library(Biobase)
library(SingleCellExperiment)

#####################
# Load Source Files #
#####################

pkgmaker::source_files('./Library', '*.R')

######################
# Specify parameters #
######################

methods<-c('MAST',
           'edgeR',
           'DESeq2',
           'CPLM',
           'ZICP',
           'scRE')

####################
# List of datasets #
####################

datasets <- c('Klein',
              'Svensson')

######################
# Loop over datasets #
######################

for (i in 1:length(datasets)){
  
  ###################
  # Extract dataset #
  ###################
  
  dataset<-datasets[i]
  cat('Running dataset:', dataset, '\n')

  ###########################
  # Load dataset one by one #
  ###########################
  
  #########
  # Klein #
  #########

  if(dataset=='Klein'){
    
    klein_data<- ReadH5AD("./Data/klein_indrops_control_GSM1599501.h5ad")
    mat5<- klein_data@assays$RNA@counts
    features <- as.data.frame(t(as.matrix(mat5)))
    metadata<-as.data.frame(rbinom(n=nrow(features), size=1, prob=0.05))
    names(metadata)<-'CellType'
    rownames(metadata)<-rownames(features)
    
  }
 
  ############
  # Svensson #
  ############
  
  if (dataset=='Svensson'){
    
    svensson_data<- ReadH5AD("./Data/svensson_chromium_control.h5ad")
    raw_counts<- svensson_data@assays$RNA@counts
    # there are two datasets, each with 2000 cells
    colnames(raw_counts) %>% stringr::str_extract("[0-9]+_") %>% table()
    mat<- raw_counts[, grepl(pattern = "20311_", x = colnames(raw_counts))] # Dataset 1
    features <- as.data.frame(t(as.matrix(mat)))
    metadata<-as.data.frame(rbinom(n=nrow(features), size=1, prob=0.05))
    names(metadata)<-'CellType'
    rownames(metadata)<-rownames(features)
    
  }
  
  ##########################################
  # Filtering and extract other parameters #
  ##########################################
  
  features <- features[, colSums(features > 0) > nrow(features)* 0.1] # Filtering
  libSize<- features %>% rowSums()
  ID<-features %>% rownames()
  
  ##############################
  # Generate Shuffled Datasets #
  ##############################
  
  simlist<-generate_mock_datasets(features = features, 
                                  metadata = metadata, 
                                  ID = ID, 
                                  libSize = libSize)
    
  #################################### 
  # Fit Per-feature Model One By One #
  ####################################
  
  for (j in seq_along(methods)){
    
    #################
    # Print message #
    #################
    
    cat('Running method:', methods[j], '\n')

    ################################################
    # Per-feature method in a parallel environment #
    ################################################
    
    outputString<-paste(workingDirectory, 'Output/Mock_', dataset, '_Results_', methods[j], '.RData', sep='')
    
    if(!file.exists(outputString)){
      
      #############################
      # Extract normalized counts #
      #############################
      
      simlist_norm<-simlist
      if (methods[j] %in% c('CPLM', 'ZICP')){
        simlist_norm<-list.SCRANnorm(simlist) 
      }

      # Set Up Clustering Environment
      start.time<-Sys.time()
      no_cores <- detectCores() - 4
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      # Run method 
      fit.function<-eval(parse(text=paste('list', methods[j], sep ='.')))
      output<-tryCatch(eval(parse(text="fit.function(
                              simlist_norm, 
                              transformation = 'NONE',
                              multiple_qvalues = FALSE)")),
                       error = function(err){NULL})
      
      # Stop the Cluster 
      stopCluster(cl)
      
      ###############
      # Save output #
      ###############
      
      stop.time<-Sys.time()
      time<-as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
      output$time<-time
      save(output, file = outputString)
    }
  }
}
