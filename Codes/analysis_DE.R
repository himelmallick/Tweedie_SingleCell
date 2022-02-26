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

datasets <- c('Kidney',
              'PBMC',
              'Petro',
              'Brain')

#############################
# Store post-filtered stats #
#############################

data_stats<-matrix(nrow = length(datasets), ncol = 4)
colnames(data_stats)<-c('nGenes', 'nCells', 'Sparsity', 'nPerCellType')
rownames(data_stats)<-datasets

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
  
  ##########
  # Kidney #
  ##########
  
  if(dataset=='Kidney'){
    
    ###########################
    # Load the Kidney dataset #
    ###########################
    
    ######################################################################################
    # Directly follow https://github.com/xuranw/MuSiC/blob/master/vignettes/vignette.Rmd #
    ######################################################################################
    
    Mousesub.eset <- readRDS(file = "./Data/Mousesubeset.rds")
    counts <- exprs(Mousesub.eset) # matrix of observations
    clusters.type = list(C1 = 'Neutro', C2 = 'Podo', C3 = c('Endo', 'CD-PC', 'LOH', 'CD-IC', 'DCT', 'PT'), C4 = c('Macro', 'Fib', 'B lymph', 'NK', 'T lymph'))
    cl.type = as.character(Mousesub.eset$cellType)
    for(cl in 1:length(clusters.type)){
      cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
    }
    pData(Mousesub.eset)$clusterType = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))
    pdata <- pData(Mousesub.eset) # data.frame of phenotypic information
    

    ###################################################
    # Subset to Epithelial (C3) and Immune (C4) cells #
    ####################################################

    keepcells<-rownames(pdata)[pdata$clusterType %in% c('C3', 'C4')]
    counts = counts[, keepcells]
    pdata = pdata[keepcells, ]
    all(colnames(counts)==rownames(pdata))
    
    ###############################
    # Force the consistent format #
    ###############################
    
    features_full<-as.data.frame(t(counts))
    features <- features_full[, colSums(features_full > 0) > nrow(features_full)* 0.2] # Following Sekula 2019
    metadata<-as.data.frame(pdata[, 'clusterType'])
    names(metadata)<-'CellType'
    metadata$CellType<-as.factor(metadata$CellType)
    metadata$CellType<-factor(metadata$CellType, levels = c('C3', 'C4'))
  
  }
  
  #########
  # Brain #
  #########

  if(dataset=='Brain'){
    
    #############
    # Load data #
    #############
    
    data(brain_scRNAseq) 

    ###############################
    # Force the consistent format #
    ###############################
    
    features<-as.data.frame(t(Y))
    metadata<-as.data.frame(design$celltype)
    names(metadata)<-'CellType'
    metadata$CellType<-factor(metadata$CellType, levels = unique(metadata$CellType))
    rownames(features)<-rownames(metadata)<-paste('cell', 1:nrow(features), sep ='')
    colnames(features)<-trimws(colnames(features)) # remove white space
    
  }
  
  ########
  # PBMC #
  ########
  
  if (dataset=='PBMC'){
    
    ######################################################################################
    # Directly follow the tutorial at https://satijalab.org/seurat/v3.0/de_vignette.html #
    ######################################################################################
    
    #########################
    # Load the PBMC dataset #
    #########################
    
    pbmc <- readRDS(file = "./Data/pbmc3k_final.rds")
    
    ####################################################################################################################
    # Directly follow https://github.com/statOmics/zinbwaveZinger/blob/master/realdata/createdata/createDataObject.Rmd #
    ####################################################################################################################
    
    keepcells = as.integer(pbmc@active.ident) %in% c(1, 2, 5)
    # keepcells = as.integer(pbmc@active.ident) %in% 1:2
    # keepgenes = pbmc@assays$RNA@var.features
    # counts = pbmc@assays$RNA@counts[keepgenes, keepcells]
    counts = pbmc@assays$RNA@counts[, keepcells]
    counts = as.matrix(counts)

    # coldata
    clusters = as.integer(pbmc@active.ident)
    # clusters = clusters[clusters %in% c(1, 2)]
    clusters = clusters[clusters %in% c(1, 2, 5)]
    cData = data.frame(seurat = clusters)
    cData$seurat<-ifelse(cData$seurat==5, 'CD8', 'CD4')
    rownames(cData) = colnames(counts)
    
    ###############################
    # Force the consistent format #
    ###############################
    
    features<-as.data.frame(t(counts))
    metadata<-cData
    names(metadata)<-'CellType'
    metadata$CellType<-as.factor(metadata$CellType)
    metadata$CellType<-factor(metadata$CellType, levels = c('CD8', 'CD4'))
  }

  #########
  # Petro #
  #########
  
  if (dataset=='Petro'){
    
    #############
    # Load data #
    #############
    
    X <- assays(petropoulos)$count
    X <- preprocess(X, min.expressed.cell = round(0.2*ncol(X))) # Following Sekula 2019
    cg <- colData(petropoulos)
    
    ###############################
    # Force the consistent format #
    ###############################
    
    features<-as.data.frame(t(X))
    metadata<-as.data.frame(cg$group)
    names(metadata)<-'CellType'
    rownames(features)<-rownames(metadata)<-paste('cell', 1:nrow(features), sep ='')
    
    #################################
    # Subset to Relevant Cell Types #
    #################################
    
    subset<-metadata %>% filter(CellType %in% c('E3', 'E4')) %>% rownames()
    features<-features[subset,]
    metadata<-as.data.frame(metadata[subset,])
    names(metadata)<-'CellType'
    metadata$CellType<-factor(metadata$CellType, levels = unique(metadata$CellType))
    all(rownames(features)==rownames(metadata))
    rownames(features)<-rownames(metadata)<-paste('cell', 1:nrow(features), sep ='')
  }
  
  ############################
  # Extract other parameters #
  ############################
  
  libSize<- features %>% rowSums()
  ID<-features %>% rownames()
  physeq<-list(features = features,
               metadata = metadata,
               libSize = libSize,
               ID = ID)
  
  
  ####################
  # Store dimensions #
  ####################
  
  data_stats[i, 'nCells']<-nrow(features)
  data_stats[i, 'nGenes']<-ncol(features)
  data_stats[i, 'Sparsity']<-round(sum(features==0)/(nrow(features)*ncol(features))*100, 1)
  data_stats[i, 'nPerCellType']<-paste(names(table(metadata)), table(metadata), collapse = ', ')
  
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
    
    outputString<-paste(workingDirectory, 'Output/', dataset, '_Results_', methods[j], '.RData', sep='')
    
    if(!file.exists(outputString)){
      
      #############################
      # Extract normalized counts #
      #############################
      
      physeq_norm<-physeq
      if (methods[j] %in% c('CPLM', 'ZICP')){
        physeq_norm<-SCRANnorm(physeq) 
      }

      # Set Up Clustering Environment
      start.time<-Sys.time()
      no_cores <- detectCores() - 4
      cl <- makeCluster(no_cores)
      registerDoParallel(cl)
      
      # Run method 
      fit.function<-eval(parse(text=paste('fit', methods[j], sep ='.')))
      output<-tryCatch(eval(parse(text="fit.function(
                              physeq_norm$features, 
                              physeq_norm$metadata, 
                              physeq_norm$libSize, 
                              physeq_norm$ID, 
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


#####################
# Post-filter stats #
#####################

View(data_stats)
