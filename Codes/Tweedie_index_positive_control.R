###################
# Clear workspace #
###################

rm(list = ls())

#########################
# Set working directory #
#########################

workingDirectory<-"<YOUR WORKING DIRECTORY>"
setwd(workingDirectory)

################
# Setup output #
################

output_phrase<-'positive_control'
outputString<-paste(workingDirectory, 'Output/', output_phrase, '.RData', sep='')

##################
# Load Libraries #
##################

library(tidyverse)
library(Seurat)
library(foreach)
library(doParallel)
library(scDatasets) # library(devtools); devtools::install_github("gongx030/scDatasets")
library(SC2P) # library(devtools); install_github("haowulab/SC2P")
library(SummarizedExperiment)
library(Biobase)
library(SingleCellExperiment)

#####################
# Load Source Files #
#####################

pkgmaker::source_files('./Library', '*.R')

##############################################
# Create a list of positive control datasets #
##############################################

list.datasets<-vector("list", length = 11) # 7 datasets from scDatasets + 4 external datasets
names(list.datasets)<-c('blakeley', 'deng', 'guo2', 'iPSC', 'loh', 
                        'petropoulos', 'pollen', 
                        'PBMC', 
                        'Kidney',
                        'Islam', 
                        'Brain')

##########################################
# Assign datasets to the list one by one #
##########################################

for (i in 1:(length(list.datasets)-4)){
  dataset<-eval(parse(text=names(list.datasets)[i]))
  X <- assays(dataset)$count
  X <- preprocess(X, min.expressed.cell = round(0.1*ncol(X)))
  list.datasets[[i]]<-as.data.frame(t(X))
}

##################################
# Add  external dataset 1 (PBMC) #
##################################

pbmc <- readRDS(file = "./Data/pbmc3k_final.rds")
counts = pbmc@assays$RNA@counts
counts = as.matrix(counts)
features<-as.data.frame(t(counts))
list.datasets[[8]]<-features

####################################
# Add  external dataset 2 (Kidney) #
####################################

Mousesub.eset <- readRDS(file = "./Data/Mousesubeset.rds")
counts <- exprs(Mousesub.eset) # matrix of observations
features<-as.data.frame(t(counts))
list.datasets[[9]]<-features

###################################
# Add  external dataset 3 (Islam) #
###################################

set.seed(49729)
load("./Data/islam.rda")
celltype <- strsplit(colnames(islam), split = "_")
celltype <- sapply(celltype, function(x) x[1])
coldata <- data.frame(condition = celltype,
                      sample = colnames(islam),
                      row.names = colnames(islam),
                      stringsAsFactors = FALSE)
rowdata <- data.frame(gene = rownames(islam),
                      row.names = rownames(islam),
                      stringsAsFactors = FALSE)
islam <- SingleCellExperiment(assays = list(exprs = islam),
                              rowData = rowdata,
                              colData = coldata)
islam <- as(islam, 'ExpressionSet')
features<-as.data.frame(t(exprs(islam)))
list.datasets[[10]]<-features

###################################
# Add  external dataset 4 (Brain) #
###################################

data(brain_scRNAseq) 
list.datasets[[11]]<-as.data.frame(t(Y))

###################################################
# Adjust so that the order matches with the paper #
###################################################

list.datasets<-list.datasets[c(8:9, 1:7, 10:11)] # PBMC first, Kidney 2nd (6th/7th in the paper)

######################
# Calculate sparsity #
######################

summary_table<-data.frame()
for (i in seq_along(names(list.datasets))){
  summary_table[i, 'nCells']<-nrow(list.datasets[[i]])
  summary_table[i, 'nGenes']<-ncol(list.datasets[[i]])
  summary_table[i, 'Sparsity']<-round(sum(list.datasets[[i]]==0)/(nrow(list.datasets[[i]])*ncol(list.datasets[[i]]))*100, 1)
    
}
rownames(summary_table)<-names(list.datasets)

###################################
# Estimate index for all datasets #
###################################

if (!file.exists(outputString)){
  
  # Setup cluster environment
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # Run method
  output<-list.estimate.tweedie.index(list.datasets)
  
  # Stop the Cluster 
  stopCluster(cl)
  
  ##########
  # Output #
  ##########
  
  names(output)<-names(list.datasets)
  save(output, file = outputString)
  
}
