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

output_phrase<-'negative_control'
outputString<-paste(workingDirectory, 'Output/', output_phrase, '.RData', sep='')

##################
# Load Libraries #
##################

library(tidyverse)
library(Seurat)
library(hdf5r)
library(foreach)
library(doParallel)

#####################
# Load Source Files #
#####################

pkgmaker::source_files('./Library', '*.R')

##############################################
# Create a list of negative control datasets #
##############################################

list.datasets<-vector("list", length = 5)
names(list.datasets)<-c('svensson1', 'svensson2', 'zheng', 'macosko', 'klein')

##########################################
# Assign datasets to the list one by one #
##########################################

#############################
# svensson chromium control #
#############################

svensson_data<- ReadH5AD("./Data/svensson_chromium_control.h5ad")
raw_counts<- svensson_data@assays$RNA@counts
# there are two datasets, each with 2000 cells
colnames(raw_counts) %>% stringr::str_extract("[0-9]+_") %>% table()
mat1<- raw_counts[, grepl(pattern = "20311_", x = colnames(raw_counts))]
mat2<- raw_counts[, grepl(pattern = "20312_", x = colnames(raw_counts))]
df1 <- as.data.frame(t(as.matrix(mat1)))
df2<- as.data.frame(t(as.matrix(mat2)))
list.datasets[[1]]<-df1
list.datasets[[2]]<-df2


#########################
# zheng gemcode control #
#########################

zheng_data<- ReadH5AD("./Data/zheng_gemcode_control.h5ad")
mat3<- zheng_data@assays$RNA@counts
df3 <- as.data.frame(t(as.matrix(mat3)))
list.datasets[[3]]<-df3 

###########################
# macosko dropseq control #
###########################

macosko_data<- ReadH5AD("./Data/macosko_dropseq_control.h5ad")
mat4<- macosko_data@assays$RNA@counts
df4 <- as.data.frame(t(as.matrix(mat4)))
list.datasets[[4]]<-df4

#########################
# klein indrops control #
#########################

klein_data<- ReadH5AD("./Data/klein_indrops_control_GSM1599501.h5ad")
mat5<- klein_data@assays$RNA@counts
df5 <- as.data.frame(t(as.matrix(mat5)))
list.datasets[[5]]<-df5

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
