######################################
# Estimate Per-Feature Tweedie Index #
######################################

##########################
# Load Required Packages #
##########################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load("cplm")

########################################
# Per-dataset Tweedie Index Estimation #
########################################

estimate.tweedie.index <- function(features){
 
 ##########################
 # Per-feature estimation #
 ##########################
 
 paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
   
   ###############################
   # Extract features one by one #
   ###############################
   
   featuresVector <- features[, x]
   data<-as.data.frame(featuresVector)
   names(data)<-'Y'
   out <- tryCatch({
      out1 <- cplm::cpglm(Y ~ 1, data = data)
      }, error = function(err){
         out1 <- try({cplm::cpglm(Y ~ 1, data = data)}) 
         return(out1)
         })
     
     #######################
     # Summarize Estimates #
     #######################
     
     if (class(out) != "try-error"){
        para<-cbind.data.frame(out$p, out$phi)
     } else{
         print(paste("Fitting problem for feature", x, "returning NA"))
         para<- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
         }
     colnames(para)<-c('index', 'dispersion')
     para$feature<-colnames(features)[x]
     return(para)
})
 
 ###################
 # Combine results #
 ###################
 
 paras<-do.call(rbind, paras)
 paras$prevalence<-apply(features, 2, function(x){mean(x>0, na.rm=T)})
 paras<-dplyr::select(paras, c('feature'), everything())
 rownames(paras)<-NULL
 return(paras)   
}


###################################################
# Tweedie Index Estimation for A List of Datasets #
###################################################

list.estimate.tweedie.index<-function(physeq){
   foreach(physeq = physeq, 
           .export = c("estimate.tweedie.index"),
           .packages = c("cplm"),
           .errorhandling = "remove") %dopar% 
      {
         start.time<-Sys.time()
         DD<-estimate.tweedie.index(physeq)
         stop.time<-Sys.time()
         time<-as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
         DD$time<-time
         return(DD)
      }
}

