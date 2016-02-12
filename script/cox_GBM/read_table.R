# Yuan Yuan 2012-12-10
# This script read the tab-separated file and handles proper rownames/colnames assignment, also deals with the end-line tab
# also remove the non-informative features
library(randomForest)

synapse.read <- function(synapse.ID) {
  synapse.data <- downloadEntity(synapse.ID)
  synapse.data <- loadEntity(synapse.data)
  fileName <- file.path(synapse.data$cacheDir, synapse.data$files[1])
  return(fileName)
}

myRead <- function(fileName) {
  start.time <- proc.time()
  data <- read.table(fileName, header=T,na.strings=c("[Pending]","[Not Available]","[Not Applicable]","null","null ","NA"),  sep="\t", quote="")
  samples <- as.character(data[,1])
  data <- data[,2:(ncol(data)-1)]
  rownames(data) <- samples
  print(paste(fileName, ": ", nrow(data), "(samples),", ncol(data), "(features)"))
  data <- na.roughfix(data) # impute the missing value in x, require randomForest library
  valid.cols <- which(apply(data, 2, sd)> 0) # remove the non-informative ones, with same values (e.g., 0) across all samples
  print(paste((ncol(data)-length(valid.cols)), "invalid cols were removed."))
  data <- data[,valid.cols]    
  runtime = proc.time() - start.time
  print(paste("========= Read", ncol(data), 'cols from', fileName, 'in',
    runtime[1], '============'))
  return (data)
}

# deals with clinical data, in which factors and numbers may be mixed but not check for flat values
myRead.simple <- function(fileName) {
  #fileName <- synapse.read(synapse.ID)
  print(paste("filename:", fileName))
  
  data <- read.table(fileName, header=T, na.strings=c("[Pending]","[Not Available]","[Not Applicable]","null","null ","NA"), sep="\t", quote="")
  samples <- as.character(data[,1])
  data <- data[,2:(ncol(data)-1)]
  rownames(data) <- samples
  print(paste(fileName, ": ", nrow(data), "(samples),", ncol(data), "(features)"))
  data <- na.roughfix(data) # impute the missing value in x, require randomForest library    
  return (data)
}
