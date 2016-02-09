#############################################################################
# Yuan Yuan, 2013-04-09
# This is the main function for predicting the clinical outcome from molecular data
# The code was modified to fetch data directly from Synapse
# Cancer: LUSC
# type: RPPA
# outcome: overall survival
#############################################################################
rm(list=ls())
library(survival)
library(glmnet)
library(foreach)
library(doMC)
registerDoMC(6)
#require(synapseClient)

#synapseLogin()

source('util.R')

##### specify the parameters #########
cancer <-  "GBM"
platform <- "miRNA+clinical"
method <- "Cox"
outdir <- "./"
exp.ID <- "syn1710368" # GBM miRNA
surv.ID <- "syn1710370" # GBM survival
train.ID <- "syn1714087"
test.ID <- "syn1714083"

clinical.ID <- "syn1715822"

print(paste("Performing survial analysis for", cancer, sep=": "))

data.dir = "/Users/daiwei89/storage/cs_proj/dap/data"


###### read in the data ##############
source("read_table.R")
exp.data <- myRead(file.path(data.dir, 'GBM_miRNA_core.txt')) # expression data
surv.data <- myRead(file.path(data.dir, 'GBM_OS_core.txt')) # survival data
mySurv <- Surv(surv.data$OS_OS, surv.data$OS_vital_status,type='right')
print(paste("Total events:", length(which(mySurv[,"status"]==1))))

clinical.data <- myRead.simple(file.path(data.dir, 'GBM_clinical_core.txt'))
clinical.data$gender <- ifelse(clinical.data$gender=="FEMALE",1, 0)


###### specify the training and test
#train.all <- read.table(synapse.read(train.ID),header=F, stringsAsFactors=F)
#test.all <- read.table(synapse.read(test.ID),header=F, stringsAsFactors=F)
train.all <- read.table(file.path(data.dir, 'GBM_train_sample_list.txt'),
header=F, stringsAsFactors=F)
test.all <- read.table(file.path(data.dir, 'GBM_test_sample_list.txt'),
header=F, stringsAsFactors=F)


result <- c()
for (seed in 1:1)
  {
    print(paste("seed =", seed))
    train.samples <- train.all[,seed]
    test.samples <- test.all[,seed]

    train.row <- match(train.samples, rownames(exp.data))
    test.row <- match(test.samples, rownames(exp.data))

    x.train <- exp.data[train.row, ]
    y.train <- mySurv[train.row]

    x.test <- exp.data[test.row, ]
    y.test <- mySurv[test.row]

    clinical.train <- clinical.data[train.row, ]
    clinical.test <- clinical.data[test.row,]

#####################################
## random survival forest model
#####################################
    source("coxcv.R")
    print("***************************************************")
    print("***************************************************")
    print("LASSO + cox: ")

    result.cox <- c()
    list[c.index, result.cox] <- try(coxcv(x.train, y.train, x.test, y.test, clinical.train, clinical.test))
    if (class(result.cox)=="try-error")
      {
        result.cox <- rep(NA, nrow(x.test))
      }
    result <- cbind(result, as.numeric(result.cox))
    c.indexes = cbind(c.indexes, c.index)
  }
print(paste('c-index. mean:', mean(c.indexes), 'std:'. sd(c.indexes),
'num bootstraps:', length(c.indexes)))

#######################################

outfile <- paste(cancer,"_", platform,"_", method,".tsv", sep="")
write(t(result),file=outfile, ncolumns=100, sep="\t")

#Save this code object into Synapse

#Save this code object into Synapse

#myCode <- createEntity(Code(list(name=paste(method, cancer, platform), parentId="syn1720423")))
#myCode <- addFile(myCode, "main.R")
#myCode <- addFile(myCode, "correlation_screen.R")
#myCode <- addFile(myCode, "coxcv.R")
#myCode <- addFile(myCode, "lasso.R")
#myCode <- addFile(myCode, "read_table.R")
#myCode <- storeEntity(myCode)
#
#myResults <- createEntity(Data(name=paste(method, cancer, platform), parentId='syn1720419'))
#used(myResults)<-list(list(entity=myCode, wasExecuted=T),
#                      list(entity=train.ID, wasExecuted=F),
#                      list(entity=test.ID, wasExecuted=F),
#                      list(entity=exp.ID, wasExecuted=F),
#                      list(entity=surv.ID, wasExecuted=F))
#myResults<-addFile(myResults, outfile)
#myResults$annotations$cancer=cancer
#myResults$annotations$dataType= platform
#myResults$annotations$method = "LASSO + cox" #'Random survival forest' # # # 
#myResults$annotations$normalization = 'None'
#myResults$annotations$featureSelection='LASSO'
#myResults$annotations$clinicalUsed = 'Yes'
#
#myResults<-storeEntity(myResults)
#evaluation<-Evaluation(synRestGET(sprintf("/evaluation/%s", '1876290')))
#submit(evaluation, myResults)  #where myResults is a file Entity
