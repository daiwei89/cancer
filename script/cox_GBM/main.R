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
require(synapseClient)

synapseLogin()

source('util.R')
source("models.R")

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
miRNA.data = scale(myRead(file.path(data.dir, 'GBM_miRNA_core.txt')))
#mRNA.data = scale(myRead(file.path(data.dir, 'GBM_mRNA_core.txt')))
#methy.data = scale(myRead(file.path(data.dir, 'GBM_methylation_core.txt')))
#CNV.data = scale(myRead(file.path(data.dir, 'GBM_CNV_core.txt')))
#exp.data = cbind(miRNA.data, mRNA.data, methy.data, CNV.data)
exp.data = cbind(miRNA.data)
print(paste('Total # of molecular features:', ncol(exp.data)))
surv.data <- myRead(file.path(data.dir, 'GBM_OS_core.txt')) # survival data
mySurv <- Surv(surv.data$OS_OS, surv.data$OS_vital_status,type='right')
print(paste("Total events:", length(which(mySurv[,"status"]==1))))

clinical.data <- myRead.simple(file.path(data.dir, 'GBM_clinical_core.txt'))
clinical.data$gender <- ifelse(clinical.data$gender=="FEMALE",1, 0)
clinical.data$performance_score = scale(clinical.data$performance_score)
clinical.data$age = scale(clinical.data$age)


###### specify the training and test
#train.all <- read.table(synapse.read(train.ID),header=F, stringsAsFactors=F)
#test.all <- read.table(synapse.read(test.ID),header=F, stringsAsFactors=F)
train.all <- read.table(file.path(data.dir, 'GBM_dev_train_sample_list.txt'),
header=F, stringsAsFactors=F)
test.all <- read.table(file.path(data.dir, 'GBM_dev_valid_sample_list.txt'),
header=F, stringsAsFactors=F)


# xgb hyper params
#eta.range = c(0.01, 0.1, 0.3, 0.5)
#lambda.range = c(2, 4, 6, 8, 10)
#num.rounds.range = c(1, 2, 4)
#alpha.range = c(0, 0.1, 1, 2, 4, 6, 8)
#max.depth.range = c(1, 3, 6, 9)

eta.range = c(0.1, 0.3)
lambda.range = c(0.1)
num.rounds.range = c(3)
alpha.range = c(0.1)
max.depth.range = c(3)
hyper.params = expand.grid(eta.range, lambda.range, num.rounds.range,
  alpha.range, max.depth.range)
print(paste('Explorint', nrow(hyper.params), 'hyper params'))

result <- c()
c.indexes.test <- c()
c.indexes.train <- c()
c.indexes.test.clinical <- c()
c.indexes.train.clinical <- c()
c.indexes.train.xgb = c()
c.indexes.test.xgb = c()
for (seed in 1:100) {
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

  print("***************************************************")
  print("***************************************************")
  print("LASSO + cox: ")

  result.cox <- c()
  list[c.index.train.clinical, c.index.test.clinical,
       c.index.train, c.index.test, result.cox] <- try(
    coxcv(x.train, y.train, x.test, y.test, clinical.train, clinical.test))
  if (class(result.cox)=="try-error") {
    result.cox <- rep(NA, nrow(x.test))
  }
  result <- cbind(result, as.numeric(result.cox))
  c.indexes.test = cbind(c.indexes.test, c.index.test)
  c.indexes.train = cbind(c.indexes.train, c.index.train)
  c.indexes.test.clinical = cbind(c.indexes.test.clinical, c.index.test.clinical)
  c.indexes.train.clinical = cbind(c.indexes.train.clinical, c.index.train.clinical)

  # xgb survival
  list[c.index.train.xgb, c.index.test.xgb] =
    xgb.surv(x.train, y.train, x.test, y.test, clinical.train, clinical.test,
    hyper.params)
  c.indexes.train.xgb = rbind(c.indexes.train.xgb, c.index.train.xgb)
  c.indexes.test.xgb = rbind(c.indexes.test.xgb, c.index.test.xgb)

  train.xgb.mean = colMeans(c.indexes.train.xgb)
  test.xgb.mean = colMeans(c.indexes.test.xgb)
  test.xgb.mean.max = max(test.xgb.mean, na.rm = T)
  print(paste('c.index.test.xgb.mean.max:', test.xgb.mean.max))
  print('======= best hyper param ==========')
  print(hyper.params[which(test.xgb.mean == test.xgb.mean.max),])
  print(paste('c.index.test.clinical:', mean(c.indexes.test.clinical)))
  print(paste('c.index.test.both:', mean(c.indexes.test)))
}
c.indexes.test = signif(c.indexes.test, 3)
c.indexes.train = signif(c.indexes.train, 3)
c.indexes.test.clinical = signif(c.indexes.test.clinical, 3)
c.indexes.train.clinical = signif(c.indexes.train.clinical, 3)

#print('dim of c.indexes.test.xgb')
#print(dim(c.indexes.test.xgb))
#print(colMeans(c.indexes.test.xgb))
# Find best hyperparam.
train.xgb.mean = colMeans(c.indexes.train.xgb)
test.xgb.mean = colMeans(c.indexes.test.xgb)
test.xgb.mean.max = max(test.xgb.mean)
print(paste('c.index.test.xgb.mean.max:', test.xgb.mean.max))
print('======= best hyper param ==========')
print(hyper.params[which(test.xgb.mean == test.xgb.mean.max),])

print(paste('c-index on test (clinical + molecular). mean:', mean(c.indexes.test),
  'median:', median(c.indexes.test),
  'Quantiles:', paste(quantile(c.indexes.test), collapse=", "),
  'std:', sd(c.indexes.test), 'num bootstraps:', length(c.indexes.test)))
print(paste('c-index on train (clinical + molecular). mean:', mean(c.indexes.train),
  'median:', median(c.indexes.train), 'Quantiles:',
  'Quantiles:', paste(quantile(c.indexes.train), collapse=", "),
  'std:', sd(c.indexes.train), 'num bootstraps:', length(c.indexes.train)))
print(paste('c-index on test (clinical). mean:', mean(c.indexes.test.clinical),
  'median:', median(c.indexes.test.clinical),
  'Quantiles:', paste(quantile(c.indexes.test.clinical), collapse=", "),
  'std:', sd(c.indexes.test.clinical),
  'num bootstraps:', length(c.indexes.test.clinical)))
print(paste('c-index on train (clinical). mean:', mean(c.indexes.train.clinical),
  'median:', median(c.indexes.train.clinical),
  'Quantiles:', paste(quantile(c.indexes.train.clinical), collapse=", "),
  'std:', sd(c.indexes.train.clinical),
  'num bootstraps:', length(c.indexes.train.clinical)))

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
