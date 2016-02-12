##############################################################################
# Runs the cox proportional hazards model using nfolds cross validation and
# returns the predicted 
#
# nfolds = number of folds in the cross validation
# dir = directory where the scripts needed will be found (lasso.R)
#
# ** requires glmnet packages
##############################################################################
library(foreach)
library(doMC)
registerDoMC(6)
library(xgboost)
library(survcomp)
source("lasso.R")
source("correlation_screen.R")

# Return (X, score) restricted to uncensored data.
# x.train.all should include clinical + molecular 
prepare.rank <- function(X.train.both, y.train) {
  death = y.train[,2]  # 1 means uncensored (dead)
  uncensored.idx = which(death != 0, arr.ind = T)
  print(paste('Got', length(uncensored.idx), 'uncensored out of total',
    length(death)))
  list(X.train.both[uncensored.idx,], y.train[uncensored.idx,1])
}

compute.auc = function(pred, scores) {
  # Compute AUC
  N = length(pred)
  stopifnot(N == length(scores))
  num_pairs = 0
  num_correct = 0
  for (i in 1:(N-2)) {
    for (j in (i+1):(N-1)) {
      true_sign = scores[i] - scores[j]
      pred_sign = pred[i] - pred[j]
      if (true_sign * pred_sign >= 0) {
        num_correct = num_correct + 1
      }
      num_pairs = num_pairs + 1
    }
  }
  return (num_correct / num_pairs)
}

xgb.surv <- function(x.train,y.train, x.test, y.test, clinical.train,
  clinical.test, hyper.params) {
  start.time <- proc.time()
  all.train <- cbind(x.train, clinical.train)
  all.test <- cbind(x.test, clinical.test)
  # xgboost ranker
  all.train.mat = data.matrix(all.train)
  all.test.mat = data.matrix(all.test)
  #all.train.mat = data.matrix(clinical.train)
  #all.test.mat = data.matrix(clinical.test)
  list[X.rank, scores] = prepare.rank(all.train.mat, y.train)

  c.indexes.train.xgb = c()
  c.indexes.valid.xgb = c()
  for (h in 1:nrow(hyper.params)) { 
    eta = hyper.params[h, 1]
    lambda = hyper.params[h, 2]
    num.rounds = hyper.params[h, 3]
    alpha = hyper.params[h, 4]
    max.depth = hyper.params[h, 5]
    bst = xgboost(data = X.rank, label = scores, max.depth = max.depth,
      eta = eta, lambda = lambda, alpha = alpha,
      nround = num.rounds, nthread = 1, objective = "rank:pairwise", verbose = 0)
    pred.train.dead = predict(bst, X.rank)
    #auc = compute.auc(pred.train.dead, scores)
    # - sign cuz predict gives expected days to live, but we want who
    # dies first.
    pred.train = -predict(bst, all.train.mat)
    pred.valid = -predict(bst, all.test.mat)
    c.index.train.xgb <- concordance.index(pred.train, y.train[,1],
      y.train[,2])$c.index
    c.index.valid.xgb <- concordance.index(pred.valid, y.test[,1],
      y.test[,2])$c.index

    # dead only.
    #uncensored.idx = which(y.train[,2] != 0, arr.ind = T)
    #dead.train.X = all.train.mat[uncensored.idx,]
    #dead.train.y = y.train[uncensored.idx,]
    #dead.train.pred = predict(bst, dead.train.X)
    #c.index.train.xgb.dead <- concordance.index(dead.train.pred, dead.train.y[,1],
    #  dead.train.y[,2])$c.index

    c.indexes.train.xgb = cbind(c.indexes.train.xgb, c.index.train.xgb)
    c.indexes.valid.xgb = cbind(c.indexes.valid.xgb, c.index.valid.xgb)
    #print(paste('auc.dead = ', auc, 'c.idx.train = ', c.index.train.xgb,
    #  'c.index.train.xgb.dead:', c.index.train.xgb.dead,
    #  'c.index.valid ', c.index.valid.xgb))
    if (is.na(c.index.train.xgb) || is.na(c.index.valid.xgb)) {
      print('hyper params causing NA')
      print(hyper.params[h,])
      print('dumping to file: dump.raw.txt')
      xgb.dump(bst, "dump.raw.txt", with.stats = T)
      #print('pred.train')
      #print(pred.train)
      #print('X:')
      #print(X.rank)
      #print('score:')
      #print(scores)
    }
    #stopifnot(!is.na(c.index.train.xgb))
    #stopifnot(!is.na(c.index.valid.xgb))
  }
  c.index.train.xgb.max = max(c.indexes.train.xgb)
  c.index.valid.xgb.max = max(c.indexes.valid.xgb)
  print(paste('c.index.train.xgb.max:', c.index.train.xgb.max))
  print(paste('c.index.valid.xgb.max:', c.index.valid.xgb.max))
  #hyper.idx = which(c.indexes.valid.xgb == c.index.valid.xgb.max)
  #print('======= best hyper param from valid ===========')
  #print(hyper.params[hyper.idx,])
  runtime = proc.time() - start.time
  print(paste('xgb.surv finished in', runtime[1]))
  list(c.indexes.train.xgb, c.indexes.valid.xgb)
}

coxcv <- function(x.train,y.train, x.test, y.test, clinical.train,
  clinical.test, hyper.params, useLASSO=TRUE) {
  cox.clinical <- coxph(y.train~., data= clinical.train)
  print(paste("Number of non-clinical features:", ncol(x.train)))

  # clinical only
  cox.clinical$coefficients[is.na(cox.clinical$coefficients)]=0 
  cox.clinical.test <- as.matrix(clinical.test)%*%cox.clinical$coefficients
  cox.clinical.train <- as.matrix(clinical.train)%*%cox.clinical$coefficients
  c.index.train.clinical <- concordance.index(cox.clinical.train, y.train[,1],
    y.train[,2])$c.index
  c.index.test.clinical <- concordance.index(cox.clinical.test, y.test[,1],
    y.test[,2])$c.index
  
  # Add molecular
  clinical.resid <- residuals(cox.clinical)
  cols.include <- correlation.screen(clinical.resid, x.train, top=sum(y.train[,2]))
  
  if (length(cols.include)==0) {
    stop("No feature passed the univariate cox screen: exit.")
  }
  x.train <- x.train[,cols.include]
  x.test <- x.test[,cols.include]

  print(paste("After univariate cox screen, features remain:", length(cols.include)))

  # further shrink by LASSO, if only few features, no need to use LASSO, note
  # 5 is a quite arbitrary setting
  if(useLASSO & length(cols.include) > 5) {
    # do LASSO without cross validation to get the features to include in the model
    # change x.train and x.test to only include those features
    cols.include <- c()
    iter <- 0
    while(length(cols.include)<1) {
      set.seed(iter+1)
      cols.include <- try(lasso(x=x.train,y=clinical.resid,above=0))
      if (class(cols.include)=="try-error") {
          print("Errors occur while calculating by LASSO, recalculating...")
          cols.include <- c()
        }
      iter <- iter+1
      if(iter> 1) {
        print(paste(length(cols.include),
          " features selected. Recalculated by LASSO:", iter))
      }
      # maximum number of iterations allowed 
      if(iter>100) {
        stop("No significant features can be selected by LASSO: exit.")
      }
    }    
    print(paste("After LASSO, features remain:", length(cols.include)))
    x.train <- x.train[,cols.include]
    x.test <- x.test[,cols.include]
  }
  #print("------------------------------------------------")
  #print(colnames(x.train))
  #print("------------------------------------------------")

  # cox model for prediction
  if (length(cols.include)==1) {
    x.train <- data.frame(x.train)
    colnames(x.train) <- "genomic"
    x.test <- data.frame(x.test)
    colnames(x.test) <- "genomic"
  }

  all.train <- cbind(x.train, clinical.train)
  all.test <- cbind(x.test, clinical.test)
  # cox model from clinical + molecular
  cox.both <- coxph(y.train~., data= all.train)

  # clinical + molecular
  cox.both$coefficients[is.na(cox.both$coefficients)]=0 
  cox.both.test <- as.matrix(all.test)%*%cox.both$coefficients
  cox.both.train <- as.matrix(all.train)%*%cox.both$coefficients
  c.index.train.both <- concordance.index(cox.both.train, y.train[,1],
    y.train[,2])$c.index
  c.index.test.both <- concordance.index(cox.both.test, y.test[,1],
    y.test[,2])$c.index
	list(c.index.train.clinical, c.index.test.clinical,
       c.index.train.both, c.index.test.both, cox.both.test)
}
