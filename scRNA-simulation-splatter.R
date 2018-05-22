# Aim here is to use a supervised learning model for cell type classification
# contains: Part 1. simulation
#           Part 2. realworld data

setwd("C:/Users/kazit/Documents/USyd/COMP5703-Capstone")
library(amap)
library(mclust)
library(randomForest)
library(lemon)
suppressWarnings(options(warn=-1))


scAdaSampling <- function(data, label, seed) {
  set.seed(seed)
  X <- data
  Y <- label
  
  model <- c()
  prob.mat <- c()
  for (i in 1:5) {
    model <- randomForest(t(X), factor(Y), ntree = 500)
    prob.mat <- predict(model, newdata=t(data), type="prob")
    X <- c()
    Y <- c()
    for(j in 1:ncol(prob.mat)) {#j = soft label class
      voteClass <- prob.mat[label==names(table(label))[j],] # select the cells ONLY from class
      idx <- sample(1:nrow(voteClass), size=nrow(voteClass), replace = TRUE, prob=voteClass[,j]) # weighted sampling with replacement
      
      X <- cbind(X, data[, rownames(voteClass)[idx]])
      Y <- c(Y, label[rownames(voteClass)[idx]])
    }
  }
  
  return(model)
}

# PCA for dimension reduction
# multi = TRUE - bagging
library(gmodels)
matPCs <- function(data, top, seed=1, multi=FALSE) {
  # genes are rows, cells are cols
  pcs <- c()
  if (multi==TRUE) {
    set.seed(seed)
    n <- round(top / 5)
    size <- round(nrow(data) / 10)
   
     for (i in 1:5) {
      subData <-data[sample(1:nrow(data), size),]
      pcs <- cbind(pcs, fast.prcomp(subData, center = TRUE, scale. = TRUE)$rotation[,1:n])
    }
    return(t(pcs))
  } else {
    #pcs <- t(prcomp(data, center = TRUE, scale. = TRUE)$rotation[,1:top])
    pcs <- t(fast.prcomp(data, center = TRUE, scale. = TRUE)$rotation[,1:top])
  }  
}

library(splatter)
#### Alternatively, parameters can be set from scratch

set.seed(1)
x <- 3

de.prob_min = 0.14
de.prob_max = 0.15
de.prob_inc = 0.01

stats_to_capture = c("de.prob","binom","ARI.model.mean","ARI.model.sd","ARI.truth")
recorder <-data.frame(matrix(ncol=length(stats_to_capture),nrow=0))
colnames(recorder) <- stats_to_capture
plots <- c()
de.prob_range <- seq(de.prob_min,de.prob_max,de.prob_inc)

for (de.prob_val in de.prob_range){
  print(c("de.prob=",de.prob_val))
  params <- newSplatParams(batchCells=x*50, group.prob = rep(1/x, time=x))
  params <- setParams(params, de.prob=de.prob_val)
  sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
  plots<- c(plots,plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")[1])
  dim(counts(sim.groups))
  
  sim.data <- log2(counts(sim.groups)+1)
  cls.truth <- colData(sim.groups)[,"Group"]
  sim.selected <- matPCs(sim.data, 10)
  
  set.seed(2)
  clust <- Kmeans(t(sim.selected), centers=3, method="pearson", iter.max = 50)
  initCls <- clust$cluster
  adjustedRandIndex(cls.truth, initCls)
  
  ARI <- c()
  for(i in 1:10) {
    model <- scAdaSampling(sim.selected, initCls, seed=i)
    final <- predict(model, newdata=t(sim.selected), type="response")
    ARI <- c(ARI, adjustedRandIndex(cls.truth, final)) # after AdaSampling
  }
  binom = binom.test(sum(ARI>adjustedRandIndex(cls.truth, initCls)), length(ARI), alternative = "greater")
  row <- data.frame("de.prob" = de.prob_val,
                    "binom" = binom$estimate,
                    "ARI.model.mean" = mean(ARI),
                    "ARI.model.sd" = sd(ARI),
                    "ARI.truth" = adjustedRandIndex(cls.truth, initCls),
                    row.names=NULL)
  recorder <- rbind(recorder,row)
}

index=1
grob_list=list()
for (plott in plots){
  if (index == 1){
    temp <- qplot(PC1, PC2, data = data.frame(plott), colour = colour_by)+ ggtitle(de.prob_range[index])+ theme(legend.position="bottom")
    grob_list[[index]] <- temp
    legend <- get_legend(grob_list[[index]])
  }else{
    temp <- qplot(PC1, PC2, data = data.frame(plott), colour = colour_by)+ ggtitle(de.prob_range[index])+ theme(legend.position="none")
    grob_list[[index]] <- temp
  }
  index = index + 1
}
nCol <- floor(sqrt(length(grob_list)))
do.call("grid.arrange", c(grob_list, ncol=nCol))


#write.csv(recorder,paste("./code/",format(Sys.time(), "%Y%m%d-%H%M%S"),"rf",".csv",sep=""))

