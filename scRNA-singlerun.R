# Aim here is to use a supervised learning model for cell type classification
# contains: Part 1. simulation
#           Part 2. realworld data

setwd("C:/Users/kazit/Documents/USyd/COMP5703 - Capstone")
library(amap)
library(mclust)
library(randomForest)
library(ggplot2)
library(gridExtra)
library(cowplot)

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
    for(j in 1:ncol(prob.mat)) {
      voteClass <- prob.mat[label==names(table(label))[j],]
      idx <- sample(1:nrow(voteClass), size=nrow(voteClass), replace = TRUE, prob=voteClass[,j])
      
      X <- cbind(X, data[, rownames(voteClass)[idx]])
      Y <- c(Y, label[rownames(voteClass)[idx]])
    }
  }
  return(model)
}



# PCA for dimension reduction
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
params <- newSplatParams(batchCells=x*50, group.prob = rep(1/x, time=x))
plots <- c()

params <- setParams(params, de.prob=y)
sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)


plots <- c(plots, plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")[1])
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
ARI
binom.test(sum(ARI>adjustedRandIndex(cls.truth, initCls)), length(ARI), alternative = "greater")
