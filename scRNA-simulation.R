# Aim here is to use a supervised learning model for cell type classification
# contains: Part 1. simulation
#           Part 2. realworld data

setwd("C:/Users/kazit/Documents/USyd/COMP5703-Capstone")
library(amap)
library(mclust)
library(randomForest)
library(ggplot2)

scAdaSampling <- function(data, label, seed,num_tree,mtry_val){
    set.seed(seed)
    X <- data
    Y <- label
  
    model <- c()
    prob.mat <- c()
    for (i in 1:5) {
        model <- randomForest(t(X), factor(Y),ntree=num_tree,mtry=mtry_val)
        prob.mat <- predict(model, newdata=t(data), type="prob")
        X <- c()
        Y <- c()
        for(j in 1:ncol(prob.mat)) { #j = soft label class
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
params <- newSplatParams(batchCells=x*50, group.prob = rep(1/x, time=x))
params <- setParams(params, de.prob=0.15)
sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")
dim(counts(sim.groups))


sim.data <- log2(counts(sim.groups)+1)
cls.truth <- colData(sim.groups)[,"Group"]
cls.truth <- as.numeric(substr(cls.truth, nchar(cls.truth[1]), nchar(cls.truth[1])))

### simulation for clustering and classification
#genes <- selectGenes(sim.data, cls.noisy, top=100)
#sim.selected <- sim.data[genes,]
sim.selected <- matPCs(sim.data, 10)

set.seed(2)
clust <- Kmeans(t(sim.selected), centers=3, method="pearson", iter.max = 50)
initCls <- clust$cluster
adjustedRandIndex(cls.truth, initCls)

#visualize kmeans
library(cluster)
clusplot(t(sim.selected),clust$cluster,color=TRUE,share=TRUE,lines=0)


recorder <-data.frame(matrix(ncol=6,nrow=0))
colnames(recorder) <- c("ntree","mtry","ARI.model.mean","ARI.model.sd","ARI.truth","binom")

ntree_min = 50
ntree_max = 1000
ntree_inc = 50

num_features = dim(sim.selected)[1]
mtry_list = round(c(0.5*sqrt(num_features),sqrt(num_features),2*sqrt(num_features)))

for (mtry_val in mtry_list){
    for (ntree in seq(ntree_min,ntree_max,ntree_inc)) {
        ARI <- c()
        print(c("mtry=",mtry_val,"ntree=",ntree))
        for(i in 1:10) {
            model <- scAdaSampling(sim.selected, initCls, seed=i,ntree,mtry_val)
            final <- predict(model, newdata=t(sim.selected), type="response")
            ARI <- c(ARI, adjustedRandIndex(cls.truth, final)) # after AdaSampling
        }
        binom = binom.test(sum(ARI>adjustedRandIndex(cls.truth, initCls)), length(ARI), alternative = "greater")
        row <- data.frame("ntree"=ntree,
                          "mtry"=mtry_val,
                          "ARI.model.mean" = mean(ARI),
                          "ARI.model.sd" = sd(ARI),
                          "ARI.truth" = adjustedRandIndex(cls.truth, initCls),
                          "binom" = binom$estimate,
                          row.names=NULL)
        recorder <- rbind(recorder,row)
    }
}

write.csv(recorder,paste(format(Sys.time(), "%H%M%S"),"rf",".csv",sep=""))
