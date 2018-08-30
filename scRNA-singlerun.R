rm(list=ls(all=TRUE))

library(amap)
library(mclust)
library(randomForest)
library(lemon)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(stats)
library(splatter)
library(gmodels)

suppressWarnings(options(warn=-1))

setwd("C:/Users/kazit/Documents/USyd/COMP5703-Capstone")

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

set.seed(1)

stats_to_capture = c("de.prob","num.groups","ARI.cluster.mean","ARI.cluster.sd","ARI.Ada.mean","ARI.Ada.sd")
recorder <-data.frame(matrix(ncol=length(stats_to_capture),nrow=0))
colnames(recorder) <- stats_to_capture

num_groups <- 2
de.prob_val <- 0.3

set.seed(2)
params <- newSplatParams(batchCells=num_groups*50, group.prob = rep(1/num_groups, time=num_groups),nGenes=10000)
params <- setParams(params, de.prob=de.prob_val)
sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)

#cat("# of Genes: ",dim(counts(sim.groups))[1],"\n")
#cat("# of Cells: ",dim(counts(sim.groups))[2],"\n")

sim.data <- log2(counts(sim.groups)+1)
cls.truth <- colData(sim.groups)[,"Group"]
cls.truth <- strtoi(substr(cls.truth,nchar(cls.truth),nchar(cls.truth)))

sim.selected <- matPCs(sim.data, 10)

graph_data <-data.frame(t(sim.selected))

p <- plot_ly(graph_data, x =graph_data$PC1, y =graph_data$PC2, z =graph_data$PC3, color = sim.groups$Group) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

p

ar_ARI_Ada <- c()
ar_ARI_clu <- c()

clust <- kmeans(t(sim.selected), centers=num_groups, iter.max = 50)
initCls <- clust$cluster

# ## averaged kmeans
# softclasses <-c()
# 
# for (i in 1:10) {
# 
#   set.seed(i)
#   clust <- Kmeans(t(sim.selected), centers=num_groups, method="pearson", iter.max = 50)
#   softcl <- clust$cluster
#   softclasses[[i]] <- softcl
# 
# }
# softclass_df <- do.call(rbind, softclasses)
# initCls <- c()
# 
# for(i in 1:ncol(softclass_df)){
# 
#   t<- table(softclass_df[,i])
#   test_cls<- as.numeric(names(t)[which.max(t)])
#   initCls <- c(initCls, test_cls)
# 
# }
# 
# names(initCls)<-colnames(sim.data)
# ##
    
for (i in 1:10) {
  model <- scAdaSampling(sim.selected, initCls, seed=i)
  final <- predict(model, newdata=t(sim.selected), type="response")
  
  ar_ARI_Ada <- c(ar_ARI_Ada,adjustedRandIndex(cls.truth, final)) # after AdaSampling
  ar_ARI_clu <- c(ar_ARI_clu,adjustedRandIndex(cls.truth, initCls)) # initial cluster softlabeling
}
    
cat("de.prob=",de.prob_val,"\tnum.groups=",num_groups,"\tCluster=",format(round(mean(ar_ARI_clu),3)),"+-",format(round(sd(ar_ARI_clu),3)),"\tAda=",format(round(mean(ar_ARI_Ada),3)),"+-",format(round(sd(ar_ARI_Ada),3)),"\n")

