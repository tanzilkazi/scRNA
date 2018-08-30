library(amap)
library(mclust)
library(randomForest)
library(lemon)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(stats)
library(splatter)

suppressWarnings(options(warn=-1))

now <- Sys.time()

setwd("C:/Users/kazit/Documents/USyd/COMP5703-Capstone")

time_print <-function(now,string){
  cat(string,"\t\t\t",Sys.time()-now,"\n")
  return (Sys.time())
}

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

set.seed(1)

de.prob_min = 0
de.prob_max = 1
de.prob_inc = 0.2

num.groups_min = 3
num.groups_max = 4
num.groups_inc = 1

stats_to_capture = c("de.prob","num.groups","ARI.cluster.mean","ARI.cluster.sd","ARI.Ada.mean","ARI.Ada.sd")
recorder <-data.frame(matrix(ncol=length(stats_to_capture),nrow=0))
colnames(recorder) <- stats_to_capture
de.prob_range <- seq(de.prob_min,de.prob_max,de.prob_inc)
num.groups_range <- seq(num.groups_min,num.groups_max,num.groups_inc)

now = time_print(now, "Initialization done...")

for (num_groups in num.groups_range){
  for (de.prob_val in de.prob_range){
    set.seed(2)
    sim.groups<-c()
    params <- newSplatParams(batchCells=num_groups*250, group.prob = rep(1/num_groups, time=num_groups,nGenes=100))
    params <- setParams(params, de.prob=de.prob_val)
    sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
    
    now = time_print(now, "Data simulation done...")
    
    #dim(counts(sim.groups))
    
    sim.data <- log2(counts(sim.groups)+1)
    cls.truth <- colData(sim.groups)[,"Group"]
    cls.truth <- strtoi(substr(cls.truth,nchar(cls.truth),nchar(cls.truth)))
    
    sim.selected <- matPCs(sim.data, 2)

    now = time_print(now, "PCA done...")
    
    ar_ARI_Ada <- c()
    ar_ARI_clu <- c()
    
    clust <- Kmeans(t(sim.selected), centers=num_groups, method="pearson", iter.max = 50)
    initCls <- clust$cluster
    
    now = time_print(now, "kmeans done...")
    
    for (i in 1:10) {
      model <- scAdaSampling(sim.selected, initCls, seed=i)
      final <- predict(model, newdata=t(sim.selected), type="response")
      ar_ARI_Ada <- c(ar_ARI_Ada,adjustedRandIndex(cls.truth, final)) # after AdaSampling
      ar_ARI_clu <- c(ar_ARI_clu,adjustedRandIndex(cls.truth, initCls)) # initial cluster softlabeling
    }
    now = time_print(now, "AdaSampling and prediction done...")
    
    row <- data.frame("de.prob" = de.prob_val,
                      "num.groups" = num_groups,
                      "ARI.cluster.mean" = format(round(mean(ar_ARI_clu),3)),
                      "ARI.cluster.sd" = format(round(sd(ar_ARI_clu),3)),
                      "ARI.Ada.mean" = format(round(mean(ar_ARI_Ada),3)),
                      "ARI.Ada.sd" = format(round(sd(ar_ARI_Ada),3)),
                      row.names=NULL)
    recorder <- rbind(recorder,row)
    cat("\t\tde.prob=",de.prob_val,"\tnum.groups=",num_groups,"\tCluster=",format(round(mean(ar_ARI_clu),3)),"+-",format(round(sd(ar_ARI_clu),3)),"\tAda=",format(round(mean(ar_ARI_Ada),3)),"+-",format(round(sd(ar_ARI_Ada),3)),"\n")
  }
}

clu<-as.numeric(as.character(recorder$ARI.cluster.mean)) # convert column in dataframe to list of simple arrays. factors are creaated if column created directly from dataframe
ada<-as.numeric(as.character(recorder$ARI.Ada.mean))

t.test(clu,ada,alternative="less",paired=TRUE)
save_string <- paste(format(Sys.time(),"%Y%m%d-%H%M%S"),"_Gr",num_groups,"_Gn",dim(counts(sim.groups))[1],sep="")
write.csv(recorder,paste("./code/",save_string,".csv",sep=""))

