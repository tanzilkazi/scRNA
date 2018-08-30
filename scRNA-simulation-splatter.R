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

de.prob_min = 0
de.prob_max = 1
de.prob_inc = 0.1

num.groups_min = 2
num.groups_max = 4
num.groups_inc = 1

stats_to_capture = c("de.prob","num.groups","ARI.cluster.mean","ARI.cluster.sd","ARI.Ada.mean","ARI.Ada.sd")
deprob_num_iter <- 1000
result_recorder <-data.frame(matrix(ncol=length(stats_to_capture),nrow=0))
colnames(result_recorder) <- stats_to_capture
de.prob_range <- seq(de.prob_min,de.prob_max,de.prob_inc)
num.groups_range <- seq(num.groups_min,num.groups_max,num.groups_inc)

seed_counter = 0
cat("deprob range:",de.prob_range,"\n")
cat("groups range:",num.groups_range,"\n")

for (num_groups in num.groups_range){
  for (de.prob_val in de.prob_range){
    avg_ARI_ada <- c()
    avg_ARI_clu <- c()
    cat ("deporb: ",de.prob_val,"\niteration:")
    temp <-0
    for (iter in 1:deprob_num_iter){
      cat(" ",iter)
      params <- newSplatParams(batchCells=num_groups*50, group.prob = rep(1/num_groups, time=num_groups),nGenes=10000,seed=seed_counter)
      params <- setParams(params, de.prob=de.prob_val)
      sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
      
      #cat("# of Genes: ",dim(counts(sim.groups))[1],"\n")
      #cat("# of Cells: ",dim(counts(sim.groups))[2],"\n")
      
      sim.data <- log2(counts(sim.groups)+1)
      cls.truth <- colData(sim.groups)[,"Group"]
      cls.truth <- strtoi(substr(cls.truth,nchar(cls.truth),nchar(cls.truth)))
      
      sim.selected <- matPCs(sim.data, 10)
      
      ar_ARI_Ada <- c()
      ar_ARI_clu <- c()
      
      clust <- Kmeans(t(sim.selected), centers=num_groups, method="pearson", iter.max = 50)
      initCls <- clust$cluster

      ## averaged kmeans
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
        
        #cat(i, "ada:", table(final==cls.truth), "clus", table(initCls==cls.truth),"\n")
        
        ar_ARI_Ada <- c(ar_ARI_Ada,adjustedRandIndex(cls.truth, final)) # after AdaSampling
        ar_ARI_clu <- c(ar_ARI_clu,adjustedRandIndex(cls.truth, initCls)) # initial cluster softlabeling
        #binom = binom.test(sum(ARI_AS>ARI_cluster), length(ARI_AS), alternative = "greater")
      }
      # cat ("ada: ",ar_ARI_Ada,"\n")
      # cat ("clu: ",ar_ARI_clu,"\n")

      avg_ARI_ada <- c(avg_ARI_ada,mean(ar_ARI_Ada))
      avg_ARI_clu <- c(avg_ARI_clu,mean(ar_ARI_clu))
      # cat ("avg_ada: ",avg_ARI_ada,"\n")
      # cat ("avg_clu: ",avg_ARI_clu,"\n")
      seed_counter<-seed_counter + 1
    } 

    row <- data.frame("de.prob" = de.prob_val,
                      "num.groups" = num_groups,
                      "ARI.cluster.mean" = format(round(mean(avg_ARI_clu),3)),
                      "ARI.cluster.sd" = format(round(sd(avg_ARI_clu),3)),
                      "ARI.Ada.mean" = format(round(mean(avg_ARI_ada),3)),
                      "ARI.Ada.sd" = format(round(sd(avg_ARI_ada),3)),
                      row.names=NULL)
    result_recorder <- rbind(result_recorder,row)
    cat("\nde.prob=",de.prob_val,"\tnum.groups=",num_groups,"\tCluster=",format(round(mean(avg_ARI_clu),3)),"+-",format(round(sd(avg_ARI_clu),3)),"\tAda=",format(round(mean(avg_ARI_ada),3)),"+-",format(round(sd(avg_ARI_ada),3)),"\n")
  }
}

clu<-as.numeric(as.character(result_recorder$ARI.cluster.mean)) # convert column in dataframe to list of simple arrays. factors are creaated if column created directly from dataframe
ada<-as.numeric(as.character(result_recorder$ARI.Ada.mean))

t.test(clu,ada,alternative="less",paired=TRUE)

save_string <- paste(format(Sys.time(),"%Y%m%d-%H%M%S"),"_Gr",num_groups,"_1000deprob_iter",sep="")
write.csv(deprob_iter_recorder,paste("./code/",save_string,".csv",sep=""))

save_string <- paste(format(Sys.time(),"%Y%m%d-%H%M%S"),"_Gr",num_groups,"_1000deprob_results",sep="")
write.csv(result_recorder,paste("./code/",save_string,".csv",sep=""))
