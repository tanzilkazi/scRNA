# Aim here is to use a supervised learning model for cell type classification
# contains: Part 1. simulation
#           Part 2. realworld data

# don't check mislabelling if occuring
# increase number of iterations
# perform binom 
# check the variances of the differences or improvements

setwd("C:/Users/kazit/Documents/USyd/COMP5703-Capstone")
library(amap)
library(mclust)
library(randomForest)
library(lemon)
library(cowplot)
library(ggplot2)
library(gridExtra)
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
num_groups <- 3

de.prob_min = 0.05
de.prob_max = 0.95
de.prob_inc = 0.05

num.groups_min = 2
num.groups_max = 4
num.groups_inc = 1

stats_to_capture = c("de.prob","num.groups","ARI.cluster.mean","ARI.cluster.sd","ARI.Ada.mean","ARI.Ada.sd")
recorder <-data.frame(matrix(ncol=length(stats_to_capture),nrow=0))
colnames(recorder) <- stats_to_capture
plots <- c()
de.prob_range <- seq(de.prob_min,de.prob_max,de.prob_inc)
num.groups_range <- seq(num.groups_min,num.groups_max,num.groups_inc)

for (num_groups in num.groups_range){
  for (de.prob_val in de.prob_range){
    #de.prob_val <- 0.5
    set.seed(2)
    params <- newSplatParams(batchCells=num_groups*50, group.prob = rep(1/num_groups, time=num_groups))
    params <- setParams(params, de.prob=de.prob_val)
    sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
    plots<- c(plots,plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")[1])
    dim(counts(sim.groups))
    
    sim.data <- log2(counts(sim.groups)+1)
    cls.truth <- colData(sim.groups)[,"Group"]
    cls.truth <- strtoi(substr(cls.truth,nchar(cls.truth),nchar(cls.truth)))
    
    sim.selected <- matPCs(sim.data, 10)
    
    ar_ARI_Ada <- c()
    ar_ARI_clu <- c()
    
    for (i in 1:10) {
      clust <- Kmeans(t(sim.selected), centers=num_groups, method="pearson", iter.max = 50)
      initCls <- clust$cluster
    
      model <- scAdaSampling(sim.selected, initCls, seed=i)
      final <- predict(model, newdata=t(sim.selected), type="response")
      
      #cat(i, "ada:", table(final==cls.truth), "clus", table(initCls==cls.truth),"\n")
      
      ar_ARI_Ada <- c(ar_ARI_Ada,adjustedRandIndex(cls.truth, final)) # after AdaSampling
      ar_ARI_clu <- c(ar_ARI_clu,adjustedRandIndex(cls.truth, initCls)) # initial cluster softlabeling
      #binom = binom.test(sum(ARI_AS>ARI_cluster), length(ARI_AS), alternative = "greater")
    }
    
    #print(ar_ARI_Ada)
    #print(ar_ARI_clu)
    
    row <- data.frame("de.prob" = de.prob_val,
                      "num.groups" = num_groups,
                      "ARI.cluster.mean" = format(round(mean(ar_ARI_clu),3)),
                      "ARI.cluster.sd" = format(round(sd(ar_ARI_clu),3)),
                      "ARI.Ada.mean" = format(round(mean(ar_ARI_Ada),3)),
                      "ARI.Ada.sd" = format(round(sd(ar_ARI_Ada),3)),
                      row.names=NULL)
    recorder <- rbind(recorder,row)
    cat("de.prob=",de.prob_val,"\tnum.groups=",num_groups,"\tCluster=",format(round(mean(ar_ARI_clu),3)),"+-",format(round(sd(ar_ARI_clu),3)),"\tAda=",format(round(mean(ar_ARI_Ada),3)),"+-",format(round(sd(ar_ARI_Ada),3)),"\n")
  }
}

recorder

# index=1
# grob_list=list()
# for (plott in plots){
#   if (index == 1){
#     legend <- get_legend(qplot(PC1, PC2, data = data.frame(plott), colour = colour_by)+ ggtitle(de.prob_range[index])+ theme(legend.position='bottom'))
#     temp <- qplot(PC1, PC2, data = data.frame(plott), colour = colour_by)+ ggtitle(de.prob_range[index])+ theme(legend.position='hidden') + theme(text=element_text(size=8),axis.text=element_text(size=8),plot.title=element_text(size=12))
#     grob_list[[index]] <- temp
#     }else{
#     temp <- qplot(PC1, PC2, data = data.frame(plott), colour = colour_by)+ ggtitle(de.prob_range[index])+ theme(legend.position='hidden') + theme(text=element_text(size=8),axis.text=element_text(size=8),plot.title=element_text(size=12))
#     grob_list[[index]] <- temp
#   }
#   index = index + 1
# }
# nCol <- floor(sqrt(length(grob_list)))
# all_plots <- grid.arrange(grobs=grob_list,ncol=nCol,legend)
# 
# recorder
# plot(all_plots)
# 
# save_string <- paste(format(Sys.time(),"%Y%m%d-%H%M%S"),"_G",num_groups,sep="")
# ggsave(paste("./code/",save_string,"_plots",".png",sep=""),plot=all_plots)
# write.csv(recorder,paste("./code/",save_string,"_perf",".csv",sep=""))

