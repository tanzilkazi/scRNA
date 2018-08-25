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

num_groups=7
de.prob_val = 0.5

params <- newSplatParams(batchCells=num_groups*50, group.prob = rep(1/num_groups, time=num_groups))
params <- setParams(params, de.prob=de.prob_val)
sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
dim(counts(sim.groups))

sim.data <- log2(counts(sim.groups)+1)
cls.truth <- colData(sim.groups)[,"Group"]
cls.truth <- strtoi(substr(cls.truth,nchar(cls.truth),nchar(cls.truth)))

sim.selected <- matPCs(sim.data, 3)
graph_data <-data.frame(t(sim.selected))

# 2D plot
# pl <- ggplot(graph_data,aes(x=graph_data$PC1, y=graph_data$PC2)) + geom_point(aes(col=sim.groups$Group))
# plot(pl)



library(plotly)

mtcars$am[which(mtcars$am == 0)] <- 'Automatic'
mtcars$am[which(mtcars$am == 1)] <- 'Manual'
mtcars$am <- as.factor(mtcars$am)

p <- plot_ly(graph_data, x =graph_data$PC1, y =graph_data$PC2, z =graph_data$PC3, color = sim.groups$Group) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

p
