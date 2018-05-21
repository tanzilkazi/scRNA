library(splatter)
set.seed(1)
x <- 5
params <- newSplatParams(batchCells=x*50, group.prob = rep(1/x, time=x))
params <- setParams(params, de.prob=0.5)
sim.groups <- splatSimulate(params, method = "groups", verbose = FALSE)
plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")
dim(counts(sim.groups))
head(sim.groups,10)
