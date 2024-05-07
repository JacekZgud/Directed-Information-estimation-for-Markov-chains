# example
library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source("./structs/stationary_distribution.R")
source("./structs/markov_simulation_src.R")
source("./structs/simulation_marginalized.R")
source('./structs/information_calculation.R')
node_dim = 2
nodes = 2
set.seed(1)
#define parent structure
ParentStructure = matrix(nrow = nodes, ncol = nodes, data = 0)
work_names = c("X", "Y")
rownames(ParentStructure) = colnames(ParentStructure) = work_names
#ParentStructure[1, 2] = 1
ParentStructure[2, 1] = 1
ParentStructure[1,1] = 1

prob=seq(from=1/2,to=0.99999,length.out=6)
info_estimator = function(a,par_struct){
  process = markov_process_init(node_dim, nodes, par_struct, work_names)
  process@trans_prob$X$prob_1 = c(1-a, a)
  process@trans_prob$X$prob_0 = 1 - process@trans_prob$X$prob_1
  
  process@trans_prob$Y$prob_1 = c(1-a,a)
  process@trans_prob$Y$prob_0 = 1 - process@trans_prob$Y$prob_1
  
  
  n_2=1000
  process = trans_entropy(process, c('Y'), n_2=n_2)
  print(a)
  return(process@trans_entropy)
}
Vectorize(info_estimator,'a')(a=prob,ParentStructure)

