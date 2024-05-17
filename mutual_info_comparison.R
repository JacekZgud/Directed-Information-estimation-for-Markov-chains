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

prob=seq(from=1/2,to=0.99999,length.out=20)
info_estimator = function(a,par_struct,n_2=100){
  proc = markov_process_init(node_dim, nodes, par_struct, work_names)
  proc@trans_prob$X$prob_1 = c(1-a, a)
  proc@trans_prob$X$prob_0 = 1 - proc@trans_prob$X$prob_1
  
  proc@trans_prob$Y$prob_1 = c(1-a,a)
  proc@trans_prob$Y$prob_0 = 1 - proc@trans_prob$Y$prob_1
  print(proc@trans_prob)
  
  proc = trans_entropy(proc, c('Y'), n_2=n_2)

  return(proc@trans_entropy)
}
infos = Vectorize(info_estimator,'a')(a=prob,ParentStructure)
plot(prob,infos,ylim = c(0,1))
info_estimator(0.999999,ParentStructure,n_2=1000)


h = function(b){
  1-b*log(b/(1/2),base=2) - (1 - b)*log((1 - b)/(1/2),base=2)
}

mutual_info=function(x){
  a=x
  b=x
  c=a*(b^2) + a*((1 - b)^2) + (1 - a)*(b^2) + 2*(1 - a)*b*(1 - b) 
  c=1- ((1 - a)*b^2 + (1 - a)*(1 - b)^2 + 2*a*b*(1 - b))
  h(c)- h(x)
}
mutual_info(0.99)
dir_pol = function(x){
  1-h(x)
}
curve(dir_pol,from=0.5,to=1,add=TRUE,col='blue')
curve(mutual_info,from=0.5,to=1,add=TRUE,col='red')
Vectorize(mutual_info)(prob) > infos
