# This is a R file that represents all functionalities of the repository.

library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source('./structs/information_calculation.R')


# parametrise the code
node_dim = 3
nodes = 3
work_names = tail(LETTERS, nodes)

#define parent structure
ParentStructure = matrix(nrow = nodes, ncol = nodes, data = 0)
rownames(ParentStructure) = colnames(ParentStructure) = work_names

# define the parent structure of nodes
diag(ParentStructure) = 1
ParentStructure[2, 3] = 1
ParentStructure[3, 2] = 1
ParentStructure[1, 2] = 1
ParentStructure[2, 1] = 1

#------------------------------------------------------------------------
#initialize class for markov_simulations
process = markov_process_init(node_dim, nodes, ParentStructure, work_names)
process@trans_prob

#example of quick manual conditional probability changes:
process@trans_prob$X
process@trans_prob$X$prob_0 = c(0.05, 0.05, 0.90, 0.90, 0.05, 0.05, 0.90, 0.90, 0.1)
process@trans_prob$X$prob_1 = c(0.90, 0.90, 0.05, 0.05, 0.90, 0.90, 0.05, 0.05, 0.1)
process@trans_prob$X$prob_2 = 1 - process@trans_prob$X$prob_1 -  process@trans_prob$X$prob_0
process@trans_prob$X



#------------------------------------------------------------------------
#transition matrix
#takes some time to calculate, approx 30 sec for n=3 and d=4.

process@statio_prob = stationary_probability(process)

process@trans_matrix_list = trans_matrix(process, TRUE)


#------------------------------------------------------------------------
#markov process simulation
#m=10^7
m = 10 ^ 2

process = simulation_runner(process, m)
process@simulation

#-------------------------------------------------------------------------
# simulate marginalized markov process
n_2 = 10 ^ 2

process = marginalized_runner(process, c('Y', 'Z'), n = n_2)

process@marg_sim

#calculate transfer entropy from V/target -----> target
process = trans_entropy(process, c('Y', 'Z'), n = n_2)
