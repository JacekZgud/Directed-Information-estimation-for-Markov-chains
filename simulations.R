# Chain X->Z->Y
library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source('./structs/information_calculation.R')

node_dim = 3
nodes = 3
set.seed(1)
#define parent structure
ParentStructure = matrix(nrow = nodes, ncol = nodes, data = 0)
diag(ParentStructure) = 1
work_names = c("X", "Y")
work_names = tail(LETTERS, nodes)
rownames(ParentStructure) = colnames(ParentStructure) = work_names

ParentStructure[2, 3] = 1
#ParentStructure[3, 1] = 1
ParentStructure[3, 2] = 1
# ParentStructure[2,5] = 1
# ParentStructure[3,2] = 1
# ParentStructure[1,2] = 1
ParentStructure[1, 2] = 1
ParentStructure[2, 1] = 1

#------------------------------------------------------------------------
#initialize class for markov_simulations
process = markov_process_init(node_dim, nodes, ParentStructure, work_names)
process@trans_prob

#examples:
process@trans_prob$X$prob_1 = c(0.05, 0.05, 0.95, 0.95)
process@trans_prob$X$prob_0 = 1 - process@trans_prob$X$prob_1

process@trans_prob$Y$prob_1 = c(0.45, 0.45, 0.55, 0.55)
process@trans_prob$Y$prob_0 = 1 - process@trans_prob$Y$prob_1

process@trans_prob$X$prob_1 = c(0, 0.5, 0.5, 1)
process@trans_prob$X$prob_0 = 1 - process@trans_prob$X$prob_1

process@trans_prob$Y$prob_1 = c(0.5, 0.5, 0.5, 0.5)
process@trans_prob$Y$prob_0 = 1 - process@trans_prob$Y$prob_1

#------------------------------------------------------------------------
#transition matrix
#takes some time to calculate, approx 30 sec for n=3 and d=4.

process@statio_prob = stationary_probability(process)
class(process@statio_prob)
process@trans_matrix_list = trans_matrix(process, TRUE)

#------------------------------------------------------------------------
#markov process simulation
#m=10^7
m = 10 ^ 4
#attr(process, "simulation_trial") = markov_sim(process, m)

process = simulation_runner(process, m)
process@simulation

X = process@simulation[, "X"]
Z = process@simulation[, "Z"]
Y = process@simulation[, "Y"]
W = process@simulation[, "W"]


process@trans_prob

#-------------------------------------------------------------------------
# simulate marginalized markov process
n_2 = 10 ^ 3

process = marginalized_runner(process, c('Y', 'Z'), n = n_2)

process@trans_prob

#calculate transfer entropy from V/target -----> target
process = trans_entropy(process, c('Y', 'Z'), n = n_2)
