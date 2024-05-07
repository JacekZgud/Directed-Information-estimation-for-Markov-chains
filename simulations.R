# Chain X->Z->Y
library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source("./structs/stationary_distribution.R")
source("./structs/markov_simulation_src.R")
source("./structs/simulation_marginalized.R")
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
ParentStructure[3, 1] = 1
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
process@trans_prob$X$prob_1 = c(0.05, 0.05, 0.95, 0.95)
process@trans_prob$X$prob_0 = 1 - process@trans_prob$X$prob_1

process@trans_prob$Y$prob_1 = c(1, 0, 1, 0)
process@trans_prob$Y$prob_0 = 1 - process@trans_prob$Y$prob_1

process@trans_prob$X$prob_1 = c(0, 0.5, 0.5, 1)
process@trans_prob$X$prob_0 = 1 - process@trans_prob$X$prob_1

process@trans_prob$Y$prob_1 = c(0.5, 0.5, 0.5, 0.5)
process@trans_prob$Y$prob_0 = 1 - process@trans_prob$Y$prob_1

#------------------------------------------------------------------------
#transition matrix
#takes some time to calculate, approx 30 sec for n=3 and d=4.

process@statio_prob = stationary_probability(process)
trans_matrix(process, TRUE)
class(process@statio_prob)
process@trans_matrix_list
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
# Stationary distribution of X

table(Z, Y) / m


XY.Z0 = table(X[Z == 0], Y[Z == 0])
XY.Z0
mosaicplot(XY.Z0)
chisq.test(XY.Z0)$expected
XY.Z1 = table(X[Z == 1], Y[Z == 1])
XY.Z1
mosaicplot(XY.Z1)
chisq.test(XY.Z1)$expected

XY. = table(X, Y)
XY.
mosaicplot(XY.)
chisq.test(XY.)$expected


# simulate marginalized markov process
n_2 = 10 ^ 2

process = marginalized_runner(process, c('X'), n_2)

process@trans_prob
process = trans_entropy(process, c('X'), n_2=n_2)

process@trans_prob

plot(c(1:n_2), out$ft[1, -c(1:3)], pch = ".")
plot(c(1:n_2), out$mt[1, -1], pch = ".")
hist(out$ft[out$y == 0], breaks = seq(0, 1, length.out = 10))
hist(out$ft[out$y == 1], breaks = seq(0, 1, length.out = 10))
out$ft[, c(1:10)]
out$ft[(out$ft == 2)]

to = markov_sim_Y(process, 1000, c('Y', 'X', 'W'))
table(to$y[, 'W'])
data.table(process@statio_prob)[, sum(statio_prob), by = eval("W")]
process@trans_prob
plot(c(1:1000), to$ft[1, -c(1:2)], pch = ".")
to$ft[, 1:20]
process@trans_prob[['Y']]
data.table(process@statio_prob)[, sum(statio_prob), by = eval(c('Y', 'Z'))]
process@trans_prob
