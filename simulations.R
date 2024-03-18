# Chain X->Z->Y
library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source("./structs/stationary_distribution.R")
source("./structs/markov_simulation_src.R")
n = 2
d = 3
set.seed(1)
#define parent structure
ParentStructure = matrix(nrow = d, ncol = d, data = 0)
diag(ParentStructure) = 1
rownames(ParentStructure) = colnames(ParentStructure) = tail(LETTERS,d)
ParentStructure[2, 3] = 1
ParentStructure[3, 1] = 1
#ParentStructure[2,5] = 1
#ParentStructure[3,2] = 1
NoParents = rowSums(ParentStructure)

#-----------------------------------------------------------------------
#initialize class for markov_simulations
process = proc_init(n,d,ParentStructure)

process@trans_prob

#------------------------------------------------------------------------
#transition matrix
#takes some time to calculate, approx 30 sec for n=3 and d=4. 

attr(process,"trans_matrix") = trans_matrix(process)
attr(process,"statio_prob") = stationary_probability(process)

process@statio_prob
#------------------------------------------------------------------------
#markov process simulation
#m=10^7
m=10^3
attr(process,"simulation_trial") = markov_sim(process,m)

process@simulation_trial

X = process@simulation_trial[,"X"]
Z = process@simulation_trial[,"Z"]
Y = process@simulation_trial[,"Y"]

process@trans_prob

#------------------------------------------------------------------------------
# Stationary distribution of X
table(X,Y,Z) / m
PX = c(1 - process@trans_prob$X[2, "prob_1"], process@trans_prob$X[1, "prob_1"])
piX = PX / sum(PX)
piX

table(X,Y,Z)/m
sum(process@statio_prob[process@statio_prob["Y"]==1,]["statio_prob"])
sum(process@statio_prob[process@statio_prob["Y"]==0,]["statio_prob"])


sum(process@simulation_trial[,"X"] == 1 & process@simulation_trial[,"Y"] == 0 & process@simulation_trial[,"Z"] == 0 )
                             
                             
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

