# Chain X->Z->Y
library(dplyr)
library(data.table)
#setwd("C:/Users/Jacek/Desktop/Symulacje Magisterka")
source("./structs/functions.R")
n = 3
d = 4
#set vertices names
NodeNames = LETTERS[1:d]

#define relation structure
ParentStructure = matrix(nrow = d, ncol = d, data = 0)
diag(ParentStructure) = 1
rownames(ParentStructure) = colnames(ParentStructure) = NodeNames
ParentStructure[2, 1] = 1
ParentStructure[2, 4] = 1
#ParentStructure[2,5] = 1
#ParentStructure[3,2] = 1

NoParents = rowSums(ParentStructure)

TransitionProbabilities = Trans_prob(Nodenames, ParentStructure, n, d)

TransitionProbabilities[["A"]]

#########################################
#transition matrix
#takes some time to calculate, approx 30 sec for n=3 and d=4. 


Trans = trans_matrix(n, d)
stationary_probability(Trans,n,d)



#######
#markov process symulation

XZY = markov_sim(TransitionProbabilities,10^3,d,n)

X = XZY[, "A"]
Z = XZY[, "B"]
Y = XZY[, "C"]

hist(X)
# Stationary distribution of X
table(X) / m
PX = c(1 - TransitionProbabilities$X[2, "prob1"], TransitionProbabilities$X[1, "prob1"])
piX = PX / sum(PX)
piX

table(X, Y, Z)

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
