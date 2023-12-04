# Chain X->Z->Y
library(dplyr)
library(data.table)
#setwd("C:/Users/Jacek/Desktop/Symulacje Magisterka")
source("./structs/functions.R")
n = 3
d = 4
NodeNames = LETTERS[1:d]
ParentStructure = matrix(nrow = d, ncol = d, data = 0)
diag(ParentStructure) = 1
rownames(ParentStructure) = colnames(ParentStructure) = NodeNames
ParentStructure[2, 1] = 1
ParentStructure[2, 4] = 1
#ParentStructure[2,5] = 1
#ParentStructure[3,2] = 1

NoParents = rowSums(ParentStructure)
TransitionProbabilities = vector("list", c(d))
names(TransitionProbabilities) = NodeNames

Trans_prob = function(Nodenames, ParentStructure, n, d) {
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  for (node in NodeNames) {
    parents = NodeNames[which(ParentStructure[node,] == 1)]
    n0_parents = NoParents[parents]
    temp = expand.grid(rep(list(0:c(n - 1)), length(n0_parents)))
    colnames(temp) = parents
    #temp['prob1'] = runif(nrow(temp))
    #print(runif(n*nrow(temp)))
    prob_matr = matrix(
      runif(n * nrow(temp)),
      nrow = nrow(temp),
      ncol = n,
      byrow = TRUE
    )
    colnames(prob_matr) = prob_vector_names
    prob_matr = t(apply(prob_matr, 1, function(x)
      x / sum(x)))
    TransitionProbabilities[[node]] = data.table(cbind(temp, prob_matr))
    #print(as.matrix(cbind(temp,prob_matr)))
    #print(data.table(cbind(temp,prob_matr)))
    setkeyv(TransitionProbabilities[[node]], parents)
  }
  return(TransitionProbabilities)
  
}
TransitionProbabilities = Trans_prob(Nodenames, ParentStructure, n, d)

TransitionProbabilities[["A"]]

#########################################
#Macierz przejśćia



start = Sys.time()
Trans = trans_matrix(n, d)
end = Sys.time()
end - start

rowSums(Trans)
A = t(Trans - diag(ncol = n ^ d, nrow = n ^ d))
A = rbind(A, rep(1, n ^ d))
b = c(rep(0, n ^ d), 1)
res_statio = qr.solve(A, b)
res_statio = cbind(expand.grid(rep(list(0:c(
  n - 1
)), d)), res_statio)
colnames(res_statio) = c(NodeNames, 'statio_prob')
sum(res_statio$statio_prob)
res_statio

# for (vertex in 1:d)
# {TransitionProbabilities[[vertex]] = numeric(2^NoParents[vertex])}
#
# TransitionProbabilities[[1]]=runif(2)
# TransitionProbabilities[[2]]=runif(4)
# TransitionProbabilities[[3]]=runif(4)

#TransitionProbabilities[[1]]=c(0.01,0.99)
#TransitionProbabilities[[2]]=c(0.1,0.1,0.9,0.9)
#TransitionProbabilities[[3]]=c(0.1,0.1,0.9,0.9)

m = 10 ^ 4
XZY = matrix(nrow = m, ncol = d)
colnames(XZY) = NodeNames
State = rep(0, d)
for (i in 1:m) {
  State = Step(State)
  XZY[i, ] = State
}

X = XZY[, "X"]
Z = XZY[, "Z"]
Y = XZY[, "Y"]
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
###########################################
#Generowanie macierzy przejśćia
#Znalezienie rozkładu stacjonarnego [układ zależmy,pominąć jedno równanie i zastąpić je sumowaniem do jedynki]
#
#https://mimuw.edu.pl/~wniem/Sym_Stoch/SyStoMC.pdf str 145
#
