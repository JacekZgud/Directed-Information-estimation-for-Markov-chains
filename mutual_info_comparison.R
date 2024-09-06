# example
library(dplyr)
library(data.table)
source("./structs/markov_process.R")
source("./structs/stationary_distribution.R")
source("./structs/markov_simulation_src.R")
source("./structs/simulation_marginalized.R")
source('./structs/information_calculation.R')

# parametrise the code
node_dim = 2
nodes = 2
work_names = c("X", "Y")


#define parent structure
ParentStructure = matrix(nrow = nodes, ncol = nodes, data = 0)
rownames(ParentStructure) = colnames(ParentStructure) = work_names
ParentStructure[2, 1] = 1
ParentStructure[1, 1] = 1

prob = seq(from = 1 / 2,
           to = 0.99999,
           length.out = 50)


# define function calculating transfer_entropy estimator

info_estimator = function(a, par_struct, n_2 = 1000) {
  proc = markov_process_init(node_dim, nodes, par_struct, work_names)
  proc@trans_prob$X$prob_1 = c(1 - a, a)
  proc@trans_prob$X$prob_0 = 1 - proc@trans_prob$X$prob_1
  
  proc@trans_prob$Y$prob_1 = c(1 - a, a)
  proc@trans_prob$Y$prob_0 = 1 - proc@trans_prob$Y$prob_1

  proc = (trans_entropy(proc, c('Y'), n = n_2))
  
  return(proc@trans_entropy)
}

info_estimator(1 / 2, ParentStructure)
infos = c()
n = 10
for (i in c(1:n)) {
  infos = infos + Vectorize(info_estimator, 'a')(a = prob, ParentStructure)
  print_progress(i,n)
}
infos = infos / n


# plot the results

h = function(b) {
  1 - b * log(b / (1 / 2), base = 2) - (1 - b) * log((1 - b) / (1 / 2), base =
                                                       2)
}

mutual_info = function(x) {
  a = x
  b = x
  c = a * (b ^ 2) + a * ((1 - b) ^ 2) + (1 - a) * (b ^ 2) + 2 * (1 - a) *
    b * (1 - b)
  c = 1 - ((1 - a) * b ^ 2 + (1 - a) * (1 - b) ^ 2 + 2 * a * b * (1 - b))
  h(c) - h(x)
}

dir_pol = function(x) {
  1 - h(x)
}


plot(
  prob,
  infos,
  ylim = c(0, 1),
  type = 'line',
  col = 'green',
  xlab = 'a (a=b)',
  ylab = 'Information'
)
legend(
  "topleft",
  legend = c("Transfer entropy",
             'M.Information|(t-1)',
             'Information Flow'),
  pch = "|",
  col = c("green", "red", 'blue')
)
curve(
  dir_pol,
  from = 0.5,
  to = 1,
  n = 1000,
  add = TRUE,
  col = 'blue'
)
curve(
  mutual_info,
  from = 0.5,
  to = 1.00000000,
  n = 10000,
  add = TRUE,
  col = 'red'
)

# check if the results are correct
Vectorize(mutual_info)(prob) > infos


#########################################################################################################
dkl = function(P, Q) {
  P = as.vector(P)
  Q = as.vector(Q)
  sum(P * log(P / Q, base = 2))
}

dkl_comp = function(x, inv = FALSE) {
  x = as.numeric(x)
  if (inv) {
    P = c(1 / 2, 1 / 2)
    Q = c(x, 1 - x)
  }
  else{
    Q = c(1 / 2, 1 / 2)
    P = c(x, 1 - x)
  }
  sum(P * log(P / Q, base = 2))
}

dkl_comp_fixed_inv_vec <- Vectorize(function(x) {
  dkl_comp(x, inv = TRUE)
})
dkl_comp_vec <- Vectorize(function(x) {
  dkl_comp(x, inv = FALSE)
})

# Plot the curve using the modified vectorized function
xes = c(0.001, 0.999)
curve(
  dkl_comp_fixed_inv_vec,
  from = xes[1],
  to = xes[2],
  xlab = 'p',
  ylab = '',
  col = 'red'
)
curve(
  dkl_comp_vec,
  from = xes[1],
  to = xes[2],
  add = TRUE,
  col = 'blue'
)
legend(
  "topleft",
  legend = c("D_Q", "D_P"),
  pch = "|",
  col = c("red", 'blue')
)
