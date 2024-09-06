# define class  markov_process
setClassUnion("arrayORmatrix", c("matrix", "array"))
require(data.table)
require(dplyr)
source("./structs/stationary_distribution.R")
source("./structs/markov_simulation_src.R")
source("./structs/simulation_marginalized.R")

setClass(
  "markov_process",
  slots = c(
    dim_num = "numeric",
    node_num = "numeric",
    node_names = "vector",
    prob_cols = "vector",
    parent_struct = "arrayORmatrix",
    trans_prob = "list",
    trans_matrix = "ANY",
    trans_matrix_list = "data.frame",
    marg_sim = 'list',
    marg_sim.ft = "arrayORmatrix" ,
    marg_sim.target_name = "vector" ,
    marg_sim.sim_target = "arrayORmatrix" ,
    statio_prob = 'data.frame',
    simulation = 'ANY',
    trans_entropy_table = 'data.frame',
    trans_entropy = 'numeric'
  )
  ,
  prototype = c(
    dim_num = NA_integer_,
    node_num = NA_integer_,
    node_names = c(NA_character_),
    prob_cols = c(NA_character_),
    parent_struct = matrix(),
    trans_prob = list(),
    trans_matrix = matrix(),
    trans_matrix_list = data.frame(),
    marg_sim = list(
      ft = array(),
      target_name = c(NA_character_),
      sim_target = array()
    ),
    statio_prob = data.frame(),
    simulation = c(),
    trans_entropy_table = data.frame(),
    trans_entropy = NA_real_
  )
)

markov_process_init <- function(n = 1,
                                d = 1,
                                parent_struct = NULL,
                                nod_names = c("X"),
                                trans_probs = NULL) {
  "
  Define markov process class
  -------------------------
  Argments:
    n - number of states
    d - number of vertices
    parent_struct - (d x d) matrix defining parental relations where rows represent parents
    nod_names - vector of names for vertices
    trans_probs - transition probabilities
  Returns:
    markov process class
  "
  NodeNames = nod_names
  prob_col = paste("prob", c(0:(n - 1)), sep = "_")
  if (is.null(parent_struct)) {
    parent_struct = matrix(nrow = d,
                           ncol = d,
                           data = 0)
    diag(parent_struct) = 1
    rownames(parent_struct) = colnames(parent_struct) = NodeNames
    
  }
  if (is.null(trans_probs)) {
    trans_probs = Trans_prob(NodeNames, parent_struct, n, d)
  }
  return(
    new(
      "markov_process",
      dim_num = n,
      node_num = d,
      node_names = c(NodeNames),
      parent_struct = parent_struct,
      trans_prob = trans_probs,
      prob_cols = c(prob_col)
    )
  )
}






# define structure for markov chain

Trans_prob <- function(Nodenames, parent_struct, n, d) {
  "
  Assign conditional probability for each vertice
  -----------------------
  Arguments:
    Nodenames - names of vertices
    ParentStructure - (d x d) matrix defining parental relations where ones in row represent parents
    n - number of states
    d - number of vertices
  Returns:
    3 dimensional array of conditional probabilities
  "
  TransitionProbabilities = vector("list", c(d))
  names(TransitionProbabilities) = Nodenames
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  NoParents = rowSums(parent_struct)
  for (node in Nodenames) {
    parents = Nodenames[which(parent_struct[node, ] == 1)]
    n0_parents = NoParents[parents]
    temp = expand.grid(rep(list(0:c(n - 1)), length(n0_parents)))
    colnames(temp) = parents
    
    if (NoParents[[node]] != 0) {
      prob_matr = matrix(
        sample.int(10 ^ 3, n * nrow(temp)),
        nrow = nrow(temp),
        ncol = n,
        byrow = TRUE
      )
      colnames(prob_matr) = prob_vector_names
      prob_matr = t(apply(prob_matr, 1, function(x)
        x / sum(x)))
      TransitionProbabilities[[node]] = data.table(cbind(temp, prob_matr))
      setkeyv(TransitionProbabilities[[node]], parents)
    }
    else{
      prob_matr = matrix(runif(n),
                         nrow = 1,
                         ncol = n,
                         byrow = TRUE)
      colnames(prob_matr) = prob_vector_names
      prob_matr = t(apply(prob_matr, 1, function(x)
        x / sum(x)))
      TransitionProbabilities[[node]] = data.table(prob_matr)
      
    }
  }
  return(TransitionProbabilities)
}


marginalized_runner <-
  function(obj,
           target = c(obj@node_names[1]),
           n = 1000,
           printer = FALSE) {
    if (nrow(obj@trans_matrix_list) == 0) {
      obj@trans_matrix_list = trans_matrix(obj, list_form = TRUE)
    }
    if (nrow(obj@statio_prob) == 0) {
      obj@statio_prob = stationary_probability(obj)
    }
    
    obj@marg_sim = markov_sim_Y(obj, n, target)
    if (printer) {
      for (i in target) {
        message(table(obj@marg_sim$sim_target[, i]) / n)
        message(data.table(obj@statio_prob)[, sum(statio_prob), by = eval(i)])
      }
    }
    return(obj)
  }

simulation_runner <- function(obj, m = 1000, printer = FALSE) {
  if (nrow(obj@statio_prob) == 0)
    obj@statio_prob = stationary_probability(obj)
  
  obj@simulation = markov_sim(obj, m)
  message('DONE')
  if (printer) {
    for (i in obj@node_names) {
      message(table(obj@simulation[, i]) / m)
      message(data.table(obj@statio_prob)[, sum(statio_prob), by = eval(i)])
    }
  }
  
  return(obj)
}
