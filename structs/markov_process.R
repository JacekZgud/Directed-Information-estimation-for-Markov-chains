#define class  markov_process

source("./structs/helpers.r")

setClass(
  "markov_process",
  slots = c(
    dim_num = "numeric",
    node_num = "numeric",
    node_names = "vector",
    parent_struct = "matrix",
    trans_prob = "list"
  )
)


proc_init <- function(n, d, parent_struct) {
  "
  Define markov process class
  -------------------------
  Argments:
    n - number of states
    d - number of vertices
    parent_struct - (d x d) matrix defining parental relations where rows represent parents
  Returns:
    markov process class
  "
  NodeNames = tail(LETTERS, d)
  TransitionProbabilities = Trans_prob(NodeNames, ParentStructure, n, d)
  return(
    new(
      "markov_process",
      dim_num = n,
      node_num = d,
      node_names = c(NodeNames),
      parent_struct = parent_struct,
      trans_prob = TransitionProbabilities
    )
  )
  
}

# define structure for markov chain

Trans_prob = function(Nodenames, ParentStructure, n, d) {
  "
  Assign conditional probability for each vertice
  -----------------------
  Arguments:
    Nodenames - names of vertices
    ParentStructure - (d x d) matrix defining parental relations where rows represent parents
    n - number of states
    d - number of vertices
  Returns:
    3 dimensional array of conditional probabilities
  "
  set.seed(1)
  TransitionProbabilities = vector("list", c(d))
  names(TransitionProbabilities) = Nodenames
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  for (node in Nodenames) {
    parents = Nodenames[which(ParentStructure[node, ] == 1)]
    n0_parents = NoParents[parents]
    temp = expand.grid(rep(list(0:c(n - 1)), length(n0_parents)))
    colnames(temp) = parents
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
    setkeyv(TransitionProbabilities[[node]], parents)
  }
  return(TransitionProbabilities)
}