# define class  markov_process
setClassUnion("arrayORmatrix", c("matrix", "array"))


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
    trans_matrix_list = "data.frame"
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
    trans_matrix_list = data.frame()
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
  if (is.null(trans_probs)) {
    trans_probs = Trans_prob(NodeNames, parent_struct, n, d)
  }
  if (is.null(parent_struct)) {
    parent_struct = matrix(nrow = d,
                           ncol = d,
                           data = 0)
    diag(parent_struct) = 1
    rownames(parent_struct) = colnames(parent_struct) = NodeNames
    
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

Trans_prob = function(Nodenames, parent_struct, n, d) {
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
  set.seed(1)
  TransitionProbabilities = vector("list", c(d))
  names(TransitionProbabilities) = Nodenames
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  NoParents = rowSums(parent_struct)
  for (node in Nodenames) {
    parents = Nodenames[which(parent_struct[node, ] == 1)]
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