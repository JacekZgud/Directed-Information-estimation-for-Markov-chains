setClassUnion("arrayORmatrix", c("matrix", "array"))
require(data.table)
source("./R/helpers.r")
#' Definition of a MarkovProcess class.
#' @exportClass MarkovProcess

setClass(
  "MarkovProcess",
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
    simulation = 'vector',
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
    trans_matrix =  data.frame(),
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
  ),
  validity = check_input
)

#' Initialize MarkovProcess class
#'
#' @param n A number of states.
#' @param d A number of vertices.
#' @param parent_struct A (d x d) matrix defining parental relations where rows represent parents.
#' @param nod_names A vector of names for vertices.
#' @param trans_probs A transition probabilities list.
#' @returns MarkovProcess class
#' @export
#' @examples
#' MarkovProcess(2,2,node_names=c("X","Y"))

MarkovProcess <- function(n = 1,
                          d = 1,
                          parent_struct = NULL,
                          node_names = c("X"),
                          trans_probs = NULL) {
  prob_col = paste("prob", c(0:(n - 1)), sep = "_")
  if (is.null(parent_struct)) {
    parent_struct <- matrix(nrow = d,
                            ncol = d,
                            data = 0)
    diag(parent_struct) <- 1
    rownames(parent_struct) <-
      colnames(parent_struct) <-  node_names

  }
  if (is.null(trans_probs)) {
    trans_probs <- Trans_prob(node_names, parent_struct, n, d)
  }
  return(
    new(
      "MarkovProcess",
      dim_num = n,
      node_num = d,
      node_names = c(node_names),
      parent_struct = parent_struct,
      trans_prob = trans_probs,
      prob_cols = c(prob_col)
    )
  )
}







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
  TransitionProbabilities <- vector("list", c(d))
  names(TransitionProbabilities) <-  Nodenames
  prob.vector.names <- paste("prob", c(0:(n - 1)), sep = '_')
  parents.number <- rowSums(parent_struct)
  for (node in Nodenames) {
    parents <- Nodenames[which(parent_struct[node,] == 1)]
    n0_parents <- parents.number[parents]
    temp <- expand.grid(rep(list(0:c(n - 1)), length(n0_parents)))
    colnames(temp) <- parents

    if (parents.number[[node]] != 0) {
      prob_matr <- matrix(
        sample.int(10 ^ 3, n * nrow(temp)),
        nrow = nrow(temp),
        ncol = n,
        byrow = TRUE,
        dimnames = list(c(), prob.vector.names)
      )
      prob_matr <- t(apply(prob_matr, 1, function(x)
        x / sum(x)))
      TransitionProbabilities[[node]] = data.table::data.table(cbind(temp, prob_matr))
      data.table::setkeyv(TransitionProbabilities[[node]], parents)
    }
    else{
      prob_matr <- matrix(
        runif(n),
        nrow = 1,
        ncol = n,
        byrow = TRUE,
        dimnames = list(c(), prob.vector.names)
      )
      prob_matr <- t(apply(prob_matr, 1, function(x)
        x / sum(x)))
      TransitionProbabilities[[node]] <-
        data.table::data.table(prob_matr)

    }
  }
  return(TransitionProbabilities)
}
