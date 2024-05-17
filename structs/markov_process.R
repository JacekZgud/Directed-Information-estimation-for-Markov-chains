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
    trans_matrix_list = "data.frame",
    marg_sim = 'list',
    statio_prob = 'data.frame',
    simulation = 'ANY'
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
    marg_sim = vector('list'),
    statio_prob = data.frame(),
    simulation = c()
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
    parents = Nodenames[which(parent_struct[node,] == 1)]
    n0_parents = NoParents[parents]
    temp = expand.grid(rep(list(0:c(n - 1)), length(n0_parents)))
    colnames(temp) = parents
    
    if(NoParents[[node]] !=0){
      prob_matr = matrix(
        sample.int(10^3,n * nrow(temp)),
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
      prob_matr = matrix(
        runif(n),
        nrow = 1,
        ncol = n,
        byrow = TRUE
      )
      colnames(prob_matr) = prob_vector_names
      prob_matr = t(apply(prob_matr, 1, function(x)
        x / sum(x)))
      TransitionProbabilities[[node]] = data.table( prob_matr)
      
    }
  }
  return(TransitionProbabilities)
}


marginalized_runner <- function(obj, target = c(obj@node_names[1]), n_2) {
  if (nrow(obj@trans_matrix_list)==0){
    obj@trans_matrix_list = trans_matrix(obj, list_form = TRUE)
  }
  if (nrow(obj@statio_prob)==0){
    obj@statio_prob = stationary_probability(obj)
  }
  
  obj@marg_sim = markov_sim_Y(obj, n_2, target)
  cat('\n','DONE')
  
  for (i in target) {
    print(table(obj@marg_sim$sim_target[, i])/n_2)
    print(data.table(obj@statio_prob)[, sum(statio_prob), by = eval(i)])
  }
  return(obj)
}

simulation_runner <- function(obj, m) {
  if (nrow(obj@statio_prob)==0)
    obj@statio_prob = stationary_probability(obj)
  
  obj@simulation = markov_sim(obj, m)
  cat('\n','DONE','\n')
  
  for (i in obj@node_names) {
    print(table(obj@simulation[, i])/m)
    print(data.table(obj@statio_prob)[, sum(statio_prob), by = eval(i)])
  }
  return(obj)
}
