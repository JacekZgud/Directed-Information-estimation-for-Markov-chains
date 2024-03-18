# This file contains all functions needed to calculate stationary distribution i.e. 
# - trans_matrix()
# - prob_transition()
# - stationary_probability()


trans_matrix = function(obj) {
  "
  Calculates transition matrix for markov process
  ----------
  Arguments:
    obj - markov_process object
  Returns:
    transition matrix
  "
  n = obj@dim_num
  d = obj@node_num
  names = obj@node_names
  Trans = expand.grid(rep(list(0:c(n - 1)), 2 * d))
  colnames(Trans) = c(paste(obj@node_names, "(t)", sep = ""),
                      paste(obj@node_names, "(t-1)", sep = ""))
  Trans['prob'] = apply(Trans, 1, function(x)
    prob_transition(x, names, obj@parent_struct, obj@trans_prob))
  return(matrix(
    as.vector(Trans['prob'])$prob,
    ncol = n ^ d,
    nrow = n ^ d,
    byrow = TRUE
  ))
}
#calculate transition probability for given configuration

prob_transition = function(x,
                           names,
                           ParentStructure,
                           TransitionProbabilities) {
  "
  Cacluate probability of given configuration at time (t) given configuration at (t-1)
  ---------------
  Arguments:
  x - vector consisting of stacked single configurations at time t and t-1
  names - Node names
  ParentStructure - (d x d) matrix defining parental relations where rows represent parents
  TransitionProbabilities - 3 dimensional array of conditional probabilities
  "
  vec = c()
  for (node in names) {
    parents = names[which(ParentStructure[node,] == 1)]
    Parents = paste(parents, "(t-1)", sep = "")
    
    ParentState = x[Parents]
    node_state = as.double(x[paste(node, "(t)", sep = "")])
    prob_col = paste("prob", node_state, sep = "_")
    prob = setDT(TransitionProbabilities[[node]])[as.list(ParentState), ..prob_col]
    vec[node] = as.double(prob)
  }
  
  return(as.double(prod(vec)))
  
}

stationary_probability = function(obj) {
  "
  Calculate stationary distribution for markov_process object
  -----------------
  Argument:
    obj - markov_process object
  Returns:
    stationary distribution
  "
  n = process@dim_num
  d = process@node_num
  A = t(process@trans_matrix - diag(ncol = n ^ d, nrow = n ^ d))
  A = rbind(A, rep(1, n ^ d))
  b = c(rep(0, n ^ d), 1)
  res_statio = qr.solve(A, b)
  res_statio = cbind(expand.grid(rep(list(0:c(
    n - 1
  )), d)), res_statio)
  colnames(res_statio) = c(obj@node_names, 'statio_prob')
  return(res_statio)
}