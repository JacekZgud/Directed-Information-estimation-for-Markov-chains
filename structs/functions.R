#define class for markov process
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
  NodeNames = tail(LETTERS,d)
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
  TransitionProbabilities = vector("list", c(d))
  names(TransitionProbabilities) = Nodenames
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  for (node in Nodenames) {
    parents = Nodenames[which(ParentStructure[node,] == 1)]
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


# One step of Markov Chain
#
Step = function(State, TransitionProbabilities,names,ParentStructure) {
  NewState = numeric()
  prob_columns = paste("prob", 1:(n - 1), sep = "_")
  for (vertex in names) {
    Parents = (ParentStructure[vertex, ] == 1)
    ParentState = State[Parents]
    Prob = TransitionProbabilities[[vertex]]
    Prob = Prob[as.list(ParentState), ..prob_columns]
    NewState[vertex] = sum(runif(1) > cumsum(as.vector(Prob)))
  }
  return(NewState)
}

#waterfall for simulation
#markov process simulation

markov_sim = function(obj,m) {
  XZY = matrix(nrow = m, ncol = obj@node_num)
  colnames(XZY) = obj@node_names
  State = rep(0, obj@node_num)
  for (i in 1:m) {
    State = Step(State, obj@trans_prob,obj@node_names,obj@parent_struct)
    XZY[i,] = State
    print_progress(i,m)
  }
  return(XZY)
}


#calculate transition matrix
trans_matrix = function(obj) {
  n=obj@dim_num
  d=obj@node_num
  names = obj@node_names
  Trans = expand.grid(rep(list(0:c(n - 1)), 2 * d))
  colnames(Trans) = c(paste(obj@node_names, "(t)", sep = ""),
                      paste(obj@node_names, "(t-1)", sep = ""))
  Trans['prob'] = apply(Trans, 1, function(x)
    prob_transition(x,names,obj@parent_struct,obj@trans_prob))
  Trans2 = matrix(
    as.vector(Trans['prob'])$prob,
    ncol = n ^ d,
    nrow = n ^ d,
    byrow = TRUE
  )
  return(Trans2)
}
#calculate transition probability for given configuration

prob_transition = function(x,names,ParentStructure,TransitionProbabilities) {
  vec = c()
  for (node in names) {
    parents = names[which(ParentStructure[node, ] == 1)]
    Parents = paste(parents, "(t-1)", sep = "")
    
    ParentState = x[Parents]
    node_state = as.double(x[paste(node, "(t)", sep = "")])
    prob_col = paste("prob", node_state, sep = "_")
    prob = setDT(TransitionProbabilities[[node]])[as.list(ParentState), ..prob_col]
    vec[node] = as.double(prob)
  }
  
  return(as.double(prod(vec)))
  
}

#solve stationary equation to obtain stationary probability distribution
stationary_probability = function(obj) {
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

#visuals
print_progress <- function(iteration, total) {
  percent_complete <- round((iteration / total) * 100 / 5)
  cat(
    "\r[",
    paste(rep("=", percent_complete), collapse = ""),
    paste(rep(" ", 20 - percent_complete), collapse = ""),
    "] ",
    percent_complete * 5,
    "%",
    " (",
    iteration,
    "/",
    total,
    ")",
    sep = ""
  )
  flush.console()
}

