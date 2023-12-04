# define structure for markov chain

Trans_prob = function(Nodenames, ParentStructure, n, d) {
  TransitionProbabilities = vector("list", c(d))
  names(TransitionProbabilities) = NodeNames
  prob_vector_names = paste("prob", c(0:(n - 1)), sep = '_')
  for (node in NodeNames) {
    parents = NodeNames[which(ParentStructure[node, ] == 1)]
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
Step = function(State,TransitionProbabilities) {
  NewState = c()
  prob_columns = paste("prob", c(1:(n - 1)), sep = "_")
  for (vertex in NodeNames) {
    Parents = (ParentStructure[vertex,] == 1)
    ParentState = State[Parents]
    Prob = TransitionProbabilities[[vertex]][as.list(ParentState), ..prob_columns]
    NewState[vertex] = sum(runif(1) > cumsum(as.vector(Prob)))
  }
  return(NewState)
}

#waterfall for simulation
#markov process symulation

markov_sim = function(structure,m,d,n){
  XZY = matrix(nrow = m, ncol = d)
  colnames(XZY) = NodeNames
  State = rep(0, d)
  for (i in 1:m) {
    State = Step(State,structure)
    XZY[i, ] = State
  }
  return(XZY)
}

#calculate transition matrix
trans_matrix = function(n, d) {
  Trans = expand.grid(rep(list(0:c(n - 1)), 2 * d))
  colnames(Trans) = c(paste(NodeNames, "(t)", sep = ""),
                      paste(NodeNames, "(t-1)", sep = ""))
  Trans['prob'] = apply(Trans, 1, function(x)
    prob_transition(x))
  Trans2 = matrix(
    as.vector(Trans['prob'])$prob,
    ncol = n ^ d,
    nrow = n ^ d,
    byrow = TRUE
  )
  return(Trans2)
}
#solve stationary equation to obtain stationary probability distribution
stationary_probability = function(transition_matrix,n,d){
  A = t(transition_matrix - diag(ncol = n ^ d, nrow = n ^ d))
  A = rbind(A, rep(1, n ^ d))
  b = c(rep(0, n ^ d), 1)
  res_statio = qr.solve(A, b)
  res_statio = cbind(expand.grid(rep(list(0:c(
    n - 1
  )), d)), res_statio)
  colnames(res_statio) = c(NodeNames, 'statio_prob')
  return(res_statio)
}


#calculate transition probability for given configuration

prob_transition = function(x) {
  vec = c()
  for (node in NodeNames) {
    parents = NodeNames[which(ParentStructure[node,] == 1)]
    Parents = paste(parents, "(t-1)", sep = "")
    
    ParentState = x[Parents]
    node_state = as.double(x[paste(node, "(t)", sep = "")])
    prob_col = paste("prob", node_state, sep = "_")
    
    # loc =  apply(TransitionProbabilities[[node]][parents], 1, function(x) all(c(x) == c(ParentState) ))
    # prob = TransitionProbabilities[[node]][loc,paste("prob",node_state,sep="_")]
    
    prob = setDT(TransitionProbabilities[[node]])[as.list(ParentState), ..prob_col]
    vec[node] = as.double(prob)
  }

  return(as.double(prod(vec)))
  
}
