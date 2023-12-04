# One step of Markov Chain
#
Step = function(State) {
  NewState = c()
  prob_columns = paste("prob", c(1:(n - 1)), sep = "_")
  for (vertex in NodeNames) {
    Parents = (ParentStructure[vertex, ] == 1)
    ParentState = State[Parents]
    #print(as.list(ParentState))
    #loc =  apply(TransitionProbabilities[[vertex]][NodeNames[Parents]], 1, function(x) all(c(x) == ParentState ))
    Prob = TransitionProbabilities[[vertex]][as.list(ParentState), ..prob_columns]
    #Prob = TransitionProbabilities[[vertex]][loc,prob_columns]
    NewState[vertex] = sum(runif(1) > cumsum(as.vector(Prob)))
  }
  return(NewState)
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


#calculate transition probability for given configuration

prob_transition = function(x) {
  vec = c()
  for (node in NodeNames) {
    parents = NodeNames[which(ParentStructure[node, ] == 1)]
    Parents = paste(parents, "(t-1)", sep = "")
    
    ParentState = x[Parents]
    node_state = as.double(x[paste(node, "(t)", sep = "")])
    prob_col = paste("prob", node_state, sep = "_")
    
    # loc =  apply(TransitionProbabilities[[node]][parents], 1, function(x) all(c(x) == c(ParentState) ))
    # prob = TransitionProbabilities[[node]][loc,paste("prob",node_state,sep="_")]
    
    prob = setDT(TransitionProbabilities[[node]])[as.list(ParentState), ..prob_col]
    #print(prob)
    vec[node] = as.double(prob)
  }
  #print(as.double(prod(vec)))
  
  return(as.double(prod(vec)))
  
}

