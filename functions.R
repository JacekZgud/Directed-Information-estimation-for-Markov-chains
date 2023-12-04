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
    Prob = setDT(TransitionProbabilities[[vertex]])[as.list(ParentState), ..prob_columns]
    #Prob = TransitionProbabilities[[vertex]][loc,prob_columns]
    NewState[vertex] = sum(runif(1) > cumsum(as.vector(Prob)))
  }
  return(NewState)
}
# Switch between binary expansion and usual (decimal) representation of integers
#
# Binary expansion of number x (sequence of l digits)
Dec2Bin = function(x, l)
{
  cyfra = 1
  b = rep(0, l)
  
  while (x > 1) {
    r = x %% 2
    x = x %/% 2
    b[cyfra] = r
    cyfra = cyfra + 1
  }
  b[cyfra] = x
  b = rev(b)
  return(b)
}

# Convert a sequence b of 0/1 to a number
Bin2Dec <- function(b)
{
  l = length(b)
  b = rev(b)
  x = sum(b * 2 ^ (1:l - 1))
  return(x)
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

