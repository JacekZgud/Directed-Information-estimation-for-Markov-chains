# here, all functions for simulation may be found i.e.
# - markov_sim()
# - Step()
#


#waterfall for simulation
#markov process simulation

markov_sim = function(obj, m) {
  "
  Calculate simulation of markov process, given relevant markov class object and length of simulation
  --------------------
  Arguments:
    obj - markov_process object
    m - length of process to be simulated
  Returns:
    Simulation of markov process of lenght m.
  "
  timer = Sys.time()
  XZY = matrix(nrow = m, ncol = obj@node_num)
  colnames(XZY) = obj@node_names
  State = rep(0, obj@node_num)
  for (i in 1:m) {
    State = Step(State, obj@trans_prob, obj@node_names, obj@parent_struct)
    XZY[i,] = State
    print_progress(i, m, timer)
  }
  return(XZY)
}


Step = function(State,
                TransitionProbabilities,
                names,
                ParentStructure) {
  "
  Calculate step for markov simulation, that is simulate process for time t given t-1 state.
  -----------------
  Arguments:
    State - vector of (t-1) state
    TransitionProbabilities - 3 dimensional array of conditional probabilities
    names - names of vertices
    ParentStructure - (d x d) matrix defining parental relations where rows represent parents
  Returns:
    simulated realization of t time process
  "
  NewState = numeric()
  prob_columns = paste("prob", 0:(n - 1), sep = "_")
  for (vertex in names) {
    Parents = (ParentStructure[vertex, ] == 1)
    ParentState = State[Parents]
    Prob = TransitionProbabilities[[vertex]]
    Prob = Prob[as.list(ParentState), ..prob_columns]
    NewState[vertex] =  sample.int(n, 1, prob = Prob) - 1
  }
  return(NewState)
}