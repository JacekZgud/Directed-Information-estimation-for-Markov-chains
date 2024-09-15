#'  Calculate simulation of markov process, given relevant MarkovProcess class object.
#'  @param obj A MarkovClass object.
#'  @param sim.length A Length of process to be simulated.
#'  @returns  Simulation of markov process of length m.
#'  @export
#'

markov.sim <- function(obj, sim.length) {
  timer <- Sys.time()
  sim.trajectory <-
    matrix(
      nrow = sim.length,
      ncol = obj@node_num,
      dimnames = list(c(), obj@node_names)
    )

  State <- rep(0, obj@node_num)
  prob_columns <- paste("prob", 0:(obj@dim_num - 1), sep = "_")

  for (step in 1:sim.length) {
    State <-
      sim.step(obj, State)
    sim.trajectory[step, ] <- State
    print_progress(step, sim.length, timer)
  }
  return(sim.trajectory)
}


sim.step <- function(obj, State) {
  NewState <- numeric()
  for (vertex in obj@node_names) {
    Parents <- (obj@parent_struct[vertex,] == 1)
    ParentState <- as.list(State[Parents])
    if (as.double(sum(Parents)) > 0) {
      Prob <- unlist(obj@trans_prob[[vertex]][ParentState, ])[obj@prob_cols]
    }
    else{
      Prob <- unlist(obj@trans_prob[[vertex]])[obj@prob_cols]
    }
    NewState[vertex] <- sample.int(obj@dim_num, 1, prob = Prob) - 1
  }
  return(NewState)
}
