#' Simulate Markov filtering procedure.
#' @param obj A MarkovProcess object.
#' @param target A vector of node names to simulate
#' @param sim.length Length of simulation to perform.
#' @return A MarkovProcess object.
#' @export

simulate.marginalized <-
  function(obj,
           target = c(obj@node_names[1]),
           sim.length = 1000) {
    if (nrow(obj@trans_matrix_list) == 0) {
      obj <- trans.matrix(obj)
    }
    if (nrow(obj@statio_prob) == 0) {
      obj <- stationary.probability(obj)
    }

    obj@marg_sim <- markov.sim.filt(obj, sim.length, target)
    return(obj)
  }
