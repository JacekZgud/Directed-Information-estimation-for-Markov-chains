#' Simulate Markov process
#' @param obj A MarkovProcess class.
#' @param sim.length Length of simulation.
#' @export

simulate <- function(obj, sim.length = 1000) {
  if (nrow(obj@statio_prob) == 0)
    obj <- stationary.probability(obj)

  obj@simulation <- markov.sim(obj, sim.length)
  return(obj)
}
