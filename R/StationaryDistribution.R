require("data.table")
#' Calculates transition matrix for MarkovProcess class
#' @param obj A MarkovProcess object.
#' @returns A MarkovProcess object with some slots populated.
#' @export
#'
#'

trans.matrix <- function(obj) {
  message('Calculating transition matrix... ')
  n <- obj@dim_num
  d <- obj@node_num
  Trans <- expand.grid(rep(list(0:c(n - 1)), 2 * d))
  colnames(Trans) <- c(paste(obj@node_names, "(t)", sep = ""),
                       paste(obj@node_names, "(t-1)", sep = ""))
  Trans['prob'] <- apply(Trans, 1, function(x)
    prob_transition(x, obj@node_names, obj@parent_struct, obj@trans_prob))
  message('DONE', '\n')
  obj@trans_matrix_list <- Trans
  obj@trans_matrix <- matrix(
    as.vector(Trans['prob'])$prob,
    ncol = n ^ d,
    nrow = n ^ d,
    byrow = TRUE
  )
  return(obj)

}


#' Calculate stationary distribution for MarkovProcess object.
#' @param obj A MarkovProcess object.
#' @returns A MarkovProcess object with populated slots assosiated with stationary distribution.
#' @export



stationary.probability <- function(obj) {
  if (is.null(nrow(obj@trans_matrix))) {
    obj <- trans.matrix(obj)
  }

  n <- obj@dim_num
  d <- obj@node_num
  A <- t(obj@trans_matrix - diag(ncol = n ^ d, nrow = n ^ d))
  A <- rbind(A, rep(1, n ^ d))
  b <- c(rep(0, n ^ d), 1)

  tryCatch({
    res_statio <- qr.solve(A, b)
  },
  error = function(e) {
    print(e)
    stop()
  },

  warning = function(w) {
    print("There was a warning message.")
  },

  finally = {
    res_statio <- cbind(expand.grid(rep(list(0:c(
      n - 1
    )), d)), res_statio)
    colnames(res_statio) <- c(obj@node_names, 'statio_prob')
    obj@statio_prob = res_statio
    return(obj)
  })


}

#' Calculate probability of given configuration at time (t) given configuration at (t-1).
#' @param x A vector consisting of stacked single configurations at time (t) and (t-1).
#' @param names A vector containing node names.
#' @param ParentStructure A (d x d) matrix defining parental relations where rows represent parents.
#' @param TransitionProbabilities A 3 dimensional array of conditional probabilities.
#' @return Probability of transition.
#' @import data.table



prob_transition <- function(x,
                            names,
                            ParentStructure,
                            TransitionProbabilities) {
  vec <- c()
  for (node in names) {
    parents <- names[which(ParentStructure[node,] == 1)]
    parents <- paste(parents, "(t-1)", sep = "")
    ParentState <- as.list(x[parents])
    node_state <- as.double(x[paste(node, "(t)", sep = "")])
    prob_column_index <- c(length(parents) + node_state + 1)

    if (length(parents) > 0) {
      prob <-
        TransitionProbabilities[[node]][ParentState,]
      vec[node] <-
        as.double(unlist(prob)[prob_column_index])
    }
    else{
      vec[node] <-
        as.double(unlist(TransitionProbabilities[[node]])[prob_column_index])
    }
  }
  return(as.double(prod(vec)))

}
