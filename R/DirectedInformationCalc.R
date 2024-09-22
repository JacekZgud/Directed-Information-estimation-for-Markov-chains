#' Calculate directed information between the nodes.
#'
#' @param obj A MarkovProcess class.
#' @param target A node name of information target.
#' @param cond A list of nodes to condition on (not implemented yet).
#' @param sim.length Length of simulation to be performed.
#' @export

trans_entropy = function(obj,
                         target = 'Y',
                         cond = NULL,
                         sim.length = 1000) {

  # perform marginalized simulation for target variable or check if it has been performed earlier
  if ((length(obj@marg_sim) == 0) |
      !all(target %in% obj@marg_sim$target_name)) {
    obj <-  simulate.marginalized(obj, target, sim.length)
  }

  ft <- data.table(obj@marg_sim$ft)
  ys <- obj@marg_sim$sim_target

  # prepare containers for relevant column names
  target <- obj@marg_sim$target_name

  nodes_without_target_vector <- setdiff(obj@node_names, target)

  prob_cols <- obj@prob_cols

  column_names <- c(paste(target, c(
    rep("(t)", length(target)), rep("(t-1)", length(target))
  ), sep = ""),
  paste(nodes_without_target_vector, "(t-1)", sep = ""))



  # marginal p(node_names)
  origin_prob <- data.table(obj@statio_prob)[, sum(statio_prob), keyby =
                                              eval(c(obj@node_names))]

  # define P_target - p(target(t)|all_nodes(t-1))
  P_target <- setDT(obj@trans_matrix_list)[, sum(prob), by = column_names]
  colnames(P_target) <- c(column_names, 'prob')
  data.table::setkeyv(P_target, paste(target, rep("(t-1)", length(target)), sep = ""))

  # calculate entropy given all previous states of markov chain
  target_entropy <- sum(P_target[, entropy(prob), keyby = c(paste(obj@node_names, '(t-1)', sep =
                                                                   ""))]$V1 * origin_prob$V1)

  # start estimation of entropy p(target(t)|target(<t))
  message('Calculating entropies...')

  end = obj@marg_sim$sim_length
  time = Sys.time()

  # entropy estimator p(y(t)|y(<t))
  entropy_target_calc <- function(index) {
    print_progress(index, end, time)
    entropy(P_target[as.list(ys[index-1, ]), sum(prob * ft[[as.character(index)]]), by =
                       eval(paste(target,
                                  rep("(t)", length(target)), sep = ""))]$V1)
  }
  entropy_only_target <- Vectorize(entropy_target_calc)(c(1:obj@marg_sim$sim_length))

  info_niem <- sum(entropy_only_target) / obj@marg_sim$sim_length

  message("DONE")

  # final transfer entropy value

  table <- data.frame(
    target = toString(target),
    origin = c(paste(
      nodes_without_target_vector, collapse = ' '
    )),
    `transfer entropy` = info_niem - target_entropy
  )

  # assign transfer entropy to relevant attributes of markov_process class
  obj@trans_entropy <- info_niem - target_entropy
  if (nrow(obj@trans_entropy_table) > 0) {
    obj@trans_entropy_table <- rbind(obj@trans_entropy_table, table)
  }
  else{
    obj@trans_entropy_table <- table
  }
  return(obj)
}

"
  message(
    '\n###################################################\n',
    '\n              Transfer entropy results             \n',
    '\n---------------------------------------------------\n',
    'Target:',
    obj@marg_sim$target_name,
    '\n',
    'Origin nodes: ',
    setdiff(obj@node_names, obj@marg_sim$target_name),
    '\n',
    'Transfer_entropy:',
    obj@trans_entropy
  )
"
