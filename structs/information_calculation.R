# here all functions used to calculate information may be found
# - entropy
# - mutual info
# - transfer entropy
require('data.table')
source('./structs/helpers.R')

entropy = function(vector, b = 2) {
  "
  Compute an entropy for vector of probabilities
  ----------------------------------
  Args:
    vector - vector of probabilities
  "
  #if (sum(vector) != 1){
  #  cat('Vector is not a probability vector')
  #}
  vector = vector[(vector != 0)]
  out = -sum(vector * log(vector, base = b))
  
  return(out)
}

trans_entropy = function(obj,
                         target = 'Y',
                         cond = NULL,
                         n = 1000) {
  "
    obj - 'markov_process' class object
    "
  
  # perform marginalized simulation for target variable or check if it has been performed earlier
  
  if ((length(obj@marg_sim) == 0) |
      (length(setdiff(obj@marg_sim$target_name, target) > 0) |
       length(setdiff(target, obj@marg_sim$target_name) > 0))) {
    obj =  marginalized_runner(obj, target, n)
  }
  
  ft = data.table(obj@marg_sim$ft)
  ys = obj@marg_sim$sim_target
  
  # prepare containers for relevant column names
  target = obj@marg_sim$target_name
  
  nodes_without_target_vector = setdiff(obj@node_names, target)
  
  prob_cols = obj@prob_cols
  
  column_names = c(paste(target, c(
    rep("(t)", length(target)), rep("(t-1)", length(target))
  ), sep = ""),
  paste(nodes_without_target_vector, "(t-1)", sep = ""))
  
  
  
  # marginal p(node_names)
  origin_prob = data.table(obj@statio_prob)[, sum(statio_prob), keyby =
                                              eval(c(obj@node_names))]
  
  # define P_target - p(target(t)|all_nodes(t-1))
  P_target = setDT(obj@trans_matrix_list)[, sum(prob), by = column_names]
  colnames(P_target) = c(column_names, 'prob')
  setkeyv(P_target, paste(target, rep("(t-1)", length(target)), sep = ""))
  
  # calculate entropy given all previous states of markov chain
  target_entropy = sum(P_target[, entropy(prob), keyby = c(paste(obj@node_names, '(t-1)', sep =
                                                                   ""))]$V1 * origin_prob$V1)
  
  # start estimation of entropy p(target(t)|target(<t))
  cat('Calculating entropies...\n')
  
  end = obj@marg_sim$sim_length
  time = Sys.time()
  
  # entropy estimator p(y(t)|y(<t))
  entropy_target_calc = function(index) {
    print_progress(index, end, time)
    entropy(P_target[as.list(ys[index, ]), sum(prob * ft[[as.character(index)]]), by =
                       eval(paste(target,
                                  rep("(t)", length(target)), sep = ""))]$V1)
  }
  entropy_only_target = Vectorize(entropy_target_calc)(c(1:obj@marg_sim$sim_length))
  
  info_niem = sum(entropy_only_target) / obj@marg_sim$sim_length
  
  
  
  # final transfer entropy value
  
  table = data.frame(
    target = target,
    origin = c(paste(
      nodes_without_target_vector, collapse = ' '
    )),
    `transfer entropy` = info_niem - target_entropy
  )
  
  # assign transfer entropy to relevant attributes of markov_process class
  attr(obj, 'trans_entropy') = info_niem - target_entropy
  if (nrow(obj@trans_entropy_table) > 0) {
    attr(obj, 'trans_entropy_table') = rbind(obj@trans_entropy_table, table)
  }
  else{
    attr(obj, 'trans_entropy_table') = table
  }
  
  
  
  cat(
    '\n---------------------------------------------------\n',
    'Target:',
    target,
    '\n',
    'Origin nodes: ',
    nodes_without_target_vector,
    '\n',
    'Transfer_entropy:',
    obj@trans_entropy,
    '\n',
    'Entropy given all past states:',
    target_entropy,
    "\n",
    'Entropy given only target :',
    info_niem,
    '\n'
    
  )
  return(obj)
}


if (sys.nframe() == 0L) {
  n_2 = 1000
  node = c('X')
  process = trans_entropy(process, node, n = n_2)
}
