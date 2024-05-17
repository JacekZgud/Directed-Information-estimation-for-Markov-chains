# here all functions used to calculate information may be found
# - entropy
# - mutual info
# - transfer entropy

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
                         n_2 = 1000) {
  "
    obj - 'markov_process' class object
    "
  if ((length(obj@marg_sim) == 0) |
      (length(setdiff(obj@marg_sim$target_name, target) > 0) |
       length(setdiff(target, obj@marg_sim$target_name) > 0))) {
    obj =  marginalized_runner(obj, target, n_2)
  }
  target = obj@marg_sim$target_name
  target_parent_nodes = vector('list')
  nodes_not_parent = vector('list')
  nodes_without_target_vector = setdiff(obj@node_names, target)
  for (i in target) {
    target_parent_nodes[[i]] = names(which(obj@parent_struct[i, ] == 1))
    
    nodes_not_parent[[i]] = setdiff(obj@node_names, target_parent_nodes[[i]])
  }
  origin_prob = data.table(obj@statio_prob)[, sum(statio_prob), keyby =
                                              eval(c(obj@node_names))]
  prob_cols = obj@prob_cols
  
  entropies_Py = rep(0, obj@dim_num ^ obj@node_num)
  Py = vector('list')
  for (i in target) {
    Py[[i]] = obj@trans_prob[[i]]
    if (!(length(nodes_not_parent[[i]]) == 0)) {
      conf = expand.grid(rep(list(0:c(
        obj@dim_num - 1
      )), length(nodes_not_parent[[i]])))
      colnames(conf) = nodes_not_parent[[i]]
      Py[[i]] = data.table(merge(conf, Py[[i]]))
    }
    ordering = c(obj@node_names, prob_cols)
    Py[[i]] = Py[[i]][, ..ordering]
    setkeyv(Py[[i]], obj@node_names)
    entropies_Py = entropies_Py + apply(data.frame(Py[[i]])[, prob_cols], 1, entropy) # by conditional independence of target
  }
  column_names = c(paste(target, c(
    rep("(t)", length(target)), rep("(t-1)", length(target))
  ), sep = ""),
  paste(nodes_without_target_vector, "(t-1)", sep = ""))
  
  P_target = setDT(obj@trans_matrix_list)[, sum(prob), by = column_names]
  colnames(P_target) = c(column_names, 'prob')
  setkeyv(P_target, paste(target, rep("(t-1)", length(target)), sep = ""))
  
  target_entropy = sum(entropies_Py * origin_prob$V1)
  cat('Entropy given all past states:', target_entropy, "\n")
  cat(
    'Entropy given all past states:',
    sum(P_target[, entropy(prob), keyby = c(paste(obj@node_names, '(t-1)', sep =
                                                    ""))]$V1 * origin_prob$V1),
    "(calculated differently)\n"
  )
  
  mt = data.table(obj@marg_sim$mt)
  ft = data.table(obj@marg_sim$ft)
  ys = obj@marg_sim$sim_target
  colnames(mt) = c(target, c(1:nrow(obj@marg_sim$sim_target)))
  setkeyv(mt, target)
  
  cat('Calculating trajectory probabilities...\n')
  
  end = nrow(ys)
  time = Sys.time()
  
  get_value = function(i) {
    y = as.vector(ys[i, ])
    j = i + length(target)
    print_progress(i, end, time)
    mt[as.list(c(y)), ..j]
  }
  saver = Vectorize(get_value)(c(1:nrow(ys)))
  cat('\nDONE\n')
  saver = as.vector(unlist(saver))

  
  entropy_only_target = apply(mt[,-..target], 2, entropy)
  #print(entropy_only_target)
  
  combinations = function(x,n){
    probs = na.omit(Reduce(`*`, shift(x, 0:(n), type="lead"),right = FALSE))
    c(sum(entropy_only_target[-c(1:n)]*probs)/sum(probs),sum(probs))
  }
  len = ceiling(0.6*length(entropy_only_target))
  
  if(len > 250){
    len=250
  }
  
  smth = Vectorize(combinations,'n')(x=saver,n=c(1:(len)))
  ind = which(smth[1,] >0 )
  info = sum(smth[1,ind]*smth[2,ind])/sum(smth[2,ind])
  print(info)
  
  
  hist(entropy_only_target)
  abline(v = target_entropy, col = 'red')
  abline(v =  info, col = 'orange')
  legend(
    "bottomright",
    legend = c(
      "Entropy | X_V(t-1))",
      'Entropy | X_Y(t...)'
    ),
    pch = "|",
    col = c("red", "orange")
  )
  attr(obj, 'trans_entropy') = info - target_entropy
  
  cat(
    '\nTarget:',
    target,
    '\n',
    'Other nodes: ',
    nodes_without_target_vector,
    '\n',
    'Transfer_entropy:',
    obj@trans_entropy,
    '\n'
  )
  return(obj)
}
n_2 = 100
process = marginalized_runner(process, c('X'), 10000)
process = trans_entropy(process, c('X'), n_2 = 10000)
plot(c(1:(length(process@marg_sim$mt[1,-2]))),process@marg_sim$mt[1,-2])
# dla mt liczymy końcowo uśrednionie przez prawdopodobienstwo trajektorii
# rozpisać na współrzędne wzorek o liczeniu informacji wzajemnej pod warunkiem cąłości historii
# rozpisać jak liczę entropie targetu pod warunkiem historii targetu
