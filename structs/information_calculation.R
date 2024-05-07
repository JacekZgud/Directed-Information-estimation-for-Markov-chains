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
  origin = setdiff(obj@node_names, target)
  origin_prob = data.table(obj@statio_prob)[, sum(statio_prob), keyby =
                                              eval(c(obj@node_names))]
  
  target_parent_nodes = vector('list')
  nodes_not_parent = vector('list')
  nodes_without_target_vector = setdiff(obj@node_names, target)
  for (i in target) {
    target_parent_nodes[[i]] = names(which(obj@parent_struct[i, ] == 1))
    
    nodes_not_parent[[i]] = setdiff(obj@node_names, target_parent_nodes[[i]])
  }
  prob_cols = obj@prob_cols
  
  Py = vector('list')
  entropies_Py = rep(0, obj@dim_num ^ obj@node_num)
  for (i in target) {
    Py[[i]] = obj@trans_prob[[i]]
    if (!(length(nodes_not_parent[[i]]) == 0)) {
      conf = expand.grid(rep(list(0:c(
        obj@dim_num - 1
      )), length(nodes_not_parent[[i]])))
      colnames(conf) = nodes_not_parent[[i]]
      Py[[i]] = data.table(merge(conf, Py[[i]]))
    }
    ordering = c(obj@node_names,prob_cols)
    Py[[i]] = Py[[i]][,..ordering]
    setkeyv(Py[[i]], obj@node_names)
    entropies_Py = entropies_Py + apply(data.frame(Py[[i]])[, prob_cols], 1, entropy) # by conditional independence of target
  }
  
  print(Py)
  print(origin_prob)
  target_entropy = sum(entropies_Py * origin_prob$V1)
  print(target_entropy)  
  
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
  #time = Sys.time()
  #for (i in 1:n) {
  #  y = as.vector(ys[i, ])
  #  j = i + length(target)
  #  saver[i] = mt[as.list(c(y)), ..j]
  #  print_progress(i, n, time)
  #}
  cat('\nDONE\n')
  probs = cumprod(as.vector(saver))
  print(probs[c(1:10)])
  entropy_only_target = apply(mt[,-..target], 2, entropy)
  #print(entropy_only_target)
  subset = c(100:length(entropy_only_target))
  attr(obj, 'trans_entropy') = sum(entropy_only_target[subset] * probs[subset]) / sum(probs[subset]) - target_entropy
  cat(
    '\nTarget:',
    target,
    '\n',
    'Other nodes: ',
    nodes_without_target_vector,
    '\n',
    'Transfer_entropy:',
    sum(entropy_only_target * probs)/sum(probs) - target_entropy,
    '\n'
  )
  return(obj)
}
n_2=10000
process = trans_entropy(process, c('X'), n_2 = 10000)
#process@trans_prob
# dla mt liczymy końcowo uśrednionie przez prawdopodobienstwo trajektorii
# rozpisać na współrzędne wzorek o liczeniu informacji wzajemnej pod warunkiem cąłości historii
# rozpisać jak liczę entropie targetu pod warunkiem historii targetu
