# here all functions for simulation with marginalized variables may be found
# - markov_sim_y
# - stepY
source("./structs/helpers.r")
source("./structs/stationary_distribution.R")

markov_sim_Y <- function(obj,
                         n,
                         target = "Y",
                         condition_set = NULL) {
  "
  Simulate process with marginalized variables
  ----------------------
  Arguments:
    obj - markov_process object
    n - length of simulation
    target - name of variable to be simulated
    condition_set - vector of nodes to condition by, NULL if no conditioning is applied
  "
  numb = obj@dim_num
  target_parent_nodes = names(which(obj@parent_struct[target, ] == 1))
  parent_without_target = target_parent_nodes[target_parent_nodes != target]
  nodes_without_target = obj@node_names[obj@node_names != target]
  nodes_not_parent = setdiff(nodes_without_target, parent_without_target)
  prob_cols = obj@prob_cols
  if (sum(obj@node_names == target) < 1) {
    print("Target out of scope")
    return(NULL)
  }
  if (is.null(condition_set)) {
    #######################################
    # preparation of relevant distributions
    
    Py = obj@trans_prob[[target]]
    conf = expand.grid(rep(list(0:c(
      process@dim_num - 1
    )), length(nodes_not_parent)))
    colnames(conf) = nodes_not_parent
    Py = data.table(merge(conf, Py))
    Py = relocate(Py, target, .before = 1)
    setkeyv(Py, c(target, nodes_without_target))
    
    if (nrow(obj@trans_matrix_list) == 0)
      P = setDT(trans_matrix(obj, list_form = TRUE))
    else
      P = setDT(obj@trans_matrix_list)
    
    setkeyv(P, c(
      paste(target, c("(t)", "(t-1)"), sep = ""),
      paste(nodes_without_target, "(t-1)", sep = "")
    ))
    
    #######################################
    # preparation for simulation vector
    fty = vector("list", c(3))
    names(fty) = c("y", "ft", "mt")
    
    fty$y = 0 # y(t)
    
    configs = expand.grid(rep(list(0:c(
      process@dim_num - 1
    )), length(nodes_without_target)))
    colnames(configs) = nodes_without_target
    values = matrix(1 / nrow(configs),
                    ncol = 1,
                    nrow = row(configs))
    colnames(values) = "prob"
    
    fty$ft = data.table(cbind(configs, values)) # P(X(t) |Y(t)=y,Y(<t))
    setkeyv(fty$ft, nodes_without_target)
    fty$ft = data.frame(fty$ft)
    
    fty$mt = c(1:obj@dim_num)    #P(Y(t+1) |Y(t)=y,Y(<t))
    
    ########################################
    # preparation for memory vectors
    Ys = rep(NULL, n) # history of target
    Fts = vector("list", n) # estimates of fts
    Mts = matrix(NA, ncol = numb, nrow = n) # estimates of mt
    colnames(Mts) = prob_cols
    
    #setup for timer
    timer = Sys.time()
    
    for (t in 1:n) {
      fty = stepY(fty, configs, Py, target, P, prob_cols)
      Ys[t] = fty$y
      Fts[[t]] = fty$ft
      Mts[t,] = fty$mt
      print_progress(t, n, timer)
      
    }
  }
  else{
    # Not Implemented Yet
    print("Conditioning not implemented yet")
    fty = NULL
  }
  
  out = fty
  out$y = Ys
  out$ft = Fts
  out$mt = Mts
  return(out)
}



stepY <- function(fr, configs, Py, target, P, prob_cols) {
  "
  Simulate step of markov process with marginalized variables without conditioning
  --------------------
  Arguments:
    fr - simulation vector consisting of y,ft,mt
      - y: Y(t) simulated
      - ft: P(X(t) |Y(t)=y,Y(<t)): probs of configurations of X
      - mt: P(Y(t+1)|Y(t)=y,Y(<t)): vector of length n
  Returns:
    fr - simulation vector
  "
  y = fr$y
  ft = fr$ft # probability of P(X(t) |Y(t)=y,Y(<t)): probs of configurations of X
  mt = fr$mt # P(Y(t+1)|Y(t)=y,Y(<t)): vector of length n
  for (i in 1:n) {
    colmn = prob_cols[i]
    
    mt[i] = sum(as.vector(unlist(Py[as.list(c(y)), ..colmn])) * as.vector(unlist(ft["prob"])))
  }
  Y = sample.int(n, 1, prob = mt) - 1
  
  ft_1 = P[as.list(c(Y, y)), sum(prob * ft$prob), by = eval(paste(colnames(configs), "(t)", sep =
                                                                    ""))]
  #print(ft_1)
  colnames(ft_1) = c(colnames(configs), "prob")
  ft$prob = ft_1$prob / sum(ft_1$prob)
  
  fr$y = Y
  fr$ft = data.frame(ft)
  fr$mt = mt
  return(fr)
}


step_Y_cond <- function(fr) {
  "
  Simulate step of markov process with marginalized variables with conditioning
  --------------------
  Arguments:
    fr - simulation vector
  Returns:
    fr - simulation vector
  "
  print("Not done yet")
  return(fr)
}
