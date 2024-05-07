# here are all functions for simulation with marginalized variables may be found
# - markov_sim_y
# - stepY
require('dplyr')
source("./structs/helpers.r")
source("./structs/stationary_distribution.R")

markov_sim_Y <- function(obj,
                         n,
                         target = c("Y"),
                         condition_set = NULL) {
  "
  Simulate obj with marginalized variables
  ----------------------
  Arguments:
    obj - markov_obj object
    n - length of simulation
    target - vector of variable names to be simulated
    condition_set - vector of nodes to condition by, NULL if no conditioning is applied
  "
  cat('Preparation for marginalized simulation...')
  numb = obj@dim_num
  
  if (length(intersect(obj@node_names, target)) < 1) {
    print("Target out of scope")
    return(NULL)
  }
  target_parent_nodes = vector('list')
  parent_without_target = vector('list')
  nodes_without_target = vector('list')
  nodes_not_parent = vector('list')
  nodes_without_target_vector = setdiff(obj@node_names, target)
  for (i in target) {
    target_parent_nodes[[i]] = names(which(obj@parent_struct[i, ] == 1))
    parent_without_target[[i]] = setdiff(target_parent_nodes[[i]], i)
    nodes_without_target[[i]] = setdiff(obj@node_names, i)
    nodes_not_parent[[i]] = setdiff(obj@node_names, target_parent_nodes[[i]])
  }
  prob_cols = obj@prob_cols
  if (is.null(condition_set)) {
    #######################################
    # preparation of relevant distributions
    
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
      ordering = c(target, nodes_without_target_vector,prob_cols)
      Py[[i]] = Py[[i]][,..ordering]
      setkeyv(Py[[i]], c(target, nodes_without_target_vector))
    }
    
    if (nrow(obj@trans_matrix_list) == 0){
      P = setDT(trans_matrix(obj, list_form = TRUE))
    }
    else{
      P = setDT(obj@trans_matrix_list)
    }
    
    setkeyv(P, c(
      paste(target, c(
        rep("(t)", length(target)), rep("(t-1)", length(target))
      ), sep = ""),
      paste(nodes_without_target_vector, "(t-1)", sep = "")
    ))
    
    
    #######################################
    # preparation for simulation vector
    fty = vector("list", c(3))
    names(fty) = c("target", "ft", "mt")
    
    #y
    fty$target = rep(0, length(target))  # y(t)
    names(fty$target) = target
    #ft
    configs = expand.grid(rep(list(0:c(
      obj@dim_num - 1
    )), length(nodes_without_target_vector)))
    colnames(configs) = nodes_without_target_vector
    values = matrix(runif(nrow(configs)),
                    ncol = 1,
                    nrow = nrow(configs))
    values[, 1] = values[, 1] / sum(values)
    colnames(values) = "prob"
    
    fty$ft = data.table(cbind(configs, values)) # P(X(t) |Y(t)=y,Y(<t))
    setkeyv(fty$ft, nodes_without_target_vector)
    fty$ft = data.frame(fty$ft)
    #mt
    fty$mt = expand.grid(rep(list(0:c(
      obj@dim_num - 1
    )), length(target)))
    colnames(fty$mt) = target    #P(Y(t+1) |Y(t)=y,Y(<t))
    values = matrix(0,
                    ncol = 1,
                    nrow = nrow(fty$mt))
    colnames(values) = "prob"
    fty$mt['prob'] = values
    
    ########################################
    # preparation for memory vectors
    Ys = matrix(NA, ncol = length(target), nrow = n) # history of target
    colnames(Ys) = target
    
    Fts = matrix(NA, ncol = n, nrow = nrow(fty$ft)) # estimates of fts
    Fts = cbind(fty$ft[, nodes_without_target_vector], Fts)
    colnames(Fts) = c(nodes_without_target_vector , as.character(1:n))
    
    Mts = matrix(NA, ncol = n, nrow = nrow(fty$mt)) # estimates of mt
    colnames(Mts) = as.character(1:n)
    Mts = cbind(fty$mt[, target], Mts)
    
    #setup for timer
    timer = Sys.time()
    cat(' DONE','\n')
    
    stepY <- function(fr, Py, P, configs, target, prob_cols) {
      "
      Simulate step of markov obj with marginalized variables without conditioning
      --------------------
      Arguments:
        fr - simulation vector consisting of y,ft,mt
          - y: Y(t) simulated
          - ft: P(X(t) |Y(t)=y,Y(<t)): probs of configurations of X
          - mt: P(Y(t+1)|Y(t)=y,Y(<t)): vector of length n
      Returns:
        fr - simulation vector
      "
      y = fr$target
      ft = fr$ft # probability of P(X(t) |Y(t)=y,Y(<t)): probs of configurations of X
      mt = fr$mt[target] # P(Y(t+1)|Y(t)=y,Y(<t)): vector of length n
      Y = y
      
      work_mt = mt
      for (var in target) {
        #setkeyv(ft,)

        temp = data.frame(work = c(0:(length(prob_cols) - 1)),
                          prob = as.vector(apply(Py[[var]][as.list(y),
                                                           ..prob_cols] * ft$prob, 2, sum)))
        colnames(temp) = c(var, 'prob')
        Y[var] = sample.int(length(prob_cols), 1, prob = temp[['prob']]) - 1
        work_mt = setDT(work_mt)[temp, on = eval(var)]
      }
      mt[, 'prob'] = apply(work_mt[, !..target], 1, prod)
      mt = data.frame(mt)
      ft_1 = P[as.list(c(Y, y)), sum(prob * ft$prob), by = eval(paste(configs, "(t)", sep =
                                                                        ""))]
      #print(ft_1)
      colnames(ft_1) = c(configs, "prob")
      ft$prob = ft_1$prob / sum(ft_1$prob)
      
      
      fr$target = Y
      fr$ft = ft
      fr$mt = mt
      return(fr)
    }
    
    #simulation loop
    for (t in 1:n) {
      fty = stepY(fty,
                  Py,
                  P,
                  nodes_without_target_vector ,
                  target,
                  prob_cols)
      Ys[t, ] = fty$target
      Fts[, as.character(t)] = fty$ft$prob
      Mts[, as.character(t)] = fty$mt$prob
      print_progress(t, n, timer)
      
    }
    cat('\n','DONE')
    
  }
  else{
    # Not Implemented Yet
    step_Y_cond <- function(fr) {
      "
      Simulate step of markov obj with marginalized variables with conditioning
      --------------------
      Arguments:
        fr - simulation vector
      Returns:
        fr - simulation vector
      "
      print("Not done yet")
      return(fr)
    }
    
    print("Conditioning not implemented yet")
    fty = NULL
  }
  out = vector('list')
  out$sim_target = Ys
  out$ft = Fts
  out$mt = Mts
  out$target_name = target
  return(out)
}


