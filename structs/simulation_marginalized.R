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

  if (length(intersect(obj@node_names, target)) < 1) {
    print("Target out of scope")
    return(NULL)
  }

  nodes_without_target_vector = setdiff(obj@node_names, target)
  prob_cols = obj@prob_cols
  if (is.null(condition_set)) {
    #######################################
    # preparation of relevant distributions
    
    if (nrow(obj@trans_matrix_list) == 0) {
      P = setDT(trans_matrix(obj, list_form = TRUE))
    }
    else{
      P = setDT(obj@trans_matrix_list)
    }
    
    column_names = c(paste(target, c(
      rep("(t)", length(target)), rep("(t-1)", length(target))
    ), sep = ""),
    paste(nodes_without_target_vector, "(t-1)", sep = ""))
    P_target = P[, sum(prob), by = column_names]
    colnames(P_target) = c(column_names, 'prob')
    setkeyv(P_target, paste(target, rep("(t-1)", length(target)), sep =
                              ""))
    
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
    fty$target = sample.int(obj@dim_num,length(target),replace = TRUE)-1  # y(t)
    names(fty$target) = target
    #ft
    not_target_configurations = expand.grid(rep(list(0:c(obj@dim_num - 1)), length(nodes_without_target_vector)))
    colnames(not_target_configurations) = nodes_without_target_vector
    values = matrix(runif(nrow(not_target_configurations)),
                    ncol = 1,
                    nrow = nrow(not_target_configurations))
    values[, 1] = values[, 1] / sum(values)
    colnames(values) = "prob"
    
    fty$ft = data.table(cbind(not_target_configurations, values)) # P(X(t) |Y(t)=y,Y(<t))
    setkeyv(fty$ft, nodes_without_target_vector)
    fty$ft = data.frame(fty$ft)
    #mt
    fty$mt = expand.grid(rep(list(0:c(obj@dim_num - 1)), length(target)))
    colnames(fty$mt) = target    #P(Y(t+1) |Y(t)=y,Y(<t))
    values = matrix(0,
                    ncol = 1,
                    nrow = nrow(fty$mt))
    values[,1] = sample.int(obj@dim_num,nrow(fty$mt),replace = TRUE)-1
    colnames(values) = "prob"
    fty$mt['prob'] = values
    fty$mt = as.data.table(fty$mt)
    setkeyv(fty$mt, target)
    
    ########################################
    # preparation for memory vectors
    Ys = matrix(NA, ncol = length(target), nrow = n) # history of target
    colnames(Ys) = target
    
    Fts = matrix(NA, ncol = n, nrow = nrow(fty$ft)) # estimates of fts
    Fts = cbind(fty$ft[, nodes_without_target_vector], Fts)
    colnames(Fts) = c(nodes_without_target_vector , as.character(1:n))
    
    Mts = matrix(NA, ncol = n, nrow = nrow(fty$mt)) # estimates of mt
    colnames(Mts) = as.character(1:n)
    Mts = cbind(as.data.frame(fty$mt[, ..target],column_names=target), Mts)
    
    #setup for timer
    timer = Sys.time()
    cat(' DONE', '\n')
    
    
    #work
    target_cols = paste(target, '(t)', sep = "")
    stepY <- function(fr, Py, P, not_target_configurations, target, prob_cols) {
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
      mt = fr$mt[, ..target] # P(Y(t+1)|Y(t)=y,Y(<t)): vector of length n
      work = P_target[as.list(y), sum(prob * ft$prob), keyby = eval(target_cols)]
      index = sample.int(nrow(work), 1, prob = work$V1)
      mt[, 'prob'] = work$V1
      Y = unlist(as.vector(work[index, ..target_cols]))
      ft_1 = P[as.list(c(Y, y)), sum(prob * ft$prob), by = eval(paste(not_target_configurations, "(t)", sep =""))]
      colnames(ft_1) = c(not_target_configurations, "prob")
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
      Ys[t,] = fty$target
      Fts[, as.character(t)] = fty$ft$prob
      Mts[, as.character(t)] = fty$mt$prob
      print_progress(t, n, timer)
      
    }
    cat(' ', 'DONE','\n')
    
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
  out$sim_length = n
  return(out)
}
if (sys.nframe() == 0L) {
  process = marginalized_runner(process, c('Y','X'), 1000)
  
  
}



