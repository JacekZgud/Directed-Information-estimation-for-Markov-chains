# here all functions used to calculate information may be found
# - entropy
# - mutual info
# - transfer entropy 




entropy = function(vector,b=2){
  "
  Compute an entropy for vector of probabilities
  ----------------------------------
  Args:
    vector - vector of probabilities
  "
  return(-sum(vector*log(vector,base=b)))
}
entropy(out$mt[1,])

trans_entropy = function(obj,target='Y',cond = NULL){
  "
  x-mtlike vector
  "
  parents = setdiff(colnames(obj@trans_prob[[target]]), obj@prob_cols)
  print(parents)
  parents_prob = data.table(process@statio_prob)[,sum(statio_prob),by=eval(parents)]
  
  setkeyv(parents_prob,parents)
  print(parents_prob)
  print(obj@trans_prob[[target]])
  mt = obj@marg_sim$mt
  sum(apply(mt,1,entropy)) - sum(apply(data.frame(obj@trans_prob[[target]])[,obj@prob_cols],1, entropy)*parents_prob)
  
}
trans_entropy(process)


# dla mt liczymy końcowo uśrednionie przez prawdopodobienstwo trajektorii 
