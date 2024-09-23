# Here all helpers functions may be found,
# - print_progress

print_progress <-
  function(iteration, total, start_time = Sys.time()) {
    percent_complete_num <- floor((iteration / total) * 100)
    percent_complete <- round((iteration / total) * 100 / 5)
    message(
      "\r[",
      paste(rep("=", percent_complete), collapse = ""),
      paste(rep(" ", 20 - percent_complete), collapse = ""),
      "] ",
      percent_complete_num,
      "%",
      " (",
      iteration,
      "/",
      total,
      "), ",
      round(difftime(Sys.time(), start_time, units = "secs"),2),
      appendLF = FALSE
    )
    flush.console()
  }


entropy <- function(vector, b = 2) {
  "
  Compute an entropy for vector of probabilities
  ----------------------------------
  Args:
    vector - vector of probabilities
  "
  vector <- vector[(vector != 0)]
  out <- -sum(vector * log(vector, base = b))

  return(out)
}

check_input <- function(object) {
  errors <- character()
  if (length(object@node_names) != object@node_num) {
    msg = paste(
      "Number of nodes",
      (object@node_names) ,
      " is different from length of node names provided",
      (object@node_num)
    )
    errors = c(errors, msg)
  }
  if ((nrow(object@parent_struct) != object@node_num)|(ncol(object@parent_struct) != object@node_num)) {
    msg = paste(
      "Parent structure matrix has wrong dimension ",
      "expected: (",
      object@node_num,
      "," ,
      object@node_num,
      "), got:",
      dim(object@parent_struct)
    )
    errors = c(errors, msg)
  }

  if (length(errors) == 0)
    TRUE
  else
    errors
}
