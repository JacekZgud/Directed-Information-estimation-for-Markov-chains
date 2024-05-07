# Here all helpers functions may be found, 
# - print_progress

print_progress <- function(iteration, total, start_time=Sys.time()) {
  percent_complete_num <- floor((iteration / total) * 100)
  percent_complete <- round((iteration / total) * 100 / 5)
  cat(
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
    difftime(Sys.time(),start_time,units = "secs" ),
    sep = ""
  )
  flush.console()
}
