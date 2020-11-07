#' open multiple threads process
#' @param nThreads The number of threads
#' @param maxSize The max memory size in global, default 500MB, the unit is MB
#' @import future
#' @export
#' @examples
#' \dotrun{
#'   enableParallel()
#' }
enableParallel <- function(nThreads = NULL, maxSize = 500) {
  nCores <- future::availableCores()
  options(future.globals.maxSize = maxSize*1024^2)
  message("The maxSize is ", maxSize, " Mb.")
  if (is.null(nThreads)) {
    if (nCores < 4)
      nThreads = nCores
    else nThreads = nCores - 3
  }
  if (!is.numeric(nThreads) || nThreads < 2)
    stop("nThreads must be numeric and at least 2.")
  if (nThreads > nCores) {
    message("Requested number of threads is higher than number of available processors (or cores)")
    message(paste("The max working processes is", nCores))
  }
  message(paste("Allowing parallel execution with up to", nThreads, "working processes."))
  future::plan("multicore", workers = nThreads)
  invisible(nThreads)
}
