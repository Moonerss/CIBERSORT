#' open multiple threads process
#' @param nThreads The number of threads
#' @param maxSize The max memory size in global, default 500MB, the unit is MB
#' @param verbose output other useful information
#' @import future
#' @importFrom future availableCores plan
#' @export
#' @examples
#' \dontrun{
#'   enableParallel()
#' }
enableParallel <- function(nThreads = NULL, maxSize = 500, verbose = FALSE) {
  nCores <- future::availableCores()
  options(future.globals.maxSize = maxSize*1024^2)
  if (verbose) message("The maxSize is ", maxSize, " Mb.")
  if (is.null(nThreads)) {
    if (nCores < 4) {
      nThreads <- nCores
    } else {
      nThreads <- nCores - 2
    }
  }
  if (!is.numeric(nThreads) || nThreads < 2)
    stop("nThreads must be numeric and at least 2.")
  if (nThreads > nCores) {
    if (verbose) {
      message("Requested number of threads is higher than number of available processors (or cores)")
      message(paste("The max working processes is", nCores))
    }
  }
  if (verbose) message(paste("Allowing parallel execution with up to", nThreads, "working processes."))
  future::plan("multicore", workers = nThreads)
  invisible(nThreads)
}
