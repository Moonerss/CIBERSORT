#' do permutations
#'
#' Do the permutations analysis
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @importFrom purrr reduce map
#' @export

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  itorect <- function(Ylist, X) {
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    return(mix_r)
  }

  if (perm == 1) {
    dist <- itorect(Ylist = Ylist, X = X)
  } else {
    dist <- purrr::map(1:perm, ~ itorect(Ylist = Ylist, X = X)) %>%
      purrr::reduce(rbind)
  }

  newList <- list("dist" = dist)
}
