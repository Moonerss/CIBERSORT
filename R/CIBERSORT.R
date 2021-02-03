#' Main functions
#'
#' The Main function of CIBERSORT
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @import utils
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#'   mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#'   results <- cibersort(sig_matrix, mixture_file)
#'   ## example 2
#'   data(LM22)
#'   data(mixed_expr)
#'   results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
#' }
cibersort <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE){

  #read in data
  if (is.character(sig_matrix)) {
    X <- read.delim(sig_matrix, header=T, sep="\t", row.names=1, check.names = F)
    X <- data.matrix(X)
  } else {
    X <- sig_matrix
  }

  if (is.character(mixture_file)) {
    Y <- read.delim(mixture_file, header=T, sep="\t", row.names=1, check.names = F)
    Y <- data.matrix(Y)
  } else {
    Y <- mixed_expr
  }


  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  }

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while (itor <= mixtures) {

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1

  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}
