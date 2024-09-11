#' Core algorithm
#'
#' The core algorithm of CIBERSORT which is used svm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @param maxSize maximum size for the computation, to be passed to the future.global.maxSize. Default to 500 (in MB) 
#' @importFrom furrr future_map
#' @importFrom future availableCores
#'@importFrom stats cor
#' @import e1071
#' @export
CoreAlg <- function(X, y, maxSize=500){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  enableParallel(maxSize=maxSize)

  if (Sys.info()['sysname'] == 'Windows') {
    out <- future_map(1:svn_itor, res)
  } else {
    if (svn_itor <= availableCores() - 2) {
      enableParallel(nThreads = svn_itor, maxSize=maxSize)
    } else {
      enableParallel(maxSize=maxSize)
    }
    out <- future_map(1:svn_itor, res)
  }

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w <- weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}
