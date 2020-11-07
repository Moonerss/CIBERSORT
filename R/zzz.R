##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("Welcome to 'CIBERSORT' package!
=======================================================================
", "You are using ", pkgname, " version ", version,
"=======================================================================")
  packageStartupMessage(msg)
}
