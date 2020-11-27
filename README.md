
# CIBERSORT

<!-- badges: start -->
<!-- badges: end -->

The goal of CIBERSORT is to run the CIBERSORT flow.

## Installation

You can install the released version of CIBERSORT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("devtools")
devtools::install_github("Moonerss/CIBERSORT")
```

## Example

``` r
library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, mixture_file)
```

