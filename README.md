README: `rpom` package
================

# Installation

The `rpom` package can be installed with the following script:

``` r
#install.packages("remotes")
remotes::install_github("aterui/rpom")
```

# Overview

The R package `rpom` provides functions to yield analytical/numerical
solutions for food chain length in branching river networks.

# Example

`fcl()` returns analytical solution of food chain length with specified
values of ecosystem size (`size`) and branching rate (`lambda`). This
function does not allow any top down effects.

``` r
## generate a food web with preferential prey model
fwb <- mcbrnet::ppm(n_species = 10, n_basal = 3, l = 10, theta = 1)

## fcl
rpom::fcl(fwb, lambda = 0.1, size = 10)
```

    ## [1] 3
    ## attr(,"tp")
    ##  [1] 1 1 1 2 3 3 2 2 3 3
    ## attr(,"p_hat")
    ##  [1] 0.9478343 0.9478343 0.9478343 0.9440926 0.9438083 0.9438083 0.9440926
    ##  [8] 0.9720463 0.9438083 0.9458792

`nfcl()` returns numerical solution of food chain length. This function
allows top down effects.

``` r
## fcl
rpom::nfcl(fwb, lambda = 0.1, size = 10)
```

    ## [1] 2
    ## attr(,"tp")
    ## [1] 1 1 1 2
    ## attr(,"p_hat")
    ##         1         2         3         4         5         6         7         8 
    ## 0.4320875 0.4783428 0.4320875 0.0000000 0.0000000 0.0000000 0.0000000 0.2923892 
    ##         9        10 
    ## 0.0000000 0.0000000
