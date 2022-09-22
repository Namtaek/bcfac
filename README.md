
# bcfac


The repository is almost a copy of Yoo Young-hoon's work, who graduated from Sungkyunkwan University with a master's degree in statistics. There is a problem with compiling Rcpp in a local environment, so it is temporarily used in public mode.

## Installation

install the development version from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("Namtaek/bcfac")
```

### For **Mac** Users Only

bartcs supports multi-threading by OpenMP. If you

- Have OpenMP installed
- And want to use OpenMP multi-threading

then install package from source with:

```r
# build from source
install.packages("bcfac", type = "source")

# count OpemMP thread
bartcs::count_omp_thread() # this should be greater than 1
```
