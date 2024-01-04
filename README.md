# RTMB - R bindings for TMB

<!-- badges: start -->
[![cran version](http://www.r-pkg.org/badges/version/RTMB)](https://cran.r-project.org/package=RTMB)
<!-- badges: end -->

## Description

The package `RTMB` provides a native R interface for *a substantial subset of* `TMB` so you can avoid coding in C++.

`RTMB` only affects the `TMB` function `MakeADFun` that builds the objective function. Once `MakeADFun` has been invoked, everything else is *exactly the same* and *models run as fast* as if coded in C++.

PROS:

- Fast to change and re-run models because no compilation is needed.
- Debugging can be performed using the normal R debugger rather than gdb.
- Most common TMB features are supported. See [current list of working examples](./tmb_examples).
- Simplified interface can automatically simulate from the model and calculate OSA residuals.

More information, including vignettes ([introduction](https://kaskr.r-universe.dev/articles/RTMB/RTMB-introduction.html) / [advanced](https://kaskr.r-universe.dev/articles/RTMB/RTMB-advanced.html)) and [documentation](https://kaskr.r-universe.dev/RTMB#reference), can be found on the [RTMB universe page](https://kaskr.r-universe.dev/RTMB).

## Install

### CRAN version

```r
install.packages('RTMB')
```

### Development version (Rtools/compiler *not* needed)

```r
install.packages('RTMB', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
```

### Experimental parallel version available as a separate package 'RTMBp'

```r
install.packages('RTMBp', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
```

### Source

NOTE: Requires at least `TMB-1.9.7`.

```r
remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB")
```

## Test the package

```r
tinytest::test_package("RTMB")
```

runs the [tests](./RTMB/inst/tinytest).

## Known issues

- If you get a segfault while installing the package, please re-install `Rcpp` from source and try again (see https://github.com/kaskr/RTMB/issues/5).
