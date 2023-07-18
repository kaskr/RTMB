# RTMB - R bindings for TMB

## Description

The package `RTMB` provides a native R interface for *a substantial subset of* `TMB` so you can avoid coding in C++.

`RTMB` only affects the `TMB` function `MakeADFun` that builds the objective function. Once `MakeADFun` has been invoked, everything else is *exactly the same* and *models run as fast* as if coded in C++.

PROS:

- Fast to change and re-run models because no compilation is needed.
- Debugging can be performed using the normal R debugger rather than gdb.
- Most common TMB features are supported. See [current list of working examples](./tmb_examples).

More information, including vignettes and documentation, can be found on the [RTMB universe page](https://kaskr.r-universe.dev/RTMB).

## Install

### Binary (Windows/Mac without compiler)

Run from R:

```r
install.packages('RTMB', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
```

**A parallel version is available as a separate package 'RTMBp'**

### Source

NOTE: Requires at least `TMB-1.9.3`.

```r
remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB")
```

## Known issues

- If you get a segfault while installing the package, please re-install `Rcpp` from source and try again (see https://github.com/kaskr/RTMB/issues/5).
