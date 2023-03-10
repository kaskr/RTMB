# RTMB - R bindings for TMB

## Motivation

The package `RTMB` provides a native R interface for *a subset of* `TMB` so you can avoid coding in C++.

`RTMB` only affects the `TMB` function `MakeADFun` that builds the objective function. Once `MakeADFun` has been invoked, everything else is *exactly the same* and *models run as fast* as if coded in C++.

PROS:

- Fast to change and re-run models because no compilation is needed.
- Debugging can be performed using the normal R debugger rather than gdb.

CONS:

- Not all features are supported. See [current list of working examples](./tmb_examples).

## Install

NOTE: Requires at least `TMB-1.9.2`.

```r
remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB")
```

## Known issues

- If you get a segfault while installing the package, please re-install `Rcpp` from source and try again (see https://github.com/kaskr/RTMB/issues/5).
