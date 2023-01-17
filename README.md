# RTMB - R bindings for TMB

## Motivation

It would be nice if we could provide a native R interface for TMB, i.e allow something like:

```r
f <- function (data, parameters) {
  a <- parameters$a
  b <- parameters$b
  sd <- exp(parameters$logsd)
  ADREPORT(sd)
  x <- data$x
  y <- data$y
  mu <- a * x + b
  -sum(dnorm(y, mu, sd, log=TRUE))
}
obj <- MakeADFun(f, random="b")
```



## Design considerations

For this to work we need to represent a vector of AD variables natively on the R side. There are at least two options:

- Represent an AD vector on the R side as 'externalptr' to the corresponding C++ structure i.e., something like `structure(new("externalptr"), class="ad")`. This gives the 'illusion' that we are working natively in R while actually the objects live in C++. The file `example.R` shows this idea.
- Represent an AD scalar on the R side by an atomic R type wide enough to encode it. In our case, and AD scalar is 128 bit which fits inside the 'RComplex' type. We can then reuse the native R vector machinery using a representation like `structure(complex(), class="ad")`. The file `complex.R` shows the idea.

The example implementations demonstrate that with a fairly small effort (~100 lines of R code) we can enable a surprisingly rich set of R features. For example, both implementations would allow the following snippet:

```r
x <- advector(1:6)
a <- array(x, c(3,2))
t(a)
z <- outer(x, x, "+")
s <- apply(z, 1, sum)
z <- sin(s)
```

## Caveats

### Sub assignment

R has special mechanisms to avoid copying unless strictly necessary. For example consider the loop:

```r
x <- advector(1:1e4)
for (i in 1:length(x)) x[i] <- x[i]+1
```

Under standard R semantics `x` is modified inplace because no other copy of `x` exists. Our 'complex' implementation gets this mechanism for free because it reuses R's vector machinery. In contrast the 'externalptr' implementation must deep copy the entire `x` vector on each modification!

### Losing the class attribute

There's always a potential risk of losing the class attribute. For example `base::diff.default` has

```
r <- unclass(x)
... ## Calculate and modify r
class(r) <- oldClass(x)
r
```

This construct is catastrophic for the 'complex' implementation: Assume x is an AD type. When unclassing 'x' we get a normal complex number which will be modified using standard complex arithmetic. The result no longer have valid encoding so when the class is restored to AD we get an illegal object. Segfault is unavoidable.

In contrast, the 'externalptr' implementation will error out soon after unclassing as no meaningful arithmetic exists for external pointers.

In both cases the problem can easily be solved using an overload, i.e.

```r
diff.advector <- function (x) head(x,-1) - tail(x,-1) 
```

### Sparse matrices

The Matrix package uses S4 classes to represent objects and C code for most sparse matrix manipulatins. This design is not easy to take advantage of by either of our AD implementations.

## Conclusion so far

- 'complex' implementation has notable benefits over the 'externalptr' implementation. However, it is more risky approach.
 
