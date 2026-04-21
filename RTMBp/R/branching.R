## Like a scalar 'ifelse'
branch <- function(test, yes, no = NULL) {
  old <- TapeConfig()
  TapeConfig(comparison="tape")
  on.exit(TapeConfig(old))
  ## Fall back to non AD ?
  if (length(test) != 1) stop("'test' must have length one")
  if (!inherits(test, "advector") || !getVariables(test)) {
    return (if (test) yes else no)
  }
  ## Make tape of both choices
  F <- MakeTape(function(x) c(yes, no), numeric())
  F$simplify()
  ## Let all references to outer variables become function arguments
  vars <- resolve_refs(.pointer(environment(F)$mod))
  ## Split and simplify the choices
  F1 <- MakeTape(function(x)F(x)[[1]], F$par())
  F2 <- MakeTape(function(x)F(x)[[2]], F$par())
  F1$simplify()
  F2$simplify()
  ## Prepend test condition to function argument
  vars <- c(test, vars)
  fp <- .pointer(environment(F1)$mod)
  gp <- .pointer(environment(F2)$mod)
  ## Call atomic branch operator
  Branch(fp, gp, vars)
}
