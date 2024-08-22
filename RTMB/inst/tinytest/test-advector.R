## Test indirect advector methods

library(RTMB)

x <- advector(1:3)

## Method call: outer -> Ops with re-cycling
expect_true(is(outer(x,x,"+"),   "advector"))
expect_true(is(outer(1:3,x,"+"), "advector"))
expect_true(is(outer(x,1:3,"+"), "advector"))

## Method call: outer -> tcrossprod
expect_true(is(outer(x,x),   "advector"))
expect_true(is(outer(1:3,x), "advector"))
expect_true(is(outer(x,1:3), "advector"))

## Method call: kronecker -> outer -> tcrossprod
expect_true(is(kronecker(x,x),   "advector"))
expect_true(is(kronecker(1:3,x), "advector"))
expect_true(is(kronecker(x,1:3), "advector"))
