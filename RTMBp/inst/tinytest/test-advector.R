## Test indirect advector methods

library(RTMBp)

x <- advector(1:3)

## Method call: outer -> Ops with re-cycling
expect_true(is(outer(x,x,"+"),   "advector"))
expect_true(is(outer(1:3,x,"+"), "advector"))
expect_true(is(outer(x,1:3,"+"), "advector"))

## See ?crossprod: **These are generic functions since R 4.4.0**
if (getRversion() >= "4.4.0") {
    ## Method call: outer -> tcrossprod
    expect_true(is(outer(x,x),   "advector"))
    expect_true(is(outer(1:3,x), "advector"))
    expect_true(is(outer(x,1:3), "advector"))

    ## Method call: kronecker -> outer -> tcrossprod
    expect_true(is(kronecker(x,x),   "advector"))
    expect_true(is(kronecker(1:3,x), "advector"))
    expect_true(is(kronecker(x,1:3), "advector"))
}
