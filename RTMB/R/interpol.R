interpol2Dfun <- function(z, xlim=c(1,nrow(z)), ylim=c(1,ncol(z)), ...) {
    con <- list(...)
    if (is.null(con$R))
        con$R <- 2
    ptr <- RTMB:::ip2D(z, xlim, ylim, con)
    function(x, y) {
        if (inherits(x, "advector") || inherits(y, "advector")) {
            RTMB:::ip2D_eval_ad(ptr, advector(x), advector(y))
        } else {
            RTMB:::ip2D_eval_num(ptr, x, y)
        }
    }
}
