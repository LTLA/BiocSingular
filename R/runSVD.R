#' @export
setMethod("runSVD", "missing", function(..., BSPARAM) {
    runSVD(..., BSPARAM=ExactParam())
})

#' @export
setMethod("runSVD", "ExactParam", function(..., BSPARAM) {
    runExactSVD(..., deferred=bsp_deferred(BSPARAM), fold=bsp_fold(BSPARAM))
})

#' @export
setMethod("runSVD", "IrlbaParam", function(..., BSPARAM) {
    do.call(runIrlbaSVD, c(list(..., deferred=bsp_deferred(BSPARAM), fold=bsp_fold(BSPARAM), extra.work=ip_extra(BSPARAM)), bsp_args(BSPARAM)))
})

#' @export
setMethod("runSVD", "RandomParam", function(..., BSPARAM) {
    do.call(runRandomSVD, c(list(..., deferred=bsp_deferred(BSPARAM), fold=bsp_fold(BSPARAM)), bsp_args(BSPARAM)))
})
