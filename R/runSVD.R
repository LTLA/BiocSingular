#' @export
setMethod("runSVD", "missing", function(..., BSPARAM) {
    runSVD(..., BSPARAM=BSPARAM)
})

#' @export
setMethod("runSVD", "ExactParam", function(..., BSPARAM) {
    runExactSVD(..., deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM))
})

#' @export
setMethod("runSVD", "IrlbaParam", function(..., BSPARAM) {
    do.call(runIrlbaSVD, c(list(..., deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM), extra.work=ip_extra(BSPARAM)), bsargs(BSPARAM)))
})

#' @export
setMethod("runSVD", "RandomParam", function(..., BSPARAM) {
    do.call(runRandomSVD, c(list(..., deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM)), bsargs(BSPARAM)))
})
