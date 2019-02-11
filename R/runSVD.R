#' @export
#' @importFrom BiocParallel SerialParam
setMethod("runSVD", "missing", function(x, k, nu=k, nv=k, center=NULL, scale=NULL, BPPARAM=SerialParam(), ..., BSPARAM) {
    runSVD(x=x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ..., BSPARAM=BSPARAM)
})

#' @export
#' @importFrom BiocParallel SerialParam
setMethod("runSVD", "ExactParam", function(x, k, nu=k, nv=k, center=NULL, scale=NULL, BPPARAM=SerialParam(), ..., BSPARAM) {
    runExactSVD(x=x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ..., deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM))
})

#' @export
#' @importFrom BiocParallel SerialParam
setMethod("runSVD", "IrlbaParam", function(x, k, nu=k, nv=k, center=NULL, scale=NULL, BPPARAM=SerialParam(), ..., BSPARAM) {
    do.call(runIrlbaSVD, 
        c(
            list(
                x=x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ...,
                deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM), extra.work=ip_extra(BSPARAM)
            ), 
            bsargs(BSPARAM)
         )
    )
})

#' @export
#' @importFrom BiocParallel SerialParam
setMethod("runSVD", "RandomParam", function(x, k, nu=k, nv=k, center=NULL, scale=NULL, BPPARAM=SerialParam(), ..., BSPARAM) {
    do.call(runRandomSVD, 
        c(
            list(
                x=x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ...,
                deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM)
            ), 
            bsargs(BSPARAM)
        )
    )
})
