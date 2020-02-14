#' @importFrom BiocParallel SerialParam
.FUN_GENERATOR <- function(FUN, ARGS) 
# Generates the function to be used in the method,
# to avoid re-typing the whole signature.
{
    function(x, k, nu=k, nv=k, center=FALSE, scale=FALSE, BPPARAM=SerialParam(), ..., BSPARAM) {
        do.call(FUN, 
            c(
                list(x=x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ...),
                ARGS(BSPARAM)
             )
        )
    }
}

#' @export
# Uses the default ExactParam() in the generic.
setMethod("runSVD", "missing", .FUN_GENERATOR(runSVD, 
    function(BSPARAM) list(BSPARAM=BSPARAM) 
))

#' @export
setMethod("runSVD", "ExactParam", .FUN_GENERATOR(runExactSVD, 
    function(BSPARAM) list(deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM))
))

#' @export
setMethod("runSVD", "IrlbaParam", .FUN_GENERATOR(runIrlbaSVD, 
    function(BSPARAM) {
        c(list(deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM), extra.work=ip_extra(BSPARAM)), bsargs(BSPARAM))
    }
))

#' @export
setMethod("runSVD", "RandomParam", .FUN_GENERATOR(runRandomSVD, 
    function(BSPARAM) {            
        c(list(deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM)), bsargs(BSPARAM))
    }
))

#' @export
#' @importFrom methods is hasMethod
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importFrom DelayedArray seed
setMethod("runSVD", "FastAutoParam", function(x, k, nu=k, nv=k, center=FALSE, scale=FALSE, 
    BPPARAM=SerialParam(), ..., BSPARAM) 
{
    BSPARAMFUN <- IrlbaParam
    if (is(x, "DelayedMatrix")) {
        if (class(x)[1]=="DelayedMatrix" || !hasMethod("%*%", class(x))) { 
            BSPARAMFUN <- RandomParam
        }
    }
    runSVD(x, k=k, nu=nu, nv=nv, center=center, scale=scale, BPPARAM=BPPARAM, ...,
        BSPARAM=BSPARAMFUN(deferred=bsdeferred(BSPARAM), fold=bsfold(BSPARAM)))
})
