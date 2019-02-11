#' @export
#' @importFrom BiocParallel SerialParam
setGeneric("runSVD", signature=c("BSPARAM"), 
    function(x, k, nu=k, nv=k, center=NULL, scale=NULL, BPPARAM=SerialParam(), ..., BSPARAM=ExactParam()) 
        standardGeneric("runSVD")
)

#' @export
setGeneric("runPCA", signature=c("x"), function(x, ...) standardGeneric("runPCA"))
