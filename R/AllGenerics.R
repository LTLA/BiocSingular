#' @export
#' @importFrom BiocParallel SerialParam
setGeneric("runSVD", signature=c("BSPARAM"), 
    function(x, k, nu=k, nv=k, center=FALSE, scale=FALSE, BPPARAM=SerialParam(), ..., BSPARAM=ExactParam()) 
        standardGeneric("runSVD")
)

#' @export
setGeneric("runPCA", signature=c("x"), function(x, ...) standardGeneric("runPCA"))

########################################
# Getters/setters for linearEmbedding.

#' @export
setGeneric("sampleFactors", function(x, ...) standardGeneric("sampleFactors"))

#' @export
setGeneric("sampleFactors<-", function(x, ..., value) standardGeneric("sampleFactors<-"))

#' @export
setGeneric("featureLoadings", function(x, ..., withDimnames=TRUE) standardGeneric("featureLoadings"))

#' @export
setGeneric("featureLoadings<-", function(x, ..., value) standardGeneric("featureLoadings<-"))

#' @export
setGeneric("factorData", function(x, ..., withDimnames=TRUE) standardGeneric("factorData"))

#' @export
setGeneric("factorData<-", function(x, ..., value) standardGeneric("factorData<-"))
