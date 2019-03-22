#' @export
#' @import methods
setClass("BiocSingularParam", contains="VIRTUAL", slots=c(deferred="logical", fold="numeric"))

#' @export
setClass("ExactParam", contains="BiocSingularParam")

#' @export
setClass("IrlbaParam", contains="BiocSingularParam", slots=c(extra.work="integer", args="list"))

#' @export
setClass("RandomParam", contains="BiocSingularParam", slots=c(args="list"))

#' @export
setClass("DeferredMatrixSeed", slots=c(.matrix="ANY", center="numeric", scale="numeric", use_center="logical", use_scale="logical", transposed="logical"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("DeferredMatrix",
    contains="DelayedMatrix",
    representation(seed="DeferredMatrixSeed")
)

#' @export
setClass("LowRankMatrixSeed", slots=c(rotation="ANY", components="ANY"))

#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
setClass("LowRankMatrix",
    contains="DelayedMatrix",
    representation(seed="LowRankMatrixSeed")
)
