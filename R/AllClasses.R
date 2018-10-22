#' @export
setClass("BiocSingularParam", contains="VIRTUAL", slots=c(fold="numeric"))

#' @export
setClass("ExactParam", contains="BiocSingularParam")

#' @export
setClass("IrlbaParam", contains="BiocSingularParam", slots=c(extra.work="integer", args="list"))

#' @export
setClass("RandomParam", contains="BiocSingularParam", slots=c(args="list"))
