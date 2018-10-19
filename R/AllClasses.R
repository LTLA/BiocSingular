#' @export
setClass("BiocSingularParam", contains="VIRTUAL")

#' @export
setClass("ExactParam", contains="BiocSingularParam")

#' @export
setClass("IrlbaParam", contains="BiocSingularParam", slots=c(extra.work="integer", maxit="integer", tol="numeric"))
