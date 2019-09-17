#' @export
bsparam <- function() {
    getOption("BiocSingularParam.default", FastAutoParam())
}
