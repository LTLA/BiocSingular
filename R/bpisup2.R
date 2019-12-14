#' @importFrom methods is
#' @importClassesFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bpisup
.bpisup2 <- function(BPPARAM) {
    bpisup(BPPARAM) || is(BPPARAM, "MulticoreParam")
}
