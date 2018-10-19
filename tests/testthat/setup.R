expect_equal_besides_sign <- function(left, right, ...) {
    ratio <- left/right
    right2 <- sweep(right, 2, sign(ratio[1,]), FUN="*")
    expect_equal(left, right2, ...)
}

expect_equal_svd <- function(left, right, ...) {
    expect_equal_besides_sign(left$u, right$u, ...)
    expect_equal_besides_sign(left$v, right$v, ...)
    expect_equal(left$d, right$d, ...)
}
