# Tests the DeferredMatrix implementation.
# library(testthat); library(BiocSingular); source("test-deferred.R")

scale_and_center <- function(y, ref, code) {
    center <- scale <- NULL

    if (code==1L) {
        center <- colMeans(ref)
        scale <- runif(ncol(ref))
        ref <- scale(ref, center=center, scale=scale)
    } else if (code==2L) {
        center <- rnorm(ncol(ref))
        ref <- scale(ref, center=center, scale=FALSE)
    } else if (code==3L) {
        scale <- runif(ncol(ref))
        ref <- scale(ref, center=FALSE, scale=scale)
    }

    # Getting rid of excess attributes.
    attr(ref, "scaled:center") <- NULL
    attr(ref, "scaled:scale") <- NULL

    def <- DeferredMatrix(y, center=center, scale=scale)
    list(def=def, ref=ref)
}

spawn_scenarios_basic <- function(NR, NC, CREATOR, REALIZER) {
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:4) {
            if (trans) { 
                # Ensure output matrix has NR rows and NC columns after t().
                y <- CREATOR(NC, NR)
            } else {
                y <- CREATOR(NR, NC)
            }
            ref <- REALIZER(y) 

            adjusted <- scale_and_center(y, ref, it) 
            if (trans) {
                adjusted$def <- t(adjusted$def)
                adjusted$ref <- t(adjusted$ref)
            }
            
            output[[counter]] <- adjusted
            counter <- counter+1L
        }
    }
    output
}

spawn_scenarios <- function(NR=50, NC=20) {
    c(
        spawn_scenarios_basic(NR, NC,
            CREATOR=function(r, c) {
                matrix(rnorm(r*c), ncol=c)
            }, 
            REALIZER=identity
        ),
        spawn_scenarios_basic(NR, NC, 
            CREATOR=function(r, c) {
                Matrix::rsparsematrix(r, c, 0.1)
            },
            REALIZER=as.matrix
        )
    )
}

##########################

set.seed(100001)
test_that("DeferredMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_s4_class(test$def, "DeferredMatrix")
        expect_identical(test$def, DeferredMatrix(DelayedArray::seed(test$def)))

        expect_identical(dim(test$def), dim(test$ref))
        expect_identical(extract_array(test$def, list(1:10, 1:10)), test$ref[1:10, 1:10])
        expect_identical(extract_array(test$def, list(1:10, NULL)), test$ref[1:10,])
        expect_identical(extract_array(test$def, list(NULL, 1:10)), test$ref[,1:10])
        expect_identical(as.matrix(test$def), test$ref)

        expect_equal(rowSums(test$def), rowSums(test$ref))
        expect_equal(colSums(test$def), colSums(test$ref))
        expect_equal(rowMeans(test$def), rowMeans(test$ref))
        expect_equal(colMeans(test$def), colMeans(test$ref))

        tdef <- t(test$def)
        expect_s4_class(tdef, "DeferredMatrix") # still a DefMat!
        expect_identical(t(tdef), test$def)
        expect_identical(as.matrix(tdef), t(test$ref))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        expect_identical(spawn_names, colnames(test$def))
        expect_s4_class(test$def, "DeferredMatrix") # still a DefMat!
    }
})

set.seed(10000101)
test_that("DeferredMatrix silly inputs work as expected", {
    default <- DeferredMatrix()
    expect_identical(dim(default), c(0L, 0L))
    val <- as.matrix(default)
    dimnames(val) <- NULL
    expect_identical(val, matrix(0,0,0))

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(DeferredMatrix(y, center=1), "length of 'center' must equal")
    expect_error(DeferredMatrix(y, scale=1), "length of 'scale' must equal")
})

set.seed(1000011)
test_that("DeferredMatrix subsetting works as expected", {
    expect_identical_and_defmat <- function(x, y) {
        expect_s4_class(x, "DeferredMatrix") # class is correctly preserved by direct seed modification.
        expect_identical(as.matrix(x), y)
    }

    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$def))
        j <- sample(ncol(test$def))
        expect_identical_and_defmat(test$def[i,], test$ref[i,])
        expect_identical_and_defmat(test$def[,j], test$ref[,j])
        expect_identical_and_defmat(test$def[i,j], test$ref[i,j])
        
        # Works with zero dimensions.
        expect_identical_and_defmat(test$def[0,], test$ref[0,])
        expect_identical_and_defmat(test$def[,0], test$ref[,0])
        expect_identical_and_defmat(test$def[0,0], test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_identical(test$def[i[1],], test$ref[i[1],])
        expect_identical(test$def[,j[1]], test$ref[,j[1]])
        expect_identical_and_defmat(test$def[i[1],drop=FALSE], test$ref[i[1],,drop=FALSE])
        expect_identical_and_defmat(test$def[,j[1],drop=FALSE], test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$def)
        expect_identical(t(alt[,i]), test$def[i,])
        expect_identical(t(alt[j,]), test$def[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_identical_and_defmat(test$def[,ch], test$ref[,ch])
    }
})

##########################
# Defining a class that can't do anything but get multiplied.
# This checks that there isn't any hidden DelayedArray realization 
# happening, which would give the same results but slower.

setClass("CrippledMatrix", slots=c(x="matrix"))

setMethod("dim", c("CrippledMatrix"), function(x) dim(x@x))

setMethod("colSums", c("CrippledMatrix"), function(x) colSums(x@x))

setMethod("rowSums", c("CrippledMatrix"), function(x) rowSums(x@x))

setMethod("sweep", c("CrippledMatrix"), function (x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
    sweep(x@x, MARGIN, STATS, FUN, check.margin, ...)
})

setMethod("%*%", c("CrippledMatrix", "ANY"), function(x, y) x@x %*% y)

setMethod("%*%", c("ANY", "CrippledMatrix"), function(x, y) x %*% y@x)

setMethod("crossprod", c("CrippledMatrix", "missing"), function(x, y) crossprod(x@x))

setMethod("crossprod", c("CrippledMatrix", "ANY"), function(x, y) crossprod(x@x, y))

setMethod("crossprod", c("ANY", "CrippledMatrix"), function(x, y) crossprod(x, y@x))

setMethod("tcrossprod", c("CrippledMatrix", "missing"), function(x, y) tcrossprod(x@x))

setMethod("tcrossprod", c("CrippledMatrix", "ANY"), function(x, y) tcrossprod(x@x, y))

setMethod("tcrossprod", c("ANY", "CrippledMatrix"), function(x, y) tcrossprod(x, y@x))

spawn_extra_scenarios <- function(NR=50, NC=20) {
    c(
        spawn_scenarios(NR, NC),
        spawn_scenarios_basic(NR, NC, 
            CREATOR=function(r, c) { 
                new("CrippledMatrix", x=matrix(runif(NR*NC), ncol=NC))
            }, 
            REALIZER=function(x) x@x
        )
    )
}

expect_equal_product <- function(x, y) {
    expect_s4_class(x, "DelayedMatrix")
    X <- as.matrix(x)

    # standardize NULL dimnames.
    if (all(lengths(dimnames(X))==0L)) dimnames(X) <- NULL 
    if (all(lengths(dimnames(y))==0L)) dimnames(y) <- NULL 
    expect_equal(X, y)
}

##########################

test_that("DeferredMatrix right multiplication works as expected", {
    possibles <- spawn_extra_scenarios(100, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), ncol=10)
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)
    }
})

test_that("DeferredMatrix left multiplication works as expected", {
    possibles <- spawn_extra_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), nrow=10)
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)
    }
})

test_that("DeferredMatrix dual multiplication works as expected", {
    # Not using the CrippledMatrix here; some scaling of the inner matrix is unavoidable
    # when the inner matrix is _not_ a DeferredMatrix but is being multiplied by one.
    possibles1 <- spawn_scenarios(10, 20)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(test1$def %*% test2$def, test1$ref %*% test2$ref)

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(test1$def[0,] %*% test2$def, test1$ref[0,] %*% test2$ref)
            expect_equal_product(test1$def %*% test2$def[,0], test1$ref %*% test2$ref[,0])
            expect_equal_product(test1$def[,0] %*% test2$def[0,], test1$ref[,0] %*% test2$ref[0,])
            expect_equal_product(test1$def[0,] %*% test2$def[,0], test1$ref[0,] %*% test2$ref[,0])
        }
    }
})

##########################

test_that("DeferredMatrix lonely crossproduct works as expected", {
    possibles <- spawn_extra_scenarios(90, 30)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal_product(crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("DeferredMatrix crossproduct from right works as expected", {
    possibles <- spawn_extra_scenarios(60, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))
    }
})

test_that("DeferredMatrix crossproduct from left works as expected", {
    possibles <- spawn_extra_scenarios(40, 100)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))
    }
})

test_that("DeferredMatrix dual crossprod works as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- DeferredMatrix(y, NULL, NULL)
    expect_identical(crossprod(bs.y, bs.y), crossprod(bs.y))
})

##########################

test_that("DeferredMatrix lonely tcrossproduct works as expected", {
    possibles <- spawn_extra_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal_product(tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("DeferredMatrix tcrossproduct from right works as expected", {
    possibles <- spawn_extra_scenarios(60, 70)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector (this doesn't work).
        z <- rnorm(ncol(ref.y))
        expect_error(tcrossprod(bs.y, z), "non-conformable")
        expect_error(tcrossprod(ref.y, z), "non-conformable")

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))
    }
})

test_that("DeferredMatrix tcrossproduct from left works as expected", {
    possibles <- spawn_extra_scenarios(80, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))
    }
})

test_that("DeferredMatrix dual tcrossprod works as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- DeferredMatrix(y, NULL, NULL)
    expect_identical(tcrossprod(bs.y, bs.y), tcrossprod(bs.y))
})

##########################

wrap_in_DefMat <- function(input, reference) 
# Wrapping an input matrix in a DeferredMatrix.
{
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:4) {
            if (trans) { 
                y <- t(input)
                ref <- t(reference)
            } else {
                ref <- reference
                y <- input
            }

            adjusted <- scale_and_center(y, ref, it)
            if (trans) {
                adjusted$def <- t(adjusted$def)
                adjusted$ref <- t(adjusted$ref)
            }

            output[[counter]] <- adjusted
            counter <- counter+1L
        }
    }
    output
}

test_that("nested DeferredMatrix works as expected", {
    basic <- matrix(rnorm(400), ncol=20)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested DefMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_DefMat(current$def, current$ref)
        }

        # Testing each one of the newly created DefMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]

            # Coercion works.
            expect_equal(as.matrix(test$def), test$ref)

            # Basic stats work.
            expect_equal(rowSums(test$ref), rowSums(test$def))
            expect_equal(colSums(test$ref), colSums(test$def))

            # Multiplication works.
            y <- matrix(rnorm(20*2), ncol=2)
            expect_equal_product(test$def %*% y, test$ref %*% y)
            expect_equal_product(t(y) %*% test$def, t(y) %*% test$ref)

            # Cross product.
            y <- matrix(rnorm(20*2), ncol=2)
            expect_equal_product(crossprod(test$def), crossprod(test$ref))
            expect_equal_product(crossprod(test$def, y), crossprod(test$ref, y))
            expect_equal_product(crossprod(y, test$def), crossprod(y, test$ref))

            # Transposed cross product.
            y <- matrix(rnorm(20*2), nrow=2) 
            expect_equal_product(tcrossprod(test$def), tcrossprod(test$ref))
            expect_equal_product(tcrossprod(test$def, y), tcrossprod(test$ref, y))
            expect_equal_product(tcrossprod(y, test$def), tcrossprod(y, test$ref))
        }
    }
})

set.seed(1200001)
test_that("deep testing of tcrossproduct internals: special mult", {
    NR <- 20
    NC <- 10
    basic <- matrix(rnorm(NC*NR), ncol=NC)
    c <- runif(NC)
    s <- runif(NR)

    ref <- t(matrix(c, NR, NC, byrow=TRUE)) %*% (basic/s^2)
    out <- BiocSingular:::.internal_mult_special(c, s, basic)
    expect_equal(ref, out)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested DefMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_DefMat(current$def, current$ref)
        }

        # Testing each one of the newly created nested DefMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]
            ref <- t(matrix(c, NR, NC, byrow=TRUE)) %*% (test$ref/s^2)
            out <- BiocSingular:::.internal_mult_special(c, s, test$def)
            expect_equal(ref, out)
        }
    }
})

set.seed(1200002)
test_that("deep testing of tcrossproduct internals: scaled tcrossprod", {
    NC <- 30
    NR <- 15
    s <- runif(NC) 
    basic <- matrix(rnorm(NC*NR), ncol=NC)

    ref <- crossprod(t(basic)/s)
    out <- BiocSingular:::.internal_tcrossprod(basic, s)
    expect_equal(ref, out)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested DefMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_DefMat(current$def, current$ref)
        }

        # Testing each one of the newly created nested DefMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]
            ref <- crossprod(t(test$ref)/s)
            out <- BiocSingular:::.internal_tcrossprod(test$def, s)
            expect_equal(ref, out)
        }
    }
})

##########################

test_that("DelayedMatrix wrapping works", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        expect_equal_product(test$def+1, test$ref+1)

        v <- rnorm(nrow(test$def))
        expect_equal_product(test$def+v, test$ref+v)
        expect_equal_product(test$def*v, test$ref*v)

        w <- rnorm(ncol(test$def))
        expect_equal_product(sweep(test$def, 2, w, "*"), sweep(test$ref, 2, w, "*"))
    }
})
