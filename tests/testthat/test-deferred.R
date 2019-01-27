# Tests the DeferredMatrix implementation.
# library(testthat); library(BiocSingular); source("test-deferred.R")

library(BiocSingular)
spawn_scenarios <- function(NR=50, NC=20) {
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:4) {
            y <- matrix(rnorm(NR*NC), ncol=NC)
            center <- scale <- NULL
    
            if (it==1L) {
                center <- colMeans(y)
                scale <- runif(ncol(y))
                ref.y <- scale(y, center=center, scale=scale)
            } else if (it==2L) {
                center <- rnorm(ncol(y))
                ref.y <- scale(y, center=center, scale=FALSE)
            } else if (it==3L) {
                scale <- runif(ncol(y))
                ref.y <- scale(y, center=FALSE, scale=scale)
            } else {
                ref.y <- y
            }
  
            # Getting rid of excess attributes.
            attr(ref.y, "scaled:center") <- NULL
            attr(ref.y, "scaled:scale") <- NULL

            def <- DeferredMatrix(y, center=center, scale=scale)
            if (trans) {
                def <- t(def)
                ref.y <- t(ref.y)
            }
            
            output[[counter]] <- list(ref=ref.y, def=def)
            counter <- counter+1L
        }
    }
    output
}

##########################

set.seed(100001)
test_that("DeferredMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_identical(dim(test$def), dim(test$ref))
        expect_identical(length(test$def), length(test$ref))
        expect_identical(as.matrix(test$def), test$ref)

        # Matrix subsetting works as expected.
        i <- sample(nrow(test$def))
        j <- sample(ncol(test$def))
        expect_identical(as.matrix(test$def[i,]), test$ref[i,])
        expect_identical(as.matrix(test$def[,j]), test$ref[,j])
        expect_identical(as.matrix(test$def[i,j]), test$ref[i,j])
        
        # Dimension dropping works as expected.
        expect_identical(test$def[i[1],], test$ref[i[1],])
        expect_identical(test$def[,j[1]], test$ref[,j[1]])
        expect_identical(as.matrix(test$def[i[1],drop=FALSE]), test$ref[i[1],,drop=FALSE])
        expect_identical(as.matrix(test$def[,j[1],drop=FALSE]), test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$def)
        expect_identical(t(alt[,i]), test$def[i,])
        expect_identical(t(alt[j,]), test$def[,j])
    }

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(DeferredMatrix(y, center=1), "length of 'center' must equal")
    expect_error(DeferredMatrix(y, scale=1), "length of 'scale' must equal")
})

##########################

test_that("DeferredMatrix right multiplication works as expected", {
    possibles <- spawn_scenarios(100, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal(bs.y %*% z, ref.y %*% z)

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), ncol=10)
        expect_equal(bs.y %*% z, ref.y %*% z)

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=ncol(ref.y))
        expect_equal(bs.y %*% z, ref.y %*% z)
    }
})

test_that("DeferredMatrix left multiplication works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal(z %*% bs.y, z %*% ref.y)

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), nrow=10)
        expect_equal(z %*% bs.y, z %*% ref.y)

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=nrow(ref.y))
        expect_equal(z %*% bs.y, z %*% ref.y)
    }
})

test_that("DeferredMatrix dual multiplication fails as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- DeferredMatrix(y, NULL, NULL)
    expect_error(bs.y %*% bs.y, "not yet supported")
})

##########################

test_that("DeferredMatrix lonely crossproduct works as expected", {
    possibles <- spawn_scenarios(90, 30)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal(crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("DeferredMatrix left crossproduct works as expected", {
    possibles <- spawn_scenarios(60, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal(crossprod(bs.y, z), crossprod(ref.y, z))
    }
})

test_that("DeferredMatrix right crossproduct works as expected", {
    possibles <- spawn_scenarios(40, 100)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal(crossprod(z, bs.y), crossprod(z, ref.y))
    }
})

##########################

test_that("DeferredMatrix lonely tcrossproduct works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal(tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("DeferredMatrix left tcrossproduct works as expected", {
    possibles <- spawn_scenarios(60, 70)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector (this doesn't work).
        z <- rnorm(ncol(ref.y))
        expect_error(tcrossprod(bs.y, z), "non-conformable")
        expect_error(tcrossprod(ref.y, z), "non-conformable")

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal(tcrossprod(bs.y, z), tcrossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal(tcrossprod(bs.y, z), tcrossprod(ref.y, z))
    }
})

test_that("DeferredMatrix right crossproduct works as expected", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal(tcrossprod(z, bs.y), tcrossprod(z, ref.y))
    }
})
