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

        expect_equal(rowSums(test$def), rowSums(test$ref))
        expect_equal(colSums(test$def), colSums(test$ref))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        expect_identical(spawn_names, colnames(test$def))
    }

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(DeferredMatrix(y, center=1), "length of 'center' must equal")
    expect_error(DeferredMatrix(y, scale=1), "length of 'scale' must equal")
})

set.seed(1000011)
test_that("DeferredMatrix subsetting works as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$def))
        j <- sample(ncol(test$def))
        expect_identical(as.matrix(test$def[i,]), test$ref[i,])
        expect_identical(as.matrix(test$def[,j]), test$ref[,j])
        expect_identical(as.matrix(test$def[i,j]), test$ref[i,j])
        
        # Works with zero dimensions.
        expect_identical(as.matrix(test$def[0,]), test$ref[0,])
        expect_identical(as.matrix(test$def[,0]), test$ref[,0])
        expect_identical(as.matrix(test$def[0,0]), test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_identical(test$def[i[1],], test$ref[i[1],])
        expect_identical(test$def[,j[1]], test$ref[,j[1]])
        expect_identical(as.matrix(test$def[i[1],drop=FALSE]), test$ref[i[1],,drop=FALSE])
        expect_identical(as.matrix(test$def[,j[1],drop=FALSE]), test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$def)
        expect_identical(t(alt[,i]), test$def[i,])
        expect_identical(t(alt[j,]), test$def[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_identical(as.matrix(test$def[,ch]), test$ref[,ch])
    }
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

test_that("DeferredMatrix dual crossprod fails as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- DeferredMatrix(y, NULL, NULL)
    expect_error(crossprod(bs.y, bs.y), "not yet supported")
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

test_that("DeferredMatrix dual tcrossprod fails as expected", {
    y <- matrix(rnorm(400), ncol=20)
    bs.y <- DeferredMatrix(y, NULL, NULL)
    expect_error(tcrossprod(bs.y, bs.y), "not yet supported")
})

##########################

test_that("nested DeferredMatrix works as expected", {
    for (it in 1:4) {
        # Setting up the scenario, with and without transposition.
        a1 <- matrix(rnorm(400), ncol=20)

        c1 <- rnorm(20)
        s1 <- runif(20)
        r1 <- DeferredMatrix(a1, c1, s1)
        ref1 <- scale(a1, c1, s1)
        if (it%%2L==0L) {
            r1 <- t(r1)
            ref1 <- t(ref1)
        }
        
        c2 <- rnorm(20)
        s2 <- runif(20)
        r2 <- DeferredMatrix(r1, c2, s2)
        ref2 <- scale(ref1, c2, s2)
        if (it > 2L) {
            r2 <- t(r2)
            ref2 <- t(ref2)
        }

        attr(ref2, "scaled:center") <- NULL
        attr(ref2, "scaled:scale") <- NULL

        # Coercion works.
        expect_equal(ref2, as.matrix(r2))

        # Basic stats work.
        expect_equal(rowSums(ref2), rowSums(r2))
        expect_equal(colSums(ref2), colSums(r2))
       
        # Multiplication works.        
        y <- matrix(rnorm(20*2), ncol=2)        
        expect_equal(ref2 %*% y, r2 %*% y)
        expect_equal(t(y) %*% ref2, t(y) %*% r2)

        # Cross product.
        y <- matrix(rnorm(20*2), ncol=2)        
        expect_equal(crossprod(ref2), crossprod(r2))
        expect_equal(crossprod(ref2, y), crossprod(r2, y))
        expect_equal(crossprod(y, ref2), crossprod(y, r2))

        # Transposed cross product.
        y <- matrix(rnorm(20*2), nrow=2) 
        expect_equal(tcrossprod(ref2), tcrossprod(r2))
        expect_equal(tcrossprod(ref2, y), tcrossprod(r2, y))
        expect_equal(tcrossprod(y, ref2), tcrossprod(y, r2))
    }
})

set.seed(1200001)
test_that("deep testing of tcrossproduct internals: special mult", {
    NC <- 20
    NR <- 10
    c <- runif(NC)
    s <- runif(NR) # NOT the scale for 'r', but for the parent DeferredMatrix, which has transposed dimensions.
    r <- matrix(rnorm(NC*NR), ncol=NC)

    ref <- t(matrix(c, NR, NC, byrow=TRUE)) %*% (r/s^2)
    out <- BiocSingular:::.internal_mult_special(c, s, r)
    expect_equal(ref, out)

    # Trying with DeferredMatrices.
    c1 <- runif(NC)
    s1 <- runif(NC)
    r1 <- DeferredMatrix(r, c1, s1)
    out <- BiocSingular:::.internal_mult_special(c, s, r1)
    ref <- BiocSingular:::.internal_mult_special(c, s, as.matrix(r1))
    expect_equal(ref, out)

    # Now with transposition.
    c2 <- runif(NR)
    s2 <- runif(NR)
    r2 <- DeferredMatrix(t(r), c2, s2)
    r2 <- t(r2)
    out <- BiocSingular:::.internal_mult_special(c, s, r2)
    ref <- BiocSingular:::.internal_mult_special(c, s, as.matrix(r2))
    expect_equal(ref, out)
})

set.seed(1200002)
test_that("deep testing of tcrossproduct internals: scaled tcrossprod", {
    NC <- 30
    NR <- 15
    s <- runif(NC) 
    r <- matrix(rnorm(NC*NR), ncol=NC)

    ref <- crossprod(t(r)/s)
    out <- BiocSingular:::.internal_tcrossprod(r, s)
    expect_equal(ref, out)

    # Trying with DeferredMatrices.
    c1 <- runif(NC)
    s1 <- runif(NC)
    r1 <- DeferredMatrix(r, c1, s1)
    out <- BiocSingular:::.internal_tcrossprod(r1, s)
    ref <- BiocSingular:::.internal_tcrossprod(as.matrix(r1), s)
    expect_equal(ref, out)

    # With transposition.
    c2 <- runif(NR)
    s2 <- runif(NR)
    r2 <- DeferredMatrix(t(r), c2, s2)
    r2 <- t(r2)
    out <- BiocSingular:::.internal_tcrossprod(r2, s)
    ref <- BiocSingular:::.internal_tcrossprod(as.matrix(r2), s)
    expect_equal(ref, out)

    # Now with nested DeferredMatrices.
    r3 <- DeferredMatrix(r1, rnorm(ncol(r1)), runif(ncol(r1)))
    out <- BiocSingular:::.internal_tcrossprod(r3, s)
    ref <- BiocSingular:::.internal_tcrossprod(as.matrix(r3), s)
    expect_equal(ref, out)

    r4 <- DeferredMatrix(r2, rnorm(ncol(r2)), runif(ncol(r2)))
    out <- BiocSingular:::.internal_tcrossprod(r4, s)
    ref <- BiocSingular:::.internal_tcrossprod(as.matrix(r4), s)
    expect_equal(ref, out)

    r5 <- t(DeferredMatrix(t(r1), rnorm(nrow(r1)), runif(nrow(r1), 10, 20)))
    out <- BiocSingular:::.internal_tcrossprod(r5, s)
    ref <- BiocSingular:::.internal_tcrossprod(as.matrix(r5), s)
    expect_equal(ref, out)
})

##########################

library(DelayedArray)
test_that("DelayedMatrix wrapping works", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        X <- DelayedArray(test$def)
        expect_equal(as.matrix(X^2), test$ref^2)
        expect_equal(as.matrix(X + 1), test$ref + 1)

        v <- nrow(X)
        expect_equal(as.matrix(X + v), test$ref + v)
        expect_equal(as.matrix(X * v), test$ref * v)

        expect_identical(type(X), "double")
    }
})
