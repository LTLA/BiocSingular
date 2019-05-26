# Tests the ResidualMatrix implementation.
# library(testthat); library(BiocSingular); source("setup.R"); source("test-residual.R")

generate_design <- function(nobs, code) {
    if (code==1L) {
        design <- NULL
    } else if (code==2L) {
        cov <- rnorm(nobs)
        design <- model.matrix(~cov)
    } else if (code==3L) {
        g <- factor(rep(1:3, length.out=nobs))
        design <- model.matrix(~0 + g)
    } else if (code==4L) {
        cov <- rnorm(nobs)
        g <- factor(rep(1:2, length.out=nobs))
        design <- model.matrix(~0 + g + cov)
    } else if (code==5L) {
        design <- cbind(rnorm(nobs))
    }
    design
}

spawn_scenarios_basic <- function(NR, NC, CREATOR, REALIZER) {
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:5) {
            if (trans) {
                # Ensure output matrix has NR rows and NC columns after t().
                y <- CREATOR(NC, NR)
            } else {
                y <- CREATOR(NR, NC)
            }
            ref <- REALIZER(y)

            # Run through a host of different design matrices.
            design <- generate_design(nrow(y), it)

            res <- ResidualMatrix(y, design)
            if (is.null(design)) {
                ref <- as.matrix(y)
                ref <- sweep(ref, 2, colMeans(ref), "-")
            } else {
                ref <- lm.fit(x=design, y=as.matrix(y))$residuals
            }

            if (trans) {
                res <- t(res)
                ref <- t(ref)
            }

            output[[counter]] <- list(res=res, ref=ref)
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
test_that("ResidualMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_s4_class(test$res, "ResidualMatrix")
        expect_identical(test$res, ResidualMatrix(DelayedArray::seed(test$res)))

        expect_identical(dim(test$res), dim(test$ref))
        expect_equal(extract_array(test$res, list(1:10, 1:10)), test$ref[1:10, 1:10])
        expect_equal(extract_array(test$res, list(1:10, NULL)), test$ref[1:10,])
        expect_equal(extract_array(test$res, list(NULL, 1:10)), test$ref[,1:10])
        expect_equal(as.matrix(test$res), test$ref)

        expect_equal(rowSums(test$res), rowSums(test$ref))
        expect_equal(colSums(test$res), colSums(test$ref))
        expect_equal(rowMeans(test$res), rowMeans(test$ref))
        expect_equal(colMeans(test$res), colMeans(test$ref))

        tdef <- t(test$res)
        expect_s4_class(tdef, "ResidualMatrix") # still a ResMat!
        expect_equal(t(tdef), test$res)
        expect_equal(as.matrix(tdef), t(test$ref))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$res)))
        colnames(test$res) <- spawn_names
        expect_identical(spawn_names, colnames(test$res))
        expect_s4_class(test$res, "ResidualMatrix") # still a ResMat!
    }
})

set.seed(10000101)
test_that("ResidualMatrix silly inputs work as expected", {
    default <- ResidualMatrix()
    expect_identical(dim(default), c(0L, 0L))
    val <- as.matrix(default)
    dimnames(val) <- NULL
    expect_identical(val, matrix(0,0,0))

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(ResidualMatrix(y, design=cbind(1)), "nrow.* should be equal")
})

set.seed(1000011)
test_that("ResidualMatrix subsetting works as expected", {
    expect_equal_and_resmat <- function(x, y) {
        expect_s4_class(x, "ResidualMatrix") # class is correctly preserved by direct seed modification.
        expect_equal(as.matrix(x), y)
    }

    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$res))
        j <- sample(ncol(test$res))
        expect_equal_and_resmat(test$res[i,], test$ref[i,])
        expect_equal_and_resmat(test$res[,j], test$ref[,j])
        expect_equal_and_resmat(test$res[i,j], test$ref[i,j])

        # Works with zero dimensions.
        expect_equal_and_resmat(test$res[0,], test$ref[0,])
        expect_equal_and_resmat(test$res[,0], test$ref[,0])
        expect_equal_and_resmat(test$res[0,0], test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_equal(test$res[i[1],], test$ref[i[1],])
        expect_equal(test$res[,j[1]], test$ref[,j[1]])
        expect_equal_and_resmat(test$res[i[1],drop=FALSE], test$ref[i[1],,drop=FALSE])
        expect_equal_and_resmat(test$res[,j[1],drop=FALSE], test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$res)
        expect_equal(t(alt[,i]), test$res[i,])
        expect_equal(t(alt[j,]), test$res[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$res)))
        colnames(test$res) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_equal_and_resmat(test$res[,ch], test$ref[,ch])
    }
})

set.seed(1000012)
test_that("ResidualMatrix centering is preserved or lost correctly", {
    NR <- 10
    NC <- 15
    for (it in 1:4) {
        design <- generate_design(NR, it)
        y <- ResidualMatrix(matrix(rnorm(NR*NC), ncol=NC), design)

        expect_true(BiocSingular:::is_centered(DelayedArray::seed(y)))
        expect_false(BiocSingular:::is_centered(DelayedArray::seed(y[1:5,])))
        expect_true(BiocSingular:::is_centered(DelayedArray::seed(y[,1:5])))
        expect_true(BiocSingular:::is_centered(DelayedArray::seed(t(y))))
    }
})

##########################
# Matrix multiplication.

test_that("ResidualMatrix right multiplication works as expected", {
    possibles <- spawn_scenarios(100, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix left multiplication works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix dual multiplication works as expected", {
    possibles1 <- spawn_scenarios(10, 20)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(test1$res %*% test2$res, test1$ref %*% test2$ref)

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(test1$res[0,] %*% test2$res, test1$ref[0,] %*% test2$ref)
            expect_equal_product(test1$res %*% test2$res[,0], test1$ref %*% test2$ref[,0])
            expect_equal_product(test1$res[,0] %*% test2$res[0,], test1$ref[,0] %*% test2$ref[0,])
            expect_equal_product(test1$res[0,] %*% test2$res[,0], test1$ref[0,] %*% test2$ref[,0])
        }
    }
})

##########################

test_that("ResidualMatrix lonely crossproduct works as expected", {
    possibles <- spawn_scenarios(90, 30)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res
        expect_equal_product(crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("ResidualMatrix crossproduct from right works as expected", {
    possibles <- spawn_scenarios(60, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix crossproduct from left works as expected", {
    possibles <- spawn_scenarios(40, 100)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix dual crossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(crossprod(test1$res, test2$res), crossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(crossprod(test1$res[,0], test2$res), crossprod(test1$ref[,0], test2$ref))
            expect_equal_product(crossprod(test1$res, test2$res[,0]), crossprod(test1$ref, test2$ref[,0]))
            expect_equal_product(crossprod(test1$res[0,], test2$res[0,]), crossprod(test1$ref[0,], test2$ref[0,]))
        }
    }
})

##########################

test_that("ResidualMatrix lonely tcrossproduct works as expected", {
    possibles <- spawn_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res
        expect_equal_product(tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("ResidualMatrix tcrossproduct from right works as expected", {
    possibles <- spawn_scenarios(60, 70)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix tcrossproduct from left works as expected", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$res

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

test_that("ResidualMatrix dual tcrossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(25, 50)
        for (test2 in possibles2) {

            expect_equal_product(tcrossprod(test1$res, test2$res), tcrossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(tcrossprod(test1$res[0,], test2$res), tcrossprod(test1$ref[0,], test2$ref))
            expect_equal_product(tcrossprod(test1$res, test2$res[0,]), tcrossprod(test1$ref, test2$ref[0,]))
            expect_equal_product(tcrossprod(test1$res[,0], test2$res[,0]), tcrossprod(test1$ref[,0], test2$ref[,0]))
        }
    }
})
