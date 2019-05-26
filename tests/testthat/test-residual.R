# Tests the ResidualMatrix implementation.
# library(testthat); library(BiocSingular); source("test-residual.R")

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

            # Run through a host of different design matrices.
            if (it==1L) {
                design <- NULL
            } else if (it==2L) {
                cov <- rnorm(nrow(y))
                design <- model.matrix(~cov)
            } else if (it==3L) {
                g <- factor(rep(1:3, length.out=nrow(y)))
                design <- model.matrix(~0 + g)
            } else if (it==4L) {
                cov <- rnorm(nrow(y))
                g <- factor(rep(1:2, length.out=nrow(y)))
                design <- model.matrix(~0 + g + cov)
            } else if (it==5L) {
                design <- cbind(rnorm(nrow(y)))
            }

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


