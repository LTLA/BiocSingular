# Tests the LowRankMatrix implementation.
# library(testthat); library(BiocSingular); source("test-lowrank.R")

spawn_scenarios <- function(NR=50, NC=20, NP=10) {
    collected <- list()

    for (i in 1:2) {
        if (i==1L) {
            comp <- matrix(rnorm(NC*NP), ncol=NP)
        } else {
            comp <- Matrix::rsparsematrix(NC, NP, 0.1)
        }
        
        for (j in 1:2) {
            if (j==1L) {
                rot <- matrix(rnorm(NR*NP), ncol=NP)
            } else {
                rot <- Matrix::rsparsematrix(NR, NP, 0.1)
            }

            ref <- as.matrix(rot %*% t(comp))
            lrm <- LowRankMatrix(rot, comp)
            collected <- c(collected, list(list(ref=ref, lrm=lrm)))
        }
    }

    collected
}

##########################

expect_identical_unnamed <- function(x, y) {
    if (all(lengths(dimnames(x))==0L)) dimnames(x) <- NULL
    if (all(lengths(dimnames(y))==0L)) dimnames(y) <- NULL
    expect_identical(x, y)
}

set.seed(500001)
test_that("LowRankMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_s4_class(test$lrm, "LowRankMatrix")
        expect_identical(test$lrm, LowRankMatrix(DelayedArray::seed(test$lrm)))

        expect_identical(dim(test$lrm), dim(test$ref))
        expect_identical(extract_array(test$lrm, list(1:10, 1:10)), test$ref[1:10, 1:10])
        expect_identical(extract_array(test$lrm, list(1:10, NULL)), test$ref[1:10,])
        expect_identical(extract_array(test$lrm, list(NULL, 1:10)), test$ref[,1:10])
        expect_identical_unnamed(as.matrix(test$lrm), test$ref)

        ttest <- t(test$lrm)
        expect_s4_class(ttest, "LowRankMatrix") # still a LRMat!
        expect_identical(t(ttest), test$lrm)
        expect_identical_unnamed(as.matrix(ttest), t(test$ref))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$lrm)))
        colnames(test$lrm) <- spawn_names
        expect_identical(spawn_names, colnames(test$lrm))
        expect_s4_class(test$lrm, "LowRankMatrix") # still a LRMat!
    }

    # Checking erronious inputs.
    expect_error(LowRankMatrix(1, 1), "must be matrix-like")
    expect_error(LowRankMatrix(cbind(1:5), cbind(1:5, 2:6)), "must be the same")
})

set.seed(500002)
test_that("LowRankMatrix subsetting works as expected", {
    expect_identical_and_lrmmat <- function(x, y) {
        expect_s4_class(x, "LowRankMatrix") # class is correctly preserved by direct seed modification.
        expect_identical_unnamed(as.matrix(x), y)
    }

    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$lrm))
        j <- sample(ncol(test$lrm))
        expect_identical_and_lrmmat(test$lrm[i,], test$ref[i,])
        expect_identical_and_lrmmat(test$lrm[,j], test$ref[,j])
        expect_identical_and_lrmmat(test$lrm[i,j], test$ref[i,j])
        
        # Works with zero dimensions.
        expect_identical_and_lrmmat(test$lrm[0,], test$ref[0,])
        expect_identical_and_lrmmat(test$lrm[,0], test$ref[,0])
        expect_identical_and_lrmmat(test$lrm[0,0], test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_identical(test$lrm[i[1],], test$ref[i[1],])
        expect_identical(test$lrm[,j[1]], test$ref[,j[1]])
        expect_identical_and_lrmmat(test$lrm[i[1],drop=FALSE], test$ref[i[1],,drop=FALSE])
        expect_identical_and_lrmmat(test$lrm[,j[1],drop=FALSE], test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$lrm)
        expect_identical(t(alt[,i]), test$lrm[i,])
        expect_identical(t(alt[j,]), test$lrm[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$lrm)))
        colnames(test$lrm) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_identical_and_lrmmat(test$lrm[,ch], test$ref[,ch])
    }
})

##########################

library(DelayedArray)
test_that("DelayedMatrix wrapping works", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        expect_identical_unnamed(as.matrix(test$lrm+1), test$ref+1)

        v <- rnorm(nrow(test$lrm))
        expect_identical_unnamed(as.matrix(test$lrm+v), test$ref+v)
        expect_identical_unnamed(as.matrix(test$lrm*v), test$ref*v)

        w <- rnorm(ncol(test$lrm))
        expect_identical_unnamed(as.matrix(sweep(test$lrm, 2, w, "*")), sweep(test$ref, 2, w, "*"))
    }
})
