# Checks out a LinearEmbeddingMatrix.
# library(BiocSingular); library(testthat); source("test-lem.R")

set.seed(1000)
ncells <- 100

library(S4Vectors)
factors <- matrix(rnorm(1000), ncol=10)
loadings <- matrix(runif(10000), ncol=10)
lem <- LinearEmbeddingMatrix(factors, loadings)

test_that("LEM construction works correctly", {
    expect_equivalent(sampleFactors(lem, withDimnames=FALSE), factors)
    expect_identical(featureLoadings(lem), loadings)
    expect_identical(nrow(factorData(lem)), ncol(lem))
    expect_identical(ncol(factorData(lem)), 0L)

    # Trying again with an initialized factorData.
    fd <- DataFrame(YAY=seq_len(ncol(factors)))
    lem2 <- LinearEmbeddingMatrix(factors, loadings, factorData=fd)
    expect_identical(factorData(lem2), fd)
    expect_identical(lem2$YAY, fd$YAY)

    # Should throw errors if it doesn't make sense.
    expect_error(LinearEmbeddingMatrix(factors, loadings[,1:2]), "must have the same number of\n *columns")
    expect_error(LinearEmbeddingMatrix(factors, loadings, DataFrame(YAY=2)), "one row per factor")

    # Add metadata.
    md <- list(a = 1, b = "two")
    lem3 <- LinearEmbeddingMatrix(factors, loadings, factorData=fd, metadata=md)
    expect_identical(md, metadata(lem3))
})

test_that("LEM setters work correctly", {
    lem2 <- lem
    sampleFactors(lem2) <- factors * 2
    expect_identical(sampleFactors(lem2, withDimnames=FALSE), factors * 2)

    lem2 <- lem
    featureLoadings(lem2) <- -loadings
    expect_identical(featureLoadings(lem2), -loadings)

    to.add <- LETTERS[seq_len(ncol(factors))]
    lem2 <- lem
    factorData(lem2)$whee <- to.add
    expect_identical(factorData(lem2)$whee, to.add)

    # Factor Data settters.
    lem2$whee <- NULL
    expect_identical(lem2$whee, NULL)
    lem2$whee <- 42
    expect_identical(lem2$whee, rep(42, ncol(lem2)))

    # Dimnames setters.
    lem2 <- lem
    colnames(lem2) <- to.add
    expect_identical(colnames(lem2), to.add)
    cell.names <- paste0("Cell", seq_len(nrow(factors)))
    rownames(lem2) <- cell.names
    expect_identical(rownames(lem2), cell.names)

    # Should throw errors if it doesn't make sense.
    expect_error(sampleFactors(lem) <- factors[,1], "must have the same number of.*columns")
    expect_error(featureLoadings(lem) <- loadings[,1], "must have the same number of.*columns")
    expect_error(factorData(lem) <- 2, "DataFrame")
    expect_error(factorData(lem) <- DataFrame(YAY=2), "one row per factor")
})

test_that("getting and setting with names makes sense", {
    rn <- paste0("Cell_", seq_len(nrow(lem)))
    rownames(lem) <- rn
    cn <- paste0("Factor_", seq_len(ncol(lem)))
    colnames(lem) <- cn

    expect_identical(rn, rownames(lem))
    expect_identical(cn, colnames(lem))

    expect_identical(rownames(sampleFactors(lem, withDimnames=FALSE)), rn)
    expect_identical(colnames(sampleFactors(lem, withDimnames=FALSE)), NULL)
    expect_identical(rownames(sampleFactors(lem)), rn)
    expect_identical(colnames(sampleFactors(lem)), cn)

    expect_identical(colnames(featureLoadings(lem, withDimnames=FALSE)), NULL)
    expect_identical(colnames(featureLoadings(lem)), cn)

    expect_identical(rownames(factorData(lem)), cn) # Doesn't make a difference, as this is the reference location.

    # This will delete the column names.
    lem2 <- lem
    factorData(lem2) <- new("DFrame", nrows=ncol(lem))
    expect_equal(ncol(factorData(lem2)), 0L)
    expect_identical(rownames(factorData(lem2)), NULL)
})

library(Matrix)
test_that("as.matrix works as expected", {
    expect_equal(as.matrix(lem), factors)

    alt <- rsparsematrix(nrow(lem), ncol(lem), density=0.1)
    sampleFactors(lem) <- alt
    expect_identical(sampleFactors(lem, withDimnames=FALSE), alt)
    expect_equal(as.matrix(lem), as.matrix(alt))
})

test_that("manipulation of metadata is correct", {
    metadata(lem)$yay <- 1
    expect_equal(metadata(lem)$yay, 1)
    metadata(lem)$yay <- "stuff"
    expect_identical(metadata(lem)$yay, "stuff")
})

#############################
#############################

test_that("rbind works correctly", {
    shuffled <- sample(nrow(lem))
    lem.alt <- lem[shuffled,]
    samish <- rbind(lem.alt)
    expect_identical(sampleFactors(lem.alt), sampleFactors(samish))
    expect_identical(featureLoadings(lem.alt), featureLoadings(samish))
    expect_identical(factorData(lem.alt), factorData(samish))

    lem2 <- rbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2, withDimnames=FALSE),
        rbind(sampleFactors(lem, withDimnames=FALSE), sampleFactors(lem.alt, withDimnames=FALSE)))
    expect_identical(featureLoadings(lem2), featureLoadings(lem))
    expect_identical(factorData(lem2), factorData(lem))
    expect_identical(rownames(lem2), NULL)
    expect_identical(colnames(lem2), NULL)

    # Works correctly with names.
    unnamed <- lem
    rownames(lem) <- paste0("CELL", seq_len(nrow(lem)))
    colnames(lem) <- paste0("FACTOR", seq_len(ncol(lem)))

    lem3 <- rbind(lem, lem[shuffled,])
    expect_identical(rownames(lem3), c(rownames(lem), rownames(lem)[shuffled]))
    expect_identical(colnames(lem3), colnames(lem))
    expect_equivalent(sampleFactors(lem3), sampleFactors(lem2))

    # Throws errors correctly.
    lem.alt <- lem[,1:2]
    expect_error(rbind(lem, lem.alt), "number of columns of matrices must match")

    lem.alt <- lem
    featureLoadings(lem.alt) <- featureLoadings(lem.alt) + 1
    expect_error(rbind(lem, lem.alt), "not identical")

    lem.alt <- lem
    lem.alt$X <- 1
    expect_error(rbind(lem, lem.alt), "not identical")
})

test_that("cbind works correctly", {
    shuffled <- sample(ncol(factors))
    lem.alt <- lem[,shuffled]
    samish <- cbind(lem.alt)
    expect_identical(sampleFactors(lem.alt), sampleFactors(samish))
    expect_identical(featureLoadings(lem.alt), featureLoadings(samish))
    expect_identical(factorData(lem.alt), factorData(samish))

    lem2 <- cbind(lem, lem.alt)
    expect_identical(sampleFactors(lem2, withDimnames=FALSE),
        cbind(sampleFactors(lem, withDimnames=FALSE), sampleFactors(lem.alt, withDimnames=FALSE)))
    expect_identical(featureLoadings(lem2), cbind(featureLoadings(lem), featureLoadings(lem.alt)))
    expect_identical(factorData(lem2), rbind(factorData(lem), factorData(lem.alt)))
    expect_identical(rownames(lem2), NULL)
    expect_identical(colnames(lem2), NULL)

    # Works correctly with names.
    unnamed <- lem
    rownames(lem) <- paste0("CELL", seq_len(nrow(lem)))
    colnames(lem) <- paste0("FACTOR", seq_len(ncol(lem)))

    lem3<- cbind(lem, lem[,shuffled])
    expect_identical(rownames(lem3), rownames(lem))
    expect_identical(colnames(lem3), c(colnames(lem), colnames(lem)[shuffled]))
    expect_equivalent(sampleFactors(lem3), sampleFactors(lem2))

    # Throws errors correctly.
    lem.alt <- lem[1:5,]
    expect_error(cbind(lem, lem.alt), "number of rows")

    lem.alt <- lem
    factorData(lem.alt)$WHEE <- 1
    expect_error(cbind(lem, lem.alt), "number of columns")
})
    
#############################
#############################

factorData(lem) <- fdata <- DataFrame(WHEE=sample(LETTERS, ncol(factors)))

test_that("subsetting works correctly for different index types", {
    keep.dimnames <- FALSE
    for (x in 1:3) {
        if (x==1) {
            by.row <- sample(nrow(factors), nrow(factors)/2)
            by.col <- sample(ncol(factors), ncol(factors)/2)
        } else if (x==2) {
            by.row <- rbinom(nrow(factors), 1, 0.5)==1
            by.col <- rbinom(ncol(factors), 1, 0.5)==1
        } else {
            colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
            rownames(lem) <- paste0("Cell_", seq_len(nrow(factors)))
            dimnames(factors) <- dimnames(lem)
            colnames(loadings) <- rownames(fdata) <- colnames(lem)

            by.row <- sample(rownames(lem), nrow(factors)/2)
            by.col <- sample(colnames(lem), ncol(factors)/2)
            keep.dimnames <- TRUE
        }

        # By row.
        lem.alt <- lem[by.row,]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[by.row,])
        expect_identical(featureLoadings(lem.alt), loadings)
        expect_identical(factorData(lem.alt), fdata)

        # By column.
        lem.alt <- lem[,by.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[,by.col])
        expect_identical(featureLoadings(lem.alt), loadings[,by.col])
        expect_identical(factorData(lem.alt), fdata[by.col,,drop=FALSE])

        # By row and column.
        lem.alt <- lem[by.row, by.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), factors[by.row,by.col])
        expect_identical(featureLoadings(lem.alt), loadings[,by.col])
        expect_identical(factorData(lem.alt), fdata[by.col,,drop=FALSE])
    }
})

test_that("subsetting works correctly with drop=TRUE", {
    # By row, with and without drop.
    keeper <- lem[1,]
    expect_identical(keeper, factors[1,])

    nodrop <- lem[1,,drop=FALSE]
    expect_identical(sampleFactors(nodrop, withDimnames=FALSE), factors[1,,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings)
    expect_identical(factorData(nodrop), fdata)

    # By column, with and without drop.
    keeper <- lem[,1]
    expect_identical(keeper, factors[,1])

    nodrop <- lem[,1,drop=FALSE]
    expect_identical(sampleFactors(nodrop, withDimnames=FALSE), factors[,1,drop=FALSE])
    expect_identical(featureLoadings(nodrop), loadings[,1,drop=FALSE])
    expect_identical(factorData(nodrop), fdata[1,,drop=FALSE])

    # Handles names correctly.
    colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
    rownames(lem) <- paste0("Cell_", seq_len(nrow(factors)))
    expect_identical(lem[,1], sampleFactors(lem)[,1])
    expect_identical(lem[2,], sampleFactors(lem)[2,])
    expect_identical(lem[1,2], sampleFactors(lem)[1,2])

    # Throws errors correctly.
    expect_error(lem[nrow(lem)+1,], "subscript out of bounds", fixed=TRUE)
    expect_error(lem["A",], "subscript out of bounds", fixed=TRUE)
})

test_that("subsetting assignment works correctly", {
    keep.dimnames <- FALSE

    for (x in 1:3) {
        if (x==1) {
            dest.row <- sample(nrow(factors), nrow(factors)/2)
            src.row <- sample(nrow(factors), nrow(factors)/2)
            dest.col <- sample(ncol(factors), ncol(factors)/2)
            src.col <- sample(ncol(factors), ncol(factors)/2)
        } else if (x==2) {
            dest.row <- seq_len(nrow(factors)) %in% sample(nrow(factors), nrow(factors)/2)
            dest.col <- seq_len(ncol(factors)) %in% sample(ncol(factors), ncol(factors)/2)
            src.row <- sample(dest.row)
            src.col <- sample(dest.col)
        } else {
            colnames(lem) <- paste0("Factor_", seq_len(ncol(factors)))
            rownames(lem) <- paste0("Gene_", seq_len(nrow(factors)))
            dimnames(factors) <- dimnames(lem)
            colnames(loadings) <- rownames(fdata) <- colnames(lem)

            dest.row <- sample(rownames(lem), nrow(factors)/2)
            src.row <- sample(rownames(lem), nrow(factors)/2)
            dest.col <- sample(colnames(lem), ncol(factors)/2)
            src.col <- sample(colnames(lem), ncol(factors)/2)
            keep.dimnames <- TRUE
        }

        # By row.
        lem.alt <- lem
        lem.alt[dest.row,] <- lem[src.row,]

        ref <- factors
        ref[dest.row,] <- ref[src.row,]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref)
        expect_identical(featureLoadings(lem.alt), loadings)
        expect_identical(factorData(lem.alt), fdata)

        # By column.
        lem.alt <- lem
        lem.alt[,dest.col] <- lem[,src.col]

        ref_sf <- factors
        ref_sf[,dest.col] <- ref_sf[,src.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref_sf)

        ref_fl <- loadings
        ref_fl[,dest.col] <- ref_fl[,src.col]
        expect_identical(featureLoadings(lem.alt), ref_fl)

        ref_fd <- fdata
        ref_fd[dest.col,] <- fdata[src.col,]
        expect_identical(factorData(lem.alt), ref_fd)

        # By row and column.
        lem.alt <- lem
        lem.alt[dest.row,dest.col] <- lem[src.row,src.col]

        ref_sf <- factors
        ref_sf[dest.row,dest.col] <- ref_sf[src.row,src.col]
        expect_identical(sampleFactors(lem.alt, withDimnames=keep.dimnames), ref_sf)

        ref_fl <- loadings
        ref_fl[,dest.col] <- ref_fl[,src.col]
        expect_identical(featureLoadings(lem.alt), ref_fl)

        ref_fd <- fdata
        ref_fd[dest.col,] <- fdata[src.col,]
        expect_identical(factorData(lem.alt), ref_fd)
    }
})
