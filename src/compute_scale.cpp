#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include <cmath>

//' @useDynLib BiocSingular
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_scale(Rcpp::RObject mat, Rcpp::RObject centering) {
    auto dptr=beachmat::create_numeric_matrix(mat);
    size_t ncols=dptr->get_ncol();
    size_t nrows=dptr->get_nrow();

    if (nrows<=1) {
        return Rcpp::NumericVector(ncols, R_NaReal);
    }

    const bool do_center=!centering.isNULL();
    Rcpp::NumericVector numeric_centers;
    if (do_center) {
        numeric_centers=Rcpp::NumericVector(centering);
        if (numeric_centers.size()!=ncols) {
            throw std::runtime_error("length of centering vector should be equal to number of columns in 'mat'");
        }
    }

    Rcpp::NumericVector output(ncols);
    Rcpp::NumericVector tmp(nrows);

    for (size_t i=0; i<ncols; ++i) {
        dptr->get_col(i, tmp.begin());

        double& current=output[i];
        for (size_t j=0; j<nrows; ++j) {
            double val=tmp[j];
            if (do_center) {
                val-=numeric_centers[i];
            }
            current+=val*val;
        }

        current/=nrows-1;
        current=std::sqrt(current);
    }
    
    return output;
}
