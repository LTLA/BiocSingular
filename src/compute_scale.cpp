#include "Rcpp.h"
#include "beachmat3/beachmat.h"
#include <cmath>

//' @useDynLib BiocSingular
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_scale(Rcpp::RObject mat, Rcpp::RObject centering) {
    auto ptr = beachmat::read_lin_block(mat);
    size_t ncols = ptr->get_ncol();
    size_t nrows = ptr->get_nrow();

    if (nrows <= 1) {
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

    if (ptr->is_sparse()) {
        auto sptr = beachmat::promote_to_sparse(ptr);
        std::vector<double> work_x(nrows);
        std::vector<int> work_i(nrows);

        for (size_t i=0; i<ncols; ++i) {
            auto idx = sptr->get_col(i, work_x.data(), work_i.data());

            double& current=output[i];
            for (size_t j = 0; j < idx.n; ++j, ++idx.x) {
                double val = *(idx.x);
                if (do_center) {
                    val -= numeric_centers[i];
                }
                current += val*val;
            }

            // Adding back the contribution of the zeros.
            if (do_center) {
                current += (nrows - idx.n) * (numeric_centers[i] * numeric_centers[i]);
            }
        }

    } else {
        std::vector<double> workspace(nrows);
        for (size_t i = 0; i < ncols; ++i) {
            auto vals = ptr->get_col(i, workspace.data());
            double& current = output[i];

            for (size_t j=0; j < nrows; ++j, ++vals) {
                double val=*vals;
                if (do_center) {
                    val-=numeric_centers[i];
                }
                current += val*val;
            }
        }
    }

    for (auto& current : output) {
        current /= nrows-1;
        current = std::sqrt(current);
    }
    
    return output;
}
