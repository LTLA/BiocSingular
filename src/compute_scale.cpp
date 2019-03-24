#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/utils/const_column.h"
#include <cmath>

template<class M>
Rcpp::NumericVector compute_scale_internal(Rcpp::RObject mat, Rcpp::RObject centering) {
    auto ptr=beachmat::create_matrix<M>(mat);
    size_t ncols=ptr->get_ncol();
    size_t nrows=ptr->get_nrow();

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
    beachmat::const_column<M> col_holder(ptr.get());

    for (size_t i=0; i<ncols; ++i) {
        col_holder.fill(i);
        auto n=col_holder.get_n();
        auto vals=col_holder.get_values();

        double& current=output[i];
        for (size_t j=0; j<n; ++j, ++vals) {
            double val=*vals;
            if (do_center) {
                val-=numeric_centers[i];
            }
            current+=val*val;
        }

        if (do_center && col_holder.is_sparse()) { // adding the contribution of the zeroes.
            current += (nrows - n) * (numeric_centers[i] * numeric_centers[i]);
        }

        current/=nrows-1;
        current=std::sqrt(current);
    }
    
    return output;
}

//' @useDynLib BiocSingular
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_scale(Rcpp::RObject mat, Rcpp::RObject centering) {
    auto rtype=beachmat::find_sexp_type(mat);
    if (rtype==INTSXP) {
        return compute_scale_internal<beachmat::integer_matrix>(mat, centering);
    } else {
        return compute_scale_internal<beachmat::numeric_matrix>(mat, centering);
    }
}
