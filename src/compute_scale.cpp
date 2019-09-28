#include "Rcpp.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"

#include <cmath>
#include <omp.h>

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

    int nth=omp_get_max_threads();
    Rcpp::NumericVector battery(nrows*nth);
    Rprintf("%i\n", nth);

    #pragma omp parallel 
    {
        auto holder=battery.begin() + omp_get_thread_num() * nrows;
        decltype(ptr) pptr=NULL;

        #pragma omp critical
        if (nth==1) {
            pptr=std::move(ptr);
        } else {
            pptr=ptr->clone();
        }

        #pragma omp for schedule(static)
        for (size_t i=0; i<ncols; ++i) {
            #pragma omp critical
            {
                pptr->get_col(i, holder);
            }

            double& current=output[i];
            auto vals=holder;
            for (size_t j=0; j<nrows; ++j, ++vals) {
                double val=*vals;
                if (do_center) {
                    val-=numeric_centers[i];
                }
                current+=val*val;
            }

            current/=nrows-1;
            current=std::sqrt(current);
        }
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

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector set_omp_threads(Rcpp::IntegerVector nthreads)
{
    if (nthreads.size()!=1L) {
        Rf_error("'nthreads' must be integer(1)");
    }
    omp_set_num_threads(nthreads[0]);
    return nthreads;
}
