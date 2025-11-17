#include "Rtatami.h"
#include "tatami_stats/tatami_stats.hpp"
#include <cmath>
#include <algorithm>
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
SEXP set_executor(SEXP ptr) {
    Rtatami::set_executor(ptr);
    return R_NilValue;
}

//' @useDynLib BiocSingular
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_center(Rcpp::RObject mat, int nthreads) {
    Rtatami::BoundNumericPointer bound(mat);
    const auto& ptr = bound->ptr;
    Rcpp::NumericVector output(ptr->ncol());
    double NR = ptr->nrow();

    if (NR == 0) {
        std::fill(output.begin(), output.end(), R_NaReal);
    } else {
        tatami_stats::sums::Options opt;
        opt.num_threads = nthreads;
        auto row_sums = tatami_stats::sums::by_column(ptr.get(), opt);
        for (int c = 0, cend = ptr->ncol(); c < cend; ++c) {
            output[c] = row_sums[c] / NR;
        }
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::List compute_center_and_scale(Rcpp::RObject mat, int nthreads) {
    Rtatami::BoundNumericPointer bound(mat);
    const auto& ptr = bound->ptr;
    auto NR = ptr->nrow();
    auto NC = ptr->ncol();

    Rcpp::NumericVector center(NC), scale(NC);
    double* cptr = static_cast<double*>(center.begin());
    double* sptr = static_cast<double*>(scale.begin());

    // Handling edge cases.
    if (NR <= 1) {
        if (NR == 0) {
            std::fill(center.begin(), center.end(), R_NaReal);
        } else {
            auto iptr = ptr->dense_row()->fetch(0, cptr);
            tatami::copy_n(iptr, NC, cptr);
        }
        std::fill(scale.begin(), scale.end(), R_NaReal);
        return Rcpp::List::create(
            Rcpp::Named("center") = center, 
            Rcpp::Named("scale") = scale
        );
    }

    if (ptr->prefer_rows()) {
        if (ptr->sparse()) {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<true>(ptr.get(), true, 0, NR, start, len);
                std::vector<double> vbuffer(len);
                std::vector<int> ibuffer(len);

                std::vector<double> tmp_means(len), tmp_vars(len);
                tatami_stats::variances::RunningSparse<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false, start);
                for (int r = 0; r < NR; ++r) {
                    auto range = ext->fetch(r, vbuffer.data(), ibuffer.data());
                    runner.add(range.value, range.index, range.number);
                }
                runner.finish();

                std::copy(tmp_means.begin(), tmp_means.end(), cptr + start);
                for (auto& v : tmp_vars) {
                    v = std::sqrt(v);
                }
                std::copy(tmp_vars.begin(), tmp_vars.end(), sptr + start);
            }, NC, nthreads);

        } else {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<false>(ptr.get(), true, 0, NR, start, len);
                std::vector<double> buffer(len);

                std::vector<double> tmp_means(len), tmp_vars(len);
                tatami_stats::variances::RunningDense<double, double, int> runner(len, tmp_means.data(), tmp_vars.data(), false);
                for (int r = 0; r < NR; ++r) {
                    auto ptr = ext->fetch(r, buffer.data());
                    runner.add(ptr);
                }
                runner.finish();

                std::copy(tmp_means.begin(), tmp_means.end(), cptr + start);
                for (auto& v : tmp_vars) {
                    v = std::sqrt(v);
                }
                std::copy(tmp_vars.begin(), tmp_vars.end(), sptr + start);
            }, NC, nthreads);
        }

    } else {
        if (ptr->sparse()) {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                tatami::Options opt;
                opt.sparse_extract_index = false;
                auto ext = tatami::consecutive_extractor<true>(ptr.get(), false, start, len, opt);
                std::vector<double> vbuffer(NR);
                for (int c = start, end = start + len; c < end; ++c) {
                    auto range = ext->fetch(c, vbuffer.data(), NULL);
                    auto paired = tatami_stats::variances::direct(range.value, range.number, NR, false);
                    cptr[c] = paired.first;
                    sptr[c] = std::sqrt(paired.second);
                }
            }, NC, nthreads);

        } else {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<false>(ptr.get(), false, start, len);
                std::vector<double> buffer(NR);
                for (int c = start, end = start + len; c < end; ++c) {
                    auto ptr = ext->fetch(c, buffer.data());
                    auto paired = tatami_stats::variances::direct(ptr, NR, false);
                    cptr[c] = paired.first;
                    sptr[c] = std::sqrt(paired.second);
                }
            }, NC, nthreads);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("center") = center, 
        Rcpp::Named("scale") = scale
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector compute_scale(Rcpp::RObject mat, Rcpp::NumericVector centers, int nthreads) {
    Rtatami::BoundNumericPointer bound(mat);
    const auto& ptr = bound->ptr;
    auto NR = ptr->nrow();
    auto NC = ptr->ncol();

    if (NC != static_cast<decltype(NC)>(centers.size())) {
        throw std::runtime_error("'centers' should be equal to the number of columns in 'mat'");
    }

    Rcpp::NumericVector output(NC);
    double* optr = static_cast<double*>(output.begin());
    const double* cptr = static_cast<const double*>(centers.begin());

    // Handling edge cases.
    if (NR <= 1) {
        std::fill(output.begin(), output.end(), R_NaReal);
        return output;
    }

    if (ptr->prefer_rows()) {
        if (ptr->sparse()) {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<true>(ptr.get(), true, 0, NR, start, len);

                std::vector<double> vbuffer(len);
                std::vector<int> ibuffer(len);
                std::vector<double> tmp_vars(len);
                std::vector<int> nonzeros(len);

                for (int r = 0; r < NR; ++r) {
                    auto range = ext->fetch(r, vbuffer.data(), ibuffer.data());
                    for (int i = 0; i < range.number; ++i) {
                        double diff = range.value[i] - cptr[range.index[i]];
                        auto offset = range.index[i] - start;
                        tmp_vars[offset] += diff * diff;
                        ++nonzeros[offset];
                    }
                }

                for (int i = 0; i < len; ++i) {
                    double center = cptr[i + start];
                    double v = tmp_vars[i] + (NR - nonzeros[i]) * center * center;
                    optr[start + i] = std::sqrt(v / static_cast<double>(NR - 1));
                }
            }, NC, nthreads);

        } else {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<false>(ptr.get(), true, 0, NR, start, len);

                std::vector<double> buffer(len);
                std::vector<double> tmp_vars(len);
                std::vector<int> nonzeros(len);

                for (int r = 0; r < NR; ++r) {
                    auto ptr = ext->fetch(r, buffer.data());
                    for (int i = 0; i < len; ++i) {
                        double diff = ptr[i] - cptr[i + start];
                        tmp_vars[i] += diff * diff;
                    }
                }

                for (auto& v : tmp_vars) {
                    v = std::sqrt(v / static_cast<double>(NR - 1));
                }
                std::copy(tmp_vars.begin(), tmp_vars.end(), optr + start);
            }, NC, nthreads);
        }

    } else {
        if (ptr->sparse()) {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                tatami::Options opt;
                opt.sparse_extract_index = false;
                auto ext = tatami::consecutive_extractor<true>(ptr.get(), false, start, len, opt);
                std::vector<double> vbuffer(NR);

                for (int c = start, end = start + len; c < end; ++c) {
                    auto range = ext->fetch(c, vbuffer.data(), NULL);
                    double center = cptr[c];

                    double tmp = 0;
                    for (int i = 0; i < range.number; ++i) {
                        double diff = range.value[i] - center;
                        tmp += diff * diff;
                    }

                    tmp += (NR - range.number) * center * center;
                    optr[c] = std::sqrt(tmp / static_cast<double>(NR - 1));
                }
            }, NC, nthreads);

        } else {
            tatami::parallelize([&](size_t, int start, int len) -> void {
                auto ext = tatami::consecutive_extractor<false>(ptr.get(), false, start, len);
                std::vector<double> buffer(NR);
                for (int c = start, end = start + len; c < end; ++c) {
                    auto ptr = ext->fetch(c, buffer.data());
                    double center = cptr[c];

                    double tmp = 0;
                    for (int r = 0; r < NR; ++r) {
                        double diff = ptr[r] - center;
                        tmp += diff * diff;
                    }
                    optr[c] = std::sqrt(tmp / static_cast<double>(NR - 1));
                }
            }, NC, nthreads);
        }
    }

    return output;
}
