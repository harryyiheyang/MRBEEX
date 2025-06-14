// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// construct_sparse_blockwise_LD_cpp
List construct_sparse_blockwise_LD_cpp(const arma::sp_mat& LD, const IntegerVector& cluster_index, const IntegerVector& cluster_sampling, double admm_rho);
RcppExport SEXP _MRBEEX_construct_sparse_blockwise_LD_cpp(SEXP LDSEXP, SEXP cluster_indexSEXP, SEXP cluster_samplingSEXP, SEXP admm_rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type LD(LDSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cluster_index(cluster_indexSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cluster_sampling(cluster_samplingSEXP);
    Rcpp::traits::input_parameter< double >::type admm_rho(admm_rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(construct_sparse_blockwise_LD_cpp(LD, cluster_index, cluster_sampling, admm_rho));
    return rcpp_result_gen;
END_RCPP
}
// get_matrix_info_cpp
List get_matrix_info_cpp(const arma::sp_mat& mat);
RcppExport SEXP _MRBEEX_get_matrix_info_cpp(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(get_matrix_info_cpp(mat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRBEEX_construct_sparse_blockwise_LD_cpp", (DL_FUNC) &_MRBEEX_construct_sparse_blockwise_LD_cpp, 4},
    {"_MRBEEX_get_matrix_info_cpp", (DL_FUNC) &_MRBEEX_get_matrix_info_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRBEEX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
