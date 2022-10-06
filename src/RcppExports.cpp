// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// count_omp_thread
int count_omp_thread();
RcppExport SEXP _bcfac_count_omp_thread() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(count_omp_thread());
    return rcpp_result_gen;
END_RCPP
}
// fit_bcf_hetero
void fit_bcf_hetero(NumericVector& Y1, NumericVector& Y0, NumericMatrix& var_count, NumericVector& var_prob, double& sigma2_exp, NumericVector& sigma2_out_hist, NumericVector& dir_alpha_hist, const NumericVector& Y, const NumericVector& trt, const NumericMatrix& X, const double trt_treated, const double trt_control, const int chain_idx, const int num_chain, const int total_iter, const int num_burn_in, const int num_thin, const int num_post_sample, const int num_tree, const NumericVector& step_prob, const double alpha, const double beta, const double nu, const double alpha2, const double beta2, const double nu2, const int num_tree_mod, const double lambda_exp, const double lambda_out, const double lambda_mod, const int boot_size, const bool is_binary_trt, const bool parallel, const bool verbose);
RcppExport SEXP _bcfac_fit_bcf_hetero(SEXP Y1SEXP, SEXP Y0SEXP, SEXP var_countSEXP, SEXP var_probSEXP, SEXP sigma2_expSEXP, SEXP sigma2_out_histSEXP, SEXP dir_alpha_histSEXP, SEXP YSEXP, SEXP trtSEXP, SEXP XSEXP, SEXP trt_treatedSEXP, SEXP trt_controlSEXP, SEXP chain_idxSEXP, SEXP num_chainSEXP, SEXP total_iterSEXP, SEXP num_burn_inSEXP, SEXP num_thinSEXP, SEXP num_post_sampleSEXP, SEXP num_treeSEXP, SEXP step_probSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nuSEXP, SEXP alpha2SEXP, SEXP beta2SEXP, SEXP nu2SEXP, SEXP num_tree_modSEXP, SEXP lambda_expSEXP, SEXP lambda_outSEXP, SEXP lambda_modSEXP, SEXP boot_sizeSEXP, SEXP is_binary_trtSEXP, SEXP parallelSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type var_count(var_countSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type var_prob(var_probSEXP);
    Rcpp::traits::input_parameter< double& >::type sigma2_exp(sigma2_expSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sigma2_out_hist(sigma2_out_histSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type dir_alpha_hist(dir_alpha_histSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type trt(trtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type trt_treated(trt_treatedSEXP);
    Rcpp::traits::input_parameter< const double >::type trt_control(trt_controlSEXP);
    Rcpp::traits::input_parameter< const int >::type chain_idx(chain_idxSEXP);
    Rcpp::traits::input_parameter< const int >::type num_chain(num_chainSEXP);
    Rcpp::traits::input_parameter< const int >::type total_iter(total_iterSEXP);
    Rcpp::traits::input_parameter< const int >::type num_burn_in(num_burn_inSEXP);
    Rcpp::traits::input_parameter< const int >::type num_thin(num_thinSEXP);
    Rcpp::traits::input_parameter< const int >::type num_post_sample(num_post_sampleSEXP);
    Rcpp::traits::input_parameter< const int >::type num_tree(num_treeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type step_prob(step_probSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha2(alpha2SEXP);
    Rcpp::traits::input_parameter< const double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< const double >::type nu2(nu2SEXP);
    Rcpp::traits::input_parameter< const int >::type num_tree_mod(num_tree_modSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_exp(lambda_expSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_out(lambda_outSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_mod(lambda_modSEXP);
    Rcpp::traits::input_parameter< const int >::type boot_size(boot_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_binary_trt(is_binary_trtSEXP);
    Rcpp::traits::input_parameter< const bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    fit_bcf_hetero(Y1, Y0, var_count, var_prob, sigma2_exp, sigma2_out_hist, dir_alpha_hist, Y, trt, X, trt_treated, trt_control, chain_idx, num_chain, total_iter, num_burn_in, num_thin, num_post_sample, num_tree, step_prob, alpha, beta, nu, alpha2, beta2, nu2, num_tree_mod, lambda_exp, lambda_out, lambda_mod, boot_size, is_binary_trt, parallel, verbose);
    return R_NilValue;
END_RCPP
}
// fit_mbart
void fit_mbart(NumericVector& Y1, NumericVector& Y0, NumericMatrix& var_count, NumericVector& var_prob, double& sigma2_exp, NumericVector& sigma2_out_hist, NumericVector& dir_alpha_hist, const NumericVector& Y, const NumericVector& trt, const NumericMatrix& X, const double trt_treated, const double trt_control, const int chain_idx, const int num_chain, const int total_iter, const int num_burn_in, const int num_thin, const int num_post_sample, const int num_tree, const NumericVector& step_prob, const double alpha, const double beta, const double nu, const double lambda_exp, const double lambda_out, const int boot_size, const bool is_binary_trt, const bool parallel, const bool verbose);
RcppExport SEXP _bcfac_fit_mbart(SEXP Y1SEXP, SEXP Y0SEXP, SEXP var_countSEXP, SEXP var_probSEXP, SEXP sigma2_expSEXP, SEXP sigma2_out_histSEXP, SEXP dir_alpha_histSEXP, SEXP YSEXP, SEXP trtSEXP, SEXP XSEXP, SEXP trt_treatedSEXP, SEXP trt_controlSEXP, SEXP chain_idxSEXP, SEXP num_chainSEXP, SEXP total_iterSEXP, SEXP num_burn_inSEXP, SEXP num_thinSEXP, SEXP num_post_sampleSEXP, SEXP num_treeSEXP, SEXP step_probSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nuSEXP, SEXP lambda_expSEXP, SEXP lambda_outSEXP, SEXP boot_sizeSEXP, SEXP is_binary_trtSEXP, SEXP parallelSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type var_count(var_countSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type var_prob(var_probSEXP);
    Rcpp::traits::input_parameter< double& >::type sigma2_exp(sigma2_expSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sigma2_out_hist(sigma2_out_histSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type dir_alpha_hist(dir_alpha_histSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type trt(trtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double >::type trt_treated(trt_treatedSEXP);
    Rcpp::traits::input_parameter< const double >::type trt_control(trt_controlSEXP);
    Rcpp::traits::input_parameter< const int >::type chain_idx(chain_idxSEXP);
    Rcpp::traits::input_parameter< const int >::type num_chain(num_chainSEXP);
    Rcpp::traits::input_parameter< const int >::type total_iter(total_iterSEXP);
    Rcpp::traits::input_parameter< const int >::type num_burn_in(num_burn_inSEXP);
    Rcpp::traits::input_parameter< const int >::type num_thin(num_thinSEXP);
    Rcpp::traits::input_parameter< const int >::type num_post_sample(num_post_sampleSEXP);
    Rcpp::traits::input_parameter< const int >::type num_tree(num_treeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type step_prob(step_probSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_exp(lambda_expSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_out(lambda_outSEXP);
    Rcpp::traits::input_parameter< const int >::type boot_size(boot_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_binary_trt(is_binary_trtSEXP);
    Rcpp::traits::input_parameter< const bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    fit_mbart(Y1, Y0, var_count, var_prob, sigma2_exp, sigma2_out_hist, dir_alpha_hist, Y, trt, X, trt_treated, trt_control, chain_idx, num_chain, total_iter, num_burn_in, num_thin, num_post_sample, num_tree, step_prob, alpha, beta, nu, lambda_exp, lambda_out, boot_size, is_binary_trt, parallel, verbose);
    return R_NilValue;
END_RCPP
}
// fit_sbart
void fit_sbart(NumericVector& Y1, NumericVector& Y0, NumericMatrix& var_count, NumericVector& var_prob, NumericVector& sigma2_out1_hist, NumericVector& sigma2_out0_hist, NumericVector& dir_alpha_hist, const NumericVector& Y_treated, const NumericVector& Y_control, const NumericVector& trt, const NumericMatrix& X, const NumericMatrix& X_treated, const NumericMatrix& X_control, const int chain_idx, const int num_chain, const int total_iter, const int num_burn_in, const int num_thin, const int num_post_sample, const int num_tree, const NumericVector step_prob, const double alpha, const double beta, const double nu, const double lambda_out1, const double lambda_out0, const int boot_size, const bool is_binary_trt, const bool parallel, const bool verbose);
RcppExport SEXP _bcfac_fit_sbart(SEXP Y1SEXP, SEXP Y0SEXP, SEXP var_countSEXP, SEXP var_probSEXP, SEXP sigma2_out1_histSEXP, SEXP sigma2_out0_histSEXP, SEXP dir_alpha_histSEXP, SEXP Y_treatedSEXP, SEXP Y_controlSEXP, SEXP trtSEXP, SEXP XSEXP, SEXP X_treatedSEXP, SEXP X_controlSEXP, SEXP chain_idxSEXP, SEXP num_chainSEXP, SEXP total_iterSEXP, SEXP num_burn_inSEXP, SEXP num_thinSEXP, SEXP num_post_sampleSEXP, SEXP num_treeSEXP, SEXP step_probSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nuSEXP, SEXP lambda_out1SEXP, SEXP lambda_out0SEXP, SEXP boot_sizeSEXP, SEXP is_binary_trtSEXP, SEXP parallelSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y0(Y0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type var_count(var_countSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type var_prob(var_probSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sigma2_out1_hist(sigma2_out1_histSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type sigma2_out0_hist(sigma2_out0_histSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type dir_alpha_hist(dir_alpha_histSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Y_treated(Y_treatedSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Y_control(Y_controlSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type trt(trtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_treated(X_treatedSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X_control(X_controlSEXP);
    Rcpp::traits::input_parameter< const int >::type chain_idx(chain_idxSEXP);
    Rcpp::traits::input_parameter< const int >::type num_chain(num_chainSEXP);
    Rcpp::traits::input_parameter< const int >::type total_iter(total_iterSEXP);
    Rcpp::traits::input_parameter< const int >::type num_burn_in(num_burn_inSEXP);
    Rcpp::traits::input_parameter< const int >::type num_thin(num_thinSEXP);
    Rcpp::traits::input_parameter< const int >::type num_post_sample(num_post_sampleSEXP);
    Rcpp::traits::input_parameter< const int >::type num_tree(num_treeSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type step_prob(step_probSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_out1(lambda_out1SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda_out0(lambda_out0SEXP);
    Rcpp::traits::input_parameter< const int >::type boot_size(boot_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_binary_trt(is_binary_trtSEXP);
    Rcpp::traits::input_parameter< const bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    fit_sbart(Y1, Y0, var_count, var_prob, sigma2_out1_hist, sigma2_out0_hist, dir_alpha_hist, Y_treated, Y_control, trt, X, X_treated, X_control, chain_idx, num_chain, total_iter, num_burn_in, num_thin, num_post_sample, num_tree, step_prob, alpha, beta, nu, lambda_out1, lambda_out0, boot_size, is_binary_trt, parallel, verbose);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bcfac_count_omp_thread", (DL_FUNC) &_bcfac_count_omp_thread, 0},
    {"_bcfac_fit_bcf_hetero", (DL_FUNC) &_bcfac_fit_bcf_hetero, 34},
    {"_bcfac_fit_mbart", (DL_FUNC) &_bcfac_fit_mbart, 29},
    {"_bcfac_fit_sbart", (DL_FUNC) &_bcfac_fit_sbart, 30},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcfac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
