#include <RcppArmadillo.h>

using namespace Rcpp;

//compareKgroups

List compareKgroups(int ngroup, LogicalMatrix y, LogicalMatrix missing, int demleader, int repleader, IntegerVector group, int nsamples, int burn, int thin, int printevery, int samplezeta, int startconfig, double varalpprop);

RcppExport SEXP _compareKgroups(SEXP ngroupSEXP, SEXP ySEXP, SEXP missingSEXP, SEXP demleaderSEXP, SEXP repleaderSEXP, SEXP groupSEXP, SEXP nsamplesSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP printeverySEXP, SEXP samplezetaSEXP, SEXP startconfigSEXP, SEXP varalppropSEXP)
{
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope Rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ngroup(ngroupSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type missing(missingSEXP);
    Rcpp::traits::input_parameter< int >::type demleader(demleaderSEXP);
    Rcpp::traits::input_parameter< int >::type repleader(repleaderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type nsamples(nsamplesSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< int >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< int >::type samplezeta(samplezetaSEXP);
    Rcpp::traits::input_parameter< int >::type startconfig(startconfigSEXP);
    Rcpp::traits::input_parameter< double >::type varalpprop(varalppropSEXP);
    rcpp_result_gen = Rcpp::wrap(compareKgroups(ngroup, y, missing, demleader, repleader, group, nsamples, burn, thin, printevery, samplezeta, startconfig, varalpprop));
    return rcpp_result_gen;
    END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_compareKgroups", (DL_FUNC) _compareKgroups, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_compareK(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
