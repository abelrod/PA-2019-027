#ifndef KGROUPSAMPLERS_HPP
#define KGROUPSAMPLERS_HPP

void normalize1(int II,int JJ, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, int repleader);

void normalize2(int II, int JJ,int KK, Rcpp::IntegerVector &gamma, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, Rcpp::NumericVector &mu, int demleader, int repleader);

void remove_and_relabel(int entry, int KK, Rcpp::IntegerVector &vector);

int maximum_vector(int KK, Rcpp::IntegerVector &vector);

void renormalize_logprobs(int LL, Rcpp::NumericVector &q);

int sample_list(int LL, Rcpp::NumericVector &q);

double logprior(int LL,int KK, Rcpp::IntegerVector &w, double alp);

double loglik(int i,int II,int JJ,int KK,Rcpp::IntegerVector &gamma, Rcpp::IntegerVector &w, Rcpp::NumericVector &z, Rcpp::NumericVector &mu ,Rcpp::NumericVector &alpha,double rho,double tau2);

void sample_omega(int II, int JJ, Rcpp::IntegerMatrix &omega, double alp, Rcpp::NumericVector &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double rho,double tau2, int KK);

double logalp_dist(double sum, int KK, int II, double alp);

void sample_alp(int II, int KK, Rcpp::IntegerMatrix &omega, double *alp, double varalpprop);

void sample_beta(int II, int JJ,int KK, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, Rcpp::NumericVector &mu, Rcpp::NumericMatrix &z, double tau2, double rho, Rcpp::IntegerVector &gamma, Rcpp::IntegerMatrix &omega);

void sample_pialpha(int II, Rcpp::NumericVector alpha, double *pialpha);

void sample_alpha(int II,int JJ, int KK, Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double w2,double pialpha);

void sample_mu(int II,int JJ, int KK,Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double eta,double kappa2);

double sample_cond_norm (double eta, double sigma);

void sample_z(int II,int JJ,int KK, Rcpp::LogicalMatrix &y, Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma);

void sample_pi(int JJ, double *p, Rcpp::NumericVector &parm);

void sample_eta(int JJ, double *eta, double kappa2, Rcpp::NumericVector &mu);


void sample_kappa2(int JJ,double *kappa2,double eta, Rcpp::NumericVector &mu);

void sample_w2(int JJ, double *w2, Rcpp::NumericVector &alpha);

void sample_rho(int II,int KK, double *rho, double tau2, Rcpp::NumericMatrix &beta, Rcpp::IntegerMatrix &omega);

void sample_tau2(int II,int KK, double *tau2, Rcpp::NumericMatrix &beta, double rho, Rcpp::IntegerMatrix &omega);

void sample_missings(int II,int JJ,int KK, Rcpp::LogicalMatrix &y, Rcpp::LogicalMatrix &missing, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma);

#endif

