#include <RcppArmadillo.h>
#include "samplers.hpp"
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
//////////////////////////////////////////////////////////
List compareKgroups(int ngroup, LogicalMatrix y, LogicalMatrix missing, int demleader, int repleader, IntegerVector group, int nsamples, int burn, int thin, int printevery, int samplezeta, int startconfig, double varalpprop)
{
    int II = y.nrow();
    int JJ = y.ncol();
    int KK = ngroup;
    Rcout << "I: " << II << " J: " << JJ << "\n" << "Iterations: " << nsamples << " Burn: " << burn << " Thin: " << thin << "\n" << "Print every: " << printevery << " \nDemLeader: " << demleader << " RepLeader: " << repleader << "\n" << std::endl;
    
    /* Initialization */
    NumericMatrix z(II, JJ);
    NumericVector mu(JJ);
    NumericVector alpha(JJ);
    NumericMatrix beta(II, KK);
    IntegerMatrix omega(II, ngroup);
    for(int i=0; i<II; i++){
        for(int k=0; k<KK; k++){
            beta(i, k)=0;
            if(startconfig==0){
                omega(i, k)=0;
            }else{
                omega(i, k)=k;
            }
        }
        for(int j=0; j<JJ; j++){
            z(i, j)=0;
        }
    }
    for(int j=0; j<JJ; j++){
        mu[j]=0;
        alpha[j]=0;
    }
    double alp     = 0.1;
    double pialpha = 0.9;
    double eta     = 0.0;
    double kappa2  = 1.0;
    double w2      = 1.0;
    double rho     = 0.0;
    double tau2    = 1.0;
    
    /* Creating variables to store samples obtained after burn */
    arma::cube betaOut = arma::zeros<arma::cube>(II, KK, nsamples-burn);
    arma::cube omegaOut = arma::zeros<arma::cube>(II, KK, nsamples-burn);
    NumericMatrix alphaOut(nsamples-burn, JJ);
    NumericMatrix muOut(nsamples-burn, JJ);
    NumericVector phiOut(nsamples-burn);
    NumericVector pialphaOut(nsamples-burn);

    /*mcmc starts here*/
    for(int m = 0; m<nsamples; m++){
        for(int ss = 0; ss<thin; ss++){
    
            R_CheckUserInterrupt();
            
            sample_alpha(II,JJ,KK,z,mu,alpha,beta,group,w2,pialpha);
            sample_pialpha(II, alpha, &pialpha);
            sample_mu(II,JJ,KK,z,mu,alpha,beta,group,eta,kappa2);
            sample_beta(II,JJ,KK,beta,alpha,mu,z,tau2,rho,group,omega);
            normalize2(II,JJ,KK,group,beta,alpha,mu,demleader,repleader);
            sample_z(II,JJ,KK,y,z,mu,alpha,beta,group);
            if(samplezeta==1){
                sample_omega(II,JJ,omega,alp,z,mu,alpha,beta,group,rho,tau2,KK);
                sample_alp(II, KK, omega, &alp, varalpprop);
            }
            sample_eta(JJ, &eta, kappa2, mu);
            sample_kappa2(JJ,&kappa2,eta,mu);
            sample_w2(JJ,&w2,alpha);
            sample_rho(II,KK,&rho,tau2,beta,omega);
            sample_tau2(II,KK,&tau2,beta,rho,omega);
            sample_pi(JJ,&pialpha,alpha);
            sample_missings(II,JJ,KK,y,missing,mu,alpha,beta,group);
            
            R_CheckUserInterrupt();
        }

        if(m % printevery == 0){
            Rcout << "Iterations completed: " << m; // << "\n";
            Rcout << "        phi (alp) = " << alp << "\n";
            /* printf("Beta0_98: %f,Beta1_98: %f, Alpha98: %f, Mu98: %f, z1010: %lf \n", beta[98][0],beta[98][1],alpha[98],mu[98], z[10][10]);
             if(KK==3) printf("Beta2_98:  %f\n",beta[98][2]);
             */
            /*Rcout << "beta0_0:  " << beta(0, 0) << " beta1_0:  " << beta(0, 1) << " beta2_0:  " << beta(0, 2) << "\n beta0_1:  " << beta(1, 0) << " beta1_1:  " << beta(1, 1) << " beta2_1: " << beta(1, 2)  << "\n ";*/
            /*  printf("Z0631: %f, Z0632: %f, Z10: %f, Z11: %f, Z12: %f \n", z[0][631],z[0][632],z[1][0],z[1][1],z[1][2]);
             printf("Eta: %f, Kappa2: %f, W2: %f, \n Rho: %f, Tau2: %f, PiAlpha: %f Alp: %f \n", eta, kappa2, w2, rho, tau2,pialpha, alp);
             */
        }

        if(m>=burn){
            for(int i=0; i<II; i++){
                for(int k=0; k<KK; k++){
                    betaOut(i, k, m-burn) = beta(i, k);
                    omegaOut(i, k, m-burn) = omega(i, k);
                    R_CheckUserInterrupt();
                }
            }
            
            alphaOut(m-burn,_) = alpha;
            muOut(m-burn,_) = mu;
            phiOut(m-burn) = alp;
            pialphaOut(m-burn) = pialpha;

            R_CheckUserInterrupt();
        }
    }
    
    return List::create(Named("beta") = betaOut, Named("alpha") = alphaOut, Named("mu") = muOut, Named("omega") = omegaOut, Named("phi") = phiOut, Named("pialpha") = pialphaOut);
    
}




