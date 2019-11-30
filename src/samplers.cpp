#define pi 3.1415926

#include <Rmath.h>
#include <Rcpp.h>
#include <RcppTN.h>

using namespace Rcpp;

int max(int a , int b)
{
    if(a > b) return a;
    else return b;
}

double max(double a , double b)
{
    if(a > b) return a;
    else return b;
}

//normalize

void normalize1(int II,int JJ, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, int repleader){
    int i, j;
    if(beta(repleader, 0)<0){
        for(i=0; i<II; i++){
            beta(i, 0) = -beta(i, 0);
            beta(i, 1) = -beta(i, 1);
        }
        for(j = 0; j < JJ; j++){
            alpha[j] = -alpha[j];
        }
    }
}

void normalize2(int II, int JJ,int KK, Rcpp::IntegerVector &gamma, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, Rcpp::NumericVector &mu, int demleader, int repleader){
    double aaa, bbb, ccc, ddd;
    int i, j, k;
    aaa = beta(demleader-1, 0);  //C starts counting at zero!!
    bbb = beta(repleader-1, 0);
    ccc = (bbb+aaa)/2.0;
    ddd = (bbb-aaa)/2.0;
    //printf("aaa:%f, bbb:%f, ccc:%f, ddd:%f \n",aaa,bbb,ccc,ddd);
    for(j=0; j<JJ; j++){
        mu[j] = mu[j] + alpha[j]*ccc;
        alpha[j] = alpha[j] * ddd;
    }
    for(i=0; i<II; i++){
        for(k=0; k<KK; k++){
            beta(i, k) = (beta(i, k) - ccc)/ddd;
        }
    }
}



///////////////////////////////////////////////////////////
//samplers





void remove_and_relabel(int entry, int KK, Rcpp::IntegerVector &vector){
    int k;
    int nentry = 0;
    int centry = vector[entry];
    
    //printf("In remove and relabel \n");
    //display_vector_int(KK,vector);
    
    for(k=0; k<KK && nentry<2; k++){
        if(vector[k]==centry){
            nentry++;
        }
    }
    //printf("KK = %d   entry = %d    centry = %d   nentry = %d", KK, entry, centry, nentry);
    vector[entry] = -1;
    if(nentry==1){
        for(k=0; k<KK; k++){
            if(vector[k]>centry){
                vector[k]--;
            }
        }
    }
}



int maximum_vector(int KK, Rcpp::IntegerVector &vector){
    int k;
    int maxval = vector[0];
    for(k=1; k<KK; k++){
        maxval = max(maxval, vector[k]);
    }
    return maxval;
}



void renormalize_logprobs(int LL, Rcpp::NumericVector &q){
    int l;
    double maxq, sumq;
    maxq = q[0];
    for(l=1; l<LL; l++){
        maxq = max(maxq, q[l]);
    }
    sumq = 0;
    for(l=0; l<LL; l++){
        q[l] = exp(q[l] - maxq);
        sumq += q[l];
    }
    for(l=0; l<LL; l++){
        q[l] = q[l]/sumq;
    }
}


int sample_list(int LL, Rcpp::NumericVector &q){
    int l;
    double u = R::runif(0.0, 1.0);
    
    double cumsum;
    l = 0;
    cumsum = q[l];
    while(u>cumsum && l<LL){
        l++;
        cumsum += q[l];
    }
    return l;
}


double logprior(int LL,int KK, Rcpp::IntegerVector &w, double alp){
    double lp;
    int nn[LL];
    int l,k;
    
    lp = lgamma(alp) - lgamma(alp+KK) + LL*log(alp);  //LL = number of groups, KK = number items in w
    for(l=0; l<LL; l++){
        nn[l] = 0;
        for(k=0; k<KK; k++){
            if(w[k] == l) nn[l] += 1;  //number in each group
        }
        lp += lgamma(nn[l]);
    }
    return lp;
}

/* double loglik(int i,int II,int JJ,int KK,int gamma[JJ],int w[KK],double z[II][JJ],double mu[JJ],double alpha[JJ],double rho,double tau2){
 int j,k,l;
 double lp = 0;
 int LL;
 LL = maximum_vector(KK,w)+1;
 double nn[LL], aa[LL], za[LL], zz[LL];
 double zstar;
 
 // for(l=0; l<LL; l++){
 // nn[l] = 0.0;
 // aa[l] = 0.0;
 // za[l] = 0.0;
 // zz[l] = 0.0;
 // }
 
 for(l=0; l<LL; l++){
 nn[l] = 0.0;        //sum is over all j such that.....
 aa[l] = 0.0;
 zz[l] = 0.0;
 za[l] = 0.0;
 //printf("l = %i  aa = %f  za = %f  zz = %f  \n",l,aa[l],za[l],zz[l]);
 
 for(k=0; k<KK; k++){  //each item of w
 //printf("in loglik: w: %i, k: %i, l: %i \n",w[k],k,l);
 if(w[k] == l){     //each item that takes on that value
 
 for(j=0; j<JJ; j++){
 if(gamma[j] == k){   //gamma from that item - j s.t. gamma_j = k where k s.t. w_k = l
 nn[l] += 1.0;
 aa[l] += alpha[j]*alpha[j];
 za[l] += (z[i][j] - mu[j])*(alpha[j]);            //need the parentheses around alpha. WHY?!?!
 zz[l] += (z[i][j] - mu[j])*(z[i][j] - mu[j]);
 //zstar = z[i][j] - mu[j];
 //zz[l] += zstar*zstar;
 //za[l] += zstar*alpha[j];
 }
 }
 //printf("l = %i  aa = %f  za = %f  zz = %f  \n",l,aa[l],za[l],zz[l]);
 //printf("aa0 %f, za0 %f, zz0 %f, aa1 %f, zz1 %f, za1 %f \n",aa[0], za[0],zz[0], aa[1], zz[1], za[1]);
 }
 }
 //-nn[l]*log(2.0*pi)/2.0 at the beginning?
 lp +=  - log(sqrt(tau2))/2.0 - log(aa[l] + 1.0/tau2)/2.0 - (zz[l] + rho*rho/tau2 - pow(za[l]+rho/tau2,2)/(aa[l]+1.0/tau2))/2.0;  //loglikelihood
 //printf("l: %i, lp: %f \n",l,lp);
 }
 //printf("aa0 %f, za0 %f, zz0 %f \n", aa[0], za[0], zz[0]);
 //printf("aa1 %f, za1 %f, zz1 %f \n", aa[1], za[1], zz[1]);
 //printf("end lp: %f\n",lp);
 return lp;
 
 } */


double loglik(int i, int II, int JJ, int KK, Rcpp::IntegerVector &gamma, Rcpp::IntegerVector &w, Rcpp::NumericVector &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha,double rho,double tau2){
    
    int j,l;
    double lp = 0;
    int LL;
    LL = maximum_vector(KK,w)+1;
    double nn[LL], aa[LL], za[LL], zz[LL];
    
    // for(l=0; l<LL; l++){
    // nn[l] = 0.0;
    // aa[l] = 0.0;
    // za[l] = 0.0;
    // zz[l] = 0.0;
    // }
    
    for(l=0; l<LL; l++){
        nn[l] = 0.0;        //sum is over all j such that.....
        aa[l] = 0.0;
        zz[l] = 0.0;
        za[l] = 0.0;
        //printf("l = %i  aa = %f  za = %f  zz = %f  \n",l,aa[l],za[l],zz[l]);
    }
    
    for(j=0; j<JJ; j++){
        nn[w[gamma[j]]] += 1.0;
        aa[w[gamma[j]]] += alpha[j]*alpha[j];
        za[w[gamma[j]]] += (z(i, j) - mu[j])*(alpha[j]);
        zz[w[gamma[j]]] += (z(i, j) - mu[j])*(z(i, j) - mu[j]);
        //zstar = z[i][j] - mu[j];
        //zz[l] += zstar*zstar;
        //za[l] += zstar*alpha[j];
    }
    
    for(l=0; l<LL; l++){
        lp +=  - log(sqrt(tau2))/2.0 - log(aa[l] + 1.0/tau2)/2.0 - (zz[l] + rho*rho/tau2 - pow(za[l]+rho/tau2,2)/(aa[l]+1.0/tau2))/2.0;   //The term log(sqrt(tau2)) is missing a 1/2
        //printf("l = %i  aa = %f  za = %f  zz = %f  \n",l,aa[l],za[l],zz[l]);
        //printf("aa0 %f, za0 %f, zz0 %f, aa1 %f, zz1 %f, za1 %f \n",aa[0], za[0],zz[0], aa[1], zz[1], za[1]);
        //-nn[l]*log(2.0*pi)/2.0 at the beginning?
        //printf("l: %i, lp: %f \n",l,lp);
    }
    
    //printf("aa0 %f, za0 %f, zz0 %f \n", aa[0], za[0], zz[0]);
    //printf("aa1 %f, za1 %f, zz1 %f \n", aa[1], za[1], zz[1]);
    //printf("end lp: %f\n",lp);
    return lp;
    
}

void sample_omega(int II, int JJ, Rcpp::IntegerMatrix &omega, double alp, Rcpp::NumericVector &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double rho,double tau2, int KK){
    int i,k,l;
    int LL,newLL;
    IntegerVector w(KK);
    int nbridge=0;
    int sum;
    
    for(i=0; i<II; i++){
        sum = 0;
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
            sum = sum + w[k];
        }
        if(sum == 0 ) nbridge += 1;
        
        //omega[i][0] = 0;    //fix the first item - then all remaining are either in this group or a new one.
        
        //printf("original w: ");
        //display_vector_int(KK,w);
        
        for(k=0; k<KK; k++){   //each item of w
            remove_and_relabel(k,KK,w);
            //printf("k = %i, after remove and relabel w: ",k);
            //display_vector_int(KK,w);
            LL = maximum_vector(KK,w)+1;    //number of groups, c starts at 0 so max of 0 means 1 group, max of 1 means 2 groups, etc.
            //printf("LL: %i\n logpr\n",LL);
            NumericVector logpr(LL+1);  // +1 for creating new group
            //display_vector_double(LL+1,logpr);
            
            for(l=0; l<LL+1; l++){    //all possible values - each group and a new one
                w[k] = l;
                //printf("\nw\n");
                //display_vector_int(KK,w);
                newLL = maximum_vector(KK,w)+1;    //how many groups are there now
                logpr[l]  = logprior(newLL,KK,w,alp);
                //printf("l: %i, newLL: %i, logprior: %f\n",l,newLL,logpr[l]);
                logpr[l] += loglik(i,II,JJ,KK,gamma,w,z,mu,alpha,rho,tau2);
                //printf("l: %i, +loglik: %f\n",l,logpr[l]);
            }
            //printf("logpr unnormalized: \n");
            //display_vector_double(LL+1,logpr);
            
            renormalize_logprobs(LL+1,logpr);
            w[k] = sample_list(LL+1,logpr);
            //omega[i][k] = w[k];
            
            /*            printf("logpr normalized: \n");
             display_vector_double(LL+1,logpr);
             printf("w: %i \n",w[k]);  */
            
            /*         if(i<5) {
             printf("LL: %i\n",LL);
             display_vector_double(LL+1,logpr);
             } */
        }
        
        for(k=0; k<KK; k++){
            //if(nbridge < KK-1) w[k] = 0.0;  //if not enough bridges, make entire vector 0,0,0 -> bridge
            if(nbridge < 2) w[k] = 0.0;  //if not enough bridges (D+1) D=dim of policy space, make entire vector 0,0,0 -> bridge
            omega(i, k) = w[k];
        }
        
        //display_vector_int(KK,w);
        //printf("\n");
    }
}

double logalp_dist(double sum, int KK, int II, double alp){
    int k;
    double theta, dist, temp;

    theta = exp(lgamma(KK) + lgamma(alp+1.0) - lgamma(KK+alp));
    dist =  II*lgamma(alp) - II*lgamma(alp+KK) + sum*log(alp) - R::pbinom(1,II,theta,0,1);
    dist += lgamma(alp+1.0) - lgamma(alp+KK);
    temp = 0.0;
    for(k=0; k<KK-1; k++){
        temp += 1/(alp+k+1);
    }
    dist += log(temp);
    return dist;
}

void sample_alp(int II, int KK, Rcpp::IntegerMatrix &omega, double *alp, double varalpprop){
    double suma = 0.0;
    double tempalp,star,praccept,uni;
    
    IntegerVector w(KK);
    int i, k;
    int LL;
    
    /*for(i=0; i<II; i++){
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
        }
        LL = maximum_vector(KK,w)+1;  //number of groups - C starts at 0 so add 1
        
        suma += LL;
        sumb += log(latent[i]);
        
        //printf("LL %i suma %f sumb %f latenti %f log %f \n",LL,suma,sumb, 1.0-latent[i], log(1.0-latent[i]));
    }
    
    newa = aalp + suma + 1.0;
    newb = balp - sumb;
    //printf("a %f suma %f newa %f, b%f sumb %f newb %f \n",a,suma,newa,b,sumb,newb);
    // *alp = rgamma(newa,1.0/newb);
     */
    //logalp = log(*alp);
    for(i=0; i<II; i++){
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
        }
        LL = maximum_vector(KK,w)+1;  //number of groups - C starts at 0 so add 1
        suma += LL;
    }
    
    tempalp = *alp;
    star = R::rlnorm(log(tempalp), varalpprop);
    praccept = logalp_dist(suma,KK,II,star) - logalp_dist(suma,KK,II,tempalp) + log(star) - log(tempalp);

    //Rcout << "alp" << tempalp << "       alpprop: " << star << "       log aceptance prob" << praccept << "\n";

    uni = R::runif(0,1);
    if(praccept >= log(uni)) *alp = star;
    
    
    
    //printf("alp %f \n",*alp);
}



void sample_beta(int II, int JJ,int KK, Rcpp::NumericMatrix &beta, Rcpp::NumericVector &alpha, Rcpp::NumericVector &mu, Rcpp::NumericMatrix &z, double tau2, double rho, Rcpp::IntegerVector &gamma, Rcpp::IntegerMatrix &omega){
    int i,j,k,l;
    double mean,var,sd;
    IntegerVector w(KK);
    double zstar;
    int LL;
    for(i=0; i<II; i++){
        
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
        }
        LL = maximum_vector(KK,w)+1;
        double bet[LL];
        double aa[LL];
        double za[LL];
        //display_vector_int(KK,w);
        for(l=0; l<LL; l++){    //each possible value
            aa[l] = 0;
            za[l] = 0;
            for(k=0; k<KK; k++){    //each item
                if(omega(i, k)==l){    //does that item equal the possible value
                    for(j=0; j<JJ; j++){
                        if(gamma[j]==k){     // {j : gamma_j = {k : w_k = l} }  j when gamma_j is the position of item that equals that possible value
                            aa[l] += alpha[j]*alpha[j];
                            //za[l] += (z[i][j]-mu[j])*(alpha[j]);
                            zstar = z(i, j) - mu[j];
                            za[l] += zstar * alpha[j];
                        }
                    }
                }
            }
            var = 1.0/(aa[l] + 1.0/tau2);
            sd = sqrt(var);
            mean = var * (za[l] + rho/tau2);
            bet[l] = R::rnorm(mean,sd);     //one beta for each possible value - repeats of that value will have the same beta
            //if(i<11) printf("l: %i, bet: %f \n",l,bet[l]);
            //for(k=0; k<KK; k++){
            //if(w[k] == l) beta[i][k] = bet[w[k]];
            //beta[i][k] = bet[w[k]];
            //printf("beta: %f  ",beta[i][k]);
            //}
            //printf("\n");
            //if(i<11) printf("beta %i mean %f sd %f \n",l,mean,sd);
            /*             for(k=0; k<KK; k++){
             if(w[k]==l) beta[i][k] = bet[k];
             if(i>94) printf("k %i wk %i l %i beta %i: %f  \n",k,w[k],l,k,beta[i][k]);
             }  */
        }
        //if(i<11) display_vector_double(LL,bet);
        //printf("\n");
        for(l=0; l<LL; l++){
            for(k=0; k<KK; k++){
                if(w[k]==l) beta(i, k) = bet[l];
                //if(i<11) printf("k %i wk %i l %i beta %i: %f  \n",k,w[k],l,k,beta[i][k]);
            } }
    }
    
}


void sample_pialpha(int II, Rcpp::NumericVector alpha, double *pialpha){
    int i;
    int numzeros = 0;
    for(i=0; i<II; i++){
        if(alpha[i]==0) numzeros++;
    }
    *pialpha = R::rbeta(1.0 + II - numzeros, 1.0 + numzeros);
}


void sample_alpha(int II,int JJ, int KK, Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double w2, double pialpha){
    int i,j;
    double var, sd, mean;
    double logscore, pr0;
    double zstar,v;
    double u;
    
    for(j=0; j<JJ; j++){
        alpha[j]=0;
        zstar=0;
        v=0;
        for(i=0; i<II; i++){
            zstar += (z(i, j)-mu[j]) * beta(i, gamma[j]);
            v += beta(i, gamma[j]) * beta(i, gamma[j]);
        }
        var = 1.0/(v+1.0/w2);
        sd = sqrt(var);
        mean = var*zstar;
        logscore = log(pialpha) - log(1.0-pialpha) - 0.5*log(w2) + 0.5*log(var) + 0.5*mean*mean/var;
        pr0 = 1.0/(1.0+exp(logscore));   //prob that alpha = 0
        u = R::runif(0.0,1.0);
        if(u > pr0) alpha[j] = R::rnorm(mean,sd);
    }
    
}

void sample_mu(int II,int JJ, int KK, Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma,double eta,double kappa2){
    int i,j;
    double var, sd, mean;
    double zstar;
    
    for(j=0; j<JJ; j++){
        zstar=0;
        for(i=0; i<II; i++){
            zstar += z(i, j) - alpha[j]*beta(i, gamma[j]);
            //if(zstar !=0) printf("zstar: %f  , z: %f , alphabeta: %f \n",zstar,z[i][j], alpha[j]*beta[i][gamma[j]]);
        }
        var = 1.0/(II+1.0/kappa2);
        sd = sqrt(var);
        mean = var * (zstar + eta/kappa2);
        //mu[j] = gsl_ran_gaussian(r,sd) + mean;
        mu[j] = R::rnorm(mean,sd);
        //mew = mu[j];
    }
    
    /*   //printf("mean %f, sd %f \n",mean, sd);
     printf("Mu gsl: %f \n", mew);
     //
     mew  = rnorm(mean,sd);
     //
     printf("Mu R: %f \n", mew);
     printf("%%%%%%%%%%%% \n"); */
    //display_matrix_double(II,JJ,z);
}


double sample_cond_norm (double eta, double sigma)        //AR version
{
    double retvalue = R::rnorm(0.0, sigma);
    while (retvalue < eta){
        retvalue = R::rnorm(0.0, sigma);
    }
    return retvalue;
}


void sample_z(int II,int JJ,int KK, Rcpp::LogicalMatrix &y, Rcpp::NumericMatrix &z, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma){
    int i, j;
    double mean;
    for(i=0; i<II; i++){
        for(j=0; j<JJ; j++){
            mean = mu[j] + alpha[j]*beta(i, gamma[j]);
            /*         printf("mu: %f  alpha: %f  beta: %f  gamma: %i\n",mu[j],alpha[j],beta[i][gamma[j]],gamma[j]);
             printf("mean: %f\n",mean); */
            if(y(i, j)==0){
                /* if(mean >= 0.0){
                 z[i][j] = mean - gsl_ran_gaussian_tail(r, mean, 1.0);
                 }else{
                 z[i][j] = mean - sample_cond_norm(r, mean, 1.0);
                 } */
                z(i, j) = RcppTN::rtn1(mean, 1.0, R_NegInf, 0.0);
            }else{
                /* if(mean < 0.0){
                 z[i][j] = mean + gsl_ran_gaussian_tail(r, -mean, 1.0);
                 }else{
                 z[i][j] = mean + sample_cond_norm(r, -mean, 1.0);
                 } */
                z(i, j) = RcppTN::rtn1(mean, 1.0, 0.0, R_PosInf);
            }
            //printf("GSL:    mean: %f    z: %f     y: %d \n", mean, z[i][j], y[i][j]);
            //printf("RTNORM: mean: %f    z: %f     y: %d \n", mean, zztop, y[i][j]);
        }
        //fprintf(file_ptr,"\r\n");
    }
}

void sample_pi(int JJ, double *p, Rcpp::NumericVector &parm)
{
    int j;
    double n=0;
    double a,b;
    double x,y;
    
    for(j=0; j<JJ; j++){
        if(parm[j]!=0) n = n + 1.0;  //how many are diff from 0
    }
    a = 1.0 + n;
    b = 1.0 + JJ - n;  //how many are zero
    //*p = gsl_ran_beta(r,a,b);
    x = R::rgamma(a,1.0);
    y = R::rgamma(b,1.0);
    *p = x/(x+y);
    
    /*   printf("GSL pialpha: %f \n", *p);
     pialpha = rbeta(a,b);
     printf("R pialpha: %f \n", pialpha); */
}

void sample_eta(int JJ, double *eta, double kappa2, Rcpp::NumericVector &mu){
    int j;
    double var, sd, mean;
    double m=0;
    
    for(j=0; j<JJ; j++){
        m += mu[j];
    }
    var = 1.0/(1.0+JJ/kappa2);
    sd = sqrt(var);
    mean = var * m/kappa2;
    //*eta = gsl_ran_gaussian(r,sd) + mean;
    *eta = R::rnorm(mean,sd);
    
    /*   printf("GSL eta: %f \n", *eta);
     et = rnorm(mean,sd);
     printf("R eta: %f \n", et);*/
    //printf("eta mean %f, sd %f \n",mean, sd);
}

void sample_kappa2(int JJ,double *kappa2,double eta, Rcpp::NumericVector &mu){
    int j;
    double m=0;
    double a,b;
    
    for(j=0; j<JJ; j++){
        m += (mu[j] - eta)*(mu[j] - eta);
    }
    a = 2.0 + JJ/2.0;
    b = 1.0 + m/2.0;
    //*kappa2 = 1.0/gsl_ran_gamma(r,a,1.0/b);
    *kappa2 = 1.0/R::rgamma(a,1.0/b);
    
    /*   printf("GSL k2: %f \n", *kappa2);
     kapp = 1.0/rgamma(a,1.0/b);
     printf("R k2: %f \n", kapp);
     //printf("a=%f, b=%f \n",a,b);
     //printf("eta from kappa2 sampler: %f \n",eta); */
    //printf("kappa2 a=%f, b=%f \n",a,b);
}

void sample_w2(int JJ, double *w2, Rcpp::NumericVector &alpha){
    int j;
    double aa = 0;
    double a, b;
    
    for(j=0; j<JJ; j++){
        aa += alpha[j]*alpha[j];
    }
    a = 2.0 + JJ/2.0;
    b = 1.0 + aa/2.0;
    //*w2 = 1.0/gsl_ran_gamma(r,a,1.0/b);
    *w2 = 1.0/R::rgamma(a,1.0/b);
    
    /*   printf("GSL w2: %f \n", *w2);
     dub = 1.0/rgamma(a,1.0/b);
     printf("R w2: %f \n", dub);
     //printf("a=%f, b=%f \n",a,b); */
}

void sample_rho(int II,int KK, double *rho, double tau2, Rcpp::NumericMatrix &beta, Rcpp::IntegerMatrix &omega){
    int i,k,l;
    double mean, sd, var;
    IntegerVector w(KK);
    int LL=0;
    double bsum=0, omsum=0;
    int first;
    
    
    for(i=0; i<II; i++){
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
        }
        LL = maximum_vector(KK,w)+1;
        double b[LL];
        double om[LL];
        
        for(l=0; l<LL; l++){
            b[l] = 0;
            om[l] = 0;
            first = 1;
            for(k=0; k<KK; k++){
                if(w[k] == l && first==1){    //don't add redundant betas
                    first = 0;
                    om[l] += 1.0;
                    b[l] += beta(i, k);
                    omsum += 1.0;
                    bsum += beta(i, k);
                }
            }
        }
        
        /*     om = om + omega[i];
         if(omega[i]==1) b1 += beta[i][0];
         if(omega[i]==0) b0 += beta[i][0] + beta[i][1]; */
        
    }
    
    //printf("omsum  %f  \n",omsum);
    var = 1.0/(1.0+omsum/tau2);    //om = number of times each type occured, used to be II + #times diff groups for 2
    sd = sqrt(var);
    mean = var * bsum/tau2;
    //*rho = gsl_ran_gaussian(r,sd) + mean;
    *rho = R::rnorm(mean,sd);
    
    //printf("rho mean %f, sd %f \n",mean, sd);
    /*   printf("GSL rho: %f \n", *rho);
     row = rnorm(mean,sd);
     printf("R rho: %f \n", row);
     //printf("mean %f, sd %f \n",mean, sd); */
}

void sample_tau2(int II,int KK, double *tau2, Rcpp::NumericMatrix &beta, double rho, Rcpp::IntegerMatrix &omega){
    int i,k,l;
    double a, b;
    double bpsum=0;
    double omsum=0;
    IntegerVector w(KK);
    int LL;
    int first;
    
    /*     for(i=0; i<II; i++){
     for(k=0; k<KK; k++){
     w[k] = omega[i][k];
     }
     LL = max(LL,maximum_vector(KK,w)+1);
     }
     double bp[LL];
     double om[LL];
     for(l=0; l<LL; l++){
     bp[l] = 0;
     om[l] = 0;
     } */
    
    for(i=0; i<II; i++){
        
        for(k=0; k<KK; k++){
            w[k] = omega(i, k);
        }
        LL = maximum_vector(KK,w)+1;
        double bp[LL];
        double om[LL];
        
        //display_vector_int(KK,w);
        
        for(l=0; l<LL; l++){
            om[l] = 0;
            bp[l] = 0;
            first = 1;
            for(k=0; k<KK; k++){
                if(w[k] == l && first==1){
                    om[l] += 1.0;
                    bp[l] += (beta(i,k)-rho)*(beta(i, k)-rho);
                    first = 0;
                    omsum += 1.0;
                    bpsum += (beta(i, k)-rho)*(beta(i, k)-rho);
                }
            }
            //printf("om  %f  omsum  %f\n",om[l],omsum);
        }
        
        /*     om += omega[i];
         if(omega[i]==0) bp0 += (beta[i][0]-rho)*(beta[i][0]-rho) + (beta[i][1]-rho)*(beta[i][1]-rho);
         if(omega[i]==1) bp1 += (beta[i][0]-rho)*(beta[i][0]-rho); */
        //printf("beta0: %f, rho: %f , bp: %f \n",beta[i][0], rho, bp);
        
    }
    
    
    a = 2.0 + omsum/2.0;
    b = 1.0 + bpsum/2.0;
    //*tau2 = 1.0/gsl_ran_gamma(r,a,1.0/b);
    *tau2 = 1.0/R::rgamma(a,1.0/b);
    
    //printf("tau2 a=%f, b=%f \n",a,b);
    /*   printf("GSL t2: %f \n", *tau2);
     taw = 1.0/rgamma(a,1.0/b);
     printf("R t2: %f \n", taw);
     //printf("bp: %f, zeta2: %f \n",bp,*zeta2);
     //printf("a=%f, b=%f \n",a,b); */
}

void sample_missings(int II,int JJ,int KK, Rcpp::LogicalMatrix &y, Rcpp::LogicalMatrix &missing, Rcpp::NumericVector &mu, Rcpp::NumericVector &alpha, Rcpp::NumericMatrix &beta, Rcpp::IntegerVector &gamma){
    int i, j;
    double x, p;
    
    for(i=0; i<II; i++){
        for(j=0; j<JJ; j++){
            if(missing(i, j) == 1){
                x = mu[j] + alpha[j]*beta(i, gamma[j]);
                //p = gsl_cdf_ugaussian_P(x);
                //y[i][j] = gsl_ran_bernoulli(r,p);
                //printf("p: %f, y: %i \n", p, y[i][j]);
                p = R::pnorm(x,0.0,1.0,TRUE,FALSE);
                y(i, j) = R::rbinom(1,p);
                //printf("p: %f, y: %i \n", p, why);
            }
        }
    }
}

