#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;    

// This is the E step of EM algorithm, formating for psi and p should be done in R
// psi should be a col vector while p should be a matrix with row as site and column as period
// [[Rcpp::export]]
double Estep_cpp( const arma::vec & psi, const arma::mat & p, // current psi and p, parameter
              const arma::vec & psi_t, const arma::mat & p_t,// old psi and p in expectation step
              const arma::mat & Y, // just the data, change all missing to 0 and record missing using missing
              const arma::mat & missing, // a matrix indicate which 0s in Y is actually missing 
              const arma::uvec & non_det){// those without detection, we do not want to calculate this each iteration 
    
    // probability of being occupied given no detection, its P(Z_i=1|Y_ij,\theta_t), only used when Y_ij=0 for all j, i.e. no detection at all
    arma::vec pPlus = exp(sum(log(1-p_t+1e-10),1) + log(psi_t+1e-10)); 
    // the full log likelihood when there exist at least one detection l(theta|Z_i,Y_ij) : 
    arma::vec p_det_Occu = log(psi+1e-10) + // occupancy rate
                            sum(((1-missing) % Y) % log(p+1e-10),1) +  // detection 
                            sum(((1-missing) % (1-Y)) % log(1-p+1e-10),1); // non detections
    arma::vec p_nOccu = log(1-psi+1e-10); // probability of absense

    arma::vec Q = p_det_Occu;
    Q(non_det) = pPlus(non_det) % p_det_Occu(non_det) +  // key part of expectation, this is the P(Z=1|Y_ij=0,\theta_t)*l(\theta|Z_i=1,Y_ij=0) 
                (1-pPlus(non_det)) % p_nOccu(non_det);   // similarly P(Z_i=0|Y_ij=0,\theta_t)*l(\theta|Z_i=0,Y_ij=0) 

    return(accu(Q));
}


// calculate the min-linear logistic psi
// [[Rcpp::export]]
arma::vec getpsi(const List & Desings, // a list of designs for each groups
                const List & betas, // list of regression coefficients 
                const arma::vec & alpha, // the intercept, should be size G
                int n_site, // number of sites, basecally number of rows of Designs[[1]]
                int G){// number of groups

    arma::mat linPred(n_site, G, arma::fill::zeros);// the linear predictors (plus intercept)
    //arma::mat design_temp;
    //arma::vec beta_temp;
    for (int i = 0; i < G; i++)
    {
        arma::mat design_temp = Desings[i];
        arma::vec beta_temp = betas[i];
        linPred.col(i) = design_temp * beta_temp;
    }

    linPred.each_row() += alpha.t(); 

    arma::vec psi = min(linPred,1);// take minimum 
    psi = 1/(1+exp(-psi));
    return(psi);
}



// calculate the min-linear logistic p
// [[Rcpp::export]]
arma::mat getp(const List & Designs_det,// should be a list of list, each element for a period
                const List & betas_det,// just a list, each element is a vector for each group
                const arma::vec & alpha,// similarly intercept
                int n_period,// number of periods 
                int n_site, // number of sites
                int G // number of groups
                ){
    arma::mat pmat(n_site,n_period, arma::fill::zeros);
    List design_temp;

    for (int i = 0; i < n_period; i++)
    {
        design_temp = Designs_det[i];
        pmat.col(i) = getpsi(design_temp, betas_det, alpha, n_site, G);
    }
    
    return(pmat);
}
