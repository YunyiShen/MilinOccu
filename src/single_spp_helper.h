#ifndef SINGLE_SPP_HELPER_H
#define SINGLE_SPP_HELPER_H

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma; 

double logLik_cpp( const arma::vec & psi, const arma::mat & p, // current psi and p, parameter
              const arma::mat & Y, // just the data, change all missing to 0 and record missing using missing
              const arma::mat & missing, // a matrix indicate which 0s in Y is actually missing 
              const arma::uvec & non_det);

arma::vec getpsi(const List & Designs, // a list of designs for each groups
                const List & betas, // list of regression coefficients 
                const arma::vec & alpha, // the intercept, should be size G
                int n_site, // number of sites, basecally number of rows of Designs[[1]]
                int G);

arma::mat getp(const List & Designs_det,// should be a list of list, each element for a period
                const List & betas_det,// just a list, each element is a vector for each group
                const arma::vec & alpha,// similarly intercept
                int n_period,// number of periods 
                int n_site, // number of sites
                int G // number of groups
                );

#endif