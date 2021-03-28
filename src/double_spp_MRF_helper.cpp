#include "single_spp_helper.h"

// we will let user specify coding. 
// This calculate the probability of 4 possible states in occupancy. 
// [[Rcpp::export]]
arma::mat getpsi_MRF(const List & Designs, // a list of designs for each groups
                const List & betas_1, // list of regression coefficients 
                const arma::vec & alpha_1, // the intercept, should be size G
                const List & betas_2, //spp2
                const arma::vec & alpha_2, // alpha for spp2
                const double & gamma, // associations between two spps
                const arma::vec & coding, // +-1 or 01
                int n_site, // number of sites, basecally number of rows of Designs[[1]]
                int G){// number of groups

    arma::vec thr1 = getpsi(Designs, betas_1, alpha_1, n_site, G);
    thr1 = logit<arma::vec>(thr1);
    arma::vec thr2 = getpsi(Designs, betas_2, alpha_2, n_site, G);
    thr2 = logit<arma::vec>(thr2);
    arma::mat psi(n_site, 4, arma::fill::zeros); // probability for 11, 10, 01, 00 
    
    // both presence
    psi.col(0) = exp(thr1 * coding(1) + 
                     thr2 * coding(1) + 
                     gamma * coding(1) * coding(1));
    
    // first presence
    psi.col(1) = exp(thr1 * coding(1) + 
                     thr2 * coding(0) + 
                     gamma * coding(1) * coding(0));
    
    // second presence
    psi.col(2) = exp(thr1 * coding(0) + 
                     thr2 * coding(1) + 
                     gamma * coding(0) * coding(1));
    
    // both absence
    psi.col(3) = exp(thr1 * coding(0) + 
                     thr2 * coding(0) + 
                     gamma * coding(0) * coding(0));

    psi.each_col() /= arma::sum(psi,1); // normalize
    return(psi);
}

// [[Rcpp::export]]
arma::vec getpartition_MRF(const List & Designs, // a list of designs for each groups
                const List & betas_1, // list of regression coefficients 
                const arma::vec & alpha_1, // the intercept, should be size G
                const List & betas_2, //spp2
                const arma::vec & alpha_2, // alpha for spp2
                const double & gamma, // associations between two spps
                const arma::vec & coding, // +-1 or 01
                int n_site, // number of sites, basecally number of rows of Designs[[1]]
                int G){// number of groups

    arma::vec thr1 = getpsi(Designs, betas_1, alpha_1, n_site, G);
    thr1 = logit<arma::vec>(thr1);
    arma::vec thr2 = getpsi(Designs, betas_2, alpha_2, n_site, G);
    thr2 = logit<arma::vec>(thr2);
    arma::mat psi(n_site, 4, arma::fill::zeros); // probability for 11, 10, 01, 00 
    
    // both presence
    psi.col(0) = exp(thr1 * coding(1) + 
                     thr2 * coding(1) + 
                     gamma * coding(1) * coding(1));
    
    // first presence
    psi.col(1) = exp(thr1 * coding(1) + 
                     thr2 * coding(0) + 
                     gamma * coding(1) * coding(0));
    
    // second presence
    psi.col(2) = exp(thr1 * coding(0) + 
                     thr2 * coding(1) + 
                     gamma * coding(0) * coding(1));
    
    // both absence
    psi.col(3) = exp(thr1 * coding(0) + 
                     thr2 * coding(0) + 
                     gamma * coding(0) * coding(0));

    return(arma::sum(psi,1));
}


// calculate detection rate when both species exist. 
// [[Rcpp::export]]
arma::vec getp_11(const List & Designs_det,// should be a list of list, each element for a period
                const arma::mat & Y_1, // detection of spp1
                const arma::mat & Y_2, // detection of spp1
                const arma::mat & missing, // missing values
                const List & betas_det_1,// just a list, each element is a vector for each group
                const arma::vec & alpha_det_1,// similarly intercept
                const List & betas_det_2,// just a list, each element is a vector for each group
                const arma::vec & alpha_det_2,// similarly intercept
                const double & gamma_d, // associations between two spps
                const arma::vec & coding, // just the coding 
                int n_period,// number of periods 
                int n_site, // number of sites
                int G){ // number of groups

    arma::mat thr1 = getp(Designs_det, betas_det_1, 
                     alpha_det_1, n_period, n_site, G);

    thr1 = logit<arma::mat>(thr1); 

    arma::mat thr2 = getp(Designs_det, betas_det_2, 
                     alpha_det_2, n_period, n_site, G);

    thr2 = logit<arma::mat>(thr2); 

    arma::mat pmat(n_site,n_period, arma::fill::zeros);
    arma::vec p_temp(n_site, arma::fill::zeros);
    arma::vec partition(n_site, arma::fill::zeros); 
    List design_temp;

    for(int i = 0; i < n_period; i++){
        design_temp = Designs_det[i];
        // Hamiltonian:
        p_temp = thr1.col(i).t() * Y_1.col(i) + 
                thr2.col(i).t() * Y_2.col(i) + 
                gamma_d * (Y_1.col(i).t() * Y_2.col(i));
        // partition function 
        partition = getpartition_MRF(design_temp, 
                betas_det_1, alpha_det_1, 
                betas_det_2, alpha_det_2,
                gamma_d, coding, n_site, G);
        pmat.col(i) = log(exp(p_temp)/partition) * (1 - missing.col(i));
    }
    arma::vec res = sum(pmat,1);
    return(exp(res));

}

// calculate detection rate when only one of the two species exist. 
// [[Rcpp::export]]
arma::vec getp_01(const List & Designs_det,// should be a list of list, each element for a period
                const arma::mat & Y_pre,// Y for the present species
                const arma::mat & missing, 
                const List & betas_det_pre,// just a list, each element is a vector for each group
                const arma::vec & alpha_det_pre,// similarly intercept
                int n_period,// number of periods 
                int n_site, // number of sites
                int G // number of groups
                ){
    arma::mat p_temp = getp(Designs_det,// should be a list of list, each element for a period
                betas_det_pre,// just a list, each element is a vector for each group
                alpha_det_pre,// similarly intercept
                n_period,// number of periods 
                n_site, // number of sites
                G // number of groups
                );
    arma::mat pmat = log(p_temp) * Y_pre + log(1-p_temp) * (1-Y_pre);
    pmat *= (1-missing);
    arma::vec res = sum(pmat, 1);
    return(exp(res));
}

// [[Rcpp::export]]
double logLik_2spp_MRF(const arma::mat & psi_MRF, // order: 11, 10, 01, 00
              const arma::vec & p_11, // current psi and p, parameter
              const arma::vec & p_10, 
              const arma::vec & p_01, // current psi and p, parameter
              const arma::uvec & non_det_1, // 01
              const arma::uvec & non_det_2, // 10
              const arma::uvec & non_det_both){ // 00

    arma::vec p_site; 
    arma::vec psi01, psi10, psi00;
    psi01 = psi_MRF.col(2);
    psi10 = psi_MRF.col(1);
    psi00 = psi_MRF.col(3);
    p_site = psi_MRF.col(0) % p_11;
    p_site(non_det_1) += psi01(non_det_1) % p_01(non_det_1);
    p_site(non_det_2) += psi10(non_det_2) % p_10(non_det_2);
    p_site(non_det_both) += (psi01(non_det_both) % p_01(non_det_both)+
                            psi10(non_det_both) % p_10(non_det_both) + 
                            psi00(non_det_both));

    return(accu(log(p_site)));
}
