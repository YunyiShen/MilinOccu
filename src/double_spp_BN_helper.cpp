#include "single_spp_helper.h"

// occupancy rate of the subordinate species
arma::vec getpsi_sub(const List & Designs, // a list of designs for each groups
                const List & betas, // list of regression coefficients 
                const arma::vec & alpha, // the intercept, should be size G
                const arma::vec & omega, // two dimensional competing parameter, intercept and slope
                int n_site, // number of sites, basecally number of rows of Designs[[1]]
                int G){// number of groups

    arma::mat linPred(n_site, G+1, arma::fill::zeros);// the linear predictors (plus intercept)
    //arma::mat design_temp;
    //arma::vec beta_temp;
    for (int i = 0; i < G; i++)
    {
        arma::mat design_temp = Designs[i];
        arma::vec beta_temp = betas[i];
        linPred.col(i) = design_temp * beta_temp;
    }

    linPred.cols(0,G-1).each_row() += alpha.t(); 
    arma::mat linPred_ab = linPred; 
    linPred.col(G) += omega(0) + omega(1); // when dominate species exist
    linPred_ab.col(G) += omega(0); //when it does not exist

    arma::mat psi(n_site,2);
    psi.col(1) = arma::min(linPred,1);// take minimum 
    psi.col(0) = arma::min(linPred_ab,1); 
    psi = 1/(1+exp(-psi));
    return(psi);
}

// this calculate the probability of the detection rate when subordinate species exist, two columns are whether the dominate exist (first column absent and second exist)
arma::mat p_det_sub(const List & Designs_det,// should be a list of list, each element for a period
                const arma::mat & Y_sub, // detection 
                const arma::mat & missing, // missing values
                const List & betas_det,// just a list, each element is a vector for each group
                const arma::vec & alpha_det,// similarly intercept
                const arma::vec & omega, // effect on detection 
                int n_period,// number of periods 
                int n_site, // number of sites
                int G // number of groups
                ){
    arma::mat pmat_z10(n_site,n_period, arma::fill::zeros);
    arma::mat pmat_z11(n_site,n_period, arma::fill::zeros);
    arma::mat ptemp;
    List design_temp;
    arma::mat res(n_site, 2);

    for (int i = 0; i < n_period; i++)
    {
        design_temp = Designs_det[i];
        ptemp = getpsi_sub(design_temp, betas_det, alpha_det, omega,n_site, G);
        pmat_z10.col(i) = ptemp.col(0); // p when dominate absent
        pmat_z11.col(i) = ptemp.col(1);
    }
    res.col(0) = arma::sum(((1-missing) % Y_sub) % log(pmat_z10+1e-10),1) + 
                 arma::sum(((1-missing) % (1-Y_sub)) % log((1-pmat_z10)+1e-10),1);
    res.col(1) = arma::sum(((1-missing) % Y_sub) % log(pmat_z11+1e-10),1) + 
                 arma::sum(((1-missing) % (1-Y_sub)) % log((1-pmat_z11)+1e-10),1);
    return(exp(res));
}


arma::mat p_det_dom(const List & Designs_det,// should be a list of list, each element for a period
                const arma::mat & Y_dom, // detection 
                const arma::mat & missing, // missing values
                const List & betas_det,// just a list, each element is a vector for each group
                const arma::vec & alpha_det,// similarly intercept
                int n_period,// number of periods 
                int n_site, // number of sites
                int G // number of groups
                ){
    arma::mat p = getp(Designs_det,// should be a list of list, each element for a period
                betas_det,// just a list, each element is a vector for each group
                alpha_det,// similarly intercept
                n_period,// number of periods 
                n_site, // number of sites
                G // number of groups
                );

    arma::vec res = arma::sum(((1-missing) % Y_dom) % log(p+1e-10),1) + 
                 arma::sum(((1-missing) % (1-Y_dom)) % log((1-p)+1e-10),1);

    return(exp(res));

}


double logLik_2spp_BN(const arma::vec & psi_dom, const arma::mat & p_dom, // current psi and p, parameter
              const arma::mat & psi_sub, const arma::mat & p_sub, // current psi and p, parameter
              const arma::uvec & non_det_dom, 
              const arma::uvec & non_det_sub,
              const arma::uvec & non_det_both){
    // two species all absent, detection is not a problem
    arma::vec p00 = (1-psi_dom) % (1-psi_sub.col(0));
    // dominate exist but not subordinate:
    arma::vec p10 = psi_dom % p_dom % (1-psi_sub.col(1));
    // subordinate exist but not dominate
    arma::vec p01 = (1 - psi_dom) % (psi_sub.col(1)) % p_sub.col(1);
    // both exist
    arma::vec p11 = psi_dom % p_dom % (psi_sub.col(1)) % p_sub.col(1);

    // if see both, can only be p11
    arma::vec p_site = p11;

    // did not see dominate but see subordinate can be p11 or p01
    p_site(non_det_dom-1) = p11(non_det_dom-1) + p01(non_det_dom-1);
    // did not see sub but se dom, 11 or 10
    p_site(non_det_sub-1) = p11(non_det_sub-1) + p10(non_det_sub-1);
    // did not see both, can be anywhere
    p_site(non_det_both-1) = p11(non_det_both-1) + p01(non_det_both-1) + 
                           p10(non_det_both-1) + p00(non_det_both-1);

    return(accu(log(p_site)));

}

