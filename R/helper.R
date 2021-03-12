# we need a convention on par in optimization algorithms, 
#  here we will first have occupancy alpha, then beta, then det alpha, det beta
#  this code is used to formate the par so that we can calculate 
#  the Q(\theta|\theta_t)
#  used in EM algorith written in C++
#
# Conventions on formating, 
#   par: is a vector, parameters, should have deminsion: 
#       sum(p_occu) + sum(p_det) + G_occu + G_det 
#   p_occu: is a vector, with G_occu elements, 
#       are number of predictors in each group
#   p_det: a vector with G_occu elements
#   G_occu: integer, number of groups of occupancy
#   G_det : integer number of groups of detections

formatpara <- function(par, p_occu, p_det, G_occu, G_det){
    beta_det <- list()
    beta_occu <- list()
    g <- 0
    alpha_occu <- par[1:G_occu]
    g <- g + G_occu
    for(i in 1:G_occu){
        beta_occu[[i]] <- par[1:p_occu[i] + g]
        g <- g + p_occu[i]
    }
    alpha_det <- par[1:G_det + g]
    for(i in 1:G_det){
        beta_det[[i]] <- par[1:p_det[i] + g]
        g <- g + p_det[i]
    }

    return(list(alpha_occu = alpha_occu, beta_occu = beta_occu, 
           alpha_det = alpha_det, beta_det = beta_det))
}

# par_t is the formatted par_t after formatpara()
#  Design_occu is a list with G_occu design matrix of each group of predictors,
#    They should all have same number of rows
#  Design_det is a list whose element has same format as Design_occu, each element
#    is for one period. 


Q_theta_theta_t <- function(par, par_t, 
                            Y, missing, # data and flag of missing
                            non_det, # indeces of non detections
                            Designs_occu, Designs_det, ## design lists, 
                            p_occu, p_det, G_occu, G_det, # dimension informations
                            n_site, n_period){ # repeat informations

    par_new <- formatpara(par, p_occu, p_det, G_occu, G_det)
    psi_t <- getpsi(Designs_occu, par_t$beta_occu, par_t$alpha_occu,
                    n_site, G_occu)
    psi <- getpsi(Designs_occu, par_new$beta_occu, par_new$alpha_occu,
                    n_site, G_occu)
    p_t <- getp(Designs_det, par_t$beta_det, par_t$alpha_det,
                    n_period, n_site, G_det)
    p <- getp(Designs_det, par_new$beta_det, par_new$alpha_det,
                    n_period, n_site, G_det)
    #browser()
    #cat("still alive\n")
    -Estep_cpp(psi,p,psi_t,p_t,Y, missing, non_det)
}


logLikR <- function(par,  
                            Y, missing, # data and flag of missing
                            non_det, # indeces of non detections
                            Designs_occu, Designs_det, ## design lists, 
                            p_occu, p_det, G_occu, G_det, # dimension informations
                            n_site, n_period){ # repeat informations

    par_new <- formatpara(par, p_occu, p_det, G_occu, G_det)
    psi <- getpsi(Designs_occu, par_new$beta_occu, par_new$alpha_occu,
                    n_site, G_occu)
    p <- getp(Designs_det, par_new$beta_det, par_new$alpha_det,
                    n_period, n_site, G_det)
    #browser()
    #cat("still alive\n")
    -logLik_cpp(psi,p,Y, missing, non_det)
}



# simulate detection history, used in 
simuY <- function(formatedpar, Designs_occu, Designs_det, ## design lists, 
                  p_occu, p_det, G_occu, G_det, # dimension informations
                  n_site, n_period){

    psi <- getpsi(Designs_occu, formatedpar$beta_occu, 
                  formatedpar$alpha_occu,
                  n_site, G_occu)

    p <- getp(Designs_det, formatedpar$beta_det, 
              formatedpar$alpha_det,
              n_period, n_site, G_det)

    occu <- runif(n_site)<=psi
    det <- matrix(runif(n_site * n_period),n_site) <= p

    det * (occu %*% matrix(1,nrow = 1, ncol = n_period))
}

