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