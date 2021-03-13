#' Single Season Occupancy Model with Min-Linear Logistic Regression
#' 
#' @description Fit a single season occupancy model based on min-linear Logistic regression.
#'
#' @param Y The detection history matrix, with row as sites and column as repeats, missing data should be NA
#' @param Designs_occu A list with each element is a design matrix without intercept, different element can have overlap columns but should not be identical
#' @param Designs_det A list of design matrix list. Each element represent the design matrix list of repeat, and the design matrix list should have the same structure like Design_occu. Structure should be consistant among repeats. 
#' @param ... Extra arguments send to optim 
#' @return A list,`$formatted` is the formatted parameter estimation and `$optim` is the raw output of optim. 
#' 
#' @details Instead of logistic regression, this model fits a min-linear logistic based occupancy, i.e. logit(psi) and logit(p) are not linear to predictors but is the minimum among several linear predictors, which represent the worst environment.  

minlinoccu <- function(Y, Designs_occu, Designs_det, ...) {
    G_occu <- length(Designs_occu) # groups
    G_det <- length(Designs_det[[1]])
    n_period <- length(Designs_det) # period
    n_site <- nrow(Designs_occu[[1]]) # sites

    p_occu <- sapply(Designs_occu, ncol ) # number of predictors in each group
    p_det <- sapply(Designs_det[[1]], ncol )

    missing <- 1 * (is.na(Y)) # which is missing
    Y[is.na(Y)] <- 0 # make all NAs 0 to avoid numeric issue
    non_det <- which (rowSums(Y)==0) # no detections
    par_t <- rnorm(G_occu + G_det + sum(p_occu) + sum(p_det))
    par_opt <- optim(par_t, logLikR,  
                            Y=Y, missing = missing, # data and flag of missing
                            non_det = non_det, # indeces of non detections
                            Designs_occu = Designs_occu, 
                            Designs_det = Designs_det, ## design lists, 
                            p_occu = p_occu, 
                            p_det = p_det, 
                            G_occu = G_occu, 
                            G_det = G_det, # dimension informations
                            n_site = n_site, 
                            n_period = n_period, ...)
        
    par_res <- formatpara(par_opt$par, p_occu, p_det, G_occu, G_det)
    return(list(formatted = par_res, optim = par_opt))
}


boot_helper <- function(i, Y,formattedpar, Designs_occu, Designs_det,
                    p_occu, p_det, G_occu, G_det, # dimension informations
                  n_site, n_period, ...){

    boot_data <- simuY(formattedpar, Designs_occu, Designs_det, ## design lists, 
                  p_occu, p_det, G_occu, G_det, # dimension informations
                  n_site, n_period)
    boot_data$det[is.na(Y)] <- NA
    minlinoccu(boot_data$det, Designs_occu, Designs_det,...)

}


#' Parametric Bootstrap for milinoccu
#' 
#' This function run n paramatric bootstrap to help construct confidence intervals
#' 
#' @param n Bootstrap sample size
#' @param Y The original data, used to determine missing
#' @param formattedpar The formatted parameter estimation 
#' @param Designs_occu Design list for occupancy, see milinoccu
#' @param Designs_det Design list for detection, see milinoccu
#' @param ... Extra arguments send to optim 
#' 
#' @return A list, each element is an output list of milinoccu

milinbootstrap <- function(n=100, Y,formattedpar, Designs_occu, Designs_det, ...){
    G_occu <- length(Designs_occu) # groups
    G_det <- length(Designs_det[[1]])
    n_period <- length(Designs_det) # period
    n_site <- nrow(Designs_occu[[1]]) # sites

    p_occu <- sapply(Designs_occu, ncol ) # number of predictors in each group
    p_det <- sapply(Designs_det[[1]], ncol )

    lapply(1:n, boot_helper, Y,formattedpar, Designs_occu, Designs_det,
                    p_occu, p_det, G_occu, G_det, # dimension informations
                    n_site, n_period, ...)
}