minlinoccu_dir <- function(Y, Designs_occu, Designs_det, ...) {
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


boot_helper <- function(i,formattedpar, Designs_occu, Designs_det,
                    p_occu, p_det, G_occu, G_det, # dimension informations
                  n_site, n_period){

    boot_data <- simuY(formattedpar, Designs_occu, Designs_det, ## design lists, 
                  p_occu, p_det, G_occu, G_det, # dimension informations
                  n_site, n_period)
    minlinoccu_dir(boot_data$det, Designs_occu, Designs_det)

}

parabootstrap <- function(n=100, formattedpar, Designs_occu, Designs_det){
    G_occu <- length(Designs_occu) # groups
    G_det <- length(Designs_det[[1]])
    n_period <- length(Designs_det) # period
    n_site <- nrow(Designs_occu[[1]]) # sites

    p_occu <- sapply(Designs_occu, ncol ) # number of predictors in each group
    p_det <- sapply(Designs_det[[1]], ncol )

    lapply(1:n, boot_helper, formattedpar, Designs_occu, Designs_det,
                    p_occu, p_det, G_occu, G_det, # dimension informations
                    n_site, n_period)
}