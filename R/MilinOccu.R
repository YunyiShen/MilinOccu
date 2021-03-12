minlinoccu_EM <- function(Y, Designs_occu, Designs_det, 
                max_iter=10000, tol=1e-5, ...) {
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

    for(i in 1:max_iter){
        form_par_t <- formatpara(par_t, p_occu, p_det, G_occu, G_det) # formatted
        print(form_par_t)
        par_opt <- optim(par_t, Q_theta_theta_t, par_t = form_par_t, 
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
        L2_diff <- sum((par_opt$par-par_t)^2)
        if(L2_diff<=tol) {
            break
        }
        par_t <- par_opt$par
        cat(i,"th iter with L2 diff: ", signif(L2_diff,3),"\n")
    }

    par_res <- formatpara(par_opt$par, p_occu, p_det, G_occu, G_det)

    return(list(formatted = par_res, optim = par_opt))
}

minlinoccu_dir <- function(Y, Designs_occu, Designs_det, 
                max_iter=10000, tol=1e-5, ...) {
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
        
    par_t <- par_opt$par
    par_res <- formatpara(par_opt$par, p_occu, p_det, G_occu, G_det)
    return(list(formatted = par_res, optim = par_opt))
}