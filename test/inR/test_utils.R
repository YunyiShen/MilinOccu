library(Rcpp)
library(RcppArmadillo)
library(Matrix)

source("./R//helper.R")
sourceCpp("./src/EM_helper.cpp")

# test par formating
G_occu <- 3
G_det <- 2
p_occu <- 1:G_occu
p_det <- 1:G_det

par <- rnorm(G_occu + G_det + sum(p_occu) + sum(p_det), 1,1)

fooo <- formatpara(par,p_occu,p_det,G_occu, G_det)

## test geting psi
n_site <- 100
Designs <- list(matrix(rnorm(n_site*p_det[1]),n_site),
    matrix(rnorm(n_site * p_det[2]),n_site))



Designs_occu <- list(matrix(rnorm(n_site*p_occu[1]),n_site),
    matrix(rnorm(n_site * p_occu[2]),n_site),
    matrix(rnorm(n_site * p_occu[3]),n_site)
    )


psi <- getpsi(Designs_occu, fooo$beta_occu, fooo$alpha_occu, n_site, G_occu)

# get p
n_period <- 10
Designs_det <- lapply(1:n_period, function(i,w){w}, Designs)

p <- getp(Designs_det,fooo$beta_det, fooo$alpha_det, n_period,n_site, G_det)


# Q(\theta|\theta_t)
Y <- simuY(fooo, Designs_occu,Designs_det,
      p_occu, p_det, G_occu, G_det, # dimension informations
      n_site, n_period)
missing <- 0*Y
non_det <- which(rowSums(Y)==0)


Estep_cpp(psi,p,psi,p, Y, missing, non_det)

Q_theta_theta_t(par, fooo, 
                            Y, missing, # data and flag of missing
                            non_det, # indeces of non detections
                            Designs_occu, Designs_det, ## design lists, 
                            p_occu, p_det, G_occu, G_det, # dimension informations
                            n_site, n_period)


logLik_cpp(psi,p, Y, missing, non_det)



