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
n_site <- 10
Desings <- list(matrix(rnorm(n_site*p_det[1]),n_site),
    matrix(rnorm(n_site * p_det[2]),n_site))



Desings_occu <- list(matrix(rnorm(n_site*p_occu[1]),n_site),
    matrix(rnorm(n_site * p_occu[2]),n_site),
    matrix(rnorm(n_site * p_occu[3]),n_site)
    )


psi <- getpsi(Desings_occu, fooo$beta_occu, fooo$alpha_occu, n_site, G_occu)

# get p
n_period <- 10
Desings_det <- lapply(1:n_period, function(i,w){w}, Desings)

p <- getp(Desings_det,fooo$beta_det, fooo$alpha_det, n_period,n_site, G_det)


# Q(\theta|\theta_t)
Y <- 1 * (matrix(rnorm(n_site*n_period),n_site)>0)
Y[1,] <- 0
missing <- 0*Y
non_det <-c(1) 


Estep_cpp(psi,p,psi,p, Y, missing, non_det)

Q_theta_theta_t(par, fooo, 
                            Y, missing, # data and flag of missing
                            non_det, # indeces of non detections
                            Desings_occu, Designs_det, ## design lists, 
                            p_occu, p_det, G_occu, G_det, # dimension informations
                            n_site, n_period)



simuY(fooo, Desings_occu,Desings_det,
      p_occu, p_det, G_occu, G_det, # dimension informations
      n_site, n_period)
