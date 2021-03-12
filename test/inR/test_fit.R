library(Rcpp)
library(RcppArmadillo)
library(Matrix)

source("./R//helper.R")
sourceCpp("./src/EM_helper.cpp")
source("./R/MilinOccu.R")


# simulate data
# par formating
G_occu <- 3
G_det <- 2
p_occu <- 1:G_occu
p_det <- 1:G_det

set.seed(42)
par <- rnorm(G_occu + G_det + sum(p_occu) + sum(p_det), 3,2)

theta <- formatpara(par,p_occu,p_det,G_occu, G_det)


n_site <- 500
n_period <- 15

Designs <- list(matrix(rnorm(n_site*p_det[1]),n_site),
    matrix(rnorm(n_site * p_det[2]),n_site))

Designs_occu <- list(matrix(rnorm(n_site*p_occu[1]),n_site),
    matrix(rnorm(n_site * p_occu[2]),n_site),
    matrix(rnorm(n_site * p_occu[3]),n_site)
    )


Designs_det <- lapply(1:n_period, function(i,w){w}, Designs)

psi <- getpsi(Designs_occu, theta$beta_occu, theta$alpha_occu, n_site, G_occu)
hist(psi)

Y <- simuY(theta, Designs_occu,Designs_det,
      p_occu, p_det, G_occu, G_det, # dimension informations
      n_site, n_period)

tryres <- minlinoccu_dir(Y, Designs_occu, Designs_det, 
         control = list(trace = 1, maxit = 10000))

