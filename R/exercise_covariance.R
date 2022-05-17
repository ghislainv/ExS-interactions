# Load libraries
library("Rcpp")
library("RcppArmadillo")
library("here")

# Source C++ code for function rmvnArma
sourceCpp(here("src/rmvnArma.cpp"))

# Param
E <- 100 # Environmental dimensions
K <- 100 # Number of sites
S <- 20  # Number of species
varBeta <- 0.1
varX <- 0.1

# Species parameters
beta <- matrix(rnorm(E * S, 0, sqrt(varBeta)), E, S)     # X, B, V
# Get covariance matrix V for environmental dimensions X
# Must be symmetric, positive semi-definite
V_E <- cov(rmvnArma(E + 1, rep(0, E), varX * diag(E)))
# Draw environmental dimensions for each site
x <- rmvnArma(K, rep(0, E), V_E)
V <- cov(x)
# Species covariance
C <- t(beta) %*% V %*% beta                             # C = B'VB
xb <- x %*% beta 

print(diag(V))