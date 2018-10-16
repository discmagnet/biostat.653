# Problem 2
library(tidyverse)
library(reshape2)
library(foreign)

set.seed(20181001)
m <- 8 # 8 subjects.
n <- 3 # 3 measurements per subject.
sigma_squared <- 4 # variance.
rho <- 0.4 # true correlation.
V0 <- matrix(c(1,rho,rho^2,
               rho,1,rho,
               rho^2,rho,1),
             nrow=3,ncol=3,byrow=TRUE) 


# design matrix:
trt_sequence <- expand.grid(c(0, 1), c(0, 1), c(0, 1))
X <- cbind(1, c(t(trt_sequence)))
beta0 <- 1
beta1 <- 2
mu_all <- X %*% matrix(c(beta0, beta1), ncol = 1)

# generate random Gaussian variables
Y_vec <- vector()
for (i in 1:(n*m)){
  Y_vec[i] <- rnorm(1,mu_all[i],sqrt(sigma_squared))
}
Y_mat <- matrix(Y_vec,nrow = 8,ncol = 3,byrow = T)


res_svd <- svd(V0) #singular value decomposition.
# take square root of the matrix in the middle:
V0_sqrt <- res_svd$u%*%diag(sqrt(res_svd$d))%*%t(res_svd$v) 
round(V0_sqrt%*%V0_sqrt,2) #round so that very small numbers are rounded to zero.
# pre-multiply each vector of independent Gaussian variables Y_mat[,i]
# to produce correlated measurements:
Y_mat_correlated <- Y_mat%*%V0_sqrt
# response variables from all m=8 subjects:
Y <- c(Y_mat_correlated) # stacking by column.
I_m_by_m <-  diag(rep(1,m))
V <- I_m_by_m%x%V0 # variance covariance matrix for all subjects.

# function to calculate weighted least squares:
wls <- function(y,X,W){
  A <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
  A%*%y
}

# function to calculate ordinary least squares:
ols <- function(y,X){
  A <- solve(t(X)%*%X)%*%t(X)
  A%*%y
}

# function to compute the sampling variance of wls estimator:
var_wls<- function(X,W,true_V,true_sigma_squared){
  A <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
  A%*%(true_sigma_squared*true_V)%*%t(A)
}

# function to compute the sampling variance of ols estimator:
var_ols<- function(X,true_V,true_sigma_squared){
  A <- solve(t(X)%*%X)%*%t(X)
  A%*%(true_sigma_squared*true_V)%*%t(A)
}

W <- solve(V)
ORLv <- var_ols(X,V,sigma_squared)
WLSv <- var_wls(X,W,V,sigma_squared)
