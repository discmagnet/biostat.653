# Problem 2
library(tidyverse)
library(reshape2)
library(foreign)

set.seed(20181001)
m <- 8             # 8 subjects.
n <- 3             # 3 measurements per subject.
sigma_squared <- 4 # variance.
rho <- c(1:9)/10   # true correlation.
eB0 <- vector()
eB1 <- vector()
for(j in 1:length(rho)){
  V0 <- matrix(c(1,        rho[j], rho[j]^2,
                 rho[j],   1,      rho[j],
                 rho[j]^2, rho[j], 1),
               nrow=3,ncol=3,byrow=TRUE) 
  
  # Create the design matrix, X
  trt_sequence <- expand.grid(c(0, 1), c(0, 1), c(0, 1))
  X <- cbind(1, c(t(trt_sequence)))
  beta0 <- 1
  beta1 <- 2
  mu_all <- X %*% matrix(c(beta0, beta1), ncol = 1)
  
  # Generate random Gaussian variables
  Y_vec <- vector()
  for (i in 1:(n*m)){
    Y_vec[i] <- rnorm(1,mu_all[i],sqrt(sigma_squared))
  }
  Y_mat <- matrix(Y_vec,nrow = 8,ncol = 3,byrow = T)
  
  # Determine Y and V (true variance)
  # Copied and pasted from class example
  res_svd <- svd(V0)
  V0_sqrt <- res_svd$u%*%diag(sqrt(res_svd$d))%*%t(res_svd$v)
  Y_mat_correlated <- Y_mat%*%V0_sqrt
  Y <- c(Y_mat_correlated)
  I_m_by_m <-  diag(rep(1,m))
  V <- I_m_by_m%x%V0
  
  # Functions to calculate, respectively:
  #  - weighted least squares
  #  - ordinary least squares
  #  - sampling variance of wls estimator
  #  - sampling variance of ols estimator
  wls <- function(y,X,W){
    A <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
    A%*%y
  }
  ols <- function(y,X){
    A <- solve(t(X)%*%X)%*%t(X)
    A%*%y
  }
  var_wls<- function(X,W,true_V,true_sigma_squared){
    A <- solve(t(X)%*%W%*%X)%*%t(X)%*%W
    A%*%(true_sigma_squared*true_V)%*%t(A)
  }
  var_ols<- function(X,true_V,true_sigma_squared){
    A <- solve(t(X)%*%X)%*%t(X)
    A%*%(true_sigma_squared*true_V)%*%t(A)
  }
  
  # Calculate Relative Efficiencies
  W <- solve(V)
  OLSv <- var_ols(X,V,sigma_squared)
  WLSv <- var_wls(X,W,V,sigma_squared)
  eB0[j] <- WLSv[1,1]/OLSv[1,1]
  eB1[j] <- WLSv[2,2]/OLSv[2,2]
}

# Plot Relative Efficiencies
# eB0
plot01 <- plot(rho,eB0,col="green")
plot01 <- lines(rho,eB0,col="green")
plot01
# eB1
plot02 <- plot(rho,eB1,col="red")
plot02 <- lines(rho,eB1,col="red")
plot02