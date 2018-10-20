# --------------------------------------------------------------
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

# --------------------------------------------------------------
# Problem 3

set.seed(20181001)
m <- 40            # 40 subjects.
n <- 3             # 3 measurements per subject.
sigma_squared <- 4 # variance.
rho <- 0.5         # true correlation.
V0 <- matrix(c(1,     rho, rho^2,
               rho,   1,   rho,
               rho^2, rho, 1),
             nrow=3,ncol=3,byrow=TRUE) 

# Create the design matrix, X
trt_sequence <- expand.grid(c(0, 1), c(0, 1), c(0, 1))
X <- cbind(1, c(t(trt_sequence)))
X <- rbind(X,X,X,X,X)
beta0 <- 1
beta1 <- 2
mu_all <- X %*% matrix(c(beta0, beta1), ncol = 1)

WLSv1 <- vector()
WLSv2 <- vector()
WLSv3 <- vector()
WLSv4 <- vector()
WLSv5 <- vector()
WLSv6 <- vector()
NREP <- 1000
for(k in 1:NREP){
  # Generate random Gaussian variables
  Y_vec <- vector()
  for (i in 1:(n*m)){
    Y_vec[i] <- rnorm(1,mu_all[i],sqrt(sigma_squared))
  }
  Y_mat <- matrix(Y_vec,nrow = m,ncol = n,byrow = T)
  
  # Determine Y and V (true variance)
  # Copied and pasted from class example
  res_svd <- svd(V0)
  V0_sqrt <- res_svd$u%*%diag(sqrt(res_svd$d))%*%t(res_svd$v)
  Y_mat_correlated <- Y_mat%*%V0_sqrt
  Y <- c(Y_mat_correlated)
  I_m_by_m <-  diag(rep(1,m))
  V_e <- var(Y_mat_correlated)
  V_emp <- I_m_by_m%x%V_e
  
  # Weight Matrix 1:
  rho12 <- 0.5
  rho13 <- 0.8
  rho23 <- 0.2
  V1 <- matrix(c(1,rho12,rho13,
                 rho12,1,rho23,
                 rho13,rho23,1),
               nrow=3,ncol=3,byrow=TRUE)
  
  W1 <- solve(I_m_by_m%x%V1)
  
  # Weight Matrix 2:
  
  W2 <- diag(n*m) # inverse of identity is still identity
  
  # Weight Matrix 3:
  rho <- 0.5
  V3 <- matrix(c(1,     rho, rho^2,
                 rho,   1,   rho,
                 rho^2, rho, 1),
               nrow=3,ncol=3,byrow=TRUE)
  
  W3 <- solve(I_m_by_m%x%V3)
  
  # Calculate WLS Estimates
  WLSv1[k] <- var_wls(X,W1,V_emp,1)[1,1]
  WLSv2[k] <- var_wls(X,W2,V_emp,1)[1,1]
  WLSv3[k] <- var_wls(X,W3,V_emp,1)[1,1]
  WLSv4[k] <- var_wls(X,W1,V_emp,1)[2,2]
  WLSv5[k] <- var_wls(X,W2,V_emp,1)[2,2]
  WLSv6[k] <- var_wls(X,W3,V_emp,1)[2,2]
}

# Sandwich Estimates
var_w1_b0_emp <- mean(WLSv1)
var_w2_b0_emp <- mean(WLSv2)
var_w3_b0_emp <- mean(WLSv3)
var_w1_b1_emp <- mean(WLSv4)
var_w2_b1_emp <- mean(WLSv5)
var_w3_b1_emp <- mean(WLSv6)

V <- I_m_by_m%x%V0

# Theoretical Estimates
var_w1_b0_real <- var_wls(X,W1,V,sigma_squared)[1,1]
var_w2_b0_real <- var_wls(X,W2,V,sigma_squared)[1,1]
var_w3_b0_real <- var_wls(X,W3,V,sigma_squared)[1,1]
var_w1_b1_real <- var_wls(X,W1,V,sigma_squared)[2,2]
var_w2_b1_real <- var_wls(X,W2,V,sigma_squared)[2,2]
var_w3_b1_real <- var_wls(X,W3,V,sigma_squared)[2,2]

# --------------------------------------------------------------
# Problem 5

# 5.1.1 - Read the data
library(haven)
chol <- read_dta("~/WORKING_DIRECTORIES/biostat.653/cholesterol.dta")

# 5.1.2 - Calculate the sample means, standard deviations,
#         and variances of the serum cholesterol levels at each
#         occasion for each treatment group.

grp1 <- subset(chol, group == 1)
grp2 <- subset(chol, group == 2)

library(dplyr)
summary01 <- summarise(grp1, mean_y1 = mean(y1,na.rm = T),
                       mean_y2 = mean(y2,na.rm = T),
                       mean_y3 = mean(y3,na.rm = T),
                       mean_y4 = mean(y4,na.rm = T),
                       mean_y5 = mean(y5,na.rm = T),
                       var_y1 = var(y1,na.rm = T),
                       var_y2 = var(y2,na.rm = T),
                       var_y3 = var(y3,na.rm = T),
                       var_y4 = var(y4,na.rm = T),
                       var_y5 = var(y5,na.rm = T),
                       std_y1 = sqrt(var_y1),
                       std_y2 = sqrt(var_y2),
                       std_y3 = sqrt(var_y3),
                       std_y4 = sqrt(var_y4),
                       std_y5 = sqrt(var_y5))
summary02 <- summarise(grp2, mean_y1 = mean(y1,na.rm = T),
                       mean_y2 = mean(y2,na.rm = T),
                       mean_y3 = mean(y3,na.rm = T),
                       mean_y4 = mean(y4,na.rm = T),
                       mean_y5 = mean(y5,na.rm = T),
                       var_y1 = var(y1,na.rm = T),
                       var_y2 = var(y2,na.rm = T),
                       var_y3 = var(y3,na.rm = T),
                       var_y4 = var(y4,na.rm = T),
                       var_y5 = var(y5,na.rm = T),
                       std_y1 = sqrt(var_y1),
                       std_y2 = sqrt(var_y2),
                       std_y3 = sqrt(var_y3),
                       std_y4 = sqrt(var_y4),
                       std_y5 = sqrt(var_y5))

# 5.1.3 - On a single graph, construct a time plot that displays
#         the mean serum cholesterol versus time (in months) for
#         the two treatment groups. Describe the general chacter-
#         istics of the time trends for the two groups.

library(ggplot2)
data513 <- data.frame(rbind(t(summary01[,1:5]),t(summary02[,1:5])),
                      c(rep(1,5),rep(2,5)),c(0:4,0:4)*6)
colnames(data513) <- c("mean","group","index")
plot513 <- ggplot(data = data513,
                  aes(x = index, y = mean, color = as.factor(group))) +
  geom_line() +
  xlab("Follow-up Time (in months)") +
  ylab("Mean Serum Cholesterol")
plot513

# 5.1.4 - Put the data in "long" format, with 5 records per subject

library(reshape2)
data514 <- melt(chol[3:7])
data514 <- mutate(data514, group = rep((c(rep(1,62),rep(2,41))),5))
data514 <- mutate(data514, id = rep(c(1:103),5))
data514 <- arrange(data514,id)

# 5.1.5 - Assuming an unstructured covariance matrix, conduct an
#         analysis of response profiles. Determine whether the patterns
#         of change over time differ in the two treatment groups.

library(nlme)
res.gls <- gls(model = value ~ factor(variable) * factor(group, levels = c(2,1)),
               data = data514,
               # Covariance structure: unstructured
               correlation = corSymm(form = ~ as.numeric(factor(variable))|id),
               weights = varIdent(form = ~ 1 | factor(variable)),
               na.action = "na.omit",
               method = "ML")
summary(res.gls)

drop1(res.gls, test = "Chisq")

# 5.1.6 - Display the estimated 5 x 5 covariance and correlation matrices
#         for the five repeated measurements of serum cholesterol.

getVarCov(res.gls) # covariance matrix
# correlation matrix in summary(res.gls)